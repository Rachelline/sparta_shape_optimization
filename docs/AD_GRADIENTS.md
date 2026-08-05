# AD gradients in SPARTA: the flux-measure bias, and what's corrected

This document explains why forward-mode automatic differentiation (AD) of
a DSMC surface tally is *biased* — not noisy, not wrong-by-a-bug, but
systematically off by a specific, derivable factor — and what the runtime
`--score-correction` switch does and does not fix. It distills a longer
local investigation log (`docs/ad_phase_c_investigation/FINDINGS.md`,
gitignored, kept as a working-notes archive) into something meant for an
external reader who wants to *use* the AD build correctly, not re-derive
the finding.

If you only read one section, read "Coverage: what's fixed and what
isn't" below before trusting any AD gradient this codebase produces.

## The setup

SPARTA's AD build replaces `sfloat` (normally `double`) with
`Sacado::Fad::SFad<double, N>`, a forward-mode dual number that carries
`N` derivative directions alongside its value. Seed a design parameter's
derivative at the point it enters the simulation (e.g. a surface point's
coordinates, via `src/read_surf.cpp`'s `SPARTA_AD_SEED_JACFILE` hook), and
ordinary arithmetic propagates `d(anything)/d(alpha)` through the entire
DSMC pipeline — particle moves, collisions, surface tallies — via the
chain rule, automatically and exactly, *for anything that's a
deterministic function of the seeded quantities*.

The problem is that a DSMC surface tally is **not** a deterministic
function of the design parameters alone. It's an expectation over a random
process — which particles happen to hit the surface — and the *set of
particles that hit* depends on the design parameters too, through a
mechanism forward-mode AD cannot see.

## The mechanism: a missing score-function term

Consider a per-hit surface tally like pressure: each specular hit deposits
normal momentum `2m|u|`, where `u = v·n` is the incoming normal velocity
component. But hits aren't drawn uniformly from the freestream — a surface
element is struck at a rate proportional to the incoming flux `|u|`, so
the *population* of particles that register a hit at all carries density
`∝ |u| f(v)`. Both the sampling weight `|u|` and the per-hit tally `2m|u|`
depend on `alpha` (through the surface normal `n(alpha)`), but in a single
concrete simulation run, the set of particles that hit is already drawn
and frozen — AD only sees the second dependency, differentiating the
tally value with the hit set held fixed.

This is exactly the missing term in the standard identity for
differentiating an expectation under a parametrized measure:

```
d/dα E_p(α)[g(α)]  =  E[dg/dα]  +  E[g · dlog p/dα]
                       ^^^^^^^^ what pathwise/forward AD computes
                                  ^^^^^^^^^^^^^^^^^^ the score-function term
                                                       AD is structurally
                                                       blind to
```

For specular reflection, both the sampling-weight term and the tally term
turn out to be **exactly equal** (the integrand is quadratic in `|u|`,
each factor linear), so AD returns exactly **half** the true derivative —
confirmed independently via Monte Carlo (`measure_theory_check.py`, no
SPARTA needed, ratio 0.5000 to four decimals at every angle tested) and
against multiple SPARTA finite-difference/analytic cross-checks (worst
observed error after applying the known 2× factor: 2.4%, within DSMC
noise). Diffuse-wall reflection carries a different, larger bias (thermal
re-emission adds an `|u|`-independent piece to the tally, breaking the
exact-half symmetry) — not yet pinned to a single clean constant the way
specular reflection is.

**This is not a bug to "fix" in the sense of finding a wrong line of
code.** Forward-mode AD is mathematically incapable of seeing a
score-function term; only the tally-side term is representable as an
ordinary derivative through deterministic arithmetic. The fix has to
*inject* the missing term some other way.

## The fix: `u / spval(u)`

`src/compute_surf.cpp` multiplies each per-hit `fluxscale` factor by
`u / spval(u)`, where `u = MathExtra::dot3(vorig, norm)` and `spval()`
strips the derivative, returning a plain value:

```cpp
if (score_correction && iorig) {
  sfloat u = MathExtra::dot3(vorig, norm);
  fluxscale *= u / spval(u);
}
```

`u / spval(u)` is **numerically exactly 1.0** — every tallied *value* is
untouched, bit-for-bit — but Sacado's ordinary product rule on a factor
whose value is 1 and whose derivative slot is `dlog|u|/dα` (since
`d(u/spval(u))/dα = du/dα / u = dlog|u|/dα` at the point `spval(u)=u`)
means the product rule applied to every `fluxscale`-weighted tally
downstream now yields `dg/dα + g·dlog|u|/dα` — the corrected, score-
function-complete derivative — with no manual expectation bookkeeping
anywhere else in the code.

Toggle: `SPARTA_AD_SCORE_CORRECTION` env var, read once per `compute surf`
command in `ComputeSurf::init()`. `ad_src/`'s drivers wrap this in
`AdCorrections` (`ad_src/shape_case.h`) and expose it as
`--score-correction`/`--no-score-correction` — see `ad_src/README.md`'s
"Flipping the AD correction switch" section for the CLI-level view; this
document only covers the underlying mechanism.

## Coverage: what's fixed and what isn't

| Tally family | Scale factor | Corrected? | Used by |
|---|---|---|---|
| Pressure, shear, heat flux | `fluxscale` | **Yes** | `HeatFluxObjective` (`etot`) |
| Force / torque (`FX/FY/FZ`, `TX/TY/TZ`) | `nfactor_inverse` | **No — not yet implemented** | `DragObjective` (`fx`) |
| Diffuse-wall reflection resampling | (re-emission velocity re-draw) | Not analyzed | any diffuse-wall run |
| Multi-collision compounding | (repeated resampling across `collide vss` events) | Not analyzed | any collisional (non-free-molecular) run |

The practical consequence, confirmed by this codebase's own golden
captures: **`--score-correction` is a documented no-op for
`--objective drag`**, because `DragObjective` tallies `fx`, which is
scaled by `nfactor_inverse` in `compute_surf.cpp` — a code path the
`u/spval(u)` insertion never touches. It has a real, measurable effect for
`--objective heatflux` (`etot`, scaled by `fluxscale`) — running the same
alpha with and without the flag produces visibly different gradient
components, including sign flips on individual components (see
`tools/ad_verify/compare_modes.py`'s output for a live comparison).

Extending the correction to `nfactor_inverse`-weighted tallies is
tracked as deferred work — see `KNOWN_GAPS.md`. The mechanism would be
analogous (a local `u/spval(u)`-style factor near `compute_surf.cpp`'s
`nfactor_inverse` multiplications), but hasn't been derived or verified
yet, and per an earlier explicit scoping decision, changing `src/`
(non-`ad_src/`) code is handled as separate, later work.

## Practical guidance

- **Never assume an AD gradient's sign or magnitude is trustworthy without
  checking which tally it came from.** `DragObjective` gradients are
  always uncorrected, regardless of the flag. `HeatFluxObjective`
  gradients are correctable, but only when the flag is actually set —
  it's off by default.
- **The bias is a multiplicative, not additive, effect on a per-tally
  basis** — for specular pressure specifically it's a clean 2× factor,
  but do not generalize that exact number to other tallies, wall models,
  or reflection types. Diffuse walls and shear tallies carry different,
  less-pinned-down factors. Verify with `compare_modes.py` or a direct
  FD comparison before trusting a new configuration's AD gradient
  quantitatively.
- **The uncorrected gradient's *direction* is usually still useful.**
  Since the missing term is a same-sign multiplier in the cases analyzed
  so far (not a sign flip or an unrelated random error), gradient-descent-
  style optimization using the raw, uncorrected AD gradient still
  generally moves downhill — this codebase's own AD-driven `shape_opt_ad`
  optimization runs converge to sensible optima using uncorrected
  gradients. Don't expect the same trajectory or iteration count as a
  finite-difference-driven run, though; a uniformly-scaled (or
  differently-biased-per-component) gradient changes IPOPT's internal
  step-size and line-search behavior even when the direction is right.
- **Values are never affected**, corrected or not — `u/spval(u)`'s value
  is always exactly 1.0. Only derivatives change.

## Where the underlying investigation lives

`docs/ad_phase_c_investigation/FINDINGS.md` (gitignored, local working
notes) has the full investigation trail: the disproven hypotheses (grid
topology, finite-plate edge effects, statistical noise, corner
ambiguity), the independent hit-count-derivative cross-check, and an
earlier, narrower, now-fixed bug (`FINDING 1`: `fix_emit_face.cpp`
truncating a derivative-carrying particle count to an integer before AD
could see it). If you're extending the AD correction to a new tally
family, that file's methodology (isolate with a rotating/deterministic
geometry, cross-check against a zero-noise analytic reference, verify
with an independent measurement of the missing term) is worth following
even though the file itself isn't committed.
