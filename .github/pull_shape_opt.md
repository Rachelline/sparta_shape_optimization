## Purpose

_Briefly describe the new feature(s), enhancement(s), or bugfix(es) included in this pull request. If this addresses an open GitHub Issue, mention the issue number, e.g. with `fixes #221` or `closes #135`, so that issue will be automatically closed when the pull request is merged_

## Author(s)

_Please state name and affiliation of the author or authors that should be credited with the changes in this pull request_

## Backward Compatibility

_Please state whether any changes in the pull request break backward compatibility for inputs, and - if yes - explain what has been changed and why_

## Implementation Notes

_Provide any relevant details about how the changes are implemented, how correctness was verified, how other features - if any - in SPARTA are affected_

## ⚠️ AD / Sacado / Kokkos Build Warning

_This fork adds an automatic-differentiation build (`-DSPARTA_ENABLE_AD=ON`, `sfloat = Sacado::Fad::SFad`). If your change touches the AD build, the Kokkos build, or the Sacado/Trilinos dependency, read and confirm the following:_

- **One Kokkos only.** Sacado's Kokkos and SPARTA's `lib/kokkos` **must be the same version**. Two different Kokkos versions in one binary is an ODR violation → silent wrong results (especially on GPU). `SFad` is header-only, so it binds to whichever Kokkos is on SPARTA's include path at compile time — the AD+Kokkos config must therefore use **one** Kokkos (`USE_EXTERNAL_KOKKOS=ON` + `Kokkos_ROOT=<trilinos-install>`), never a second one.
- **Version pin is load-bearing.** Sacado must come from the Trilinos release whose bundled Kokkos == SPARTA's `lib/kokkos` version. As of this writing: SPARTA `lib/kokkos` = **5.0.2** → pin **`trilinos-release-17-0-0`** (its Kokkos is 5.0.2). Trilinos 16.2 ships Kokkos 4.7.4 (too old); 17.1.x ships 5.1.1 (newer). **Rule: match Sacado to Kokkos — bump Sacado, never downgrade Kokkos.**
- **If you bump `lib/kokkos`,** you MUST move the Trilinos pin to the release carrying the matching Kokkos, or the AD+Kokkos build breaks. The configure-time version guard in `cmake/common/process/sparta_build_options.cmake` should catch a mismatch and error out — do not disable it.
- **Minimal Trilinos only.** Build Sacado with `Trilinos_ENABLE_ALL_PACKAGES=OFF` + `Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF` + `Trilinos_ENABLE_Sacado=ON` (+ `Trilinos_ENABLE_Kokkos=ON` for the AD+Kokkos config). Never depend on the full Trilinos suite.
- **C++ standard.** Trilinos 17.0.0 builds at **C++20**; the AD build's minimum standard is at least C++17 (Kokkos 5.0.2 floor). Confirm the standard bump does not regress the non-AD build.
- **Confirm before merging:**
  - [ ] AD build (`-DSPARTA_ENABLE_AD=ON`) reproduces stock (non-AD) results within regression tolerance (`1.e-7`)
  - [ ] If AD+Kokkos: Sacado's Kokkos version == SPARTA's `lib/kokkos` version, and only one Kokkos is linked
  - [ ] The Trilinos pin / Kokkos version guard is still correct after this change

## Post Submission Checklist

_Please check the fields below as they are completed_
- [ ] The feature or features in this pull request is complete
- [ ] Suitable new documentation files and/or updates to the existing docs are included
- [ ] One or more example input decks are included
- [ ] The source code follows the SPARTA formatting guidelines

## Further Information, Files, and Links

_Put any additional information here, attach relevant text or image files, and URLs to external sites (e.g. DOIs or webpages)_
