/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   sparta_util: thin wrappers around SPARTA's public C library interface
   (src/library.h), shared by every ad_src driver/objective/runner.
   Replaces 4 near-identical private copies of cmd()/the open-argv-dance
   that previously lived in shape_case.cpp, power_law_main.cpp,
   drag_objective.cpp, and heat_flux_objective.cpp.
------------------------------------------------------------------------- */

#ifndef SPARTA_AD_SRC_UTIL_H
#define SPARTA_AD_SRC_UTIL_H

#include <string>

// Opens a no-MPI SPARTA instance. verbose=false runs with
// "-log none -screen none" (silent); verbose=true keeps SPARTA's normal
// stdout logging. Caller must close_sparta() the returned handle.
void *open_sparta(bool verbose);
void close_sparta(void *spa);

// Issues one SPARTA input-script command. Truncates at 1023 chars
// (matches sparta_command()'s own line-length assumption).
void cmd(void *spa, const char *str);

// sparta_extract_compute() takes a non-const char* for the compute
// name even though it never writes through it; this hides the
// const_cast every call site otherwise has to repeat.
void *extract_compute(void *spa, const char *name, int style, int type);

// Throws std::runtime_error(msg). Replaces the old per-file die()
// helpers, which each printed a message and called exit(1) directly --
// that skips destructors for any RAII objects on the stack (std::vector
// et al.), and hard-kills the process even from inside a context (an
// IPOPT line-search callback) where a caller further up might want to
// recover. Callers that can't safely let a C++ exception escape into
// non-exception-safe code (IPOPT's own call stack, notably --
// see shape_tnlp.cpp's eval_f/eval_grad_f) must catch it locally.
[[noreturn]] void die(const std::string &msg);

#endif
