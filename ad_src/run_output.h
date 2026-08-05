/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   run_output: shared experiment-output-directory naming, replacing two
   near-verbatim copies (opt_main.cpp, power_law_main.cpp). Keeps the
   stricter of the two originals: checks mkdir's return and errno rather
   than silently proceeding into failing ofstreams.
------------------------------------------------------------------------- */

#ifndef SPARTA_AD_SRC_RUN_OUTPUT_H
#define SPARTA_AD_SRC_RUN_OUTPUT_H

#include <string>

// Creates and returns "output/<prefix>_<YYYY-MM-DD>_<experiment>", or
// that path with a "_2", "_3", ... suffix appended if it already exists
// (same-day reruns of the same experiment name). Throws
// std::runtime_error if "output/" or the final directory can't be
// created.
std::string make_output_dir(const std::string &prefix,
                            const std::string &experiment);

#endif
