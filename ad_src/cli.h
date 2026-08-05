/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   cli: shared CLI-token parsing. Replaces two disagreeing copies
   (shape_main.cpp's atof-based parser, which silently returned 0.0 on
   garbage; opt_main.cpp's istringstream-based parser, which silently
   stopped at the first bad token) with one strict parser that reports
   an error instead of guessing.
------------------------------------------------------------------------- */

#ifndef SPARTA_AD_SRC_CLI_H
#define SPARTA_AD_SRC_CLI_H

#include <string>
#include <vector>

std::string trim(const std::string &s);

// Comma- or whitespace-separated tokens. Throws std::invalid_argument
// (with the offending token and the whole input string) on anything
// that doesn't parse cleanly as a double/int.
std::vector<double> parse_doubles(const std::string &s);
std::vector<int> parse_ints(const std::string &s);

#endif
