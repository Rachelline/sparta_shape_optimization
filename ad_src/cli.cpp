/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   cli.cpp -- see cli.h
------------------------------------------------------------------------- */

#include "cli.h"

#include <cctype>
#include <sstream>
#include <stdexcept>

std::string trim(const std::string &s)
{
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char) s[a])) a++;
  while (b > a && std::isspace((unsigned char) s[b - 1])) b--;
  return s.substr(a, b - a);
}

namespace {

// Splits on commas (turned into spaces) and/or whitespace, then hands
// each token to `parse_one`. Throws std::invalid_argument, naming the
// bad token and the original string, on the first token that `parse_one`
// can't consume in full.
template <typename T, typename ParseOne>
std::vector<T> parse_tokens(const std::string &s, const char *what,
                            ParseOne parse_one)
{
  std::vector<T> out;
  std::string tmp = s;
  for (char &ch : tmp)
    if (ch == ',') ch = ' ';

  std::istringstream ss(tmp);
  std::string tok;
  while (ss >> tok) {
    size_t pos = 0;
    T val;
    try {
      val = parse_one(tok, pos);
    } catch (const std::exception &) {
      throw std::invalid_argument(std::string("not ") + what + ": '" + tok +
                                  "' in '" + s + "'");
    }
    if (pos != tok.size())
      throw std::invalid_argument(std::string("trailing garbage after ") +
                                  what + " '" + tok + "' (from '" + s + "')");
    out.push_back(val);
  }
  return out;
}

}  // namespace

std::vector<double> parse_doubles(const std::string &s)
{
  return parse_tokens<double>(s, "a number",
      [](const std::string &tok, size_t &pos) { return std::stod(tok, &pos); });
}

std::vector<int> parse_ints(const std::string &s)
{
  return parse_tokens<int>(s, "an integer",
      [](const std::string &tok, size_t &pos) { return std::stoi(tok, &pos); });
}
