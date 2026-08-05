/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   sparta_util.cpp -- see sparta_util.h
------------------------------------------------------------------------- */

#include "sparta_util.h"
#include "library.h"

#include <cstdio>
#include <stdexcept>

void *open_sparta(bool verbose)
{
  void *spa;
  char a0[] = "sparta", a1[] = "-log", a2[] = "none", a3[] = "-screen";
  char *argv_verbose[] = {a0, a1, a2};
  char *argv_quiet[]   = {a0, a1, a2, a3, a2};
  if (verbose) sparta_open_no_mpi(3, argv_verbose, &spa);
  else         sparta_open_no_mpi(5, argv_quiet, &spa);
  return spa;
}

void close_sparta(void *spa)
{
  sparta_close(spa);
}

void cmd(void *spa, const char *str)
{
  char buf[1024];
  std::snprintf(buf, sizeof(buf), "%s", str);
  sparta_command(spa, buf);
}

void *extract_compute(void *spa, const char *name, int style, int type)
{
  return sparta_extract_compute(spa, const_cast<char *>(name), style, type);
}

void die(const std::string &msg)
{
  throw std::runtime_error(msg);
}
