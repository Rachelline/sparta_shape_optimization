/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   run_output.cpp -- see run_output.h
------------------------------------------------------------------------- */

#include "run_output.h"

#include <cerrno>
#include <cstring>
#include <ctime>
#include <stdexcept>
#include <sys/stat.h>

namespace {

bool dir_exists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

}  // namespace

std::string make_output_dir(const std::string &prefix,
                            const std::string &experiment)
{
  if (mkdir("output", 0755) != 0 && errno != EEXIST)
    throw std::runtime_error(std::string("cannot create 'output': ") +
                             std::strerror(errno));

  char datestr[32];
  time_t t = time(nullptr);
  strftime(datestr, sizeof(datestr), "%Y-%m-%d", localtime(&t));

  std::string base = "output/" + prefix + "_" + datestr + "_" + experiment;
  std::string dir = base;
  for (int suffix = 2; dir_exists(dir); suffix++)
    dir = base + "_" + std::to_string(suffix);

  if (mkdir(dir.c_str(), 0755) != 0 && errno != EEXIST)
    throw std::runtime_error("cannot create output dir '" + dir + "': " +
                             std::strerror(errno));

  return dir;
}
