// ----------------------------------------------------------------------
/*!
 * \brief Parsed command line options
 */
// ----------------------------------------------------------------------

#pragma once

#include <time.h>
#include <string>

struct Options
{
  bool verbose = false;
  int jobs = 0;
  std::string outfile;
  std::string infile;
  std::string configfile;
  time_t time;
};
