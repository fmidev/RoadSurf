// ----------------------------------------------------------------------
/*!
 * \brief Parsed command line options
 */
// ----------------------------------------------------------------------

#pragma once

#include <macgyver/DateTime.h>
#include <optional>
#include <string>

struct Options
{
  bool verbose = false;
  int jobs = 0;
  std::string outfile;
  std::string infile;
  std::string configfile;
  std::string archive;
  std::optional<Fmi::DateTime> time;
};
