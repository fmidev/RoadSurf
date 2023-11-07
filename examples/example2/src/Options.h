// ----------------------------------------------------------------------
/*!
 * \brief Parsed command line options
 */
// ----------------------------------------------------------------------

#pragma once

#include <macgyver/DateTime.h>
#include <boost/optional.hpp>
#include <string>

struct Options
{
  bool verbose = false;
  int jobs = 0;
  std::string outfile;
  std::string infile;
  std::string configfile;
  std::string archive;
  boost::optional<Fmi::DateTime> time;
};
