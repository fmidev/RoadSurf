// ----------------------------------------------------------------------
/*!
 * \brief Parsed command line options
 */
// ----------------------------------------------------------------------

#pragma once

#include <boost/date_time/posix_time/ptime.hpp>
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
  boost::optional<boost::posix_time::ptime> time;
};
