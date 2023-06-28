#include "FileTools.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <fstream>
#include <iostream>

// ----------------------------------------------------------------------
/*!
 * \brief Read file containing a list of stations
 *
 * The file consists of rows of form "ID,lon,lat,name":
 *
 * Sample:
 *
 * 1001,25.0,65.0,Helsinki Kaisaniemi
 * 1002,25.0,65.1,Helsinki Vantaa
 *
 */
//------------------------------------------------------------------------

std::list<NFmiStation> read_filename_locations(const std::string& pFilename)
{
  std::list<NFmiStation> locations;

  std::ifstream in(pFilename);
  if (!in)
    throw std::runtime_error("Failed to open '" + pFilename + "' for reading");

  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty())
      continue;
    std::vector<std::string> parts;
    boost::algorithm::split(parts, line, boost::algorithm::is_any_of(","));
    if (parts.size() != 4)
      throw std::runtime_error("File '" + pFilename + "' line '" + line +
                               "' does not contain 4 comma separated fields");

    auto id = std::stoi(parts[0]);
    auto name = parts[3];
    auto lon = std::stod(parts[1]);
    auto lat = std::stod(parts[2]);

    locations.emplace_back(NFmiStation(id, name, lon, lat));
  }

  return locations;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read simulation points based on JSON filenames
 * Sample configuration:
 *
 * "points":
 * {
 *       "filename": "foo.txt";
 *       "filename": ["foo1.txt", "foo2.txt"];
 * }
 */
// ----------------------------------------------------------------------

std::list<NFmiStation> read_filename_locations(const Json::Value& pJson)
{
  Json::Value nulljson;
  // get points.keyword settings
  auto j = pJson.get("points", nulljson);
  if (!j.isNull())
    j = j.get("filename", nulljson);
  if (j.isNull())
    return {};

  std::list<NFmiStation> locations;
  if (j.isString())
    locations.splice(locations.end(), read_filename_locations(j.asString()));
  else if (j.isArray())
  {
    for (const auto& filename : j)
      locations.splice(locations.end(), read_filename_locations(filename.asString()));
  }
  else
    throw std::runtime_error("points.filename must be a string or an array of strings");

  return locations;
}
