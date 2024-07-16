#include "RoadSurfSource.h"
#include "FileTools.h"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/timer/timer.hpp>
#include <cpr/cpr.h>
#include <json/json.h>
#include <json/reader.h>
#include <json/writer.h>
#include <smartmet/macgyver/NearTree.h>
#include <smartmet/macgyver/NearTreeLatLon.h>
#include <smartmet/macgyver/StringConversion.h>
#include <smartmet/macgyver/TimeParser.h>
#include <smartmet/newbase/NFmiPoint.h>
#include <smartmet/newbase/NFmiStringTools.h>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Observations = std::vector<Weather>;    // Observations sorted by time
using DataMap = std::map<int, Observations>;  // observations for each point

}  // namespace
// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

class RoadSurfSource::Impl
{
 public:
  Impl(const NFmiMetTime& pWallClock,
       const NFmiMetTime& pStartTime,
       const NFmiMetTime& pEndTime,
       const Json::Value& pSource,
       const Json::Value& pConfig);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon) const;

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const;

 private:
  const NFmiMetTime mWallClock;
  std::vector<std::string> mParams;
  NFmiMetTime mStartTime;
  NFmiMetTime mEndTime;
  int timeStep = -1;
  bool is_newTimeStep = true;
  NFmiMetTime lastTime = NFmiMetTime(2000, 1, 1, 0, 0, 0);  // initialize with some time
  NFmiMetTime time = NFmiMetTime(2000, 1, 1, 0, 0, 0);      // initialize with some time

  mutable DataMap mData;
  using Station = Fmi::NearTreeLatLon<int>;
  using StationTree = Fmi::NearTree<Station, Fmi::NearTreeLatLonDistance<Station>>;
  mutable StationTree mStationTree;
  std::vector<double> grid_coord;
  std::vector<double> grid_size;

  NFmiPoint mLatLon;
  double mRadius;
  void check_params();

  static double interpolate(const Observations& pObservations,
                            const NFmiMetTime& pTime,
                            int pPos,
                            int pMaxTimeGap,
                            double Weather::*pPtr);

  void read_file(const std::string& mFilename);
  void read_obs_list(const std::string& mFilename, const std::string& path);

};  // class AsciiSource::Impl

// ----------------------------------------------------------------------
/*!
 * \brief Read data from a file
 */
// ----------------------------------------------------------------------

void RoadSurfSource::Impl::read_file(const std::string& mFilename)
{
  // Read input file
  std::ifstream in(mFilename.c_str());
  if (!in)
    throw std::runtime_error("Failed to open '" + mFilename + "' for reading");

  std::string line;
  // Test by reading first line
  if (!std::getline(in, line))
    throw std::runtime_error("Failed to read line 1 from '" + mFilename + "'");

  using boost::algorithm::contains;

  // Data rows: yy mm dd hh tair rh vz rr1h rform srad lrad tsurf

  int yy = 0;
  int mm = 0;
  int dd = 0;
  int hh = 0;
  int mi = 0;
  int ss = 0;
  int pnum = 0;
  int dataRow = 0;
  int dataColumn = 0;
  int rownum = 1;
  int pointNum = 0;
  auto last_pos = mData.end();  // Pointer to Observation class in datamap

  // Go trough rows until the end of file
  while (in.good())
  {
    // read line to "line" variable
    std::getline(in, line);
    // if timestamp row
    if (line.rfind("AIKA:", 0) == 0)
    {
      // Parse timestamp
      std::string timestamp = line.substr(5, 12);
      yy = std::stoi(timestamp.substr(0, 4));
      mm = std::stoi(timestamp.substr(4, 2));
      dd = std::stoi(timestamp.substr(6, 2));
      hh = std::stoi(timestamp.substr(8, 2));
      mi = std::stoi(timestamp.substr(10, 2));
      ss = 0;          // seconds
      dataRow = 0;     // used to calculate lat
      dataColumn = 0;  // used to calculate lon
      pointNum = 0;
      // create time variable
      time = NFmiMetTime(yy, mm, dd, hh, mi, ss);
      // if different than last timesamp,
      if (lastTime != time)
      {
        is_newTimeStep = true;
        timeStep++;
      }
      else
      {
        is_newTimeStep = false;
      }
      lastTime = time;
    }
    // Get parameter number
    else if (line.rfind("PARAMETRI:", 0) == 0)
    {
      // string stream goes through words in line
      std::string temp;
      std::stringstream words;
      words << line;
      while (!words.eof())
      {
        words >> temp;
        // Check if integer
        if (std::stringstream(temp) >> pnum)
          break;
      }
    }
    else
    {
      // first word in line
      std::string number_string = line.substr(0, line.find(' '));
      bool is_valid;

      // check if first parameter in line is float
      try
      {
        static_cast<void>(std::stod(number_string));
        is_valid = true;
      }
      catch (std::exception& ia)
      {
        is_valid = false;
      }
      // read data if float was found
      if (is_valid)
      {
        // isstringstream goes trough line word by word
        std::istringstream iss(line);
        std::string value_string;
        while (iss >> value_string)
        {
          // convert string to double
          auto value = std::stod(value_string);

          int id = pointNum;
          // make map with location id as key and Observation class as value
          // get iterator to the element
          last_pos = mData.find(id);
          // if key was not in the map yet, insert it there
          if (last_pos == mData.end())
          {
            auto pos_bool = mData.insert(std::make_pair(id, Observations()));
            last_pos = pos_bool.first;
            if (grid_coord.empty())
              throw std::runtime_error("Grid not set");

            double lat = grid_coord[1] + dataRow * (grid_coord[3] - grid_coord[1]) / grid_size[1];
            double lon =
                grid_coord[0] + dataColumn * (grid_coord[2] - grid_coord[0]) / grid_size[0];
            mStationTree.insert(Station{lon, lat, id});
          }
          // if new timestep, generate Weather class where to insert data
          if (is_newTimeStep)
          {
            Weather data{time};
            last_pos->second.push_back(data);
          }
          // if value is not missing, store it to the corresponding paramter slot
          if (value > -9999)
          {
            // get the observation class object for current time and location
            auto& obs = last_pos->second;

            if (pnum == 21)
              obs[timeStep].vz = value;
            if (pnum == 4)
              obs[timeStep].t2m = value;
            if (pnum == 13)
              obs[timeStep].rh2m = value;
            if (pnum == 317)
              obs[timeStep].swdn = value;
            if (pnum == 315)
              obs[timeStep].lwdn = value;
            if (pnum == 353)
              obs[timeStep].simuprec = value;
            if (pnum == 354)
              obs[timeStep].simuprec = value / 3;  // 3 hour precipitation
            if (pnum == 5)
              obs[timeStep].troad = value;
            if (pnum == 57)
              obs[timeStep].phase = value;
          }
          dataColumn++;
          pointNum++;
        }
        dataRow++;
        dataColumn = 0;
      }
      rownum++;
    }
    if (in.eof())
      break;
  }
  if (in.bad())
    throw std::runtime_error("Failed to read row " + std::to_string(rownum) + " from '" +
                             mFilename + "'");
}
// ----------------------------------------------------------------------
/*!
 * \brief Read data from a files given in a list
 */
// ----------------------------------------------------------------------

void RoadSurfSource::Impl::read_obs_list(const std::string& mFilename, const std::string& path)
{
  // Read file list
  std::cout << "Reading file " << mFilename << std::endl;
  std::ifstream in(mFilename.c_str());
  if (!in)
    throw std::runtime_error("Failed to open '" + mFilename + "' for reading");

  // To sort files in cronological order, create vector containing pairs,
  // where first is timestamp and the second filename.
  std::string fileName;
  std::vector<std::pair<std::string, std::string>> fileNames;
  int rownum = 1;
  // Go trough rows in file name list until the end of file
  while (in.good())
  {
    // read line to "line" variable
    std::getline(in, fileName);
    if (!fileName.empty() and !std::all_of(fileName.begin(), fileName.end(), isspace))
    {  // full path to data file
      std::string full_name = path + fileName;
      // open datafile to get timestamp
      std::ifstream intime(full_name.c_str());
      if (!intime)
        throw std::runtime_error("Failed to open '" + full_name + "' for reading");

      std::string line;
      if (!std::getline(intime, line))
        throw std::runtime_error("Failed to read line 1 from '" + full_name + "'");

      // read line to "line" variable
      std::getline(intime, line);
      std::string ending;
      // if timestamp row
      if (line.rfind("AIKA:", 0) == 0)
      {
        // Parse timestamp
        std::string time_stamp = line.substr(5, 12);
        // as filenames will be organized alphabetically, put b at the end of radar observation
        // so that they will override "ilmasto_hila" when available
        if (fileName.rfind("Tutkasade") == 0)
        {
          ending = "b";
        }
        else
        {
          ending = "a";
        }

        fileNames.emplace_back(make_pair(time_stamp + ending, fileName));
      }
      intime.close();
    }
    rownum++;
    if (in.eof())
      break;
  }
  if (in.bad())
    throw std::runtime_error("Failed to read row " + std::to_string(rownum) + " from '" + path +
                             mFilename + "'");
  // sort vector containing pairs of timetamps an filenames
  std::sort(fileNames.begin(), fileNames.end());

  // read each file
  for (const auto& name : fileNames)
    read_file(path + name.second);
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct source from ASCII file
 *
 * Sample definition:
 *
 *        {
 *            "name":        "ascii",
 *            "filename":    "RoadRunner.txt"
 *        }
 */
// ----------------------------------------------------------------------

RoadSurfSource::Impl::Impl(const NFmiMetTime& pWallClock,
                           const NFmiMetTime& pStartTime,
                           const NFmiMetTime& pEndTime,
                           const Json::Value& pSource,
                           const Json::Value& pConfig)
    : mWallClock(pWallClock), mStartTime(pStartTime), mEndTime(pEndTime)
{
#if 0
  Json::StyledWriter writer;
  std::cout << "JSON:\n" << writer.write(theSource) << std::endl;
#endif

  Json::Value nulljson;

  // Extract name for error messages only. The name is known to exist
  // since DataSourceFactory depends on it.

  auto name = pSource.get("name", nulljson).asString();
  mRadius = 2.0;
  // Parse resource location

  auto jfile = pSource.get("filename", nulljson);
  if (jfile.isNull())
    throw std::runtime_error("filename setting missing for ascii source");

  // Parameters

  auto json = pSource.get("params", nulljson);
  if (json.isNull())
    std::cerr << "Warning: Data source '" + name + "' has no active parameters" << std::endl;
  else
  {
    if (!json.isArray())
      throw std::runtime_error("Parameter setting for source '" + name + "' must be an array");
    for (const auto& p : json)  // NOLINT(cppcheck-useStlAlgorithm)
      mParams.push_back(p.asString());
    check_params();
  }

  // Get projection It should be given in json as "projection":"latlon:
  // west_longitude,south_latitude,east_longitude,north_latitude"
  auto jproj = Json::Path(".points.projection").resolve(pConfig);
  if (!jproj.isNull())
  {
    // read projection setting as string
    const std::string projection = jproj.asString();
    const char* separator = ":";
    // make vector by separating string by separator
    const auto parts = NFmiStringTools::Split<std::vector<std::string>>(projection, separator);
    // make vectors for coordinate points and number of grid points
    const auto aparts1 = NFmiStringTools::Split<std::vector<std::string>>(parts[1], "/");
    const auto aparts2 = NFmiStringTools::Split<std::vector<std::string>>(parts[2], "/");
    grid_coord = NFmiStringTools::Split<std::vector<double>>(aparts1[0]);
    grid_size = NFmiStringTools::Split<std::vector<double>>(aparts2[0]);
  }

  // Read data
  std::string mFilename = jfile.asString();
  std::string list = "ascii_obs_list";
  // check if list of files
  if (name == list)
  {
    // Read observation data from obs_list
    std::string path = pSource.get("path", nulljson).asString();
    if (path.empty())
      throw std::runtime_error("file path missing from ascii source");
    read_obs_list(mFilename, path);
  }
  else
    read_file(mFilename);
}

// ----------------------------------------------------------------------
/*!
 * \brief Check that the required querydata parameters exist
 */
// ----------------------------------------------------------------------

void RoadSurfSource::Impl::check_params()
{
  static std::set<std::string> known_params = {{"roadtemperature"},
                                               {"airtemperature"},
                                               {"dewpoint"},
                                               {"humidity"},
                                               {"windspeed"},
                                               {"longwaveradiation"},
                                               {"shortwaveradiation"},
                                               {"precipitation"},
                                               {"precipitationrate"},
                                               {"precipitationform"}};

  // {"precipitationform"}};

  for (const auto& param : mParams)
  {
    auto pos = known_params.find(param);
    if (pos == known_params.end())
      throw std::runtime_error("Unsupported SmartMet parameter '" + param + "'");
  }
}
// ----------------------------------------------------------------------
/*!
 * \brief Intepolate a single value
 *
 * The given index is such that theObservations[thePos].date >= theTime
 */
// ----------------------------------------------------------------------

double RoadSurfSource::Impl::interpolate(const Observations& pObservations,
                                         const NFmiMetTime& pTime,
                                         int pPos,
                                         int pMaxTimeGap,
                                         double Weather::*pPtr)
{
  if (pObservations[pPos].date == pTime)
  {
    double value = pObservations[pPos].*pPtr;
    if (!is_missing(value))
      return value;
  }

  if (pPos == 0)
    return MISSING;

  // Search for first value after current time with valid value

  unsigned int pos2;
  double value2 = MISSING;
  for (pos2 = pPos; pos2 < pObservations.size(); ++pos2)
  {
    value2 = pObservations[pos2].*pPtr;
    if (!is_missing(value2))
      break;
  }

  if (is_missing(value2))
    return value2;

  // Search for first value before current time with valid value

  unsigned int pos1;
  double value1 = MISSING;
  for (pos1 = pPos - 1;; --pos1)
  {
    value1 = pObservations[pos1].*pPtr;
    if (pos1 == 0 || !is_missing(value1))
      break;
  }

  if (is_missing(value1))
    return value1;

  // Do time interpolation if the gap is not too long

  const auto& t1 = pObservations[pos1].date;
  const auto& t2 = pObservations[pos2].date;

  const int gap = t2.DifferenceInMinutes(t1);
  if (gap > pMaxTimeGap)
    return MISSING;

  // Do time interpolation

  const int gap1 = pTime.DifferenceInMinutes(t1);

  double value = ((gap - gap1) * value1 + gap1 * value2) / gap;
  return value;
}

// ----------------------------------------------------------------------
/*!
 * \brief Get the weather for the given point
 */
// ----------------------------------------------------------------------

void RoadSurfSource::Impl::GetWeather(InputData& pData,
                                      const SimulationTimes& pTimes,
                                      const NFmiPoint& pLatLon) const
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return;

  // No station near enough?
  auto nearest_station =
      mStationTree.nearest(Station{pLatLon.X(), pLatLon.Y()}, Station::ChordLength(mRadius));
  if (!nearest_station)
    return;
  // No data for the station?!?!
  auto myData = mData.find(nearest_station->ID());
  if (myData == mData.end())
    return;
  const auto& obs = myData->second;
  if (obs.empty())
    return;
  for (unsigned int pos = 0; pos < obs.size(); ++pos)
  {
  }
  // Process all times
  for (std::size_t i = 0; i < pTimes.size(); i++)
  {
    const auto& pt = pTimes[i].pt;
    if (pt < mStartTime || pt > mEndTime)
      continue;
    // Not in the time interval?
    if (pt < obs.front().date || pt > obs.back().date)
      continue;

    // Find first position where t >= time. This is guaranteed to succeed
    // after the above tests.

    unsigned int pos;
    for (pos = 0; pos < obs.size(); ++pos)
      if (obs[pos].date.PosixTime() >= pt)
        break;
    const int max_time_gap = 180;

    double Weather::*ptr = nullptr;  // pointer to double member of struct

    Weather data{pt};

    for (const auto& p : mParams)
    {
      if (p == "roadtemperature")
        ptr = &Weather::troad;
      else if (p == "airtemperature")
        ptr = &Weather::t2m;
      else if (p == "dewpoint")
        ptr = &Weather::tdew2m;
      else if (p == "humidity")
        ptr = &Weather::rh2m;
      else if (p == "windspeed")
        ptr = &Weather::vz;
      else if (p == "longwaveradiation")
        ptr = &Weather::lwdn;
      else if (p == "shortwaveradiation")
        ptr = &Weather::swdn;
      else if (p == "precipitation" || p == "precipitationrate")
        ptr = &Weather::simuprec;
      else if (p == "precipitationform")
        ptr = &Weather::phase;

      if (ptr)
        data.*ptr = interpolate(obs, pt, pos, max_time_gap, ptr);
    }

    // clamp rh to 0-100
    if (!is_missing(data.rh2m))
      data.rh2m = std::max(0.0, std::min(100.0, data.rh2m));

    // sanity check on precipitation amount
    if (data.simuprec > 100)
      data.simuprec = MISSING;

    // TODO: tdew????

    if (!is_missing(data.troad))
      pData.TSurfObs[i] = data.troad;
    if (!is_missing(data.t2m))
      pData.tair[i] = data.t2m;
    if (!is_missing(data.rh2m))
      pData.Rhz[i] = data.rh2m;
    if (!is_missing(data.vz))
      pData.VZ[i] = data.vz;
    if (!is_missing(data.lwdn))
      pData.LW[i] = data.lwdn;
    if (!is_missing(data.swdn))
      pData.SW[i] = data.swdn;
    if (!is_missing(data.simuprec))
      pData.prec[i] = data.simuprec;
    if (!is_missing(data.phase))
      pData.PrecPhase[i] = data.phase;
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish latest available road observation time for setting coupling end time
 */
// ----------------------------------------------------------------------

std::optional<NFmiMetTime> RoadSurfSource::Impl::GetLatestObsTime(
    const NFmiPoint& pLatLon, const std::string& variable) const
{
  // Not in this source if parameter is not defined in the configuration
  auto pos = find(mParams.begin(), mParams.end(), variable);
  if (pos == mParams.end())
    return {};

  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return {};

  // No station near enough?
  auto nearest_station =
      mStationTree.nearest(Station{pLatLon.X(), pLatLon.Y()}, Station::ChordLength(mRadius));
  if (!nearest_station)
    return {};

  // No data for the station?!?!
  auto station_pos = mData.find(nearest_station->ID());
  if (station_pos == mData.end())
    return {};

  const auto& obs = station_pos->second;

  for (std::size_t i = obs.size(); i > 0; i--)
  {
    if (!is_missing(obs[i - 1].troad))
      return {obs[i - 1].date};
  }
  return {};
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct a SmartMet source for weather
 */
// ----------------------------------------------------------------------

RoadSurfSource::RoadSurfSource(const NFmiMetTime& pWallClock,
                               const NFmiMetTime& pStartTime,
                               const NFmiMetTime& pEndTime,
                               const Json::Value& pSource,
                               const Json::Value& pConfig)
    : impl(new RoadSurfSource::Impl(pWallClock, pStartTime, pEndTime, pSource, pConfig))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Get weather for the given point and time
 */
// ----------------------------------------------------------------------

void RoadSurfSource::GetWeather(InputData& pData,
                                const SimulationTimes& pTimes,
                                const NFmiPoint& pLatLon) const
{
  return impl->GetWeather(pData, pTimes, pLatLon);
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish latest available road observation time for setting coupling end time
 */
// ----------------------------------------------------------------------

std::optional<NFmiMetTime> RoadSurfSource::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                              const std::string& variable) const
{
  return impl->GetLatestObsTime(pLatLon, variable);
}
