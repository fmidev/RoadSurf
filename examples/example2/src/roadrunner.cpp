#include "DataManager.h"
#include "InputData.h"
#include "InputSettings.h"
#include "InputParameters.h"
#include "JsonTools.h"
#include "LocalParameters.h"
#include "Mutex.h"
#include "OutputData.h"
#include "PointMode.h"
#include "QueryDataSymbols.h"
#include "QueryDataTools.h"
#include "SimulationTime.h"
#include "SkyView.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/asio.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/vector.hpp>
#include <json/json.h>
#include <smartmet/macgyver/StringConversion.h>
#include <smartmet/macgyver/TimeParser.h>
#include <smartmet/newbase/NFmiFastQueryInfo.h>
#include <smartmet/newbase/NFmiMetTime.h>
#include <smartmet/newbase/NFmiQueryData.h>
#include <atomic>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

// The fortran API
extern "C"
{
  void runsimulation(OutputPointers* pOutputPointers,
                     const InputPointers* pInputPointers,
                     const InputSettings* pSettings,
                     const InputParameters* pInputParams,
                     const LocalParameters* lParameters);
 
}
namespace
{
MutexType read_mutex;
}
using namespace boost::archive;
using namespace boost::posix_time;
using namespace boost::filesystem;
// Parsed command line options
Options options;

// ----------------------------------------------------------------------
/*!
 * \brief Parse command line options
 * \return True on success, false if execution is to be stopped
 */
// ----------------------------------------------------------------------

bool parse_options(int argc, char* argv[])
{
  namespace po = boost::program_options;

  std::string stime;

  po::options_description desc("Available options");
  // clang-format off
  desc.add_options()("help,h", "print help message")
      ("verbose,v", po::bool_switch(&options.verbose), "verbose mode")
      ("jobs,j", po::value(&options.jobs)->implicit_value(-1), "number of parallel jobs")
      ("outfile,o", po::value(&options.outfile), "output filename")
      ("infile,i", po::value(&options.infile), "additional comma separated input datafiles (mode=ascii) to be used last")
      ("time,t", po::value(&stime), "simulation time")
      ("config,c", po::value(&options.configfile), "configuration file");
  // clang-format on

  po::positional_options_description p;
  p.add("config", 1);

  po::variables_map opt;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), opt);

  po::notify(opt);

  if (opt.count("help") != 0)
  {
    std::cout << "Usage: roadrunner [options] [configfile]\n\n"
                 "Run road model using the given configuration.\n\n"
              << desc << std::endl;
    return false;
  }

  if (options.configfile.empty())
    throw std::runtime_error("Configuration file not given");

  if (!boost::filesystem::exists(options.configfile))
    throw std::runtime_error("Configuration file '" + options.configfile + "' missing");

  // Plain -j means use all cores
  if (options.jobs < 0)
    options.jobs = std::thread::hardware_concurrency();

  if (opt.count("time") != 0)
    options.time = Fmi::TimeParser::parse(stime);

  return true;
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct the simulation timesteps once for speed
 */
// ----------------------------------------------------------------------

SimulationTimes get_times(const InputSettings& pSettings)
{
  SimulationTimes times;

  const auto dt = boost::posix_time::seconds(static_cast<long>(pSettings.DTSecs));
  auto t = pSettings.start_time;
  for (int i = 0; i < pSettings.SimLen; i++)
  {
    times.emplace_back(SimulationTime{t});
    t += dt;
  }

  return times;
}


// ----------------------------------------------------------------------
/*!
 * \brief Read data required for simulation
 */
// ----------------------------------------------------------------------

boost::optional<InputData> read_input(const NFmiPoint& pLonLat,
                                      const SimulationTimes& pTimes,
                                      InputSettings& pSettings,
                                      const DataManager& pDataManager,
                                      const InputParameters& parameters,
                                      LocalParameters& lParameters,
                                      SkyView& skyview)
{
  InputData data{pSettings.SimLen};
  // Establish simulation times

  for (int i = 0; i < pSettings.SimLen; i++)
  {
    data.year[i] = pTimes[i].year;
    data.month[i] = pTimes[i].month;
    data.day[i] = pTimes[i].day;
    data.hour[i] = pTimes[i].hour;
    data.minute[i] = pTimes[i].minute;
    data.second[i] = pTimes[i].second;
  }
  // Set Initialization length used by relaxation
  const auto init_duration = pSettings.forecast_time - pSettings.start_time;
  const auto init_secs = init_duration.total_seconds();
  lParameters.InitLenI = 1 + static_cast<int>(init_secs / pSettings.DTSecs);

  std::string time_string = Fmi::to_iso_string(pSettings.forecast_time);
  std::string coord_string = std::to_string(pLonLat.X()) + "_" + std::to_string(pLonLat.Y());
  pDataManager.GetWeather(data, pTimes, pLonLat);
  skyview.GetSkyVariables(lParameters, data, pLonLat);
  // Verify the model has all required parameters
  lParameters.lat = pLonLat.Y();
  lParameters.lon = pLonLat.X();
  for (int i = 0; i < pSettings.SimLen; i++)
  {
    const auto& pt = pTimes[i].pt;

    if (is_missing(data.tair[i]))
    {
      std::cout << "2m temperature missing at " + to_iso_extended_string(pt) << " " << pLonLat.X()
                << " " << pLonLat.Y() << std::endl;
      return {};
    }
    if (is_missing(data.Rhz[i]))
    {
      std::cout << "RH missing at " + to_iso_extended_string(pt) << " " << pLonLat.X() << " "
                << pLonLat.Y() << std::endl;
      return {};
    }
    if (is_missing(data.prec[i]))
    {
      std::cout << "Precipitation missing at " + to_iso_extended_string(pt) << " " << pLonLat.X()
                << " " << pLonLat.Y() << std::endl;
      return {};
    }
    if (is_missing(data.SW[i]))
    {
      std::cout << "SW missing at " + to_iso_extended_string(pt) << " " << pLonLat.X() << " "
                << pLonLat.Y() << std::endl;
      return {};
    }
    if (is_missing(data.LW[i]))
    {
      std::cout << "LW missing at " + to_iso_extended_string(pt) << " " << pLonLat.X() << " "
                << pLonLat.Y() << std::endl;
      return {};
    }
    if (is_missing(data.VZ[i]))
    {
      std::cout << "VZ missing at  " + to_iso_extended_string(pt) << " " << pLonLat.X() << " "
                << pLonLat.Y() << std::endl;
      return {};
    }

#if 0
    if (!is_missing(weather.troad))
      data.TSurfObs[i] = weather.troad;
#endif
    }
    // If relaxation is used, set first forecasted values
    //(Should be forecasted values at the time of the latest observation, but this is
    // good enough with small timesteps)
    if (pSettings.use_relaxation == 1)
    {
      lParameters.tair_relax = -9999.9;
      lParameters.VZ_relax = -9999.9;
      lParameters.RH_relax = -9999.9;
      auto lastTairObsTime = pDataManager.GetLatestObsTime(pLonLat, "airtemperature");
      if (lastTairObsTime)
      {
        boost::posix_time::ptime obstime = *lastTairObsTime;
        const auto duration = obstime - pSettings.start_time;
        const auto total_secs = duration.total_seconds();
        lParameters.InitLenI = static_cast<int>(total_secs / pSettings.DTSecs) + 1;
        lParameters.tair_relax = data.tair[lParameters.InitLenI];
        lParameters.VZ_relax = data.VZ[lParameters.InitLenI];
        lParameters.RH_relax = data.Rhz[lParameters.InitLenI];
      }
    }
    // if coupling is used, determine coupling index and latest surface temperature observations
    if (pSettings.use_coupling == 1)
    {
      lParameters.couplingTsurf = -9999.9;
      lParameters.couplingIndexI = -9999;
      int i = data.TSurfObs.size() - 1;
      // Find the latest surface temeprature observation
      while (i >= 0 and (is_missing(data.TSurfObs[i]) || data.TSurfObs[i] < -100))
      {
        i = i - 1;
      }
      // If i is at least as large as the copuling length
      if (i >= static_cast<int>(pSettings.coupling_minutes * 60 / pSettings.DTSecs))
      {
        lParameters.couplingTsurf = data.TSurfObs[i];
        lParameters.couplingIndexI = i;
        int couplingStartI =
            i - static_cast<int>(pSettings.coupling_minutes * 60 / pSettings.DTSecs);
        // During coupling there should not be surface temperature values as input
        for (int j = i; j > couplingStartI; j--)
        {
          data.TSurfObs[j] = -9999.9;
        }
      }
    }

  return data;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read mask from querydata file
 */
// ----------------------------------------------------------------------

std::vector<bool> read_querydata_mask(const Json::Value& mask, NFmiFastQueryInfo* qi)
{
  Json::Value nulljson;

  auto path = mask.get("path", nulljson).asString();
  if (path.empty())
    throw std::runtime_error("file path missing from mask");

  auto formula = mask.get("enable", nulljson).asString();
  if (formula.empty())
    throw std::runtime_error("Formula for querydata mask missing");

  // Read mask data
  NFmiQueryData mask_qd(path);
  NFmiFastQueryInfo mask_qi(&mask_qd);

  // Set up the expression evaluator

  stx::ParseTree expr = stx::parseExpression(formula);

  QueryDataSymbols symbols;
  symbols.setData(mask_qi);

  int enabled = 0;
  int disabled = 0;

  auto nlocations = qi->SizeLocations();
  std::vector<bool> point_mask(nlocations, false);

  std::size_t pos = 0;
  for (qi->ResetLocation(); qi->NextLocation(); ++pos)
  {
    symbols.setLocation(qi->LatLon());
    stx::AnyScalar val = expr.evaluate(symbols);

    if (!val.isBooleanType())
      throw std::runtime_error("Expression " + formula + " value must be boolean");

    bool flag = val.getBoolean();
    point_mask[pos] = flag;
    if (flag)
      ++enabled;
    else
      ++disabled;
  }

  std::cout << "Using querydata mask " << path << std::endl
            << "\tenabled  " << enabled << " points" << std::endl
            << "\tdisabled " << disabled << " points" << std::endl;

  return point_mask;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read mask from ASCII file
 */
// ----------------------------------------------------------------------

std::vector<bool> read_ascii_mask(const Json::Value& mask, NFmiFastQueryInfo* qi)
{
  Json::Value nulljson;

  auto path = mask.get("path", nulljson).asString();
  if (path.empty())
    throw std::runtime_error("file path missing from mask");

  // Read the mask data
  std::cout << "Reading ASCII mask " << path << std::endl;
  std::ifstream in(path.c_str());
  if (!in)
    throw std::runtime_error("Failed to read mask");

  // create vector holding vectors of mask values (inner vectors represent west-east grid rows)
  // false means that the grid point is skipped

  auto nlocations = qi->SizeLocations();
  std::vector<bool> point_mask(nlocations, false);

  std::vector<std::vector<bool> > mask_values;

  std::string line;
  char no_data = '-';

  // read mask into vector of bool vectors

  while (in.good())
  {
    std::getline(in, line);
    std::vector<bool> row;
    if (!line.empty())
    {
      // Go trough characters is row
      for (char mask_value : line)
      {
        if (mask_value == no_data)
          row.push_back(false);
        else
          row.push_back(true);
      }
      mask_values.push_back(row);
    }
    if (in.eof())
      break;
  }
  if (in.bad())
    throw std::runtime_error("Failed to read row from " + path);

  // Read some grid parameters from NFmiFastQueryInfo

  const NFmiArea* area = qi->Area();
  const NFmiGrid* grid = qi->Grid();
  auto top = area->TopLeftLatLon().Y();
  auto left = area->TopLeftLatLon().X();
  auto right = area->BottomRightLatLon().X();
  auto bottom = area->BottomRightLatLon().Y();
  auto xnum = grid->XNumber();
  auto ynum = grid->YNumber();
  auto xwidth = (right - left) / (xnum - 1);
  auto ywidth = (top - bottom) / (ynum - 1);

  // Map from 2D grid into 1D vector of locations
  for (std::size_t loc_index = 0; loc_index < nlocations; ++loc_index)
  {
    qi->LocationIndex(loc_index);
    auto lonlat = qi->LatLon();
    double x = lonlat.X();
    double y = lonlat.Y();
    // determine location of the mask value for current point based on grid values
    long cnum = std::round((x - left) / xwidth);
    long rnum = std::round((y - bottom) / ywidth);
    // add mask value to location
    point_mask[loc_index] = mask_values[rnum][cnum];
  }

  return point_mask;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read mask to skip certain points
 * \return Vector of mask values in same order as locations are in NFmiFastQueryInfo.
 */
// ----------------------------------------------------------------------

std::vector<bool> read_mask(const Json::Value& mask, NFmiFastQueryInfo* qi)
{
  auto nlocations = qi->SizeLocations();

  // enable all points by default

  if (mask.isNull())
    return std::vector<bool>(nlocations, true);

  // get mask type

  Json::Value nulljson;
  std::string type = mask.get("type", nulljson).asString();

  if (type == "ascii")
    return read_ascii_mask(mask, qi);

  if (type == "querydata")
    return read_querydata_mask(mask, qi);

  throw std::runtime_error("Unknown mask type '" + type + "'");
}

// ----------------------------------------------------------------------
/*!
 * \brief Write output to file for single coordinate
 */
// ----------------------------------------------------------------------

void put_file_contents(const std::string& outfile, const OutputData& output)
{
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("Failed to open " + outfile + " for writing");
  out << output;
  out.close();
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for a single coordinate
 */
// ----------------------------------------------------------------------

void run_coordinate(const Json::Value& pJson,
                    InputSettings& pSettings,
                    const DataManager& pDataManager)
{
  // Read data for a single point

  auto jlon = Json::Path(".points.longitude").resolve(pJson);
  auto jlat = Json::Path(".points.latitude").resolve(pJson);

  if (jlon.isNull() || jlat.isNull())
    throw std::logic_error("points.longitude and/or points.latitude not set");

  // run road weather model

  Json::Value nulljson;
  InputParameters parameters{pSettings, pJson.get("parameters", nulljson)};
  LocalParameters lParameters(pSettings);

  SkyView skyview(pJson.get("parameters", nulljson));
  NFmiPoint lonlat(jlon.asDouble(), jlat.asDouble());
  const auto times = get_times(pSettings);
  // Get input data for point. Returns nothing if missing value is found in the data
  const auto input =
      read_input(lonlat, times, pSettings, pDataManager, parameters, lParameters, skyview);

  if (input)
  {
    // initialize output vectors with missing values
    OutputData output(pSettings.SimLen);
    const InputData& inputT = *input;
    // Pass the data to Fortran via pointers
    auto outPointers = output.pointers();
    auto inPointers = inputT.pointers();
    runsimulation(
        &outPointers, &inPointers, &pSettings,
        &parameters,&lParameters);
    if (options.outfile.empty() || options.outfile == "-")
    {
      std::cout << output << std::endl;
    }
    else
    {
      put_file_contents(options.outfile, output);
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for querydata locations
 */
// ----------------------------------------------------------------------

void run_querydata_locations_sync(NFmiQueryData& pQ,
                                  const Json::Value& pJson,
                                  InputSettings& pSettings,
                                  const DataManager& pDataManager)
{
  // Early error checks
  auto jfilename = Json::Path(".output.filename").resolve(pJson);
  if (jfilename.isNull() && options.outfile.empty())
    throw std::runtime_error("Output filename not set");

  Json::Value nulljson;
  InputParameters parameters{pSettings, pJson.get("parameters", nulljson)};
  SkyView skyview(pJson.get("parameters", nulljson));
  const auto times = get_times(pSettings);
  NFmiFastQueryInfo qi(&pQ);
  long nlocations = qi.SizeLocations();
  auto stride = get_write_stride(qi, pSettings);

  // Read mask if set in json file
  auto mask = Json::Path(".points.mask").resolve(pJson);
  const auto point_mask = read_mask(mask, &qi);

  // set ratio of allowed missing data (e.g 0.1 means 10 %)
  auto jallowed_missing_ratio = Json::Path(".points.allow_missing_ratio").resolve(pJson);

  double allowed_missing_ratio;
  if (jallowed_missing_ratio.isNull())
    allowed_missing_ratio = 0.1;
  else
    allowed_missing_ratio = jallowed_missing_ratio.asDouble();

  // for tracking how many points are skipped due to missing data
  int missing_data_points = 0;

  for (long loc_index = 0; loc_index < nlocations; ++loc_index)
  {
    if (loc_index % 1000 == 0)
     std::cout << "\t" << loc_index << " / " << nlocations << " points" << std::endl;
    if (point_mask.at(loc_index))
    {
      qi.LocationIndex(loc_index);
      LocalParameters lParameters(pSettings);
      auto lonlat = qi.LatLon();
      // Get input data for point. Returns nothing if missing value is found in the data
      const auto input =
          read_input(lonlat, times, pSettings, pDataManager, parameters, lParameters, skyview); 
     if (input)
      {
        // initialize output vectors with missing values
        OutputData output(pSettings.SimLen);
        const InputData& inputT = *input;
        // Pass the data to Fortran via pointers
        auto outPointers = output.pointers();
        auto inPointers = inputT.pointers();
        runsimulation(
            &outPointers, &inPointers, &pSettings,
            &parameters, &lParameters);
 //       std::cout<<"here"<<std::endl;
  //      runsimulationExample(&outPointers, &inPointers, &pSettings,
   //                       &parameters,&lParameters);
        write_timeseries(qi, loc_index, inputT, output, stride);
      }
      else
      {
        missing_data_points++;
        // If there are too many missing data points, throw error
        if ((double)missing_data_points / nlocations > allowed_missing_ratio)
          throw std::runtime_error("Too much missing data");
      }
    }
  }
  std::string filename = (options.outfile.empty() ? jfilename.asString() : options.outfile);
  pQ.Write(filename);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for querydata locations
 */
// ----------------------------------------------------------------------

void run_querydata_locations_async(NFmiQueryData& pQ,
                                   const Json::Value& pJson,
                                   InputSettings& pSettings,
                                   const DataManager& pDataManager)
{
  // Early error checks

  auto jfilename = Json::Path(".output.filename").resolve(pJson);
  if (jfilename.isNull())
    throw std::runtime_error("Output filename not set");

  Json::Value nulljson;
  InputParameters parameters{pSettings, pJson.get("parameters", nulljson)};
  SkyView skyview(pJson.get("parameters", nulljson));

  const auto times = get_times(pSettings);

  std::cout << "Number of jobs = " << options.jobs << std::endl;

  boost::asio::thread_pool pool(options.jobs);

  long nlocations = pQ.Info()->SizeLocations();

  // for tracking how many points are skipped due to missing data, atomic_int is thread safe
  std::atomic_int missing_data_points{0};

  // set ratio of allowed missing data (e.g 0.1 means 10 %)
  auto jallowed_missing_ratio = Json::Path(".points.allow_missing_ratio").resolve(pJson);

  double allowed_missing_ratio;
  if (jallowed_missing_ratio.isNull())
    allowed_missing_ratio = 0.1;
  else
    allowed_missing_ratio = jallowed_missing_ratio.asDouble();

  // Establish common write stride only once since NFmiMetTime construction
  // is expensive due to libc TZ mutexes

  NFmiFastQueryInfo tmpinfo(&pQ);
  auto stride = get_write_stride(tmpinfo, pSettings);

  // Static thread locals don't need capture, and plain thread locals cannot be captured.
  // Hence a static thread_local here, and no "&qi" in the lambda capture list. Tried to
  // make write_timeseries thread safe using a single NFmiFastQueryInfo, but failed. All
  // methods used seem thread safe, but I must have missed something. Having separate
  // copies for each thread is especially expensive if the data is not gridded, and
  // thousands of station structs have to be copied.

  static thread_local std::unique_ptr<NFmiFastQueryInfo> qi;

  // Make sure NFmiQueryData latlon cache is initialized single threaded
  qi.reset(new NFmiFastQueryInfo(&pQ));
  qi->LatLon(0);
  auto mask = Json::Path(".points.mask").resolve(pJson);
  const auto point_mask = read_mask(mask, &tmpinfo);

  std::atomic_int simulation_counter{0};
  const auto true_location_count = std::count(point_mask.begin(), point_mask.end(), true);

  for (long loc_index = 0; loc_index < nlocations; ++loc_index)
  {
    if (!point_mask.at(loc_index))
      continue;

    boost::asio::post(
        pool,
        [loc_index,
         nlocations,
         allowed_missing_ratio,
         &missing_data_points,
         &simulation_counter,
         &true_location_count,
         &pQ,
         &stride,
         &times,
         &pSettings,
         &parameters,
         &skyview,
         &pDataManager]()
        {
          if (!qi)
            qi.reset(new NFmiFastQueryInfo(&pQ));

          LocalParameters lParameters(pSettings);
          auto lonlat = qi->LatLon(loc_index);

          // Get input data for point. Returns nothing if missing value is found in the data
          auto input =
              read_input(lonlat, times, pSettings, pDataManager, parameters, lParameters, skyview);

          if (input)
          {
            // initialize output vectors with missing values
            OutputData output(pSettings.SimLen);
            InputData inputT = *input;
            // Pass the data to Fortran via pointers
            auto outPointers = output.pointers();
            auto inPointers = inputT.pointers();
            runsimulation(&outPointers,
                          &inPointers,
                          &pSettings,
                          &parameters,
                          &lParameters);
            write_timeseries(*qi, loc_index, inputT, output, stride);
          }
          else
          {
            missing_data_points++;
            // If there are too many missing data points, throw error
            if (1.0 * missing_data_points / nlocations > allowed_missing_ratio)
              throw std::runtime_error("Too much missing data");
          }

          ++simulation_counter;
          if (simulation_counter % 1000 == 0)
            std::cout << "\tcompleted " << simulation_counter << " / " << true_location_count
                      << " points" << std::endl;
        });
  }

  pool.join();

  std::string filename = (options.outfile.empty() ? jfilename.asString() : options.outfile);
  pQ.Write(filename);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for querydata locations
 */
// ----------------------------------------------------------------------

void run_querydata_locations(NFmiQueryData& pQ,
                             const Json::Value& pJson,
                             InputSettings& pSettings,
                             const DataManager& pDataManager)
{
  if (options.jobs > 1)
    run_querydata_locations_async(pQ, pJson, pSettings, pDataManager);
  else
    run_querydata_locations_sync(pQ, pJson, pSettings, pDataManager);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for several locations
 */
// ----------------------------------------------------------------------

void run_coordinates(const Json::Value& pJson,
                     InputSettings& pSettings,
                     const DataManager& pDataManager)
{
  std::cout << "Running multiple coordinate simulation" << std::endl;

  // Establish the simulation points: latlon, numeric id and name

  auto qd = create_coordinates_querydata(pSettings, pJson);

  run_querydata_locations(*qd, pJson, pSettings, pDataManager);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model for grid
 */
// ----------------------------------------------------------------------

void run_grid(const Json::Value& pJson,
              InputSettings& pSettings,
              const DataManager& pDataManager)
{
  std::cout << "Running gridded simulation" << std::endl;

  auto qd = create_grid_querydata(pSettings, pJson);
  run_querydata_locations(*qd, pJson, pSettings ,pDataManager);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model
 */
// ----------------------------------------------------------------------

void run(const Json::Value& pJson, InputSettings& pSettings,  const DataManager& pDataManager)
{
  switch (pointmode(pJson))
  {
    case PointMode::Coordinate:
      run_coordinate(pJson, pSettings, pDataManager);
      break;
    case PointMode::Coordinates:
      run_coordinates(pJson, pSettings, pDataManager);
      break;
    case PointMode::Grid:
      run_grid(pJson, pSettings, pDataManager);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Add --infile option to JSON config
 */
// ----------------------------------------------------------------------

void add_option_infile(Json::Value& pJson)
{
  if (options.infile.empty())
    return;

  // We accept a comma separated list of files
  std::list<std::string> infiles;
  boost::algorithm::split(infiles, options.infile, boost::algorithm::is_any_of(","));

  int counter = 1;
  for (const auto& infile : infiles)
  {
    Json::Value jsource(Json::objectValue);
    jsource["name"] = "command line option " + Fmi::to_string(counter++);
    jsource["type"] = "ascii";
    jsource["filename"] = infile;

    if (!pJson.isMember("input"))
    {
      Json::Value jarray(Json::arrayValue);
      jarray.append(jsource);
      pJson["input"] = jarray;
    }
    else if (!pJson["input"].isArray())
      throw std::runtime_error("'input' object must be an array in the JSON configuration");
    else
      pJson["input"].append(jsource);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Add command line options to JSON configuration
 */
// ----------------------------------------------------------------------

void add_options(Json::Value& pJson)
{
  add_option_infile(pJson);
}

// Return true if there is a file containing string in the target directory
bool fileContainingStringExists(std::string target_directory, InputSettings& settings)
{
  std::string itrFile;
  std::string target_substring = Fmi::to_iso_string(settings.forecast_time);
  if (!exists(target_directory))
    return false;
  directory_iterator end_itr;
  for (directory_iterator itr(target_directory); itr != end_itr; ++itr)

  {
    itrFile = itr->path().filename().string();
    if (itrFile.find(target_substring) != std::string::npos)
      return true;
  }
  return false;
}
// ----------------------------------------------------------------------
/*!
 * \brief Main program
 */
// ----------------------------------------------------------------------
int main(int argc, char* argv[])
try  // NOLINT(cppcheck-syntaxError)
{
  if (!parse_options(argc, argv))
    return 0;

  // Read the JSON configurationi
  const Json::Value nulljson;
  auto json = read_json(options.configfile);
  add_options(json);
  // Initialize settings to default values

  InputSettings settings{json, options};
  DataManager datamanager;

  // If input data is not archived
  if (options.archive.empty() || !fileContainingStringExists(options.archive, settings))
  {
    // Start the data manager with wall clock time

    datamanager.init(settings.forecast_time, settings.start_time, settings.end_time, json);
  }
  // Various run-modes:
  //
  // - Single coordinate, ascii output, mostly for development purposes
  // - Multiple unstructured coordinates, 1D querydata output
  // - Gridded points, 2D querydata output

  run(json, settings, datamanager);

  return 0;
}
catch (const std::exception& e)
{
  std::cerr << "Error: " << e.what() << std::endl;
  return 1;
}
