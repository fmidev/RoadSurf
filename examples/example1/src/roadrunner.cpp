#include "DataHandler.h"
#include "InputData.h"
#include "InputSettings.h"
#include "InputParameters.h"
#include "JsonTools.h"
#include "LocalParameters.h"
#include "OutputData.h"
#include "WorkQueue.h"
#include "SkyView.h"
#include <json/json.h>
#include <json/writer.h>
#include <atomic>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <time.h>
#include <iterator>
#include <unistd.h>
#include <sys/stat.h>
// The fortran API
extern "C"
{
  void runsimulation(OutputPointers* pOutputPointers,
                     const InputPointers* pInputPointers,
                     const InputSettings* pSettings,
                     const InputParameters* pInputParams,
                     const LocalParameters* lParameters); 
}
// Parsed command line options
Options options;

// Check if file exists
// std::filesystem::exist not available in C++11
inline bool file_exists (const std::string& name)
{
   struct stat buffer;
   return (stat (name.c_str(), &buffer) == 0);
}

//Check if value is missing
inline bool is_missing(double pValue)
{
  return (std::isnan(pValue) || pValue < -9000);
}
//Args structure required by multi thread run
struct Args{
   int loc_index;
};

// ----------------------------------------------------------------------
/*!
 * \brief Parse command line options
 * \return True on success, false if execution is to be stopped
 */
// ----------------------------------------------------------------------

bool parse_options(int argc, char* argv[])
{
  //forecast start time
  struct tm forecast_start_time;
  forecast_start_time.tm_isdst=-1;
  std::string stime;
  //Option descriptions
  std::vector<std::string> desc;
  desc.push_back("-h print help message");
  desc.push_back("-j number of parallel jobs");
  desc.push_back("-t simulation time");
  desc.push_back("-c configuration file");
  int descLen = static_cast<int>(desc.size());

  int opt; 
  //Go trough options
  while ((opt = getopt(argc,argv,"hj:t:c:")) !=-1)
  {
     switch (opt)
     {
     //print help message
     case 'h':
     
        std::cout << "Usage: roadrunner [options] [configfile]\n\n"
                 "Run road model using the given configuration.\n"
               <<std::endl;
        for (int j = 1;j<descLen;++j)
        {
              std::cout<< desc[j] << std::endl;
        }
        return false;
     //parse forecast start time
     case 't':
     
       stime=optarg;
       strptime(stime.c_str(),"%Y%m%dT%H%M",&forecast_start_time);
       forecast_start_time.tm_sec=0;
       options.time = mktime(&forecast_start_time);
       break;
     //number of cores used in multithread run
     case 'j':
        options.jobs = atoi(optarg);
	//use all cores if smaller than zero
        if (options.jobs < 0)
           options.jobs = std::thread::hardware_concurrency();
        break;
     //Configuration file
     case 'c': 
         options.configfile=optarg;
         break;
     default:
        std::cout << "Usage: roadrunner [options] [configfile]\n"<<std::endl;
        return false;
        
     }   
  }
  //Read configuration file as last argument if not already given
  if (options.configfile.empty() and optind<argc)
  {
     options.configfile=argv[optind];
  }
  //Throw error if configuration file is not given
  if (options.configfile.empty())
     throw std::runtime_error("Configuration file not given");

  //Check that config file exists
  if (!file_exists(options.configfile))
    throw std::runtime_error("Configuration file '" + options.configfile + "' missing");
  return true;
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct the simulation timesteps
 */
// ----------------------------------------------------------------------

std::vector<time_t> get_times(const InputSettings& pSettings)
{
  std::vector<time_t> times;

  const auto dt = pSettings.DTSecs; //time step
  auto t = pSettings.start_time;    //start time
  for (int i = 0; i < pSettings.SimLen; i++)
  {
    times.emplace_back(t);
    t += dt;
  }

  return times;
}


// ----------------------------------------------------------------------
/*!
 * \brief Read data required for simulation
 */
// ----------------------------------------------------------------------

InputData read_input(InputSettings& pSettings,
                                      const std::vector<time_t> pTimes,
                                      const DataHandler& pDataHandler,
                                      const InputParameters& parameters,
                                      LocalParameters& lParameters,
                                      SkyView& skyview,
                                      int pointID,
                                      bool& success)
{
    //Class for storing input data, initialized with missing values
    InputData data{pSettings.SimLen};
    // Set Initialization length used by relaxation
    const auto init_secs = pSettings.forecast_time - pSettings.start_time;
    lParameters.InitLenI = 1 + static_cast<int>(init_secs / pSettings.DTSecs);

    //Get data for the given point
    pDataHandler.GetWeather(data, pointID);
   
    //Get sky view factor and local horizon angles for given point
    skyview.GetSkyVariables(lParameters, data, pointID);
    //For convertin time to string
    const std::size_t MAXLEN=20;
    std::array<char,MAXLEN> timeString;
    
    std::string fmt="%Y-%m-%dT%H:%M";
    // Verify the model has all required parameters
    for (int i = 0; i < pSettings.SimLen; i++)
    {
      //convert time to string
      auto pTime=gmtime(&pTimes[i]);
      std::strftime(timeString.data(),MAXLEN,fmt.c_str(),pTime);
      if (is_missing(data.tair[i]))
      {
        std::cout << "2m temperature missing at " << timeString.data() << " " << lParameters.lat
                  << " " << lParameters.lon << std::endl;
        success=false;
        return data;
      }
      if (is_missing(data.Rhz[i]))
      {
        std::cout << "RH missing at " << timeString.data() << " " << lParameters.lat << " "
                  << lParameters.lon << std::endl;
        success=false;
        return data;
      }
      if (is_missing(data.prec[i]))
      {
        std::cout << "Precipitation missing at " << timeString.data() << " " << lParameters.lat
                  << " " << lParameters.lon << std::endl;
        success=false;
        return data;
      }
      if (is_missing(data.SW[i]))
      {
        std::cout << "SW missing at " << timeString.data() << " " << lParameters.lat << " "
                  << lParameters.lon << std::endl;
        success=false;
        return data;
      }
      if (is_missing(data.LW[i]))
      {
        std::cout << "LW missing at " << timeString.data() << " " << lParameters.lat << " "
                  << lParameters.lon << std::endl;
        success=false;
        return data;
      }
      if (is_missing(data.VZ[i]))
      {
        std::cout << "VZ missing at  " << timeString.data() << " " << lParameters.lat << " "
                  << lParameters.lon << std::endl;
        success=false;
        return data;
      }

    }
    // If relaxation is used, set first forecasted values
    //(Should be forecasted values at the time of the latest observation, but this is
    // good enough with small timesteps)
    if (pSettings.use_relaxation == 1)
    {
      lParameters.tair_relax = -9999.9;
      lParameters.VZ_relax = -9999.9;
      lParameters.RH_relax = -9999.9;
      //Get index of last observed air temperature value
      int lastTairObsIndex = pDataHandler.GetLatestObsIndex(pointID);
      if (lastTairObsIndex>-1)
      {
        lParameters.InitLenI = lastTairObsIndex; //length of initialization period
        lParameters.tair_relax = data.tair[lastTairObsIndex];
        lParameters.VZ_relax = data.VZ[lastTairObsIndex];
        lParameters.RH_relax = data.Rhz[lastTairObsIndex];
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
  success=true;
  return data;
}

// ----------------------------------------------------------------------
/*!
 * \brief Save output to json object
 */
// ----------------------------------------------------------------------
void save_output(OutputData& output, InputData input, LocalParameters& lParameters, 
                 int pointID, InputSettings& pSettings,std::vector<time_t> times,
                 int loc_index,Json::Value& forecast)
{
   //Output step
   int step=pSettings.outputStep*60/pSettings.DTSecs;
   //Values for time to string conversion
   const std::size_t MAXLEN=20;
   std::array<char,MAXLEN> timeString; 
   std::string fmt="%Y-%m-%dT%H:%M";
   //output data vectors
   Json::Value timeVec(Json::arrayValue);
   Json::Value tsurf(Json::arrayValue);
   Json::Value water(Json::arrayValue);
   Json::Value ice(Json::arrayValue);
   Json::Value snow(Json::arrayValue);
   Json::Value deposit(Json::arrayValue);
   //Go trhough data in steps determined by outptu step
   for (int i = 0; i < pSettings.SimLen; i+=step)
   {
      //time tos string
      auto pTime=gmtime(&times[i]);
      pTime->tm_isdst=-1;
      std::strftime(timeString.data(),MAXLEN,fmt.c_str(),pTime);
      //add data to vectors
      timeVec.append(timeString.data());
      tsurf.append(Json::Value(output.TsurfOut[i]));
      water.append(Json::Value(output.WaterOut[i]));
      ice.append(Json::Value(output.IceOut[i]));
      snow.append(Json::Value(output.SnowOut[i]));
      deposit.append(Json::Value(output.DepositOut[i]));
   }
   //Set data to json object
   forecast[loc_index]["statId"]=pointID;
   forecast[loc_index]["lat"]=lParameters.lat;
   forecast[loc_index]["lon"]=lParameters.lon;
   forecast[loc_index]["time"]=timeVec;
   forecast[loc_index]["RoadTemperature"]=tsurf;
   forecast[loc_index]["Water"]=water;
   forecast[loc_index]["Ice"]=ice;
   forecast[loc_index]["Snow"]=snow;
   forecast[loc_index]["Deposit"]=deposit;
}
// ----------------------------------------------------------------------
/*!
 * \brief Write json output
 */
// ----------------------------------------------------------------------

void write_output(std::string outputFileName, const Json::Value& forecast)
{
  std::ofstream ofs(outputFileName);
  Json::StreamWriterBuilder builder;
  //Some settings for writing
  builder["indentation"]="   ";
  builder.settings_["precision"]=7;
  //set precision to decimal, does not seem to work
  builder.settings_["precisionType"]="decimal";
  std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
  //Write json data
  writer->write(forecast,&ofs);
  ofs.close();
}

// ----------------------------------------------------------------------
/*!
 * \brief Run the points one by one
 */
// ----------------------------------------------------------------------

void run_locations_sync(const Json::Value& pJson,
                                  InputSettings& pSettings,
                                  const DataHandler& pDataHandler)
{
  // Check that output file is given
  auto jfilename = Json::Path(".output.filename").resolve(pJson);
  if (jfilename.isNull() && options.outfile.empty())
    throw std::runtime_error("Output filename not set");
  std::string outputFileName=jfilename.asString();
  Json::Value nulljson;
  //Set parameters
  InputParameters parameters{pSettings, pJson.get("parameters", nulljson)};
  //Read sky view factors and local horizon angels
  SkyView skyview(pJson.get("parameters", nulljson));
  //Get simulation times
  const auto times = get_times(pSettings);
  //Get vector of location information
  auto locations=pDataHandler.GetLocations();
  //Get vector of point ids
  auto pointIDs=pDataHandler.GetPointIDs();
  //Variable to determine wheter input data was succesfully read
  bool success=false;
   // Number of forecast locations
  int nlocations = static_cast<int>(locations.size());
  //Json object to store output forecast
  Json::Value forecast;  

  //Go trough all forecast locations
  for (int loc_index = 0; loc_index < nlocations; ++loc_index)
  {
     //Point location
     LocalParameters lParameters=locations[loc_index];
     success=false;
     //Get input data for point
     const auto input =
        read_input( pSettings, times, pDataHandler, parameters, lParameters, 
                  skyview, pointIDs[loc_index], success);
     //If input data was succesfully obtained
     if (success)
     {
       //Report progress every 1000 points
       if (loc_index % 1000 == 0)
          std::cout << "\t" << loc_index << " / " << nlocations << " points" << std::endl;
        // initialize output vectors with missing values
       OutputData output(pSettings.SimLen);
        // Pass the data to Fortran via pointers
       auto outPointers = output.pointers();
       auto inPointers = input.pointers();
       //Run the road weather model
       runsimulation(
            &outPointers, &inPointers, &pSettings,
            &parameters, &lParameters);
       //Save output to json object
       save_output(output, input, lParameters,pointIDs[loc_index], pSettings,times,
                 loc_index,forecast);

      }
  }
  //write output to file
  write_output(outputFileName,forecast);
}

// ----------------------------------------------------------------------
/*!
 * \brief Run points in parallel
 */
// ----------------------------------------------------------------------

void run_locations_async(const Json::Value& pJson,
                                   InputSettings& pSettings,
                                   const DataHandler& pDataHandler)
{
  // Check that output file is givne
  auto jfilename = Json::Path(".output.filename").resolve(pJson);
  if (jfilename.isNull())
    throw std::runtime_error("Output filename not set");
  std::string outputFileName=jfilename.asString();

  Json::Value nulljson;
  //Set input parameters
  InputParameters parameters{pSettings, pJson.get("parameters", nulljson)};
  //Get sky view factor and local horizon angles
  SkyView skyview(pJson.get("parameters", nulljson));
  //Get simulation times
  const auto times = get_times(pSettings);
  //Stae number of jobs done in parallel
  std::cout << "Number of jobs = " << options.jobs << std::endl;

  //Get vector with location information
  auto locations=pDataHandler.GetLocations();
  //Get vector with point ids
  auto pointIDs=pDataHandler.GetPointIDs();
  //number of forecast locations
  int nlocations = static_cast<int>(locations.size());
  // atomic int is thread safe
  std::atomic_int simulation_counter{0};
  Json::Value forecast;  

  //Function for running the road weather model
  const auto run_location=[nlocations,&simulation_counter,&times,&pSettings,
                     &parameters,&skyview,&pDataHandler,&locations,&pointIDs,
                     &forecast](Args args)
  {
       bool success=false;
       int pointID=pointIDs[args.loc_index];
       //Get location info for point
       LocalParameters lParameters=locations[args.loc_index];
       //Get input for point
       const auto input =
          read_input(pSettings, times, pDataHandler, parameters, lParameters, 
                    skyview,pointID, success);
       //If obtaining input was success
       if (success)
       {
          // initialize output vectors with missing values
         OutputData output(pSettings.SimLen);
          // Pass the data to Fortran via pointers
         auto outPointers = output.pointers();
         auto inPointers = input.pointers();
	 //Run the road weather model
         runsimulation(
              &outPointers, &inPointers, &pSettings,
              &parameters, &lParameters);
	 //save output to json object
         save_output(output, input, lParameters,pointIDs[args.loc_index], pSettings,times,
                   args.loc_index,forecast);
  
            ++simulation_counter;
	    //Report progress every 1000 points
            if (simulation_counter % 1000 == 0)
              std::cout << "\tcompleted " << simulation_counter << " / " << nlocations
                        << " points" << std::endl;
         }
  };
  //Work queue for running the simulations in parallel
  WorkQueue<Args> workqueue(run_location,options.jobs);
  //Go trouhg points
  for (int loc_index = 0; loc_index < nlocations; ++loc_index)
  {
       Args args{loc_index};
       workqueue(args); //send to workqueue
  }
  workqueue.join_all();
  //Write output to file
  write_output(outputFileName,forecast);

}

// ----------------------------------------------------------------------
/*!
 * \brief Run the model
 */
// ----------------------------------------------------------------------

void run(const Json::Value& pJson, InputSettings& pSettings, 
         const DataHandler& pDataHandler)
{
  //Run simulations in parallel if number of jobs is larger than one
  if (options.jobs > 1)
    run_locations_async(pJson, pSettings, pDataHandler);
  //Ohterwise run point by point
  else
    run_locations_sync(pJson, pSettings, pDataHandler);
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
  // Initialize settings to default values
  InputSettings settings{json, options};
  //Class for handling input data
  DataHandler datahandler;

  // Start the data handler
  datahandler.init(settings.forecast_time, settings.start_time, settings.end_time,
                     json,settings.SimLen,settings.DTSecs);
  //run the model
  run(json, settings, datahandler);

  return 0;
}
catch (const std::exception& e)
{
  std::cerr << "Error: " << e.what() << std::endl;
  return 1;
}
