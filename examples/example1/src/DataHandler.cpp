#include "DataHandler.h"
#include "GenericSourceFactory.h"

// ----------------------------------------------------------------------
/*!
 * \brief Construct configured data sources
 *
 * Sample settings:
 *
 *   "time":
 *   {
 *       "analysis": 48,
 *       "forecast": 120
 *   },
 *   "input":
 *   [
 *     {
 *        "name": "RWSObs",
 *        "path": "example_json",
 *        "type": "json",
 *        "source": "observations"
          "params":["Temperature 2m","Humidity","WindSpeed","Precipitation"]
 *     },
 *     ....
 *   ]
 *
 */
// ----------------------------------------------------------------------

//Destructor
DataHandler::~DataHandler() { for(auto * ptr : mDataSources) delete ptr; }

//Read in data from different data sources given in configuration file
void DataHandler::init(const time_t& pWallClock,
                       const time_t& pStartTime,
                       const time_t& pEndTime,
                       const Json::Value& pConfig,
                       const int SimLen,
                       const int DTSecs)
{
  mForecastTime = pWallClock; //Time when forecast phase starts
  mStartTime = pStartTime;    //Simulation start time
  mEndTime = pEndTime;        //Simulation end time

  // Read Data source settings
  Json::Value nulljson;

  auto json = pConfig.get("input", nulljson);
  if (json.isNull())
    throw std::runtime_error("Config variable 'input' must be set");
  if (!json.isArray())
    throw std::runtime_error(
        "Config variable 'input' must be an array of JSON objects defining data sources");
  //Go trough data sources
  for (const auto& source : json)
  {
    //Create Data source, the input data is read when the source is created
    mDataSources.push_back(
        GenericSourceFactory::create(mForecastTime, mStartTime, mEndTime,
	                            source, pConfig,SimLen,DTSecs));
    //Check if source is set as observation data (Used for determining 
    //last observation time)
    if (source.get("source", nulljson) == "observations")
      mDataSources.back()->is_observation = true;
  }
}
// ----------------------------------------------------------------------
/*!
 * \brief Get weather data for point
 */
// ----------------------------------------------------------------------

void DataHandler::GetWeather(InputData& pData, int pointID) const
{
  //Go trough all data sources and overwrite pData with the last read data
  for (const auto& source : mDataSources)
  {
    source->GetWeather(pData, pointID);
    //std::cout<<pData<<std::endl;
  }
  //  std::exit(EXIT_FAILURE);
}
// ----------------------------------------------------------------------
/*!
 * \brief Get local paramters object for each point
 */
// ----------------------------------------------------------------------

std::vector<LocalParameters> DataHandler::GetLocations() const
{
  return mDataSources[0]->GetLocations();
}
// ----------------------------------------------------------------------
/*!
 * \brief Get point ids
 */
// ----------------------------------------------------------------------

std::vector<int> DataHandler::GetPointIDs() const
{
  return mDataSources[0]->GetPointIDs();
}
// ----------------------------------------------------------------------
/*!
 * \brief Get index of latest observation available based on air temperature
 */
// ----------------------------------------------------------------------

int DataHandler::GetLatestObsIndex(int pointID) const
{
  int maxIndex=-1;
  //Go trough data sources
  for (const auto& source : mDataSources)
  {
    //Check if observation
    if (source->is_observation)
    {
      //Get index of latest observation
      int tmp = source->GetLatestObsIndex(pointID);
      //remplace maxIndex if new index is bigger
      if (tmp>-1)
      {
        if (maxIndex<0 || (tmp > maxIndex))
          maxIndex = tmp;
      }
    }
  }

  return maxIndex;
}
