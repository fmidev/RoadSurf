#include <cmath>
#include <cstdlib>
#include "JsonSource.h"
#include "JsonTools.h"
#include "InputData.h"
#include "MeteorologyTools.h"
#include "LocalParameters.h"
#include <json/json.h>
#include <json/reader.h>
#include <time.h>
#include <algorithm>

//Class for reading JSon input
class JsonSource::Impl
{
 public:
  Impl(const time_t& pWallClock,
       const time_t& pStartTime,
       const time_t& pEndTime,
       const Json::Value& pSource,
       const Json::Value& pConfig,
       int SimLen,
       int DTSecs);

  //Get input data for point
  void GetWeather(InputData& pData, int pointID) const;
  //Get vector of class containing lat and lon
  std::vector<LocalParameters> GetLocations() const;
  //Get vector containing point ids
  std::vector<int> GetPointIDs() const;
  //Get index of latest observation
  int GetLatestObsIndex(int pointID) const;

 private:
   std::vector<int> pointIDs;
   std::vector<InputData> weatherData;
   std::vector<LocalParameters> locations;
   void interpolate(InputData rawData,
                 InputData& interpolatedData,
                 std::vector<time_t> rawtime,
                 std::vector<time_t> simtime);
};

JsonSource::~JsonSource() = default;

namespace
{
// Copy the value if it is above the missing value limit
template <typename T>
void copy_if_above(
    const std::vector<T>& src, std::vector<T>& dst, int srcPos, int dstPos, double limit)
{
  if (src[srcPos] > limit)
    dst[dstPos] = src[srcPos];
}

// Interpolate linearly in time if neither raw value is missing
void interpolate_if_valid(const std::vector<double>& raw,
                          std::vector<double>& interpolated,
                          int rawPos,
                          int simPos,
                          const std::vector<time_t>& rawtime,
                          const std::vector<time_t>& simtime,
                          double limit)
{
  if (raw[rawPos] > limit and raw[rawPos + 1] > limit)
    interpolated[simPos] = raw[rawPos] + (simtime[simPos] - rawtime[rawPos]) *
                                             (raw[rawPos + 1] - raw[rawPos]) /
                                             (rawtime[rawPos + 1] - rawtime[rawPos]);
}

// Copy the raw data at an equal time if the values are not missing
void copy_raw_values(const InputData& rawData, InputData& interpolatedData, int rawPos, int simPos)
{
  copy_if_above(rawData.tair, interpolatedData.tair, rawPos, simPos, -100.0);
  copy_if_above(rawData.tdew, interpolatedData.tdew, rawPos, simPos, -100.0);
  copy_if_above(rawData.VZ, interpolatedData.VZ, rawPos, simPos, -100.0);
  copy_if_above(rawData.Rhz, interpolatedData.Rhz, rawPos, simPos, -100.0);
  copy_if_above(rawData.prec, interpolatedData.prec, rawPos, simPos, -100.0);
  copy_if_above(rawData.SW, interpolatedData.SW, rawPos, simPos, -100.0);
  copy_if_above(rawData.LW, interpolatedData.LW, rawPos, simPos, -100.0);
  copy_if_above(rawData.SW_dir, interpolatedData.SW_dir, rawPos, simPos, -100.0);
  copy_if_above(rawData.LW_net, interpolatedData.LW_net, rawPos, simPos, -1000.0);
  copy_if_above(rawData.TSurfObs, interpolatedData.TSurfObs, rawPos, simPos, -100.0);
  copy_if_above(rawData.PrecPhase, interpolatedData.PrecPhase, rawPos, simPos, -100.0);
}

// Interpolate the raw data between two raw times
void interpolate_raw_values(const InputData& rawData,
                            InputData& interpolatedData,
                            int rawPos,
                            int simPos,
                            const std::vector<time_t>& rawtime,
                            const std::vector<time_t>& simtime)
{
  interpolate_if_valid(rawData.tair, interpolatedData.tair, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.tdew, interpolatedData.tdew, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.VZ, interpolatedData.VZ, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.Rhz, interpolatedData.Rhz, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.prec, interpolatedData.prec, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.SW, interpolatedData.SW, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(rawData.LW, interpolatedData.LW, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(
      rawData.SW_dir, interpolatedData.SW_dir, rawPos, simPos, rawtime, simtime, -100.0);
  interpolate_if_valid(
      rawData.LW_net, interpolatedData.LW_net, rawPos, simPos, rawtime, simtime, -1000.0);
  interpolate_if_valid(
      rawData.TSurfObs, interpolatedData.TSurfObs, rawPos, simPos, rawtime, simtime, -100.0);
  copy_if_above(rawData.PrecPhase, interpolatedData.PrecPhase, rawPos + 1, simPos, -100.0);
}

// Find the first raw and simulation positions to process
void find_start_positions(const std::vector<time_t>& rawtime,
                          const std::vector<time_t>& simtime,
                          int& rawPos,
                          int& simPos)
{
  //Original data length
  int rawLen = static_cast<int>(rawtime.size());
  //Simulation length
  int simLen = static_cast<int>(simtime.size());
  //if simulation time starts after raw data time
  if (rawtime[0]<simtime[0])
  {
     // Search for first raw time value before simTime
     for (rawPos=0; rawPos < rawLen; ++rawPos)
     {
        if (rawtime[rawPos]>=simtime[0])
           break;
     }
     rawPos=rawPos-1;
     simPos=0;
  //If simulation time starts before raw time
  } else if (simtime[0]<rawtime[0])
  {
     //Search for first simulation value before raw time
     for (simPos=0; simPos < simLen; ++simPos)
     {
        if (simtime[simPos]>=rawtime[0])
           break;
     }
     rawPos=0;
  }
}

// Read the times of the raw data
std::vector<time_t> read_raw_times(const Json::Value& dataArray, int dataLen)
{
  std::vector<time_t> rawtime;
  //Read times to vector from the data
  for (int index=0; index<dataLen; index++)
  {
     struct tm timeStruct{};
     std::string timeString=dataArray[index].asString();
     strptime(timeString.c_str(),"%Y-%m-%d %H:%M",&timeStruct);
     timeStruct.tm_sec=0;
     timeStruct.tm_isdst=-1;
     rawtime.push_back(mktime(&timeStruct));
  }
  return rawtime;
}

// Store a value of the variable with the given index in the variable name table
void store_raw_value(InputData& rawData, int var, int index, double value)
{
  switch (var)
  {
    case 0:
      rawData.tair[index] = value;
      break;
    case 1:
      rawData.Rhz[index] = value;
      break;
    case 2:
      rawData.tdew[index] = value;
      break;
    case 3:
      rawData.VZ[index] = value;
      break;
    case 4:
      rawData.PrecPhase[index] = value;
      break;
    case 5:
      rawData.prec[index] = value;
      break;
    case 6:
      rawData.LW_net[index] = value;
      break;
    case 7:
      rawData.LW[index] = value;
      break;
    case 8:
      rawData.SW[index] = value;
      break;
    case 9:
      rawData.SW_dir[index] = value;
      break;
    case 10:
      rawData.TSurfObs[index] = value;
      break;
    default:
      break;
  }
}

// Read the data for each variable
void read_raw_values(const Json::Value& j,
                     InputData& rawData,
                     int dataLen,
                     const std::string* variableNames,
                     int Nvar)
{
  for (int var=0; var<Nvar; var++)
  {
     Json::Value nulljson;
     //Get data array from json
     auto dataArray=j.get(variableNames[var],nulljson);
     if (!dataArray.isNull())
     {
        //Reach data for each time step
        for (int index=0; index<dataLen; index++)
          store_raw_value(rawData, var, index, dataArray[index].asDouble());
     }
  }
}

// Calculate tdew from rh or other way around if only one was in the input data
void calculate_humidity(InputData& rawData, int dataLen)
{
  for (int index=0; index<dataLen; index++)
  {
      if (rawData.tdew[index]<-100 and rawData.Rhz[index]>-100 and rawData.tair[index]>-100)
         rawData.tdew[index] = CalcTdewOrRH(rawData.tair[index], -9999.9, rawData.Rhz[index]);
      if (rawData.Rhz[index]<-100 and rawData.tdew[index]>-100 and rawData.tair[index]>-100)
         rawData.Rhz[index] = CalcTdewOrRH(rawData.tair[index], rawData.tdew[index], -9999.9);
  }
}

// Fill the time arrays in interpolated data
void fill_simulation_times(InputData& interpolatedData,
                           const std::vector<time_t>& simTime,
                           int SimLen)
{
  struct tm *tt = nullptr;
  for (int simPos=0; simPos < SimLen; ++simPos)
  {
      struct ::tm ttBuffer{};  // localtime_r is thread safe
      tt=localtime_r(&simTime[simPos], &ttBuffer);
      tt->tm_year+=1900;
      tt->tm_mon+=1;
      interpolatedData.year[simPos]=tt->tm_year;
      interpolatedData.month[simPos]=tt->tm_mon;
      interpolatedData.day[simPos]=tt->tm_mday;
      interpolatedData.hour[simPos]=tt->tm_hour;
      interpolatedData.minute[simPos]=tt->tm_min;
      interpolatedData.second[simPos]=tt->tm_sec;
  }
}

// Copy the data of a station to the input data if not missing
void copy_station_values(const InputData& w, InputData& pData)
{
  //Lenght of simulation
  int SimLen = static_cast<int>(w.tair.size());
  int simPos=0; //Index in data
  //Add data to pData point by point if not missing
  for (simPos=0; simPos < SimLen; ++simPos)
  {
    copy_if_above(w.tair, pData.tair, simPos, simPos, -100.0);
    copy_if_above(w.tdew, pData.tdew, simPos, simPos, -100.0);
    copy_if_above(w.VZ, pData.VZ, simPos, simPos, -100.0);
    copy_if_above(w.Rhz, pData.Rhz, simPos, simPos, -100.0);
    copy_if_above(w.prec, pData.prec, simPos, simPos, -100.0);
    copy_if_above(w.SW, pData.SW, simPos, simPos, -100.0);
    copy_if_above(w.LW, pData.LW, simPos, simPos, -100.0);
    copy_if_above(w.SW_dir, pData.SW_dir, simPos, simPos, -100.0);
    copy_if_above(w.LW_net, pData.LW_net, simPos, simPos, -1000.0);
    copy_if_above(w.TSurfObs, pData.TSurfObs, simPos, simPos, -100.0);

    copy_if_above(w.year, pData.year, simPos, simPos, -100.0);
    copy_if_above(w.month, pData.month, simPos, simPos, -100.0);
    copy_if_above(w.day, pData.day, simPos, simPos, -100.0);
    copy_if_above(w.hour, pData.hour, simPos, simPos, -100.0);
    copy_if_above(w.minute, pData.minute, simPos, simPos, -100.0);
    copy_if_above(w.second, pData.second, simPos, simPos, -100.0);
  }
}
}  // namespace


// ----------------------------------------------------------------------
/*!
 * \brief Interpolates raw data for each time step in simulation
 */
// ----------------------------------------------------------------------
void JsonSource::Impl::interpolate(InputData rawData,
                 InputData& interpolatedData,
                 std::vector<time_t> rawtime,
                 std::vector<time_t> simtime)
{
  //Original data length
  int rawLen = static_cast<int>(rawtime.size());
  //Simulation length
  int simLen = static_cast<int>(simtime.size());
  int rawPos=0;
  int simPos=0;
  find_start_positions(rawtime, simtime, rawPos, simPos);
  //Go trouhg points until raw data ends or end of simulation
  //is reached
  while (rawPos+1<rawLen and simPos<simLen)
  {  //Check if times are equal
     if (std::abs(simtime[simPos]-rawtime[rawPos])<0.01)
     {
         //If times are equal, insert raw data to interpolatd data
         //if value is not missing
         copy_raw_values(rawData, interpolatedData, rawPos, simPos);
         simPos++;
     //If simtime equals next time in rawtime, advance positin in raw data
     }else if (std::abs(simtime[simPos]-rawtime[rawPos+1])<0.01)
     {
         rawPos++;
     }else
     {
         //If neither raw data point is missing, do the interpolation
         interpolate_raw_values(rawData, interpolatedData, rawPos, simPos, rawtime, simtime);
         simPos++;
     }
  }
}
// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

JsonSource::Impl::Impl(const time_t& /*pWallClock*/,
                           const time_t& pStartTime,
                           const time_t& /*pEndTime*/,
                           const Json::Value& pSource,
                           const Json::Value& /*pConfig*/,
                           const int SimLen,
                           const int DTSecs)
{
  //Variable names in json data
  const int Nvar=11;
  std::string variableNames[Nvar]={"Temperature 2m","Humidity","DewPoint","WindSpeed",
                "PrecipitationForm", "Precipitation","RadiationNetSurfaceLW",
                "RadiationLW","RadiationGlobal","RadiationDirectSW","RoadTemperature"};
  //Initialize vector containing simulation times
  std::vector<time_t> simTime;
  auto t = pStartTime; //start time
  for (int i = 0; i<SimLen; i++)
  {
    simTime.emplace_back(t);
    t += DTSecs;
  }
  //Read input file path from json file and read the file
  Json::Value nulljson;
  auto path = pSource.get("path", nulljson).asString();
  auto json = read_json(path);

  //Go trough locations
  for (const auto& j : json)
  {
    //Station identification number. It is used to identify the
    //right input data for each point
    int id = j["statId"].asInt();
    pointIDs.push_back(id);

    //Class containing latitude and longitude
    LocalParameters localParams;
    double lon = j["lon"].asDouble();
    double lat = j["lat"].asDouble();
    localParams.lat=lat;
    localParams.lon=lon;
    locations.push_back(localParams);
    //Initialize interpolated data
    InputData interpolatedData{SimLen};

    //Determine input data length
    auto dataArray=j.get("time",nulljson);
    int dataLen = static_cast<int>(dataArray.size()); 
    
    if (dataLen==0) 
    {
       weatherData.push_back(interpolatedData);
       continue;
    }
    //Initialize input data class
    InputData rawData{dataLen};
    //vecotr for input data times
    std::vector<time_t> rawtime = read_raw_times(dataArray, dataLen);
    //read data for each variables
    read_raw_values(j, rawData, dataLen, variableNames, Nvar);
    //Calculate tdew from rh or other way around if only one was in the input data
    calculate_humidity(rawData, dataLen);
    //Fill the time arrays in interpolated data
    fill_simulation_times(interpolatedData, simTime, SimLen);
    //Do interpolation from original data times to simulation times
    interpolate(rawData,interpolatedData,rawtime,simTime);
    //Add data to vector
    weatherData.push_back(interpolatedData);
   }

}
// ----------------------------------------------------------------------
/*!
 * \brief Get data data for point
 */
// ----------------------------------------------------------------------
void JsonSource::Impl::GetWeather(InputData& pData,int pointID) const
{
  //Find index in the pointID vector with the given pointID
  auto index=find(pointIDs.begin(),pointIDs.end(),pointID);
  //If index is found
  if (index!=pointIDs.end())
  {
     int statPos=index-pointIDs.begin();//This gives the right index for point
     copy_station_values(weatherData[statPos], pData);
  }
}
// ----------------------------------------------------------------------
/*!
 * \brief Get vector with location data for each point
 */
// ----------------------------------------------------------------------
std::vector<LocalParameters> JsonSource::Impl::GetLocations() const
{
  return locations;
}
// ----------------------------------------------------------------------
/*!
 * \brief Get vector of pointIDs
 */
// ----------------------------------------------------------------------
std::vector<int> JsonSource::Impl::GetPointIDs() const
{
  return pointIDs;
}

// ----------------------------------------------------------------------
/*!
 * \brief Get latest available road observation time 
 */
// ----------------------------------------------------------------------

int JsonSource::Impl::GetLatestObsIndex(int pointID) const
{
  //Gind inidex of statID in the vector
  auto index=find(pointIDs.begin(),pointIDs.end(),pointID);
  if (index!=pointIDs.end())
  {
     int statPos=index-pointIDs.begin(); //Get the right index
     int SimLen = static_cast<int>(weatherData[statPos].tair.size()); 
     //Go trough data backwards until first observation is found
     for (int i = SimLen; i > 0; i--)
     {
       if (weatherData[statPos].tair[i - 1]>-100)
          return i;
     }
     return -9999;
   }
   return -9999;
}
// ----------------------------------------------------------------------
/*!
 * \brief Construct a JsonSource source for Json data
 */
// ----------------------------------------------------------------------

JsonSource::JsonSource(const time_t& pWallClock,
                               const time_t& pStartTime,
                               const time_t& pEndTime,
                               const Json::Value& pSource,
                               const Json::Value& pConfig,
                               const int SimLen,
                               const int DTSecs)
    : impl(new JsonSource::Impl(pWallClock, pStartTime, pEndTime, pSource, pConfig,SimLen,DTSecs))
{
}
// ----------------------------------------------------------------------
/*!
 * \brief Get input data for the given point 
 */
// ----------------------------------------------------------------------

void JsonSource::GetWeather(InputData& pData, int pointID) const

{
  impl->GetWeather(pData, pointID);
}
// ----------------------------------------------------------------------
/*!
 * \brief Get point locations
 */
// ----------------------------------------------------------------------

std::vector<LocalParameters> JsonSource::GetLocations() const

{
  return impl->GetLocations();
}
// ----------------------------------------------------------------------
/*!
 * \brief point ids
 */
// ----------------------------------------------------------------------
std::vector<int> JsonSource::GetPointIDs() const
{
  return impl->GetPointIDs();
}
// ----------------------------------------------------------------------
/*!
 * \brief Get index of last observation
 */
// ----------------------------------------------------------------------

int JsonSource::GetLatestObsIndex(int pointID) const
{
   return impl->GetLatestObsIndex(pointID);
}
