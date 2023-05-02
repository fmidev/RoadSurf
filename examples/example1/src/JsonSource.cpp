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
       const int SimLen,
       const int DTSecs);

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
  //Go trouhg points until raw data ends or end of simulation
  //is reached
  while (rawPos+1<rawLen and simPos<simLen)
  {  //Check if times are equal
     if (abs(simtime[simPos]-rawtime[rawPos])<0.01)
     {
         //If times are equal, insert raw data to interpolatd data
         //if value is not missing
         if (rawData.tair[rawPos]>-100.0)
            interpolatedData.tair[simPos]=rawData.tair[rawPos];
         if (rawData.tdew[rawPos]>-100.0)
            interpolatedData.tdew[simPos]=rawData.tdew[rawPos];
         if (rawData.VZ[rawPos]>-100.0)
            interpolatedData.VZ[simPos]=rawData.VZ[rawPos];
         if (rawData.Rhz[rawPos]>-100.0)
            interpolatedData.Rhz[simPos]=rawData.Rhz[rawPos];
         if (rawData.prec[rawPos]>-100.0)
            interpolatedData.prec[simPos]=rawData.prec[rawPos];
         if (rawData.SW[rawPos]>-100.0)
            interpolatedData.SW[simPos]=rawData.SW[rawPos];
         if (rawData.LW[rawPos]>-100.0)
            interpolatedData.LW[simPos]=rawData.LW[rawPos];
         if (rawData.SW_dir[rawPos]>-100.0)
            interpolatedData.SW_dir[simPos]=rawData.SW_dir[rawPos];
         if (rawData.LW_net[rawPos]>-1000.0)
            interpolatedData.LW_net[simPos]=rawData.LW_net[rawPos];
         if (rawData.TSurfObs[rawPos]>-100.0)
            interpolatedData.TSurfObs[simPos]=rawData.TSurfObs[rawPos];
         if (rawData.PrecPhase[rawPos]>-100.0)
            interpolatedData.PrecPhase[simPos]=rawData.PrecPhase[rawPos];
         simPos++;
     //If simtime equals next time in rawtime, advance positin in raw data
     }else if (abs(simtime[simPos]-rawtime[rawPos+1])<0.01)
     {
         rawPos++;
     }else 
     {
         //If neither raw data point is missing, do the interpolation
         if (rawData.tair[rawPos]>-100.0 and rawData.tair[rawPos+1]>-100.0)
            interpolatedData.tair[simPos] = rawData.tair[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.tair[rawPos+1]-rawData.tair[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.tdew[rawPos]>-100.0 and rawData.tdew[rawPos+1]>-100.0)
            interpolatedData.tdew[simPos] = rawData.tdew[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.tdew[rawPos+1]-rawData.tdew[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.VZ[rawPos]>-100.0 and rawData.VZ[rawPos+1]>-100.0)
            interpolatedData.VZ[simPos] = rawData.VZ[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.VZ[rawPos+1]-rawData.VZ[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.Rhz[rawPos]>-100.0 and rawData.Rhz[rawPos+1]>-100.0)
            interpolatedData.Rhz[simPos] = rawData.Rhz[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.Rhz[rawPos+1]-rawData.Rhz[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.prec[rawPos]>-100.0 and rawData.prec[rawPos+1]>-100.0)
            interpolatedData.prec[simPos] = rawData.prec[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.prec[rawPos+1]-rawData.prec[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.SW[rawPos]>-100.0 and rawData.SW[rawPos+1]>-100.0)
            interpolatedData.SW[simPos] = rawData.SW[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.SW[rawPos+1]-rawData.SW[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.LW[rawPos]>-100.0 and rawData.LW[rawPos+1]>-100.0)
            interpolatedData.LW[simPos] = rawData.LW[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.LW[rawPos+1]-rawData.LW[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.SW_dir[rawPos]>-100.0 and rawData.SW_dir[rawPos+1]>-100.0)
            interpolatedData.SW_dir[simPos] = rawData.SW_dir[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.SW_dir[rawPos+1]-rawData.SW_dir[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.LW_net[rawPos]>-1000.0 and rawData.LW_net[rawPos+1]>-1000.0)
            interpolatedData.LW_net[simPos] = rawData.LW_net[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.LW_net[rawPos+1]-rawData.LW_net[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.TSurfObs[rawPos]>-100.0 and rawData.TSurfObs[rawPos+1]>-100.0)
            interpolatedData.TSurfObs[simPos] = rawData.TSurfObs[rawPos]+
                (simtime[simPos]-rawtime[rawPos])*
                (rawData.TSurfObs[rawPos+1]-rawData.TSurfObs[rawPos])/
                (rawtime[rawPos+1]-rawtime[rawPos]);
         if (rawData.PrecPhase[rawPos+1]>-100.0) 
            interpolatedData.PrecPhase[simPos] = rawData.PrecPhase[rawPos+1];

         simPos++;
     }
  }
}
// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

JsonSource::Impl::Impl(const time_t& pWallClock,
                           const time_t& pStartTime,
                           const time_t& pEndTime,
                           const Json::Value& pSource,
                           const Json::Value& pConfig,
                           const int SimLen,
                           const int DTSecs)
{
  //Variable names in json data
  const int Nvar=11;
  std::string variableNames[Nvar]={"Temperature 2m","Humidity","DewPoint","WindSpeed",
                "PrecipitationForm", "Precipitation","RadiationNetSurfaceLW",
                "RadiationLW","RadiationGlobal","RadiationDirectSW","RoadTemperature"};
  int simPos=0;
  //Initialize vector containing simulation times
  std::vector<time_t> simTime;
  auto t = pStartTime; //start time
  struct tm *tt;
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
    std::vector<time_t> rawtime;
    //Read times to vector from the data
    for (int index=0; index<dataLen; index++)
    {
       struct tm timeStruct;
       std::string timeString=dataArray[index].asString();
       strptime(timeString.c_str(),"%Y-%m-%d %H:%M",&timeStruct);
       timeStruct.tm_sec=0;
       timeStruct.tm_isdst=-1;
       rawtime.push_back(mktime(&timeStruct));
    }
    //read data for each variables
    for (int var=0; var<Nvar; var++)
    {
       Json::Value nulljson;
       //Get data array from json
       auto dataArray=j.get(variableNames[var],nulljson);
       if (!dataArray.isNull())
       { 
          //Reach data for each time step
          for (int index=0; index<dataLen; index++)
          {
              double value=dataArray[index].asDouble();
              if (var==0)
                 rawData.tair[index]=value;
              if (var==1)
                 rawData.Rhz[index]=value;
              if (var==2)
                 rawData.tdew[index]=value;
              if (var==3)
                 rawData.VZ[index]=value;
              if (var==4)
                 rawData.PrecPhase[index]=value;
              if (var==5)
                 rawData.prec[index]=value;
              if (var==6)
                 rawData.LW_net[index]=value;
              if (var==7)
                 rawData.LW[index]=value;
              if (var==8)
                 rawData.SW[index]=value;
              if (var==9)
                 rawData.SW_dir[index]=value;
              if (var==10)
                 rawData.TSurfObs[index]=value;
          }
       }
     }
    //Calculate tdew from rh or other way around if only one was in the input data
    for (int index=0; index<dataLen; index++)
    {
        if (rawData.tdew[index]<-100 and rawData.Rhz[index]>-100 and rawData.tair[index]>-100)
           rawData.tdew[index] = CalcTdewOrRH(rawData.tair[index], -9999.9, rawData.Rhz[index]);
        if (rawData.Rhz[index]<-100 and rawData.tdew[index]>-100 and rawData.tair[index]>-100)
           rawData.Rhz[index] = CalcTdewOrRH(rawData.tair[index], rawData.tdew[index], -9999.9);
    }
    //Fill the time arrays in interpolated data
    for (simPos=0; simPos < SimLen; ++simPos)
    {
        tt=localtime(&simTime[simPos]);
        tt->tm_year+=1900;
        tt->tm_mon+=1;
        interpolatedData.year[simPos]=tt->tm_year;
        interpolatedData.month[simPos]=tt->tm_mon;
        interpolatedData.day[simPos]=tt->tm_mday;
        interpolatedData.hour[simPos]=tt->tm_hour;
        interpolatedData.minute[simPos]=tt->tm_min;
        interpolatedData.second[simPos]=tt->tm_sec;
    }
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
     //Lenght of simulation
     int SimLen = static_cast<int>(weatherData[statPos].tair.size()); 
     int simPos=0; //Index in data
     //Add data to pData point by point if not missing
     for (simPos=0; simPos < SimLen; ++simPos)
     {
        if (weatherData[statPos].tair[simPos]>-100.0)
           pData.tair[simPos]=weatherData[statPos].tair[simPos];
        if (weatherData[statPos].tdew[simPos]>-100.0)
           pData.tdew[simPos]=weatherData[statPos].tdew[simPos];
        if (weatherData[statPos].VZ[simPos]>-100.0)
           pData.VZ[simPos]=weatherData[statPos].VZ[simPos];
        if (weatherData[statPos].Rhz[simPos]>-100.0)
           pData.Rhz[simPos]=weatherData[statPos].Rhz[simPos];
        if (weatherData[statPos].prec[simPos]>-100.0)
           pData.prec[simPos]=weatherData[statPos].prec[simPos];
        if (weatherData[statPos].SW[simPos]>-100.0)
           pData.SW[simPos]=weatherData[statPos].SW[simPos];
        if (weatherData[statPos].LW[simPos]>-100.0)
           pData.LW[simPos]=weatherData[statPos].LW[simPos];
        if (weatherData[statPos].SW_dir[simPos]>-100.0)
           pData.SW_dir[simPos]=weatherData[statPos].SW_dir[simPos];
        if (weatherData[statPos].LW_net[simPos]>-1000.0)
           pData.LW_net[simPos]=weatherData[statPos].LW_net[simPos];
        if (weatherData[statPos].TSurfObs[simPos]>-100.0)
           pData.TSurfObs[simPos]=weatherData[statPos].TSurfObs[simPos];
   
        if (weatherData[statPos].year[simPos]>-100.0)
           pData.year[simPos]=weatherData[statPos].year[simPos];
        if (weatherData[statPos].month[simPos]>-100.0)
           pData.month[simPos]=weatherData[statPos].month[simPos];
        if (weatherData[statPos].day[simPos]>-100.0)
           pData.day[simPos]=weatherData[statPos].day[simPos];
        if (weatherData[statPos].hour[simPos]>-100.0)
           pData.hour[simPos]=weatherData[statPos].hour[simPos];
        if (weatherData[statPos].minute[simPos]>-100.0)
           pData.minute[simPos]=weatherData[statPos].minute[simPos];
        if (weatherData[statPos].second[simPos]>-100.0)
           pData.second[simPos]=weatherData[statPos].second[simPos];
     }
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
