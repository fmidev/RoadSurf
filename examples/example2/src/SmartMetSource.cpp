#include "SmartMetSource.h"
#include "FileTools.h"
#include "MeteorologyTools.h"
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
#include <numeric>
namespace
{
using Observations = std::vector<Weather>;    // observations sorted by time
using DataMap = std::map<int, Observations>;  // observations for each station

// Last searched location and the data for it is cached for speed
thread_local std::optional<NFmiPoint> myLatLon;
thread_local DataMap::const_iterator myData;

std::string format_smartmet_time(NFmiMetTime pTime, int pTimeMargin)
{
  pTime.ChangeByMinutes(pTimeMargin);

  std::string ret = pTime.ToStr(kYYYYMMDD).CharPtr();
  ret += 'T';
  ret += pTime.ToStr(kHour + kMinute).CharPtr();
  return ret;
}
}  // namespace

// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

class SmartMetSource::Impl
{
 public:
  Impl(const NFmiMetTime& pWallClock,
       const NFmiMetTime& pStartTime,
       const NFmiMetTime& pEndTime,
       const Json::Value& pSource,
       const Json::Value& pConfig);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon);

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& p) const;

 private:
  const NFmiMetTime mWallClock;
  std::vector<std::string> mParams;
  NFmiMetTime mStartTime;
  NFmiMetTime mEndTime;
  int mTimeMargin = 0;

  // We may use lazy initialization for coordinate input
  mutable bool mInitialized = false;

  mutable DataMap mData;

  using Station = Fmi::NearTreeLatLon<int>;
  using StationTree = Fmi::NearTree<Station, Fmi::NearTreeLatLonDistance<Station>>;
  mutable StationTree mStationTree;

  double mRadius;

  std::string mProtocol;
  std::string mHost;
  std::string mPlugin;
  std::optional<std::string> mKeyword;
  std::string mFmisid;
  std::string mProducer = "road";

  std::string mParameters;
  std::string mNumberF;
  std::string mDateF;
  std::string mTroadF;
  std::string mTroadSF;
  std::string mT2mF;
  std::string mTdewF;
  std::string mRhF;
  std::string mWspdF;
  std::string mLradF;
  std::string mSradF;
  std::string mSradDirF;
  std::string mLradNetF;
  std::string mRr1hF;
  std::string mRrF;
  std::string mRformF;

  void check_params();

  static double interpolate(const Observations& pObservations,
                            const NFmiMetTime& pTime,
                            int pPos,
                            int pMaxTimeGap,
                            double Weather::*pPtr);

  void read_keyword(const Json::Value& pConfig);
  void read_ids(const Json::Value& pConfig);
  void read_latlon(const NFmiPoint& pLatLon,const Json::Value& pConfig);
  void read_data(const cpr::Url& pUrl, const cpr::Parameters& pParams,const Json::Value& pConfig);
  bool uses_param(const std::string& pParam) const;
  void SameValueCheck(Observations& pObs,int sameLimit);
  void RemoveTooManySameValues( Observations& pObs,int sameLimit, double Weather::*pPtr);
  void RemoveSpikes(Observations& pObs,double spikeSize, std::size_t maxSpikeDuration,
		double Weather::*pPtr);
  void RemoveOutliers();
  void SpikeCheck(Observations& pObs,double spikeSize,std::size_t maxSpikeDuration);
  void MaxMinCheck(Observations& pObs);
};  // class SmartMetSource::Impl

// ----------------------------------------------------------------------
/*!
 * \brief Read observations for a keyword from the server
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::read_keyword(const Json::Value& pConfig)
{
  if (mInitialized || !mKeyword)
    return;

  std::cout << "Reading keyword " << *mKeyword << std::endl;

  // Date conditions

  auto mydate1 = format_smartmet_time(mStartTime, -mTimeMargin);
  //  auto mydate2 = format_smartmet_time(mEndTime, +mTimeMargin);
  auto mydate2 = format_smartmet_time(mWallClock, +mTimeMargin);
  cpr::Parameters params({{"keyword", *mKeyword},
                          {"param", mParameters},
                          {"format", "json"},
                          {"lang", "fi"},
                          {"starttime", mydate1},
                          {"endtime", mydate2},
                          {"producer", mProducer},
                          {"precision", "full"},
                          {"tz", "UTC"}});

  cpr::Url url{mProtocol + "://" + mHost + "/" + mPlugin};

  read_data(url, params,pConfig);
}

// ----------------------------------------------------------------------
/*!
 * \brief Read observations for certain stations from the server
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::read_ids(const Json::Value& pConfig)
{
  if (mInitialized)
    return;

  std::cout << "Reading data by id " << std::endl;

  // Date conditions
  auto locations = read_filename_locations(pConfig);
  std::string idList;
  std::list<NFmiStation>::iterator it;
  for (it = locations.begin(); it != locations.end(); ++it)
  {
    idList.append(std::to_string(it->GetIdent()));
    idList.append(",");
  }
  idList.pop_back();  // remove last ","
  // idList.pop_back();
  auto mydate1 = format_smartmet_time(mStartTime, -mTimeMargin);
//  auto mydate2 = format_smartmet_time(mEndTime, +mTimeMargin);
  auto mydate2 = format_smartmet_time(mWallClock, +mTimeMargin);
  cpr::Parameters params({{"fmisid", idList},
                          {"param", mParameters},
                          {"format", "json"},
                          {"lang", "fi"},
                          {"starttime", mydate1},
                          {"endtime", mydate2},
                          {"producer", mProducer},
                          {"precision", "full"},
                          {"tz", "UTC"}});

  cpr::Url url{mProtocol + "://" + mHost + "/" + mPlugin};

  read_data(url, params,pConfig);
}
// ----------------------------------------------------------------------
/*!
 * \brief Read observations for a single coordinate from the server
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::read_latlon(const NFmiPoint& pLatLon,const Json::Value& pConfig)
{
  if (mInitialized)
    return;

  std::cout << "Reading " << pLatLon;

  // Date conditions

  auto mydate1 = format_smartmet_time(mStartTime, -mTimeMargin);
  auto mydate2 = format_smartmet_time(mEndTime, +mTimeMargin);

  cpr::Parameters params(
      {{"param", mParameters},
       {"format", "json"},
       {"lang", "fi"},
       {"starttime", mydate1},
       {"endtime", mydate2},
       {"producer", mProducer},
       {"precision", "full"},
       {"tz", "UTC"},
       {"lonlat", Fmi::to_string(pLatLon.X()) + "," + Fmi::to_string(pLatLon.Y())}});

  cpr::Url url{mProtocol + "://" + mHost + "/" + mPlugin};

  read_data(url, params,pConfig);
}

// ----------------------------------------------------------------------
/*!
 * \brief Read data from the server
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::read_data(const cpr::Url& pUrl, const cpr::Parameters& pParams,
		const Json::Value& pConfig)
{
  if (mInitialized)
    throw std::runtime_error("Cannot read SmartMet data twice, logic error?");

  mInitialized = true;

  // FIXME: does not compile
  // std::cout << "Reading data from " << pUrl << "?" << pParams.content << std::endl;

  cpr::Response response;

  {
    boost::timer::auto_cpu_timer timer(2, "\tFetching SmartMet data took %t CPU, %w sec real\n");
    response = cpr::Get(pUrl, pParams);
  }

  if (response.status_code != 200)
  {
    std::string msg = "Failed to retrieve data from SmartMet server. Status code = " +
                      std::to_string(response.status_code) + ".Headers:";
    for (const auto& name_value : response.header)  // NOLINT(cppcheck-useStlAlgorithm)
      msg += "\n\t" + name_value.first + " : " + name_value.second;

    throw std::runtime_error(msg);
  }

  if (response.text.empty())
  {
    std::cerr << "Warning: No observations available from " << pUrl << std::endl;
    return;
  }

  // Parse JSON

  Json::Value json;
  Json::Reader reader;
  bool ok = reader.parse(response.text, json);
  if (!ok)
  {
    std::string msg = reader.getFormattedErrorMessages();
    std::cerr << "Failed to parse SmartMet server response:\n"
              << msg << "\n\nJSON text content:\n"
              << response.text << "\n";

    throw std::runtime_error("Failed to read SmartMet server data");
  }

#if 0
  Json::StyledWriter writer;
  std::cout << "JSON:\n" << writer.write(json) << std::endl;
#endif

  // Store all the data

  int last_id = -1;             // does not exist
  auto last_pos = mData.end();  // not inserting right now

  std::cout << "Processing " << json.size() << " SmartMet rows" << std::endl;

  for (const auto& j : json)
  {
    int id = j["fmisid"].asInt();
    double lon = j["longitude"].asDouble();
    double lat = j["latitude"].asDouble();
    auto time = Fmi::TimeParser::parse_iso(j["time"].asString());

    // Update current observation vector and station number
    if (id != last_id)
    {
      auto pos_bool = mData.insert(std::make_pair(id, Observations()));
      if (!pos_bool.second)
        throw std::runtime_error("Failed to start reading observations for station " +
                                 std::to_string(id));
      last_pos = pos_bool.first;
      last_id = id;
      mStationTree.insert(Station{lon, lat, id});
    }

    Weather data{time};

    Json::Value nulljson;

#if 0
    Json::StyledWriter writer;
    std::cout << "JSON:\n" << writer.write(j) << std::endl;
#endif
    auto jval = j.get(mTroadF, nulljson);
    if (!jval.isNull())
    {
      data.troad = jval.asDouble();
    }
    else
    {
      jval = j.get(mTroadSF, nulljson);
      if (!jval.isNull())
      {
        data.troad = jval.asDouble();
      }
    }
    jval = j.get(mT2mF, nulljson);
    if (!jval.isNull())
      data.t2m = jval.asDouble();

    jval = j.get(mTdewF, nulljson);
    if (!jval.isNull())
      data.tdew2m = jval.asDouble();

    jval = j.get(mRhF, nulljson);
    if (!jval.isNull())
      data.rh2m = jval.asDouble();

    jval = j.get(mWspdF, nulljson);
    if (!jval.isNull())
      data.vz = jval.asDouble();

    jval = j.get(mLradF, nulljson);
    if (!jval.isNull())
      data.lwdn = jval.asDouble();

    jval = j.get(mSradF, nulljson);
    if (!jval.isNull())
      data.swdn = jval.asDouble();

    jval = j.get(mLradNetF, nulljson);
    if (!jval.isNull())
      data.lwnet = jval.asDouble();

    jval = j.get(mSradDirF, nulljson);
    if (!jval.isNull())
      data.swdir = jval.asDouble();
    jval = j.get(mRr1hF, nulljson);
    if (!jval.isNull())
      data.simuprec = jval.asDouble();

    jval = j.get(mRrF, nulljson);
    if (!jval.isNull())
      data.simuprec = jval.asDouble();

    jval = j.get(mRformF, nulljson);
    if (!jval.isNull())
      data.phase = jval.asDouble();

    if (is_missing(data.tdew2m) and !is_missing(data.rh2m) and !is_missing(data.t2m))
      data.tdew2m = CalcTdewOrRH(data.t2m, MISSING, data.rh2m);

    if (is_missing(data.rh2m) and !is_missing(data.tdew2m) and !is_missing(data.t2m))
      data.rh2m = CalcTdewOrRH(data.t2m, data.tdew2m, MISSING);
    last_pos->second.push_back(data);
  }

  auto obs_qc = Json::Path(".model.obs_qc").resolve(pConfig, 0).asInt();
  if (obs_qc > 0)
  {
     auto maxSameTolerance = Json::Path(".model.maxSameTolerance").resolve(pConfig, 0).asInt();
     std::map<int, Observations>::iterator it;
     auto spikeSize = Json::Path(".model.minSpikeSize").resolve(pConfig, 0).asDouble();
     auto maxSpikeDuration = Json::Path(".model.maxSpikeDuration").resolve(pConfig, 0).asInt();
     for (it = mData.begin(); it != mData.end();it++) {
        auto obs=it->second;

        SameValueCheck(obs,maxSameTolerance);
        SpikeCheck(obs,spikeSize,maxSpikeDuration);
        MaxMinCheck(obs);
 	it->second=obs;
     }
  }
//Outlier check needs further studies
/*  if (obs_qc == 2)
  {
     RemoveOutliers();
  }*/


  mStationTree.flush();
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct querydata source for weather
 *
 * Sample source:
 *
 *        {
 *            "name":        "smartmet",
 *            "protocol":    "http",
 *            "host":        "smartmet.fmi.fi",
 *            "plugin",      "timeseries",
 *            "maxdistance": 2.0,
 *            "params":      ["airtemperature","roadtemperature","dewpoint","humidity"],
 *            "timemargin":  0,
 *            "fields":
 *            {
 *                "number":          "fmisid",
 *                "date":            "time",
 *                "airtemperature":  "Temperature",
 *                "roadtemperature": "RoadTemperature",
 *                "dewpoint":        "DewPoint",
 *                "humidity":        "Humidity",
 *            },
 *            "stations",
 *            {
 *                "keyword":   "tiesääasemat",
 *                "producer":  "road",
 *                "number":    "fmisid",
 *                "longitude": "longitude",
 *                "latitude":  "latitude"
 *            }
 *        }
 */
// ----------------------------------------------------------------------

SmartMetSource::Impl::Impl(const NFmiMetTime& pWallClock,
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

  // Parse resource location

  mProtocol = pSource.get("protocol", "http").asString();
  mHost = pSource.get("host", "smartmet.fmi.fi").asString();
  mPlugin = pSource.get("plugin", "timeseries").asString();

  // Maximum distance for observations

  mRadius = pSource.get("maxdistance", 2.0).asDouble();

  if (mRadius > 10)
    throw std::runtime_error("Effective maxdistance for stations > 10km is not permitted");

  // How many minutes to expand the request time interval

  mTimeMargin = pSource.get("timemargin", 0).asInt();

  // Parameters

  auto json = pSource.get("params", nulljson);
  if (json.isNull())
    std::cerr << "Warning: Data source '" + name + "' has no active parameters" << std::endl;
  else
  {
    if (!json.isArray())
      throw std::runtime_error("Parameter setting for source '" + name + "' must be an array");
    for (const auto& p : json)
      mParams.push_back(p.asString());
    check_params();
  }

  // Station settings

  json = pSource.get("stations", nulljson);
  if (json.isNull())
    throw std::runtime_error("Data source '" + name + "' must have a stations setting");
  if (!json.isObject())
    throw std::runtime_error("Data source '" + name + "' stations setting must be a JSON object");

  // Keyword member is set only when specified
  auto j_keyword = json.get("keyword", nulljson);
  if (!j_keyword.isNull())
    mKeyword = j_keyword.asString();

  auto j_fmisid = json.get("fmisid_points", nulljson);
  if (!j_fmisid.isNull())
    mFmisid = j_fmisid.asString();

  mProducer = json.get("producer", "road").asString();

  auto s_number = json.get("number", "fmisid").asString();
  auto s_lon = json.get("longitude", "longitude").asString();
  auto s_lat = json.get("latitude", "latitude").asString();

  // Column names for parameters to be fetched from the server

  json = pSource.get("fields", nulljson);
  if (json.isNull())
    throw std::runtime_error("Data source '" + name + "' must have a fields setting");
  if (!json.isObject())
    throw std::runtime_error("Data source '" + name + "' fields setting must be a JSON object");

  // These are obligatory
  mNumberF = json.get("number", "fmisid").asString();
  mDateF = json.get("date", "time").asString();
  mTroadF = json.get("roadtemperature", "").asString();
  mT2mF = json.get("airtemperature", "").asString();
  mTdewF = json.get("dewpoint", "").asString();
  mRhF = json.get("humidity", "").asString();
  mWspdF = json.get("windspeed", "").asString();
  mLradF = json.get("longwaveradiation", "").asString();
  mSradF = json.get("shortwaveradiation", "").asString();
  mLradNetF = json.get("netlongwaveradiation", "").asString();
  mSradDirF = json.get("directshortwaveradiation", "").asString();
  mRr1hF = json.get("precipitation", "").asString();
  mRrF = json.get("precipitationrate", "").asString();
  mRformF = json.get("precipitationform", "").asString();

  // Check if road surface temperature sensor is specifiedy by using RoadTemperature(:1)
  // syntax, in this case data is obtained from server with label RoadTemperature_#1
  std::string colon = ":";
  size_t col_pos = mTroadF.find(colon);
  if (col_pos != std::string::npos)
  {
    std::string valName = mTroadF.substr(0, col_pos - 1);
    std::string sensorN = mTroadF.substr(col_pos + 1, 1);
    std::string combiner = "_#";
    mTroadSF = valName + combiner + sensorN;
  }

  // order: fmisid,time,lon,lat,troad?,t2m?,tdew?,rh?,lrad?,srad?,rr1h?,rform?

  mParameters = s_number + "," + mDateF + "," + s_lon + "," + s_lat;

  if (!mTroadF.empty())
    mParameters += "," + mTroadF;
  if (!mT2mF.empty())
    mParameters += "," + mT2mF;
  if (!mTdewF.empty())
    mParameters += "," + mTdewF;
  if (!mRhF.empty())
    mParameters += "," + mRhF;
  if (!mWspdF.empty())
    mParameters += "," + mWspdF;
  if (!mLradF.empty())
    mParameters += "," + mLradF;
  if (!mSradF.empty())
    mParameters += "," + mSradF;
  if (!mLradNetF.empty())
    mParameters += "," + mLradNetF;
  if (!mSradDirF.empty())
    mParameters += "," + mSradDirF;
  if (!mRr1hF.empty())
    mParameters += "," + mRr1hF;
  if (!mRrF.empty())
    mParameters += "," + mRrF;
  if (!mRformF.empty())
    mParameters += "," + mRformF;

  // If a keyword is set, read the data now
  //
  std::string tr = "true";
  if (mFmisid == tr)
    read_ids(pConfig);
  else
    read_keyword(pConfig);

  // Add dew point/relative humidity to parameters if only the other is found, as the
  // other is calculated from the other
  if (uses_param("humidity") and uses_param("airtemperature") and !uses_param("dewpoint"))
  {
    std::string dewpoint = "dewpoint";
    mParams.push_back(dewpoint);
  }
  else if (uses_param("dewpoint") and uses_param("airtemperature") and !uses_param("humidity"))
  {
    std::string humidity = "humidity";
    mParams.push_back(humidity);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Check that the required querydata parameters exist
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::check_params()
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
 * \brief Return true if the given parameter is wanted for simulation
 */
// ----------------------------------------------------------------------

bool SmartMetSource::Impl::uses_param(const std::string& pParam) const
{
  return (std::find(mParams.begin(), mParams.end(), pParam) != mParams.end());
}
// ----------------------------------------------------------------------
/*!
 * \brief Intepolate a single value
 *
 * The given index is such that theObservations[thePos].date >= theTime
 */
// ----------------------------------------------------------------------

double SmartMetSource::Impl::interpolate(const Observations& pObservations,
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
//Got Help from ChatGPT
void SmartMetSource::Impl::RemoveTooManySameValues(Observations& pObs,int sameLimit, double Weather::*pPtr){
   int count=1;
   double previous= pObs[0].*pPtr;
   for (std::size_t i = 1; i < pObs.size(); ++i) {
      if (abs(pObs[i].*pPtr-previous)<0.001){
         ++count;
      }else{
         if (count>=sameLimit){
            for (std::size_t j = i-count; j<i;++j){
               pObs[j].*pPtr=MISSING;
             }
         }
	 count=1;
      }
      previous=pObs[i].*pPtr;
   }
   //Check the last values
   if (count>=sameLimit){
      for (std::size_t j =pObs.size()-count; j<pObs.size(); ++j){
          pObs[j].*pPtr=MISSING;
      }

   }

}
void SmartMetSource::Impl::SameValueCheck(Observations& pObs,int sameLimit){
    double Weather::*ptr = nullptr;  // pointer to double member of struct

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

      if (ptr)
        RemoveTooManySameValues(pObs,sameLimit, ptr);
    }
}
void SmartMetSource::Impl::MaxMinCheck(Observations& pObs){
    double Weather::*ptr = nullptr;  // pointer to double member of struct

    for (const auto& p : mParams)
    {
      if (p == "roadtemperature")
        ptr = &Weather::troad;
      else if (p == "airtemperature")
        ptr = &Weather::t2m;
      else if (p == "dewpoint")
        ptr = &Weather::tdew2m;

      if (ptr)
         for (std::size_t i = 0; i < pObs.size(); i++) {
            if (pObs[i].*ptr>50.0 || pObs[i].*ptr<-50.0){
               pObs[i].*ptr=MISSING;
	    }
	 }
    }
}
void SmartMetSource::Impl::RemoveSpikes(Observations& pObs,double spikeSize, std::size_t maxSpikeDuration,
		double Weather::*pPtr){
   double previous= pObs[0].*pPtr;
   for (std::size_t i = 1; i < pObs.size(); i++) {
      if (abs(pObs[i].*pPtr-previous)>spikeSize){
	std::size_t checkLimit=i+maxSpikeDuration;
        auto checkEnd=std::min(pObs.size(),checkLimit);
        for (std::size_t j=i+1; j<checkEnd; j++){
           if (abs(pObs[i].*pPtr-pObs[j].*pPtr)>spikeSize){
             for (std::size_t k =i; k<j; ++k){
                pObs[k].*pPtr=MISSING;
             }
           }
        }
     }
     previous=pObs[i].*pPtr;
   }
}
void SmartMetSource::Impl::SpikeCheck(Observations& pObs,double spikeSize,std::size_t spikeDuration){
    double Weather::*ptr = nullptr;  // pointer to double member of struct

    for (const auto& p : mParams)
    {
      if (p == "roadtemperature")
        ptr = &Weather::troad;
      else if (p == "airtemperature")
        ptr = &Weather::t2m;
      else if (p == "dewpoint")
        ptr = &Weather::tdew2m;

      if (ptr)
        RemoveSpikes(pObs,spikeSize, spikeDuration,ptr);
    }
}
void SmartMetSource::Impl::RemoveOutliers()
{
     std::size_t minSize=50; //minimum datapoin limit
     const auto dt = Fmi::Seconds(10*60);
     auto t = mStartTime.PosixTime();
     auto endTime=mWallClock.PosixTime();
     auto tNext=t+dt;
     while (t<= endTime){
	std::vector<double> values;
	std::vector<int> pos_vector;
	std::vector<int> it_vector;
        std::map<int, Observations>::iterator it;
        for (it = mData.begin(); it != mData.end();it++) {
           auto obs=it->second;
           for (std::size_t i = 0; i < obs.size(); i++) {
               auto pos_time=obs[i].date.PosixTime();
               if (pos_time>=t and pos_time<tNext){
		   if (!is_missing(obs[i].t2m)){
                      it_vector.push_back(it->first);
		      pos_vector.push_back(i);
                      values.push_back(obs[i].t2m);
		   }
               }else if (pos_time>=tNext){
		   break;
	       }
	   }
	}
        if (values.size()>=minSize){
          double mean = std::accumulate(values.begin(),values.end(),0.0) / values.size();
          double sum2=std::inner_product(values.begin(),values.end(),values.begin(),0.0);
          double stdDev=std::sqrt(sum2/values.size()-mean*mean);
          double treshold=3*stdDev;
          std::cout<<mean<<" "<<stdDev<<" "<<treshold<<std::endl;
          for (std::size_t j =0; j<values.size();j++){
             std::cout<<values[j]<<" "<<pos_vector[j]<<" "<<it_vector[j]<<std::endl;
             if(abs(values[j]-mean)>=treshold){
          	  std::cout<<"outlier detected"<<std::endl;
             }
          }
	}
        t+=dt;
	tNext+=dt;

     }

}

// ----------------------------------------------------------------------
/*!
 * \brief Get the weather for the given point
 */
// ----------------------------------------------------------------------

void SmartMetSource::Impl::GetWeather(InputData& pData,
                                      const SimulationTimes& pTimes,
                                      const NFmiPoint& pLatLon)
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return;

//if (!mInitialized)
  //  read_latlon(pLatLon);

  if (!myLatLon || *myLatLon != pLatLon)
  {
    // Update new active location
    myLatLon = pLatLon;
    // No station near enough?
    auto nearest_station =
        mStationTree.nearest(Station{pLatLon.X(), pLatLon.Y()}, Station::ChordLength(mRadius));
    if (!nearest_station)
      return;
    // No data for the station?!?!
    myData = mData.find(nearest_station->ID());
    if (myData == mData.end())
      return;
  }

  const auto& obs = myData->second;
  if (obs.empty())
    return;


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
    if (!is_missing(data.tdew2m))
      pData.tdew[i] = data.tdew2m;
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

std::optional<NFmiMetTime> SmartMetSource::Impl::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                                    const std::string& p) const
{
  double Weather::*ptr = nullptr;
  // Not in this source if parameter is not defined in the configuration
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
  if (!ptr)
    return {};

  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return {};

//if (!mInitialized)
 //   read_latlon(pLatLon);

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
    if (!is_missing(obs[i - 1].*ptr))
      return obs[i - 1].date;
  }
  return {};
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct a SmartMet source for weather
 */
// ----------------------------------------------------------------------

SmartMetSource::SmartMetSource(const NFmiMetTime& pWallClock,
                               const NFmiMetTime& pStartTime,
                               const NFmiMetTime& pEndTime,
                               const Json::Value& pSource,
                               const Json::Value& pConfig)
    : impl(new SmartMetSource::Impl(pWallClock, pStartTime, pEndTime, pSource, pConfig))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Get weather for the given point and time
 */
// ----------------------------------------------------------------------

void SmartMetSource::GetWeather(InputData& pData,
                                const SimulationTimes& pTimes,
                                const NFmiPoint& pLatLon) const

{
  impl->GetWeather(pData, pTimes, pLatLon);
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish latest available road observation time for setting coupling end time
 */
// ----------------------------------------------------------------------

std::optional<NFmiMetTime> SmartMetSource::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                              const std::string& variable) const
{
  return impl->GetLatestObsTime(pLatLon, variable);
}
