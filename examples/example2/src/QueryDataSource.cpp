#include "QueryDataSource.h"
#include "MeteorologyTools.h"
#include <smartmet/newbase/NFmiMultiQueryInfo.h>
#include <smartmet/newbase/NFmiParameterName.h>
#include <smartmet/newbase/NFmiStringTools.h>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace
{
using Observations = std::vector<Weather>;  // observations sorted by time

// Thread local status variables
thread_local std::map<std::string, int> myTimeIndex;
thread_local std::map<std::string, Observations> myData;
thread_local std::map<std::string, std::unique_ptr<NFmiMultiQueryInfo>> myInfo;
}  // namespace

// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

class QueryDataSource::Impl
{
 public:
  Impl(const NFmiMetTime& pWallClock,
       const NFmiMetTime& pStartTime,
       const NFmiMetTime& pEndTime,
       const Json::Value& pSource);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon) const;

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const;

 private:
  using CompareValues = std::vector<long>;  // NFmiTime::GetCompareValue
  using Times = std::vector<NFmiMetTime>;

  void error(const std::string& pMessage) const;
  void check_params();
  bool uses_param(FmiParameterName pParam) const;
  void update_timeseries(const NFmiPoint& pLatLon) const;
  void update_pointdata_timeseries(const NFmiPoint& pLatLon) const;
  void update_griddata_timeseries(const NFmiPoint& pLatLon) const;

  // Access thread local instance of multiqueryinfo
  NFmiMultiQueryInfo& info() const
  {
    auto pos = myInfo.find(mName);
#ifdef MULTIQUERYINFO_CAN_BE_COPIED
    if (pos == myInfo.end())
      myInfo.insert(std::make_pair(
          mName, std::unique_ptr<NFmiMultiQueryInfo>(new NFmiMultiQueryInfo(*mPrivateInfo))));
#else
    // This takes the risk that different threads might see different data in the directory
    if (pos == myInfo.end())
      myInfo.insert(std::make_pair(
          mName, std::unique_ptr<NFmiMultiQueryInfo>(new NFmiMultiQueryInfo(mPath))));
#endif
    return *myInfo.at(mName);
  }

  Observations& data() const
  {
    auto pos = myData.find(mName);
    if (pos != myData.end())
      return pos->second;
    myData.insert(std::make_pair(mName, Observations()));
    return myData.at(mName);
  }

  int& timeindex() const
  {
    auto pos = myTimeIndex.find(mName);
    if (pos != myTimeIndex.end())
      return pos->second;
    myTimeIndex.insert(std::make_pair(mName, 0));
    return myTimeIndex.at(mName);
  }

  double interpolate(const Observations& pObservations,
                     long pCompareValue,
                     int pPos,
                     int pMaxTimeGap,
                     double Weather::*pPtr) const;

  // nearest time version
  double nearest(const Observations& pObservations,
                 long pCompareValue,
                 int pPos,
                 int pMaxTimeGap,
                 double Weather::*pPtr) const;

  const Json::Value mJson;
  long mWallClockCompareValue;

  std::string mName;
  std::string mPath;
  std::vector<std::string> mParams;
  std::vector<FmiParameterName> mParamNames;

  const NFmiMetTime mWallClock;
  NFmiMetTime mStartTime;
  NFmiMetTime mEndTime;
  long mStartTimeCompareValue = 0;
  long mEndTimeCompareValue = 0;

  bool mObservationFlag = false;

  // Cached time comparison values and the respective time index range

  CompareValues mCompareValues;
  Times mTimes;
  long mFirstTimeIndex = -1;
  long mLastTimeIndex = -1;

  // Main thread private info to be used in constructor only or to be copied for thread local use
  std::unique_ptr<NFmiMultiQueryInfo> mPrivateInfo;

  std::string precConversion;
  int time_shift = 0;  // minutes, used to shift time stamp
};                     // class QueryDataSource::Impl


QueryDataSource::~QueryDataSource() = default;

///-----------------------------------------------------------------------
/*!
 * return kfmiparameter of the corresponding varaible
 */
//-----------------------------------------------------------------------
FmiParameterName parse_parameter_name(const std::string& pName)
{
  if (pName == "roadtemperature")
    return kFmiRoadTemperature;
  if (pName == "airtemperature")
    return kFmiTemperature;
  if (pName == "humidity")
    return kFmiHumidity;
  if (pName == "dewpoint")
    return kFmiDewPoint;
  if (pName == "windspeed")
    return kFmiWindSpeedMS;
  if (pName == "depth")
    return cf_depth;
  throw std::runtime_error("Unknown parameter name '" + pName + "'");
}

// ----------------------------------------------------------------------
/*!
 * \brief Error handler
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::error(const std::string& pMessage) const
{
  // name and type are known to be set at this time since DataSourceFactory checks both

  throw std::logic_error("Data source '" + mName + "' of type '" + mJson["type"].asString() + ' ' +
                         pMessage);
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct querydata source for weather
 */
// ----------------------------------------------------------------------

QueryDataSource::Impl::Impl(const NFmiMetTime& pWallClock,
                            const NFmiMetTime& pStartTime,
                            const NFmiMetTime& pEndTime,
                            const Json::Value& pSource)
    : mJson(pSource), mWallClock(pWallClock), mStartTime(pStartTime), mEndTime(pEndTime)
{
  // Read settings
  Json::Value nulljson;
  mName = mJson.get("name", nulljson).asString();
  auto jpath = mJson.get("path", nulljson);
  auto jparams = mJson.get("params", nulljson);
  auto precConv = mJson.get("precConversion", nulljson);
  auto time_shift_param = mJson.get("timeShift", nulljson);

  // Validate settings
  if (jpath.isNull())
    error("has no path setting");

  if (jparams.isNull())
    error("has no params setting");

  if (!jparams.isArray())
    error("params setting must be an array of parameter names");

  if (!precConv.isNull())
    precConversion = precConv.asString();
  if (!time_shift_param.isNull())
    time_shift = time_shift_param.asInt();
  // Store settings
  mPath = jpath.asString();

  for (const auto& jparam : jparams)  // NOLINT(cppcheck-useStlAlgorithm) simpler this way
    mParams.push_back(jparam.asString());

  // Initialize private info
  mPrivateInfo.reset(new NFmiMultiQueryInfo(mPath));

  check_params();

  // The modeling time period

  mWallClockCompareValue = mWallClock.GetCompareValue();
  mStartTimeCompareValue = mStartTime.GetCompareValue();
  mEndTimeCompareValue = mEndTime.GetCompareValue();

  // The data is of observation/analysis type if its last valid time is <= its creation time
  mPrivateInfo->LastTime();
  mObservationFlag = (mPrivateInfo->ValidTime() <= mPrivateInfo->OriginTime());
  mLastTimeIndex = mPrivateInfo->TimeIndex();

  // Initialize time comparevalues

  bool first = true;
  for (mPrivateInfo->ResetTime(); mPrivateInfo->NextTime();)
  {
    auto t = mPrivateInfo->ValidTime();
    if (time_shift != 0)
    {
      if (time_shift < 0)
        t.PreviousMetTime(abs(time_shift));
      else
        t.NextMetTime(time_shift);
    }
    auto tc = t.GetCompareValue();

    if (tc < mStartTimeCompareValue)
      continue;

    // If this is the first time we insert something, wind
    // back to the previous time which should be < theTime

    if (first)
    {
      if (mPrivateInfo->TimeIndex() > 0)
        mPrivateInfo->PreviousTime();
      t = mPrivateInfo->ValidTime();
      if (time_shift != 0)
      {
        if (time_shift < 0)
          t.PreviousMetTime(abs(time_shift));
        else
          t.NextMetTime(time_shift);
      }
      tc = t.GetCompareValue();
      mFirstTimeIndex = mPrivateInfo->TimeIndex();
      first = false;
    }

    mTimes.emplace_back(t);
    mCompareValues.push_back(tc);

    // Abort loop if we got one time beyond the interval
    if (tc > mEndTimeCompareValue)
    {
      mLastTimeIndex = mPrivateInfo->TimeIndex();
      break;
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Check that the required querydata parameters exist
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::check_params()
{
  static std::map<std::string, FmiParameterName> mapping = {
      {"roadtemperature", kFmiRoadTemperature},
      {"airtemperature", kFmiTemperature},
      {"dewpoint", kFmiDewPoint},
      {"humidity", kFmiHumidity},
      {"windspeed", kFmiWindSpeedMS},
      {"longwaveradiation", kFmiRadiationLW},
      {"shortwaveradiation", kFmiRadiationGlobal},
      {"directshortwaveradiation", kFmiRadiationSW},
      {"netlongwaveradiation", kFmiRadiationNetSurfaceLW},
      {"precipitation", kFmiPrecipitation1h},
      {"precipitationrate", kFmiPrecipitationRate},
      {"precipitationform", kFmiPrecipitationForm},
      {"depth", cf_depth}};

  for (const auto& param : mParams)
  {
    auto pos = mapping.find(param);
    if (pos == mapping.end())
      error("parameter '" + param + "' is unknown");

    if (!mPrivateInfo->Param(pos->second))
      error("parameter '" + param + "' is not available in '" + mPath + "'");

    mParamNames.push_back(pos->second);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Return true if the given parameter is wanted for simulation
 */
// ----------------------------------------------------------------------

bool QueryDataSource::Impl::uses_param(FmiParameterName pParam) const
{
  return (std::find(mParamNames.begin(), mParamNames.end(), pParam) != mParamNames.end());
}

// ----------------------------------------------------------------------
/*!
 * \brief Intepolate a single value
 *
 * The given index is such that theObservations[thePos].date >= theTime
 */
// ----------------------------------------------------------------------

double QueryDataSource::Impl::interpolate(const Observations& pObservations,
                                          long pCompareValue,
                                          int pPos,
                                          int pMaxTimeGap,
                                          double Weather::*pPtr) const
{
  if (mCompareValues[pPos] == pCompareValue)
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

  const int gap = mCompareValues[pos2] - mCompareValues[pos1];
  if (gap > pMaxTimeGap)
    return MISSING;

  // Do time interpolation

  const int gap1 = pCompareValue - mCompareValues[pos1];

  double value = ((gap - gap1) * value1 + gap1 * value2) / gap;
  return value;
}

// ----------------------------------------------------------------------
/*!
 * \brief Intepolate a single value using nearest time
 *
 * The given index is such that pObservations[pPos].date >= pTime
 */
// ----------------------------------------------------------------------

double QueryDataSource::Impl::nearest(const Observations& pObservations,
                                      long pCompareValue,
                                      int pPos,
                                      int pMaxTimeGap,
                                      double Weather::*pPtr) const
{
  if (mCompareValues[pPos] == pCompareValue)
  {
    double value = pObservations[pPos].*pPtr;
    if (!is_missing(value))
      return value;
  }

  if (pPos == 0)
    return MISSING;

  // Candidate positions
  unsigned int pos1 = pPos - 1;
  unsigned int pos2 = pPos;

  const int gap1 = pCompareValue - mCompareValues[pos1];
  const int gap2 = mCompareValues[pos2] - pCompareValue;

  if (std::min(gap1, gap2) > pMaxTimeGap)
    return MISSING;

  if (gap1 < gap2)
    return pObservations[pos1].*pPtr;

  return pObservations[pos2].*pPtr;
}

// ----------------------------------------------------------------------
/*!
 * \brief Extract native timestep values from point data
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::update_pointdata_timeseries(const NFmiPoint& pLatLon) const
{
  // Check against missing data
  if (mFirstTimeIndex < 0)
    return;

  // Get ref to thread local queryinfo
  auto& q = info();

  if (!q.NearestPoint(pLatLon))
    return;

  // We want to extract one timestep before and one timestep after
  // the interval we are modeling

  // First extract the times
  auto& d = data();

  for (long time_index = mFirstTimeIndex; time_index <= mLastTimeIndex; ++time_index)
  {
    const auto& t = mTimes[time_index - mFirstTimeIndex];
    Weather data{t};
    d.emplace_back(data);
  }

  // Then handle one parameter at a time for speed. Should be rewritten when C++14 lambdas
  // are available using a lambda which is given a lambda for storing the value.

  float value;

  if (uses_param(kFmiRoadTemperature) && q.Param(kFmiRoadTemperature))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].troad = value;

  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].t2m = value;

  if (uses_param(kFmiDewPoint) && q.Param(kFmiDewPoint))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].tdew2m = value;

  if (uses_param(kFmiHumidity) && q.Param(kFmiHumidity))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].rh2m = std::max(0.0F, std::min(100.0F, value));

  if (uses_param(kFmiWindSpeedMS) && q.Param(kFmiWindSpeedMS))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].vz = value;

  if (uses_param(kFmiRadiationLW) && q.Param(kFmiRadiationLW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].lwdn = value;

  if (uses_param(kFmiRadiationGlobal) && q.Param(kFmiRadiationGlobal))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].swdn = value;

  if (uses_param(kFmiRadiationNetSurfaceLW) && q.Param(kFmiRadiationNetSurfaceLW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].lwnet = value;

  if (uses_param(kFmiRadiationSW) && q.Param(kFmiRadiationSW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].swdir = value;

  if (uses_param(kFmiPrecipitation1h) && q.Param(kFmiPrecipitation1h))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].simuprec = value;

  if (uses_param(kFmiPrecipitationRate) && q.Param(kFmiPrecipitationRate))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].simuprec = value;

  if (uses_param(kFmiPrecipitationForm) && q.Param(kFmiPrecipitationForm))
  {
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      if ((value = q.FloatValue()) != kFloatMissing)
      {
        if (!precConversion.empty())
          value = ConvertPrecPhase(precConversion, value);
        d[idx - mFirstTimeIndex].phase = value;
      }
    }
  }
  if (uses_param(cf_depth) && q.Param(cf_depth))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.FloatValue()) != kFloatMissing)
        d[idx - mFirstTimeIndex].depth = value;

  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature) && uses_param(kFmiHumidity) &&
      q.Param(kFmiHumidity) && !uses_param(kFmiDewPoint))
  {
    // Calculate Tdew frm T and RH
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      q.Param(kFmiTemperature);
      auto t = q.FloatValue();
      q.Param(kFmiHumidity);
      auto rh = q.FloatValue();
      if (t != kFloatMissing && rh != kFloatMissing)
        d[idx - mFirstTimeIndex].tdew2m = CalcTdewOrRH(t, MISSING, rh);
    }
  }
  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature) && uses_param(kFmiDewPoint) &&
      q.Param(kFmiDewPoint) && !uses_param(kFmiHumidity))
  {
    // Calculate RH frm T and TDEW
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      q.Param(kFmiTemperature);
      auto t = q.FloatValue();
      q.Param(kFmiDewPoint);
      auto tdew = q.FloatValue();
      if (t != kFloatMissing && tdew != kFloatMissing)
        d[idx - mFirstTimeIndex].rh2m = CalcTdewOrRH(t, tdew, MISSING);
    }
  }
 // std::ofstream outputFile("/smartmet/data/virve/txt/atm_forecast.txt",std::ios::app);
 // if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature)){
 //    outputFile<<pLatLon.Y()<<" "<<pLatLon.X()<<std::endl;
 //    outputFile<<"date air_temperature humidity wind_speed prec_phase"<<std::endl;
 //    for (std::size_t i = 0; i < d.size(); i++){
 //       outputFile<<std::fixed<<std::setprecision(2)<<d[i].date<<" "<<d[i].t2m<<" "<<d[i].rh2m<<" "<<d[i].vz<<" "<<d[i].phase<<std::endl;
 //    }
 // }

 // if (uses_param(kFmiRadiationLW) && q.Param(kFmiRadiationLW)){
 //    outputFile<<"date SW LW prec"<<std::endl;
 //    for (std::size_t i = 0; i < d.size(); i++){
 //       outputFile<<std::fixed<<std::setprecision(2)<<d[i].date<<" "<<d[i].swdn<<" "<<d[i].lwdn<<" ";
 //       outputFile<<std::fixed<<std::setprecision(4)<<d[i].simuprec<<std::endl;
 //    }
 // }
 // outputFile.close();
}
// ----------------------------------------------------------------------
/*!
 * \brief Extract native timestep values from grid data
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::update_griddata_timeseries(const NFmiPoint& pLatLon) const
{
  // Check against missing data
  if (mFirstTimeIndex < 0)
    return;

  // Get ref to thread local queryinfo
  auto& q = info();

  auto width = q.Grid()->XNumber();
  auto height = q.Grid()->YNumber();
  NFmiLocationCache loc_cache = q.CalcLocationCache(pLatLon, width, height);

  // First extract the times

  auto& d = data();

  for (long time_index = mFirstTimeIndex; time_index <= mLastTimeIndex; ++time_index)
  {
    const auto& t = mTimes[time_index - mFirstTimeIndex];
    Weather data{t};
    d.emplace_back(data);
  }

  // Then handle one parameter at a time for speed. Should be rewritten when C++14 lambdas
  // are available using a lambda which is given a lambda for storing the value.

  float value;

  if (uses_param(kFmiRoadTemperature) && q.Param(kFmiRoadTemperature))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].troad = value;

  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].t2m = value;

  if (uses_param(kFmiDewPoint) && q.Param(kFmiDewPoint))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].tdew2m = value;

  if (uses_param(kFmiHumidity) && q.Param(kFmiHumidity))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].rh2m = std::max(0.0F, std::min(100.0F, value));

  if (uses_param(kFmiWindSpeedMS) && q.Param(kFmiWindSpeedMS))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].vz = value;

  if (uses_param(kFmiRadiationLW) && q.Param(kFmiRadiationLW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].lwdn = value;

  if (uses_param(kFmiRadiationGlobal) && q.Param(kFmiRadiationGlobal))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].swdn = value;

  if (uses_param(kFmiRadiationNetSurfaceLW) && q.Param(kFmiRadiationNetSurfaceLW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].lwnet = value;

  if (uses_param(kFmiRadiationSW) && q.Param(kFmiRadiationSW))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].swdir = value;

  if (uses_param(kFmiPrecipitation1h) && q.Param(kFmiPrecipitation1h))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].simuprec = value;

  if (uses_param(kFmiPrecipitationRate) && q.Param(kFmiPrecipitationRate))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].simuprec = value;

  if (uses_param(kFmiPrecipitationForm) && q.Param(kFmiPrecipitationForm))
  {
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
      {
        if (!precConversion.empty())
          value = ConvertPrecPhase(precConversion, value);
        d[idx - mFirstTimeIndex].phase = value;
      }
    }
  }
  if (uses_param(cf_depth) && q.Param(cf_depth))
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
      if ((value = q.CachedInterpolation(loc_cache)) != kFloatMissing)
        d[idx - mFirstTimeIndex].depth = value;

  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature) && uses_param(kFmiHumidity) &&
      q.Param(kFmiHumidity) && !uses_param(kFmiDewPoint))
  {
    // Calculate Tdew frm T and RH
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      q.Param(kFmiTemperature);
      auto t = q.FloatValue();
      q.Param(kFmiHumidity);
      auto rh = q.FloatValue();
      if (t != kFloatMissing && rh != kFloatMissing)
        d[idx - mFirstTimeIndex].tdew2m = CalcTdewOrRH(t, MISSING, rh);
    }
  }
  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature) && uses_param(kFmiDewPoint) &&
      q.Param(kFmiDewPoint) && !uses_param(kFmiHumidity))
  {
    // Calculate RH frm T and TDEW
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      q.Param(kFmiTemperature);
      auto t = q.CachedInterpolation(loc_cache);
      q.Param(kFmiDewPoint);
      auto tdew = q.CachedInterpolation(loc_cache);
      if (t != kFloatMissing && tdew != kFloatMissing)
        d[idx - mFirstTimeIndex].rh2m = CalcTdewOrRH(t, tdew, MISSING);
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Extract native timestep values for the given coordinate
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::update_timeseries(const NFmiPoint& pLatLon) const
{
  data().clear();

  if (!info().IsGrid())
    update_pointdata_timeseries(pLatLon);
  else
    update_griddata_timeseries(pLatLon);
}

// ----------------------------------------------------------------------
/*!
 * \brief Get the weather for the given point in time
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::GetWeather(InputData& pData,
                                       const SimulationTimes& pTimes,
                                       const NFmiPoint& pLatLon) const
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return;

  // Check against missing data
  if (mFirstTimeIndex < 0)
    return;

  update_timeseries(pLatLon);

  auto& d = data();

  if (d.empty())
    return;

  // Process all times

  for (std::size_t i = 0; i < pTimes.size(); i++)
  {
    // const auto& pt = pTimes[i].pt;
    const auto& t = pTimes[i].t;

    auto tcmp = t.GetCompareValue();

    // No data outside the time range
    if (tcmp < mStartTimeCompareValue || tcmp > mEndTimeCompareValue)
      continue;

    // Do not provide observations if simulation time would be exceeded.
    if (mObservationFlag && tcmp > mWallClockCompareValue)
      continue;

    // Search the first time point >= theTime. We are likely to only advance
    // in time, hence we are not likely to take any steps back in this
    // first loop.

    std::size_t pos = timeindex();

    while (pos > 0)
    {
      if (mCompareValues[pos - 1] < tcmp)
        break;
      --pos;
    }

    // This is also an unlikely path. In the first iteration it may perform
    // several loops, after that 0 is much more likely than one.
    for (; pos < d.size(); ++pos)
    {
      if (mCompareValues[pos] >= tcmp)
        break;
    }

    if (pos >= d.size())
      continue;

    timeindex() = pos;

    Weather data{t};  // Note: NFmiMetTime instead of posix_time to avoid conversions

    if (mCompareValues[pos] == tcmp)
      data = d[pos];
    else
    {
      // Must interpolate

      const int max_time_gap = 180;

      double Weather::*ptr = nullptr;  // pointer to double member of struct

      for (const auto& p : mParamNames)
      {
        // If there is dewpoint data and temperature data but not relative
        // humidity data, interpolate the RH values calculated in
        // update_timeseries
        if (p == kFmiDewPoint && uses_param(kFmiTemperature) && !uses_param(kFmiHumidity))
        {
          ptr = &Weather::rh2m;
          data.*ptr = interpolate(d, tcmp, pos, max_time_gap, ptr);
        }
        else if (p == kFmiHumidity && uses_param(kFmiTemperature) && !uses_param(kFmiDewPoint))
        {
          ptr = &Weather::tdew2m;
          data.*ptr = interpolate(d, tcmp, pos, max_time_gap, ptr);
        }
        if (p == kFmiPrecipitationForm)
        {
          ptr = &Weather::phase;
          data.*ptr = nearest(d, tcmp, pos, max_time_gap, ptr);
        }
        else
        {
          if (p == kFmiRoadTemperature)
            ptr = &Weather::troad;
          else if (p == kFmiTemperature)
            ptr = &Weather::t2m;
          else if (p == kFmiDewPoint)
            ptr = &Weather::tdew2m;
          else if (p == kFmiHumidity)
            ptr = &Weather::rh2m;
          else if (p == kFmiWindSpeedMS)
            ptr = &Weather::vz;
          else if (p == kFmiRadiationLW)
            ptr = &Weather::lwdn;
          else if (p == kFmiRadiationGlobal)
            ptr = &Weather::swdn;
          else if (p == kFmiRadiationNetSurfaceLW)
            ptr = &Weather::lwnet;
          else if (p == kFmiRadiationSW)
            ptr = &Weather::swdir;
          else if (p == kFmiPrecipitation1h || p == kFmiPrecipitationRate)
            ptr = &Weather::simuprec;
          else if (p == cf_depth)
            ptr = &Weather::depth;

          if (ptr)
            data.*ptr = interpolate(d, tcmp, pos, max_time_gap, ptr);
        }
      }
    }

    // clamp rh to 0-100
    if (!std::isnan(data.rh2m))
      data.rh2m = std::max(0.0, std::min(100.0, data.rh2m));

    // sanity check on precipitation amount
    if (data.simuprec > 100)
      data.simuprec = MISSING;

    // TODO: tdew???

    if (!is_missing(data.troad))
      pData.TSurfObs[i] = data.troad;
    if (!is_missing(data.t2m))
      pData.tair[i] = data.t2m;
    if (!is_missing(data.rh2m))
      pData.Rhz[i] = data.rh2m;
    if (!is_missing(data.vz))
      pData.VZ[i] = data.vz;
    if (!is_missing(data.tdew2m))
      pData.tdew[i] = data.tdew2m;
    if (!is_missing(data.lwdn))
      pData.LW[i] = data.lwdn;
    if (!is_missing(data.swdn))
      pData.SW[i] = data.swdn;
    if (!is_missing(data.lwnet))
      pData.LW_net[i] = data.lwnet;
    if (!is_missing(data.swdir))
      pData.SW_dir[i] = data.swdir;
    if (!is_missing(data.simuprec))
      pData.prec[i] = data.simuprec;
    if (!is_missing(data.phase))
      pData.PrecPhase[i] = data.phase;
    if (!is_missing(data.depth))
      pData.Depth[i] = data.depth;
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish latest available road observation time for setting coupling end time
 */
// ----------------------------------------------------------------------

std::optional<NFmiMetTime> QueryDataSource::Impl::GetLatestObsTime(
    const NFmiPoint& pLatLon, const std::string& variable) const
{
  // Check against missing data
  if (mFirstTimeIndex < 0)
    return {};

  FmiParameterName parameterName = parse_parameter_name(variable);
  // Not in this source if parameter is not defined in the configuration
  auto pos = find(mParamNames.begin(), mParamNames.end(), parameterName);
  if (pos == mParamNames.end())
    return {};

  // Get ref to thread local queryinfo
  auto& q = info();

  q.Param(parameterName);  // check_params guarantees success

  std::optional<NFmiMetTime> ret;

  for (long time_index = mFirstTimeIndex; time_index <= mLastTimeIndex; ++time_index)
  {
    double tmp = q.InterpolatedValue(pLatLon);
    if (!is_missing(tmp))
      ret = mTimes[time_index - mFirstTimeIndex];
  }

  return ret;
}

/// ----------------------------------------------------------------------
/*!
 * \brief Construct a querydata source for weather
 */
// ----------------------------------------------------------------------

QueryDataSource::QueryDataSource(const NFmiMetTime& pWallClock,
                                 const NFmiMetTime& pStartTime,
                                 const NFmiMetTime& pEndTime,
                                 const Json::Value& pSource)
    : impl(new QueryDataSource::Impl(pWallClock, pStartTime, pEndTime, pSource))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Get weather for the given point and time
 */
// ----------------------------------------------------------------------

void QueryDataSource::GetWeather(InputData& pData,
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

std::optional<NFmiMetTime> QueryDataSource::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                               const std::string& variable) const
{
  return impl->GetLatestObsTime(pLatLon, variable);
}
