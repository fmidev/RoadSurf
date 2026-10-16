#include "QueryDataSource.h"
#include "MeteorologyTools.h"
#include <cmath>
#include <smartmet/newbase/NFmiMultiQueryInfo.h>
#include <smartmet/newbase/NFmiParameterName.h>
#include <smartmet/newbase/NFmiStringTools.h>
#include <cmath>
#include <string>
#include <utility>
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
       Json::Value  pSource);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon) const;

  std::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const;

 private:
  using CompareValues = std::vector<long>;  // NFmiTime::GetCompareValue
  using Times = std::vector<NFmiMetTime>;

  [[noreturn]] void error(const std::string& pMessage) const;
  void check_params();
  bool uses_param(FmiParameterName pParam) const;
  void update_timeseries(const NFmiPoint& pLatLon) const;
  void update_pointdata_timeseries(const NFmiPoint& pLatLon) const;
  void update_griddata_timeseries(const NFmiPoint& pLatLon) const;
  void init_times();
  void init_observation_times(Observations& d) const;

  template <typename Getter, typename Converter>
  void store_param(NFmiMultiQueryInfo& q,
                   Observations& d,
                   FmiParameterName pParam,
                   double Weather::*pPtr,
                   Getter get_value,
                   Converter convert) const;

  template <typename Getter>
  void extract_timeseries(NFmiMultiQueryInfo& q, Observations& d, Getter get_value) const;

  std::optional<std::size_t> find_time_position(long tcmp, std::size_t size) const;
  void interpolate_weather(Weather& data,
                           const Observations& d,
                           long tcmp,
                           std::size_t pos) const;

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

namespace
{
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

// Shift the time by the given number of minutes
void apply_time_shift(NFmiMetTime& t, int time_shift)
{
  if (time_shift != 0)
  {
    if (time_shift < 0)
      t.PreviousMetTime(abs(time_shift));
    else
      t.NextMetTime(time_shift);
  }
}

// The member of the weather struct for the parameter, or nullptr if there is none
double Weather::*weather_member(FmiParameterName p)
{
  switch (p)
  {
    case kFmiRoadTemperature:
      return &Weather::troad;
    case kFmiTemperature:
      return &Weather::t2m;
    case kFmiDewPoint:
      return &Weather::tdew2m;
    case kFmiHumidity:
      return &Weather::rh2m;
    case kFmiWindSpeedMS:
      return &Weather::vz;
    case kFmiRadiationLW:
      return &Weather::lwdn;
    case kFmiRadiationGlobal:
      return &Weather::swdn;
    case kFmiRadiationNetSurfaceLW:
      return &Weather::lwnet;
    case kFmiRadiationSW:
      return &Weather::swdir;
    case kFmiPrecipitation1h:
    case kFmiPrecipitationRate:
      return &Weather::simuprec;
    case cf_depth:
      return &Weather::depth;
    default:
      return nullptr;
  }
}

// Store the valid values of the weather into the input data at the given position
void store_weather(InputData& pData, std::size_t i, const Weather& data)
{
  set_if_valid(pData.TSurfObs[i], data.troad);
  set_if_valid(pData.tair[i], data.t2m);
  set_if_valid(pData.Rhz[i], data.rh2m);
  set_if_valid(pData.VZ[i], data.vz);
  set_if_valid(pData.tdew[i], data.tdew2m);
  set_if_valid(pData.LW[i], data.lwdn);
  set_if_valid(pData.SW[i], data.swdn);
  set_if_valid(pData.LW_net[i], data.lwnet);
  set_if_valid(pData.SW_dir[i], data.swdir);
  set_if_valid(pData.prec[i], data.simuprec);
  set_if_valid(pData.PrecPhase[i], data.phase);
  set_if_valid(pData.Depth[i], data.depth);
}
}  // namespace

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
                            Json::Value  pSource)
    : mJson(std::move(pSource)), mWallClock(pWallClock), mStartTime(pStartTime), mEndTime(pEndTime)
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

  init_times();
}

// ----------------------------------------------------------------------
/*!
 * \brief Initialize the times and their compare values for the modeling period
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::init_times()
{
  bool first = true;
  for (mPrivateInfo->ResetTime(); mPrivateInfo->NextTime();)
  {
    auto t = mPrivateInfo->ValidTime();
    apply_time_shift(t, time_shift);
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
      apply_time_shift(t, time_shift);
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

  unsigned int pos2 = 0;
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

  unsigned int pos1 = 0;
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
 * \brief Initialize the observations with the times of the modeling period
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::init_observation_times(Observations& d) const
{
  for (long time_index = mFirstTimeIndex; time_index <= mLastTimeIndex; ++time_index)
  {
    const auto& t = mTimes[time_index - mFirstTimeIndex];
    Weather data{t};
    d.emplace_back(data);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Store the valid values of a parameter into the observations
 */
// ----------------------------------------------------------------------

template <typename Getter, typename Converter>
void QueryDataSource::Impl::store_param(NFmiMultiQueryInfo& q,
                                        Observations& d,
                                        FmiParameterName pParam,
                                        double Weather::*pPtr,
                                        Getter get_value,
                                        Converter convert) const
{
  if (!uses_param(pParam) || !q.Param(pParam))
    return;

  for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
  {
    float value = get_value();
    if (value != kFloatMissing)
      d[idx - mFirstTimeIndex].*pPtr = convert(value);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Extract the parameter values using the given value getter
 */
// ----------------------------------------------------------------------

template <typename Getter>
void QueryDataSource::Impl::extract_timeseries(NFmiMultiQueryInfo& q,
                                               Observations& d,
                                               Getter get_value) const
{
  // Then handle one parameter at a time for speed.

  const auto as_is = [](float value) { return value; };

  store_param(q, d, kFmiRoadTemperature, &Weather::troad, get_value, as_is);
  store_param(q, d, kFmiTemperature, &Weather::t2m, get_value, as_is);
  store_param(q, d, kFmiDewPoint, &Weather::tdew2m, get_value, as_is);
  store_param(q, d, kFmiHumidity, &Weather::rh2m, get_value,
              [](float value) { return std::max(0.0F, std::min(100.0F, value)); });
  store_param(q, d, kFmiWindSpeedMS, &Weather::vz, get_value, as_is);
  store_param(q, d, kFmiRadiationLW, &Weather::lwdn, get_value, as_is);
  store_param(q, d, kFmiRadiationGlobal, &Weather::swdn, get_value, as_is);
  store_param(q, d, kFmiRadiationNetSurfaceLW, &Weather::lwnet, get_value, as_is);
  store_param(q, d, kFmiRadiationSW, &Weather::swdir, get_value, as_is);
  store_param(q, d, kFmiPrecipitation1h, &Weather::simuprec, get_value, as_is);
  store_param(q, d, kFmiPrecipitationRate, &Weather::simuprec, get_value, as_is);
  store_param(q, d, kFmiPrecipitationForm, &Weather::phase, get_value,
              [this](float value)
              {
                if (!precConversion.empty())
                  value = ConvertPrecPhase(precConversion, value);
                return value;
              });
  store_param(q, d, cf_depth, &Weather::depth, get_value, as_is);

  if (uses_param(kFmiTemperature) && q.Param(kFmiTemperature) && uses_param(kFmiHumidity) &&
      q.Param(kFmiHumidity) && !uses_param(kFmiDewPoint))
  {
    // Calculate Tdew frm T and RH
    for (long idx = mFirstTimeIndex; idx <= mLastTimeIndex && q.TimeIndex(idx); ++idx)
    {
      q.Param(kFmiTemperature);
      auto t = get_value();
      q.Param(kFmiHumidity);
      auto rh = get_value();
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
      auto t = get_value();
      q.Param(kFmiDewPoint);
      auto tdew = get_value();
      if (t != kFloatMissing && tdew != kFloatMissing)
        d[idx - mFirstTimeIndex].rh2m = CalcTdewOrRH(t, tdew, MISSING);
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Extract native timestep values from point data
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
  init_observation_times(d);

  extract_timeseries(q, d, [&q]() { return q.FloatValue(); });
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
 * \brief Extract native timestep values from grid data
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
  init_observation_times(d);

  extract_timeseries(q, d, [&q, &loc_cache]() { return q.CachedInterpolation(loc_cache); });
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
 * \brief Find the first time position >= the given time, and remember it
 */
// ----------------------------------------------------------------------

std::optional<std::size_t> QueryDataSource::Impl::find_time_position(long tcmp,
                                                                      std::size_t size) const
{
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
  for (; pos < size; ++pos)
  {
    if (mCompareValues[pos] >= tcmp)
      break;
  }

  if (pos >= size)
    return {};

  timeindex() = pos;
  return pos;
}

// ----------------------------------------------------------------------
/*!
 * \brief Interpolate the weather parameters in time
 */
// ----------------------------------------------------------------------

void QueryDataSource::Impl::interpolate_weather(Weather& data,
                                                const Observations& d,
                                                long tcmp,
                                                std::size_t pos) const
{
  // Must interpolate

  const int max_time_gap = 180;

  for (const auto& p : mParamNames)
  {
    // If there is dewpoint data and temperature data but not relative
    // humidity data, interpolate the RH values calculated in
    // update_timeseries
    if (p == kFmiDewPoint && uses_param(kFmiTemperature) && !uses_param(kFmiHumidity))
      data.rh2m = interpolate(d, tcmp, pos, max_time_gap, &Weather::rh2m);
    else if (p == kFmiHumidity && uses_param(kFmiTemperature) && !uses_param(kFmiDewPoint))
      data.tdew2m = interpolate(d, tcmp, pos, max_time_gap, &Weather::tdew2m);

    if (p == kFmiPrecipitationForm)
    {
      data.phase = nearest(d, tcmp, pos, max_time_gap, &Weather::phase);
    }
    else
    {
      double Weather::*ptr = weather_member(p);  // pointer to double member of struct
      if (ptr)
        data.*ptr = interpolate(d, tcmp, pos, max_time_gap, ptr);
    }
  }
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

    auto position = find_time_position(tcmp, d.size());
    if (!position)
      continue;

    std::size_t pos = *position;

    Weather data{t};  // Note: NFmiMetTime instead of posix_time to avoid conversions

    if (mCompareValues[pos] == tcmp)
      data = d[pos];
    else
      interpolate_weather(data, d, tcmp, pos);

    // clamp rh to 0-100
    if (!std::isnan(data.rh2m))
      data.rh2m = std::max(0.0, std::min(100.0, data.rh2m));

    // sanity check on precipitation amount
    if (data.simuprec > 100)
      data.simuprec = MISSING;

    // TODO: tdew???

    store_weather(pData, i, data);
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
