#include "AsciiSource.h"
#include <boost/algorithm/string/predicate.hpp>
#include <json/json.h>
#include <json/reader.h>
#include <json/writer.h>
#include <smartmet/macgyver/StringConversion.h>
#include <smartmet/newbase/NFmiPoint.h>
#include <fstream>

// ----------------------------------------------------------------------
/*!
 * \brief Implementation details
 */
// ----------------------------------------------------------------------

class AsciiSource::Impl
{
 public:
  Impl(const NFmiMetTime& pWallClock,
       const NFmiMetTime& pStartTime,
       const NFmiMetTime& pEndTime,
       const Json::Value& pSource);

  void GetWeather(InputData& pData, const SimulationTimes& pTimes, const NFmiPoint& pLatLon) const;

  boost::optional<NFmiMetTime> GetLatestObsTime(const NFmiPoint& pLatLon,
                                                const std::string& variable) const;

 private:
  const NFmiMetTime mWallClock;
  NFmiMetTime mStartTime;
  NFmiMetTime mEndTime;

  using Observations = std::vector<Weather>;  // observations sorted by time
  Observations mData;

  std::string mFilename;
  NFmiPoint mLatLon;

  void check_params();

  double interpolate(const NFmiMetTime& pTime,
                     int pPos,
                     int pMaxTimeGap,
                     double Weather::*pPtr) const;

  void read_file();

};  // class AsciiSource::Impl

// ----------------------------------------------------------------------
/*!
 * \brief Read data from a file
 */
// ----------------------------------------------------------------------

void AsciiSource::Impl::read_file()
{
  std::cout << "Reading file " << mFilename << std::endl;

  std::ifstream in(mFilename.c_str());
  if (!in)
    throw std::runtime_error("Failed to open '" + mFilename + "' for reading");

  // First line: settings, for example "coupling" and/or "relaxation"

  int rownum = 1;
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("Failed to read line 1 from '" + mFilename + "'");

  using boost::algorithm::contains;

  bool coupling_on = contains(line, "coupling");
  bool relaxation_on = contains(line, "relaxation");

  // Next line is for coupling settings:
  // YY MM DD HH MI SS couplingIndexI couplingIndexJ troad (for coupling)

  if (coupling_on)
  {
    if (!std::getline(in, line))
      throw std::runtime_error("Failed to read coupling line from '" + mFilename + "'");
    ++rownum;
    // ignored
  }

  // Next line is for relaxation:
  // tair_relax vz_relax rh_relax

  if (relaxation_on)
  {
    if (!std::getline(in, line))
      throw std::runtime_error("Failed to read relaxation line from '" + mFilename + "'");
    ++rownum;
    // ignored
  }

  // Next line: simulen initleni initlenj lat lon n/b/t/l

  int simulen = 0;
  int initleni = 0;
  int initlenj = 0;
  double lon = 0;
  double lat = 0;
  std::string type;

  in >> simulen >> initleni >> initlenj >> lon >> lat >> type;

  if (!in)
    throw std::runtime_error("Failed to read simulation settings line from '" + mFilename + "'");

  ++rownum;

  mLatLon = NFmiPoint(lon, lat);

  // Data rows: yy mm dd hh tair rh vz rr1h rform srad lrad tsurf

  int yy;
  int mm;
  int dd;
  int hh;
  int mi;
  int ss;
  int rform;
  double tair;
  double rh;
  double vz;
  double rr1h;
  double srad;
  double lrad;
  double tsurf;

  while (in.good())
  {
    in >> yy >> mm >> dd >> hh >> mi >> ss >> tair >> rh >> vz >> rr1h >> rform >> srad >> lrad >>
        tsurf;

    if (in.bad())
      throw std::runtime_error("Failed to read row " + std::to_string(rownum) + " from '" +
                               mFilename + "'");

    if (in.eof())
      break;

    ++rownum;

    Weather w{NFmiMetTime(yy, mm, dd, hh, mi, ss)};

    if (tsurf > -9999)
      w.troad = tsurf;
    if (tair > -9999)
      w.t2m = tair;
    if (rh > -9999)
      w.rh2m = rh;
    if (vz > -9999)
      w.vz = vz;
    if (rr1h > -9999)
      w.simuprec = rr1h;
    if (rform != -99)
      w.phase = rform;
    if (srad > -9999)
      w.swdn = srad;
    if (lrad > -9999)
      w.lwdn = lrad;

    mData.emplace_back(w);
  }
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

AsciiSource::Impl::Impl(const NFmiMetTime& pWallClock,
                        const NFmiMetTime& pStartTime,
                        const NFmiMetTime& pEndTime,
                        const Json::Value& pSource)
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

  auto jfile = pSource.get("filename", nulljson);
  if (jfile.isNull())
    throw std::runtime_error("filename setting missing for ascii source '" + name + "'");

  mFilename = jfile.asString();

  read_file();
}

// ----------------------------------------------------------------------
/*!
 * \brief Intepolate a single value
 *
 * The given index is such that theObservations[thePos].date >= theTime
 */
// ----------------------------------------------------------------------

double AsciiSource::Impl::interpolate(const NFmiMetTime& pTime,
                                      int pPos,
                                      int pMaxTimeGap,
                                      double Weather::*pPtr) const
{
  if (mData[pPos].date == pTime)
  {
    double value = mData[pPos].*pPtr;
    if (!is_missing(value))
      return value;
  }

  if (pPos == 0)
    return MISSING;

  // Search for first value after current time with valid value

  unsigned int pos2;
  double value2 = MISSING;
  for (pos2 = pPos; pos2 < mData.size(); ++pos2)
  {
    value2 = mData[pos2].*pPtr;
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
    value1 = mData[pos1].*pPtr;
    if (pos1 == 0 || !is_missing(value1))
      break;
  }

  if (is_missing(value1))
    return value1;

  // Do time interpolation if the gap is not too long

  const auto& t1 = mData[pos1].date;
  const auto& t2 = mData[pos2].date;

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
 * \brief Get the weather for the given point in time
 */
// ----------------------------------------------------------------------

void AsciiSource::Impl::GetWeather(InputData& pData,
                                   const SimulationTimes& pTimes,
                                   const NFmiPoint& pLatLon) const
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return;

  if (pLatLon != mLatLon)
    return;

  // Process all times

  for (std::size_t i = 0; i < pTimes.size(); i++)
  {
    const auto& pt = pTimes[i].pt;

    if (pt < mStartTime || pt > mEndTime)
      continue;

    // Not in the time interval?
    if (pt < mData.front().date || pt > mData.back().date)
      continue;

    // Find first position where t >= time. This is guaranteed to succeed
    // after the above tests.

    unsigned int pos;
    for (pos = 0; pos < mData.size(); ++pos)
      if (mData[pos].date >= pt)
        break;

    Weather data{pt};

    const int max_time_gap = 180;

    data.troad = interpolate(pt, pos, max_time_gap, &Weather::troad);
    data.t2m = interpolate(pt, pos, max_time_gap, &Weather::t2m);
    data.rh2m = interpolate(pt, pos, max_time_gap, &Weather::rh2m);
    data.vz = interpolate(pt, pos, max_time_gap, &Weather::vz);
    data.lwdn = interpolate(pt, pos, max_time_gap, &Weather::lwdn);
    data.swdn = interpolate(pt, pos, max_time_gap, &Weather::swdn);
    data.simuprec = interpolate(pt, pos, max_time_gap, &Weather::simuprec);
    data.phase = interpolate(pt, pos, max_time_gap, &Weather::phase);

    // clamp rh to 0-100
    if (!is_missing(data.rh2m))
      data.rh2m = std::max(0.0, std::min(100.0, data.rh2m));

    // sanity check on precipitation amount
    if (data.simuprec > 100)
      data.simuprec = MISSING;

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

boost::optional<NFmiMetTime> AsciiSource::Impl::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                                 const std::string& variable) const
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return {};

  if (pLatLon != mLatLon)
    return {};

  for (std::size_t i = mData.size(); i > 0; i--)
  {
    if (variable == "roadtemperature")
    {
      if (!is_missing(mData[i - 1].troad))
        return mData[i - 1].date;
    }
    else if (variable == "airtemperature")
    {
      if (!is_missing(mData[i - 1].t2m))
        return mData[i - 1].date;
    }
    else if (variable == "humidity")
    {
      if (!is_missing(mData[i - 1].rh2m))
        return mData[i - 1].date;
    }
    else if (variable == "windspeed")
    {
      if (!is_missing(mData[i - 1].vz))
        return mData[i - 1].date;
    }
  }
  return {};
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct a SmartMet source for weather
 */
// ----------------------------------------------------------------------

AsciiSource::AsciiSource(const NFmiMetTime& pWallClock,
                         const NFmiMetTime& pStartTime,
                         const NFmiMetTime& pEndTime,
                         const Json::Value& pSource)
    : impl(new AsciiSource::Impl(pWallClock, pStartTime, pEndTime, pSource))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Get weather for the given point and time
 */
// ----------------------------------------------------------------------

void AsciiSource::GetWeather(InputData& pData,
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

boost::optional<NFmiMetTime> AsciiSource::GetLatestObsTime(const NFmiPoint& pLatLon,
                                                           const std::string& variable) const
{
  return impl->GetLatestObsTime(pLatLon, variable);
}
