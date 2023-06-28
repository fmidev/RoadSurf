#include "QueryDataTools.h"
#include "FileTools.h"
#include "InputData.h"
#include "InputSettings.h"
#include "OutputData.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <smartmet/locus/Query.h>
#include <smartmet/newbase/NFmiAreaFactory.h>
#include <smartmet/newbase/NFmiFastQueryInfo.h>
#include <smartmet/newbase/NFmiHPlaceDescriptor.h>
#include <smartmet/newbase/NFmiParamBag.h>
#include <smartmet/newbase/NFmiParamDescriptor.h>
#include <smartmet/newbase/NFmiParameterName.h>
#include <smartmet/newbase/NFmiPoint.h>
#include <smartmet/newbase/NFmiQueryData.h>
#include <smartmet/newbase/NFmiQueryDataUtil.h>
#include <smartmet/newbase/NFmiTimeDescriptor.h>
#include <smartmet/newbase/NFmiVPlaceDescriptor.h>

// ----------------------------------------------------------------------
/*!
 * \brief Read simulation points based on JSON configuration for keywords
 *
 * Sample configuration:
 *
 * "points":
 * {
 *       "keyword":
 *       {
 *             "host": "host.fmi.fi",
 *             "database": "fminames",
 *             "user": "username",
 *             "password": "pass",
 *             "keyword": "roadstations",
 *             "features": "SYNOP",
 *             "countries": "FI"
 *       }
 * }

 */
// ----------------------------------------------------------------------

std::list<NFmiStation> read_keyword_locations(const Json::Value& pJson)
{
  Json::Value nulljson;

  // get points.keyword settings
  auto j = pJson.get("points", nulljson);
  if (!j.isNull())
    j = j.get("keyword", nulljson);
  if (j.isNull())
    return {};

  auto jhost = j.get("host", nulljson);
  auto jdatabase = j.get("database", nulljson);
  auto juser = j.get("user", nulljson);
  auto jpassword = j.get("password", nulljson);
  auto jkeyword = j.get("keyword", nulljson);
  auto jfeatures = j.get("features", "SYNOP");
  auto jcountries = j.get("countries", "");  // default = all

  Locus::QueryOptions options;
  options.SetFeatures(jfeatures.asString());
  options.SetCountries(jcountries.asString());

  Locus::Query query(
      jhost.asString(), juser.asString(), jpassword.asString(), jdatabase.asString());
  auto places = query.FetchByKeyword(options, jkeyword.asString());

  if (places.empty())
    throw std::logic_error("Failed to read any locations from database for keyword '" +
                           jkeyword.asString() + "'");

  std::list<NFmiStation> locations;
  for (const auto& location : places)  // NOLINT(cppcheck-useStlAlgorithm)
    locations.emplace_back(NFmiStation(location.id, location.name, location.lon, location.lat));

  return locations;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read file containing a list of stations
 *
 * The file consists of rows of form "ID,lon,lat,name":
 *
 * Sample:
 *
 * 1001,25.0,65.0,Helsinki Kaisaniemi
 * 1002,25.0,65.1,Helsinki Vantaa
 *
 */
// ----------------------------------------------------------------------
/*!
 * \brief Helper function for creating a parameter
 */
// ----------------------------------------------------------------------

NFmiParam create_param(FmiParameterName pNumber,
                       const std::string& pName,
                       FmiInterpolationMethod pInterpolation,
                       const std::string& pPrecision = "",
                       float pMinValue = kFloatMissing,
                       float pMaxValue = kFloatMissing)
{
  NFmiParam p(pNumber, pName);
  p.InterpolationMethod(pInterpolation);
  if (!pPrecision.empty())
    p.Precision(pPrecision);
  if (pMinValue != kFloatMissing)
    p.MinValue(pMinValue);
  if (pMaxValue != kFloatMissing)
    p.MaxValue(pMaxValue);
  return p;
}

// ----------------------------------------------------------------------
/*!
 * \brief Create querydata parameter descriptor
 */
// ----------------------------------------------------------------------

NFmiParamDescriptor create_param_descriptor(const InputSettings& /* pSettings */,
                                            const Json::Value& /* pJson */)
{
  // Parameters to be written.

  auto p_tsurf = create_param(kFmiRoadTemperature, "Surface temperature", kLinearly);
  auto p_tair = create_param(kFmiTemperature, "Temperature", kLinearly);
  auto p_tdew = create_param(kFmiDewPoint, "Dew point temperature", kLinearly);
  auto p_tdewdiff = create_param(kFmiDewPointDeficit, "Dew point deficit", kLinearly);
  auto p_snow = create_param(kFmiRoadSnowCover, "Snow cover", kLinearly, "%.2f", 0);
  auto p_water = create_param(kFmiWaterStorage, "Water cover", kLinearly, "%.2f", 0);
  auto p_ice = create_param(kFmiRoadIceCover, "Ice cover", kLinearly, "%.2f", 0);
  auto p_deposit = create_param(kFmiRoadFrostCover, "Frost cover", kLinearly, "%.2f", 0);
  auto p_ice2 = create_param(
      kFmiRoadIceCoverOnVehiclePath, "Ice cover on vehicle path", kLinearly, "%.2f", 0);
  // Collect the parameters. Note that the order has to be the same as in write_timeseries
  NFmiParamBag pbag;

  pbag.Add(NFmiDataIdent(p_tsurf));
  pbag.Add(NFmiDataIdent(p_tair));
  pbag.Add(NFmiDataIdent(p_tdew));
  pbag.Add(NFmiDataIdent(p_tdewdiff));
  pbag.Add(NFmiDataIdent(p_snow));
  pbag.Add(NFmiDataIdent(p_water));
  pbag.Add(NFmiDataIdent(p_ice));
  pbag.Add(NFmiDataIdent(p_deposit));
  pbag.Add(NFmiDataIdent(p_ice2));

  return {pbag};
}

// ----------------------------------------------------------------------
/*!
 * \brief Create querydata time descriptor
 */
// ----------------------------------------------------------------------

NFmiTimeDescriptor create_time_descriptor(const InputSettings& pSettings,
                                          const Json::Value& pJson)
{
  const int default_saved_start_time = 0;
  //  const int default_saved_time_step = 60;

  NFmiMetTime time1 = pSettings.forecast_time;
  NFmiMetTime time2 = pSettings.end_time;
  NFmiMetTime origintime = pSettings.forecast_time;

  auto start = Json::Path(".output.start").resolve(pJson, default_saved_start_time).asInt();
  //  auto step = Json::Path(".output.step").resolve(pJson, default_saved_time_step).asInt();

  if (pSettings.outputStep <= 0)
    throw std::logic_error("Output timestep must be > 0 minutes");

  if (pSettings.outputStep * 60 % static_cast<int>(pSettings.DTSecs) != 0)
    throw std::logic_error("Output timestep must be a multiple of simulation timestep");

  // First time to be saved must be >= time1 and a multiple of the chosen timestep

  time1.ChangeByMinutes(start);
  if (time1.IsLessThan(NFmiMetTime(pSettings.start_time)))
    time1 = pSettings.start_time;

  time1.SetTimeStep(pSettings.outputStep, true, kForward);

  // Last time to be saved must be <= time2 and a multiple of the chosen timestep

  time2.SetTimeStep(pSettings.outputStep, true, kBackward);

  NFmiTimeBag tbag(time1, time2, pSettings.outputStep);
  return {origintime, tbag};
}

// ----------------------------------------------------------------------
/*!
 * \brief Create querydata place descriptor
 */
// ----------------------------------------------------------------------

NFmiHPlaceDescriptor create_place_descriptor(const std::list<NFmiStation>& pStations)
{
  NFmiLocationBag lbag;
  for (const auto& station : pStations)
  {
    lbag.AddLocation(station);
  }
  return {lbag};
}

// ----------------------------------------------------------------------
/*!
 * \brief Create empty output querydata
 */
// ----------------------------------------------------------------------

std::unique_ptr<NFmiQueryData> create_coordinates_querydata(const InputSettings& pSettings,
                                                            const Json::Value& pJson)
{
  auto locations = read_keyword_locations(pJson);
  locations.splice(locations.end(), read_filename_locations(pJson));

  auto hdesc = create_place_descriptor(locations);
  auto pdesc = create_param_descriptor(pSettings, pJson);
  auto tdesc = create_time_descriptor(pSettings, pJson);

  NFmiQueryInfo qi(pdesc, tdesc, hdesc, NFmiVPlaceDescriptor());
  return std::unique_ptr<NFmiQueryData>{NFmiQueryDataUtil::CreateEmptyData(qi)};
}

// ----------------------------------------------------------------------
/*!
 * \brief Create grid from json settings
 */
// ----------------------------------------------------------------------

NFmiGrid create_grid(const Json::Value& pJson)
{
  auto jproj = Json::Path(".points.projection").resolve(pJson);
  if (jproj.isNull())
    throw std::logic_error("Simulation projection not set");

  auto area = NFmiAreaFactory::Create(jproj.asString());
  NFmiGrid grid(area.get(),
                static_cast<int>(round(area->Width() + 1)),
                static_cast<int>(round(area->Height() + 1)));
  return grid;
}

// ----------------------------------------------------------------------
/*!
 * \brief Create empty output querydata
 */
// ----------------------------------------------------------------------

std::unique_ptr<NFmiQueryData> create_grid_querydata(const InputSettings& pSettings,
                                                     const Json::Value& pJson)
{
  auto grid = create_grid(pJson);

  auto pdesc = create_param_descriptor(pSettings, pJson);
  auto tdesc = create_time_descriptor(pSettings, pJson);
  auto hdesc = NFmiHPlaceDescriptor(grid);

  NFmiQueryInfo qi(pdesc, tdesc, hdesc, NFmiVPlaceDescriptor());
  return std::unique_ptr<NFmiQueryData>{NFmiQueryDataUtil::CreateEmptyData(qi)};
}

WriteStride get_write_stride(NFmiFastQueryInfo& qi, const InputSettings& pSettings)
{
  qi.FirstTime();
  auto time1 = qi.ValidTime();
  auto diff = time1.DifferenceInMinutes(NFmiMetTime(pSettings.start_time));
  std::size_t startpos = 60 * diff / static_cast<int>(pSettings.DTSecs);

  qi.NextTime();
  auto time2 = qi.ValidTime();
  auto diff2 = time2.DifferenceInMinutes(time1);
  std::size_t stride = 60 * diff2 / static_cast<int>(pSettings.DTSecs);

  return {startpos, stride};
}

std::vector<double> calc_difference(std::vector<double> vec1, std::vector<double> vec2)
{
  std::vector<double> dif;
  for (std::size_t i = 0; i < vec1.size(); i++)
  {
    if (!std::isnan(vec1[i]) and vec1[i] > -9000 and !std::isnan(vec2[i]) and vec2[i] > -9000)
      dif.push_back(vec1[i] - vec2[i]);
    else
      dif.push_back(-9999);
  }
  return dif;
}

template <typename T>
void write_timeseries(NFmiFastQueryInfo& pQ,
                      long pLocIndex,
                      long pParamIndex,
                      const std::vector<T>& pValues,
                      const WriteStride& pStride)
{
  long level_index = 0;

  auto pos = pStride.start;
  const auto step = pStride.step;

  long ntimes = pQ.SizeTimes();

  for (long time_index = 0; time_index < ntimes; ++time_index)
  {
    if (pValues.at(pos) != -9999)
    {
      auto idx = pQ.Index(pParamIndex, pLocIndex, level_index, time_index);
      pQ.IndexFloatValue(idx, pValues[pos]);
    }
    pos += step;
  }
}

void write_timeseries(NFmiFastQueryInfo& pQ,
                      long pLocIndex,
                      const InputData& pInput,
                      const OutputData& pOutput,
                      const WriteStride& pStride)
{
  // Store the values into querydata. Note that the order of parameters
  // has to be the same as when creating the location bag in create_param_descriptor

  long param_index = 0;
  std::vector<double> dewPointDeficit = calc_difference(pOutput.TsurfOut, pInput.tdew);

  write_timeseries(pQ, pLocIndex, param_index++, pOutput.TsurfOut, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pInput.tair, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pInput.tdew, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, dewPointDeficit, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pOutput.SnowOut, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pOutput.WaterOut, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pOutput.IceOut, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pOutput.DepositOut, pStride);
  write_timeseries(pQ, pLocIndex, param_index++, pOutput.Ice2Out, pStride);

  // kFmiRoadConditionSeverity  same as TrafficIndex???
  // kFmiRoadNotification???
}
