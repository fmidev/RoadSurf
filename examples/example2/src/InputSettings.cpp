#include "InputSettings.h"
#include "JsonTools.h"
#include <smartmet/macgyver/TimeParser.h>

// ----------------------------------------------------------------------
/*!
 * \brief Establish forecast time
 */
// ----------------------------------------------------------------------

Fmi::DateTime get_forecast_time(const Json::Value &pJson,
                                           const boost::optional<Fmi::DateTime> &pTime)
{
  // Command line option overrides config options
  if (pTime)
    return *pTime;

  Json::Value nulljson;
  auto tmp = Json::Path(".time.now").resolve(pJson);
  if (!tmp.isNull())
    return Fmi::TimeParser::parse(tmp.asString());

  auto now = Fmi::SecondClock::universal_time();

  // Round down to to full minutes

  auto tday = now.time_of_day();
  return {now.date(), Fmi::TimeDuration(tday.hours(), tday.minutes(), 0)};
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish model start time
 */
// ----------------------------------------------------------------------

Fmi::DateTime get_model_start_time(const Fmi::DateTime &pWallClock,
                                              const Json::Value &pJson)
{
  auto hours = Json::Path(".time.analysis").resolve(pJson, 24).asInt();
  return pWallClock - Fmi::Hours(hours);
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish model end time
 */
// ----------------------------------------------------------------------

Fmi::DateTime get_model_end_time(const Fmi::DateTime &pWallClock,
                                            const Json::Value &pJson)
{
  auto hours = Json::Path(".time.forecast").resolve(pJson, 48).asInt();
  return pWallClock + Fmi::Hours(hours);
}

// ----------------------------------------------------------------------
/*!
 * \brief Constuct from JSON settings
 */
// ----------------------------------------------------------------------

InputSettings::InputSettings(const Json::Value &pJson, const Options &pOptions)
    : forecast_time(get_forecast_time(pJson, pOptions.time))
{
  // Establish model start and end times

  start_time = get_model_start_time(forecast_time, pJson);
  end_time = get_model_end_time(forecast_time, pJson);

  // Override some defaults

  const Json::Value nulljson;
  const auto json = pJson.get("model", nulljson);
  if (!json.isNull())
  {
    override(&use_coupling, json, "use_coupling");
    override(&use_relaxation, json, "use_relaxation");
    override(&DTSecs, json, "DTSecs");
    override(&tsurfOutputDepth, json, "tsurfOutputDepth");
    override(&NLayers, json, "NLayers");
    override(&couplingEffectReduction, json, "couplingEffectReduction");
    override(&use_obs_qc, json, "use_obs_qc");
  }
  const auto json2 = pJson.get("output", nulljson);
  if (!json.isNull())
  {
    override(&outputStep, json2, "step");
  }
  // Establish required vector size for simulation.
  const auto duration = end_time - start_time;
  const auto total_secs = duration.total_seconds();

  SimLen = 1 + static_cast<int>(total_secs / DTSecs);

  int coupling_minutes_input = Json::Path(".time.coupling_minutes").resolve(pJson, 0).asInt();
  if (coupling_minutes_input > 0)
    coupling_minutes = coupling_minutes_input;
}
