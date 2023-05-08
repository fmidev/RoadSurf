#include "InputSettings.h"
#include "JsonTools.h"
#include <iostream>
// ----------------------------------------------------------------------
/*!
 * \brief Establish forecast start time
 */
// ----------------------------------------------------------------------

time_t get_forecast_time(const Json::Value &pJson,
                         const time_t &pTime)
{
  // Command line option overrides config options
  if (pTime)
    return pTime;

  struct tm forecast_start_time;
  Json::Value nulljson;
  //get forecast start time from config file
  auto tmp = Json::Path(".time.now").resolve(pJson).asString();
  if (!tmp.empty())
    {
       strptime(tmp.c_str(),"%Y%m%dT%H%M",&forecast_start_time);
       forecast_start_time.tm_isdst=-1;
       return  mktime(&forecast_start_time);
    }
  //If forecast start time is not given in command line or config file,
  //use curretn time
  time_t now=time(0);
  tm *forecast_start = gmtime(&now);
  // Round down to to full minutes
  forecast_start->tm_sec=0;
  return mktime(forecast_start);

}

// ----------------------------------------------------------------------
/*!
 * \brief Establish model start time
 */
// ----------------------------------------------------------------------

time_t get_model_start_time(const time_t &pWallClock,
                            const Json::Value &pJson)
{
  auto hours = Json::Path(".time.analysis").resolve(pJson, 24).asInt();
  return pWallClock - hours*3600;
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish model end time
 */
// ----------------------------------------------------------------------

time_t get_model_end_time(const time_t &pWallClock,
                          const Json::Value &pJson)
{
  auto hours = Json::Path(".time.forecast").resolve(pJson, 48).asInt();
  return pWallClock + hours*3600;
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
    override(&couplingEffectReduction, json, "couplingEffectReduction");
  }
  const auto json2 = pJson.get("output", nulljson);
  if (!json.isNull())
  {
    override(&outputStep, json2, "step");
  }
  // Establish required vector size for simulation.
  const auto total_secs = end_time - start_time;

  SimLen = 1 + static_cast<int>(total_secs / DTSecs);

  //Lenght of coupling period
  int coupling_minutes_input = Json::Path(".time.coupling_minutes").resolve(pJson, 0).asInt();
  if (coupling_minutes_input > 0)
    coupling_minutes = coupling_minutes_input;
}
