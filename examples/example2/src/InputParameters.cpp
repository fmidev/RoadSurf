#include "InputParameters.h"
#include "InputSettings.h"
#include "JsonTools.h"

// ----------------------------------------------------------------------
/*!
 * \brief Construct based on input model settings only
 */
// ----------------------------------------------------------------------

InputParameters::InputParameters(const InputSettings& pSettings)
{
  MinPrecmm = 0.05 * pSettings.DTSecs / 3600.0;
  MinWatmms = 0.01 * pSettings.DTSecs / 3600.0;
  MinSnowmms = 0.1 * pSettings.DTSecs / 3600.0;
  MaxWatmms = MaxPormms + MaxExtmms;
  WDampLim = 0.1 * MaxPormms;
  WWetLim = 0.9 * MaxPormms;
  WWearLim = 0.1 * MaxPormms;
  MinDepmms = 0.01 * pSettings.DTSecs / 3600.0;
  MinIcemms = 0.05 * pSettings.DTSecs / 3600.0;
}

// ----------------------------------------------------------------------
/*!
 * \brief Construct based on input model settings and JSON overrides
 */
// ----------------------------------------------------------------------

InputParameters::InputParameters(const InputSettings& pSettings, const Json::Value& pJson)
    : InputParameters(pSettings)
{
  // Defaults are now set

  if (pJson.isNull())
    return;

  // Override some defaults

  // time dependent variables
  override(&NightOn, pJson, "NightOn");
  override(&NightOff, pJson, "NightOff");
  override(&CalmLimDay, pJson, "CalmLimDay");
  override(&CalmLimNgt, pJson, "CalmLimNgt");
  override(&TrfFricNgt, pJson, "TrfFricNgt");
  override(&TrFfricDay, pJson, "TrFfricDay");

  // physical constants
  override(&Grav, pJson, "Grav");
  override(&SB_Const, pJson, "SB_Const");
  override(&VK_Const, pJson, "VK_Const");
  override(&LVap, pJson, "LVap");
  override(&LFus, pJson, "LFus");
  override(&WatDens, pJson, "WatDens");
  override(&SnowDens, pJson, "SnowDens");
  override(&IceDens, pJson, "IceDens");
  override(&DepDens, pJson, "DepDens");
  override(&WatMHeat, pJson, "WatMHeat");
  override(&PorEvaF, pJson, "PorEvaF");

  // physical properties related to road point
  override(&ZRefW, pJson, "ZRefW");
  override(&ZRefT, pJson, "ZRefT");
  override(&ZeroDisp, pJson, "ZeroDisp");
  override(&ZMom, pJson, "ZMom");
  override(&ZHeat, pJson, "ZHeat");
  override(&Emiss, pJson, "Emiss");
  override(&Albedo, pJson, "Albedo");
  override(&Albedo_surroundings, pJson, "Albedo_Surroundings");
  override(&MaxPormms, pJson, "MaxPormms");
  override(&TClimG, pJson, "TClimG");
  override(&DampDpth, pJson, "DampDpth");
  override(&Omega, pJson, "Omega");
  override(&AZ, pJson, "AZ");
  override(&DampWearF, pJson, "DampWearF");
  override(&AlbDry, pJson, "AlbDry");
  override(&AlbSnow, pJson, "AlbSnow");
  override(&vsh1, pJson, "vsh1");
  override(&vsh2, pJson, "vsh2");
  override(&Poro1, pJson, "Poro1");
  override(&Poro2, pJson, "Poro2");
  override(&RhoB1, pJson, "RhoB1");
  override(&RhoB2, pJson, "RhoB2");
  override(&Silt1, pJson, "Silt1");
  override(&Silt2, pJson, "Silt2");

  // limit values

  override(&freezing_limit_normal, pJson, "freezing_limit_normal");
  override(&snow_melting_limit_normal, pJson, "snow_melting_limit_normal");
  override(&ice_melting_limit_normal, pJson, "ice_melting_limit_normal");
  override(&frost_melting_limit_normal, pJson, "frost_melting_limit_normal");
  override(&frost_formation_limit_normal, pJson, "frost_formation_limit_normal");
  override(&T4Melt_normal, pJson, "T4Melt_normal");

  override(&TLimColdH, pJson, "TLimColdH");
  override(&TLimColdL, pJson, "TLimColdL");
  override(&WetSnowFormR, pJson, "WetSnowFormR");
  override(&WetSnowMeltR, pJson, "WetSnowMeltR");

  override(&PLimSnow, pJson, "PLimSnow");
  override(&PLimRain, pJson, "PLimRain");
  override(&MaxSnowmms, pJson, "MaxSnowmms");
  override(&MaxDepmms, pJson, "MaxDepmms");
  override(&MaxIcemms, pJson, "MaxIcemms");
  override(&MaxExtmms, pJson, "MaxExtmms");

  override(&Snow2IceFac, pJson, "Snow2IceFac");
}
