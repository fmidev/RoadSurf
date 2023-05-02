#pragma once
#include "Constants.h"
#include <json/json.h>
#include <cmath>

struct InputModelSettings;

//!> input parameters given by modelRunner.cpp

struct InputParameters
{
  InputParameters() = delete;

  InputParameters(const InputModelSettings& pSettings, const Json::Value& pJson);
  explicit InputParameters(const InputModelSettings& pSettings);

  // time dependent variables
  double NightOn = 19.0;     //!< Beginning hour of night traffic (UTC)
  double NightOff = 4.0;     //!< Ending hour of night traffic (UTC)
  double CalmLimDay = 1.5;   //!< Minimum wind speed day (m/s)
  double CalmLimNgt = 0.4;   //!< Minimum wind speed night (m/s)
  double TrfFricNgt = 5.0;   //!< Traffic induced friction night (W/m2)
  double TrFfricDay = 10.0;  //!< Traffic induced friction day (W/m2)

  // physical constants
  double Grav = 9.81;          //!< Gravitational acceleration (m/s2)
  double SB_Const = 5.67e-8;   //!< Stefan-boltzman constant (W/m2K4)
  double VK_Const = 0.4;       //!< Von Karman's konstant
  double LVap = 2.452E6;       //!< Latent heat of water vaporisation (J/kg)
  double LFus = 0.334E6;       //!< Latent heat of fusion (constant,not calculated)
  double WatDens = 999.87;     //!< Density of water at 0 C
  double SnowDens = 100.0;     //!< Density of snow (fresh,Oke p.44)
  double IceDens = 920.0;      //!< Density of ice (Oke p.44)
  double DepDens = 920.0;      //!< Density of deposit
  double WatMHeat = 333000.0;  //!< Heat of ablation for water (J/kg)
  double PorEvaF = 1.0;        //!< Pore resistance factor for evaporation

  // physical properties related to road point
  double ZRefW = 10.0;                //!< Wind reference height (m)
  double ZRefT = 2.0;                 //!< Wind reference height (m)
  double ZeroDisp = 0.0;              //!< Zero displacement height (m)
  double ZMom = 0.4000;               //!< Roughness factor for momentum (m)
  double ZHeat = 0.0010;              //!< Roughness factor for heat (m)
  double Emiss = 0.95;                //!< Emissivity constant of the surface
  double Albedo = 0.10;               //!< Dry ground albedo
  double Albedo_surroundings = 0.15;  //!< albedo of surroundings
  double MaxPormms = 1.0;             //!< maximum porosity
  double TClimG = 6.4;                //!< Climatological temperature at the bottom layer
  double DampDpth = 2.7;              //!< Damping depth
  double Omega = 2.0 * M_PI / 365.0;  //!< constant to calculate bottom layer temperature
  double AZ = 0.6;                    //!< constant to calculate bottom layer temperature
  double DampWearF = 0.5;             //!< Damp surface (poer) wear factor
  double AlbDry = 0.1;                //!> Asphalth alpedo
  double AlbSnow = 0.6;               //!< Snow albedo
  double vsh1 = 1.94e+06;             //!< Heat capacity of surface layers (dry)
  double vsh2 = 1.28e+06;             //!< Heat capacity of deep ground layers (dry)
  double Poro1 = 0.1;                 //!< Porosity for surface layers
  double Poro2 = 0.4;                 //!< Porosoty for deep ground layers
  double RhoB1 = 2.11;                //!< Bulk density for surface layers
  double RhoB2 = 1.6;                 //!< Bulk density for ground layers
  double Silt1 = 0.1;                 //!< Clay fraction for surface layers
  double Silt2 = 0.8;                 //!< Clay fraction for deep ground layers

  // limit values

  double freezing_limit_normal = -0.25;        //!< Freezing limit (C)
  double snow_melting_limit_normal = 0.25;     //!< Melting limit for snow (C)
  double ice_melting_limit_normal = 0.25;      //!< Melting limit for ice (C)
  double frost_melting_limit_normal = 1.25;    //!< Melting limit for deposit (C)
  double frost_formation_limit_normal = 0.25;  //!< Dew/deposit formation limit (C)
  double T4Melt_normal = 0.25;                 //!< Normal melting limit (C)

  double TLimColdH = -19.0;   //!< Higher limit for cold ground temperature
  double TLimColdL = -21.0;   //!< Lower limit for cold ground temperature
  double WetSnowFormR = 0.1;  //!< Water to snow ratio for wet snow formation
  double WetSnowMeltR = 0.6;  //!< Water to snow ratio for wet snow melting

  double PLimSnow = 0.3;      //!< Precipitation interpretation : snow limit
  double PLimRain = 0.7;      //!< Precipitation interpretation : rain limit
  double MaxSnowmms = 100.0;  //!< MAX snow storage : mm (plowed away if above)
  double MaxDepmms = 2.0;     //!< MAX deposit storage : mm
  double MaxIcemms = 50.0;    //!< MAX ice storage : mm
  double MaxExtmms = 1.0;     //!< MAX extra (surface) water content (mm)


  // missing values
  double MissValI = -9999;
  double MissValR = -99.99;

  // Traffic index : precipitation and wind speed classification limits
  // - numerical values for precipitation as per hour, converted to per time step


  double Snow2IceFac = 0.5;  //!< Snow to ice transition factor
                                 //!< less than this limit
  // Derived values

  double MinPrecmm;   //!< MIN precipitation  mm/hour
  double MinWatmms;   //!< MIN water storage : mm
  double MinSnowmms;  //!< MIN snow storage : mm (accounts for snow "drift")

  double MaxWatmms;

  double WDampLim;  //!< Dry/Damp limit for water
  double WWetLim;   //!< Damp/Wet limit for water

  double WWearLim;  //!< Wear limit for water (below this only evaporation)

  double MinDepmms;  //!< MIN deposit storage : mm
  double MinIcemms;  //!< MIN ice storage : mm (MinPrec => condens won't freeze(?))

};
