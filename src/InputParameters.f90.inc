!MIT License
!Copyright (c) 2023 FMI Open Development
!> input parameters given by modelRunner.cpp
Type, Bind(C) :: InputParameters

      REAL(C_DOUBLE) :: NightOn         !< Beginning hour of night traffic (UTC)
      REAL(C_DOUBLE) :: NightOff        !< Ending hour of night traffic (UTC)
      REAL(C_DOUBLE) :: CalmLimDay      !< Minimum wind speed day (m/s)
      REAL(C_DOUBLE) :: CalmLimNgt      !< Minimum wind speed night (m/s)
      REAL(C_DOUBLE) :: TrfFricNgt      !< Traffic induced friction (night) (W/m2)
      REAL(C_DOUBLE) :: TrFfricDay      !< Traffic induced friction day (W/m2)

      REAL(C_DOUBLE) :: Grav            !< Gravitational acceleration (m/s2)
      REAL(C_DOUBLE) :: SB_Const        !< Stefan-boltzman constant (W/m2K4)
      REAL(C_DOUBLE) :: VK_Const        !< Von Karman's konstant
      REAL(C_DOUBLE) :: LVap            !< Latent heat of water vaporisation (J/kg)
      REAL(C_DOUBLE) :: LFus            !< Latent heat of fusion (constant,not calculated)
      REAL(C_DOUBLE) :: WatDens         !< Density of water at 0 C
      REAL(C_DOUBLE) :: SnowDens        !< Density of snow (fresh,Oke p.44)
      REAL(C_DOUBLE) :: IceDens         !< Density of ice (Oke p.44)
      REAL(C_DOUBLE) :: DepDens         !< Density of deposit
      REAL(C_DOUBLE) :: WatMHeat        !< Heat of ablation for water (J/kg)
      REAL(C_DOUBLE) :: PorEvaF         !< Pore resistance factor for evaporation

      REAL(C_DOUBLE) :: ZRefW           !< Wind reference height (m)
      REAL(C_DOUBLE) :: ZRefT           !< Wind reference height (m)
      REAL(C_DOUBLE) :: ZeroDisp        !< Zero displacement height (m)
      REAL(C_DOUBLE) :: ZMom            !< Roughness factor for momentum (m)
      REAL(C_DOUBLE) :: ZHeat           !< Roughness factor for heat (m)
      REAL(C_DOUBLE) :: Emiss           !< Emissivity constant of the surface
      REAL(C_DOUBLE) :: Albedo          !< Dry ground albedo
      REAL(C_DOUBLE) :: Albedo_surroundings !< albedo of surroundings
      REAL(C_DOUBLE) :: MaxPormms       !< maximum porosity
      REAL(C_DOUBLE) :: TClimG          !< Climatological temperature at the bottom layer
      REAL(C_DOUBLE) :: DampDpth        !< Damping depth
      REAL(C_DOUBLE) :: Omega           !<constant to calculate bottom layer temperature
      REAL(C_DOUBLE) :: AZ              !<constant to calculate bottom layer temperature
      REAL(C_DOUBLE) :: DampWearF       !< Damp surface (poer) wear factor
      REAL(C_DOUBLE) :: AlbDry          !< Asphalt albedo
      REAL(C_DOUBLE) :: AlbSnow         !< Snow albedo
      REAL(C_DOUBLE) :: vsh1            !< Heat capacity of surface layers (dry)
      REAL(C_DOUBLE) :: vsh2            !< Heat capacity of deep ground layers (dry)
      REAL(C_DOUBLE) :: Poro1           !< Porosity for surface layers
      REAL(C_DOUBLE) :: Poro2           !< Porosoty for deep ground layers
      REAL(C_DOUBLE) :: RhoB1           !< Bulk density for surface layers
      REAL(C_DOUBLE) :: RhoB2           !< Bulk density for ground layers
      REAL(C_DOUBLE) :: Silt1           !< Clay fraction for surface layers
      REAL(C_DOUBLE) :: Silt2           !< Clay fraction for deep ground layers
       
      REAL(C_DOUBLE) :: freezing_limit_normal         !< Freezing limit (C)
      REAL(C_DOUBLE) :: snow_melting_limit_normal     !< Melting limit for snow (C)
      REAL(C_DOUBLE) :: ice_melting_limit_normal      !< Melting limit for ice (C)
      REAL(C_DOUBLE) :: frost_melting_limit_normal    !< Melting limit for deposit (C)
      REAL(C_DOUBLE) :: frost_formation_limit_normal  !< Dew/deposit formation limit (C)
      REAL(C_DOUBLE) :: T4Melt_normal                 !< Normal melting limit (C)
                                                      !< surface class interpretation

      REAL(C_DOUBLE) :: TLimColdH;      !< Higher limit for cold ground temperature
      REAL(C_DOUBLE) :: TLimColdL;      !< Lower limit for cold ground temperature
      REAL(C_DOUBLE) :: WetSnowFormR;   !< Water to snow ratio for wet snow formation
      REAL(C_DOUBLE) :: WetSnowMeltR;   !< Water to snow ratio for wet snow melting

      REAL(C_DOUBLE) :: PLimSnow;       !< Precipitation interpretation : snow limit
      REAL(C_DOUBLE) :: PLimRain;       !< Precipitation interpretation : rain limit
      REAL(C_DOUBLE) :: MaxSnowmms;     !< MAX snow storage : mm (plowed away if above)
      REAL(C_DOUBLE) :: MaxDepmms;      !< MAX deposit storage : mm
      REAL(C_DOUBLE) :: MaxIcemms;      !< MAX ice storage : mm
      REAL(C_DOUBLE) :: MaxExtmms;      !< MAX extra (surface)
                                        !< water content (mm)

      REAL(C_DOUBLE) :: MissValI;       !< Missing value, integer
      REAL(C_DOUBLE) :: MissValR;       !< Missing value, real
                                                                                                          
      REAL(C_DOUBLE) :: Snow2IceFac;  !< Snow to ice transition factor

      REAL(C_DOUBLE) :: MinPrecmm;      !< MIN precipitation  mm/hour
      REAL(C_DOUBLE) :: MinWatmms;      !< MIN water storage : mm
      REAL(C_DOUBLE) :: MinSnowmms;     !< MIN snow storage : mm
                                        !< (accounts for snow "drift")
      REAL(C_DOUBLE) :: MaxWatmms;
      REAL(C_DOUBLE) :: WDampLim;       !< Dry/Damp limit for water
      REAL(C_DOUBLE) :: WWetLim;        !< Damp/Wet limit for water 

      REAL(C_DOUBLE) :: WWearLim;       !< Wear limit for water (below this only
                                        !< evaporation)

      REAL(C_DOUBLE) :: MinDepmms;      !< MIN deposit storage : mm
      REAL(C_DOUBLE) :: MinIcemms;      !< MIN ice storage : mm (MinPrec => condens
                                        !< won't freeze(?))

end type InputParameters
