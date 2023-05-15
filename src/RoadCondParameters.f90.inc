!MIT License
!Copyright (c) 2023 FMI Open Development
!> Parameters to determine storage terms and road condition
type RoadCondParameters
      Real(8) :: WDampLim                      !< Dry/Damp limit for water
      Real(8) :: WWetLim                       !< Damp/Wet limit for water
      Real(8) :: WWearLim                      !< Wear limit for water (below this
                                               !< only evaporation)
      Real(8) :: Snow2IceFac                   !< Snow to ice transition factor
      Real(8) :: SnowIceRat                    !< Snow-Ice ratio
      Real(8) :: MissValI                      !< Missing interger value
      Real(8) :: MissValR                      !< Missing real value
      real(8) :: MinPrecmm                     !< MIN precipitation : mm/hour
      real(8) :: MinWatmms                     !< MIN water storage : mm
      real(8) :: MinSnowmms                    !< MIN snow storage : mm   
                                               !< (accounts for snow "drift")
      real(8) :: MinDepmms                     !< MIN deposit storage : mm
      real(8) :: MinIcemms                     !< MIN ice storage : mm (MinPrec =>
                                               !< condens won't freeze(?))
      real(8) :: MaxSnowmms                    !< MAX snow storage : mm (plowed
                                               !< away if above)
      real(8) :: MaxDepmms                     !< MAX deposit storage : mm
      real(8) :: MaxIcemms                     !< MAX ice storage : mm
      real(8) :: MaxExtmms                     !< MAX extra (surface) water content (mm)
      real(8) :: MaxWatmms                     !< MAX water content
      real(8) :: AlbDry                        !< Dry asphalt albedo
      real(8) :: AlbSnow                       !< Snow albedo
      Real(8) :: WatDens                       !< Density of water at 0 C
      Real(8) :: SnowDens                      !< Density of snow (fresh, Oke p.44)
      Real(8) :: IceDens                       !< Density of ice (Oke p.44)
      Real(8) :: DepDens                       !< Density of deposit
      Real(8) :: WatMHeat                      !< Heat of ablation for water (J/kg
      Real(8) :: PorEvaF                       !< Pore resistance factor for evaporation
      Real(8) :: DampWearF                     !< Damp surface (poer) wear factor
      Real(8) :: TLimFreeze                    !< Freezing limit (C)
      Real(8) :: TLimMeltSnow                  !< Melting limit for snow (C)
      Real(8) :: TLimMeltIce                   !< Melting limit for ice (C)
      Real(8) :: TLimMeltDep                   !< Melting limit for deposit (C)
      Real(8) :: TLimDew                       !< Dew/deposit formation limit (C)
      Real(8) :: TLimColdH                     !< Higher limit for cold ground
                                               !< temperature
      Real(8) :: TLimColdL                     !< Lower limit for cold ground
                                               !< temperature
      Real(8) :: WetSnowFormR                  !< Water to snow ratio for wet snow
                                               !< formation
      Real(8) :: WetSnowMeltR                  !< Water to snow ratio for wet snow
                                               !< melting
      Real(8) :: PLimSnow                      !< Precipitation interpretation :
                                               !< snow limit
      Real(8) :: PLimRain                      !< Precipitation interpretation :
                                               !< rain limit

      Logical :: WetSnowFrozen                 !< True if wet snow is frozen
      Real(8) :: freezing_limit_normal         !< Freezing limit
      Real(8) :: snow_melting_limit_normal     !< Snow melting limit
      Real(8) :: ice_melting_limit_normal      !< ice melting limit
      Real(8) :: frost_melting_limit_normal    !< frost melting limit
      Real(8) :: frost_formation_limit_normal  !< frost formation limit
      Real(8) :: T4Melt_normal                 !< limit for melting normal
      Logical :: forceIceMelting               !< True if ice melting is forced
      Logical :: forceSnowMelting              !< True if snow melting is forced
      Logical :: CanMeltingChangeTemperature   !< True if surface temperature is changed
                                               !< in melting
end type RoadCondParameters
