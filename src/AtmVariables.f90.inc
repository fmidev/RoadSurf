!> Variables for atmospheric properties
type atmVariables
      real(8) :: Tair                  !< air temperature (C)
      real(8) :: VZ                    !< wind speed (m/s)
      real(8) :: Tdew                  !<dew point temperature (C)
      real(8) :: RHz                   !< relative humidity (%)
      real(8) :: PrecInTStep           !< precipitation in model time step
      real(8) :: TairInitEnd           !< Air temperature at the end of
                                       !< initialization (C)
      real(8) :: VZInitEnd             !< wind speed at the end of initialization
                                       !< (m/s)
      real(8) :: RhzInitEnd            !< relative humidity at the end of
                                       !< initialization (%)
      real(8) :: BLCond                !< boundary layer conductance
      real(8) :: RNet                  !< Net radiation (W/m2)
      real(8) :: LE_Flux               !< latent heat flux (W/m2)
      real(8) :: RainIntensity         !< rain intensity (mm/h)
      real(8) :: SnowIntensity         !< snow intensity (mm/h)
      real(8) :: TairR                 !< Air temperature from forecast at
                                       !< the end of initialization (C)
      real(8) :: VZR                   !< wind speed from forecastat the
                                       !< end of initialization (m/s)
      real(8) :: RhzR                  !< relative humidity from forecast
                                       !< at the end of initialization (%)
      real(8) :: CalmLim               !< Minimum wind speed
      real(8) :: SensibleHeatFlux      !< sensible heat flux (W/m2)
      real(8) :: RainmmTS              !< amount of rain in time step
      real(8) :: SnowmmTS              !< amount of snow in time step
      Integer :: SnowType              !< snow type (wet or dry)
      Integer :: PrecType              !< Precipitation type 
                                       !< (1=Water,2=Sleet,3=Snow)
end type atmVariables
