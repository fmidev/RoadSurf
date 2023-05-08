!> Variables for model settings
type ModelSettings
      integer:: InitLenI                    !< The length of initialization period,
                                            !< input timesteps
      integer:: SimLen                      !< Lenght of simulation
      logical :: use_coupling               !< coupling is used
      logical :: use_relaxation             !< true if relaxation is used
                                                
      integer:: NLayers                     !< number of ground layers
      real(8):: DTSecs                        !< time step in seconds
      real(8) :: tsurfOutputDepth              !< depth to interpolate output surface temperature
      logical :: simulation_failed          !< true if simulation has been failed
      real(8) :: Tph                           !< time steps per hour
      real(8) :: NightOn                    !< Beginning hour of night traffic (UTC)
      real(8) :: NightOff                   !< Ending hour of night traffic (UTC)
      Real(8) :: CalmLimDay                   !< Minimum wind speed (m/s)
      Real(8)  :: CalmLimNgt                   !< Minimum wind speed (m/s)
      Real(8)  :: TrfFricNgt                   !< Traffic induced friction (W/m2)
      Real(8)  :: TrFfricDay                   !< Traffic induced friction (W/m2)
      Integer :: coupling_minutes           !< Lenght of the coupling period in minutes
      Real(8) :: couplingEffectReduction       !< Parameter used to calculate radiation
                                            !< coefficient after coupling

      Integer :: outputStep                 !< Frequency of output in minutes
end type ModelSettings
