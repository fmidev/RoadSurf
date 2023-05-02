!> model settings given as input by modelRunner.cpp
Type, Bind(C) :: InputModelSettings
      INTEGER(C_INT) :: SimLen                  !< Lenght of simulation
      INTEGER(C_INT) :: use_coupling            !< 1 if coupling is used
      INTEGER(C_INT) :: use_relaxation          !< 1 if relaxation is used
      REAL(C_DOUBLE) :: DTSecs                  !< time step in seconds
      REAL(C_DOUBLE) :: tsurfOutputDepth        !< depth to interpolate surface temperature (m)
      INTEGER(C_INT) :: coupling_minutes        !< coupling length in minutes
      REAL(C_DOUBLE) :: couplingEffectReduction !< Parameter used to calculate radiation
                                                !< coefficient after coupling

      INTEGER(C_INT) :: outputStep              !< model output frequency in minutes
end type InputModelSettings
