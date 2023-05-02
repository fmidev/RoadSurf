!> local parameters given by modelRunner.cpp
Type, Bind(C) :: LocalParameters
      REAL(C_DOUBLE) :: tair_relax;     !< tair for relaxation
      REAL(C_DOUBLE) :: VZ_relax;       !< wind speed for relaxation
      REAL(C_DOUBLE) :: RH_relax;       !< relative humidity for relaxation
      INTEGER(C_INT) :: couplingIndexI; !< index in input data where coupling starts
      REAL(C_DOUBLE) :: couplingTsurf;  !< surface temperature value used in coupling
      REAL(C_DOUBLE) :: lat;            !< latitude
      REAL(C_DOUBLE) :: lon;            !< longitude
      REAL(C_DOUBLE) :: sky_view        !< Sky view factor
      INTEGER(C_INT) :: InitLenI        !< The length of initialization
                                        !< period, input timesteps
end type LocalParameters
