!MIT License
!Copyright (c) 2023 FMI Open Development
!> input data arrays
Type :: InputArrays
      INTEGER(C_INT), POINTER :: timeForFortran(:) !< time as 1D array with
                                                   !< form(year1,month1,day1,hour1,
                                                   !< minute1,second1, year2,month2...)

      REAL(C_DOUBLE), POINTER :: Tair(:)        !< air temperature (C)
      REAL(C_DOUBLE), POINTER :: Tdew(:)        !< Dew point temperature (C)
      REAL(C_DOUBLE), POINTER :: VZ(:)          !< wind speed (m/s)
      REAL(C_DOUBLE), POINTER :: Rhz(:)         !< relative humidity (%)
      REAL(C_DOUBLE), POINTER :: prec(:)        !< precipitation (mm/h)
      REAL(C_DOUBLE), POINTER :: SW(:)          !< incoming short wave radiation
                                                !< (W/m2)
      REAL(C_DOUBLE), POINTER :: LW(:)          !< incoming long wave radiation
                                                !< (W/m2)
      REAL(C_DOUBLE), POINTER :: SW_dir(:)      !< direct short wave radiation
                                                !< (W/m2)
      REAL(C_DOUBLE), POINTER :: LW_net(:)      !< net long wave radiation
                                                !< (W/m2)
      REAL(C_DOUBLE), POINTER :: TSurfObs(:)    !< surface temperature  (C)
      INTEGER(C_INT), POINTER :: PrecPhase(:)   !< precipitation phase (Hail=6;
                                                !< FreezingRain=5; FreezingDrizzle=4;
                                                !< Snow=3; Sleet=2; Rain=1; Drizzle=0;
                                                
      REAL(C_DOUBLE),POINTER  :: local_horizons(:) !< local horizon angles
      REAL(C_DOUBLE), POINTER :: depth(:)       !< depth to calculate output temperature (m)
      INTEGER(C_INT), POINTER :: year(:)        !< year
      INTEGER(C_INT), POINTER :: month(:)       !< month
      INTEGER(C_INT), POINTER :: day(:)         !< day
      INTEGER(C_INT), POINTER :: hour(:)        !< hour
      INTEGER(C_INT), POINTER :: minute(:)      !< minute
      INTEGER(C_INT), POINTER :: second(:)      !< second
      

end type InputArrays
