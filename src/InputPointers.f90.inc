!MIT License
!Copyright (c) 2023 FMI Open Development
!> pointers to input data arrays given by modelRunner.cpp
Type, Bind(C) :: InputPointers
      INTEGER(C_INT) :: inputLen        !<input data length

      TYPE(C_PTR) :: c_tair             !< pointer to input air temperature array
      TYPE(C_PTR) :: c_tdew             !< pointer to input dew point temperature array
      TYPE(C_PTR) :: c_VZ               !< pointer to input wind speed array
      TYPE(C_PTR) :: c_Rhz              !< pointer to input relative humidity array
      TYPE(C_PTR) :: c_prec             !< pointer to input precipitation array
      TYPE(C_PTR) :: c_SW               !< pointer to input short wave radiation array
      TYPE(C_PTR) :: c_LW               !< pointer to input long wave array
      TYPE(C_PTR) :: c_SW_dir           !< pointer to input direct short wave radiation array
      TYPE(C_PTR) :: c_LW_net           !< pointer to input net long wave array
      TYPE(C_PTR) :: c_TSurfObs         !< pointer to input tsurf obs array
      TYPE(C_PTR) :: c_PrecPhase        !< pointer to input prec phase array
      TYPE(C_PTR) :: c_local_horizons   !< pointer to local horizon angle array
      TYPE(C_PTR) :: c_Depth            !< pointer to input depth array
      TYPE(C_PTR) :: c_year             !< pointer to input year array
      TYPE(C_PTR) :: c_month            !< pointer to input month array
      TYPE(C_PTR) :: c_day              !< pointer to input day array
      TYPE(C_PTR) :: c_hour             !< pointer to input hour array
      TYPE(C_PTR) :: c_minute           !< pointer to input minute array
      TYPE(C_PTR) :: c_second           !< pointer to input second array

end type InputPointers

