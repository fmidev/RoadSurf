!> Arrays for model output data
type outputArrays
      REAL(C_DOUBLE), POINTER :: TsurfOut(:)            !< Road surface temperature (C)
      REAL(C_DOUBLE), POINTER :: SnowOut(:)             !< snow storage (water
                                                        !< equivalent mm)
      REAL(C_DOUBLE), POINTER :: WaterOut(:)            !< water storage (mm)
      REAL(C_DOUBLE), POINTER :: IceOut(:)              !< ice storage (water
                                                        !< equivalent mm)
      REAL(C_DOUBLE), POINTER :: DepositOut(:)          !< deposit storage
                                                        !< (water equivalent mm)
      REAL(C_DOUBLE), POINTER :: Ice2Out(:)             !< secondary ice storage


end type outputArrays