!> pointers to output data arrays given by modelRunner.cpp
Type, Bind(C) :: OutputDataPointers

      INTEGER(C_INT) :: outputLen               !<output data length

      TYPE(C_PTR) :: c_TsurfOut;                !< pointer to output surface
                                                !< temperature array
      TYPE(C_PTR) :: c_SnowOut;                 !< pointer to output snow array
      TYPE(C_PTR) :: c_WaterOut;                !< pointer to output water array
      TYPE(C_PTR) :: c_IceOut;                  !< pointer to output ice array
      TYPE(C_PTR) :: c_DepositOut;              !< pointer to output deposit array
      TYPE(C_PTR) :: c_Ice2Out;                 !< pointer to output secondary
                                                !< ice array

end type OutputDataPointers
