Submodule (RoadSurf) ConnectArrays

   Implicit None
   contains
      !> Connect C pointers to fortran arrays
      module Subroutine ConnectFortran2Carrays(inPointers,modelInput,&
                 outPointers,modelOutput)
         use RoadSurfVariables
      
         type(InputPointers), intent(IN) :: inPointers      !< pointers to input
                                                              !< data arrays given by
                                                              !< modelRunner.cpp
         type(OutputPointers), intent(INOUT) :: outPointers !< pointers to
                                                              !< output data arrays
         type(InputArrays),intent(OUT) :: modelInput          !< Arrays for model input data
                                                              !< input data
         type(OutputArrays), intent(OUT) :: modelOutput       !< Arrays for model output data
                                                              !< output data arrays
         !Set data from C++ to fortran arrays
         call SetData2InputArrays(inPointers, modelInput)
         !Connect pointers to OutputArrays given by C++ side
         call ConnectOutputArrays2Pointers(outPointers, modelOutput) 
      
      end subroutine ConnectFortran2Carrays

end submodule ConnectArrays
      !> saves input data from modelRunner.cpp to fortran arrays
Subroutine SetData2InputArrays(inPointers, modelInput)
   USE, INTRINSIC :: ISO_C_BINDING
   use RoadSurfVariables
   implicit none
   type(InputPointers), intent(IN) :: inPointers !< pointers to input data
                                                   !< arrays given by modelRunner.cpp
   type(InputArrays), intent(OUT) :: modelInput !< Arrays for model input data

   CALL C_F_POINTER(inPointers%c_tair, modelInput%Tair, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_tdew, modelInput%Tdew, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_VZ, modelInput%VZ, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_Rhz, modelInput%Rhz, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_prec, modelInput%prec, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_SW, modelInput%SW, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_LW, modelInput%LW, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_SW_dir, modelInput%SW_dir, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_LW_net, modelInput%LW_net, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_TSurfObs, modelInput%TSurfObs,&
    (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_PrecPhase, modelInput%PrecPhase,&
    (/inPointers%inputLen/))

   CALL C_F_POINTER(inPointers%c_local_horizons, modelInput%local_horizons, (/360/))
   CALL C_F_POINTER(inPointers%c_Depth, modelInput%depth,&
    (/inPointers%inputLen/))

   CALL C_F_POINTER(inPointers%c_year, modelInput%year, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_month, modelInput%month, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_day, modelInput%day, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_hour, modelInput%hour, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_minute, modelInput%minute, (/inPointers%inputLen/))
   CALL C_F_POINTER(inPointers%c_second, modelInput%second, (/inPointers%inputLen/))
end subroutine

!> connect output arrays to pointers given by modelRunner.cpp
Subroutine ConnectOutputArrays2Pointers(outPointers, modelOutput)
   USE, INTRINSIC :: ISO_C_BINDING
   use RoadSurfVariables
   implicit none
   type(OutputPointers), intent(INOUT) :: outPointers !< pointers to
                                                             !< output data arrays
   type(OutputArrays), intent(OUT) :: modelOutput !< Arrays for model output data

   CALL C_F_POINTER(outPointers%c_TsurfOut, modelOutput%TsurfOut, &
                    (/outPointers%outputLen/))
   CALL C_F_POINTER(outPointers%c_SnowOut, modelOutput%SnowOut, &
                    (/outPointers%outputLen/))
   CALL C_F_POINTER(outPointers%c_WaterOut, modelOutput%WaterOut, &
                    (/outPointers%outputLen/))
   CALL C_F_POINTER(outPointers%c_IceOut, modelOutput%IceOut, &
                    (/outPointers%outputLen/))
   CALL C_F_POINTER(outPointers%c_DepositOut, modelOutput%DepositOut, &
                    (/outPointers%outputLen/))
   CALL C_F_POINTER(outPointers%c_Ice2Out, modelOutput%Ice2Out, &
                    (/outPointers%outputLen/))

end subroutine
