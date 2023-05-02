Submodule (RoadSurf) ConnectArrays

   Implicit None
   contains
      !> Connect C pointers to fortran arrays
      module Subroutine ConnectFortran2Carrays(inputPointers,modelInput,&
                 outputPointers,modelOutput)
         use RoadSurfVariables
      
         type(DataPointers), intent(IN) :: inputPointers      !< pointers to input
                                                              !< data arrays given by
                                                              !< modelRunner.cpp
         type(OutputDataPointers), intent(INOUT) :: outputPointers !< pointers to
                                                              !< output data arrays
         type(InputArrays),intent(OUT) :: modelInput          !< Arrays for model input data
                                                              !< input data
         type(outputArrays), intent(OUT) :: modelOutput       !< Arrays for model output data
                                                              !< output data arrays
         !Set data from C++ to fortran arrays
         call SetData2InputArrays(inputPointers, modelInput)
         !Connect pointers to outputArrays given by C++ side
         call ConnectOutputArrays2Pointers(outputPointers, modelOutput) 
      
      end subroutine ConnectFortran2Carrays

end submodule ConnectArrays
      !> saves input data from modelRunner.cpp to fortran arrays
Subroutine SetData2InputArrays(inputPointers, modelInput)
   USE, INTRINSIC :: ISO_C_BINDING
   use RoadSurfVariables
   implicit none
   type(DataPointers), intent(IN) :: inputPointers !< pointers to input data
                                                   !< arrays given by modelRunner.cpp
   type(InputArrays), intent(OUT) :: modelInput !< Arrays for model input data

   CALL C_F_POINTER(inputPointers%c_tair, modelInput%Tair, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_tdew, modelInput%Tdew, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_VZ, modelInput%VZ, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_Rhz, modelInput%Rhz, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_prec, modelInput%prec, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_SW, modelInput%SW, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_LW, modelInput%LW, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_SW_dir, modelInput%SW_dir, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_LW_net, modelInput%LW_net, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_TSurfObs, modelInput%TSurfObs,&
    (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_PrecPhase, modelInput%PrecPhase,&
    (/inputPointers%inputLen/))

   CALL C_F_POINTER(inputPointers%c_local_horizons, modelInput%local_horizons, (/360/))
   CALL C_F_POINTER(inputPointers%c_Depth, modelInput%depth,&
    (/inputPointers%inputLen/))

   CALL C_F_POINTER(inputPointers%c_year, modelInput%year, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_month, modelInput%month, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_day, modelInput%day, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_hour, modelInput%hour, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_minute, modelInput%minute, (/inputPointers%inputLen/))
   CALL C_F_POINTER(inputPointers%c_second, modelInput%second, (/inputPointers%inputLen/))
end subroutine

!> connect output arrays to pointers given by modelRunner.cpp
Subroutine ConnectOutputArrays2Pointers(outputPointers, modelOutput)
   USE, INTRINSIC :: ISO_C_BINDING
   use RoadSurfVariables
   implicit none
   type(OutputDataPointers), intent(INOUT) :: outputPointers !< pointers to
                                                             !< output data arrays
   type(outputArrays), intent(OUT) :: modelOutput !< Arrays for model output data

   CALL C_F_POINTER(outputPointers%c_TsurfOut, modelOutput%TsurfOut, &
                    (/outputPointers%outputLen/))
   CALL C_F_POINTER(outputPointers%c_SnowOut, modelOutput%SnowOut, &
                    (/outputPointers%outputLen/))
   CALL C_F_POINTER(outputPointers%c_WaterOut, modelOutput%WaterOut, &
                    (/outputPointers%outputLen/))
   CALL C_F_POINTER(outputPointers%c_IceOut, modelOutput%IceOut, &
                    (/outputPointers%outputLen/))
   CALL C_F_POINTER(outputPointers%c_DepositOut, modelOutput%DepositOut, &
                    (/outputPointers%outputLen/))
   CALL C_F_POINTER(outputPointers%c_Ice2Out, modelOutput%Ice2Out, &
                    (/outputPointers%outputLen/))

end subroutine
