#include "Constants.h"
#include "ConstantsExtra.h"

!>Runs road weather model simulation
SUBROUTINE runsimulation(outputPointers, inputPointers,&
                         inputSettings, inputParam,&
                         localParam) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING
   use RoadSurfVariables
   use RoadSurf
   Implicit None

   type(OutputDataPointers), intent(INOUT) :: outputPointers !< pointers to
                                                             !< output data arrays

   type(DataPointers), intent(IN) :: inputPointers          !< pointers to input
                                                            !< data arrays given
                                                            !< by modelRunner.cpp

   type(InputModelSettings), intent(IN) :: inputSettings    !< model settings
                                                            !< given as input by
                                                            !< modelRunner.cpp

   type(InputParameters), intent(IN) :: inputParam          !< input parameters
                                                            !< given by modelRunner.cpp
   type(LocalParameters), intent(IN) :: localParam          !< Local parameters
                                                            !< given by modelrunner.cpp
   integer:: i                                  !< index of inputdata time steps

   type(InputArrays) :: modelInput              !< Arrays for model input data
   type(outputArrays) :: modelOutput            !< Arrays for model output data
   type(physicalParameters) :: phy              !< Physical paremeters used in
                                                !< the model
   type(groundVariables) :: ground              !< Varibales for ground properties
   type(surfaceVariables) :: surf               !< Variables for surface properties
                                                !< outside of car tyre tracks
   type(atmVariables) :: atm                    !< Variables for atmospheric
                                                !< properties
                                                !< properties
   type(couplingVariables) :: coupling          !< variables used in coupling(adjusting
                                                !< radiation to fit observed
                                                !< surface temperature)

   type(modelSettings) :: settings              !< Variables for model settings
                                                !< conditions
   type(roadCondParameters) :: condParam        !< Parameters to determine
                                                !< storage terms and road condition

!---------INITIALIZE----------------------------------------------------------
   call ConnectFortran2Carrays(inputPointers,modelInput,&
           outputPointers,modelOutput)

   call Initialization(modelInput, inputSettings, settings, &
                       modelOutput, atm, surf, inputParam, localParam,&
                       coupling, phy, ground, condParam) 
!---------START SIMULATION----------------------------------
   !Start temperature profile simulation, go trough input data values
   i = 1
   Do while (i < settings%SimLen .and. (settings%Simulation_Failed .eqv. .false.))
      Call checkValues(modelInput, i, settings, surf,localParam)

      !Check if coupling is on
      if (settings%use_coupling) Then
         !Check if in coupling phase, save variables if at the start of the
         ! coupling phase,
         !Load saved variables if coupling is started again.
         !After coupling, calculate radiation coefficients
         call CouplingOperations1(i, coupling, surf, settings, ground, modelInput,&
                                 CondParam,localParam)

      end if

      !set current values to atm%Tair etc
      call setCurrentValues(i, ground%Tmp, modelInput, atm, settings, surf, coupling,&
          ground)

      !If relaxation is used
      if (settings%use_relaxation) Then
         !Smooth t2m, rh and wind values when moving from initialization phase
         !to forecasting phase
         call relaxationOperations(i, atm, settings,ground%Tmp)

      end if
      !Calculate temperature profile and storage values one timestep forward
      call roadModelOneStep(i, phy, ground, surf, atm,&
                            settings, coupling, modelInput, condParam,&
                            inputParam, localParam)
      
      !Save output
      call saveOutput(modelOutput, i, surf)

      !Coupling control if at the end of the coupling period
      call checkEndCoupling(i, settings, coupling, surf)

      i = i + 1
   end do

   !If simulation is not failed, make calculations for last value
   !and save output for last time step if index matches wanted
   !output time step
   if (settings%Simulation_Failed .eqv. .false.) Then

      !Make still calculation for i=SimLen (the last value)

      !set last input values as interpolated values
      call lastValues(modelInput, atm, settings, ground, surf)

      !Calculate temperature profile and storage values one time step forward
      call roadModelOneStep(settings%SimLen, phy, ground, surf, atm,&
                            settings, coupling, modelInput, condParam,&
                            inputParam,localParam)

      !Save output
      call saveOutput(modelOutput, i, surf)
!--------SIMULATION END----------------------------
   end if

End Subroutine

!>Executes one time step of heat balance model and road condition calculation
Subroutine roadModelOneStep(input_idxI, phy, ground, surf, atm,&
                            settings, coupling, modelInput, condParam,&
                            inputParam,localParam)

   use RoadSurfVariables
   use RoadSurf

   Implicit None
   integer, intent(IN) :: input_idxI                !< current index I in input data
   type(inputArrays), intent(INOUT) :: modelInput      !< Arrays for model input data
   type(couplingVariables), intent(IN) :: coupling  !< variables used in coupling
   type(modelSettings), intent(INOUT) :: settings   !< Variables for model settings
   type(physicalParameters), intent(INOUT) :: phy   !< Physical paremeters used
                                                    !< in the model
   type(groundVariables), intent(INOUT) :: ground   !< Varibales for ground
                                                    !< properties
   type(surfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                    !< properties
   type(atmVariables), intent(INOUT) :: atm         !< Variables for atmospheric
                                                    !< properties

   type(roadCondParameters), intent(INOUT) :: condParam     !< Parameters to
                                                            !< determine storage
                                                            !< terms and road
                                                            !< condition
   type(inputParameters), intent(IN) :: inputParam  !< model input parameters
   type(LocalParameters), intent(IN) :: localParam  !< Local parameters
                                                    !< given by modelrunner.cpp
   type(wearingFactors) :: wearF   !< wearing factors

   !Determine wheter the precipitation is rain or snow and add to storage
   call PrecipitationToStorage(settings,condParam,modelInput%PrecPhase(input_idxI),&
                               atm,surf)
   !Make radiation corrections based on sky view and local horizon angles
   if (localParam%sky_view<1.0 .and. localParam%sky_view>-0.01)then
      call modRadiationBySurroundings(modelInput,inputParam,localParam,input_idxI)
   end if
   !Calculate temperature profile one time step forward
   !Checks also for melting (can affect temperature)
   call balanceModelOneStep(modelInput%SW(input_idxI), &
                            modelInput%LW(input_idxI), &
                            phy, ground, surf, atm, settings, coupling, modelInput,&
                            input_idxI,condParam)
   call wearFactors(condParam%Snow2IceFac, settings%Tph, surf, wearF)
  ! ************* WEAR FACTORS
   !Calculate storage terms and road condition
   call roadCond(phy%MaxPormms, surf, atm, settings, &
                 condParam,wearF)

   ! *************  ALBEDO
   call calcAlbedo(ground%Albedo, surf, condParam)
 
end Subroutine
