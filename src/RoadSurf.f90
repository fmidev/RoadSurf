!MIT License
!Copyright (c) 2023 FMI Open Development
Module RoadSurf
   implicit none

   Interface
      !PUBLIC
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
      end subroutine ConnectFortran2Carrays

      !> Initializes model variables and parameters
      module Subroutine Initialization(modelInput, inSettings, settings, &
                                modelOutput, atm, surf, inputParam, localParam,&
                                coupling, phy, ground,condParam)
      
         USE, INTRINSIC :: ISO_C_BINDING
         use RoadSurfVariables
      
         type(InputSettings), intent(IN) :: inSettings    !< model settings
                                                                  !< given as input by
                                                                  !< modelRunner.cpp
      
         type(InputParameters), intent(IN) :: inputParam      !< input parameters
                                                              !< given by modelRunner.cpp
         type(LocalParameters), intent(IN) :: localParam      !< Local parameters
                                                               !< given by modelrunner.cpp
      
         type(inputArrays), intent(OUT) :: modelInput         !< Arrays for model
                                                              !< input data
         type(OutputArrays), intent(OUT) :: modelOutput       !< Arrays for model
                                                              !< output data
         type(AtmVariables), intent(OUT) :: atm               !< Variables for
                                                              !< atmospheric properties
         type(CouplingVariables), intent(OUT) :: coupling     !< variables used in
                                                              !< coupling(adjusting
                                                              !< radiation to fit
                                                              !< observed surface
                                                              !< temperature)
      
         type(ModelSettings), intent(OUT) :: settings         !< Variables for model
                                                              !< settings
      
         type(PhysicalParameters), intent(OUT) :: phy         !< Physical paremeters
                                                              !< used in the model
         type(GroundVariables), intent(OUT) :: ground         !< Varibales for ground
                                                              !< properties
         type(SurfaceVariables), intent(OUT) :: surf          !< Variables for surface
                                                              !< properties
         type(RoadCondParameters), intent(OUT) :: condParam   !< Parameters to determine storage
                                                              !< terms and road condition
      end subroutine Initialization      

     !> Checks input data for abnormal values
      module Subroutine CheckValues(modelInput, i, settings, surf,localParam)
         use RoadSurfVariables
      
         type(inputArrays), intent(INOUT) :: modelInput  !< Arrays for model input data
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(ModelSettings), intent(INOUT) :: settings !< Variables for model settings
         type(LocalParameters), intent(IN) :: localParam  !< local parameters given by
                                                          !< modelRunner.cpp
      end subroutine CheckValues

     !>Handles start and end of coupling and radiation coefficient calculation after
     !> coupling
      module Subroutine CouplingOperations1(i, coupling, surf, settings, ground, modelInput,&
                                     CP,localParam)
         use RoadSurfVariables
      
         type(ModelSettings), intent(IN) :: settings          !< Variables for model
                                                              !< settings
         type(inputArrays), intent(INOUT) :: modelInput       !< Arrays for model
                                                              !< input data
         type(RoadCondParameters), intent(IN) :: CP           !< Parameters to
                                                              !< determine storage
                                                              !< terms and road condition
      
         integer, intent(INOUT) :: i                          !< index for point in
                                                              !< input data
         type(CouplingVariables), intent(INOUT) :: coupling   !< variables used in
                                                              !< coupling(adjusting
                                                              !< radiation to fit
                                                              !< observed surface
                                                              !< temperature)
      
         type(SurfaceVariables), intent(INOUT) :: surf        !< Variables for surface
                                                              !< properties
         type(GroundVariables), intent(INOUT) :: ground       !< Varibales for ground
                                                              !< properties
         type(localParameters),intent(IN) :: localParam       !< local parameters
      end subroutine CouplingOperations1

   !> Use relaxation to air temperature, wind speed and relative humidity
   !> after initialization phase. This is done to avoid jump when moving
   !> from observed atmospheric values to forecasted ones.
      module Subroutine RelaxationOperations(i, atm, settings,Tmp)
         use RoadSurfVariables
   
         integer, intent(IN) :: i                     !< index of inputdata time steps
         type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
         type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                      !< properties
         real(8), dimension(0:16), intent(INOUT)::Tmp    !< Temperatures for each layer
      end subroutine RelaxationOperations

      !>set values to Tair etc
      module Subroutine SetCurrentValues(i, Tmp, modelInput, atm, settings, surf,coupling,ground)
         use RoadSurfVariables
      
         integer, intent(IN) :: i                     !< index of inputdata time steps
      
         type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
         type(inputArrays), intent(IN) :: modelInput  !< Arrays for model input data
         type(CouplingVariables), intent(IN) :: coupling !< variables used in
                                                         !< coupling(adjusting
                                                         !< radiation to fit observed
                                                         !< surface temperature)
      
         type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                      !< properties
         type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
         type(GroundVariables), intent(INOUT) :: ground   !< Varibales for ground
                                                          !< properties
         real(8), dimension(0:16), intent(INOUT)::Tmp    !< Temperatures for each layer

      end subroutine SetCurrentValues
      !>Calculates values for next time step using heat balance model
      module Subroutine BalanceModelOneStep(SWi, LWi, phy, ground, surf, atm, &
                                      settings, coupling, modelInput,&
                                      inputIdx,condParam)
         use RoadSurfVariables
      
         real(8), intent(IN) :: SWi                       !< Downwelling short wave
                                                          !< radiation
         real(8), intent(IN) :: LWi                       !< Downwelling long wave
                                                          !< radiation
      
         type(PhysicalParameters), intent(INOUT) :: phy   !< Physical paremeters used
                                                          !< in the model
         type(CouplingVariables), intent(IN) :: coupling  !< variables used in coupling
         type(inputArrays), intent(IN) :: modelInput      !< Arrays for model input data
         type(GroundVariables), intent(INOUT) :: ground   !< Varibales for ground
                                                          !< properties
         type(SurfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                          !< properties
         type(AtmVariables), intent(INOUT) :: atm         !< Variables for atmospheric
                                                          !< properties
         type(ModelSettings), intent(INOUT) :: settings   !< Variables for model settings
         type(RoadCondParameters), intent(IN) :: condParam  !< Parameters to
                                                            !< determine storage
                                                            !< terms and road
                                                            !< condition
      
         integer, intent(IN)::inputIdx                    !< Index in input data
      end subroutine BalanceModelOneStep
      !>Save values to output arrays
      module Subroutine SaveOutput(modelOutput, i, surf)
      
         use RoadSurfVariables
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
                                                      !< properties
         type(OutputArrays), intent(INOUT) :: modelOutput !< Arrays for model input data
      end subroutine
      !> check if at the end of coupling period
      module Subroutine CheckEndCoupling(i, settings, coupling, surf)
         use RoadSurfVariables
         integer, intent(IN) :: i                     !<index for point in input data
      
         type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
         type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
         type(CouplingVariables), intent(INOUT) :: coupling !< variables used in
                                                            !< coupling(adjusting
                                                            !< radiation to fit
                                                            !< observed surface temperature)
      end subroutine

      !> Determine precipitation type and add to storage
      module Subroutine PrecipitationToStorage(settings,CP,PrecPhase,atm,surf)
         use RoadSurfVariables
         
         type(ModelSettings),intent(IN)::settings      !> Variables for model settings
         type(RoadCondParameters),intent(IN)::CP       !< Parameters to determine
                                                       !< storage terms and road 
                                                       !< condition
         integer, intent(IN) :: PrecPhase              !< precPhase (Hail = 6;
                                                       !< FreezingRain = 5; 
                                                       !< FreezingDrizzle= 4; 
                                                       !< Snow = 3; Sleet = 2;
                                                       !< Rain = 1; Drizzle = 0;)
      
         type(AtmVariables),Intent(INOUT)::atm         !< Variables for atmpheric properties
         type(SurfaceVariables), Intent(INOUT) :: surf !< Variables for surface properties
      end subroutine PrecipitationToStorage
      !> Uses sky view factor and local horizon angles to modify incoming
      !> radiation fluxes
      module Subroutine ModRadiationBySurroundings(modelInput,inputParam,localParam,i)
      use RoadSurfVariables
         type(inputArrays),intent(INOUT) ::modelInput !< model input arrays
         type(inputParameters),intent(IN) :: inputParam !< model input parameters
         type(LocalParameters), intent(IN) :: localParam  !< Local parameters
                                                          !< given by modelrunner.cpp
         integer, intent(IN) ::i            !< model input index
      end subroutine ModRadiationBySurroundings
      !> Determine wear factors (how much storage terms are reduced by traffic)
      module Subroutine WearFactors(Snow2IceFac, Tph, surf, wearF)
         use RoadSurfVariables
         real(8), intent(IN)    :: Tph                   !< time steps per hour
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(WearingFactors), intent(OUT) :: wearF   !< wearing factors
         Real(8), intent(INOUT)    :: Snow2IceFac        !< Snow to ice transition factor
      end subroutine WearFactors
      !> Determines road surface condition for program Simulation      !
      module Subroutine RoadCond(MaxPormms, surf, atm, settings, &
                          CP,wearF)
         use RoadSurfVariables
      
         real(8), intent(IN)    :: MaxPormms                 !< maximum water in asphalt pores
         type(ModelSettings), intent(IN) :: settings      !< Variables for model settings
         type(RoadCondParameters), intent(INOUT) :: CP    !< Parameters to determine
                                                          !< storage terms and road
                                                          !< condition
      
         type(SurfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                          !< properties
         type(AtmVariables), Intent(INOUT) :: atm         !< Variables for atmospheric
                                                          !< properties
         type(WearingFactors),intent(IN) :: wearF         !< wearing factors
      end subroutine RoadCond
      !> Calculates albedo
      module Subroutine CalcAlbedo(albedo, surf, cp)
         use RoadSurfVariables
      
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                      !< storage terms and road condition
      
         real(8), intent(INOUT) :: Albedo                !< surface albedo
      end subroutine CalcAlbedo
 end Interface

   

   public :: ConnectFortran2Carrays
   public :: Initialization  
   public :: CheckValues
   public :: CouplingOperations1
   public :: RelaxationOperations
   public :: SetCurrentValues
   public :: BalanceModelOneStep
   public :: SaveOutput
   public :: CheckEndCoupling
   public :: PrecipitationToStorage
   public :: ModRadiationBySurroundings
   public :: WearFactors
   public :: RoadCond
   public :: CalcAlbedo

end Module RoadSurf
