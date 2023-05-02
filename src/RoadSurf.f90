Module RoadSurf
   implicit none

   Interface
      !PUBLIC
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
      end subroutine ConnectFortran2Carrays

      !> Initializes model variables and parameters
      module Subroutine Initialization(modelInput, inputSettings, settings, &
                                modelOutput, atm, surf, inputParam, localParam,&
                                coupling, phy, ground,condParam)
      
         USE, INTRINSIC :: ISO_C_BINDING
         use RoadSurfVariables
      
         type(InputModelSettings), intent(IN) :: inputSettings    !< model settings
                                                                  !< given as input by
                                                                  !< modelRunner.cpp
      
         type(InputParameters), intent(IN) :: inputParam      !< input parameters
                                                              !< given by modelRunner.cpp
         type(LocalParameters), intent(IN) :: localParam      !< Local parameters
                                                               !< given by modelrunner.cpp
      
         type(inputArrays), intent(OUT) :: modelInput         !< Arrays for model
                                                              !< input data
         type(outputArrays), intent(OUT) :: modelOutput       !< Arrays for model
                                                              !< output data
         type(atmVariables), intent(OUT) :: atm               !< Variables for
                                                              !< atmospheric properties
         type(couplingVariables), intent(OUT) :: coupling     !< variables used in
                                                              !< coupling(adjusting
                                                              !< radiation to fit
                                                              !< observed surface
                                                              !< temperature)
      
         type(modelSettings), intent(OUT) :: settings         !< Variables for model
                                                              !< settings
      
         type(physicalParameters), intent(OUT) :: phy         !< Physical paremeters
                                                              !< used in the model
         type(groundVariables), intent(OUT) :: ground         !< Varibales for ground
                                                              !< properties
         type(surfaceVariables), intent(OUT) :: surf          !< Variables for surface
                                                              !< properties
         type(roadCondParameters), intent(OUT) :: condParam   !< Parameters to determine storage
                                                              !< terms and road condition
      end subroutine Initialization      

     !> Checks input data for abnormal values
      module Subroutine checkValues(modelInput, i, settings, surf,localParam)
         use RoadSurfVariables
      
         type(inputArrays), intent(INOUT) :: modelInput  !< Arrays for model input data
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(surfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(modelSettings), intent(INOUT) :: settings !< Variables for model settings
         type(LocalParameters), intent(IN) :: localParam  !< local parameters given by
                                                          !< modelRunner.cpp
      end subroutine checkValues

     !>Handles start and end of coupling and radiation coefficient calculation after
     !> coupling
      module Subroutine CouplingOperations1(i, coupling, surf, settings, ground, modelInput,&
                                     CP,localParam)
         use RoadSurfVariables
      
         type(modelSettings), intent(IN) :: settings          !< Variables for model
                                                              !< settings
         type(inputArrays), intent(INOUT) :: modelInput       !< Arrays for model
                                                              !< input data
         type(roadCondParameters), intent(IN) :: CP           !< Parameters to
                                                              !< determine storage
                                                              !< terms and road condition
      
         integer, intent(INOUT) :: i                          !< index for point in
                                                              !< input data
         type(couplingVariables), intent(INOUT) :: coupling   !< variables used in
                                                              !< coupling(adjusting
                                                              !< radiation to fit
                                                              !< observed surface
                                                              !< temperature)
      
         type(surfaceVariables), intent(INOUT) :: surf        !< Variables for surface
                                                              !< properties
         type(groundVariables), intent(INOUT) :: ground       !< Varibales for ground
                                                              !< properties
         type(localParameters),intent(IN) :: localParam       !< local parameters
      end subroutine CouplingOperations1

   !> Use relaxation to air temperature, wind speed and relative humidity
   !> after initialization phase. This is done to avoid jump when moving
   !> from observed atmospheric values to forecasted ones.
      module Subroutine relaxationOperations(i, atm, settings,Tmp)
         use RoadSurfVariables
   
         integer, intent(IN) :: i                     !< index of inputdata time steps
         type(modelSettings), intent(IN) :: settings  !< Variables for model settings
         type(atmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                      !< properties
         real(8), dimension(0:16), intent(INOUT)::Tmp    !< Temperatures for each layer
      end subroutine relaxationOperations

      !>set values to Tair etc
      module Subroutine setCurrentValues(i, Tmp, modelInput, atm, settings, surf,coupling,ground)
         use RoadSurfVariables
      
         integer, intent(IN) :: i                     !< index of inputdata time steps
      
         type(modelSettings), intent(IN) :: settings  !< Variables for model settings
         type(inputArrays), intent(IN) :: modelInput  !< Arrays for model input data
         type(couplingVariables), intent(IN) :: coupling !< variables used in
                                                         !< coupling(adjusting
                                                         !< radiation to fit observed
                                                         !< surface temperature)
      
         type(atmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                      !< properties
         type(surfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
         type(groundVariables), intent(INOUT) :: ground   !< Varibales for ground
                                                          !< properties
         real(8), dimension(0:16), intent(INOUT)::Tmp    !< Temperatures for each layer

      end subroutine setCurrentValues
      !>Calculates values for next time step using heat balance model
      module Subroutine balanceModelOneStep(SWi, LWi, phy, ground, surf, atm, &
                                      settings, coupling, modelInput,&
                                      inputIdx,condParam)
         use RoadSurfVariables
      
         real(8), intent(IN) :: SWi                       !< Downwelling short wave
                                                          !< radiation
         real(8), intent(IN) :: LWi                       !< Downwelling long wave
                                                          !< radiation
      
         type(physicalParameters), intent(INOUT) :: phy   !< Physical paremeters used
                                                          !< in the model
         type(couplingVariables), intent(IN) :: coupling  !< variables used in coupling
         type(inputArrays), intent(IN) :: modelInput      !< Arrays for model input data
         type(groundVariables), intent(INOUT) :: ground   !< Varibales for ground
                                                          !< properties
         type(surfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                          !< properties
         type(atmVariables), intent(INOUT) :: atm         !< Variables for atmospheric
                                                          !< properties
         type(modelSettings), intent(INOUT) :: settings   !< Variables for model settings
         type(roadCondParameters), intent(IN) :: condParam  !< Parameters to
                                                            !< determine storage
                                                            !< terms and road
                                                            !< condition
      
         integer, intent(IN)::inputIdx                    !< Index in input data
      end subroutine balanceModelOneStep
      !>Save values to output arrays
      module Subroutine saveOutput(modelOutput, i, surf)
      
         use RoadSurfVariables
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(surfaceVariables), intent(IN) :: surf   !< Variables for surface properties
                                                      !< properties
         type(outputArrays), intent(INOUT) :: modelOutput !< Arrays for model input data
      end subroutine
      !> check if at the end of coupling period
      module Subroutine checkEndCoupling(i, settings, coupling, surf)
         use RoadSurfVariables
         integer, intent(IN) :: i                     !<index for point in input data
      
         type(modelSettings), intent(IN) :: settings  !< Variables for model settings
         type(surfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
         type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                            !< coupling(adjusting
                                                            !< radiation to fit
                                                            !< observed surface temperature)
      end subroutine

      !> Determine precipitation type and add to storage
      module Subroutine PrecipitationToStorage(settings,CP,PrecPhase,atm,surf)
         use RoadSurfVariables
         
         type(modelSettings),intent(IN)::settings      !> Variables for model settings
         type(roadCondParameters),intent(IN)::CP       !< Parameters to determine
                                                       !< storage terms and road 
                                                       !< condition
         integer, intent(IN) :: PrecPhase              !< precPhase (Hail = 6;
                                                       !< FreezingRain = 5; 
                                                       !< FreezingDrizzle= 4; 
                                                       !< Snow = 3; Sleet = 2;
                                                       !< Rain = 1; Drizzle = 0;)
      
         type(atmVariables),Intent(INOUT)::atm         !< Variables for atmpheric properties
         type(surfaceVariables), Intent(INOUT) :: surf !< Variables for surface properties
      end subroutine PrecipitationToStorage
      !> Uses sky view factor and local horizon angles to modify incoming
      !> radiation fluxes
      module Subroutine modRadiationBySurroundings(modelInput,inputParam,localParam,i)
      use RoadSurfVariables
         type(inputArrays),intent(INOUT) ::modelInput !< model input arrays
         type(inputParameters),intent(IN) :: inputParam !< model input parameters
         type(LocalParameters), intent(IN) :: localParam  !< Local parameters
                                                          !< given by modelrunner.cpp
         integer, intent(IN) ::i            !< model input index
      end subroutine modRadiationBySurroundings
      !> Determine wear factors (how much storage terms are reduced by traffic)
      module Subroutine wearFactors(Snow2IceFac, Tph, surf, wearF)
         use RoadSurfVariables
         real(8), intent(IN)    :: Tph                   !< time steps per hour
         type(surfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(wearingFactors), intent(OUT) :: wearF   !< wearing factors
         Real(8), intent(INOUT)    :: Snow2IceFac        !< Snow to ice transition factor
      end subroutine wearFactors
      !> Determines road surface condition for program Simulation      !
      module Subroutine roadCond(MaxPormms, surf, atm, settings, &
                          CP,wearF)
         use RoadSurfVariables
      
         real(8), intent(IN)    :: MaxPormms                 !< maximum water in asphalt pores
         type(modelSettings), intent(IN) :: settings      !< Variables for model settings
         type(roadCondParameters), intent(INOUT) :: CP    !< Parameters to determine
                                                          !< storage terms and road
                                                          !< condition
      
         type(surfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                          !< properties
         type(atmVariables), Intent(INOUT) :: atm         !< Variables for atmospheric
                                                          !< properties
         type(wearingFactors),intent(IN) :: wearF         !< wearing factors
      end subroutine roadCond
      !> Calculates albedo
      module Subroutine calcAlbedo(albedo, surf, cp)
         use RoadSurfVariables
      
         type(surfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(roadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                      !< storage terms and road condition
      
         real(8), intent(INOUT) :: Albedo                !< surface albedo
      end subroutine calcAlbedo
 end Interface

   

   public :: ConnectFortran2Carrays
   public :: Initialization  
   public :: checkValues
   public :: CouplingOperations1
   public :: relaxationOperations
   public :: setCurrentValues
   public :: balanceModelOneStep
   public :: saveOutput
   public :: checkEndCoupling
   public :: PrecipitationToStorage
   public :: modRadiationBySurroundings
   public :: wearFactors
   public :: roadCond
   public :: calcAlbedo

end Module RoadSurf
