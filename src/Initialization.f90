!MIT License
!Copyright (c) 2023 FMI Open Development
#include "Constants.h"

submodule (RoadSurf) Init
   Implicit none
   contains
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
      
         !Initialize settings
         call initSettings(inSettings, settings,inputParam,localParam)
         
         !Initialize output arrays
         call initOutputArrays(settings%SimLen, modelOutput)
      
         !set input parameters
         call setInputParam(localParam, atm, coupling, settings)
        
         !Initialize variables and parameters
         call initVariablesAndParameters(modelInput, coupling, phy, ground,&
               surf, atm, settings,inputParam,condParam)
        
      end subroutine Initialization
end submodule Init
      !> Initializes model values
Subroutine initVariablesAndParameters(modelInput, &
                                      coupling, phy, ground, surf, atm, settings,&
                                      inputParam,condParam)

   use RoadSurfVariables
   implicit none
   type(InputParameters), intent(IN) :: inputParam      !< input parameters
                                                        !< given by modelRunner.cpp
   type(inputArrays), intent(INOUT) :: modelInput       !< Arrays for model input data
   type(AtmVariables), intent(INOUT) :: atm             !< Variables for atmospheric properties
   type(CouplingVariables), intent(INOUT) :: coupling   !< variables used in coupling(adjusting
                                                        !< radiation to fit observed surface
                                                        !< temperature)

   type(ModelSettings), intent(INOUT) :: settings       !< Variables for model settings

   type(PhysicalParameters), intent(OUT) :: phy         !< Physical paremeters
                                                        !< used in the model
   type(GroundVariables), intent(OUT) :: ground         !< Varibales for ground properties
   type(SurfaceVariables), intent(INOUT) :: surf        !< Variables for surface properties
   type(RoadCondParameters), intent(OUT) :: condParam   !< Parameters to determine storage
                                                        !< terms and road condition

   real(8) :: depth
   real(8) ::t_output
   settings%simulation_failed = .false.

   settings%Tph = settings%DTSecs/3600.0 !Time steps per hour

   !Initialize ground layer depths, temperature profile, physical parameters
   !and  ground physical properties
   call allocator(settings,ground,coupling)
   call initDepth(ground, settings%NLayers)
   call initSurf(surf, .true.)
   call InitParam(ground%Albedo, phy,inputParam)

   call initTemp(modelInput%TsurfOBS(1), modelInput%Tair(1), &
                modelInput%depth(1), phy, settings, ground, modelInput, surf)

   call initVariables(ground, atm,settings)
   
   call HCapValues(phy)
   call ground_prop_init(settings%NLayers, ground)
   call CalcCC(settings%NLayers, phy, ground)
   call CalcHCapHCond(settings%NLayers, settings%DTSecs, phy, ground, atm)

   call calcCapDZCondDZ(settings%NLayers, ground)
   !Initialize coupling values
   call initCoupling(coupling)
   call initCouplingTimes(coupling, settings)

  !Initialize min and max limits and other parameters for road condition calculation
   call condInit(condParam, surf,inputParam)
 
  !Make sure that wind speed is at least 0.4
   !(small value causes problems in aerodynamic resistance calculation)
   if (modelInput%VZ(1) < 0.4) Then
      modelInput%VZ(1) = 0.4
   end if

   !First interpolated values
   atm%Tair = modelInput%Tair(1)
   atm%VZ = modelInput%VZ(1)
   atm%Rhz = modelInput%RHz(1)
   depth=modelInput%depth(1)
   !Calculate temperaturea at given depth
   if (depth>=0) Then
      Call getTempAtDepth(ground,depth,t_output)
      surf%TsurfAve=t_output
   else
      surf%TSurfAve = (ground%Tmp(1) + ground%Tmp(2))/2.0 !Average temperature of 
   end if
   !Calculate boundary layer conductance and latent heat flux
   call CalcBLCondAndLE(surf%TSurfAve, surf%EvapmmTS, settings%DtSecs, surf%SrfWatmms,&
                        phy, atm)

   !Don't use coupling if not surface temperature observation
   if (coupling%lastTsurfObs < -100) Then
      coupling%Coupling_failed = .true.
   end if


end subroutine

!>Allocates model output and input arrays
Subroutine allocator(settings,ground,coupling)
   use RoadSurfVariables
   Implicit None

   type(modelSettings), intent(INOUT) :: settings   !< Variables for model settings
   type(groundVariables), intent(OUT) :: ground     !< Variables for ground
                                                    !< properties
   type(couplingVariables), intent(OUT) :: coupling     !< variables used in
                                                        !< coupling(adjusting
                                                        !< radiation to fit
                                                        !< observed surface
                                                        !< temperature)

   
   allocate (ground%condDZ(settings%NLayers+1))
   allocate (ground%capDZ(settings%NLayers+1))
   allocate (ground%Wcont(settings%NLayers+1))
   allocate (ground%VSH(settings%NLayers+1))
   allocate (ground%HS(settings%NLayers+1))
   allocate (ground%CC(settings%NLayers+1))
   allocate (ground%Tmp(0:settings%NLayers+1))
   allocate (ground%TmpNw(0:settings%NLayers+1))
   allocate (ground%DyC(settings%NLayers+1))
   allocate (ground%DyK(settings%NLayers+1))
   allocate (ground%ZDpth(settings%NLayers+1))
   allocate (ground%GCond(0:settings%NLayers+1))
   allocate (coupling%TmpSave(0:settings%NLayers+1))

end Subroutine

!>Initialization of ground heat capacity and conductivity
subroutine ground_prop_init(NLayers,ground)
   use RoadSurfVariables
   implicit none
   Integer, intent(IN) :: Nlayers           !< Number of ground layers
   type(groundVariables), intent(INOUT) :: ground     !< Variables for ground
                                                    !< properties

   Integer :: j, i
   Integer::ConductivityLayers              !<Number of layers to calcuate
                                            !< conductivity

   !Calculate depth difference between layers (mid point - mid point)
   ground%DyC(1) = (ground%ZDpth(2) - ground%ZDpth(1))/2.0
   DO j = 2, Nlayers
         ground%DyC(j) = (ground%ZDpth(j + 1) - ground%ZDpth(j - 1))/2.0
   END DO

   ConductivityLayers = NLayers

   !Calculate the thicknesses of individual layers
   Do j = 1, ConductivityLayers

      ground%DyK(j) = ground%ZDpth(j + 1) - ground%ZDpth(j)

   END DO
   !Water content
   ground%WCont(1) = 0.01
   ground%WCont(2) = 0.01
   ground%WCont(3) = 0.3
   ground%WCont(4) = 0.3
   Do I = 5, NLayers
      ground%WCont(I) = 0.3
   end do
End Subroutine ground_prop_init

!>Initialize depth array
Subroutine initDepth(ground, NLayers)
   use RoadSurfVariables
   implicit none

   integer, intent(IN)::NLayers             !< Number of ground layers
   type(groundVariables), intent(INOUT) :: ground     !< Variables for ground
                                                    !< properties
   real(8):: ZAdd
   integer ::I

   ZAdd = 0.02
   ground%ZDpth(1) = 0.0
   !Calculates layer depths so that the thicknesses increase with depth

   Do I = 1, NLayers
      ground%ZDpth(I + 1) = ground%ZDpth(I) + 0.0103*1.4**(I - 1) + ZAdd
   End do

End Subroutine

!>Initialize temperature array
Subroutine initTemp(Tsurf, Tair,depth, phy, settings, ground, modelInput, surf)
   use RoadSurfVariables
   implicit none
   type(PhysicalParameters), intent(IN) :: phy    !< Physical paremeters used in the model
   type(ModelSettings), intent(IN) :: settings    !< Variables for model settings
   Real(8), intent(IN) :: Tsurf                   !< Surface temperature
   Real(8), intent(IN) :: Tair                    !< Air temperature
   Real(8), intent(IN) :: depth                   !> Depth to calculate tsurfave
   type(inputArrays), intent(IN) :: modelInput    !< Arrays for model input data
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   type(SurfaceVariables), intent(INOUT) :: surf  !< Variables for surface properties

   integer :: i, juld
   real(8) :: t_output

   ground%Tmp(0) = Tair !0th index is air temperature
   !First four are same as observed surface temperature, if observation
   !is available
   if (Tsurf > -100) Then
      Do i = 1, 4
         ground%Tmp(i) = Tsurf
      end do
   else
      Do i = 1, 4
         ground%Tmp(i) = Tair
      end do
   end if

   Call JulDay(modelInput, juld, 1)
   ground%Tmp(settings%NLayers + 1) = phy%TClimG + phy%AZ*Sin(phy%Omega*juld + &
      phy%Omega*(-170) -(ground%ZDpth(settings%NLayers+ 1)/phy%DampDpth))


   !Temperature of the rest of the layers approaches linearly to climatological value
   Do i = 5, settings%NLayers
      ground%Tmp(i) = ground%Tmp(4) + (ground%Tmp(settings%NLayers + 1) -&
       ground%Tmp(4))/ (ground%ZDpth(settings%NLayers + 1) - &
       ground%ZDpth(4))*(ground%ZDpth(i) - ground%ZDpth(4))
   End Do

   Do i = 0, settings%NLayers + 1
      ground%TmpNW(i) = ground%Tmp(i)
   end Do
   if (depth>=0) Then
      Call getTempAtDepth(ground,depth,t_output)
      surf%TsurfAve=t_output
   else
     surf%TsurfAve = 0.5*(ground%Tmp(1) + ground%Tmp(2))
   end if
End Subroutine

!> Initialize surface values
Subroutine initSurf(surf, wearOn)
   use RoadSurfVariables
   Implicit none
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
   logical, intent(IN) :: wearOn

   surf%Q2Melt = 0.0
   surf%VeryCold = .false.
   surf%WearSurf = wearOn
   surf%TrfFric = 5.0
   surf%EvapmmTS = 0.0
   surf%TSurfObs = -99.9
   surf%SrfWatmms = 0.0
   surf%SrfSnowmms = 0.0
   surf%SrfIcemms= 0.0
   surf%SrfIce2mms = 0.0
   surf%SrfDepmms = 0.0
   
End subroutine
!>Initialize parameters
Subroutine InitParam(Albedo, phy,inputParam)
   use RoadSurfVariables
   Implicit none

   type(InputParameters), intent(IN) :: inputParam          !< input parameters
                                                            !< given by modelRunner.cpp
   type(PhysicalParameters), intent(INOUT) :: phy !< Physical paremeters used in
                                                  !< the model
   real(8), intent(OUT) :: Albedo                  !< Surface albedo
   

   phy%Grav = inputParam%Grav            !< Gravitational acceleration (m/s2)
   phy%SB_Const = inputParam%SB_Const    !< Stefan-boltzman constant (W/m2K4)
   phy%VK_Const = inputParam%VK_Const    !< Von Karman's konstant
   phy%ZRefW = inputParam%ZRefW          !< Wind reference height (m)
   phy%ZRefT = inputParam%ZRefT          !< Wind reference height (m)
   phy%ZeroDisp = inputParam%ZeroDisp    !< Zero displacement height (m)
   phy%ZMom = inputParam%ZMom            !< Roughness factor for momentum (m)
   phy%ZHeat = inputParam%ZHeat          !< Roughness factor for heat (m)
   !< used to calcualte aerodynamic resistance
   phy%logMom=Log((phy%ZRefW + phy%ZMom)/phy%ZMom)     
   phy%logHeat =Log((phy%ZRefW + phy%ZHeat)/phy%ZHeat)
   !< used to calculate boundary layer conductance
   phy%logCond=Log((phy%ZRefW - phy%ZeroDisp +&
               phy%ZHeat)/phy%ZHeat)
   !< used to calculate friction velocity
   phy%logUstar=Log((phy%ZRefW - phy%ZeroDisp + phy%ZMom)/&
               phy%ZMom)
   phy%Emiss = inputParam%Emiss          !< Emissivity constant of the surface
   Albedo = inputParam%Albedo            !< Dry ground albedo
   phy%MaxPormms = inputParam%MaxPormms

   phy%TClimG = inputParam%TClimG        !< Climatological temperature at the bottom layer
   phy%DampDpth = inputParam%DampDpth    !< Damping depth
   phy%Omega = inputParam%Omega          !<constant to calculate bottom
                                         !< layer temperature
   phy%AZ = inputParam%AZ                !<constant to calculate bottom layer temperature
   phy%LVap = inputParam%LVap            !< Latent heat of water vaporisation (J/kg)
   phy%LFus = inputParam%LFus            !< Latent heat of fusion (constant, not calculated)
   
   phy%vsh1 = inputParam%vsh1            !< Heat capacity for surface layers
   phy%vsh2 = inputParam%vsh2            !< Heat capacity for deep ground layers
   phy%Poro1 = inputParam%Poro1          !< Porosity for surface layers
   phy%Poro2 = inputParam%Poro2          !< Porosity for deep ground layers
   phy%RhoB1 = inputParam%RhoB1          !< Bulk density for surface layers
   phy%RhoB2 = inputParam%RhoB2          !< Bulk density for deep ground layers
   phy%Silt1 = inputParam%Silt1          !< Clay fraction for surface layers
   phy%Silt2 = inputParam%Silt2          !< Clay fraction for deep ground layers
End Subroutine

!> Initializes atmoshpere and surface variables
Subroutine initVariables(ground, atm,settings)
   use RoadSurfVariables
   Implicit none
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties
   type(ModelSettings), intent(OUT) :: settings         !< Variables for model
                                                              !< settings

   integer :: i

   atm%Tair = -99.9
   atm%VZ = -99.9
   atm%Rhz = -99.9
   atm%PrecInTStep = -99.9
   atm%BLCond = -99.9
   atm%TairInitEnd = -99.9
   atm%VZInitEnd = -99.9
   atm%RhzInitEnd = -99.9
   atm%RainmmTS = 0.0
   atm%SnowmmTS = 0.0
   atm%SnowType = SURFACE_SNOW_DRY
   atm%CalmLim = 0.4
   atm%SensibleHeatFlux = -9999.9
   ground%GroundFlux = -9999.9
   Do I = 1, settings%NLayers
      ground%VSH(i) = -99.9
      ground%HS(i) = -99.9
      ground%CC(i) = -99.9
      ground%Gcond(i) = -99.9
   End Do
   ground%Gcond(0) = -99.9

End Subroutine

!> Initializes output arrays with missing values
Subroutine initOutputArrays(SimLen, modelOutput)
   use RoadSurfVariables
   implicit None
   integer, intent(IN):: SimLen                     !< length of output data arrays
   type(OutputArrays), intent(INOUT) :: modelOutput !< Arrays for model output data

   integer ::i
   Do i = 1, SimLen
      modelOutput%SnowOut(i) = -9999.0
      modelOutput%WaterOut(i) = -9999.0
      modelOutput%IceOut(i) = -9999.0
      modelOutput%Ice2Out(i) = -9999.0
      modelOutput%DepositOut(i) = -9999.0
      modelOutput%TsurfOut(i) = -9999.0
   end do
end Subroutine


!> Initialize car observation array for multicopuling
!> Multicoupling is not used in current model version
subroutine initTsurfObsArrays(coupling)
   use RoadSurfVariables
   Implicit None
   integer :: i
   type(CouplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)


!initialize observations used in coupling
   Do i = 1, 48
      coupling%carObsTime(i, 1) = -99
      coupling%carObsTime(i, 2) = -99
      coupling%carObsTime(i, 3) = -99
      coupling%carObsTime(i, 4) = -99
      coupling%carObsTime(i, 5) = -99
      coupling%carObsTime(i, 6) = -99
      coupling%obsI(i) = -99
      coupling%obsTsurf(i) = -99.0
   end do

end subroutine initTsurfObsArrays

!> Initializes model setting values
subroutine initSettings(inSettings, settings,inputParam,localParam)

   use RoadSurfVariables
   Implicit None

   type(InputSettings), intent(IN) :: inSettings !< model settings given
                                                         !< as input by modelRunner.cpp
   type(InputParameters), intent(IN) :: inputParam      !< input parameters
                                                        !< given by modelRunner.cpp
   type(LocalParameters), intent(IN) :: localParam  !< Local parameters
                                                         !< given by modelrunner.cpp
   type(ModelSettings), intent(OUT) :: settings !< Variables for model settings
   logical :: int2Logical
   settings%SimLen = inSettings%SimLen
   settings%InitLenI = localParam%InitLenI
   settings%DTSecs = inSettings%DTSecs
   settings%tsurfOutputDepth = inSettings%tsurfOutputDepth
   settings%NLayers=inSettings%NLayers

   settings%NightOn = inputParam%NightOn 
   settings%NightOff =  inputParam%NightOff
   settings%CalmLimDay =  inputParam%CalmLimDay
   settings%CalmLimNgt =  inputParam%CalmLimNgt

   settings%TrfFricNgt =  inputParam%TrfFricNgt
   settings%TrFfricDay =  inputParam%TrFfricDay

   settings%use_coupling = int2Logical(inSettings%use_coupling)
   settings%use_relaxation = int2Logical(inSettings%use_relaxation)
   settings%coupling_minutes=inSettings%coupling_minutes
   settings%couplingEffectReduction=inSettings%couplingEffectReduction
   settings%outputStep=inSettings%outputStep
   settings%force_tsurf=int2Logical(inSettings%force_tsurf)
 
end subroutine initSettings

!> Initializes parameters used in road condition calculation
Subroutine condInit(condParam, surf,inputParam)
   use RoadSurfVariables
   Implicit None

   type(InputParameters), intent(IN) :: inputParam      !< input parameters
                                                        !< given by modelRunner.cpp
   type(SurfaceVariables), intent(INOUT) :: surf        !< Variables for surface
                                                        !< properties
   type(RoadCondParameters), intent(OUT) :: condParam   !< Parameters to
                                                        !< determine storage
                                                        !< terms and road condition
   ! * Wet snow density/thickness same as that of dry snow
! * All thicknesses in mm equivalent water
   condParam%WatDens = inputParam%WatDens       !< Density of water at 0 C
   condParam%SnowDens = inputParam%SnowDens     !< Density of snow (fresh,!< Oke p.44)
   condParam%IceDens = inputParam%IceDens       !< Density of ice (Oke p.44)
   condParam%DepDens = inputParam%DepDens       !< Density of deposit
   condParam%WatMHeat = inputParam%WatMHeat     !< Heat of ablation for water (J/kg)
   condParam%PorEvaF = inputParam%PorEvaF       !< Pore resistance factor for evaporation
   condParam%DampWearF = inputParam%DampWearF   !< Damp surface (poer) wear factor

   !< Freezing limit (C)
   condParam%freezing_limit_normal =inputParam%freezing_limit_normal  
   !< Melting limit for snow (C)
   condParam%snow_melting_limit_normal =inputParam%snow_melting_limit_normal
   !< Melting limit for ice (C)
   condParam%ice_melting_limit_normal = inputParam%ice_melting_limit_normal
   !< Melting limit for deposit (C)
   condParam%frost_melting_limit_normal = inputParam%frost_melting_limit_normal                                                    
   !< Dew/deposit formation limit	
   condParam%frost_formation_limit_normal =inputParam%frost_formation_limit_normal 
                                                        
   condParam%T4Melt_normal = inputParam%T4Melt_normal
  
   condParam%TLimFreeze = inputParam%freezing_limit_normal
   condParam%TLimMeltSnow = inputParam%snow_melting_limit_normal
   condParam%TLimMeltIce = inputParam%ice_melting_limit_normal
   condParam%TLimMeltDep = inputParam%frost_melting_limit_normal
   condParam%TLimDew = inputParam%frost_formation_limit_normal
   surf%T4Melt = inputParam%T4Melt_normal
  
   condParam%TLimColdH =inputParam%TLimColdH        !< Higher limit for cold ground temperature
   condParam%TLimColdL = inputParam%TLimColdL       !< Lower limit for cold ground temperature
   condParam%WetSnowFormR = inputParam%WetSnowFormR !< Water to snow ratio for wet snow formation
   condParam%WetSnowMeltR =inputParam%WetSnowMeltR  !< Water to snow ratio for wet snow melting
! * Rain type parameters
   condParam%PLimSnow =inputParam%PLimSnow          !< Precipitation interpretation : snow limit
   condParam%PLimRain = inputParam%PLimRain         !< Precipitation interpretation : rain limit
   condParam%MinPrecmm = inputParam%MinPrecmm    !< MIN precipitation  mm/hour
   condParam%MinWatmms = inputParam%MinWatmms    !< MIN water storage : mm
   condParam%MinSnowmms =inputParam%MinSnowmms   !< MIN snow storage : mm
                                                 !<    (accounts for snow "drift")
 
   condParam%MinDepmms = inputParam%MinDepmms    !< MIN deposit storage : mm
   condParam%MinIcemms = inputParam%MinIcemms    !< MIN ice storage : mm (MinPrec => condens
                                                 !< won't freeze(?))

   condParam%MaxSnowmms =inputParam%MaxSnowmms   !< MAX snow storage : mm (plowed away if above)
   condParam%MaxDepmms = inputParam%MaxDepmms    !< MAX deposit storage : mm
   condParam%MaxIcemms = inputParam%MaxIcemms    !< MAX ice storage : mm
   condParam%MaxExtmms = inputParam%MaxExtmms    !< MAX extra (surface)
                                                 !< water content (mm) 

   condParam%MaxWatmms =inputParam%MaxWatmms
   condParam%AlbDry = inputParam%AlbDry
   condParam%AlbSnow = inputParam%AlbSnow
   condParam%MissValI = inputParam%MissValI
   condParam%MissValR = inputParam%MissValR
   condParam%WDampLim = inputParam%WDampLim    !< Dry/Damp limit for water
   condParam%WWetLim = inputParam%WWetLim      !< Damp/Wet limit for water Ver.6.23: (0.5 => 0.9)
   condParam%WWearLim = inputParam%WWearLim    !< Wear limit for water (below this only
                                               !< evaporation)

   condParam%Snow2IceFac =inputParam%Snow2IceFac   !< Snow to ice transition factor
   condParam%WetSnowFrozen=.false.
   condParam%forceIceMelting=.false.
   condParam%forceSnowMelting=.false.
   condPAram%CanMeltingChangeTemperature=.true.
end subroutine condInit
!> Turn integer to locigal
!> Used because its easier to paas integers from c++ to frotran
Logical function int2Logical(intVal)
   implicit none
   integer :: intVal

   if (intVal == 1) Then
      int2Logical = .true.
   else
      int2Logical = .false.
   end if
   return
end function

