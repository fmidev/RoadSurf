!MIT License
!Copyright (c) 2023 FMI Open Development
Submodule (RoadSurf) BalanceModel
   Implicit None
   contains
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
         real(8) ::depth                                  !< depth to take temperature
         real(8) :: t_output                              !< temperature at depth
         !Set traffic friction and check that wind is not below minimum
         call SetDayDependendVariables(settings, surf, modelInput, atm, inputIdx)
      
         !Calculate boundary layer conductance
         call CalcBLCondAndLE(surf%TSurfAve, surf%EvapmmTS, settings%DtSecs, &
                              surf%SrfWatmms, phy, atm)
      
         !Calculate net radiation
         Call CalcRNet(phy%Emiss, phy%SB_Const, surf%TSurfAve, ground%Albedo, SWi, &
                       LWi, atm%RNet, coupling%SwRadCof, coupling%LwRadCof)
         !Calculate heat capacity and conductance
         call CalcHCapHCond(settings%NLayers, settings%DTSecs, phy, ground, atm)
      
         !Calculate heat capacity and conductance derivatives
         call calcCapDZCondDZ(settings%NLayers, ground)
      
         !Calculate temperature profile for the next time step
         call calcProfile(settings%Nlayers, settings%DTSecs, surf%TrfFric, &
                          ground, atm)
      
      
         !Calculate heat storage
         call calcHStor(ground)
      
         !Check if output depth is given in settings
         if (settings%tsurfOutputDepth>=0.0) Then
            depth=settings%tsurfOutputDepth
         else
            !If not, set to value given in input (-9999.9 if missing)
            depth=modelInput%depth(InputIdx)
         end if
      
         !> Calculate the heat needed to melt/freeze the whole uppermost snow/ice layer.
         call NewMeltFreezeHeat(settings%DTSecs, surf, condParam)
         !Check if melting
         call melting(coupling%inCouplingPhase, coupling%lastTsurfObs, &
                      ground, depth,surf, &
                      condParam)
      
         ground%Tmp = ground%TmpNw
        
         !If depth is set, interpolate temperature from given depth
         if (depth>=0)Then
            Call getTempAtDepth(ground,depth,t_output)
            surf%TsurfAve=t_output
         else
         !Otherwise use average of the first two layers
            surf%TSurfAve = (ground%Tmp(1) + ground%Tmp(2))/2.0
         end if
      
      End Subroutine
end submodule

!>Calculates temperature profile one time step forward
Subroutine calcProfile(Nlayers, DTSecs, TrfFric, ground, atm)
   use RoadSurfVariables
   implicit none
   Integer, intent(IN) :: NLayers               !< Number of ground layers
   real(8), intent(IN) :: DTSecs                !< time step in seconds
   real(8), intent(IN) :: TrfFric               !< Surface heating caused by traffic

   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   integer :: j
   Real(8),allocatable :: Gflux(:)                 !< Heat flux
   Integer :: NLayerTemp                        !< Number of layers to caclulate
                                                !< temperature
   !Brutsaert, W, Evaporation into the Atmosphere. Theory, History, and
   !Applications. D. Reidel Publishing Company, Dordrecht, Holland, 1984
   
   !Campbell, Gaylon S., Soil Physics wit Basic. Elsevier, 
   !the Netherlands, 1985
   allocate (Gflux(0:NLayers+1))
   NLayerTemp = NLayers

   atm%SensibleHeatFlux = atm%BlCond*(ground%Tmp(0) - ground%Tmp(1))

   !Heat flux from air to ground
   GFlux(0) = atm%RNet - atm%LE_Flux + TrfFric + atm%SensibleHeatFlux
   ground%TmpNw = ground%Tmp

   !Calculate heat flux for different layers
   do j = 1, NLayerTemp
      GFlux(j) = ground%condDZ(j)*(ground%Tmp(j + 1) - ground%Tmp(j))
   end do
   ground%GroundFlux = GFlux(3)

   !Calculate new temperatures
   do j = 1, NLayers
      ground%TmpNw(j) = ground%Tmp(j) + DTsecs*(ground%capDZ(j)*(GFlux(j) -&
       Gflux(j - 1)))
   end do
end subroutine calcProfile

!>Calculate heat capacity and conductance derivatives
Subroutine calcCapDZCondDZ(NLayers, ground)
   use RoadSurfVariables
   implicit none
   Integer, intent(IN) :: NLayers               !< number of ground layers
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties

   Integer :: NLayerTemp                        !< Number of layers to caclulate
                                                !< temperature
   integer:: j                                  !< loop integer

   NLayerTemp = NLayers

   !Calculate help variables to use in the temperature profile calcualtion
   ground%condDZ(1) = -(ground%CC(1)/ground%DYK(1))
   ground%capDZ(1) = -(1/(ground%DYC(1)*ground%VSH(1)))

   DO j = 2, NLayerTemp

      ground%condDZ(j) = -(ground%CC(j)/ground%DYK(j))
      ground%capDZ(j) = -(1/(ground%DYC(j)*ground%VSH(j)))

   END DO

End Subroutine

!> Determine some values needed in heat capacity and conductivity calculation
Subroutine HCapValues(phy)
   use RoadSurfVariables
   Implicit None
   type(PhysicalParameters), intent(INOUT) :: phy !< Physical paremeters used in
                                                  !< the model 
   !Campbell, Gaylon S., Soil Physics wit Basic. Elsevier, 
   !the Netherlands, 1985, p.32

   phy%Afc1 = 0.65-0.78*phy%RhoB1 + 0.60*phy%RhoB1*phy%RhoB1
   phy%Bfc1 = 1.06*phy%RhoB1
   If (phy%Silt1 > 0.00001) Then
      phy%Cfc1 = 1 + 2.6/sqrt(phy%Silt1)
   Else
      phy%Cfc1 = 0.
   End If
   phy%Dfc1 = 0.03+0.1*phy%RhoB1*phy%RhoB1
   phy%Efc1 = 4

   phy%Afc2 = 0.65-0.78*phy%RhoB2 + 0.60*phy%RhoB2*phy%RhoB2
   phy%Bfc2 = 1.06*phy%RhoB2
   If (phy%Silt2 > 0.00001) Then
      phy%Cfc2 = 1 + 2.6/sqrt(phy%Silt2)
   Else
      phy%Cfc2 = 0.
   End If
   phy%Dfc2 = 0.03+0.1*phy%RhoB2*phy%RhoB2
   phy%Efc2 = 4

End Subroutine

!>Calculates Ground layers heat capacity and heat conductance
Subroutine CalcHCapHCond(NLayers, DTSecs, phy, ground, atm)
   use RoadSurfVariables
   implicit none
   integer, intent(IN) ::NLayers                !< number of ground layers
   real(8), intent(IN)::DTsecs                  !< time step in seconds
   type(PhysicalParameters), intent(IN) :: phy  !< Physical paremeters used in
                                                !< the model
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   type(AtmVariables), intent(IN) :: atm        !< Variables for atmospheric
                                                !< properties
   integer :: i      !< loop integer
   real(8) :: RooWT  !< water density (kg/m3)
   real(8) :: CWT    !< specific heat capacity of water (kJ/kgK)
   real(8) :: CHWT   !< Volumetric heat capacity for water
   real(8) :: tmp2   !< surface temperature power 2 

   ! ************* HEAT CAPACITY AND HEAT CONDUCTIVITY BY LAYERS
   !Campbell, Gaylon S, An introduction to Environmental Biophysics,
   !Springer Verlag, New York 1986, Appendix
   !Oke, T.R. Boundary Layer Climates. Methuen, New Yourk 1987

   ground%GCond(0) = atm%BLCond ! Boundary layer conductance

   !Go trough ground layers
   Do I = 1, NLayers
      !if above feezing
      If (ground%TmpNw(I) >= 0) Then ! Water
         tmp2=ground%TmpNw(I)*ground%TmpNw(I)
         !< water density (kg/m3)
         RooWT = -0.0050*tmp2 + &
                 0.0079*ground%TmpNw(I) + 1000.0028 
         !specific heat capacity of water (kJ/kgK)
         CWT = 0.0000102*tmp2*tmp2- &
               0.0017169*tmp2*ground%TmpNw(I)+ &
               0.11516*tmp2 &
               - 3.4739*ground%TmpNw(I) + 4217.2
      Else ! Ice (Oke p.44)
         RooWT = 920.0
         CWT = 2100.0
      End If
      CHWT = RooWT*CWT !Volumetric heat capacity for water
      !Volumetric heat capacity for moist ground is calculate as 
      !weihghted average of the dry ground and water
      If (I <= 2) Then
         ground%VSH(I) = (1.0-phy%Poro1)*phy%VSH1 + ground%WCont(I)*CHWT
      Else
         ground%VSH(I) = (1.0-phy%Poro2)*phy%VSH2 + ground%WCont(I)*CHWT
      End If
      !Volumetric heat capacity (J/m3K) is converted to intesity units (W/m2K)
      If (I == 1) Then
         ! Surface node at ground surface
         ground%HS(I) = ground%VSH(I)*(ground%ZDpth(I + 1) - &
              ground%ZDpth(I))/(2.0*DTSecs)
          
      Else !
         ground%HS(I) = ground%VSH(I)*(ground%ZDpth(I + 1) - ground%ZDpth(I -&
          1))/(2.0*DTSecs)
      End If
      !Heat conductivity (W/mK) converted to (W/m2K)
      ground%GCond(I) = ground%CC(I)/(ground%ZDpth(I + 1) - ground%ZDpth(I))
   end do

end subroutine

!>Calculate heat conductivity
Subroutine CalcCC(NLayers, phy, ground)
   use RoadSurfVariables
   implicit none
   integer, intent(IN) ::NLayers                !< number of ground layers
   type(PhysicalParameters), intent(IN) :: phy  !< Physical paremeters used in
                                                !< the model
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   integer :: i

   ! ************* HEAT CAPACITY BY LAYERS
   !Campbell, Gaylon S., Soil Physics wit Basic. Elsevier, 
   !the Netherlands, 1985, p.32

   !Heat conductivity based on water content
   Do I = 1, NLayers

      If (I <= 2) Then
         ground%CC(I) = phy%Afc1 + phy%Bfc1*ground%WCont(I) - &
             (phy%Afc1 - phy%Dfc1)*exp(-(phy%Cfc1*ground%WCont(I))**phy%Efc1)
      Else
         ground%CC(I) = phy%Afc2 + phy%Bfc2*ground%WCont(I) - &
             (phy%Afc2 - phy%Dfc2)*exp(-(phy%Cfc2*ground%WCont(I))**phy%Efc2)
      End If
          
   end do
end subroutine

!>Calculate net radiation
Subroutine CalcRNet(Emiss, SB_Const, TSurfAve, Albedo, SW, LW, RNet, SwRadCof,&
                    LWRadCof)
   Implicit None
   real(8), intent(IN) :: Emiss    !< Emissivity constant of the surface
   real(8), intent(IN) :: SB_Const !< Stefan-boltzman constant (W/m2K4)
   real(8), intent(IN) :: TSurfAve !< average temperature of the first two layers
   real(8), intent(IN) :: Albedo   !< surface albedo
   real(8), intent(IN) :: SW       !< Downwelling short wave radiation
   real(8), intent(IN) :: LW       !< Downwelling long wave radiation
   real(8), intent(IN) :: LWRadCof !< radiation coefficient to use for long wave
                                   !< radiation
   real(8), intent(IN) :: SWRadCof !< radiation coefficient to use for short wave
                                   !< radiation
   Real(8), intent(OUT) :: RNet    !< Net radiation
   Real(8) :: RBB
   Real(8):: TsurfK,TsurfK2

   !Brutsaert, W, Evaporation into the Atmosphere. Theory, History, and
   !Applications. D. Reidel Publishing Company, Dordrecht, Holland, 1984

   TsurfK=TSurfAve+273.15
   TsurfK2=TsurfK*TsurfK
   RBB = Emiss*SB_Const*(TsurfK2*TsurfK2) ! Black body emission
   !SwRadCof and LWRadCof determined in coupling
   RNet = (1.-Albedo)*SW*SwRadCof + Emiss*LW*LWRadCof - RBB
End Subroutine

!>Calculate heat variable describing the stored heat
!!to the surface from previous time step
Subroutine calcHStor(ground)
   use RoadSurfVariables
   Implicit none

   type(GroundVariables), intent(INOUT) :: ground
   real(8) ::T1Ave, TN1Ave

   T1Ave = (ground%Tmp(1) + 3.*ground%Tmp(2))/4.
   TN1Ave = (ground%TmpNw(1) + 3.*ground%TmpNw(2))/4.

   ground%HStor = ground%HS(1)*(TN1Ave - T1Ave)
End Subroutine

!>Calculate julian day
Subroutine JulDay(modelInput, juday, i)
!---------------------------------------------------------------------------
!   * Calculate Julian day
!   * Variable leapcorr = 1 for leap year, 0 otherwise
!   * Note : Jan 1st is day 1, not day 0
   use RoadSurfVariables
   Implicit None

   type(inputArrays), intent(IN) :: modelInput  !< Model input arrays
   integer, intent(IN) :: i                     !< Index in input data
   integer, intent(OUT) ::juday                 !< Julian day

   Integer, Parameter, Dimension(24) :: &
      MonEnd = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, &
                 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/) 
                ! first array for normal year, second for leap year

   Integer :: syear, smon, sday, leapcorr

   syear = modelInput%year(i)
   smon = modelInput%month(i)
   sday = modelInput%day(i)
   leapcorr = 1 - Min(Mod(syear, 4), 1) + Min(Mod(syear, 100), 1) - Min(Mod(syear,&
    400), 1)
   juday = MonEnd(smon + leapcorr*12) + sday

End Subroutine JulDay

!>Set traffic friction and check that wind is not below minimum
Subroutine SetDayDependendVariables(settings, surf, modelInput, atm, inputIdx)
   use RoadSurfVariables
   Implicit None

   integer, intent(IN) :: inputIdx              !< Input data index
   type(inputArrays), intent(IN) :: modelInput  !< Model input arrays
   type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties

   integer :: shour

   shour = modelInput%hour(inputIdx)

   If ((shour >= settings%NightOn) .or. (shour <= settings%NightOff)) Then
      ! Night time
!       * Emulate heat effect of traffic by wind induced turbulence and friction
!         - different coefficients for day and night
      atm%CalmLim = settings%CalmLimNgt
      surf%TrfFric = settings%TrfFricNgt

   Else
      ! Daytime
      atm%CalmLim = settings%CalmLimDay
      surf%TrfFric = settings%TrfFricDay
   End If

   !Minimum wind speed
   If (atm%VZ < atm%CalmLim) Then

      atm%VZ = atm%CalmLim
   end if
End subroutine SetDayDependendVariables

!Give interpolated temperature at certain depth
Subroutine getTempAtDepth(ground,depth,output_tsurf)
   use RoadSurfVariables
   Implicit None
   type(GroundVariables), intent(INOUT) :: ground   !< Ground related variables
   Real(8),intent(IN) :: depth                      !< depth for output temperature
   Real(8),intent(OUT) :: output_tsurf              !< output temperature
   Integer :: idx,zlen

   zlen=size(ground%ZDpth)
   !if depth is 0.0 return temperature of the first layer
   if (abs(depth-0.0)<0.00001) Then
      output_tsurf=ground%Tmp(1)
   !if depth is greatern than the depth of the last layer return
   !< the temperature of the last layer
   else if (depth>ground%ZDpth(zlen)) Then
      output_tsurf=ground%Tmp(zlen)
   else
     !Determine where the given depth is in relation to the depth array
     Do idx=1,zlen-1
        if (depth>ground%ZDpth(idx) .AND. depth <= ground%ZDpth(idx+1)) Then
           EXIT
        end if
     end do 
     !Interpoilate temperature to the given depth    
     output_tsurf=ground%Tmp(idx)+(depth-ground%ZDpth(idx))*&
                  (ground%Tmp(idx+1)-ground%Tmp(idx))/(ground%ZDpth(idx+1)-ground%ZDpth(idx))
   end if
End subroutine getTempAtDepth

