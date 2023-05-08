
!>Calculate boundary layer conductance
Subroutine CalcBLCondAndLE(TSurfAve, EvapmmTS, DtSecs, SrfWatmms, phy, atm)
   use RoadSurfVariables
   Implicit none
   real(8), intent(IN) :: TSurfAve              !< average temperature of the
                                                !< first two layers
   real(8), intent(IN) :: DtSecs                !< time steps in seconds
   real(8), intent(IN) :: SrfWatmms             !< water storage

   type(PhysicalParameters), intent(IN) :: phy  !< Physical paremeters used in
                                                !< the model
   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties

   real(8), intent(INOUT) :: EvapmmTS           !< Evaporation

   Real(8), Parameter    :: ConvLim = 0.001     !< Convergence limit for BLCond
                                                !< (absolute diff)
   Integer, Parameter :: MaxIter = 40           !< Max number of iterations for
                                                !< BLCond (40)
   Real(8) :: PSIM       !< Stability factor, momentum
   Real(8) :: PSIH       !< Stability factor, heat
   Real(8) :: UStar      !< friction velocity
   Real(8) :: Stab       !< Stability parameter 
   Real(8) :: RAero      !< Aerodynamic resistance
   Real(8) :: Tair       !< Air temperature
   Real(8) :: VZ         !< wind speed
   Real(8) :: Rhz        !< relative humidity
   Real(8) :: BLCond     !< Boundary layer conductance
   Integer :: j          !< Index for loops
   Real(8) :: TaK        !< Air temperature in Kelvins
   Real(8) :: AirDens    !< Air density
   Real(8) :: AirHCap    !< Air specific heat 
   Real(8) :: AirVCap    !< Ari volumetric spedific heat
   Real(8) :: PsychC     !< Psychrometric constant [kPa/K]
   Real(8) ::  WatDen    !< Water density
   Real(8) :: BLCond_old !< previous iteration Boundary layer conductance

   !Campbell, Gaylon S., Soil Physics wit Basic. Elsevier, 
   !the Netherlands, 1985
   !Campbell, Gaylon S, An introduction to Environmental Biophysics,
   !Springer Verlag, New York 1986, Appendix

   Tair = atm%Tair
   VZ = atm%VZ
   Rhz = atm%Rhz
   BLCond = atm%BLcond

   TaK = Tair + 273.15                          !< Air temp C => K
   !* Temperature dependent
   AirDens = 100000.0/(287.05*TaK)              !< Air density [kg/m3]
   AirHCap = 1005.0+((TaK - 250.0)**2)/3364.    !< Air specific heat [J/kgK]
   AirVCap = AirHCap*AirDens                    !< Air volumetric specific heat
                                                !< [J/kgK]
   PsychC = 0.1*(0.00063*TaK + 0.47496)         !< Psychrometric constant [kPa/K]
   WatDen = -0.0050*TSurfAve*TSurfAve + 0.0079*TSurfAve + 1000.0028 !< Water
                                                                    !< density[kg/ms]

!********** CALCULATE BOUNDARY LAYER CONDUCTANCE (BLCond)
   PSIM = 0.0 ! STAB. CORR. FACTOR FOR MOMENTUM (INITIAL VALUE)
   PSIH = 0.0 ! STAB. CORR. FACTOR FOR HEAT (INITIAL VALUE)
   Stab= 0.0
   Do j = 1, MaxIter

!       SAVE OLD VALUE
      BLCond_Old = BLCond

!       FRICTION VELOCITY
      UStar = phy%VK_Const*VZ/(phy%logUstar + PSIM)
      If (UStar < 0.0) Then
         Write (*, *) ' ERROR : UStar negative,vz ', vz
         write (*, *) Tair, VZ, Rhz, BLCond, TSurfAve
      End If

!       BOUNDARY LAYER CONDUCTANCE / STABILITY PARAMETER
      BLCond = AirVCap*phy%VK_Const*UStar/(phy%logCond + PSIH)
      Stab = -phy%VK_Const*phy%ZRefT*phy%Grav*BLCond*(TSurfAve - Tair)/&
         (AirVCap*(Tair +273.15)*(UStar*UStar*UStar))
      IF (Stab > 1) Stab = 1

!       STABILITY CORRECTION FACTORS
      If (Stab .GT. 0) Then ! STABLE CONDITION
         PSIH = 4.7*Stab
         PSIM = PSIH
      Else ! UNSTABLE CONDITION
         PSIH = -2.0*Log((1.0+Sqrt(1.0-16.0*Stab))/2.0)
         PSIM = 0.6*PSIH
      End If

!       CHECK ITERATION CONVERGENCE
      If ((abs(BLCond - BLCond_Old) < ConvLim) .and. (j >= 5)) Then ! OK; stop
                                                                    !< iteration
         EXIT
      end if
   end do

   If ((abs(BLCond - BLCond_Old) > 10*ConvLim) .and. (j >= 5)) Then
      Write (*, "(' Max number of BLCond iterations (MaxIter,BLCond_Old,BLCond) :',I5,2F10.5)")&
             j, BLCond_Old, BLCond
   end if

   call calcRaero(phy%logMom,phy%logHeat, PSim, PSIH, phy%VK_Const, VZ,&
                  Raero)
   call CalcLE(TSurfAve, Tair, RHz, AirDens, AirHCap, PsychC, RAero, phy%LVap,&
               phy%LFus, WatDen, DtSecs, SrfWatmms, atm%LE_Flux, EvapmmTS)

   atm%BLcond = BLCond
End Subroutine

!> Calculates aerodynamic resistance
Subroutine calcRaero(logMom,logHeat, PSIM, PSIH, VK_Const, VZ, Raero)
   implicit none
   Real(8), intent(IN)      :: PSIM    !< Stability parameter for momentum
   Real(8), intent(IN)      :: PSIH    !< Stability parameter for heat
   Real(8), intent(IN)      :: VZ      !< Wind speed
   Real(8), intent(IN)      :: VK_Const !< Von Karman's konstant
   Real(8), intent(IN)      :: logMom  !< log((ZRefW+ZMom)/ZMom)
   Real(8), intent(IN)      :: logHeat !<log((ZRefW+ZHeat)/ZHeat)
   Real(8), intent(OUT) :: RAero       !< Aerodynamic resistance

   !Tourula T., M. Heikinheimo, 1998: Modelling evapotranspiration from a barley
   !field over the growing season, Agricultural and Forest Meteorology 91, p.237-250

   RAero = (logMom + PSim) &
           *(logHeat + PSIH)/(VK_Const*VK_Const*VZ)
   if (RAero > 30.0) Then
      RAero = 30.
   End If

End subroutine

!> Calculates latent heat flux and evaporation
Subroutine CalcLE(TSurfAve, TAmb, Rhz, AirDens, AirHCap, PsychC, RAero, LVap, LFus,&
                  WatDen, DtSecs, SrfWatmms, LE_Flux, EvapmmTS)
   implicit none

   real(8), intent(IN) :: TSurfAve     !< average temperature of the first two layers
   real(8), intent(IN) :: TAmb         !< Air temperature (C)
   real(8), intent(IN) :: Rhz          !< Relative humidity (%)
   real(8), intent(IN) :: AirDens      !< Air density
   real(8), intent(IN) :: AirHCap      !< Air heat capacity
   real(8), intent(IN) :: PsychC       !< Psychrometric constant [kPa/K]
   real(8), intent(IN) :: RAero        !< Aerodynamic resistance
   real(8), intent(IN) :: LVap         !< latent heat of water vaporization
   real(8), intent(IN) :: LFus         !< latent heat of sublimation
   real(8), intent(IN) :: WatDen       !< Density of water
   real(8), intent(IN) :: DtSecs       !< time step in seconds
   real(8), intent(IN) :: SrfWatmms    !< water storage
   real(8), intent(INOUT) :: LE_Flux   !< latent heat flux
   real(8), intent(INOUT) :: EvapmmTS  !< Evaporation
   real(8) :: ESat  !< starutated water wapor pressure
   real(8) :: ESurf !< water wapor pressure over surface
   real(8) :: EAir  !< water wapor pressure air

   !Calder, I.R., 1990: Evaporation in teh Uplands. Jhon Wiley & sons

   ! Water wapor pressure over surface
   If (TSurfAve < 0) Then
      ESat = 0.61078*Exp(21.875*TSurfAve/(TSurfAve + 265.5)) ! over ice
   Else
      ESat = 0.61078*Exp(17.269*TSurfAve/(TSurfAve + 237.3)) ! over water
   End If
   ESurf = ESat

   ! * Air water vapor pressures
   If (TAmb < 0) Then
      ESat = 0.61078*Exp(21.875*TAmb/(TAmb + 265.5)) ! over ice
   Else
      ESat = 0.61078*Exp(17.269*TAmb/(TAmb + 237.3)) ! over water
   End If
   EAir = Min((0.01*Rhz), 1.0)*ESat

   !********** LATENT HEAT FLUX (W/m2) AND EVAPORATION (mm/TimeStep)
   !** Calder : Evaporation in the Uplands p.17-18
   LE_Flux = (AirDens*AirHCap*(ESurf - EAir))/(PsychC*RAero)

   If (TSurfAve >= 0.0) Then
      EvapmmTS = (LE_Flux/(LVap*WatDen))*1000.0*DTSecs ! water
   Else
      EvapmmTS = (LE_Flux/(LFus*WatDen))*1000.0*DTSecs ! frost
   End If

   !********** CHECK WATER AVAILABILITY

   If ((LE_Flux > 0.0) .and. (SrfWatmms <= 0.0)) Then ! No water to evaporate
      LE_Flux = 0.0
      EvapmmTS = 0.0
   End If
end subroutine

