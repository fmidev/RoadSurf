
#include "Constants.h"

Submodule (RoadSurf) Coupling
   implicit none
   contains
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
         real(8) :: DTs                                          !< input time step in seconds
      
         DTs = settings%DTSecs
      
         coupling%inCouplingPhase = .false.
         !Determine if simulation in coupling phase
         if (i >= coupling%couplingStartI(coupling%CoupPhaseN) .and. &
             i <= coupling%couplingEndI(coupling%CoupPhaseN)) Then
      
            coupling%inCouplingPhase = .true.
      
         end if
         !At start of coupling period save simulated values so that
         !they can be restored later
         if (i == coupling%couplingStartI(coupling%CoupPhaseN) .and. &
             coupling%Coupling_iterations == 0) Then
            call saveDataForCoupling(i, ground%Albedo, ground%Tmp, coupling, surf, &
                                     modelInput)
            !Initialize coefficients
            coupling%SwRadCof = 1.0
            coupling%LWRadCof = 1.0
            coupling%SW_correction = 0.0
            coupling%LW_correction = 0.0
      
         end if
         !If coupling is started again, return to start of coupling period
         if (coupling%start_coupling_again) Then
            !Restore parameters at the beginning of coupling
            call uploadDataForCoupling(i, ground%Albedo, ground%Tmp, coupling, &
                                       surf, modelInput)
            coupling%start_coupling_again = .false.
            !Set coefficient for short wave radiation if it has larger value than
            !Long wave radiation or sky view factor is NOT used
            IF (modelInput%SW(i) > modelInput%LW(i).AND. &
             .NOT. (localParam%sky_view<1.0 .and. localParam%sky_view>-0.01)) Then
               coupling%SwRadCof = coupling%RadCoeff
               coupling%LWRadCof = 1.0
            Else
               !Othervise set coefficients for long wave radiation
               coupling%SwRadCof = 1.0
               coupling%LWRadCof = coupling%RadCoeff
            End if
      
         end if
      
         !Radiation coefficients after coupling
         !Returns gradually to 1
         if (i > coupling%couplingEndI(coupling%CoupPhaseN)) Then
      
            coupling%SwRadCof = 1.0+coupling%SW_correction*exp(-((DTs*i) - &
                          (DTs*coupling%couplingEndI(coupling%CoupPhaseN)))/settings%couplingEffectReduction)
            coupling%LWRadCof = 1.0+coupling%LW_correction*exp(-((DTs*i) - &
                          (DTs*coupling%couplingEndI(coupling%CoupPhaseN)))/settings%couplingEffectReduction)
         end if
         if (coupling%inCouplingPhase) Then
      
          !Check if there is a need to force snow/ice melting
          ! to prevent coupling from getting stuck
            call snowIceCheck(coupling%LastTsurfObs, surf, CP)
         end if
      
      end Subroutine
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
      
      
         !Coupling control if at the end of the coupling period
         if (settings%use_coupling .and. &
             i == coupling%CouplingEndI(coupling%CoupPhaseN) .and. &
             (coupling%coupling_failed .eqv. .false.)) Then
            !Determine new radiation coefficient if necessary
            call CouplingOperations2(surf%TSurfAve, coupling)
         end if
      
      end subroutine
end submodule

!> Determine new radiation coefficient if necessary
Subroutine CouplingOperations2(TSurfAve, coupling)
   use RoadSurfVariables
   Implicit none

   real(8), intent(INOUT) :: TsurfAve              !< average temperature of the
                                                !< first two layers
   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)

   !save temperature if at first iteration
   if (coupling%Coupling_iterations == 0) Then
      coupling%Tsurf_end_coup1 = TSurfAve
   end if
   !Change radiation coefficient if necessary
   call Coupling_control(TSurfAve, coupling)

   coupling%Coupling_iterations = coupling%Coupling_iterations + 1
end subroutine

!> Initialize coupling variables
Subroutine initCoupling(coupling)
   use RoadSurfVariables
   Implicit None

   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)
   coupling%Coupling_iterations = 0
   coupling%TsurfNearestAbove = -9999.0
   coupling%TsurfNearestBelow = -9999.0
   coupling%RadCoeff = 1.0
   coupling%RadCoefNearestAbove = -9999.0
   coupling%RadCoefNearestBelow = -9999.0
   coupling%RadCoeffPrevious = 1.0
   coupling%SwRadCof = 1.0
   coupling%LWRadCof = 1.0
   coupling%start_coupling_again = .false.
   coupling%Coupling_failed = .false.
   coupling%SW_correction = 0.0
   coupling%LW_correction = 0.0
   coupling%inCouplingPhase = .false.
   coupling%CoupPhaseN = 1
   coupling%LastTsurfObs = coupling%obsTsurf(1)

End Subroutine initCoupling

!>Save data at the beginning of the coupling period
Subroutine saveDataForCoupling(datai, Albedo, Tmp, coupling, surf,&
 modelInput)
   use RoadSurfVariables
   implicit none

   integer, intent(IN):: datai                  !<index for point in input data
   real(8), intent(IN) :: Albedo                !< Surface albedo
   real(8), dimension(0:16), intent(IN)::Tmp    !< Temperatures for each layer
   type(surfaceVariables), intent(IN) :: surf   !< Variables for surface properties
   type(inputArrays), intent(IN) :: modelInput  !< Arrays for model
   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)
   integer :: i
   integer :: couplingLen !Lenght of coupling period
   couplingLen=coupling%couplingEndI(1)-coupling%couplingStartI(1)+1

   coupling%saveDatai = datai
   coupling%TsurfAveSave = surf%TsurfAve
   coupling%srfWatmmsSave = surf%SrfWatmms
   coupling%srfIce2mmsSave = surf%SrfIce2mms
   coupling%srfIce2mmsSave = surf%SrfIce2mms
   coupling%srfDepmmsSave = surf%SrfDepmms
   coupling%srfSnowmmsSave = surf%SrfSnowmms
   coupling%AlbedoSave = Albedo
   coupling%VeryColdSave = surf%VeryCold
   do i = 0, 16
      coupling%TmpSave(i) = Tmp(i)
   end do

   do i=1,couplingLen
      coupling%SWSave(i)=modelInput%SW(coupling%couplingStartI(1)+i-1)
      coupling%SWDirSave(i)=modelInput%SW_dir(coupling%couplingStartI(1)+i-1)
      coupling%LWSave(i)=modelInput%LW(coupling%couplingStartI(1)+i-1)
   end do 

end subroutine saveDataForCoupling

!>Upload saved data if coupling is started again
Subroutine uploadDataForCoupling(datai, Albedo, Tmp, coupling, surf,&
 modelInput)
   use RoadSurfVariables
   implicit none

   integer, intent(INOUT):: datai                       !<index for point in
                                                        !< input data
   real(8), intent(INOUT) :: Albedo                        !<Surface albedo
   real(8), dimension(0:16), intent(INOUT)::Tmp            !< Temperatures for each
                                                        !< layer
   type(couplingVariables), intent(INOUT) :: coupling   !< variables used in
                                                        !< coupling(adjusting
                                                        !< radiation to fit
                                                        !< observed surface
                                                        !< temperature)

   type(surfaceVariables), intent(INOUT) :: surf        !< Variables for surface
                                                        !< properties
   type(inputArrays), intent(INOUT) :: modelInput          !< Arrays for model

   integer :: i,couplingLen

   couplingLen=coupling%couplingEndI(1)-coupling%couplingStartI(1)+1
   datai = coupling%saveDatai
   surf%TsurfAve = coupling%TsurfAveSave
   surf%srfWatmms = coupling%SrfWatmmsSave
   surf%srfIce2mms = coupling%SrfIce2mmsSave
   surf%srfIce2mms = coupling%SrfIce2mmsSave
   surf%srfDepmms = coupling%SrfDepmmsSave
   surf%srfSnowmms = coupling%SrfSnowmmsSave
   Albedo = coupling%AlbedoSave
   surf%VeryCold = coupling%VeryColdSave
   do i = 0, 16
      Tmp(i) = coupling%TmpSave(i)
   end do

   do i=1,couplingLen
      modelInput%SW(coupling%couplingStartI(1)+i-1)=coupling%SWSave(i)
      modelInput%SW_dir(coupling%couplingStartI(1)+i-1)=coupling%SWDirSave(i)
      modelInput%LW(coupling%couplingStartI(1)+i-1)=coupling%LWSave(i)
   end do 

end subroutine

!>Check if there is a need to force snow/ice melting
!! to prevent coupling from getting stuck
Subroutine snowIceCheck(LastTsurfObs, surf, CP)
   use RoadSurfVariables
   implicit none
   real(8), intent(IN) :: LastTsurfObs             !< latest surface temperature
                                                !< observation
   type(roadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(surfaceVariables), intent(INOUT) :: surf !< Variables for surface properties

   if (LastTsurfObs > CP%TLimMeltSnow .and. surf%SrfSnowmms > 0.00) Then
      !Force snow melting
      surf%SrfWatmms = surf%SrfWatmms + surf%SrfSnowmms
      surf%SrfSnowmms = 0.00
   end if
   if (LastTsurfObs > CP%TLimMeltIce .and. surf%SrfIcemms > 0.00) Then
      !Force ice melting
      surf%SrfWatmms = surf%SrfWatmms + surf%Srficemms
      surf%SrfIcemms = 0.00
   end if
   if (LastTsurfObs > CP%TLimMeltIce .and. surf%SrfIce2mms > 0.00) Then
      !Force ice melting
      surf%SrfWatmms = surf%SrfWatmms + surf%Srfice2mms
      surf%SrfIce2mms = 0.00
   end if
   if (LastTsurfObs > CP%TLimMeltDep .and. surf%SrfDepmms > 0.00) Then
      !Force ice melting
      surf%SrfWatmms = surf%SrfWatmms + surf%SrfDepmms
      surf%SrfDepmms = 0.00
   end if

end subroutine

!>Determine new coupling values if necessary
Subroutine Coupling_control(TSurfAve, coupling)
   use RoadSurfVariables
   Implicit None
   real(8), intent(INOUT) ::    TsurfAve           !< average temperature of the
                                                !< first two layers
   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)
   real(8) ::    TDifBelow !> Difference between nearest value below TsurfAve and TsurfAve
   real(8) ::    TdifAbove !> Difference between nearest value above TsurfAve and TsurfAve

   !Crevier, L. and Y. Delage, 2001: METRo: A new model for road-condition forecasting
   !in Canada. Journal of Applied Meteorology, 40(11), 2026-2037
   !
   !Karsisto, V. P. Nurmi, M. Kangas, M. Hippi, C. Fortelius, S. Niemelä and 
   !H. Järvinen, 2016: Improving road weather model forecasts by adjusting the radiation
   !input. Meteorological Applications, 23, 503-513

   coupling%start_coupling_again=.false.
   !change to Kelvins
   TsurfAve = TsurfAve + 273.16
   coupling%LastTsurfObs = coupling%LastTsurfObs + 273.16

   !If coupling is not failed
   If (.not. coupling%Coupling_failed) Then

      If (coupling%Coupling_iterations .eq. 0) Then
         coupling%Tsurf_end_coup1 = TsurfAve !Save first value (radcof=1)
      END IF

      !If too much iterations
      If (coupling%Coupling_iterations .eq. 25) Then

         ! If first guess was better than last, return to start and use radcof=1
         ! Otherwise continue from here with radcof=1
         if (abs(coupling%Tsurf_end_coup1 - coupling%LastTsurfObs) <&
          abs(TsurfAve - coupling%LastTsurfObs)) Then
            coupling%start_coupling_again = .true.
         end if

         coupling%SwRadCof = 1.0
         coupling%LWRadCof = 1.0
         coupling%SW_correction = 0.0
         coupling%LW_correction = 0.0
         coupling%RadCoeff = 1.0
         coupling%Coupling_failed = .true.

         !If surface temperature observation is missing, return to start of
         !< coupling and use radcof=1
      Else if (coupling%LastTsurfObs < -100) Then
         coupling%SwRadCof = 1.0
         coupling%LWRadCof = 1.0
         coupling%SW_correction = 0.0
         coupling%LW_correction = 0.0
         coupling%RadCoeff = 1.0
         coupling%Coupling_failed = .true.
         coupling%start_coupling_again = .true.

         !IF abnormal temperature value, return to start of coupling and use radcof=1
      Else if (TsurfAve<170.0 .or. TsurfAve>400.0 .or.&
       coupling%Coupling_failed) Then
         coupling%SwRadCof = 1.0
         coupling%LWRadCof = 1.0
         coupling%SW_correction = 0.0
         coupling%LW_correction = 0.0
         coupling%Coupling_failed = .true.
         coupling%start_coupling_again = .true.
         coupling%RadCoeff = 1.0

         !If forecasted Tsurf greater than observed
      Else if (TsurfAve - coupling%LastTsurfObs .gt. 0.1) Then

         !Save the guess first time
         if (coupling%TsurfNearestAbove < -100) Then
            coupling%TsurfNearestAbove = TsurfAve
            coupling%RadCoefNearestAbove = coupling%RadCoeff

            !compare to previous guess, if it is nearer observed temperature save it
         else if (coupling%TsurfNearestAbove - coupling%LastTsurfObs > &
                  TsurfAve - coupling%LastTsurfObs) Then
            coupling%TsurfNearestAbove = TsurfAve
            coupling%RadCoefNearestAbove = coupling%RadCoeff
         end if

         coupling%start_coupling_again = .true.

         !Calculate new guess, if there is one guess where observed temperature
         !is greater than forecasted and one where observed is smaller
         If (coupling%TsurfNearestAbove > -100 .and. &
             coupling%TsurfNearestBelow > -100) Then
            TDifAbove = coupling%TsurfNearestAbove - coupling%LastTsurfObs
            TDifBelow = coupling%LastTsurfObs - coupling%TsurfNearestBelow
            coupling%RadCoeff = coupling%RadCoefNearestAbove - &
                      TDifAbove/(TDifAbove +TDifBelow)* &
                      (coupling%RadCoefNearestAbove -coupling%RadCoefNearestBelow)

         Else
            !Decrease radcof by half
            coupling%RadCoeff = 0.5*coupling%RadCoeff 
         END IF

         !If not change, reset
         If (abs(coupling%RadCoeff - coupling%RadCoeffPrevious) < 0.00005) Then
            coupling%TsurfNearestAbove = -9999
            coupling%TsurfNearestBelow = -9999
         end if

         if (coupling%RadCoeff < 0.01) Then
            write (*, *) "coupling coefficient too small, coupling failed"
            coupling%RadCoeff =1.0
            coupling%coupling_failed = .true.
            coupling%SwRadCof = 1.0
            coupling%LWRadCof = 1.0
            coupling%SW_correction = 0.0
            coupling%LW_correction = 0.0
         end if

         coupling%RadCoeffPrevious = coupling%RadCoeff
         !If forecasted Tsurf is smaller that observed
      Else if (coupling%LastTsurfObs - TsurfAve .gt. 0.1) Then

         !Save guess the first time
         if (coupling%TsurfNearestBelow < -100) Then
            coupling%TsurfNearestBelow = TsurfAve
            coupling%RadCoefNearestBelow = coupling%RadCoeff

            !compare to previous guess, if it is nearer observed temperature save it
         else if (coupling%TsurfNearestBelow - coupling%LastTsurfObs < TsurfAve&
          - coupling%LastTsurfObs) Then
            coupling%TsurfNearestBelow = TsurfAve
            coupling%RadCoefNearestBelow = coupling%RadCoeff
         end if

         coupling%start_coupling_again = .true.

         !Calculate new guess, if there is one guess where observed temperature
         !is greater than forecasted and one where observed is smaller
         If (coupling%TsurfNearestAbove > -100 .and. coupling%TsurfNearestBelow&
          > -100) Then
            TDifAbove = coupling%TsurfNearestAbove - coupling%LastTsurfObs
            TDifBelow = coupling%LastTsurfObs - coupling%TsurfNearestBelow
            coupling%RadCoeff = coupling%RadCoefNearestAbove - &
                TDifAbove/(TDifAbove +TDifBelow)* &
                (coupling%RadCoefNearestAbove -coupling%RadCoefNearestBelow)

         Else
            !Increace radcof
            coupling%RadCoeff = 2.0*coupling%RadCoeff 
            
         End if

         !If not change, reset
         If (abs(coupling%RadCoeff - coupling%RadCoeffPrevious) < 0.00005) Then
            coupling%TsurfNearestAbove = -9999
            coupling%TsurfNearestBelow = -9999
         end if
         coupling%RadCoeffPrevious = coupling%RadCoeff
      Else
         if (coupling%RadCoeff > 3.0) Then
            write (*, *) "coupling coefficient too big, coupling failed"
            coupling%coupling_failed = .true.
            coupling%RadCoeff =1.0
            coupling%SwRadCof = 1.0
            coupling%LWRadCof = 1.0
            coupling%SW_correction = 0.0
            coupling%LW_correction = 0.0
         end if
         !Coupling was successful
         coupling%SW_correction = coupling%SwRadCof - 1.0
         coupling%LW_correction = coupling%LWRadCof - 1.0
         coupling%Coupling_failed = .false.
         coupling%Coupling_iterations = -1
         coupling%TsurfNearestAbove = -9999.0
         coupling%TsurfNearestBelow = -9999.0
         coupling%RadCoeff = 1.0
         coupling%RadCoefNearestAbove = -9999.0
         coupling%RadCoefNearestBelow = -9999.0
         coupling%RadCoeffPrevious = 1.0
         if (coupling%CoupPhaseN < coupling%NObs) Then
            coupling%CoupPhaseN = coupling%CoupPhaseN + 1
            coupling%LastTsurfObs = coupling%obsTsurf(coupling%CoupPhaseN) + 273.16
         end if
      END IF

   END IF
   !Return to Celcius
   TsurfAve = TsurfAve - 273.16
   coupling%LastTsurfObs = coupling%LastTsurfObs - 273.16
End Subroutine

!>Determine coupling times from car obs times
!>This function assumes that coupling can be done multiple times, separately for each car observation
!>However, multicoupling feature is not used in the current model version
Subroutine initCouplingTimes(coupling, settings)
   use RoadSurfVariables
   Implicit None

   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)

   type(modelSettings), intent(INOUT) :: settings !< Variables for model settings
   real(8) :: DTs                                  !< model time step in seconds

   integer:: i
   integer:: couplingLen

   DTs = settings%DTSecs

   !initialize
   Do i = 1, 48
      coupling%couplingStartI(i) = -99
      coupling%couplingEndI(i) = -99
   end do

   !If only one coupling time
   if (settings%use_coupling .and. coupling%obsI(1) > -1) Then
      coupling%couplingEndI(1) = coupling%obsI(1)
      if (coupling%obsI(1) <= settings%coupling_minutes*60/DTs) then
         coupling%couplingStartI(1) = 1
      else
         !start coupling three hours before the observation
         coupling%couplingStartI(1) = coupling%obsI(1) - &
                                      int(settings%coupling_minutes*60/DTs) 

      end if

   else
      settings%use_coupling = .false.

   end if
   couplingLen=coupling%couplingEndI(1)-coupling%couplingStartI(1)+1
   allocate(coupling%SWSave(couplingLen))
   allocate(coupling%SWDirSave(couplingLen))
   allocate(coupling%LWSave(couplingLen))
   Do i =1,couplingLen
       coupling%SWSave(i)=-9999.9
       coupling%SWDirSave(i)=-9999.9
       coupling%LWSave(i)=-9999.9
   end do
End subroutine

!>Determine radiation correction coefficient to use in the
!! model when coefficient is given as input.
!! (The coefficient approaches gradually one as the forecast advances)
subroutine couplingCofWithInputRadCof(inputRC, coupling, i, settings)
   use RoadSurfVariables
   Implicit none

   integer, intent(IN) :: i                     !<index for point in input data
   type(inputRadiationCoefficient), intent(IN) :: inputRC !< Variables used when
                                                          !< radiation
                                                          !< coefficient is
                                                          !< given as input

   type(modelSettings), intent(IN) :: settings  !< Variables for model settings
   type(couplingVariables), intent(INOUT) :: coupling !< variables used in
                                                      !< coupling(adjusting
                                                      !< radiation to fit
                                                      !< observed surface temperature)


   !If it is right time to use input radiation coefficient
   if (i >= inputRC%inputRadCofI) Then
      !Calculate radiation coefficients for the current time step
      coupling%SWRadCof = 1.0+coupling%SW_correction*exp(-(( settings%DTSecs*i) - &
                           (settings%DTSecs*inputRC%inputRadCofI))/settings%couplingEffectReduction)
      coupling%LWRadCof = 1.0+coupling%LW_correction*exp(-(( settings%DTSecs*i) - &
                           (settings%DTSecs*inputRC%inputRadCofI))/settings%couplingEffectReduction)
   end if

end subroutine couplingCofWithInputRadCof
