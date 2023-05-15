!MIT License
!Copyright (c) 2023 FMI Open Development
#include "Constants.h"

submodule (RoadSurf) Cond
   implicit none
   contains
!> Determines road surface condition for program Simulation      !
      module Subroutine RoadCond(MaxPormms, surf, atm, settings, &
                          CP,wearF)
         use RoadSurfVariables
      
      !------------------------------------------------------------------------------!
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
         real    :: SrfExtmms, SrfPormms !< water outside and insied porous media
         Real    :: Melted               !< Melted amount (metres/timestep)
      
      
      
         Melted = 0.0
         CP%SnowIceRat = 0.0 ! Snow to snow+ice ratio
         atm%SnowType = SURFACE_SNOW_DRY
         ! Between the limits, keep the old setting hysteresis : *H > *L)
         If ((surf%VeryCold) .and. (surf%TSurfAve > CP%TLimColdH)) Then
            surf%VeryCold = .false.
         end if
         If ((.not. surf%VeryCold) .and. (surf%TSurfAve < CP%TLimColdL)) Then
            surf%VeryCold= .true.
         end if
      
      ! ************* WATER STORAGE
         call WaterStorage(MaxPormms, wearF%WatWear, SrfExtmms, &
                           SrfPormms, surf, CP)
      
      ! ************* SNOW STORAGE
         call SnowStorage(SrfExtmms, Melted, settings%DTSecs, &
                          wearF, SrfPormms, MaxPormms, surf, CP, &
                           atm)
      
      ! *************  ICE STORAGE
         call IceStorage(Melted, SrfExtmms, SrfPormms, MaxPormms, &
                         settings%DTSecs, surf, CP, wearF)
      
      ! *************  DEPOSIT STORAGE
         call DepositStorage(wearF%DepWear, SrfExtmms, SrfPormms, MaxPormms, &
                             surf, CP)
      
      ! *************  WATER STORAGE LIMITS RECHECK
      ! * Checked also in the beginning of storage calculations
      
         If (surf%SrfWatmms < CP%MinWatmms) surf%SrfWatmms = 0.0 ! Stop from going negative
         If (surf%SrfWatmms > CP%MaxWatmms) surf%SrfWatmms = CP%MaxWatmms ! Overflow
      
      End Subroutine
      
      !--------------------------------------------------------------------------------------------------
      !> Determine wear factors (how much storage terms are reduced by traffic)
      module Subroutine WearFactors(Snow2IceFac, Tph, surf, wearF)
         use RoadSurfVariables
         real(8), intent(IN)    :: Tph                   !< time steps per hour
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(WearingFactors), intent(OUT) :: wearF   !< wearing factors
         Real(8), intent(INOUT)    :: Snow2IceFac        !< Snow to ice transition factor
      
         ! ************* WEAR FACTORS
      
         wearF%SnowTran = (0.2+0.25)*surf%SrfSnowmms
         wearF%SnowTran = Max(wearF%SnowTran, 0.01)
      
      ! In case of a small snow layer (< 0.2mm) wearing more effective (x3)
         If (surf%SrfSnowmms < 0.2) Then
            wearF%SnowTran = wearF%SnowTran*3
         End If
      
         Snow2IceFac = 0.25/(0.2+0.25)
         wearF%SnowTran = wearF%SnowTran*Tph ! - to mm/time step
      
         wearF%IceWear = 1.1*2.0*0.145*surf%SrfIcemms ! - car traffic
         wearF%IceWear = Max(wearF%IceWear, 0.01)
         wearF%IceWear = wearF%IceWear*Tph ! - to mm/time step
         wearF%IceWear2 = 1.1*2.0*(4.0*0.290)*surf%SrfIce2mms ! - car traffic
         wearF%IceWear2 = Max(wearF%IceWear2, 0.01)
         wearF%IceWear2 = wearF%IceWear2*Tph ! - to mm/time step
      
         wearF%DepWear = 0.5*2.0*(4.0*0.290)*surf%SrfDepmms ! - car traffic
         wearF%DepWear = Max(wearF%DepWear, 0.01)
         wearF%DepWear = wearF%DepWear*Tph ! - to mm/time step
      
         wearF%WatWear = 0.145*surf%SrfWatmms ! - car traffic
         wearF%WatWear = Max(wearF%WatWear, 0.06)
         wearF%WatWear = 10*wearF%WatWear*Tph ! - to mm/time step
      End subroutine
      !> Calculates albedo
      module Subroutine CalcAlbedo(albedo, surf, cp)
         use RoadSurfVariables
      
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                      !< storage terms and road condition
      
         real(8), intent(INOUT) :: Albedo                !< surface albedo
         real(8) :: IceSum, IceMax ! Ice thickness dependent albedo variables
      
         ! *************  SURFACE PROPERTIES
         ! * Change in emissivity assumed negligible
         ! * NOTE : surface condition dependent
      
         ! For ice, albedo varies linerily from CP%AlbDry to CP%AlbSnow
         ! between ice+depostit thickness 0-IceMax;
         If (surf%wearSurf) Then ! Only with wearing surface
            IceSum = 0.5*(surf%SrfIcemms + surf%SrfIce2mms) + surf%SrfDepmms
            If (IceSum < 0.0) IceSum = 0.0
            IceMax = 1.5
            Albedo = CP%AlbDry
      
            If (surf%SrfSnowmms>0.01 .and. surf%SrfSnowmms>surf%SrfIcemms) Then
               Albedo = CP%AlbSnow
            
            else if (surf%SrfIcemms>0.01 .or. surf%SrfDepmms>0.01) Then
               If (IceSum < Icemax) Then
                  Albedo = CP%AlbDry + (IceSum/IceMax)*(CP%AlbSnow - CP%AlbDry)
               Else
                  Albedo = CP%AlbSnow
               End If
            End If
         End If
      
      End Subroutine
end submodule
!--------------------------------------------------------------------------------------------------
!> Determine precipitation type
Subroutine CalcPrecType(PrecPhase, DTSecs, atm, CP)
   use RoadSurfVariables
   Implicit None

   integer, intent(IN) :: PrecPhase             !< precPhase (Hail = 6; 
                                                !< FreezingRain =5; 
                                                !< FreezingDrizzle = 4; Snow = 3;
                                                !< Sleet = 2; Rain = 1; Drizzle = 0;)

   real, intent(IN):: DTSecs                    !< Input time step in seconds
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties
   Logical :: UseInterpr                        !< use in-built phase
                                                !< interpretation control
   Real(8)    :: PRain, PExp                       !< Probability coefficients
                                                !< (prec interpretation)

!         ************* PRECIPTATION TYPE
! * Snow amount set by in-built interpretation if left to missing value or if forced
!   - Ref: Koistinen, Saltikoff: Experience on customer products of accumulated snow,
!  sleet and rain(Int.Seminar on Advanced Weather Radar Systems, Switzerland, March 1998)
! * atm%PrecType set to (-1 = not defined; should appear onply if PrecInTStep <MinPrecmm)
!    1 = water
!    2 = sleet
!    3 = snow

   atm%RainmmTS = 0.0
   atm%SnowmmTS = 0.0
   atm%PrecType = -1 ! Reset atm%PrecType
   UseInterpr = .true. ! Reset in-built phase interpretation control
   If (PrecPhase > CP%MissValI) Then
!        * Precipitation PHASE input (MetEd)
!          - FmiNumberOfPrecipitationForms = 7;
!          - kTHail = 6;
!          - kTFreezingRain = 5;
!          - kTFreezingDrizzle = 4;
!          - kTSnow = 3;
!          - kTSleet = 2;
!          - kTRain = 1;
!          - kTDrizzle = 0;

      UseInterpr = .false. ! By default no interpretation
      If (atm%PrecInTStep <= CP%MinPrecmm) Then
         atm%PrecInTStep = 0.0 ! Check minimum level
         atm%PrecType = -1 ! No precipitation -> No atm%PrecType
         atm%RainmmTS=0.0
         atm%SnowmmTS=0.0
      Else
         Select Case (PrecPhase)
            ! WATER (above minimum)
         Case (PRECIPITATION_NONE, PRECIPITATION_RAIN, &
               PRECIPITATION_FREEZING_DRIZZLE, PRECIPITATION_FREEZING_RAIN)
            atm%RainmmTS = atm%PrecInTStep ! - interpretation should handle
            atm%SnowmmTS = 0.0 !   freezing drizzle and rain
            atm%PrecType = 1 ! - precipitation type
            atm%SnowType = SURFACE_SNOW_WET ! - rain on snow makes it wet
         Case (PRECIPITATION_SLEET) ! SLEET
            atm%SnowmmTS = atm%PrecInTStep/2. ! - half water, half snow
            atm%RainmmTS = atm%SnowmmTS ! - atm%SnowType wet
            atm%PrecType = 2 ! - precipitation type
            atm%SnowType = SURFACE_SNOW_WET !
         Case (PRECIPITATION_SNOW, PRECIPITATION_HAIL) ! SNOW

            atm%SnowmmTS = atm%PrecInTStep ! - all prec as snow
            atm%PrecType = 3 ! - precipitation type
            atm%RainmmTS = 0.0 !
         Case Default ! INTERPRETATION MISSING
            UseInterpr = .true. ! - force own interpretation
         End Select
      End If
   End If

! * Intepretation : Used if UseInterpr set to true
!    - default = true ; set to false if snow amount obtained from input (above)
!    - NOTE : rainmm may not have missing value (such cases set to 0 in input)

   If (UseInterpr) Then
      If (atm%PrecInTStep <= CP%MinPrecmm) Then
         atm%PrecInTStep = 0.0 ! Check minimum level
         atm%PrecType = -1 ! No precipitation -> No atm%PrecType
         atm%RainmmTS=0.0
         atm%SnowmmTS=0.0
      Else 
         atm%SnowmmTS = 0.0 ! Amount of snow (mm_water/time_step)
         PExp = 22.0-2.7*atm%Tair - 0.20*atm%Rhz
         PRain = 1.0/(1.0+exp(PExp))
         If (PRain < CP%PLimSnow) Then ! SNOW
            atm%SnowmmTS = atm%PrecInTStep
            atm%PrecType = 3
         Else If (PRain > CP%PLimRain) Then ! WATER
            atm%RainmmTS = atm%PrecInTStep
            atm%SnowType = SURFACE_SNOW_WET ! Rain on snow makes it wet
            atm%PrecType = 1
         Else ! SLEET
            atm%SnowmmTS = atm%PrecInTStep/2. ! - half water, half snow
            atm%RainmmTS = atm%SnowmmTS
            atm%SnowType = SURFACE_SNOW_WET
            atm%PrecType = 2
         End If
      End if
   End If
   atm%RainIntensity = (atm%RainmmTS/DTSecs)*3600.0 ! - intensity, rain mm/h
   atm%SnowIntensity = (atm%SnowmmTS/DTSecs)*3600.0 ! - intensity, snow mm/h

End subroutine
!--------------------------------------------------------------------------------------------------
