#include "Constants.h"

Submodule (RoadSurf) prec2Storage
   implicit none
   contains
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
      
      ! ************* PRECIPTATION TYPE
         call CalcPrecType(PrecPhase, settings%DTSecs, atm, CP)
         surf%SrfWatmms = surf%SrfWatmms + atm%RainmmTS ! Rain (snow water not included)
         surf%SrfSnowmms = surf%SrfSnowmms + atm%SnowmmTS ! Snow precipitation (mm snow)
      End subroutine

end submodule
!>Calculate water storage one time step forward
Subroutine WaterStorage(MaxPormms, WatWear, SrfExtmms, &
                        SrfPormms, surf, CP)
   use RoadSurfVariables
   Implicit none

   real(8), intent(IN) ::MaxPormms                 !< maximum water in asphalt pores
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties

   real(8), intent(INOUT) :: WatWear               !< Water storage reduction
                                                !< caused by traffic
   real(8), intent(OUT) :: SrfExtmms               !< Surface water content
   real(8), intent(OUT) :: SrfPormms               !< Water in pores

   ! ************* WATER STORAGE

! * Evaporation / condensation
   !Only for bare and warm surface
   If ((surf%SrfSnowmms <= 0.0) .and. (surf%SrfIcemms <= 0.0) &
       .and. (surf%SrfDepmms <= 0.0) .and. (surf%TSurfAve > CP%TLimDew)) Then
      If (surf%SrfWatmms > MaxPormms) Then
         !surface evaporation/condensation
         surf%SrfWatmms = surf%SrfWatmms - surf%EvapmmTS 
      Else
         !pore evaporation/condensation
         surf%SrfWatmms = surf%SrfWatmms - CP%PorEvaF*surf%EvapmmTS    
  
       End If
   End If

! * Water wear by traffic
   If (surf%wearSurf .and. surf%SrfWatmms > 0.0) Then ! Wear (only above treshold)
      If (surf%SrfWatmms < CP%WWearLim) WatWear = 0.0 ! - no wear below specified limit
      If (surf%SrfWatmms > CP%WWetLim) Then
         surf%SrfWatmms = surf%SrfWatmms - WatWear ! - water amount is reduced

      Else
        !less wear below wet limit
         surf%SrfWatmms = surf%SrfWatmms - CP%DampWearF*WatWear 

      End If
   End If

! * Water storage limits checked also at the end of storage calculations
   !stop from going negative
   If (surf%SrfWatmms < CP%MinWatmms) surf%SrfWatmms = 0.0 
   If (surf%SrfWatmms > CP%MaxWatmms) surf%SrfWatmms = CP%MaxWatmms ! - overflow

   SrfExtmms = Max((surf%SrfWatmms - MaxPormms), 0.) ! Water on surface
   SrfPormms = Min(surf%SrfWatmms, MaxPormms) ! Water in pores
End subroutine

!-------------------------------------------------------------------------------------
!>Calculate snow storage one time step forward
Subroutine SnowStorage(SrfExtmms, Melted, DTSecs, wearF, SrfPormms, &
                       MaxPormms, surf, CP, atm)

   use RoadSurfVariables
   Implicit None

   real(8), intent(IN) ::DTSecs                    !< Time step in seconds
   type(WearingFactors), intent(IN) :: wearF    !< wearing factors
   real(8), intent(IN) ::MaxPormms                 !< Maximum water content in pores
   type(RoadCondParameters), intent(INOUT) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
   type(AtmVariables), Intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties

   real(8), intent(INOUT) :: SrfExtmms             !< surface water content
   real(8), intent(INOUT) :: Melted                !< Melted amount (metres/timestep)
   real(8), intent(OUT) :: SrfPormms               !< Water in pores
   real(8) :: WatSnowRat                           !< Water to snow ratio
   real(8) :: RDummy

   WatSnowRat = 0.0 ! (Surface) water to water+snow ratio
! ************* SNOW STORAGE
! - Heat needed for melting is calculated at the end of the subroutine, and
!   melting heat balance in RoadTemp


   RDummy = SrfExtmms + surf%SrfSnowmms ! Total snow+water in mm water
   If (RDummy > 0.001) Then !
      WatSnowRat = SrfExtmms/RDummy ! Surface water to snow+water ratio
   Else !
      WatSnowRat = 0.0 ! Dry snow with no water
   End If

   RDummy = surf%SrfSnowmms + surf%SrfIcemms ! Total snow+water in mm water
   If (RDummy > 0.001) Then !
      CP%SnowIceRat = surf%SrfSnowmms/RDummy ! Surface water to snow+water ratio
   Else !
      CP%SnowIceRat = 0.0 ! Interpret surface icy
   End If

   If (surf%SrfSnowmms > 0.0) Then
      !wet snow : forming (initially dry)
      If (WatSnowRat > CP%WetSnowFormR) atm%SnowType = SURFACE_SNOW_WET    
   Else ! - only dry => wet allowed
      atm%SnowType = SURFACE_SNOW_DRY ! Dry set to default value when
   End If !  snow storage gets empty

   If (surf%SrfSnowmms > 0.0) Then
      If (surf%SrfDepmms > 0.0) Then ! Deposit under snow ..
         surf%SrfIcemms = surf%SrfIcemms + surf%SrfDepmms ! ... to ice
         surf%SrfDepmms = 0.0 ! Deposit storage to zero
      End If
   End If

   If (surf%SrfSnowmms > 0.0) Then
      !force mwlting when salt on road
      if (CP%forceSnowMelting)Then
         surf%SrfWatmms = surf%srfWatmms+surf%SrfSnowmms
         surf%SrfSnowmms = 0.0
      !normal melting
      else if ((surf%Q2Melt > 0.0) .and. (surf%TSurfAve >= CP%TLimMeltSnow)) Then 
         !Melted amount (meters/timestep)   
         Melted = (surf%Q2Melt*DTSecs)/(CP%WatMHeat*CP%WatDens)
         surf%SrfSnowmms = surf%SrfSnowmms - 1000.*Melted ! Snow melts ...
         surf%SrfWatmms = surf%SrfWatmms + 1000.*Melted ! ... and forms water
      End If
   End If
   If (surf%wearSurf .and. surf%SrfSnowmms > 0.0) Then ! Wear : snow wears to ice

         !snow amount is reduced to ice
         surf%SrfSnowmms = surf%SrfSnowmms - wearF%SnowTran 
         surf%SrfIcemms = surf%SrfIcemms + CP%Snow2IceFac*wearF%SnowTran 
         surf%SrfIce2mms = surf%SrfIce2mms + CP%Snow2IceFac*wearF%SnowTran !
   End If
   !Wet snow
   If ((surf%SrfSnowmms > 0.0) .and. (atm%SnowType == SURFACE_SNOW_WET)) Then 
      If (WatSnowRat > CP%WetSnowMeltR) Then ! - Melting (high water content)
         surf%SrfWatmms = surf%SrfWatmms + surf%SrfSnowmms !   * forms water
         surf%SrfSnowmms = 0.0 !   * snow melts
         atm%SnowType = SURFACE_SNOW_DRY !   * default snow type
      End If
      If (surf%TSurfAve < CP%TLimFreeze) Then ! - Freezing
         !< all at once to ice
         surf%SrfIcemms = surf%SrfIcemms + surf%SrfSnowmms + surf%SrfWatmms 
         ! including water storage                                                  
         surf%SrfIce2mms = surf%SrfIce2mms + surf%SrfSnowmms + surf%SrfWatmms 
                                                            
         atm%SnowType = SURFACE_SNOW_DRY !   * default snow type
         If(surf%SrfSnowmms>0.5) THEN
               CP%WetSnowFrozen=.true.  !Wet snow is frozen
         END IF
         surf%SrfSnowmms = 0.0
         surf%SrfWatmms = 0.0

      End If
   End If

   SrfExtmms = Max((surf%SrfWatmms - MaxPormms), 0.) ! Water on surface
   SrfPormms = Min(surf%SrfWatmms, MaxPormms) ! Water in pores

   If (surf%SrfSnowmms < CP%MinSnowmms) surf%SrfSnowmms = 0.0 ! Stop from going
                                                              !< too small
   If (surf%SrfSnowmms > CP%MaxSnowmms) Then ! Snow "overflow"
      surf%SrfSnowmms = surf%SrfSnowmms - (CP%MaxSnowmms/2.) ! - reduce snow
                                                             !< amount to half
   End If

end subroutine
!------------------------------------------------------------------------------------
!>Calculate ice storage one time step forward
Subroutine IceStorage(Melted, SrfExtmms, SrfPormms, MaxPormms, &
                      DTSecs, surf, CP, wearF)

   use RoadSurfVariables
   Implicit None

   real(8), intent(IN) :: MaxPormms                !< Maximum water content in pores
   real(8), intent(IN) :: DTSEcs                   !< time step in seconds
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(WearingFactors), intent(IN) :: wearF    !< wearing factors
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties

   real(8), intent(INOUT) :: Melted                !< Melted amount (metres/timestep)
   real(8), intent(OUT) :: SrfExtmms               !< surface water content
   real(8), intent(OUT) :: SrfPormms               !< water content in pores

   ! *************  ICE STORAGE
! - Heat needed for melting is calculated at the end of the subroutine, and
!   melting heat balance in RoadTemp

   If (surf%TSurfAve < CP%TLimFreeze .and. surf%SrfWatmms>0.0) Then ! Freezing
      surf%SrfIcemms = surf%SrfIcemms + surf%SrfWatmms ! - all water to ice
      surf%SrfIce2mms = surf%SrfIce2mms + surf%SrfWatmms ! - no freezing for deposit
      surf%SrfWatmms = 0.0

   End If
   If ((surf%SrfSnowmms <= 0.) .and. (surf%SrfIcemms > 0.)) Then ! Melting
      if (CP%forceIceMelting)Then
         surf%SrfWatmms=surf%SrfWatmms+surf%SrfIcemms
         surf%SrfIcemms=0.0
         surf%SrfIce2mms=0.0
      !< only on snowfree ice 
      else if ((surf%Q2Melt > 0.0) .and. (surf%TSurfAve >= CP%TLimMeltIce)) Then 
        !Melted amount (meters/timestep)
         Melted = (surf%Q2Melt*DTSecs)/(CP%WatMHeat*CP%WatDens) 

         surf%SrfIcemms = surf%SrfIcemms - 1000.*Melted ! - both at same rate
         surf%SrfIce2mms = surf%SrfIce2mms - 1000.*Melted
         surf%SrfWatmms = surf%SrfWatmms + 1000.*Melted ! ... adds water storage
      End If
   End If
   If (surf%wearSurf .and. surf%SrfIcemms > 0.) Then ! Wear : also under snow
      !Ice amount is reduced to secondary ice at faster rate
      surf%SrfIcemms = surf%SrfIcemms - wearF%IceWear 
      surf%SrfIce2mms = surf%SrfIce2mms - wearF%IceWear2 

   End If

   SrfExtmms = Max((surf%SrfWatmms - MaxPormms), 0.) ! Water on surface
   SrfPormms = Min(surf%SrfWatmms, MaxPormms) ! Water in pores

   If (surf%SrfIcemms < CP%MinIcemms) surf%SrfIcemms = 0.0 ! Stop from going too
                                                           !< small
   If (surf%SrfIcemms > CP%MaxIcemms) Then ! Ice "overflow"
      surf%SrfIcemms = CP%MaxIcemms ! - reduce to maximum
   End If

   If (surf%SrfIce2mms < CP%MinIcemms) surf%SrfIce2mms = 0.0 ! Stop from going
                                                             !< too small
   If (surf%SrfIce2mms > CP%MaxIcemms) Then ! Ice "overflow"
      surf%SrfIce2mms = CP%MaxIcemms ! - reduce to maximum
   End If
   
End Subroutine

!-------------------------------------------------------------------------------------
!>Calculate deposit storage one time step forward
Subroutine DepositStorage(DepWear, SrfExtmms, SrfPormms, &
                          MaxPormms, surf, CP)

   use RoadSurfVariables
   Implicit None

   real(8), intent(IN) :: DepWear                  !< Deposit reduction caused by
                                                !< traffic
   real(8), intent(IN) :: MaxPormms                !< Maximum water content in pores
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
   real(8), intent(OUT) :: SrfExtmms               !< surface water content
   real(8), intent(OUT) :: SrfPormms               !< water content in pores

   ! *************  DEPOSIT STORAGE
! - Removed conditions for ice-free surface and TSurf <= CP%TLimDew

   If (surf%EvapmmTS < 0.0) Then ! Condensation only
      surf%SrfDepmms = surf%SrfDepmms - surf%EvapmmTS ! - no evaporation
   End If

   If (surf%TSurfAve > CP%TLimMeltDep) Then ! Melting
      surf%SrfWatmms = surf%SrfWatmms + surf%SrfDepmms !  - increase water storage
      surf%SrfDepmms = 0.0 !  - all deposit to water
   End If

   If (surf%wearSurf .and. (surf%SrfSnowmms <= 0.0) & ! Wear only on snow-free
                                                      !< surface
       .and. (surf%SrfDepmms > 0)) Then
      surf%SrfDepmms = surf%SrfDepmms - DepWear
   End If

   SrfExtmms = Max((surf%SrfWatmms - MaxPormms), 0.) ! Water on surface
   SrfPormms = Min(surf%SrfWatmms, MaxPormms) ! Water in pores
   If (surf%SrfDepmms < CP%MinDepmms) surf%SrfDepmms = 0.0 ! Stop from going too
                                                           !< small
   If (surf%SrfDepmms > CP%MaxDepmms) Then ! Deposit "overflow" ?
       !excess into water storage
      surf%SrfWatmms = surf%SrfWatmms + (surf%SrfDepmms - CP%MaxDepmms) 
      surf%SrfDepmms = CP%MaxDepmms ! - reduce to maximum
   End If

End Subroutine

!----------------------------------------------------------------------------------

!>Calculate melting
Subroutine melting(HStor, inCouplingPhase, TsurfObsLast, HS, TmpNw,ZDpth,depth,&
                  surf,CP)
   use RoadSurfVariables
   Implicit None
   real(8), intent(IN) :: HStor                        !< Descriptes stored heat to
                                                    !< the surface from previous
                                                    !< time step

   real(8), intent(IN) :: TsurfObsLast                 !< last surface temperature
                                                    !< observation
   logical, intent(IN) ::inCouplingPhase            !< true if in coupling phase
   real(8), dimension(16), intent(IN):: HS             !< Heat capacity in intensity units
                                                    !< (W/m^2K)
   type(SurfaceVariables), intent(INOUT) :: surf    !< Variables for surface
                                                    !< properties

   real(8), dimension(0:16), intent(INOUT):: TmpNw     !< Next time step layer
                                                    !< temperatures
   real(8), dimension(16), intent(INOUT):: ZDpth      !< Layer depths
   real(8) :: depth                                 !< depth to calculate outpu temp
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   real(8) :: QAvail
   real(8) :: QLeftOver
   real(8) :: t_output

!*************** SNOW/ICE  (only for wearing surface)
! *  Q2Melt = Amount of heat needed for  ice/snow (W/m2)
! *  T4Melt = Limit temperature for  ice/snow
! * NOTE : Only snow OR ice can melt;
! *  Q2Melt in : amount heat (W/m2) needed to melt all ice/snow
! *       out : actual amount of heat (W/m2) used to melt ice/snow
! * During , surface temperature kept at  T4Melt+0.01
! * Even if roadcond determined a Q2 value during the previous calculation step, no
!   calculation is performed if  conditions are not met here, and Q2 is set to zero
!    - temperature may have changed during this time step
!
!
! * Wearing surface

   If ((surf%SrfSnowmms > 0.0) .or. (surf%SrfIcemms > 0.0) &
       .or. (surf%SrfIce2mms > 0.0)) Then
      MLT: Do
         !Melting is forced in salty situation
         if (.not. CP%CanMeltingChangeTemperature) Then
            Exit MLT
         end if
         If ((HStor <= 0.00001) .or. (surf%TSurfAve <= surf%T4Melt) .or. &
            (surf%Q2Melt <=0) .or. &
            (inCouplingPhase .and. TsurfObsLast < surf%T4Melt)) Then

            
            if (surf%TSurfAve < 0.5) Then

               surf%Q2Melt = 0.0 ! No
               Exit MLT ! => exit with no
            ! No change in temerpature if it is high, would drop too much
            else if  (surf%TSurfAve>2.0) Then
               Exit MLT
            end if
         End If 
         QAvail = HS(1)*(TmpNw(1) - surf%T4Melt) ! Heat available for
         If (surf%Q2Melt >= QAvail) Then ! All available heat used (partly melts)
            surf%Q2Melt = QAvail ! - temp remains at  temp
            TmpNw(1) = surf%T4Melt + 0.01 ! - offset to guarantee freezing
            TmpNw(2) = surf%T4Melt + 0.01
            
         Else ! Only part used => no change in Q2Melt
            QLeftOver = QAvail - surf%Q2Melt ! - heat left over from ...
            TmpNw(1) = surf%T4Melt + (QLeftOver/HS(1)) !  ... increases temperature
            TmpNw(2) = surf%T4Melt + 0.01
         End If

         
         if (depth>=0)Then
            Call getTempAtDepth(TmpNw,ZDpth,depth,t_output)
            surf%TsurfAve=t_output
         else
            surf%TSurfAve = 0.5*(TmpNw(1) + TmpNw(2))
         end if
         Exit MLT
      End Do MLT
   Else
      surf%Q2Melt = 0.0 ! Set to zero value if no
   End If
!

End Subroutine

!-------------------------------------------------------------------------------

!> Calculate the heat needed to melt/freeze the whole uppermost snow/ice layer.
!!Used and updated in energy balance calculation.
!!Thickesses in equivalent water mm.
Subroutine NewMeltFreezeHeat(DTSecs, surf, CP)
   use RoadSurfVariables
   Implicit None

   real(8), intent(IN) ::DTSecs                    !< Model time step in seconds
   type(RoadCondParameters), intent(IN) :: CP   !< Parameters to determine
                                                !< storage terms and road condition
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties

   ! * Melting for snow and ice
   surf%Q2Melt = 0.0 ! Default
   If ((surf%SrfSnowmms > 0.0)) Then
      !! Snow melt heat (W/m2) 
      surf%Q2Melt = CP%WatMHeat*CP%WatDens*(surf%SrfSnowmms/1000.)/DTSecs 
      surf%T4Melt = CP%TLimMeltSnow                                        
   End If
   If ((surf%SrfSnowmms <= 0.0) .and. (surf%SrfIcemms > 0.0)) Then
      !Ice melt heat (W/m2)
      surf%Q2Melt = CP%WatMHeat*CP%WatDens*(surf%SrfIcemms/1000.)/DTSecs 

      surf%T4Melt = CP%TLimMeltIce
   End If
   If (surf%Q2Melt < 0.0) surf%Q2Melt = 0.0 ! Just to be sure ...
End subroutine

!-------------------------------------------------------------------------------
