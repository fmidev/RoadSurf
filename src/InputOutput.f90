
!>Set input parameters defined in modelRunner.cpp
Subroutine setInputParam(localParam, atm, coupling, settings)

   use RoadSurfVariables

   Implicit None

   type(LocalParameters), intent(IN) :: localParam  !< local parameters given by
                                                    !< modelRunner.cpp
   type(ModelSettings), intent(INOUT) :: settings   !< Variables for model settings
   type(AtmVariables), intent(INOUT) :: atm         !< Variables for atmospheric
                                                    !< properties
   type(CouplingVariables), intent(OUT) :: coupling !< variables used in coupling

   coupling%NObs = -99

   atm%TairR = real(localParam%tair_relax, 4)
   atm%VZR = real(localParam%VZ_relax, 4)
   atm%RhzR = real(localParam%RH_relax, 4)

   if (atm%TairR<-100.0 .or. atm%TairR>100.0 .or. atm%VZR <0.0 .or. & 
       atm%VZR>100.0 .or. atm%RhzR<0.0 .or. atm%RhzR>110)then
       settings%use_relaxation=.false.
   end if
    
   !initialize coupling obs arrays
   call initTsurfObsArrays(coupling)
   coupling%obsI(1) = localParam%couplingIndexI
   coupling%obsTsurf(1) = localParam%couplingTsurf
   coupling%lastTsurfObs = localParam%couplingTsurf
   coupling%NObs = 1
   if (localParam%couplingTsurf < -100 .or. coupling%obsI(1)<1) then
      settings%use_coupling = .false.
   end if


End Subroutine setInputParam

Submodule (RoadSurf) ValueControl
   Implicit None
   contains 
     !> Checks input data for abnormal values
      module Subroutine CheckValues(modelInput, i, settings, surf,localParam)
         use RoadSurfVariables
      
         type(inputArrays), intent(INOUT) :: modelInput  !< Arrays for model input data
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
         type(ModelSettings), intent(INOUT) :: settings !< Variables for model settings
         type(LocalParameters), intent(IN) :: localParam  !< local parameters given by
                                                          !< modelRunner.cpp
      
         if (modelInput%Tair(i) < -90.0 .or. modelInput%Tair(i) > 100.0 &
             .or. modelInput%Tdew(i) < -90 .or. modelInput%Tdew(i) > 100.0 &
             .or. modelInput%RHz(i) < -0.1 .or. modelInput%RHz(i) > 120.0 &
             .or. modelInput%VZ(i) < -1.0 .or. modelInput%VZ(i) > 100.0 &
             .or. modelInput%SW(i) < -0.1 .or. modelInput%SW(i) > 4000.0 &
             .or. modelInput%LW(i) < -0.1 .or. modelInput%LW(i) > 1000.0 &
             .or. modelInput%prec(i) < -0.1 .or. modelInput%prec(i) > 500.0) Then
      
            write (*, *) "BAD input value! ",  modelInput%Tair(i), &
               modelInput%Tdew(i),modelInput%RHz(i), modelInput%VZ(i), modelInput%SW(i),&
               modelInput%LW(i),  modelInput%prec(i)
            settings%simulation_failed = .true.
         end if
         if (localParam%sky_view<1.0 .and. localParam%sky_view>-0.01) Then
             if (modelInput%SW_dir(i)< -0.1 .or. modelInput%SW_dir(i) > 4000.0 &
                 .or. modelInput%LW_net(i) < -1000.0 .or. modelInput%LW_net(i) > 1000.0) Then
                write(*,*) "BAD input value: SW_dir,LW_net",modelInput%SW_dir(i),modelInput%LW_net(i)
                settings%simulation_failed = .true.
             end if
         end if
         if (modelInput%SW_dir(i)>modelInput%SW(i)) then
            modelInput%SW_dir(i)=modelInput%SW(i)
         end if
         if (surf%TSurfAve < -100.0 .or. surf%TsurfAve > 100.0) Then
            write (*, *) "Abnormal surface temperature", surf%TSurfAve,&
                          i,localParam%lat, localParam%lon
            settings%simulation_failed = .true.
         end if
      
      end subroutine CheckValues
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
      
         real(8) ::depth                                  !< depth to take temperature
         real(8) :: t_output                              !< temperature at depth
      
         atm%Tair = modelInput%Tair(i)
         atm%Tdew = modelInput%Tdew(i)
         atm%VZ = modelInput%VZ(i)
         atm%Rhz = modelInput%RHz(i)
         atm%PrecInTStep = modelInput%prec(i)/3600*settings%DTSecs
         Tmp(0) = atm%Tair
      
         !call CalcTDew(REAL(atm%Tair, 8), atm%TDew, REAL(atm%Rhz, 8))
         !If in initialization phase, set surface temperature to observed value
         if (i <= settings%InitLenI) Then
            !if not in coupling phase
            if (modelInput%TsurfOBS(i) > -100.0) Then
      
               if (.not. settings%use_coupling .or. &
                   i < coupling%couplingStartI(coupling%CoupPhaseN)) then
                  surf%TSurfObs = modelInput%TsurfOBS(i)
                  Tmp(1) = surf%TSurfObs
                  Tmp(2) = surf%TSurfObs
                  if (settings%tsurfOutputDepth>=0.0) Then
                     depth=settings%tsurfOutputDepth
                  else
                     !If not, set to value given in input (-9999.9 if missing)
                     depth=modelInput%depth(i)
                  end if
        
                  !interpolate temperature from given depth if if is not missihg
                  if (depth>=0)Then
                     Call getTempAtDepth(ground%Tmp,ground%ZDpth,depth,t_output)
                     surf%TsurfAve=t_output
                  else
                     surf%TSurfAve = (Tmp(1) + Tmp(2))/2.0
                  end if
      
               else
                  surf%TSurfObs = -9999.0
               end if
      
            else
               surf%TSurfObs = -9999.0
            end if
      
         end if
      End subroutine
      !>Save values to output arrays
      module Subroutine SaveOutput(modelOutput, i, surf)
      
         use RoadSurfVariables
         integer, intent(IN) ::i                      !< index of inputdata time steps
         type(SurfaceVariables), intent(IN) :: surf   !< Variables for surface properties
                                                      !< properties
         type(OutputArrays), intent(INOUT) :: modelOutput !< Arrays for model input data
      
         modelOutput%SnowOut(i) = surf%SrfSnowmms
         modelOutput%WaterOut(i) = surf%SrfWatmms
         modelOutput%IceOut(i) = surf%SrfIcemms
         modelOutput%Ice2Out(i) = surf%SrfIce2mms
         modelOutput%DepositOut(i) = surf%SrfDepmms
         modelOutput%TsurfOut(i) = surf%TsurfAve
      End Subroutine
end submodule

!>Sets last input values for interpolated values
Subroutine lastValues(modelInput, atm, settings, ground, surf)
   use RoadSurfVariables
   Implicit None

   type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
   type(inputArrays), intent(IN) :: modelInput  !< Arrays for model input data
   type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                !< properties
   type(GroundVariables), intent(INOUT) :: ground !< Varibales for ground properties
   type(SurfaceVariables), intent(INOUT) :: surf !< Variables for surface properties
   real(8) :: depth
   real(8) :: t_output

   atm%Tair = modelInput%Tair(settings%SimLen)
   atm%Tdew = modelInput%Tdew(settings%SimLen)
   atm%VZ = modelInput%VZ(settings%SimLen)
   atm%Rhz = modelInput%RHz(settings%SimLen)
   atm%PrecInTStep = modelInput%prec(settings%SimLen)/3600*settings%DTSecs
   ground%Tmp(0) = atm%Tair
   depth=modelInput%depth(settings%SimLen)
   !Calculate temperaturea at given depth
   if (depth>=0) Then
      Call getTempAtDepth(ground%Tmp,ground%ZDpth,depth,t_output)
      surf%TsurfAve=t_output
   else
      surf%TSurfAve = (ground%Tmp(1) + ground%Tmp(2))/2.0 !Average temperature of 
                                                           !< the first two layers
    end if

end subroutine

!-----------------------------------------------------------------------------
!> Calculate relative humidity from temperature and dew point
Subroutine CalcRh(T2m, Tdew, Rhz, SimLen)

   IMPLICIT NONE

   real(8), dimension(SimLen), intent(IN)::T2m  !< Air temperature (C)
   real(8), dimension(SimLen), intent(IN)::Tdew !< Dew point temperature (C)
   integer, intent(IN) :: SimLen                !< Lenght of simulation
   real(8), dimension(SimLen), intent(OUT) :: Rhz !< relative humidity (%)

   Real(8), Parameter :: AFact = 0.61078           !< e in kPa
   Real(8), Parameter :: Alphai = 21.875           !< over ice
   Real(8), Parameter :: Betai = 265.5             !< over ice
   Real(8), Parameter :: Alphaw = 17.269           !< over water
   Real(8), Parameter :: Betaw = 237.3             !< over water

   integer :: i
   Real(8)            :: Alpha, Beta, ESatT, ESatTD, MissValR

   MissValR = -9999.9
   Do i = 1, SimLen
      If ((Tdew(i) > MissValR) .and. (T2m(i) > MissValR)) Then
         If (T2m(i) >= 0.) Then
            Alpha = Alphaw
            Beta = Betaw
         Else
            Alpha = Alphai
            Beta = Betai
         End If
         ESatT = AFact*Exp(Alpha*T2m(i)/(T2m(i) + Beta))
         ESatTD = AFact*Exp(Alpha*Tdew(i)/(Tdew(i) + Beta))
         Rhz(i) = Min((ESatTD/ESatT)*100.0, 100.0)
      End If
   End Do

End Subroutine CalcRh

!> Calculate dew point from temperature and humidity (only for one value)
Subroutine CalcTDew(T2m, Tdew, Rhz)
!-------------------------------------------------------------------------

   IMPLICIT NONE

   Real(8), intent(IN) :: T2m   !< Air temperature (C)
   Real(8), intent(IN) :: Rhz   !< relative humidity (%)
   Real(8), intent(OUT):: Tdew     !< Dew point temperature (C)

   Real(8), Parameter :: AFact = 0.61078 ! e in kPa
   Real(8), Parameter :: Alphai = 21.875 ! over ice
   Real(8), Parameter :: Betai = 265.5 ! over ice
   Real(8), Parameter :: Alphaw = 17.269 ! over water
   Real(8), Parameter :: Betaw = 237.3 ! over water

   Real(8)         :: Alpha, Beta, EPr, EprSat, XX

   If (T2m >= 0.) Then
      Alpha = Alphaw
      Beta = Betaw
   Else
      Alpha = Alphai
      Beta = Betai
   End If
   EPrSat = AFact*Exp(Alpha*T2m/(T2m+Beta))
   Epr = 0.01*Rhz*EprSAt
   XX = Log(Epr/AFact)
   TDew = Beta*XX/(Alpha - XX)

End Subroutine CalcTDew
