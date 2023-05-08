Submodule (RoadSurf) RelaxModule
   Implicit none
   contains

   !>Use relaxation to air temperature, wind speed and relative humidity
   !> after initialization phase. This is done to avoid jump when moving
   !> from observed atmospheric values to forecasted ones.
   module Subroutine RelaxationOperations(i, atm, settings,Tmp)
      use RoadSurfVariables
   
      integer, intent(IN) :: i                     !< index of inputdata time steps
      type(ModelSettings), intent(IN) :: settings  !< Variables for model settings
      type(AtmVariables), intent(INOUT) :: atm     !< Variables for atmospheric
                                                   !< properties
      real(8), dimension(0:16), intent(INOUT)::Tmp    !< Temperatures for each layer
      real(8) :: DTs
      integer :: initLI
      DTs = settings%DTSecs
      initLI = settings%initLenI
      
      !If current time step is the end of the initialization period
      if (i == initLI) Then
         atm%TairInitEnd = atm%Tair
         atm%VZInitEnd = atm%VZ
         atm%RhzInitEnd = atm%Rhz
      end if
   
      !If not in the initialization period, make adjustments to
      !air temperature, wind speed and relative humidity
      if (i > initLI) Then
         atm%Tair = atm%Tair - (atm%TairR - atm%TairInitEnd) &
                    *exp(-((DTs*i) - (DTs*initLI))/(4.*3600.))
         Tmp(0)=atm%Tair
         atm%VZ = atm%VZ - (atm%VZR - atm%VZInitEnd) &
                  *exp(-((DTs*i) - (DTs*initLI))/(4.*3600.))
         atm%Rhz = atm%Rhz - (atm%RhzR - atm%RhzInitEnd) &
                   *exp(-((DTs*i) - (DTs*initLI))/(4.*3600.))
         if (atm%Rhz > 100.) Then
            atm%Rhz = 100.0
         end if
      end if
   
      call CalcTDew(atm%Tair, atm%TDew, atm%Rhz)
   End Subroutine
end submodule RelaxModule
