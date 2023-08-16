!MIT License
!Copyright (c) 2023 FMI Open Development
!> Variables for ground properties
type GroundVariables
      real(8) :: Albedo                        !< surface albedo
      real(8) :: HStor                         !< Descriptes stored heat to the
                                               !< surface from previous time step
      real(8), allocatable:: condDZ(:)         !< variable for temperature profile
                                               !< calculation
      real(8), allocatable:: capDZ(:)          !< variable for temperature profile
                                               !< calculation
      real(8), allocatable:: Wcont(:)          !< water content in ground levels
      real(8), allocatable:: VSH(:)            !< Volumetric heat capacity of
                                               !< ground layers (J/(m^3K))
      real(8), allocatable:: HS(:)             !< Heat capacity in intensity units
                                               !< (W/m^2K)
      real(8), allocatable:: CC(:)             !< ground heat conductivity
      real(8), allocatable:: Tmp(:)            !< Temperatures for each layer
      real(8), allocatable:: TmpNw(:)          !< Next time step layer temperatures
      Real(8), allocatable:: DyC(:)            !< Layers thicknesses for heat
                                               !< capacity calculation
      Real(8), allocatable:: DyK(:)            !< Layer thicknesses for heat conductivity
                                               !< calculation
      Real(8), allocatable:: ZDpth(:)          !< Depths of ground layers
      Real(8), allocatable :: GCond(:)         !< Heat conductivity divided by
                                               !< layer height
      real(8) :: GroundFlux                    !< ground heat flux (between 2nd and
                                               !< 3rd layers)
end type GroundVariables
