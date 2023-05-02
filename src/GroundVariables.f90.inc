!> Variables for ground properties
type groundVariables
      real(8) :: Albedo                        !< surface albedo
      real(8) :: HStor                         !< Descriptes stored heat to the
                                            !< surface from previous time step
      real(8), dimension(16):: condDZ          !< variable for temperature profile
                                            !< calculation
      real(8), dimension(16):: capDZ           !< variable for temperature profile
                                            !< calculation
      real(8), dimension(16):: Wcont           !< water content in ground levels
      real(8), dimension(16):: VSH             !< Volumetric heat capacity of
                                            !< ground layers (J/(m^3K))
      real(8), dimension(16):: HS              !< Heat capacity in intensity units
                                            !< (W/m^2K)
      real(8), dimension(16):: CC              !< ground heat conductivity
      real(8), dimension(0:16):: Tmp           !< Temperatures for each layer
      real(8), dimension(0:16):: TmpNw         !< Next time step layer temperatures
      Real(8), dimension(16):: DyC             !< Layers thicknesses for heat
                                            !< capacity calculation
      Real(8), dimension(16):: DyK             !< Layer thicknesses for heat conductivity
                                            !< calculation
      real(8), dimension(16):: ZDpth           !< Depths of ground layers
      Real(8), dimension(0:16) :: GCond        !< Heat conductivity divided by
                                            !< layer height
      real(8) :: GroundFlux                    !< ground heat flux (between 2nd and
                                            !< 3rd layers)
end type groundVariables
