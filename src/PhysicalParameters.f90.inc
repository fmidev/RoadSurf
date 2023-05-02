!> Physical parameters used in the model
type physicalParameters
      real(8) :: VK_Const  !< Von Karman's konstant
      real(8) :: SB_const  !< Stefan-boltzman constant (W/m2K4)
      real(8) :: ZRefW     !< Wind reference height (m)
      real(8) :: ZRefT     !< Teperature reference height (m)
      real(8) :: ZeroDisp  !< Zero displacement height (m)
      real(8) :: ZMom      !< Roughness factor for momentum (m)
      real(8) :: ZHeat     !< Roughness factor for heat (m)
      real(8) :: logMom    !< Log((ZRefW+ZMom)/ZMom)
      real(8) :: logHeat   !< Log((ZRefW+ZHeat)/ZHeat)
      real(8) :: logCond   !< Log((phy%ZRefW - phy%ZeroDisp +&
                           !< phy%ZHeat)/phy%ZHeat)
      real(8) :: logUstar  !< Log((phy%ZRefW - phy%ZeroDisp + phy%ZMom)/&
	                   !< phy%ZMom)
      real(8) :: Grav      !< Gravitational acceleration (m/s2)
      real(8) :: Emiss     !< Emissivity constant of the surface
      real(8) ::Afc1       !< Help variable used when calculating heat conductivity
      real(8) :: Bfc1      !< Help variable used when calculating heat conductivity
      real(8) :: Cfc1      !< Help variable used when calculating heat conductivity
      real(8) :: Dfc1      !< Help variable used when calculating heat conductivity
      real(8) :: Efc1      !< Help variable used when calculating heat conductivity
      real(8) :: Afc2      !< Help variable used when calculating heat conductivity
      real(8) :: Bfc2      !< Help variable used when calculating heat conductivity
      real(8) :: Cfc2      !< Help variable used when calculating heat conductivity
      real(8) :: Dfc2      !< Help variable used when calculating heat conductivity
      real(8) :: Efc2      !< Help variable used when calculating heat conductivity
      real(8) ::Poro1      !< porosity of asphalt
      real(8) :: Poro2     !< porosity of ground layers
      real(8) :: vsh1      !< volumetric heat capacity of dry asphalt
      real(8) :: vsh2      !< volumetric heat capacity of dry ground
      real(8) :: LVap      !< latent heat of water vaporization
      real(8) :: LFus      !< latent heat of sublimation
      real(8) :: TClimG    !< Climatological temperature at the bottom ground layer
      real(8) :: MaxPormms !< maximum water in asphalt pores
      real(8) :: DampDpth  !< Damping depth
      real(8) :: Omega     !< Frequency of year variation (2pi/365)
      real(8) :: AZ        !< Amplitude at z(m+1), used to calculate bottom layer
                           !< temperature
      real(8) :: Silt1     !< Clay fraction for surface layers
      real(8) :: Silt2     !< Clay fraction for deep ground layers
      real(8) :: RhoB1     !< Bulk density for surface layers
      real(8) :: RhoB2     !< Bulk density for deep ground layers
end type physicalParameters
