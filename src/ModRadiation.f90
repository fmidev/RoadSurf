
submodule (RoadSurf) RadMod
   implicit none
   contains
      !> Uses sky view factor and local horizon angles to modify incoming
      !> radiation fluxes
      module SUBROUTINE ModRadiationBySurroundings(modelInput,inputParam,localParam,i)
      use RoadSurfVariables
         type(inputArrays),intent(INOUT) ::modelInput !< model input arrays
         type(inputParameters),intent(IN) :: inputParam !< model input parameters
         type(LocalParameters), intent(IN) :: localParam  !< Local parameters
                                                          !< given by modelrunner.cpp
         integer, intent(IN) ::i            !< model input index
         
         !LOCAL VARIABLES
         REAL(8) :: sun_azim    ! local azimuth clockwise from north
         REAL(8) :: LW_surroundings  ! longwave radiation from surroundings
         REAL(8) :: sun_elevation ! sun elevation angle
         REAL(8) :: dif_SW      ! diffuse short wave radiation										   
         REAL(8) :: SW_ref      ! reflected short wave radiation by surroundings
         real(8) :: horizon_in_sun_dir !Horizon angle in the direction of the sun
         REAL(8)  :: shadow_fac  !  shadow factor
         INTEGER :: azim_idx    ! index in local horizon array that corresponds
                                ! sun azimuth angle
      
         !Calculation is based on paper: "Parametrization of orographic effects
         !on surface radiation in HIRLAM, Senkova et al. 2007"
         
         !Diffuse radiation is difference between global and direct sw radiation
         dif_SW=modelInput%SW(i)-modelInput%SW_dir(i)
         
         !Caclulate upward lw radiation from net lw radiation and downwelling lw radiation
         !This is the radiation emitted by the road surroundings
         LW_surroundings=modelInput%LW_net(i)-modelInput%LW(i)
         
         !Calculate sun position
         Call SunPosition(modelInput,localParam,i,sun_elevation,sun_azim)
      
         !Take local horizon angle that corresponds sun azimuth
         horizon_in_sun_dir=0.
         azim_idx=NINT(sun_azim)
         if (azim_idx==360) then
            azim_idx=0
         end if
         horizon_in_sun_dir=modelInput%local_horizons(azim_idx+1)
      
         !Location is in shadow if sun elevation is lower than local horizon angle   
         if (horizon_in_sun_dir>(sun_elevation))THEN   
            shadow_fac=0.0
         else
            shadow_fac=1.0
         end if
         
         !If sun is above horizon
         IF (sun_elevation>0.0) THEN  
          !Reduces direct solar radiation if the sun is behind obstacle
          modelInput%SW_dir(i) = modelInput%SW_dir(i)  * shadow_fac  
          !Reflected radiation from surroundings
          SW_ref = inputParam%Albedo_surroundings* modelInput%SW_dir(i) + &
              inputParam%Albedo_surroundings * dif_SW 
          !Total incoming diffuse short wave radiation
          dif_SW =localParam%sky_view  *  dif_SW + (1.0 - localParam%sky_view) *  SW_ref
      
          !Total incoming short wave radiaiont
          modelInput%SW(i) =dif_SW +  modelInput%SW_dir(i) 
      
          END IF  
      
          !Total incoming long wave radiation
          modelInput%LW(i) = localParam%sky_view*modelInput%LW(i)+&
                            (1.0-localParam%sky_view)*(-LW_surroundings)
      
      END SUBROUTINE
 end submodule
