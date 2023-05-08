!> Variables for surface properties
type SurfaceVariables
      real(8) :: TsurfAve            !< average temperature of the
                                     !< first two layers
      real(8) :: SrfWatmms           !< water storage (mm)
      real(8) :: SrfSnowmms          !< snow storage (water equivalent mm)
      real(8) :: SrfIcemms           !< ice storage (water equivalent mm)
      real(8) :: SrfIce2mms          !< Secondary ice storage (water
                                     !< equivalent mm)
      real(8) :: SrfDepmms           !< deposit storage (water
                                     !< equivalent mm)
      real(8) :: Q2Melt              !< Amount of heat needed for 
                                     !< ice/snow to melt (W/m2)/
                                     !< actual amount of heat (W/m2)
                                     !< used to melt ice/snow

      real(8) :: T4Melt              !< Limit temperature for 
                                     !< ice/snow to melt
      real(8) :: TrfFric             !< Surface heating caused by traffic
      real(8) :: EvapmmTS            !< Evaporation
      Logical :: VeryCold            !< true if ground temperature is cold
      Logical :: WearSurf            !< true if storage terms can be
                                     !< reduced by traffic wearing

      real(8):: TsurfOBS             !< surface temperature observation
end type SurfaceVariables
