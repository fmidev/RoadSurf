!> variables used in coupling (adjusting radiation to fit observed surface temperature)
type CouplingVariables
      real(8), dimension(0:16)::TmpSave            !< Saved layer temperatures at
                                                !< the beginning of the coupling
                                                !< period

      integer :: Coupling_iterations            !< number of coupling interations
      real(8) :: TsurfNearestAbove                 !< Suraface temperature at the
                                                !< end of couplgin that is
                                                !< nearest to the observed temperature
                                                !< but is still above it

      real(8) :: TsurfNearestBelow                 !< Suraface temperature at the
                                                !< end of couplgin that is
                                                !< nearest to the observed temperature
                                                !< but is still below it

      real(8) :: RadCoeff                          !< used when determining
                                                !< radiation coefficient
      logical :: Down                           !< forecasted temperature has
                                                !< been greater than observed
      logical :: start_coupling_again           !< if true, start coupling
                                                !< period again
      logical :: Coupling_failed                !< if ture, coupling has been failed,
                                                !< continue without coupling
      logical :: inCouplingPhase                !< if true, model is currently
                                                !< in coupling phase
      logical :: VeryColdSave                   !< Saved veryCold value at the
                                                !< start of the coupling period
      real(8) :: RadCoefNearestAbove               !< coefficient by which TsurfNearestAbove
                                                !< was achieved
      real(8) :: RadCoefNearestBelow               !< coefficient by which TsurfNearestBelow
                                                !< was achieved

      real(8) :: RadCoeffPrevious                 !< Radcoeff2 at previous iteration
      real(8) :: SWRadCof                          !< radiation coefficient to use
                                                !< for short wave radiation
      real(8) :: LWRadCof                          !< radiation coefficient to use
                                                !< for long wave radiation
      real(8) :: SW_correction                     !< final radiation coefficient
                                                !< for short wave radiation
      real(8) :: LW_correction                     !< final radiation coefficient
                                                !< for long wave radiation
      real(8) :: Tsurf_end_coup1                   !< surace temperature at the end
                                                !< of first coupling period
      integer, dimension(48):: couplingStartI   !< i index to start coupling
      integer, dimension(48):: couplingEndI     !< i index to end coupling
      integer, dimension(48):: couplingStartJ   !< j index to start coupling
      integer, dimension(48):: couplingEndJ     !< j index to end coupling
      real(8) :: TSurfAveSave                      !< Saved surface tmperature at
                                                !< the start of the coupling period
      real(8) :: SrfWatmmsSave                     !< Saved water stroage at the
                                                !< start of the coupling period
      real(8) :: SrfIcemmsSave                     !< Saved ice storage at the
                                                !< start of the coupling period
      real(8) :: SrfIce2mmsSave                    !< Saved secondary ice storage
                                                !< at the start of the coupling
                                                !< period

      real(8) :: SrfDepmmsSave                     !< Saved deposit storage at the
                                                !< start of the coupling period
      real(8) :: SrfSnowmmsSave                    !< Saved snow storage at the
                                                !< start of the coupling period
      real(8):: AlbedoSave                        !< Saved albedo at the start of
                                                !< the coupling period
      real(8) :: lastTsurfObs                      !< last surface temperature
                                                !< observation
      integer :: saveDatai                      !< i at the start of the
                                                !< coupling period
      integer :: saveDataj                      !< j index at the start of the
                                                !< coupling period
      real(8), dimension(:), allocatable:: SWSave  !< saved sw data
      real(8), dimension(:), allocatable:: SWDirSave  !< saved dir sw data
      real(8), dimension(:), allocatable:: LWSave  !< saved lw data
      integer :: NObs                           !< number of car observations
      integer :: CoupPhaseN                     !< Coupling Phase Number
      integer, dimension(48, 6):: carObsTime    !< observation time for car obs (year,
                                                !< month,day,hour,minute second)
      integer, dimension(48)::obsI              !< i index corresponding car
                                                !< observation time
      real(8), dimension(48)  :: obsTsurf          !< observed car surface temperature

end type CouplingVariables
