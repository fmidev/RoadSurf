!> wearing ratios due to traffic
type wearingFactors
      Real(8)    :: SnowTran      !< Snow transform factor : mm/timestep,
      Real(8)    :: SnowTran2     !< Snow transform factor nigth,
      Real(8)    :: SnowTranDef   !< default value for transfor
      Real(8)    :: DepWear       !< Deposit wear : mm/timestep
      Real(8)    :: IceWear       !< Ice wear : mm/timestep 
      Real(8)    :: IceWear2      !< Ice wear2 : mm/timestep
      Real(8)    :: IceWearNight  !< Ice wear night : mm/timestep
      Real(8)    :: IceWear2Night !< Ice wear2 night : mm/timestep
      Real(8)    :: IceWearSW     !< Ice wear by sun : mm/timestep
      Real(8)    :: IceWear2SW    !< Ice wear2 by sun: mm/timestep
      Real(8)    :: WatWear       !< Water wear : mm/timestep

end type wearingFactors
