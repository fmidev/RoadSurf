!>Variables used when radiation coefficient is given as input
!>NOT CURRENTLY IN USE IN THE MODEL
type InputRadiationCoefficient
      real:: inputRadCofSW                      !< Input short wave radiation coeffient
      real :: inputRadCofLW                     !< Input long waver radiation coefficient
      Integer :: inputRadCofI                   !< i index to start using input
                                                !< radiation coefficient
      integer, dimension(6)::timeToUseInputRadCof !< time to use input radiation
                                                  !< coefficient
end type InputRadiationCoefficient
