! MODULE: RoadSurfVariables

! DESCRIPTION:
!> This module contains all structures used by road weather model

Module RoadSurfVariables
   USE, INTRINSIC :: ISO_C_BINDING

   Implicit None
  
   INCLUDE 'InputPointers.f90.inc'
   INCLUDE 'OutputPointers.f90.inc'
   INCLUDE 'InputArrays.f90.inc'
   INCLUDE 'InputSettings.f90.inc'
   INCLUDE 'InputParameters.f90.inc'
   INCLUDE 'LocalParameters.f90.inc'
   INCLUDE 'OutputArrays.f90.inc'
   INCLUDE 'PhysicalParameters.f90.inc'
   INCLUDE 'GroundVariables.f90.inc'
   INCLUDE 'SurfaceVariables.f90.inc'
   INCLUDE 'AtmVariables.f90.inc'
   INCLUDE 'CouplingVariables.f90.inc'
   INCLUDE 'ModelSettings.f90.inc'
   INCLUDE 'InputRadiationCoefficient.f90.inc'
   INCLUDE 'RoadCondParameters.f90.inc'
   INCLUDE 'WearingFactors.f90.inc'


End Module RoadSurfVariables
