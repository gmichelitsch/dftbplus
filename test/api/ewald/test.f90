program test
  use, intrinsic :: iso_fortran_env, only : stdOut => output_unit
  use dftbp_api_ewald
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: nAtom
  real(dp), allocatable :: invRMat(:,:), coords(:,:)
  character(20) :: formStr

  nAtom = 2
  allocate(invRMat(2, nAtom))
  allocate(coords(3, nAtom))
  coords(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  coords(:, 2) = [1.0_dp, 1.0_dp, 1.0_dp]
  call invR(invRMat, nAtom, coords)
  write(formStr, "(A,I0,A)") "(", nAtom, "E20.12)"
  write(stdOut, "(A)") "Received 1/R matrix:"
  write(stdOut, formStr) transpose(invRMat)
  
end program test
