program test
  use, intrinsic :: iso_fortran_env, only : stdOut => output_unit
  use dftbp_api_ewald
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  type(EwaldCalculator) :: ewaldCalc
  integer :: nAtom
  real(dp) :: latVecs(3, 3)
  real(dp) :: latConst
  real(dp), allocatable :: coords(:,:), grad1(:,:), grad2(:,:)
  real(dp), allocatable :: charges(:), pot1(:), pot2(:)
  real(dp) :: alpha, tolerance
  character(20) :: formStr

  nAtom = 2
  latConst = 5.427092_dp / 0.529177249_dp
  allocate(coords(3, nAtom))
  coords(:,:) = latConst * reshape([0.0_dp, 0.0_dp, 0.0_dp,&
      & 0.25_dp, 0.25_dp, 0.25_dp], [3, nAtom])
  latVecs(:,:) = latConst * reshape([0.0_dp, 0.5_dp, 0.5_dp,&
      & 0.5_dp, 0.0_dp, 0.5_dp,&
      & 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate(charges(nAtom))
  charges(:) = [1.47548659_dp, -1.47548659_dp]
  alpha = 0.491_dp
  tolerance = 1e-9_dp

  allocate(pot1(nAtom))
  allocate(pot2(nAtom))
  allocate(grad1(3, nAtom))
  allocate(grad2(3, nAtom))

  call init(ewaldCalc, nAtom, alpha, tolerance)
  call setLatticeVectors(ewaldCalc, latVecs)
  call getPeriodicPotential(ewaldCalc, coords, charges, pot1)
  call getCentralPotential(ewaldCalc, coords, charges, pot2)
  call getPeriodicGradient(ewaldCalc, coords, charges, grad1)
  call getCentralGradient(ewaldCalc, coords, charges, grad2)

  write(stdOut, "(A)") "Full periodic electrostatic potential:"
  write(stdOut, "(E20.12)") pot1
  write(stdOut, "(A)") "Electrostatic potential contribution from central cell:"
  write(stdOut, "(E20.12)") pot2
  write(stdOut, "(A)") "Difference:"
  write(stdOut, "(E20.12)") pot1 - pot2

  write(stdOut, "(A)") "Full periodic electrostatic gradient:"
  write(stdOut, "(3E20.12)") grad1
  write(stdOut, "(A)") "Electrostatic gradient contribution from central cell:"
  write(stdOut, "(3E20.12)") grad2
  write(stdOut, "(A)") "Difference:"
  write(stdOut, "(3E20.12)") grad1 - grad2

end program test
