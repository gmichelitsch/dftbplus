!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_hambuider
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus
  use dftbp_constants, only : AA__Bohr
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 3

  integer, parameter :: nExtChrg = 2

  ! H2O coordinates
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! H2O atom types
  integer, parameter :: species(nAtom) = [1, 2, 2]

  character(100), parameter :: slakoFiles(2, 2) = reshape([character(100) :: &
      & "external/slakos/origin/mio-1-1/O-O.skf",&
      & "external/slakos/origin/mio-1-1/H-O.skf",&
      & "external/slakos/origin/mio-1-1/O-H.skf",&
      & "external/slakos/origin/mio-1-1/H-H.skf"], [2, 2])

  character(1), parameter :: maxAngNames(4) = ["s", "p", "d", "f"]

  integer, parameter :: nOrb = 6


  real(dp), parameter :: densityMtx0(nOrb, nOrb, 1) = reshape([&
      &  0.187613546023362E+01_dp, -0.422098028383389E+00_dp, -0.137815635002950E-16_dp,&
      &  0.833683185489424E-32_dp,  0.348806542202691E-01_dp,  0.348806542202691E-01_dp,&
      &  0.000000000000000E+00_dp,  0.109234265677572E+01_dp,  0.503788387448307E-15_dp,&
      &  0.898790576397687E-31_dp,  0.488617118109771E+00_dp,  0.488617118109771E+00_dp,&
      &  0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.107455133453139E+01_dp,&
      &  0.503324339152872E-15_dp,  0.589571777683362E+00_dp, -0.589571777683362E+00_dp,&
      &  0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.000000000000000E+00_dp,&
      &  0.200000000000000E+01_dp, -0.340996956212552E-15_dp,  0.340996956212552E-15_dp,&
      &  0.348806542202691E-01_dp,  0.488617118109771E+00_dp,  0.589571777683362E+00_dp,&
      & -0.340996956212552E-15_dp,  0.571252723492489E+00_dp, -0.757054438291368E-01_dp,&
      &  0.348806542202691E-01_dp,  0.488617118109771E+00_dp, -0.589571777683362E+00_dp,&
      &  0.340996956212552E-15_dp, -0.757054438291368E-01_dp,  0.571252723492489E+00_dp],&
      & [nOrb, nOrb, 1])


  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input

  real(dp), allocatable :: densityMtx(:,:,:), hamiltonian(:,:,:), overlap(:,:)
  real(dp) :: atomCharges(nAtom)
  integer :: denseShape(3)
  integer :: matDim, nSpin
  real(dp) :: coords(3, nAtom)
  integer :: devNull
  type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pParserOpts

  ! Note: setting the global standard output to /dev/null will also suppress run-time error messages
  open(newunit=devNull, file="/dev/null", action="write")
  call TDftbPlus_init(dftbp, outputUnit=devNull, calcOnlyHS=.true.)
  !call TDftbPlus_init(dftbp, calcOnlyHS=.true.)

  ! You should provide the dftb_in.hsd and skfiles as found in the
  ! test/prog/dftb+/non-scc/Si_2/ folder
  call dftbp%getEmptyInput(input)
  call input%getRootNode(pRoot)
  call setChild(pRoot, "Geometry", pGeo)
  call setChildValue(pGeo, "Periodic", .false.)
  call setChildValue(pGeo, "TypeNames", ["O", "H"])
  coords(:,:) = 0.0_dp
  call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coords)
  call setChild(pRoot, "Hamiltonian", pHam)
  call setChild(pHam, "Dftb", pDftb)
  call setChildValue(pDftb, "Scc", .true.)
  call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
  call setChildValue(pMaxAng, "O", maxAngNames(getMaxAngFromSlakoFile(slakoFiles(1, 1)) + 1))
  call setChildValue(pMaxAng, "H", maxAngNames(getMaxAngFromSlakoFile(slakoFiles(2, 2)) + 1))
  call setChild(pDftb, "SlaterKosterFiles", pSlakos)
  call setChild(pSlakos, "Type2FileNames", pType2Files)
  call setChildValue(pType2Files, "Prefix", "external/slakos/origin/mio-1-1/")
  call setChildValue(pType2Files, "Separator", "-")
  call setChildValue(pType2Files, "Suffix", ".skf")
  call setChild(pRoot, "ParserOptions", pParserOpts)
  call setChildValue(pParserOpts, "ParserVersion", 5)
  
  print "(A)", 'Input tree in HSD format:'
  call dumpHsd(input%hsdTree, output_unit)
  
  call dftbp%setupCalculator(input)
  denseShape(:) = dftbp%denseMatrixShape()
  matDim = denseShape(1)
  nSpin = denseShape(3)
  allocate(densityMtx(matDim, matDim, nSpin))
  allocate(hamiltonian(matDim, matDim, nSpin))
  allocate(overlap(matDim, matDim))

  call dftbp%getNonInteractingDensMat(densityMtx)
  print "(A)", "Non-interacting density matrix"
  print "(6F13.8)", densityMtx

  coords(:,:) = initialCoords
  call dftbp%setGeometry(coords)
  call dftbp%getRealOverlap(overlap)
  print "(A)", "Overlap obtained for current geometry"
  print "(6F13.8)", overlap

  densityMtx(:,:,:) = densityMtx0
  print "(A)", "Density matrix set by caller"
  print "(6F13.8)", densityMtx
  call dftbp%setRealDensity(densityMtx)

  call dftbp%getGrossCharges(atomCharges)
  print *, 'Gross atomic charges obtained from density matrix:'
  print "(3F13.8)", atomCharges
  
  call dftbp%getRealHamiltonian(hamiltonian)
  print "(A)", "Hamiltonian obtained for current density matrix"
  print "(6F13.8)", hamiltonian

  call TDftbPlus_destruct(dftbp)


end program test_hambuider
