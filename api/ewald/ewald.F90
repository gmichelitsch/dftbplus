!> Exports the Ewald summation routines of DFTB+.
module dftbp_api_ewald
  use accuracy, only : dp
  use constants, only : pi
  use simplealgebra, only : determinant33, invert33
  use periodic, only : getCellTranslations, getLatticePoints
  use coulomb
  implicit none
  private

  public :: EwaldCalculator
  public :: init, setLatticeVectors, getPeriodicPotential, getCentralPotential


  type :: EwaldCalculator
    private
    integer :: nAtom
    real(dp) :: cellVol
    real(dp) :: alpha
    real(dp) :: tolerance
    real(dp) :: rCutoff, gCutoff
    real(dp), allocatable :: rCellVecs(:,:), gCellVecs(:,:)
  end type EwaldCalculator
  
  interface init
    module procedure EwaldCalculator_init
  end interface init

  interface setLatticeVectors
    module procedure EwaldCalculator_setLatticeVectors
  end interface setLatticeVectors

  interface getPeriodicPotential
    module procedure EwaldCalculator_getPeriodicPotential
  end interface getPeriodicPotential

  interface getCentralPotential
    module procedure EwaldCalculator_getCentralPotential
  end interface getCentralPotential
  

contains


  !> Initialises an Ewald-sum calculator.
  subroutine EwaldCalculator_init(this, nAtom, alpha, tolerance)

    !> Instance
    type(EwaldCalculator), intent(out) :: this

    !> Number of atoms in the system
    integer, intent(in) :: nAtom

    !> Ewald summation parameter (usually betwee 0.1 and 1.0)
    real(dp), intent(in) :: alpha

    !> Maximal terms to consider in the Ewald summation
    real(dp), intent(in) :: tolerance

    this%nAtom = nAtom
    this%alpha = alpha
    this%tolerance = tolerance
    
  end subroutine EwaldCalculator_init


  !> Sets the lattice vectors for the current lattice.
  subroutine EwaldCalculator_setLatticeVectors(this, latVecs)

    !> Instance.
    type(EwaldCalculator), intent(inout) :: this

    !> Lattice vectors (x,y,z|1,2,3)
    real(dp), intent(in) :: latVecs(:,:)

    real(dp), allocatable :: dummy(:,:)
    real(dp) :: recVecs(3, 3)

    recVecs = latVecs(:,:)
    call invert33(recVecs, latVecs)
    recVecs(:,:) = 2.0_dp * pi * transpose(recVecs)
    this%cellVol = abs(determinant33(latVecs))
    this%rCutoff = getMaxREwald(this%alpha, this%tolerance)
    this%gCutoff = getMaxGEwald(this%alpha, this%cellVol, this%tolerance)
    call getCellTranslations(dummy, this%rCellVecs, latVecs, recVecs / (2.0_dp * pi), this%rCutoff)
    call getLatticePoints(this%gCellVecs, recVecs, latVecs / (2.0_dp * pi), this%gCutoff,&
        & onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    this%gCellVecs(:,:) = matmul(recVecs, this%gCellVecs)

  end subroutine EwaldCalculator_setLatticeVectors


  !> Returns the electrostatic potential contributed by the atoms in the central cell.
  subroutine EwaldCalculator_getCentralPotential(this, coords, charges, potential)

    !> Instance
    type(EwaldCalculator), intent(inout) :: this

    !> Coordinates of the atoms. Shape: (3, nAtom).
    real(dp), intent(in) :: coords(:,:)

    !> Charges of the atoms. Shape: (nAtom)
    real(dp), intent(in) :: charges(:)

    !> Resulting potential. Shape: (nAtom)
    real(dp), intent(out) :: potential(:)

    call sumInvR(potential, this%nAtom, this%nAtom, coords, coords, charges, ignoreSamePos=.true.)

  end subroutine EwaldCalculator_getCentralPotential


  !> Returns the electrostatic potential considering all periodic images.
  subroutine EwaldCalculator_getPeriodicPotential(this, coords, charges, potential)

    !> Instance
    type(EwaldCalculator), intent(inout) :: this

    !> Coordinates of the atoms. Shape: (3, nAtom).
    real(dp), intent(in) :: coords(:,:)

    !> Charges of the atoms. Shape: (nAtom)
    real(dp), intent(in) :: charges(:)

    !> Resulting potential. Shape: (nAtom)
    real(dp), intent(out) :: potential(:)

    call sumInvR(potential, this%nAtom, this%nAtom, coords, coords, charges, this%rCellVecs,&
        & this%gCellVecs, this%alpha, this%cellVol)
    
  end subroutine EwaldCalculator_getPeriodicPotential
    

end module dftbp_api_ewald
