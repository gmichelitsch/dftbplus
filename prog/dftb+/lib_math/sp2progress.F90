!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to apply the sp2 method using the progress library.
!>
!> This is a special interface to apply the sp2 solvers from the progress library for the dftb+
!> code.
!>
!> See https://github.com/lanl/qmd-progress for more information
!>
module sp2progress
  use dftbp_bml
  use dftbp_progress
  use constants
  use sparse2bml
  use periodic
  use parallelks
  use densedescr
  use orbitals
  use message
  implicit none
  private

  public :: TSp2Solver


  type :: TSp2Solver
    private
    integer :: iGenZ = 0
    logical :: tInitSp2 = .false.
    logical :: tInitZ = .false.
    type(bml_matrix_t) :: ham
    type(bml_matrix_t) :: over
    type(bml_matrix_t) :: rho
    type(bml_matrix_t) :: orthoH   
    type(bml_matrix_t) :: orthoRho
    type(bml_matrix_t) :: zk1
    type(bml_matrix_t) :: zk2
    type(bml_matrix_t) :: zk3
    type(bml_matrix_t) :: zk4
    type(bml_matrix_t) :: zk5
    type(bml_matrix_t) :: zk6
    type(bml_matrix_t) :: zMat
    type(genZSPinp) :: zsp
    type(sp2data_type) :: sp2
  contains
    procedure :: getDensity => TSp2Solver_getDensity
    procedure :: buildZMatrix => TSp2Solver_buildZMatrix
  end type TSp2Solver
  

contains

  !> This routine implements the sp2 technique with the routines from the progress library.
  subroutine TSp2Solver_getDensity(this, ham, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, parallelKS, denseDesc, orb, nEl, rhoPrim)
    class(TSp2Solver), intent(inout) :: this
    real(dp), intent(in) ::  ham(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iSparseStart(:,:)
    integer, intent(in) :: img2CentCell(:)
    type(TParallelKS), intent(in) :: parallelKS
    type(TDenseDescr), intent(in) :: denseDesc
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: nEl(:)
    real(dp), intent(out) :: rhoPrim(:,:)
    
    real(dp) ::  bandFilling
    integer ::  iSpin
    character(100) :: errorStr

    print *, "CALLIING SP2 get Density"
    
    ! Parsing sp2 input parameters.
    if (.not. this%tInitSp2)then
      call prg_parse_sp2(this%sp2, "progress.in")
      if (this%sp2%mdim < 0) then
        this%sp2%mdim = orb%norb
      end if
      this%tInitSp2 = .true.
    end if

    if (bml_get_type(this%zMat) /= this%sp2%bml_type) then
      call error("SP2 and ZSP bml types differ")
    end if

    if (bml_get_M(this%zMat) /= this%sp2%mdim) then
      call error("SP2 and ZSP M dimension for sparse matrices are different")
      write(errorStr, "(A,I0,A,I0,A)") "SP2 and ZSP M dimension are different (", this%sp2%mdim,&
          & " versus ", bml_get_M(this%zMat), ")"
    endif

    iSpin = parallelKS%localKS(2, 1)
    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, orb%norb, this%sp2%mdim, this%ham)
    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, orb%norb, this%sp2%mdim,&
        & this%orthoH)
    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, orb%norb, this%sp2%mdim,&
        & this%orthoRho)

    if (bml_get_N(this%ham) > orb%norb) then
      call bml_clear(this%ham)
    end if

    call foldToRealBml(ham(:,iSpin), neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell, this%ham)

    call prg_orthogonalize(this%ham, this%zMat, this%orthoH, this%sp2%threshold, this%sp2%bml_type,&
        & this%sp2%verbose)

    ! WARNING: Works only for spin unpolarized cases
    bandFilling = (nEl(iSpin) / 2.00_dp) / real(orb%nOrb, dp)

    ! Perform SP2 from progress
    if (this%sp2%flavor == "Basic")then
      call prg_sp2_basic(this%orthoH, this%orthoRho, this%sp2%threshold, bandFilling,&
          & this%sp2%minsp2iter, this%sp2%maxsp2iter, this%sp2%sp2conv, this%sp2%sp2tol,&
          & this%sp2%verbose)
    else if (this%sp2%flavor == "Alg1") then
      call prg_sp2_alg1(this%orthoH, this%orthoRho, this%sp2%threshold, bandFilling,&
          & this%sp2%minsp2iter, this%sp2%maxsp2iter, this%sp2%sp2conv, this%sp2%sp2tol,&
          & this%sp2%verbose)
    else if(this%sp2%flavor == "Alg2") then
      call prg_sp2_alg2(this%orthoH, this%orthoRho, this%sp2%threshold, bandFilling,&
          & this%sp2%minsp2iter, this%sp2%maxsp2iter, this%sp2%sp2conv, this%sp2%sp2tol,&
          & this%sp2%verbose)
    else
      write(errorStr, "(3A)") "Invalid SP2 flavor '", this%sp2%flavor, "'"
      call error(errorStr)
    endif

    call bml_deallocate(this%orthoh)
    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, orb%nOrb, this%sp2%mdim, this%rho)

    if (bml_get_N(this%ham) > orb%norb) then
      call bml_clear(this%rho)
    end if

    call prg_deorthogonalize(this%orthoRho, this%zMat, this%rho, this%sp2%threshold,&
        & this%sp2%bml_type, this%sp2%verbose)
    call bml_deallocate(this%orthoRho)

    ! Transforming rho from bml to sparse dftb+.
    rhoPrim(:,:) = 0.0_dp
    call unfoldFromRealBml(this%rho, neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell, rhoPrim(:,iSpin))

  end subroutine TSp2Solver_getDensity

  
  !> Computes the inverse overlap congruence transform.
  !>
  !> Computes the inverse overlap needed to orthogonalized the Hamiltonian
  !> matrix before applying the SP2 algorithm. The factors of the inverse
  !> overlap congruence transform can be constructed in with an O(N) scaling.
  !> Details of this implementation can be found in:
  !> http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00154
  !>
  subroutine TSp2Solver_buildZMatrix(this, over, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, parallelKS, denseDesc, orb)

    !> Instance.
    class(TSp2Solver), intent(inout) :: this

    !> Overlap matrix in the dftb+ compressed format.
    real(dp), intent(in) ::  over(:)

    !> Neighbor list for each atom.
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> K-points and spins.
    type(TParallelKS), intent(in) :: parallelKS

    !> Dense matrix descriptor.
    type(TDenseDescr), intent(in) :: denseDesc

    !> Orbital (basis) information.
    type(TOrbitals), intent(in) :: orb

    ! Parse input file.
    if(.not. this%tInitZ)then
      call prg_parse_zsp(this%zsp, "progress.in")
      if(this%zsp%mdim < 0) then
        this%zsp%mdim = orb%nOrb
      end if
      this%tInitZ = .true.
    end if
    
    this%iGenZ = this%iGenZ + 1

    if (bml_get_N(this%over) <= 0) then
      call bml_zero_matrix(this%zsp%bml_type, bml_element_real, dp, orb%nOrb,&
          & this%zsp%mdim, this%zMat)
      call bml_zero_matrix(this%zsp%bml_type, bml_element_real, dp, orb%nOrb, this%zsp%mdim,&
          & this%over)
    end if
    
    call bml_zero_matrix(this%zsp%bml_type, bml_element_real, dp, orb%nOrb, this%zsp%mdim,&
        & this%over)

    ! From dftb+ to bml format
    call foldToRealBml(over, neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell, this%over)

    ! Congruence transformation.
    if (this%zsp%zsp) then
      ! Sparse Z matrix building
      call prg_buildzsparse(this%over, this%zMat, this%iGenZ, this%zsp%mdim, this%zsp%bml_type,&
          & this%zk1, this%zk2, this%zk3, this%zk4, this%zk5, this%zk6, this%zsp%nfirst,&
          & this%zsp%nrefi, this%zsp%nreff, this%zsp%numthresi, this%zsp%numthresf,&
          & this%zsp%integration, this%zsp%verbose)
    else
      ! Build Z matrix using diagonalization (usual method).
      call prg_buildzdiag(this%over, this%zMat, this%zsp%numthresf, this%zsp%mdim,&
          & this%zsp%bml_type)
    end if

  end subroutine TSp2Solver_buildZMatrix
  

end module sp2progress
