!> To apply the sp2 method using the progress library.
!! \brief This is a special interface to apply the sp2 solvers from
!! the progress library for the dftb+ code.
!!
!! See https://github.com/lanl/qmd-progress for more information
!!
#:include 'common.fypp'

module sp2progress

#:if WITH_PROGRESS
  use bml
  use prg_sp2_mod
  use prg_sp2parser_mod
  use prg_genz_mod
  use prg_nonortho_mod

  use constants
  use sparse2bml
  use periodic
  use parallelks
  use densedescr
  use orbitals

  private

  public :: sp2prg, buildZprg

  integer, public             ::  igenz
  logical, public             ::  initsp2 = .false.
  logical, public             ::  initz = .false.
  type(bml_matrix_t)          ::  orthoh_bml, orthoz_bml
  type(bml_matrix_t), public  ::  ham_bml, over_bml, rho_bml, zk1_bml
  type(bml_matrix_t), public  ::  zk2_bml, zk3_bml, zk4_bml, zk5_bml
  type(bml_matrix_t), public  ::  zk6_bml, zmat_bml
  type(genZSPinp), public     ::  zsp
  type(sp2data_type), public  ::  sp2

contains

  !> This routine implements the sp2 technique with the
  !! routines from the progress library.
  !!
  subroutine sp2prg(sp2, ham, zmat_bml, rhoPrim, neighborList, nNeighbor, iSparseStart, img2CentCell,&
       & solver, parallelKS, denseDesc, orb, nEl, initsp2)
    implicit none
    integer                          ::  iSpin, mdim
    integer, intent(in)              ::  iSparseStart(:,:), img2CentCell(:), nNeighbor(:), solver
    logical                          ::  initsp2
    real(dp)                         ::  bndfil
    real(dp), intent(in)             ::  ham(:,:), nEl(:)
    type(TDenseDescr), intent(in)    ::  denseDesc
    type(TNeighborList), intent(in)  ::  neighborList
    type(TOrbitals), intent(in)      ::  orb
    type(TParallelKS), intent(in)    ::  parallelKS
    type(bml_matrix_t)               ::  ham_bml, orthoH_bml, orthorho_bml, rho_bml
    type(bml_matrix_t)               ::  zmat_bml
    type(sp2data_type)               ::  sp2
    real(dp), intent(inout)          ::  rhoPrim(:,:)

    !Parsing sp2 input parameters.
    if(.not.initsp2)then
       call prg_parse_sp2(sp2,"progress.in")
       if(sp2%mdim < 0)sp2%mdim = orb%norb
       initsp2 = .true.
    endif

    if(bml_get_type(zmat_bml) .ne. sp2%bml_type)then
       stop "ERROR: SP2 and ZSP bml types differ"
    endif

    if(bml_get_M(zmat_bml) .ne. sp2%mdim)then
       write(*,*)bml_get_M(zmat_bml),sp2%mdim
       stop "ERROR: SP2 and ZSP M dimension for sparse matrices are different"
    endif

    iSpin = parallelKS%localKS(2, 1)
    call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,orb%norb,sp2%mdim,ham_bml)
    call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,orb%norb,sp2%mdim,orthoH_bml)
    call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,orb%norb,sp2%mdim,orthorho_bml)

    if(bml_get_N(ham_bml) > orb%norb) call bml_clear(ham_bml)

    !From dftb+ to bml format
    call foldToRealBml(ham(:,iSpin), neighborList%iNeighbor, nNeighbor, orb, &
         denseDesc%iAtomStart, iSparseStart, img2CentCell, ham_bml) !, sp2%threshold)

    call prg_orthogonalize(ham_bml,zmat_bml,orthoH_bml,sp2%threshold,sp2%bml_type,sp2%verbose)

    bndfil=(nEl(iSpin)/2.00_dp)/real(orb%norb,dp) !WARNING: ONLY FOR SPIN UNPOLARIZED

    !Perform SP2 from progress
    if(sp2%flavor == "Basic")then
       call prg_sp2_basic(orthoH_bml, orthoRho_bml,sp2%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter, &
            & sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    elseif(sp2%flavor == "Alg1")then
       call prg_sp2_alg1(orthoH_bml,orthoRho_bml,sp2%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter, &
            & sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    elseif(sp2%flavor == "Alg2")then
       call prg_sp2_alg2(orthoH_bml,orthoRho_bml,sp2%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter, &
            & sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    else
       stop "No valid SP2 flavor"
    endif

    call bml_deallocate(orthoh_bml)

    call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,orb%norb,sp2%mdim,rho_bml)

    if(bml_get_N(ham_bml) > orb%norb) call bml_clear(rho_bml)

    call prg_deorthogonalize(orthoRho_bml,zmat_bml,rho_bml,sp2%threshold,sp2%bml_type,sp2%verbose)

    call bml_deallocate(orthoRho_bml)

    !Transforming rho from bml to sparse dftb+.
    rhoPrim = 0.0_dp
    call unfoldFromRealBml(rho_bml, neighborList%iNeighbor, nNeighbor, orb, &
         & denseDesc%iAtomStart, iSparseStart, img2CentCell, rhoPrim(:,iSpin))

    ! call bml_print_matrix("rho_bml", rho_bml, 1, 10, 1, 10)

  end subroutine sp2prg

  !> Routine to compute the inverse overlap congruence transform.
  !> \brief Computes the inverse overlap needed to orthogonalized the Hamiltonian
  !! matrix before applying the SP2 algorithm. The factors of the inverse
  !! overlap congruence transform can be constructed in with an O(N) scaling.
  !! Details of this implementation can be found in:
  !! http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00154
  !!
  !! \param over Overlap matrix in the dftb+ compressed format.
  !! \param over_bml Overlap matrix in the BML format (could be dense or ellpack).
  !! \param zk1_bml-zk6_bml Last six Z matrices.
  !! \param zsp Structure with settings that are read from an imput file.
  !! \param neighborList Neighbor list for each atom.
  !! \param nNeighbor Number of neighbors for each atom (incl. itself).
  !! \param iSparseStart Indexing array for the sparse Hamiltonian.
  !! \param img2CentCell Atomic mapping indexes.
  !! \param parallelKS K-points and spins.
  !! \param denseDesc Dense matrix descriptor
  !! \param orb Orbital (basis) information.
  !! \param initz Logical variable that tracks the initialization of this routine.
  !!
  subroutine buildZprg(zsp, over, over_bml, zmat_bml, zk1_bml, zk2_bml, zk3_bml, &
       & zk4_bml, zk5_bml, zk6_bml, neighborList, nNeighbor, &
       & iSparseStart, img2CentCell, parallelKS, denseDesc, orb, initz)
    implicit none
    integer                            ::  iSpin, mdim
    integer, intent(in)                ::  iSparseStart(:,:), img2CentCell(:), nNeighbor(:)
    logical                            ::  initz
    real(dp), intent(in)               ::  over(:)
    type(TDenseDescr), intent(in)      ::  denseDesc
    type(TNeighborList), intent(in)    ::  neighborList
    type(TOrbitals), intent(in)        ::  orb
    type(TParallelKS), intent(in)      ::  parallelKS
    type(bml_matrix_t)                 ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  over_bml
    type(bml_matrix_t), intent(inout)  ::  zk1_bml, zk2_bml, zk3_bml
    type(bml_matrix_t), intent(inout)  ::  zk4_bml, zk5_bml, zk6_bml, zmat_bml
    type(genZSPinp)                    ::  zsp

    !Parse input file.
    if(.not.initz)then
       call prg_parse_zsp(zsp,"progress.in")
       if(zsp%mdim < 0) zsp%mdim = orb%norb
       initz = .true.
    endif

    igenz = igenz + 1

    if(bml_get_N(over_bml).le.0)then
       call bml_zero_matrix(zsp%bml_type,bml_element_real,dp,orb%norb,zsp%mdim,zmat_bml)
       call bml_zero_matrix(zsp%bml_type,bml_element_real,dp,orb%norb,zsp%mdim,over_bml)
    endif

    iSpin = parallelKS%localKS(2, 1)
    call bml_zero_matrix(zsp%bml_type,bml_element_real,dp,orb%norb,zsp%mdim,over_bml)

    !From dftb+ to bml format
    call foldToRealBml(over, neighborList%iNeighbor, nNeighbor, orb, &
         denseDesc%iAtomStart, iSparseStart, img2CentCell, over_bml) !, sp2%threshold)

    if(zsp%zsp)then !Congruence transformation.
       !Build Z using the technique implemented in http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00154
       call prg_buildzsparse(over_bml,zmat_bml,igenz,zsp%mdim,&
            zsp%bml_type, zk1_bml,zk2_bml,zk3_bml&
            ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
            zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)
    else
       !Build Z matrix using diagonalization (usual method).
       call prg_buildzdiag(over_bml,zmat_bml,zsp%numthresf,zsp%mdim,zsp%bml_type)
    endif

  end subroutine buildZprg

#:endif

end module sp2progress
