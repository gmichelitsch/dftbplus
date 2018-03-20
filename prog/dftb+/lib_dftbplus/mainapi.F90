!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_mainapi
  use dftbp_environment, only : TEnvironment
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_main, only : processGeometry, getMullikenPopulation
  use dftbp_initprogram, only : initProgramVariables, destructProgramVariables, coord0, latVec,&
      & tCoordsChanged, tLatticeChanged, energy, derivs, TRefExtPot, refExtPot, tExtField, orb,&
      & nAtom, nSpin, q0, qOutput, qBlockOut, qiBlockOut, qInput, qBlockIn, qiBlockIn, sccCalc,&
      & tExtChrg, tForces, chrgForces, ham, over, rhoPrim, iRhoPrim, neighborList, nNeighbor,&
      & denseDesc, iSparseStart, img2CentCell
  use dftbp_spin, only : qm2ud, ud2qm
  use dftbp_sparse2dense, only : packHS, unpackHS, blockSymmetrizeHS
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges
  public :: setRealDensity, getRealHamiltonian, getRealOverlap, getNonInteractingDensMat
  public :: denseMatrixShape

  logical :: tChargesChanged = .true.

contains


  subroutine setGeometry(coords, latVecs)
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: latVecs(:,:)
    
    coord0(:,:) = coords
    tCoordsChanged = .true.
    if (present(latVecs)) then
      latVec(:,:) = latVecs
      tLatticeChanged = .true.
    else
      tLatticeChanged = .false.
    end if

  end subroutine setGeometry


  subroutine getEnergy(env, merminEnergy)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: merminEnergy

    call recalcGeometry(env)
    merminEnergy = energy%EMermin
    
  end subroutine getEnergy


  subroutine getGradients(env, gradients)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: gradients(:,:)

    call recalcGeometry(env)
    gradients(:,:) = derivs
    
  end subroutine getGradients


  subroutine getGrossCharges(atomCharges)
    real(dp), intent(out) :: atomCharges(:)

    atomCharges(:) = sum(q0(:, :, 1) - qOutput(:, :, 1), dim=1)
    
  end subroutine getGrossCharges


  !> Sets up an external electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be negative.
  !>
  subroutine setExternalPotential(atomPot, shellPot, potGrad)

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(refExtPot%atomPot)) then
        allocate(refExtPot%atomPot(nAtom, nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [nAtom]))
      refExtPot%atomPot(:,1) = -atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(refExtPot%shellPot)) then
        allocate(refExtPot%shellPot(orb%mShell, nAtom, nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, nAtom]))
      refExtPot%shellPot(:,:,1) = -shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(refExtPot%potGrad)) then
        allocate(refExtPot%potGrad(3, nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, nAtom]))
      refExtPot%potGrad(:,:) = -potGrad
    end if
    tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up external point charges
  subroutine setExternalCharges(chargeCoords, chargeQs, blurWidths)

    !> Coordiante of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    tExtChrg = .true.
    if (tForces) then
      if (.not. allocated(chrgForces)) then
        allocate(chrgForces(3, size(chargeQs)))
      end if
    end if
    call sccCalc%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)
    
  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(chargeGradients)

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(tForces .and. allocated(chrgForces))

    chargeGradients(:,:) = chrgForces
    
  end subroutine getExtChargeGradients


  subroutine setRealDensity(env, densityMatrix)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(in) :: densityMatrix(:,:,:)

    integer :: iSpin

    call recalcGeometry(env)
    rhoPrim(:, :) = 0.0_dp
    do iSpin = 1, size(densityMatrix, dim=3)
      call packHS(rhoPrim(:,iSpin), densityMatrix(:,:,iSpin), neighborlist%iNeighbor, nNeighbor,&
          & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
    end do
    call ud2qm(rhoPRim)
    tChargesChanged = .true.
    call recalcCharges(env)
    
  end subroutine setRealDensity


  subroutine getRealHamiltonian(env, hamiltonian)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: hamiltonian(:,:,:)

    integer :: iSpin

    call recalcGeometry(env)
    do iSpin = 1, size(hamiltonian, dim=3)
      call unpackHS(hamiltonian(:,:,iSpin), ham(:,iSpin), neighborList%iNeighbor, nNeighbor,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call blockSymmetrizeHS(hamiltonian(:,:,iSpin), denseDesc%iAtomStart)
    end do
    call qm2ud(hamiltonian)

  end subroutine getRealHamiltonian


  subroutine getRealOverlap(env, overlap)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: overlap(:,:)

    call recalcGeometry(env)
    call unpackHS(overlap, over, neighborList%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call blockSymmetrizeHS(overlap, denseDesc%iAtomStart)

  end subroutine getRealOverlap


  subroutine getNonInteractingDensMat(densityMatrix)
    real(dp), intent(out) :: densityMatrix(:,:,:)

    integer :: iAt, iOrb, iSpin
    integer :: ind

    densityMatrix(:,:,:) = 0.0_dp
    ind = 1
    do iSpin = 1, size(densityMatrix, dim=3)
      do iAt = 1, nAtom
        do iOrb = 1, orb%nOrbAtom(iAt)
          densityMatrix(ind, ind, iSpin) = qInput(iOrb, iAt, iSpin)
          ind = ind + 1
        end do
      end do
    end do
    call qm2ud(densityMatrix)
    
  end subroutine getNonInteractingDensMat


  function denseMatrixShape() result(matDim)
    integer :: matDim(3)

    matDim(:) = [orb%nOrb, orb%nOrb, nSpin]

  end function denseMatrixShape


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine recalcGeometry(env)
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc
    
    if (tLatticeChanged .or. tCoordsChanged) then
      call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc)
      tLatticeChanged = .false.
      tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry


  subroutine recalcCharges(env)
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc
    
    if (.not. tChargesChanged) then
      return
    end if

    call getMullikenPopulation(rhoPrim, over, orb, neighborList, nNeighbor, img2CentCell,&
        & iSparseStart, qInput, iRhoPrim=iRhoPrim, qBlock=qBlockIn, qiBlock=qiBlockIn)
    qOutput(:,:,:) = qInput
    if (allocated(qBlockIn)) then
      qBlockOut(:,:,:,:) = qBlockIn
    end if
    if (allocated(qiBlockIn)) then
      qiBlockOut(:,:,:,:) = qiBlockIn
    end if
    tChargesChanged = .false.
    call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc)
    
  end subroutine recalcCharges

  
end module dftbp_mainapi
