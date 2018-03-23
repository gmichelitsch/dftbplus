!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module kernelmixer
  use assert
  use accuracy
  implicit none
  private

  public :: TKernelMixer, TKernelMixerInp
  public :: init, reset, mix


  type :: TKernelMixerInp

    !> Mixing parameter for the mixer
    real(dp) :: mixParam

    !> Displacement for the finite difference
    real(dp) :: delta

  end type TKernelMixerInp


  type :: TKernelMixer
    private

    !> Input charge at the point where the derivative is calculated
    real(dp), allocatable :: qInp0(:)

    !> Output charge at the point where the derivative is calculated
    real(dp), allocatable :: qOut0(:)

    !> Difference betwee output and input charges at the point where deriv. is calculated
    real(dp), allocatable :: qDiff0(:)

    !> v-vector: (qOut - qInp) / ||qOut - qInp||
    real(dp), allocatable :: vv(:)

    !> Output charge when at the positiv displacement during finite difference calculation
    real(dp), allocatable :: qOutPlus(:)

    !> Norm of qDiff0
    real(dp) :: qDiff0Norm

    !> Total charges
    real(dp) :: qTot

    !> Finite difference used to calculate the derivative
    real(dp) :: delta

    !> Mixing parameter
    real(dp) :: cc

    !> Current mixing phase (0 - central, 1 - positive finite diff., -1 - negative finite diff.)
    integer :: phase

  end type TKernelMixer


  interface init
    module procedure TKernelMixer_init
  end interface init


  interface reset
    module procedure TKernelMixer_reset
  end interface reset


  interface mix
    module procedure TKernelMixer_mix
  end interface mix



contains


  !> Creates a Kernel mixer instance.
  subroutine TKernelMixer_init(this, input)

    !> an initialized Kernel mixer on exit
    type(TKernelMixer), intent(out) :: this

    !> Mixing parameter
    type(TKernelMixerInp), intent(in) :: input

    @:ASSERT(input%mixParam > 0.0_dp)
    @:ASSERT(input%delta > 0.0_dp)

    this%cc = input%mixParam
    this%delta = input%delta
    this%phase = 0
    this%qTot = -1.0_dp

  end subroutine TKernelMixer_init


  !> Makes the mixer ready for a new SCC cycle
  subroutine TKernelMixer_reset(this, nElem)

    !> Kernel mixer instance
    type(TKernelMixer), intent(inout) :: this

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    if (allocated(this%qInp0) .and. size(this%qInp0) /= nElem) then
      deallocate(this%qInp0)
      deallocate(this%qOut0)
      deallocate(this%qDiff0)
      deallocate(this%vv)
      deallocate(this%qOutPlus)
    end if
    if (.not. allocated(this%qInp0)) then
      allocate(this%qInp0(nElem))
      allocate(this%qOut0(nElem))
      allocate(this%qDiff0(nElem))
      allocate(this%vv(nElem))
      allocate(this%qOutPlus(nElem))
    end if

  end subroutine TKernelMixer_reset


  !> Mixes charges according to the modified Kernel method
  subroutine TKernelMixer_mix(this, qInpResult, qDiff)

    !> The Kernel mixer
    type(TKernelMixer), intent(inout) :: this

    !> Input charges on entry, mixed charges on exit.
    real(dp), intent(inout) :: qInpResult(:)

    !> Charge difference between output and input charges
    real(dp), intent(in) :: qDiff(:)

    real(dp), allocatable :: uu(:), qOutMinus(:), qOutDeriv(:)

    !print *, 'Entering mixer', this%phase
    !print *, 'qInpResult:', qInpResult
    !print *, 'qDiff:', qDiff

    select case (this%phase)

    case (0)
      if (this%qTot < 0.0_dp) then
        this%qTot = sum(qInpResult)
      end if
      this%qInp0(:) = qInpResult
      this%qOut0(:) = qInpResult + qDiff
      this%qDiff0(:) = qDiff
      this%qDiff0Norm = sqrt(sum(this%qDiff0**2))
      this%vv(:) = this%qDiff0 / this%qDiff0Norm
      qInpResult(:) = this%qInp0 + this%delta * this%vv
      this%phase = 1

    case (1)
      this%qOutPlus(:) = qInpResult + qDiff
      qInpResult(:) = this%qInp0 - this%delta * this%vv
      this%phase = -1

    case (-1)
      qOutMinus = qInpResult + qDiff
      qOutDeriv = (this%qOutPlus - qOutMinus) / (2.0_dp * this%delta)
      uu = qOutDeriv + (1.0_dp - this%cc) / this%cc * this%vv
      qInpResult(:) = this%qInp0 + this%cc * this%qDiff0&
          & + this%cc**2 * this%qDiff0Norm * uu&
          & / (1.0_dp - this%cc * dot_product(this%vv, uu))
      this%phase = 0

    end select

    qInpResult = qInpResult * (this%qTot / sum(qInpResult))

    !print *, 'Exiting mixer'
    !print *, 'qResult:', qInpResult
    !print *, 'Total charge:', sum(qInpResult)

  end subroutine TKernelMixer_mix


end module kernelmixer
