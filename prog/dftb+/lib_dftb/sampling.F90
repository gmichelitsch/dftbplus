#:include 'common.fypp'

module sampling
  use assert
  use accuracy
  use linkedlist
  use simplealgebra
  use message
  implicit none
  private

  public :: getSuperSampling

contains

  !!* Creates a K-points sampling, equivalent to folding of a reciprocal point
  !!*   of a super lattice.
  !!* @param coeffs Coefficients of the lattice vectors in the linear
  !!*   combination for the super lattice vectors (should be integer values)
  !!* @param shifts Shift of the grid along the three small reciprocal lattice
  !!*   vectors (between 0.0 and 1.0)
  !!* @param latVecs Lattice vector of the original grid
  !!* @param recVecs2p Reciprocal lattice vectors in 2p units (inverse(latVecs))
  !!* @param kPoints Contains the kPoints on exit.
  !!* @param kWeights Contains the weights of the kPoints on exit.
  !!* @param reduceByInversion If points equivalent by inversion should be
  !!*   reduced.
  !!* @desc The routine calculates those reciprocal lattice points of the
  !!*   super lattice, which are inside the Brillouin zone of the original
  !!*   lattice. The resulting points are then all shifted by sum(shift(i)*B(i))
  !!*   where B(i) are the reciprocal lattice vectors of the super lattice.
  !!*   Finally, points equivalent by inversion are reduced, unless specified
  !!*   otherwise.
  subroutine getSuperSampling(coeffs, shifts, kPoints, kWeights, reduceByInversion)
    real(dp), intent(in) :: coeffs(:,:)
    real(dp), intent(in) :: shifts(:)
    real(dp), allocatable, intent(out) :: kPoints(:,:)
    real(dp), allocatable, intent(out) :: kWeights(:)
    logical, intent(in), optional :: reduceByInversion

    real(dp), allocatable :: allKPoints(:,:), allKWeights(:)
    logical, allocatable :: irreducible(:)
    logical :: tReduce
    real(dp) :: invCoeffs(3,3), rr(3)
    integer :: imgRange(2,3), itmp3(3)
    integer :: nAllKPoint, nKPoint
    integer :: i1, i2, i3
    type(listRealR1) :: lr1

    real(dp), parameter :: tol = 1e-4_dp
    real(dp), parameter :: minLim = -tol, maxLim = 1.0_dp - tol

    @:ASSERT(all(shape(coeffs) == (/ 3, 3 /)))
    ! check they are integers
    @:ASSERT(all(coeffs - nint(coeffs) < epsilon(1.0_dp)))
    @:ASSERT(size(shifts) == 3)

    if (present(reduceByInversion)) then
      tReduce = reduceByInversion
    else
      tReduce = .true.
    end if

    ! Get the eight corners of the original (big) reciprocal unit cell as linear
    ! combination of the reciprocal lattice vectors (B) of the superlattice
    ! Note: b = B * N^T (b/B: rec.lat.vec. of lattice/superlattice)
    imgRange(:,:) = 0
    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          itmp3 = i1*nint(coeffs(1,:)) + i2*nint(coeffs(2,:)) &
              & + i3*nint(coeffs(3,:))
          imgRange(1,:) = min(itmp3, imgRange(1,:))
          imgRange(2,:) = max(itmp3, imgRange(2,:))
        end do
      end do
    end do
    ! Decrease by one to have the range [min, max)
    imgRange(2,:) = imgRange(2,:) - 1

    ! invCoeffs = (N^-1)^T
    call invert33(invCoeffs, coeffs)
    invCoeffs = transpose(invCoeffs)
    call init(lr1)

    do i1 = imgRange(1, 1), imgRange(2, 1)
      do i2 = imgRange(1, 2), imgRange(2, 2)
        do i3 = imgRange(1, 3), imgRange(2, 3)
          ! relative coordinate with respect to the original reciprocal lattice
          rr(:) = matmul(invCoeffs, real((/ i1, i2, i3 /), dp))
          if (all(rr >= minLim) .and. all(rr < maxLim)) then
            ! Add point + shift vector
            call append(lr1, rr + matmul(invCoeffs, shifts))
          end if
        end do
      end do
    end do

    nAllKPoint = len(lr1)
    if (abs(real(nAllKPoint,dp) - abs(determinant33(coeffs))) > tol) then
      call error("Monkhorst-Pack routine failed to find all K-points.")
    end if

    allocate(allKPoints(3, nAllKPoint))
    allocate(allKWeights(nAllKPoint))
    call asArray(lr1, allKPoints)
    call destruct(lr1)
    allKPoints = modulo(allKPoints, 1.0_dp)
    allKWeights = 1.0_dp / real(nAllKPoint, dp)

    ! Reduce by inversion if needed
    if (tReduce) then
      allocate(irreducible(nAllKPoint))
      irreducible(:) = .true.
      do i1 = 1, nAllKPoint
        if (.not. irreducible(i1)) then
          cycle
        end if
        rr(:) = modulo(-1.0_dp * allKPoints(:,i1), 1.0_dp)
        do i2 = i1 + 1, nAllKPoint
          if (.not. irreducible(i2)) then
            cycle
          end if
          if (all(abs(allKPoints(:,i2) - rr(:)) < tol)) then
            irreducible(i2) = .false.
            allKWeights(i1) = allKWeights(i1) + allKWeights(i2)
          end if
        end do
      end do
      nKPoint = count(irreducible)
      allocate(kPoints(3, nKPoint))
      allocate(kWeights(nKPoint))
      i1 = 1
      i2 = 1
      do while (i2 <= nKpoint)
        if (irreducible(i1)) then
          kPoints(:,i2) = allKPoints(:,i1)
          kWeights(i2) = allKWeights(i1)
          i2 = i2 + 1
        end if
        i1 = i1 + 1
      end do
    else
      allocate(kPoints(3, nAllKPoint))
      allocate(kWeights(nAllKPoint))
      kPoints(:,:) = allKPoints
      kWeights(:) = allKWeights
    end if

  end subroutine getSuperSampling


end module sampling
