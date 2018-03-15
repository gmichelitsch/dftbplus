!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides interface to the BML and Progress libraries.
module bmlprogress
  use bml
  implicit none
  private

  public :: bml_matrix_t
  public :: bml_set_row, bml_get_row
  public :: bml_transpose_triangle, bml_adjungate_triangle

end module bmlprogress
