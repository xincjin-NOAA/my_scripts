! Minimal stand-in for GSI's kinds module (GSI builds with r_kind = double,
! i_kind = 4-byte integer). Only what phil*.f90 and hilbert_curve.f90 use.
module kinds
  implicit none
  integer, parameter :: i_kind = selected_int_kind(8)
  integer, parameter :: sp     = selected_real_kind(6)
  integer, parameter :: dp     = selected_real_kind(15)
  integer, parameter :: r_kind = dp
end module kinds

! Minimal stand-in for GSI's constants module.
module constants
  use kinds, only: r_kind
  implicit none
  real(r_kind), parameter :: zero = 0.0_r_kind, one = 1.0_r_kind, &
                             two  = 2.0_r_kind, four = 4.0_r_kind
end module constants
