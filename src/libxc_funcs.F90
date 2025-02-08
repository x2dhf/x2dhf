module xc_f03_funcs_m
  use, intrinsic :: iso_c_binding
  implicit none

  public
       ! List of functionals
#include "libxc_inc.f90"

  ! These are old names kept for compatibility
  integer(c_int), parameter, public :: &
    XC_LDA_X_1D             =  21,     &
    XC_GGA_X_BGCP           =  38,     &
    XC_GGA_C_BGCP           =  39,     &
    XC_GGA_C_BCGP           =  39,     &
    XC_GGA_C_VPBE           =  83,     &
    XC_GGA_XC_LB            = 160,     &
    XC_MGGA_C_CC06          = 229,     &
    XC_GGA_K_ABSR1          = 506,     &
    XC_GGA_K_ABSR2          = 507,     &
    XC_LDA_C_LP_A           = 547,     &
    XC_LDA_C_LP_B           = 548,     &
    XC_MGGA_C_LP90          = 564

end module xc_f03_funcs_m
