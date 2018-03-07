! ### banner ###
!
!     Prints a banner of the program

module printBanner_m
  implicit none
contains
  subroutine printBanner
    implicit none

    write(*,1020)

01020 format(//,25x,'FINITE DIFFERENCE 2D HARTREE-FOCK  (version 2.1)',//)

  end subroutine printBanner
end module printBanner_m
