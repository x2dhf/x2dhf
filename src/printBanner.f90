! ### banner ###  
!
!     Prints a banner of the program
! 

subroutine printBanner
  implicit none
  
  write(*,1020)
  
01020 format(//,25x,'FINITE DIFFERENCE 2D HARTREE-FOCK  (version 2.3)',//)
  
end subroutine printBanner



