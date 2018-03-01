! ### vpoly1q ###
!  
!     This function uses the Horner scheme to calculate value of the polynomial 
!     stored in array a at a particular point
!
function vpoly1q (x,a)
  use params
  use commons16

  implicit none
  integer :: i
  real (PREC16) :: vpoly1q
  real (PREC16) :: x
  real (PREC16), dimension(9) :: a

  vpoly1q=0.0_PREC
  do i=iord,2,-1
     vpoly1q=(vpoly1q+a(i))*x
  enddo
  vpoly1q=vpoly1q+a(1)
  return
end function vpoly1q
