
! ### vpolyq ###
!
!     This function uses the Horner scheme to calculate value of the polynomial
!     stored in array a at a particular point

function vpolyq (mup,a)
  use params
  use commons16

  implicit none
  integer :: i,mup
  real (PREC16) :: vpolyq
  real (PREC16) :: x
  real (PREC16), dimension(kend) :: a

  x=vmuq(mup)
  vpolyq=0.0_PREC16
  do i=iord,2,-1
     vpolyq=(vpolyq+a(i))*x
  enddo
  vpolyq=vpolyq+a(1)

end function vpolyq
