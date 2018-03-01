! ### setmaxunit ###
!
!     Determines and sets maximum unit number allowed

subroutine setmaxunit

  use commons8

  implicit none

  integer :: i,iunit
  do i=99,1000 
     iunit=i
     open(iunit,status='scratch',form='unformatted',err=100)
     close(iunit)
  enddo
  maxunit=1000
  return
100 maxunit=i-1

end subroutine setmaxunit
