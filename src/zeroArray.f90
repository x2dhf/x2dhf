! ### zeroArray ###
!
!     Zeroize an array.

module zeroArray_m
  implicit none
contains
  subroutine zeroArray (n,array)

    use params

    implicit none

    integer :: i,n
    real (PREC), dimension(*) ::  array

    do i=1,n
       array(i)=0.0_PREC
    enddo

  end subroutine zeroArray
end module zeroArray_m
