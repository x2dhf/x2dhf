! ### lpderq ###
!
!     This routine calculates coefficients of the first and second
!     derivative of the polynomial stored in a

module lpderq_m
  implicit none
contains
  subroutine lpderq (a,a1,a2)

    use params
    use commons16

    implicit none

    integer :: i
    real (PREC16), dimension(9) :: a,a1,a2

    do i=1,iord-1
       a1(i)=a(i+1)*dble(i)
    enddo
    a1(iord)=0.0_PREC16

    do i=1,iord-2
       a2(i)=a1(i+1)*dble(i)
    enddo
    a2(iord-1)=0.0_PREC16
    a2(iord)  =0.0_PREC16

  end subroutine lpderq
end module lpderq_m
