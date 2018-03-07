! ### setPrecision ###

!     Calculates floating-point precision and integer/real variable
!     lengths

module setPrecision_m
  implicit none
contains
  subroutine setPrecision
    use commons8
    implicit none

    integer :: i,i4tmp1,i8tmp1,i8tmp2
    real (PREC) :: o

    !     Warning!

    !     Precision calculated by the code below depends on compiler
    !     optimization flags used. This routine must be compiled with the
    !     minimum optimization on x86 systems to avoid reporting precision
    !     of the extended IEEE 754 arithmetic.

    !     When in doubt try using function DLAMCH('e') (LAPACK auxiliary routine)

    precis=one
    do i=1,10000
       precis=precis/two
       o=one+precis
       if ( o.eq.one ) goto 100
    enddo
100 continue

    !     set default lengths of integer and real constants and variables used

    i8tmp1=2**30
    i8tmp2=100000*i8tmp1
    i4tmp1=i8tmp2
    if (i4tmp1.eq.0) then
       lengthint=4
    else
       lengthint=8
    endif

    if (precis.gt.1e-20_PREC) then
       lengthfp =8
    else
       lengthfp =16
    endif

  end subroutine setPrecision
end module setPrecision_m
