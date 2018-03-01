! ### setPrecision ###

!     Calculates floating-point precision and integer/real variable
!     lengths

subroutine setPrecision

  use commons8

  implicit none

  integer :: i,i4tmp1,i8tmp1,i8tmp2
  real (PREC) :: o

  !     pi from bc -l
  !     scale=32
  !     a(1)*4=           3.14159265358979323846264338327948
  !     scale=64
  !     a(1)*4=           3.14159265358979323846264338327950
  
  pii=atan(one)*four
  
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
