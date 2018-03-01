!  ### lpcoeffq ###
! 
!      This routine calculates coefficients of the (sub)Lagrange 
!      polynomial for a grid point k
!      \prod_{i=1,i\ne k}^{9} {(x-x_{i}) \over (x_{k}-x_{i})}
! 
!      x_{k}= vmu(mup-5+k), k=1,...,9 
! 
!      Last modification: 2.01.01
!      
subroutine lpcoeffq (mup,k,coeffq)
  use params
  use commons16

  implicit none

  integer :: i,ib,ic1,ic2,j,k,mup

  real (PREC16) :: c1,denom

  real (PREC16), dimension(9) :: coeffq
  real (PREC16), dimension(12) :: a,b

!  calculate denominator product

  denom=1.0_PREC16
  ib=mup-(iord/2+1)
  do i=kbeg,kend
     if (i.ne.k) denom=denom*(vmuq(ib+k)-vmuq(ib+i))
  enddo
  
  do i=1,12
     a(i)=0.0_PREC16
     b(i)=0.0_PREC16
  enddo
  
  !     calculate nominator product: 
  !     a(1)+a(2)*x+a(3)*x^2+...+a(9)*x^8
  !     storing coefficients in a and b
  
  !     if (k.ne.kbeg) then
  !     a(1)=-vmuq(ib+kbeg)
  !     else
  !     a(1)=-vmuq(ib+kbeg+1)	  
  !     endif
  !     a(2)=1.0_PREC	
  
  a(1)=1.0_PREC16
  
  !  multiply polynomial a by (x+c1), c1=-vmuq(ib+i)
  
  ic2=1
  do i=kbeg,kend
     if (i.ne.k) then
        ic2=ic2+1
        c1=-vmuq(ib+i)
        b(1)=a(1)*c1
        do ic1=2,ic2
           b(ic1)=a(ic1)*c1+a(ic1-1)
        enddo
        b(ic2+1)=a(ic2)
        do j=1,ic2+1
           a(j)=b(j)
        enddo
     endif
  enddo
  !     
  do i=1,iord
     coeffq(i)=a(i)/denom
  enddo

end subroutine lpcoeffq

