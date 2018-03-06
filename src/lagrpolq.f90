! ### lagrpolq ###
!
subroutine lagrpolq (dc1,dc2)
  use params
  use discret
  use commons8
  use commons16
  use vpolyq_m
  use lpcoeffq_m

  implicit none

  integer :: i,ib,ie,ig,imup,k,mup

  real (PREC16), dimension(kbeg:kend) :: coeffq,coeff1,coeff2
  real (PREC), dimension(kbeg:kend,7,kbeg:kend) :: dc1,dc2
  real(PREC16), dimension(maxmu) :: vmuq
  
!      derivative coeff. from 8th-order Lagrange interpolation formula
!      are stored in dc[1-2](ngbound,imu,k), where ngbound is the number of
!      the grid boundaries (1 for 1-2, 2 for 2-3 etc), imu is one of the 7
!      intergrid mu values and k points to a derivative coefficient.

!      the grid boundaries (1 for 1-2, 2 for 2-3 etc), imu is one of the 7
!      intergrid mu values and k points to a derivative coefficient.

!      For each value of mup (p) a Lagrange polynomial of the 8th-order
!      employing 4 grid points to the left and right of imu (imu is included)
!      is build (x_{0} corresponds to the grid boundary mu value)

!      f^{(p)}(x) = \sum_{k=p-4}^{p+4} f(x_{k})
!                   \prod_{i=p-4,i\ne k}^{p+4} {(x-x_{i}) \over (x_{k}-x_{i})}

!      Polynomials for each k have to be calculated and their first and second
!      derivatives obtained. The values of the derivative polynomials at grid
!      point x_{p} (mup) are the differentiation coefficients sought.

!      in order to quarantee the real*8 in the Lagrange coefficients
!      the quadruple precision have to be used during their generation

!      iord=9
!      kbeg=1
!      kend=9

!  FIXME if possible use a function for projecting real*8 into real*16 variables

  do i=1,mxnmu
     !        vmuq(i)=_QFLOAT_(vmu(i))
     vmuq(i)=(vmu(i))
  enddo

  !     loop over grid boundaries first (ngrids>1)

  do ig=1,ngrids-1
     ib=iemu(ig)-3
     ie=iemu(ig)+3
     !        loop over mup values (mup=ib...ie)
     !        loop over k
     imup=0
     do mup=ib,ie
        imup=imup+1
        do k=kbeg,kend
           call lpcoeffq(mup,k,coeffq,vmuq)
           call lpderq(coeffq,coeff1,coeff2)
           dc1(ig,imup,k)=vpolyq(mup,coeff1,vmuq)
           dc2(ig,imup,k)=vpolyq(mup,coeff2,vmuq)
        enddo
     enddo
  enddo

end subroutine lagrpolq
