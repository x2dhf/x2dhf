! ### testn2f ###
!   Tests nabla^2 f

subroutine testn2f(f,wk0,wk1,wk2,wk3)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ii,itest,j
  real (PREC) :: x,z
  real (PREC), dimension(*) :: f,wk0,wk1,wk2,wk3

  itest=1
  
  !   take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
  if (itest.eq.1) then
     do i=1,mxnmu
        ii=(i-1)*nni
        do j=1,nni
           !             t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
           x=(r/2.0_PREC)*vxi1(i)*veta1(j)
           z=(r/2.0_PREC)*vxi(i)*veta(j)
           f(ii+j)=1.0_PREC
           !            write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
        enddo
     enddo
     
     call n2f(f,wk0,wk1,wk2,wk3)
     call pmtx(nni,mxnmu,wk3,ione,ione,incrni,incrmu)
     
  endif
  
end subroutine testn2f
