! ### testnfng ###
!     Tests nabla f nabla g

subroutine testnfng(rhot,grhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

  use params
  use discret
  use commons8

  implicit none
  integer :: i,ii,itest,j
  real (PREC) :: t1,z
  real (PREC), dimension(*) ::  wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,grhot

  itest=1
  
  if (itest.eq.1) then
     do i=1,mxnmu
        ii=(i-1)*nni
        do j=1,nni
           rhot(ii+j) =0.d0
           grhot(ii+j)=0.d0
        enddo
     enddo
     
     call nfng (rhot,grhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
     call pmtx(nni,mxnmu,wk7,ione,ione,incrni,incrmu)
     
  endif
  
end subroutine testnfng
