subroutine nuclder1(iorb,psi,dznu,dzmu,dernu,dermu,wk2,wk3,wk4,wk5)
  use params
  use discret
  use solver
  use commons8

  implicit none
  integer :: i,i1beg,imu,inu,iorb,iprt,nmut,ngorb,ngpot

  real (PREC), dimension(nni,*) :: psi,dermu,dernu
  real (PREC), dimension(*) :: wk2,wk3,wk4,wk5
  real (PREC), dimension(6) :: dzmu
  real (PREC), dimension(5000) :: dznu

  i1beg=i1b(iorb)
  nmut=i1mu(iorb)
  ngorb=i1si(iorb)
  ngpot=i2si(iorb)	

  !   interpolate psi at (nni,1)

  psi(nni,1)= exeven(1)*psi(nni-1,1)+exeven(2)*psi(nni-2,1)+exeven(3)*psi(nni-3,1) &
       +exeven(4)*psi(nni-4,1)+exeven(5)*psi(nni-5,1)
 
  call putin (nni,nmut,isym,psi,wk2)	
  call diffnu (nmut,wk2,wk3,wk4,wk5)
  call putout (nni,nmut,dernu,wk3)

  call diffmu (nmut,wk2,wk3)
  call putout (nni,nmut,dermu,wk3)

  dernu(nni,1)= exeven(1)*dernu(nni-1,1)+exeven(2)*dernu(nni-2,1)+exeven(3)*dernu(nni-3,1) &
       +exeven(4)*dernu(nni-4,1)+exeven(5)*dernu(nni-5,1)

  dermu(nni,1)= exeven(1)*dermu(nni,2)+exeven(2)*dermu(nni,3)+exeven(3)*dermu(nni,4)       &
       +exeven(4)*dermu(nni,5)+exeven(5)*dermu(nni,6)

  !  having derivative over nu calculate the derivative over z
  !  near (0,0,-R/2)

  do inu=nni-5,nni
     dznu(inu)=(-1.0_PREC)*2.0_PREC/(r*veta1(inu))*dernu(inu,1)
  enddo

  dzmu(1)=.0_PREC
  do imu=2,6
     dzmu(imu)=(-1.0_PREC)*2.0_PREC/(r*vxi1(imu))*dermu(nni,imu)
  enddo
  
  !   calculate the derivative over z at (0,0,R/2+) by extrapolation

  dznu(nni)= exeven(1)*dznu(nni-1)+exeven(2)*dznu(nni-2)+exeven(3)*dznu(nni-3) &
       +exeven(4)*dznu(nni-4)+exeven(5)*dznu(nni-5)

  !   calculate the derivative over z at (0,0,ER/2-) by extrapolation

  dzmu(1)= exeven(1)*dzmu(2)+exeven(2)*dzmu(3)+exeven(3)*dzmu(4)+exeven(4)*dzmu(5)+exeven(5)*dzmu(6)

! FIXME
  iprt=1
  if (iprt.eq.0) then        
     write(*,*) 'psi(nni-5,...nni,1)'
     write(*,2000) (psi(i,1),i=nni-5,nni)
     write(*,*) 'dznu(nni-5,nni)'
     write(*,2000) (dznu(i),i=nni-5,nni)
     write(*,*) 'psi(5,..1)'
     write(*,2000) (psi(nni,i),i=6,1,-1)
     write(*,*) 'dzmu(6,...,1)'
     write(*,2000) (dzmu(i),i=6,1,-1)
02000 format(1x,3e25.12/1x,3e25.12)
  endif

  dznu(nni)=occ(iorb)*2.0_PREC*psi(nni,1)*dznu(nni)
  dzmu(1)  =occ(iorb)*2.0_PREC*psi(nni,1)*dzmu(1)

end subroutine nuclder1



