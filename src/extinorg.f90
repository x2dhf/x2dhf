! ### extinorg ###
subroutine extinorg(psi)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iprt

  real (PREC), dimension(nni,*) :: psi
  
  iprt=1
  if (iprt.eq.0) then        
     write(*,*) 'psi(nni-5..,1) ',(psi(i,1),i=nni-5,nni)
     write(*,*) 'psi(nni,1..)   ',(psi(nni,i),i=1,6)
  endif
  
  psi(1,1)= exeven(1)*psi(2,1)+exeven(2)*psi(3,1)+exeven(3)*psi(4,1)+exeven(4)*psi(5,1)+exeven(5)*psi(6,1)

  psi(nni,1)= exeven(1)*psi(nni-1,1)+exeven(2)*psi(nni-2,1)+exeven(3)*psi(nni-3,1) &
       +exeven(4)*psi(nni-4,1)+exeven(5)*psi(nni-5,1)

  iprt=1
  if (iprt.eq.0) then        
     write(*,*)
     write(*,*) 'psi(nni-5..,1) ',(psi(i,1),i=nni-5,nni)
     write(*,*) 'psi(nni,1..)   ',(psi(nni,i),i=1,6)
  endif

end subroutine extinorg


