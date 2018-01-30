
! ### tden ###

!    This routine initializes arrays used for differentiation needed in nuclder

subroutine tden(iorb,ngorb1,psi,wk2)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iborb1,iorb,iorb1,ngorb1
  real (PREC) :: coo
  real (PREC), dimension(*) :: psi,wk2

  do i=1,ngorb1
     wk2(i)=.0_PREC
  enddo

  call extinorg(psi)

  !      do iorb1=1,norb
  do iorb1=iorb,iorb
     iborb1=i1b(iorb1)	   	 
     ngorb1=i1si(iorb1)	 	 
     coo=occ(iorb1)
     do i=1,ngorb1
        wk2(i)=wk2(i)+coo*psi(iborb1+i-1)*psi(iborb1+i-1)
     enddo
  enddo
  
end subroutine tden
