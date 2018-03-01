! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### potAsympt ### 
!
!     determines the asymptotic values of Coulomb and exchange
!     potentials needed for a given Fock equation 

subroutine potAsympt (iorb1,pot,excp)
  use params
  use commons8

  implicit none
  integer :: i3beg,ibeg,idel,iorb1,iorb2,k,ngrid

  real (PREC), dimension(*) :: pot,excp

!   determine values of Coulomb potentials in the asymptotic region
!   from the multipole expansion

  if (itouch(iorb1).eq.1) then	
     ibeg=i2b(iorb1)
     call coulAsympt(iorb1,pot(ibeg))
  endif
  
  if (islat.eq.1) return
  
  !   initialize amul array needed to evaluate exchange potential in
  !   the asymptotic region by calling zasyx
  
  do iorb2=1,norb
     if (itouch(iorb1).eq.0.and.itouch(iorb2).eq.0) goto 10
     k=i3xk(iorb1,iorb2)
     !      check if k=0
     if (k.eq.0) goto 10
     ngrid=i3si(k)
     if (ilc(k).ne.0) then
        
        !         the exchange potential address array should already be modified by
        !         call to rfdexch routine
        
        i3beg=i3b(k)
        idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
        if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)
        
        call exchAsympt (idel,k,excp(i3beg))
        
        if (ilc(k).eq.2) then
           idel=mgx(6,iorb2)+mgx(6,iorb1)
           k=k+norb*(norb+1)/2
           i3beg=i3beg+ngrid
           call exchAsympt (idel,k,excp(i3beg))
        endif
     endif
10   continue
  enddo
  
end subroutine potAsympt
