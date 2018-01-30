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
! ### fockform ###
!
!     Print the Fock equation for every orbital.

subroutine fockform 
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,iorb,iorb1,ipc,kex,idexp,norb2,orborder		
  real (PREC) :: coo,coo0,coo1,coo2,w
  parameter (orborder=1)

  norb2=norb*norb		

  ! exchange contributions due to different pair of orbitals
  ! Ka:  sigma i - nonsigma j
  ! Kb:  nonsigma i - nonsigma i 
  ! Kc:  nonsigma i - nonsigma j 

  if (nel.eq.1) return

  if (norb.eq.1) then
     
     !      FIXME what about one non-sigma orbital            
     !      contribution from coulomb potential of one sigma orbital
     !     call copy (ngpot,pot(ibpot),ione,wk1,ione)
     write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gut(iorb)
     write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
          1.0,iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb)
     return
  endif

  if (orborder.eq.1) goto 1000

  do iorb=1,norb
     write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gut(iorb)

     !      add contributions from Coulomb potentials
        
     do iorb1=1,norb
        kex=iorb+norb*(iorb1-1)
        ipc=iorb1+norb*(iorb-1)
           
        if (iorb.le.iorb1) then
           idexp=iorb+iorb1*(iorb1-1)/2
        else
           idexp=iorb1+iorb*(iorb-1)/2
        endif
        
        coo=occ(iorb1)
        if (iorb.eq.iorb1) coo=coo-1.0_PREC

        !           write(*,'(2i5," J1  " ,1Pe12.2)') iorb,iorb1,coo
        write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
             coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1)
     enddo

     !      add contributions from exchange potentials 
     do iorb1=1,norb
        kex=iorb+norb*(iorb1-1)
        ipc=iorb1+norb*(iorb-1)
        
        if (iorb.le.iorb1) then
           idexp=iorb+iorb1*(iorb1-1)/2
        else
           idexp=iorb1+iorb*(iorb-1)/2
        endif
        
        if (iorb1.ne.iorb)  then
           coo=gec(kex)
           !              write(*,'(2i5," Ka  " ,1Pe12.2)') iorb,iorb1,-coo
           write(*,'(15x,f6.2, "  Ka(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
           
           if (ilc(idexp).gt.1) then
              coo=gec(kex+norb2)
              !              call axpy (ngexp,coo,excp(ibexp+ngexp),ione,wk2,ione)
              !                 write(*,'(2i5," Kc  " ,1Pe12.2)') iorb,iorb1,-coo
              write(*,'(15x,f6.2, "  Kc(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                   -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
              
           endif
        else
           if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
              coo=gec(kex)
              coo2=coo
              !                 write(*,'(2i5," Kb  " ,1Pe12.2)') iorb,iorb1,-coo
              write(*,'(15x,f6.2, "  Kb(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                   -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
           endif
        endif
     enddo
  enddo

1000 continue

  do iorb=norb,1,-1
     write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gut(iorb)

     !      add contributions from Coulomb potentials
        
     do iorb1=norb,1,-1
        kex=iorb+norb*(iorb1-1)
        ipc=iorb1+norb*(iorb-1)
           
        if (iorb.le.iorb1) then
           idexp=iorb+iorb1*(iorb1-1)/2
        else
           idexp=iorb1+iorb*(iorb-1)/2
        endif
        
        coo=occ(iorb1)
        
        !         FIXME               
        if (iorb.eq.iorb1) coo=coo-1.0_PREC
        coo0=coo
        
        !        call axpy (ngpot1,coo,pot(ibpot1),ione,wk1,ione)
        
        !           write(*,'(2i5," J1  " ,1Pe12.2)') iorb,iorb1,coo
        write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
             coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1)
     enddo

     !      add contributions from exchange potentials 
     do iorb1=norb,1,-1
        kex=iorb+norb*(iorb1-1)
        ipc=iorb1+norb*(iorb-1)
        
        if (iorb.le.iorb1) then
           idexp=iorb+iorb1*(iorb1-1)/2
        else
           idexp=iorb1+iorb*(iorb-1)/2
        endif
        
        coo=occ(iorb1)
        
        if (iorb.eq.iorb1) coo=coo-1.0_PREC
        coo0=coo
        
        coo1=0.0_PREC
        coo2=0.0_PREC
        
        if (iorb1.ne.iorb)  then
           coo=gec(kex)
           coo1=coo
           !           call copy (ngexp,excp(ibexp),ione,wk2,ione)
           !           call scal (ngexp,coo,wk2,ione)
           !              write(*,'(2i5," Ka  " ,1Pe12.2)') iorb,iorb1,-coo
           write(*,'(15x,f6.2, "  Ka(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
           
           
           if (ilc(idexp).gt.1) then
              coo=gec(kex+norb2)
              coo2=coo
              !              call axpy (ngexp,coo,excp(ibexp+ngexp),ione,wk2,ione)
              !                 write(*,'(2i5," Kc  " ,1Pe12.2)') iorb,iorb1,-coo
              write(*,'(15x,f6.2, "  Kc(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                   -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
              
           endif
        else
           if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
              coo=gec(kex)
              coo2=coo
              !                 write(*,'(2i5," Kb  " ,1Pe12.2)') iorb,iorb1,-coo
              write(*,'(15x,f6.2, "  Kb(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                   -coo,iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb),bond(iorb),gut(iorb)
           endif
        endif
     enddo
  enddo



  return
end subroutine fockform
