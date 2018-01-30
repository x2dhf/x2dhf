! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996,2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### wtdisknat ### 

!     Writes orbitals, potenials, Lagrange multipliers (diagonal and
!     off-diagonal) and multipole expansion coefficients to a disk file
!     in an unformatted form

subroutine wtdisknat (cw_orb,cw_coul,cw_exch)
  use params
  use discret
  use memory
  use scf
  use commons8

  implicit none
  integer :: i,ierr,iorb1,iorb2,k

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

  write(iout6,1050)
1050 format(1x,'... writing functions to disk (current precision) ...')
  
  !   write a header into the orbital output file (deprecated)

  if (inpform.eq.0) then
     write (iout21,err=1000) i1b,i2b,i3b,i1e,i2e,i3e,i1si,i2si,i3si,i1ng,i2ng,i3ng,i1mu,i2mu,i3mu	
  else
     write (iout21,formint,err=1005) i1b,i2b,i3b,i1e,i2e,i3e,i1si,i2si,i3si,i1ng,i2ng,i3ng,i1mu,i2mu,i3mu	
  endif

  !   add orbitals
  
  do i=1,norb
     call writea(iout21,i1si(i),cw_orb(i1b(i)),ierr)
     if (ierr.ne.0) then
        write(iout6,*) 'error detected when writing orbital',i
        stop 'wtdisknat'
     endif
  enddo
  
  !   append the following arrays to the output orbital file 
  write(iout21,err=1020) area
  write(iout21,err=1020) eng
  write(iout21,err=1020) engo
  write(iout21,err=1020) cmulti
  write(iout21,err=1020) excdi
  write(iout21,err=1020) excqu
  write(iout21,err=1020) excoc
  write(iout21,err=1020) exche
  write(iout21,err=1020) exc5
  write(iout21,err=1020) exc6
  write(iout21,err=1020) exc7
  write(iout21,err=1020) exc8
  
  !   write out Coulomb potentials
  do i=1,norb
     call writea(iout22,i2si(i),cw_coul(i2b(i)),ierr)
     if (ierr.ne.0) then
        write(iout6,*) 'error detected when writing coulomb potential',i
        stop 'wtdisknat'
     endif
  enddo
  
  !   write out exchange potentials
  if (imethod.eq.1) then  
     if (iform.eq.2.or.iform.eq.3) then
        do iorb1=1,norb
           do iorb2=iorb1,norb
              k=iorb1+iorb2*(iorb2-1)/2
              if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50
              call writea(iout23,i3si(k),cw_exch(i3b(k)),ierr)
              if (ierr.ne.0) then
                 write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                 stop 'wtdisknat'
              endif
              if (iorb1.eq.iorb2) goto 50
              if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50
              call writea(iout23,i3si(k),cw_exch(i3b(k)+i3si(k)),ierr)
              if (ierr.ne.0) then
                 write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                 stop 'wtdisknat'
              endif
50            continue
           enddo
        enddo
        rewind(iout23)
        
        !         if iform=0 exchange potentials do not have to be saved since 
        !         they are continually updated during the scf process
        
     elseif (iform.eq.0.or.iform.eq.1) then
        if (iform.eq.1) then
           call wtdexch1(cw_exch)
        endif
     endif
     
  elseif (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) then
     !       call writea(iout23,i3si(1),cw_exch(i3b(1)),ierr)         
     !       call writea(iout23,i3si(1),cw_exch(length3-mxsize+1),ierr)
     call writea(iout23,mxsize,cw_exch(length3-mxsize),ierr)
     
     if (ierr.ne.0) then
        write(iout6,*) 'error detected when writing local exchange potential'
        stop 'wtdisknat'
     endif
  endif
  
  rewind(iout21)
  rewind(iout22)
  rewind(iout23)
  
  
  if (iprint(121).ne.0) then
     write(*,*)
     write(*,*) 'wtdisknat: output orbitals'
     !      print orbitals 
     do i=1,norb
        write(*,*) iorn(i),' ',bond(i)
        call pmtx(nni,i1mu(i),cw_orb(i1b(i)),ione,ione,incrni,incrmu)
     enddo
  endif
  
  if(iprint(122).ne.0) then
     write(*,*)
     write(*,*) 'wtdisknat: output Coulomb potentials'
     !      print Coulomb potential
     do i=1,norb
        write(*,*) iorn(i),' ',bond(i)
        call pmtx(nni,i1mu(i),cw_coul(i2b(i)),ione,ione,incrni,incrmu)
     enddo
  endif
  
  !   print exchange potentials
  if(iprint(123).ne.0.and.nexch.ge.1) then
     if (imethod.eq.1) then  
        write(*,*)
        write(*,*) 'wtdisknat: output exchange potentials'
        
        do i=1,nexch
           write(*,*) 'exchange potential ',i 
           call pmtx(nni,i1mu(i),cw_exch(i3b(i)),ione,ione,incrni,incrmu)
        enddo
     elseif (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) then
        write(*,*)
        write(*,*) 'wtdisknat: output local exchange potential'
        call pmtx(nni,i1mu(1),cw_exch(length3-mxsize),ione,ione,incrni,incrmu)
     endif
  endif
  
  return

1000 continue
  write(*,*) 'wtdisknat: error encountered when writing integer arrays to disk (unformatted form)'
  stop 'wtdisknat'	
  
  
1005 continue
  write(*,*) 'wtdisknat: error encountered when writing integer arrays to disk (formatted form)'
  stop 'wtdisknat'	
      
1020 continue
  write(*,*) 'wtdisknat: error encountered when writing an extra real array to disk'
  write(*,*) '      '
  stop 'wtdisknat'
  
end subroutine wtdisknat








