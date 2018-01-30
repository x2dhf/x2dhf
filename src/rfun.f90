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
! c ### rfun ###	

!     Reads functions from a disk file in an unformatted form

subroutine rfun (norbt,cw_orb,cw_coul,cw_exch,wk8,wk16)
  use params
  use discret
  use memory
  use scf
  use commons8

  implicit none
  integer :: i,ierr,ii,ioffset,iorb1,iorb2,j,k,norbt
  integer, dimension(9750) :: i4tmp
  integer*8, dimension(9750) :: i8tmp

  real (PREC), dimension(*) :: wk8,cw_orb,cw_coul,cw_exch
  real (PREC), dimension(60) :: r8tmp1
  real (PREC), dimension(1200) :: r8tmp2
  real (PREC), dimension(3660) :: r8tmp3

  real (PREC16), dimension(*) :: wk16
  real (PREC16), dimension(60) :: r16tmp1
  real (PREC16), dimension(1200) :: r16tmp2
  real (PREC16), dimension(3660) :: r16tmp3

  if (inpform.eq.0) then
     if (lengthint.eq.4) then
        read (iinp11,err=1000) i4tmp
     else
        read (iinp11,err=1000) i8tmp
     endif
  else
     read (iinp11,formint,err=1000) i4tmp
  endif

  
  !     if there are no virtual orbitals defined retrieve all orbitals and
  !     Coulomb functions from the corresponding disk data files ignoring
  !     inhyd flags. The flags are restored when returning from this
  !     routine.
  
  !     ioffset is no longer used
  !      ioffset=norb-norbt
  ioffset=0
  
  if (ini4.eq.0) then
     do i=1,norb
        inhyd(i)=0
     enddo
  endif
  
  !     retrieve orbitals 
  
  !      do i=1,norbt
  !     
  do i=1,norb
     if (inhyd(i).eq.1) goto 100
     
     if (inpform.eq.0) then
        if (lengthfp.eq.8) then
           call reada8(iinp11,i1si(i),wk8,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading orbital',i
!              stop 'rfun 8'
              write(iout6,*) 'continue with crossed fingers ...'
           endif
           do j=1,i1si(i)
              cw_orb(i1b(i)+j-1)=wk8(j)
           enddo
        endif
        
        if (lengthfp.eq.16) then
           call reada16(iinp11,i1si(i),wk16,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading orbital',i
              stop 'rfun 16'
           endif
           do j=1,i1si(i)
              cw_orb(i1b(i)+j-1)=wk16(j)
           enddo
        endif
     else
        read (iinp11,formfp,err=1004) (wk8(ii),ii=1,i1si(i))
        !            call readaf(iinp11,i1si(i),wk8,ierr)
        !            if (ierr.ne.0) then
        !               write(iout6,*) 'error detected when reading orbital',i
        !               stop 'rfun 8'
        !            endif				   	
        do j=1,i1si(i)
           cw_orb(i1b(i)+j-1)=wk8(j)
        enddo
     endif
100  continue
  enddo
  
  !     retrieve the extra data from the orbital input file
  
  if (lengthfp.eq.8) then
     if (inpform.eq.0) then
        read(iinp11,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,area)
        
        read(iinp11,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,eng)
        
        read(iinp11,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,engo)
        
        read(iinp11,end=1008,err=1010) r8tmp2
        do i=1,1200
           cmulti(i)=r8tmp2(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           excdi(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           excqu(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           excoc(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           exche(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc5(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc6(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc7(i)=r8tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc8(i)=r8tmp3(i)
        enddo
     else
        read(iinp11,formfp64,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,area)

        read(iinp11,formfp64,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,eng)
        
        read(iinp11,formfp64,end=1008,err=1010) r8tmp1
        call rfunaux(r8tmp1,engo)

        read(iinp11,formfp64,end=1008,err=1010) r8tmp2
        do i=1,1200
           cmulti(i)=r8tmp2(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           excdi(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           excqu(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           excoc(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           exche(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc5(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc6(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc7(i)=r8tmp3(i)
        enddo
        read(iinp11,formfp64,end=1008,err=1010) r8tmp3
        do i=1,3660
           exc8(i)=r8tmp3(i)
        enddo
     endif
  else
     if (inpform.eq.0) then
        read(iinp11,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,area)

        read(iinp11,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,eng)

        read(iinp11,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,engo)

        read(iinp11,end=1008,err=1010) r16tmp2
        do i=1,1200
           cmulti(i)=r16tmp2(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           excdi(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           excqu(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           excoc(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           exche(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc5(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc6(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc7(i)=r16tmp3(i)
        enddo
        read(iinp11,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc8(i)=r16tmp3(i)
        enddo
     else
        read(iinp11,formfp128,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,area)

        read(iinp11,formfp128,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,eng)
        
        read(iinp11,formfp128,end=1008,err=1010) r16tmp1
        call rfunaux16(r16tmp1,engo)

        read(iinp11,formfp128,end=1008,err=1010) r16tmp2
        do i=1,1200
           cmulti(i)=r16tmp2(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           excdi(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           excqu(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           excoc(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           exche(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc5(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc6(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc7(i)=r16tmp3(i)
        enddo

        read(iinp11,formfp128,end=1008,err=1010) r16tmp3
        do i=1,3660
           exc8(i)=r16tmp3(i)
        enddo
     endif
  endif
  rewind iinp11	
  
  !     read in Coulomb potentials 
  
  do i=1,norb
     if (inhyd(i).eq.1) goto 200 
     if (inpform.eq.0) then
        if (lengthfp.eq.8) then
           call reada8(iinp12,i2si(i),wk8,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading coulomb potential',i
              stop 'rfun 8'
           endif
           do j=1,i2si(i)
              cw_coul(i2b(i)+j-1)=wk8(j)
           enddo
        endif
        if (lengthfp.eq.16) then
           call reada16(iinp12,i2si(i),wk16,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading coulomb potential',i
              stop 'rfun 16'
           endif
           do j=1,i2si(i)
              cw_coul(i2b(i)+j-1)=wk16(j)
           enddo
        endif
     else
        read (iinp12,formfp,err=1005) (wk8(ii),ii=1,i2si(i))
        !            call readaf(iinp12,i2si(i),wk8,ierr)
        do j=1,i2si(i)
           cw_coul(i2b(i)+j-1)=wk8(j)
        enddo
     endif
200  continue
  enddo
  
  !     read in exchange potentials from the exhange potential input file
  !     (only for HF calculations)
  
  if (imethod.eq.1) then
     if (iform.eq.1.or.iform.eq.3.and.ini.ne.6) then
        do iorb1=1,norbt
           do iorb2=iorb1,norbt
              k=iorb1+iorb2*(iorb2-1)/2
              if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50
              
              if (inpform.eq.0) then                  
                 if (lengthfp.eq.8) then
                    call reada8(iinp13,i3si(k),wk8,ierr)
                    if (ierr.ne.0) then
                       write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                       stop 'rfun 8'
                    endif
                    do j=1,i3si(k)
                       cw_exch(i3b(k)+j-1)=wk8(j)
                    enddo
                 endif
                 if (lengthfp.eq.16) then
                    call reada16(iinp13,i3si(k),wk16,ierr)
                    if (ierr.ne.0) then
                       write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                       stop 'rfun 16'
                    endif
                    do j=1,i3si(k)
                       cw_exch(i3b(k)+j-1)=wk16(j)
                    enddo
                 endif
              else
                 read (iinp13,formfp,err=1010) (wk8(j),j=1,i3si(k))
                 do j=1,i3si(k)
                    cw_exch(i3b(k)+j-1)=wk8(j)
                 enddo
              endif
              
              if (iorb1.eq.iorb2) goto 50
              if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50
              
              if (inpform.eq.0) then                  
                 if (lengthfp.eq.8) then
                    call reada8(iinp13,i3si(k),wk8,ierr)
                    if (ierr.ne.0) then
                       write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                       stop 'rfun 8'
                    endif
                    do j=1,i3si(k)
                       cw_exch(i3b(k)+i3si(k)+j-1)=wk8(j)
                    enddo
                 endif
                 if (lengthfp.eq.16) then
                    call reada16(iinp13,i3si(k),wk16,ierr)
                    if (ierr.ne.0) then
                       write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                       stop 'rfun 16'
                    endif
                    do j=1,i3si(k)
                       cw_exch(i3b(k)+i3si(k)+j-1)=wk16(j)
                    enddo
                 endif
              else
                 read(iinp13,formfp,err=1010) (wk8(j),j=1,i3si(k))
                 do j=1,i3si(k)
                    cw_exch(i3b(k)+i3si(k)+j-1)=wk8(j)
                 enddo
              endif
50            continue
           enddo
        enddo
        rewind iinp13
     endif
     write(*,*) '... orbitals and potentials retrieved ... '
  endif
  
  if (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) then
     if (inpform.eq.0) then                  
        if (lengthfp.eq.8) then
           call reada8(iinp13,i3si(1),wk8,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading local exchange potential'
              stop 'rfun 8'
           endif
           do j=1,i3si(1)
              cw_exch(length3-mxsize+j-1)=wk8(j)
           enddo
        endif
        if (lengthfp.eq.16) then
           call reada16(iinp13,i3si(1),wk16,ierr)
           if (ierr.ne.0) then
              write(iout6,*) 'error detected when reading local exchange potential'
              stop 'rfun 16'
           endif
           do j=1,i3si(1)
              cw_exch(length3-mxsize+j-1)=wk16(j)
           enddo
        endif
     else
        read (iinp13,formfp,err=1010) (wk8(j),j=1,i3si(1))
        do j=1,i3si(1)
           cw_exch(length3-mxsize+j-1)=wk8(j)
        enddo
     endif
  endif
  
  
  if (lengthfpin.eq.8) then
     formfp=formfp64
  else
     formfp=formfp128
  endif
  
  if(iprint(115).ne.0) then
     write(*,*)
     write(*,*) 'rfun: input orbitals'
     !        print orbitals 
     do i=1,norb
        write(*,*) iorn(i),' ',bond(i)
        call pmtx(nni,i1mu(i),cw_orb(i1b(i)),ione,ione,incrni,incrmu)
     enddo
  endif
  
  if(iprint(116).ne.0) then
     write(*,*)
     write(*,*) 'rfun: input Coulomb potentials'
     !        print Coulomb potential
     do i=1,norb
        write(*,*) iorn(i),' ',bond(i)
        call pmtx(nni,i1mu(i),cw_coul(i2b(i)),ione,ione,incrni,incrmu)
     enddo
  endif
  
  if(iprint(117).ne.0.and.nexch.ge.1) then
     write(*,*)
     write(*,*) 'rfun: input exchange potentials'
     !        print exchange potentials
     do i=1,nexch
        write(*,*) 'exchange potential ',i 
        call pmtx(nni,i1mu(i),cw_exch(i3b(i)),ione,ione,incrni,incrmu)
     enddo
  endif
  
  do i=1,norb
     inhyd(i)=inhydlcao(i)
  enddo
  
  
  return
  
1000 continue
  write(*,*) '... error detected when retrieving the disk file ... '
  stop 'rfun'
  
  
1004 continue
  write(*,*) '... error detected when retrieving orbital function', i
  stop 'rfun'
  
1005 continue
  write(*,*) '... error detected when retrieving potential function', i
  stop 'rfun'
  
1008 continue
  write(*,*) '... end of file detected when retrieving extension of the orbital input file ...'
  write(*,*) '... disk file without extension retrieved ... '
  idump=1
  goto 1100
  
1010 write(*,*) '... error detected when retrieving extension of the orbital input file ...'
  write(*,*) '... disk file without extension retrieved ... '
  idump=1
1100 continue
  
end subroutine rfun

