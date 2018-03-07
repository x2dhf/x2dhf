! ***************************************************************************
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### rfun_int ###
!
!     Reads functions from disk in unformatted form and interpolates
!     to a new grid.
!
module rfun_int_m
  implicit none
contains

  subroutine rfun_int (norb_p,cw_orb,cw_coul,cw_exch,wk8,wk16,cw_sctch)
    use params
    use discret
    use scf
    use commons8

    use dointerp_m
    use reada8_m
    use reada16_m

    implicit none

    integer :: i,ica,idel,ierr,ioffset,iorb1,iorb2,ipex,isym,j,k,norb_p

    integer, dimension(60) :: i1b_p,i2b_p,i1e_p,i2e_p,i1si_p,i2si_p,i1ng_p,i2ng_p,i1mu_p,i2mu_p
    integer, dimension(1830) :: i3b_p,i3e_p,i3si_p,i3ng_p,i3mu_p

    real (PREC), dimension(*) :: wk8,cw_orb,cw_coul,cw_exch,cw_sctch

    real (PREC16), dimension(*) :: wk16
    real (PREC), dimension(:), allocatable :: fbefore

    read (iinp11,err=1000) i1b_p,i2b_p,i3b_p,i1e_p,i2e_p,i3e_p, &
         i1si_p,i2si_p,i3si_p,i1ng_p,i2ng_p,i3ng_p,i1mu_p,i2mu_p,i3mu_p

    ioffset=norb-norb_p
    ica=1
    write(iout6,1100)
01100 format(' ... interpolating orbitals:',/,'  ',$)
01110 format(i4,$)

    allocate(fbefore(nni_p*nmu_p(1)))

    do i=1,norb_p
       write(iout6,1110) i
       ipex=mgx(6,i)
       ipex=ipex-2*(ipex/2)
       if (ipex.eq.0) then
          isym= 1
       else
          isym=-1
       endif
       if (lengthfp.eq.8) then
          call reada8(iinp11,i1si_p(i+ioffset),wk8,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading orbital',i
             stop 'rfun 8'
          endif
          do j=1,i1si_p(i+ioffset)
             fbefore(j)=wk8(j)
          enddo
       endif

       if (lengthfp.eq.16) then
          call reada16(iinp11,i1si_p(i+ioffset),wk16,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading orbital',i
             stop 'rfun 16'
          endif
          do j=1,i1si_p(i+ioffset)
             fbefore(j)=wk16(j)
          enddo
       endif

       call dointerp (ica,nmu_p(1),nmu(1),fbefore,cw_orb(i1b(i+ioffset)))

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when reading orbital',i
          stop 'inifun'
       endif
    enddo

    !     retrieve extra data from the orbital input file

    read(iinp11,end=1010,err=1010) area
    read(iinp11,end=1010,err=1010) eng
    read(iinp11,end=1010,err=1010) engo
    read(iinp11,end=1010,err=1010) cmulti
    read(iinp11,end=1010,err=1010) excdi
    read(iinp11,end=1010,err=1010) excqu
    read(iinp11,end=1010,err=1010) excoc
    read(iinp11,end=1010,err=1010) exche
    read(iinp11,end=1010,err=1010) exc5
    read(iinp11,end=1010,err=1010) exc6
    read(iinp11,end=1010,err=1010) exc7
    read(iinp11,end=1010,err=1010) exc8

    write(iout6,*)
    rewind iinp11

    ica=2
    isym=1
    write(iout6,1102)
01102 format(' ... interpolating Coulomb potentials:',/,'  ',$)

    do i=1,norb_p
       write(iout6,1110) i
       if (lengthfp.eq.8) then
          call reada8(iinp12,i2si_p(i+ioffset),wk8,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading coulomb potential',i
             stop 'rfun 8'
          endif
          do j=1,i2si_p(i+ioffset)
             fbefore(j)=wk8(j)
          enddo
       endif
       if (lengthfp.eq.16) then
          call reada16(iinp12,i2si_p(i+ioffset),wk16,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading coulomb potential',i
             stop 'rfun 16'
          endif
          do j=1,i2si_p(i+ioffset)
             fbefore(j)=wk16(j)
          enddo
       endif

       call dointerp (ica,nmu_p(1),nmu(1),fbefore,cw_coul(i2b(i+ioffset)) )

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when reading coulomb potential',i
          stop 'rfun_int'
       endif
    enddo
    write(iout6,*)
    rewind iinp12

    if (imethod.eq.2) then
       write(*,*) 'rfun_int: disk file without extension retrieved '
       write(*,*) '      '
       rewind iinp13
       deallocate(fbefore)
       return
    endif

    if (imethod.eq.1) then
       ica=3
       if (iform.eq.0.or.iform.eq.2) then
          write(*,'("rfun_int: cannot interpolate exchange potentials when they are being retrieved as " &
               & "separate files")')
          stop "rfun_int"
       endif

       write(iout6,1104)
01104  format(' ... interpolating exchange potentials for orbitals:',/,'  ',$)
       do iorb1=1,norb_p
          write(iout6,1110) iorb1

          do iorb2=iorb1,norb_p
             k=iorb1+iorb2*(iorb2-1)/2

             idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
             if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)
             ipex=idel-2*(idel/2)
             if (ipex.eq.0) then
                isym= 1
             else
                isym=-1
             endif

             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50
             if (lengthfp.eq.8) then
                call reada8(iinp13,i3si_p(k),wk8,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 8'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk8(j)
                enddo
             endif
             if (lengthfp.eq.16) then
                call reada16(iinp13,i3si_p(k),wk16,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 16'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk16(j)
                enddo
             endif

             call dointerp (ica,nmu_p(1),nmu(1),fbefore,cw_exch(i3b(k)))
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                stop 'rfun_int'
             endif
             if (iorb1.eq.iorb2) goto 50
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50
             if (lengthfp.eq.8) then
                call reada8(iinp13,i3si_p(k),wk8,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 8'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk8(j)
                enddo
             endif
             if (lengthfp.eq.16) then
                call reada16(iinp13,i3si_p(k),wk16,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 16'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk16(j)
                enddo
             endif

             call dointerp (ica,nmu_p(1),nmu(1),fbefore,cw_exch(i3b(k)+i3si(k)))

             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                stop 'rfun_int'
             endif
00050        continue
          enddo
       enddo
       write(iout6,*)
       write(iout6,*)
       rewind iinp13
    endif

    deallocate(fbefore)
    return
01000 continue
    write(*,*) '... error detected when retrieving the disk file ... '
    stop 'rfun_int'
01010 write(*,*) '... error detected when retrieving extension of the orbital input file ...'
    write(*,*) '... disk file without extension retrieved ... '
    idump=1

  end subroutine rfun_int
end module rfun_int_m
