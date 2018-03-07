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
! ### wtdisk32f ###

!     Writes orbitals, potenials, Lagrange multipliers and multipole
!     expansion coefficients to a disk file in a formatted form:
!     integers are 4 and reals are 8-byte long.

module wtdisk32f_m
  implicit none
contains
  subroutine wtdisk32f (cw_orb,cw_coul,cw_exch,i4tmp,r8l60,r8l1200,r8l3660,r8mxsize)
    use params
    use discret
    use memory
    use scf
    use commons8

    use writea32f_m
    use wtdexch1f_m

    implicit none

    integer :: i,ierr,iorb1,iorb2,ioffset,j,k

    integer :: i4tmp1
    integer, dimension(*) :: i4tmp
    real (PREC), dimension(*) :: r8l60,r8l1200,r8l3660,r8mxsize
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    write(iout6,1050)

    ioffset=0
    do i=1,maxorb
       i4tmp(ioffset+i)=i1b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp(ioffset+i)=i2b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp(ioffset+i)=i3b(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp(ioffset+i)=i1e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp(ioffset+i)=i2e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp(ioffset+i)=i3e(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp(ioffset+i)=i1si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp(ioffset+i)=i2si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp(ioffset+i)=i3si(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp(ioffset+i)=i1ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp(ioffset+i)=i2ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp(ioffset+i)=i3ng(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp(ioffset+i)=i1mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp(ioffset+i)=i2mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp(ioffset+i)=i3mu(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    write (iout21,formint,err=1000) (i4tmp(i),i=1,ioffset)

    !   add orbitals
    do i=1,norb
       i4tmp1=i1si(i)
       do j=1,i4tmp1
          r8mxsize(j)=cw_orb(i1b(i)+j-1)
       enddo
       call writea32f(iout21,i4tmp1,r8mxsize,ierr)
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisk32f'
       endif
    enddo

    !   append the following arrays to the output orbital file
    do i=1,maxorb
       r8l60(i)=area(i)
    enddo
    !   write(iout21,err=1020) area
    write(iout21,formfp64,err=1020) (r8l60(i),i=1,maxorb)

    do i=1,maxorb
       r8l60(i)=eng(i)
    enddo
    !   write(iout21,err=1020) eng
    write(iout21,formfp64,err=1020) (r8l60(i),i=1,maxorb)


    do i=1,maxorb*maxorb
       r8l3660(i)=engo(i)
    enddo
    !   write(iout21,err=1020) engo
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*maxorb)

    do i=1,20*maxorb
       r8l1200(i)=cmulti(i)
    enddo
    !    write(iout21,err=1020) cmulti
    write(iout21,formfp64,err=1020) (r8l1200(i),i=1,20*maxorb)


    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=excdi(i)
    enddo
    !   write(iout21,err=1020) excdi
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=excqu(i)
    enddo
    !   write(iout21,err=1020) excqu
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=excoc(i)
    enddo
    !    write(iout21,err=1020) excoc
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=exche(i)
    enddo
    !   write(iout21,err=1020) exche
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=exc5(i)
    enddo
    !   write(iout21,err=1020) exc5
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=exc6(i)
    enddo
    !   write(iout21,err=1020) exc6
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=exc7(i)
    enddo
    !   write(iout21,err=1020) exc7
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r8l3660(i)=exc8(i)
    enddo
    !   write(iout21,err=1020) exc8
    write(iout21,formfp64,err=1020) (r8l3660(i),i=1,maxorb*(maxorb+1))

    !   write out Coulomb potentials

    do i=1,norb
       i4tmp1=i2si(i)
       do j=1,i4tmp1
          r8mxsize(j)=cw_coul(i2b(i)+j-1)
       enddo
       call writea32f(iout22,i4tmp1,r8mxsize,ierr)
       !       call writea(iout22,i2si(i),cw_coul(i2b(i)),ierr)

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisk32f'
       endif
    enddo

    !   write out exchange potentials
    if (imethod.eq.1) then
       if (iform.eq.2.or.iform.eq.3) then
          do iorb1=1,norb
             do iorb2=iorb1,norb
                k=iorb1+iorb2*(iorb2-1)/2
                if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50

                i4tmp1=i3si(k)
                do j=1,i4tmp1
                   r8mxsize(j)=cw_exch(i3b(k)+j-1)
                enddo
                call writea32f(iout23,i4tmp1,r8mxsize,ierr)
                !                call writea(iout23,i3si(k),cw_exch(i3b(k)),ierr)

                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                   stop 'wtdisk32f'
                endif

                if (iorb1.eq.iorb2) goto 50
                if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50

                i4tmp1=i3si(k)
                do j=1,i4tmp1
                   r8mxsize(j)=cw_exch(i3b(k)+i3si(k)+j-1)
                enddo
                call writea32f(iout23,i4tmp1,r8mxsize,ierr)

                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                   stop 'wtdisk32f'
                endif

50              continue
             enddo
          enddo
          rewind(iout23)

          !         if iform=0 exchange potentials do not have to be saved since
          !         they are continually updated during the scf process

       elseif (iform.eq.0.or.iform.eq.1) then
          if (iform.eq.1) then
             ! FIXME
             call wtdexch1f(cw_exch)
          endif
       endif
    elseif (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) then
       i4tmp1=i3si(1)
       do j=1,i4tmp1
          r8mxsize(j)=cw_exch(length3-mxsize+j-1)
       enddo
       call writea32f(iout23,i4tmp1,r8mxsize,ierr)

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisk32f'
       endif
    endif

    rewind(iout21)
    rewind(iout22)
    rewind(iout23)

    return

1000 continue
    write(iout6,1070)
    stop 'wtdisk32f'

1020 continue
    write(*,*) 'wtdisk32: error encountered when writing an extension to orbital file'
    write(*,*) '      '
    stop 'wtdisk32'
1050 format(1x,'... writing functions to disk (i32f) ...')
1070 format(//1x,'error! can not write data to disk'//)
    !1080 format(//1x,'error! can not write data to disk'//)
  end subroutine wtdisk32f
end module wtdisk32f_m
