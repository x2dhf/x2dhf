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
! ### wtdisk128 ###

!     Writes orbitals, potenials, Lagrange multipliers and multipole
!     expansion coefficients to a disk file in an unformatted form
!     (integers are 4 and reals are 8-byte long)

module wtdisk128_m
  implicit none
contains
  subroutine wtdisk128 (cw_orb,cw_coul,cw_exch,i8tmp,r16l60,r16l1200,r16l3660,r16mxsize)
    use params
    use discret
    use memory
    use scf
    use commons8

    use writea128_m
    use wtdexch1_m

    implicit none
    integer :: i,iorb1,iorb2,ioffset,j,k
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    integer*8 :: i8tmp1,ierr8
    integer*8, dimension(*) :: i8tmp
    real (PREC16), dimension(*) :: r16l60,r16l1200,r16l3660,r16mxsize

    write(iout6,1050)

    ioffset=0
    do i=1,maxorb
       i8tmp(ioffset+i)=i1b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp(ioffset+i)=i2b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp(ioffset+i)=i3b(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp(ioffset+i)=i1e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp(ioffset+i)=i2e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp(ioffset+i)=i3e(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp(ioffset+i)=mxsize
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp(ioffset+i)=mxsize
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp(ioffset+i)=mxsize
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp(ioffset+i)=i1ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp(ioffset+i)=i2ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp(ioffset+i)=i3ng(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp(ioffset+i)=i1mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp(ioffset+i)=i2mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp(ioffset+i)=i3mu(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    write (i8out21,err=1000) (i8tmp(i),i=1,ioffset)

    !   add orbitals
    do i=1,norb,1
       i8tmp1=mxsize
       do j=1,mxsize
          r16mxsize(j)=cw_orb(i1b(i)+j-1)
       enddo
       call writea128(i8out21,i8tmp1,r16mxsize,ierr8)
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisk128'
       endif
    enddo

    ! FIXME
    !   append the following arrays to the output orbital file
    do i=1,maxorb
       r16l60(i)=area(i)
    enddo
    !    write(i8out21,err=1020) area
    write(i8out21,err=1020) (r16l60(i),i=1,maxorb)

    do i=1,maxorb
       r16l60(i)=eng(i)
    enddo
    !    write(i8out21,err=1020) eng
    write(i8out21,err=1020) (r16l60(i),i=1,maxorb)

    do i=1,maxorb*maxorb
       r16l3660(i)=engo(i)
    enddo
    !    write(i8out21,err=1020) engo
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*maxorb)

    do i=1,20*maxorb
       r16l1200(i)=cmulti(i)
    enddo
    !    write(i8out21,err=1020) cmulti
    write(i8out21,err=1020) (r16l1200(i),i=1,20*maxorb)

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=excdi(i)
    enddo
    !    write(i8out21,err=1020) excdi
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=excqu(i)
    enddo
    !    write(i8out21,err=1020) excqu
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=excoc(i)
    enddo
    !    write(i8out21,err=1020) excoc
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=exche(i)
    enddo
    !    write(i8out21,err=1020) exche
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=exc5(i)
    enddo
    !    write(i8out21,err=1020) exc5
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=exc6(i)
    enddo
    !    write(i8out21,err=1020) exc6
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=exc7(i)
    enddo
    !    write(i8out21,err=1020) exc7
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    do i=1,maxorb*(maxorb+1)
       r16l3660(i)=exc8(i)
    enddo
    !    write(i8out21,err=1020) exc8
    write(i8out21,err=1020) (r16l3660(i),i=1,maxorb*(maxorb+1))

    !   write out Coulomb potentials

    do i=1,norb
       do j=1,mxsize
          r16mxsize(j)=cw_coul(i2b(i)+j-1)
       enddo
       i8tmp1=mxsize
       call writea128(i8out22,i8tmp1,r16mxsize,ierr8)

       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisk128'
       endif
    enddo

    !   write out exchange potentials

    if (imethod.eq.1) then
       if (iform.eq.2.or.iform.eq.3) then

          do iorb1=1,norb
             do iorb2=iorb1,norb
                k=iorb1+iorb2*(iorb2-1)/2
                if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50

                do j=1,i3si(k)
                   r16mxsize(j)=cw_exch(i3b(k)+j-1)
                enddo
                i8tmp1=i3si(k)
                call writea128(i8out23,i8tmp1,r16mxsize,ierr8)
                !                call writea(i8out23,i3si(k),cw_exch(i3b(k)),ierr)

                if (ierr8.ne.0) then
                   write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                   stop 'wtdisk128'
                endif
                if (iorb1.eq.iorb2) goto 50
                if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50

                do j=1,i3si(k)
                   r16mxsize(j)=cw_exch(i3b(k)+i3si(k)+j-1)
                enddo
                i8tmp1=i3si(k)
                call writea128(i8out23,i8tmp1,r16mxsize,ierr8)
                if (ierr8.ne.0) then
                   write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                   stop 'wtdisk128'
                endif
50              continue
             enddo
          enddo
          rewind(i8out23)

          !         if iform=0 exchange potentials do not have to be saved since
          !         they are continually updated during the scf process

       elseif (iform.eq.0.or.iform.eq.1) then
          if (iform.eq.1) then
             ! FIXME
             call wtdexch1(cw_exch)
          endif
       endif

    elseif (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) then
       do j=1,i3si(1)
          r16mxsize(j)=cw_exch(length3-mxsize+j-1)
       enddo
       i8tmp1=i3si(1)
       call writea128(i8out23,i8tmp1,r16mxsize,ierr8)
       !       call writea(i8out23,i3si(1),cw_exch(i3b(1)),ierr)

       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisk128'
       endif
    endif

    rewind(i8out21)
    rewind(i8out22)
    rewind(i8out23)

    return

1000 continue
    write(iout6,1070)
    stop 'wtdisk128'

    ! FIXME
1020 continue
    write(*,*) 'wtdiski128: error encountered when writing an extension to orbital file'
    write(*,*) '      '
    stop 'wtdiski128'
1050 format(1x,'... writing functions to disk (r128) ...')
1070 format(//1x,'error! can not write data to disk'//)
    !1080 format(//1x,'error! can not write data to disk'//)
  end subroutine wtdisk128
end module wtdisk128_m
