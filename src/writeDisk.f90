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
! ### writeDisk ###

!     Writes orbitals, potenials, Lagrange multipliers (diagonal
!     and off-diagonal) and multipole expansion coefficients to a disk
!     file in either formatted or unformatted form

module writeDisk_m
  implicit none
contains
  subroutine writeDisk (cw_orb,cw_coul,cw_exch)
    use params
    use discret
    use scf
    use commons8

    use getdatetime_m
    use wtdisknat_m
    use wtdisk32_m
    use wtdisk32f_m
    use wtdisk64_m
    use wtdisk64f_m
    use wtdisk128_m
    use wtdisk128f_m

    implicit none
    character*80 :: datetime
    integer :: i4,i8,i16,li4tmp,li8tmp

    integer, allocatable :: i4tmp(:)
    integer*8, allocatable :: i8tmp(:)

    real (PREC), allocatable :: r8l60(:),r8l1200(:),r8l3660(:),r8mxsize(:)
    real (PREC), allocatable :: r16l60(:),r16l1200(:),r16l3660(:),r16mxsize(:)
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    data i4,i8,i16/4,8,16/

    ! fixed precision
    rewind(iout24)
    call getDateTime(datetime)
    write(iout24,'(a80)') header
    write(iout24,'(a80)') datetime
    write(iout24,formint) ngrids,nni,nmu
    write(iout24,formfp64) r,rgrid
    write(iout24,formfp64) z1,z2
    write(iout24,formint) norb,nel,nexch

    call flush(iout24)

    if (inout32.eq.1) then
       li4tmp=10*maxorb+5*maxorb*(maxorb+1)/2

       allocate(i4tmp(li4tmp))
       allocate(r8l60(maxorb))
       allocate(r8l1200(20*maxorb))
       allocate(r8l3660(maxorb*(maxorb+1)))
       allocate(r8mxsize(nni*mxnmu))
       call wtdisk32 (cw_orb,cw_coul,cw_exch,i4tmp,r8l60,r8l1200,r8l3660,r8mxsize)

       deallocate(i4tmp)
       deallocate(r8l60)
       deallocate(r8l1200)
       deallocate(r8l3660)
       deallocate(r8mxsize)

    elseif (inout32.eq.2) then
       li4tmp=10*maxorb+5*maxorb*(maxorb+1)/2
       allocate(i4tmp(li4tmp))
       allocate(r8l60(maxorb))
       allocate(r8l1200(20*maxorb))
       allocate(r8l3660(maxorb*(maxorb+1)))
       allocate(r8mxsize(nni*mxnmu))
       call wtdisk32f (cw_orb,cw_coul,cw_exch,i4tmp,r8l60,r8l1200,r8l3660,r8mxsize)

       deallocate(i4tmp)
       deallocate(r8l60)
       deallocate(r8l1200)
       deallocate(r8l3660)
       deallocate(r8mxsize)

    elseif (inout64.eq.1) then
       li8tmp=10*maxorb+5*maxorb*(maxorb+1)/2
       allocate(i8tmp(li8tmp))
       allocate(r8l60(maxorb))
       allocate(r8l1200(20*maxorb))
       allocate(r8l3660(maxorb*(maxorb+1)))
       allocate(r8mxsize(nni*mxnmu))
       call wtdisk64 (cw_orb,cw_coul,cw_exch,i8tmp,r8l60,r8l1200,r8l3660,r8mxsize)

       deallocate(i8tmp)
       deallocate(r8l60)
       deallocate(r8l1200)
       deallocate(r8l3660)
       deallocate(r8mxsize)

    elseif (inout64.eq.2) then
       li8tmp=10*maxorb+5*maxorb*(maxorb+1)/2
       allocate(i8tmp(li8tmp))
       allocate(r8l60(maxorb))
       allocate(r8l1200(20*maxorb))
       allocate(r8l3660(maxorb*(maxorb+1)))
       allocate(r8mxsize(nni*mxnmu))
       call wtdisk64f (cw_orb,cw_coul,cw_exch,i8tmp,r8l60,r8l1200,r8l3660,r8mxsize)

       deallocate(i8tmp)
       deallocate(r8l60)
       deallocate(r8l1200)
       deallocate(r8l3660)
       deallocate(r8mxsize)

    elseif (inout128.eq.1) then
       li8tmp=10*maxorb+5*maxorb*(maxorb+1)/2
       allocate(i8tmp(li8tmp))
       allocate(r16l60(maxorb))
       allocate(r16l1200(20*maxorb))
       allocate(r16l3660(maxorb*(maxorb+1)))
       allocate(r16mxsize(nni*mxnmu))

       call wtdisk128 (cw_orb,cw_coul,cw_exch,i8tmp,r16l60,r16l1200,r16l3660,r16mxsize)

       deallocate(i8tmp)
       deallocate(r16l60)
       deallocate(r16l1200)
       deallocate(r16l3660)
       deallocate(r16mxsize)


    elseif (inout128.eq.2) then
       li8tmp=10*maxorb+5*maxorb*(maxorb+1)/2
       allocate(i8tmp(li8tmp))
       allocate(r16l60(maxorb))
       allocate(r16l1200(20*maxorb))
       allocate(r16l3660(maxorb*(maxorb+1)))
       allocate(r16mxsize(nni*mxnmu))

       call wtdisk128f (cw_orb,cw_coul,cw_exch,i8tmp,r16l60,r16l1200,r16l3660,r16mxsize)

       deallocate(i8tmp)
       deallocate(r16l60)
       deallocate(r16l1200)
       deallocate(r16l3660)
       deallocate(r16mxsize)
    else
       call wtdisknat(cw_orb,cw_coul,cw_exch)
    endif

  end subroutine writeDisk
end module writeDisk_m
