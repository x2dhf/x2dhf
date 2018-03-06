! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### dointerp_mu ###
!
!     Performs interpolations of functions in mu variable.
!
subroutine dointerp_mu (nnit,nmuall_p,nmuall,fbefore,fafter)
  use params
  use discret
  use commons8
  use commons16

  implicit none
  integer :: i,idebug1,idebug2,idebug3,imu,imu_p,k, &
       nmu_first,nmu_last,nnit,nmuall_p,nmuall
  real (PREC) :: rerror

  real (PREC16) xmu
  real (PREC16), dimension(kend) :: coeffq
  real (PREC16), dimension(kend,kend) :: coeffq2
  real (PREC), dimension(nnit,nmuall_p) :: fbefore
  real (PREC), dimension(nnit,nmuall) :: fafter
  real (PREC16), external ::  vpoly1q
  rerror=2
  idebug1=0
  idebug2=0
  idebug3=0

  !     interpolation for the first iord2 points in mu variable
  do nmu_first=1,nmuall-1
     if(vmu(nmu_first).ge.vmu_p(iord2)) exit
  enddo

  do k=1,kend
     call lpcoeffq(iord2+1,k,coeffq)
     do i=1,kend
        coeffq2(i,k)=coeffq(i)
     enddo

  enddo
  do imu=1,nmu_first
     xmu=vmu(imu)
     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug1.eq.1) then
           if (abs(fafter(ini,imu)-fbefore(ini,3)).gt.abs(fbefore(ini,3))*rerror) then
              write(*,*) 'first'
              write(*,'(2i5,e15.3,4x,5e15.3)') ini,imu,fafter(ini,imu),(fbefore(ini,k),k=1,kend)
           endif
        endif

     enddo
  enddo

  !     Interpolation for the last iord2 points in mu variable.
  !     Determine the location of the tail region in the new grid

  do nmu_last=1,nmuall-1
     if(vmu(nmu_last).ge.vmu_p(nmuall_p)) exit
  enddo

  do k=1,kend
     call lpcoeffq(nmuall_p-iord2,k,coeffq)
     do i=1,kend
        coeffq2(i,k)=coeffq(i)
     enddo
  enddo
  do imu=nmu_last-iord2+1,nmuall
     xmu=vmu(imu)
     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,nmuall_p-iord2+k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug2.eq.1) then
           if (abs(fafter(ini,imu)-fbefore(ini,nmuall_p-iord2+2)).gt. &
                abs(fbefore(ini,nmuall_p-iord2+2))*rerror) then
              write(*,*) 'last'
              write(*,'(2i5,e15.3,4x,5e15.3)') ini,imu,fafter(ini,imu), &
                   (fbefore(ini,nmuall_p-iord2+k),k=1,kend)
           endif
        endif
     enddo
  enddo

  !     interpolation for the inner points in this region

  do imu=nmu_first+1,nmu_last-iord2-1
     xmu=vmu(imu)
     do imu_p=1,nmuall_p-1
        if(vmu_p(imu_p).ge.xmu) exit
     enddo

     do k=1,kend
        call lpcoeffq(imu_p,k,coeffq)
        do i=1,kend
           coeffq2(i,k)=coeffq(i)
        enddo
     enddo

     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           ! Array overflow in the second index of fbefore:
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,imu_p-iord2-1+k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug3.eq.1) then
           if (abs(fafter(ini,imu)-fbefore(ini,imu_p-iord2+1)).gt. &
                abs(fbefore(ini,imu_p-iord2+1))*rerror) then
              write(*,*) 'inside'
              write(*,'(2i5,e15.3,4x,5e15.3,2i5)') ini,imu,fafter(ini,imu), &
                   (fbefore(ini,imu_p-iord2-1+k),k=1,kend),nmu_first,nmu_last
           endif
        endif
     enddo
  enddo
end subroutine dointerp_mu
