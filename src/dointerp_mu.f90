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

  do imu=1,nmuall
     xmu=vmu(imu)
     do imu_p=1,nmuall_p-1
        if(vmu_p(imu_p).ge.xmu) exit
     end do

     ! Handle case at beginning or end of array
     if(imu_p .lt. iord2+1) then
        imu_p=iord2+1
     end if
     if(imu_p .gt. nmuall_p-iord2) then
        imu_p=nmuall_p-iord2
     end if

     do k=1,kend
        call lpcoeffq(imu_p,k,coeffq)
        do i=1,kend
           coeffq2(i,k)=coeffq(i)
        enddo
     enddo

     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,imu_p-iord2-1+k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug3.eq.1) then
           if (abs(fafter(ini,imu)-fbefore(ini,imu_p-iord2+1)).gt. &
                abs(fbefore(ini,imu_p-iord2+1))*rerror) then
              write(*,'(2i5,e15.3,4x,5e15.3,2i5)') ini,imu,fafter(ini,imu), &
                   (fbefore(ini,imu_p-iord2-1+k),k=1,kend),nmu_first,nmu_last
           endif
        endif
     enddo
  enddo
end subroutine dointerp_mu
