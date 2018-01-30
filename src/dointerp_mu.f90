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
subroutine dointerp_mu (nnit,nmumin_p,nmumax_p,nmumin,nmumax,fbefore,fafter) 
  use params
  use discret
  use commons8
  use commons16

  implicit none
  integer :: i,idebug1,idebug2,idebug3,imu,imu_first,imu_p,k, &
       nmu_first,nmu_last,nnit,nmumin_p,nmumax_p,nmumin,nmumax
  real (PREC) :: rerror

  real (PREC16) xmu
  real (PREC16), dimension(9) :: coeffq
  real (PREC16), dimension(9,9) :: coeffq2
  real (PREC), dimension(nnit,*) :: fbefore,fafter
  real (PREC16), external ::  vpoly1q  
  rerror=2
  idebug1=0
  idebug2=0
  idebug3=0
  
  !     interpolation for the first iord2 points in mu variable
  
  do imu=nmumin,nmumax
     if(vmu(imu).ge.vmu_p((nmumin_p-1)+iord2)) goto 115
  enddo
00115 continue
  nmu_first=imu
  
  !     print *,"dointerp_mu: nmumin,nmu_first",nmumin,nmu_first
  
  do k=1,kend
     call lpcoeffq((nmumin_p-1)+iord2+1,k,coeffq)
     do i=1,kend
        coeffq2(i,k)=coeffq(i)
     enddo
  enddo
  do imu=nmumin,nmu_first
     xmu=vmu(imu)            
     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,nmumin_p-1+k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug1.eq.1) then 
           if (abs(fafter(ini,imu)-fbefore(ini,nmumin_p+2)).gt.abs(fbefore(ini,nmumin_p+2))*rerror) then
              write(*,'(2i5,e15.3,4x,3e15.3)') ini,imu,fafter(ini,imu),(fbefore(ini,nmumin_p-1+k),k=1,kend)
           endif
        endif
        
     enddo
  enddo
  
  !     Interpolation for the last iord2 points in mu variable.
  !     Determine the location of the tail region in the new grid
  
  do imu=nmumin,nmumax
     if(vmu(imu).ge.vmu_p(nmumax_p)) goto 120
  enddo
00120 continue
  nmu_last=imu
  if (nmu_last.gt.nmumax) nmu_last=nmumax
  
  do k=1,kend
     call lpcoeffq(nmumax_p-iord2,k,coeffq)
     do i=1,kend
        coeffq2(i,k)=coeffq(i)
     enddo
  enddo
  do imu=nmu_last-iord2+1,nmumax
     xmu=vmu(imu)            
     do ini=1,nnit
        fafter(ini,imu)=0.0_PREC
        do k=1,kend
           ! iord or iord2 (3x)
           fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,nmumax_p-iord+k)*vpoly1q(xmu,coeffq2(1,k))
        enddo
        if (idebug2.eq.1) then 
           if (abs(fafter(ini,imu)-fbefore(ini,nmumax_p-iord+2)).gt. &
                abs(fbefore(ini,nmumax_p-iord+2))*rerror) then
              write(*,'(2i5,e15.3,4x,3e15.3)') ini,imu,fafter(ini,imu), &
                   (fbefore(ini,nmumax_p-iord+k),k=1,kend)
           endif
        endif
     enddo
  enddo
  
  !     interpolation for the inner points in this region
  
  do imu=nmu_first+1,nmu_last-iord2
     xmu=vmu(imu)            
     do i=1,nmumax_p
        if(vmu_p(i).ge.xmu) goto 130
     enddo
00130 continue
     imu_p=i
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
              write(*,'(2i5,e15.3,4x,3e15.3,2i5)') ini,imu,fafter(ini,imu), &
                   (fbefore(ini,imu_p-iord2-1+k),k=1,kend),nmu_first,nmu_last
           endif
        endif
     enddo
  enddo
end subroutine dointerp_mu
