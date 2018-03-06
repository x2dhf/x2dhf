! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### dointerp_nu ###
!
!     Performs interpolations of functions in ni variable.
!
module dointerp_nu_m
  implicit none
contains
  subroutine dointerp_nu (nmuall,fbefore,fafter,vmuq)
    use params
    use discret
    use commons8
    use commons16
    use lpcoeffq_m
    use vpoly1q_m

    implicit none
    integer :: i,imu,ini_p,k,nmuall,nni_first,nni_last

    real (PREC), dimension(nni_p,nmuall) :: fbefore
    real (PREC), dimension(nni,nmuall) :: fafter
    real (PREC16) xni
    real (PREC16), dimension(kend) :: coeffq
    real (PREC16), dimension(kend,kend) :: coeffq2
    real (PREC16), dimension(maxmu) ::  vmuq

    ! interpolation for the first iord2 points in ni variable
    do nni_first=1,nni-1
       if(vni(nni_first).ge.vni_p(iord2)) exit
    enddo

    do k=1,kend
       call lpcoeffq(iord2+1,k,coeffq,vmuq)
       do i=1,kend
          coeffq2(i,k)=coeffq(i)
       enddo
    enddo

    do ini=1,nni_first
       xni=vni(ini)
       do imu=1,nmuall
          fafter(ini,imu)=0.0_PREC
          do k=1,kend
             fafter(ini,imu)=fafter(ini,imu)+fbefore(k,imu)*vpoly1q(xni,coeffq2(1,k))
          enddo
          if (ini.eq.1.and.(imu.eq.1.or.imu.eq.nmuall)) then
             fafter(ini,imu)=fbefore(ini,imu)
          endif
       enddo
    enddo

    !  Interpolation for the last iord2 points in ni variable.
    !  Determine the location of the tail region in the new grid

    do nni_last=1,nni-1
       if(vni(nni_last).ge.vni_p(nni_p-iord2+1)) exit
    enddo

    do k=1,kend
       call lpcoeffq(nni_p-iord2,k,coeffq,vmuq)
       do i=1,kend
          coeffq2(i,k)=coeffq(i)
       enddo
    enddo

    do ini=nni_last-iord2+1,nni
       xni=vni(ini)
       do imu=1,nmuall
          fafter(ini,imu)=0.0_PREC
          do k=1,kend
             fafter(ini,imu)=fafter(ini,imu)+fbefore(nni_p-iord+k,imu)*vpoly1q(xni,coeffq2(1,k))
          enddo
          if (ini.eq.nni.and.(imu.eq.1.or.imu.eq.nmuall)) then
             fafter(nni,imu)=fbefore(nni_p,imu)
          endif
       enddo
    enddo

    !     interpolation for the inner points in this region

    do ini=nni_first+1,nni_last-iord2-1
       xni=vni(ini)
       do ini_p=1,nni_p-1
          if(vni_p(ini_p).ge.xni) exit
       enddo

       do k=1,kend
          call lpcoeffq(ini_p,k,coeffq,vmuq)
          do i=1,kend
             coeffq2(i,k)=coeffq(i)
          enddo
       enddo

       do imu=1,nmuall
          fafter(ini,imu)=0.0_PREC
          do k=1,kend
             fafter(ini,imu)=fafter(ini,imu)+fbefore(ini_p-iord2-1+k,imu)*vpoly1q(xni,coeffq2(1,k))
          enddo
       enddo
    enddo
  end subroutine dointerp_nu
end module dointerp_nu_m
