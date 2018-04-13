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
    integer :: i,imu,ini_p,k,nmuall

    real (PREC), dimension(nni_p,nmuall) :: fbefore
    real (PREC), dimension(nni,nmuall) :: fafter
    real (PREC16) xni
    real (PREC16), dimension(kend) :: coeffq
    real (PREC16), dimension(kend,kend) :: coeffq2
    real (PREC16), dimension(maxmu) ::  vmuq

    do ini=1,nni
       xni=vni(ini)
       do ini_p=1,nni_p-1
          if(vni_p(ini_p).ge.xni) exit
       enddo

       ! Handle case at beginning or end of array
       if(ini_p .lt. iord2+1) then
          ini_p=iord2+1
       end if
       if(ini_p .gt. nni_p-iord2) then
          ini_p=nni_p-iord2
       end if

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
