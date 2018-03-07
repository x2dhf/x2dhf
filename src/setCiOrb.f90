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
! ### setCiOrb ###

! sets Ci symmetry of a given orbital: the values from the [0,pi/2] c region replace
! those from [pi/2,pi] one

module setCiOrb_m
  implicit none
contains
  subroutine setCiOrb (n2,nmut,orb,ihsym)
    use params
    use discret
    use commons8

    implicit none
    integer :: ihsym,mi,mis,n2,ni,nis,nmut
    real (PREC), dimension(nni,nmut) :: orb

    if (iprint(105).eq.0) then
       !       set symmetry according to input data
       if (ihsym.eq.1) then
          do ni=1,n2-1
             do mi=1,nmut
                orb(nni+1-ni,mi)= orb(ni,mi)
             enddo
          enddo
       else
          do ni=1,n2-1
             do mi=1,nmut
                orb(nni+1-ni,mi)=-orb(ni,mi)
             enddo
          enddo
          do mi=1,nmut
             orb(n2,mi)=0.0_PREC
          enddo
       endif
    else

       mis=10
       nis=2

       ! check symmetry at a given point and impose the same symmetry at other ones

       if (orb(nis,mis)*orb(nni-nis+1,mis).gt.0.0_PREC) then
          do ni=1,n2-1
             do mi=1,nmut
                orb(nni+1-ni,mi)= orb(ni,mi)
             enddo
          enddo
          ihsym= 1
       else
          do ni=1,n2-1
             do mi=1,nmut
                orb(nni+1-ni,mi)=-orb(ni,mi)
             enddo
          enddo

          do mi=1,nmut
             orb(n2,mi)=0.0_PREC
          enddo
          ihsym=-1
       endif
    endif

  end subroutine setCiOrb
end module setCiOrb_m
