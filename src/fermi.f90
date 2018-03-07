! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fermi ###

!     Analyzes Fermi charge distribution.

module fermi_m
  implicit none
contains
  subroutine fermi
    use params
    use discret
    use commons8

    implicit none

    integer :: iprt,mu,nc1,nc2,ni
    real (PREC) :: a,at3,atw,c,facto1,fmtoau,rr,rrms,rrmsfm,t,tfm,xr2,xrr

    if (iprint(230).ne.0) then
       write(*,*) ' '
       write(*,*) '  finite nuclei with the Fermi charge distribution'
    endif

    !     calculate parameters of the finite nucleus charge distribution
    !     (Fremi distribution)

    fmtoau=1.0e-13_PREC/ainfcm

    if (z1.ne.0.0_PREC) then

       !        set atomic weight for centre A

       atw=z1atmass
       at3=atw**(1.0_PREC/3.0_PREC)
       rrmsfm = 0.8360_PREC*at3+0.5700_PREC
       tfm = 2.30_PREC

       !        change units from fm into bohr

       rrms = rrmsfm*fmtoau
       t = tfm*fmtoau
       a = t/(4.00_PREC*log(3.00_PREC))
       facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
       c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)

       ni=nni
       nc1=0
       xr2=r/2.00_PREC
       do mu=1,mxnmu
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc1=nc1+1
             if (iprint(231).ne.0) write(*,*) 'mu,xrr',mu,xrr
          endif
       enddo

       mu=1
       nc2=0
       do ni=nni,nni/2,-1
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc2=nc2+1
             if (iprint(231).ne.0) write(*,*) 'ni,xrr',ni,xrr
          endif
       enddo

       if (iprt.ne.0) then
          write(*,*) '    center A: '
          write(*,*) '    c = ',c,' bohr', '  a = ',a,' bohr'
          write(*,*) '    nu=-1:',nc1, ' points inside radius c'
          write(*,*) '    mu= 1:',nc2, ' points inside radius c'
          write(*,*) ' '
       endif
    elseif (z2.ne.0.0_PREC) then

       !           set atomic weight for centre B

       atw=z2atmass
       at3=atw**(1.0_PREC/3.0_PREC)
       rrmsfm = 0.8360_PREC*at3+0.5700_PREC
       tfm = 2.300_PREC

       !        change units from fm into bohr

       rrms = rrmsfm*fmtoau
       t = tfm*fmtoau
       a = t/(4.00_PREC*log (3.00_PREC))
       facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
       c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)

       ni=1
       nc1=0
       xr2=r/2.00_PREC
       do mu=1,mxnmu
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc1=nc1+1
             if (iprint(231).ne.0) write(*,*) 'mu,xrr',mu,xrr
          endif
       enddo

       mu=1
       nc2=0
       do ni=nni,nni/2,-1
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc2=nc2+1
             if (iprint(231).ne.0) write(*,*) 'ni,xrr',ni,xrr
          endif
       enddo

       if (iprint(230).ne.0) then
          write(*,*) '     Center B: '
          write(*,*) '     c = ',c,' bohr', '  a = ',a,' bohr'
          write(*,*) '     nu=1:',nc1, ' points inside radius c'
          write(*,*) '     mu= 1:',nc2, ' points inside radius c'
       endif
    endif

  end subroutine fermi
end module fermi_m
