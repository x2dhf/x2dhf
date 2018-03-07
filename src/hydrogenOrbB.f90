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
! ### hydrogenOrbB ###
!
!     This routine returns an arrys containing a hydrogenic orbital
!     centred on centre A for a given set of n,l,m and Z.

subroutine hydrogenOrbB (orbital,n2,l2,m2,ez2)
  use params
  use discret
  use commons8
  use factor_m
  use plaguer_m
  use plegendg_m

  implicit none
  integer :: igp,imu,in,inioff,n2,l2,m2

  real (PREC) :: costh2,ez2,fn2,r2t,rr,shn2,z

  real (PREC), dimension(*) :: orbital

!     normalization factor for Laguere polynomials

  fn2=(2.0_PREC*ez2/dble(n2))**(3.0_PREC/2.0_PREC+dble(l2))*sqrt(factor(n2+l2)/(2.0_PREC*dble(n2)*factor(n2-l2-1)))/factor(2*l2+1)

  shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*sqrt((2*l2+1)*factor(l2-m2)/factor(l2+m2))

  !     loop over grid points

  do in=1,nni
     do imu=1,mxnmu
        inioff=(imu-1)*nni
        igp=inioff+in

        !           for each grid point, e.i. for (vmu(imu),vni(ini))
        !           determine its distance |_r1| and |_r2| from the nuclei A
        !           and B and cosine of the polar angles costh1 and costh2
        !           between z axis and the vectors _r1 and _r2
        !
        !           rr=sqrt(xi(j,6)-eta(i,3)+1.0_PREC)
        !           rr=sqrt(vxi12(j)-veta12(i)+1.0_PREC)
        !           rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
        !           costh=eta(i,1)*xi(j,4)/rr

        rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
        z=(r/2.0_PREC)*vxi(imu)*veta(in)
        r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

        if (r2t.lt.precis) then
           !              jk 02/01
           costh2=-one
        else
           costh2=-(z-r/2.0_PREC)/r2t
        endif

        !           calculate radial part of the hydrogenic orbital centered on
        !           both the nuclei

        orbital(igp)=fn2*plaguer(n2,l2,ez2,r2t)*shn2*plegendg(l2,m2,costh2)
     enddo
  enddo

end subroutine hydrogenOrbB
