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
! c ### hydrogenOrbA ###
! c
! c     This routine returns an arrys containing a hydrogenic orbital
! c     centred on centre A for a given set of n,l,m and Z.

subroutine hydrogenOrbA (orbital,n1,l1,m1,ez1)
  use params
  use discret
  use commons8
  use factor_m
  use plaguer_m
  use plegendg_m

  implicit none
  integer :: igp,imu,in,inioff,n1,l1,m1

  real (PREC) :: costh1,ez1,fn1,r1t,rr,shn1,z
  real (PREC), dimension(*) :: orbital

  !     normalization factor for Laguere polynomials

  fn1=(2.0_PREC*ez1/dble(n1))**(3.0_PREC/2.0_PREC+dble(l1))*sqrt(factor(n1+l1)/(2.0_PREC*dble(n1)*factor(n1-l1-1)))/factor(2*l1+1)

  !     normalization factor for spherical harmonics

  shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))
  if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))

  !          multiplication factor fot the associate Legendre polynomials
  !     alp1=(-1.0_PREC/2.0_PREC)**m1*factor(l1+m1)/factor(l1-m1)/factor(m1)

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
        r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))

        if (r1t.lt.precis) then
           !              jk 02/01
           costh1=-1.0_PREC
        else
           costh1=(z+r/2.0_PREC)/r1t
        endif

        !           calculate radial part of the hydrogenic orbital centered on
        !           both the nuclei
        orbital(igp)=fn1*plaguer(n1,l1,ez1,r1t)*shn1*plegendg(l1,m1,costh1)
     enddo
  enddo

end subroutine hydrogenOrbA
