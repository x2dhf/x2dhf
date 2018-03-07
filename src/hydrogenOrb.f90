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

module hydrogenOrb_m
  implicit none
contains


  subroutine hydrogenOrbA (orbital,n1,l1,m1,ez1)
    use params
    use discret
    use commons8
    use factor_m
    use plaguer_m
    use plegendg_m
    use normalization

    implicit none
    integer :: igp,imu,in,inioff,n1,l1,m1

    real (PREC) :: costh1,ez1,fn1,r1t,rr,shn1,z
    real (PREC), dimension(*) :: orbital

    !     normalization factor for Laguerre polynomials
    fn1=laguerre_normalization(ez1,n1,l1)
    !     normalization factor for spherical harmonics
    shn1=sphharm_normalization(l1,m1)

    do in=1,nni
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          igp=inioff+in

          !           for each grid point, e.i. for (vmu(imu),vni(ini))
          !           determine its distance |_r1| and |_r2| from the nuclei A
          !           and B and cosine of the polar angles costh1 and costh2
          !           between z axis and the vectors _r1 and _r2

          rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
          z=(r/2.0_PREC)*vxi(imu)*veta(in)
          r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))

          if (r1t.lt.precis) then
             !              jk 02/01
             costh1=-1.0_PREC
          else
             costh1=(z+r/2.0_PREC)/r1t
             ! Sanity check
             if(costh1 < -1.0_PREC) costh1 = -1.0_PREC
             if(costh1 > 1.0_PREC) costh1 = 1.0_PREC
          endif

          ! calculate radial part of the hydrogenic orbital centered
          ! on the nucleus
          orbital(igp)=fn1*plaguer(n1,l1,ez1,r1t)*shn1*plegendg(l1,m1,costh1)
       enddo
    enddo
  end subroutine hydrogenOrbA

  subroutine hydrogenOrbB (orbital,n2,l2,m2,ez2)
    use params
    use discret
    use commons8
    use factor_m
    use plaguer_m
    use plegendg_m
    use normalization

    implicit none
    integer :: igp,imu,in,inioff,n2,l2,m2

    real (PREC) :: costh2,ez2,fn2,r2t,rr,shn2,z
    real (PREC), dimension(*) :: orbital

    !     normalization factor for Laguere polynomials
    fn2=laguerre_normalization(ez2,n2,l2)
    !     normalization factor for spherical harmonics
    shn2=sphharm_normalization(l2,m2)

    do in=1,nni
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          igp=inioff+in

          !           for each grid point, e.i. for (vmu(imu),vni(ini))
          !           determine its distance |_r1| and |_r2| from the nuclei A
          !           and B and cosine of the polar angles costh1 and costh2
          !           between z axis and the vectors _r1 and _r2
          rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
          z=(r/2.0_PREC)*vxi(imu)*veta(in)
          r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

          if (r2t.lt.precis) then
             !              jk 02/01
             costh2=-one
          else
             costh2=-(z-r/2.0_PREC)/r2t
             ! Sanity check
             if(costh2 < -1.0_PREC) costh2 = -1.0_PREC
             if(costh2 > 1.0_PREC) costh2 = 1.0_PREC
          endif

          ! calculate radial part of the hydrogenic orbital centered
          ! on the nucleus
          orbital(igp)=fn2*plaguer(n2,l2,ez2,r2t)*shn2*plegendg(l2,m2,costh2)
       enddo
    enddo

  end subroutine hydrogenOrbB
end module hydrogenOrb_m
