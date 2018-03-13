! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### initHF ###

!     This routine initializes molecular orbitals as linear combinations
!     of Hartree-Fock functions (taken from the qrhf program) on centres
!     A and B.  In the case of HF or HFS calculations Coulomb (exchange)
!     potentials are approximated as a linear combination of
!     Thomas-Fermi (1/r) potentials of the two centres. If method OED is
!     chosen the potential functions are set to zero.

module initHF_m
  implicit none
contains
  subroutine initHF (psi,pot,excp,f2,f4,wgt2,wk0)
    use params
    use discret
    use commons8

    use factor_m
    use flp_m
    use norm94_m
    use initPot_m
    use plegendg_m

    implicit none


    integer ::  igp,ihf1,ihf2,ilabel,imu,in,inioff,iord,ishift,mxmax1,mxmax2, &
         nhforb,ngridpts,nwf1,nwf2,ouf2dhf1,ouf2dhf2
    integer :: i,j,iorb,l1,m1,n1,l2,m2,n2
    real (PREC) :: costh1,costh2,psi1,psi2,psi1prev,psi2prev,shn1,shn2,r1t,r2t,rr,xnorm, &
         z,zhf1,zhf2

    parameter (nhforb=20,ngridpts=3000,iord=3,ouf2dhf1=8,ouf2dhf2=9)

    integer, dimension(nhforb) :: nhf1,lhf1,nhf2,lhf2

    real (PREC), dimension(*) :: psi,pot,excp,f2,f4,wgt2,wk0
    real (PREC), dimension(ngridpts) :: rhf1,rhf2,phf1t,phf2t
    real (PREC), dimension(nhforb) :: ehf1,qc1,ehf2,qc2
    real (PREC), dimension(nhforb,ngridpts) :: phf1,phf2

    character*8 :: atom1,term1,atom2,term2



    !     Initialization of molecular orbitals

    print *,'... initializing molecular orbitals from Hartree-Fock functions ...'

    ilabel=0

    !     read HF functions for centre A
    if (z1.ne.zero) then
       open(ouf2dhf1,file='1dhf_centreA.orb',form='formatted',status='old')

       read(ouf2dhf1,*) atom1,term1,zhf1,nwf1,mxmax1

       if (mxmax1.gt.ngridpts) then
          write(*,*) "Error: too many grid points in the 1dHF data", &
               &           " for centre A -- increase the value of ngridpts parameter"
          stop 'initHF'
       endif

       read(ouf2dhf1,'(2x,20i5)') (nhf1(i),i=1,nwf1)
       read(ouf2dhf1,'(2x,20i5)') (lhf1(i),i=1,nwf1)
       read(ouf2dhf1,'(3x,20f5.0)') (qc1(i),i=1,nwf1)
       read(ouf2dhf1,'(24x,20d24.16)') (ehf1(i),i=1,nwf1)
       do j=1,mxmax1
          read(ouf2dhf1,'(20d24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
       enddo

       if (iprint(220).ne.0) then
          write(*,'(" Orbitals on centre Z1 (Z=",f4.0,"):")') zhf1
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf1
             write(*,'(4x,3i5,e16.6)') i,nhf1(i),lhf1(i), ehf1(i)
          enddo

       endif

       if (iprint(222).ne.0) then
          do j=1,mxmax1
             write(*,'(20e24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
          enddo
       endif

    endif
    close(ouf2dhf1)
    
    !     read HF functions for centre B

    if (z2.ne.zero) then
       open(ouf2dhf2,file='1dhf_centreB.orb',form='formatted',status='old')

       read(ouf2dhf2,*) atom2,term2,zhf2,nwf2,mxmax2

       if (mxmax2.gt.ngridpts) then
          write(*,*) "Error: too many grid points in the 1dHF data", &
               " for centre B -- increase the value of ngridpts parameter"
          stop 'initHF'
       endif

       read(ouf2dhf2,'(2x,20i5)') (nhf2(i),i=1,nwf2)
       read(ouf2dhf2,'(2x,20i5)') (lhf2(i),i=1,nwf2)
       read(ouf2dhf2,'(3x,20f5.0)') (qc2(i),i=1,nwf2)
       read(ouf2dhf2,'(24x,20d24.16)') (ehf2(i),i=1,nwf2)
       do j=1,mxmax2
          read(ouf2dhf2,'(20d24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
       enddo

       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z2 (Z=",f4.0,"):")') zhf2
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf2
             write(*,'(4x,3i5,e16.6)') i,nhf2(i),lhf2(i), ehf2(i)
          enddo

       endif

       if (iprint(222).ne.0) then
          do j=1,mxmax2
             write(*,'(20e24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
          enddo
       endif
    endif
    close(ouf2dhf2)
    
    !     loop over orbitals

    do iorb=1,norb

       ishift=i1b(iorb)-1
       n1=mgx(1,iorb)
       l1=mgx(2,iorb)
       m1=mgx(3,iorb)

       !        normalization factor for spherical harmonics

       shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))
       if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))

       if (co1(iorb).ne.0.0_PREC) then
          ihf1=0
          do i=1,nwf1
             if (nhf1(i).eq.n1.and.lhf1(i).eq.l1) ihf1=i
          enddo
          if (ihf1.eq.0) then
             print *,"initHF: no proper atomic orbital for centre Z1 found"
             stop 'iniHF'
          endif

          do j=1,mxmax1
             phf1t(j)=phf1(ihf1,j)
          enddo
       endif

       n2=mgx(4,iorb)
       l2=mgx(5,iorb)
       m2=mgx(6,iorb)

       !        normalization factor for spherical harmonics

       shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*sqrt((2*l2+1)*factor(l2-m2)/factor(l2+m2))

       if (co2(iorb).ne.0.0_PREC) then
          ihf2=0
          do i=1,nwf2
             if (nhf2(i).eq.n2.and.lhf2(i).eq.l2) ihf2=i
          enddo
          if (ihf2.eq.0) then
             print *,"inihf: no proper atomic orbital for centre Z2 found"
             stop 'iniHF'
          endif

          do j=1,mxmax2
             phf2t(j)=phf2(ihf2,j)
          enddo
       endif

       if (iprint(221).ne.0) then
          print *,'inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
          print *,'inihf: n2,l2,m2,shn2,ihf2',n2,l2,m2,shn2,ihf2
       endif

       !        loop over grid points

       psi1prev=0.0_PREC
       psi2prev=0.0_PREC
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=ishift+inioff+in

             !              for each grid point, e.i. for (vmu(imu),vni(ini))
             !              determine its distance |_r1| and |_r2| from the nuclei A
             !              and B and cosine of the polar angles costh1 and costh2
             !              between z axis and the vectors _r1 and _r2

             psi1=0.0_PREC
             psi2=0.0_PREC

             rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
             z=(r/2.0_PREC)*vxi(imu)*veta(in)
             r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))
             r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

             !              calculate radial part of the hydrogenic orbital centered
             !              on both the nuclei

             if (r1t.lt.precis) then
                costh1=0.0_PREC
             else
                costh1=(z+r/2.0_PREC)/r1t
             endif
             !
             if (r2t.lt.precis) then
                costh2=0.0_PREC
             else
                costh2=-(z-r/2.0_PREC)/r2t
             endif

             !              calculate radial part of the HF orbital taken from (R)HF
             !              program as P=R*r

             if     (co1(iorb).ne.0.0_PREC) then
                if (r1t.ge.rhf1(mxmax1)) then
                   psi1=0.0_PREC
                elseif (r1t.ge.precis) then
                   psi1=flp(iord,mxmax1,rhf1,phf1t,r1t)/(r1t)*shn1*plegendg(l1,m1,costh1)
                   psi1prev=psi1
                else
                   psi1=psi1prev
                endif
             endif

             if (co2(iorb).ne.0.0_PREC) then
                if (r2t.ge.rhf2(mxmax2)) then
                   psi2=0.0_PREC
                elseif (r2t.ge.precis) then
                   psi2=flp(iord,mxmax2,rhf2,phf2t,r2t)/(r2t)*shn2*plegendg(l2,m2,costh2)
                   psi2prev=psi2
                else
                   psi2=psi2prev
                endif
             endif
             psi(igp)=co1(iorb)*psi1+co2(iorb)*psi2
          enddo
       enddo
       if (iprint(220).ne.0) then
          if (ilabel.eq.0) then
             ilabel=1
             write (*,*)
             write (*,*) 'Normalization of LCAOs'
          endif
          call norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
          write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm
1115      format(i4,1x,a8,a1,3x,e22.16,2e16.2)
       endif
    enddo

    !  initialize Coulomb and exchange potentials

    call initPot(psi,pot,excp,f2,f4,wk0)
  end subroutine initHF
end module initHF_m
