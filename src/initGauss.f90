! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### initGauss ###

!     Molecular orbitals are initialized as a linear combination of
!     Gaussian orbitals. Parameters of the orbitals and the coeeficients
!     are provided by the GAUSSIAN program.

!     Coulomb (HF) potentials are initialized as a linear combination of 
!     -ez1/r1 and -ez2/r2. For the initialization of exchange potentials 
!     see routine tfpot.

subroutine initGauss (psi,pot,excp,f2,f4,wgt2,wk0)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ic1,icen1,icen2,icen3,igauss,igp,imu,in,inioff,iorb,iorbnorm,ipnzero,ipb,iprt,ishift, &
       l1,m1,nc2,ngrid,norbt
  real (PREC) :: c1,c1sum,costh1,d1,d1r,dl1,exp_max,expt,fn1,fnorm,psi00,psi1,psi01,r1,thnorm2,thnorm3, &
       shn1,xnorm,xnorm_prv,xnorm_prev,z
  real (PREC), dimension(*) :: psi,pot,excp,f2,f4,wgt2,wk0
  real (PREC), external :: plegendg

  !     set the maximum value of an argument of the exp function 
  parameter (exp_max=700.d0)
  parameter (thnorm2=1.d-2, thnorm3=1.d-2)
  
  !     evaluate overlap matrix over gaussian basis functions
  
  igauss=0
  if (idbg(562).ne.0) igauss=1
  if (igauss.ne.0) then
!!!     call gauss_ovlap (psi,pot,excp,f2,f4,wgt2,wk0)
     stop
  endif

  !     Initialization of molecular orbitals 

  if (ini.eq.4.) then
     norbt=1
  else
     norbt=norb
  endif
  
  write(*,1114)
  icen1=1
  icen2=2
  icen3=3
  iorb=0
  iorbnorm=1

900 continue
  
  !     Gaussian normal output
  !     In the case of homonuclear molecules Gaussian program reverses 
  !     the order of the centres (with respect to the order defined  
  !     by the input data (sigma g/u and nonsigma orbitals are effected). 
  !     That is why the other ordering is tried when the norm of a numerical 
  !     orbital is off from unity by thnorm2.
  
  if     ((ini.eq.2.or.ini.eq.3).and.idbg(565).eq.0) then
     icen1=1
     icen2=2
     icen3=3
     xnorm_prv=0.0_PREC
  elseif ((ini.eq.2.or.ini.eq.3).and.idbg(565).eq.1) then
     icen1=3
     icen2=2
     icen3=1
  endif
  iorb=iorb+1

  if (iorb.gt.norbt) goto 910  

  ! FIXME
  if (iorbnorm.eq.1) then
     if(ini.eq.3.and.mm(iorb).ne.0) then 
        iprt=0
     else
        iprt=1
     endif
     if (ini.eq.2) then
        iprt=0
     endif
  endif
  
  
  psi00=0.0_PREC
  psi01=0.0_PREC
  
  if (idbg(566).ne.0) then
     i1b(iorb)=1
     ngrid= i1si(iorb)
     do i=1,ngrid
        psi(i)=0.0_PREC
     enddo
  endif
  
  ishift=i1b(iorb)-1
  ngrid= i1si(iorb)
  igp=ishift
  do i=1,ngrid
     psi(ishift+i)=0.0_PREC
  enddo
  
  !     loop over basis set functions
  
  ipnzero=0
  nc2=0
  c1sum=0.0_PREC
  do ipb=1,npbasis
     l1  =lprim(ipb)
     m1  =abs(mprim(ipb))
!     print *,iorb,ipb,l1,m1,mm(iorb)
     if (m1.ne.mm(iorb)) goto 999
     !        Gaussian long (customized) output.
     !        In most cases PX components are used to construct numerical
     !        molecular orbital. Sometimes PY components are used instead to
     !        get the norm of the orbital correct. If the norm is off by more
     !        than thnorm3 than idbg(561) is set to 1 and the orbital is 
     !        initialized afresh. 

     if (ini.eq.3.and.idbg(561).eq.0) then
        if (mprim(ipb).lt.0) goto 999 
     elseif (ini.eq.3.and.idbg(561).eq.1) then
        if (mprim(ipb).gt.0) goto 999 
     endif
     ic1 =icgau(ipb)
     d1  =primexp(ipb)
     
     !        sometimes the gaussian long output contains unprintable exponents
     !        which must be singled out
     if (ini.eq.3) then
        if (d1.gt.99999999.0_PREC) goto 999
     endif
     
     fn1 =fngau2(ipb)
     shn1=shngau(ipb)
     fnorm=fn1*shn1
     c1 =primcoef(iorb,ipb)

     if (abs(c1).gt.0.0_PREC) then
        ipnzero=ipnzero+1
 !                  write(*,1111) iorb,ipb,ic1,l1,m1,d1,c1,fn1,shn1
 !                  write(*,1112) iorb,ipb,ic1,l1,mprim(ipb),c1,d1
1111    format(5i4,f20.9,3e15.5)
1112    format(5i4,2f20.9)
        
        !           ic1=1 basis function at centre Z1
        !           ic1=2 basis function at bond centre
        !           ic1=3 basis function at centre Z2
        
        !           loop over grid points
        
        do imu=1,mxnmu
           inioff=(imu-1)*nni
           do in=1,nni
              igp=ishift+inioff+in
              
              !                 for each grid point, i.e. for (vmu(imu),vni(ini)) determine 
              !                 its distance |_r1| from the nuclei Z_1 and Z_2 and cosine 
              !                 of the polar angles costh1 and costh2 betp ween z axis 
              !                 and the vectors _r1 and _r2
              
              !                 rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)		
              
              z=(r/2.0_PREC)*vxi(imu)*veta(in)
              if (ic1.eq.icen1) then
                 !                    gaussians centred on Z1
                 r1=(r/2.0_PREC)*(vxi(imu)+veta(in))
                 if (r1.lt.precis) then
                    costh1=0.0_PREC
                 else   
                    costh1=(z+r/2.0_PREC)/r1
                 endif
              elseif (ic1.eq.icen3) then
                 !                    gaussians centred on Z2
                 r1=(r/2.0_PREC)*(vxi(imu)-veta(in))
                 if (r1.lt.precis) then
                    costh1=0.0_PREC
                 else   
                    costh1=(z-r/2.0_PREC)/r1	
                 endif
              elseif (ic1.eq.icen2) then
                 !                    gaussians centred on the bond centre
                 r1=sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)*r/2.0_PREC
                 if (r1.lt.precis) then
                    costh1=0.0_PREC
                 else   
                    costh1= (z)/r1	
                 endif
              endif
              
              !                 calculate the gaussian function centered on one of the nuclei
              
              d1r=d1*r1*r1
              if (d1r.gt.exp_max) then
                 expt=0.0_PREC
              else
                 expt=1.0_PREC/exp(d1r)
              endif
              
              if (r1.lt.precis.and.l1.eq.0) then
                 psi1=fnorm*expt*plegendg(l1,m1,costh1)
              elseif (r1.lt.precis.and.l1.gt.0) then
                 psi1=0.0_PREC
              else
                 dl1=l1
                 psi1=fnorm*r1**dl1*expt*plegendg(l1,m1,costh1)
              endif
99            continue
              psi(igp)=psi(igp)+c1*psi1
           enddo
        enddo
        !           c1.ne.0
     endif
     !        
999  continue
  enddo

  call norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
  write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm

!  if (iprt.eq.1) write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm
  
  if (iorbnorm.eq.1) then
     !         if(ini.eq.2) then
     if(ini.eq.2.and.abs(xnorm-1.0_PREC).gt.thnorm2) then
        
        !   Try reversed ordering of atomic centres. This might lead to a better orbital
        idbg(561)=1
        iorb=iorb-1
        iorbnorm=2
        xnorm_prev=xnorm
     endif
     !         if(ini.eq.3.and.mm(iorb).ne.0) then
     if(ini.eq.3.and.mm(iorb).ne.0.and.abs(xnorm-1.0_PREC).gt.thnorm3) then
        
        idbg(561)=1
        iorb=iorb-1
        iorbnorm=2
        xnorm_prev=xnorm
     endif
  elseif (iorbnorm.eq.2) then
     if (abs(xnorm-1.0_PREC).le.abs(xnorm_prev-1.0_PREC)) then
        iorbnorm=3
        idbg(561)=0
        iprt=1
     else
        iorbnorm=3
        xnorm=xnorm_prev
        iorb=iorb-1
        idbg(561)=0
        iprt=1
     endif
  endif
  goto 900
910 continue
  write(*,*)
  
  if (idbg(560).ne.0) stop 'inigauss'
  
  !     initialize Coulomb and exchange potentials
  
  call initPot(psi,pot,excp,f2,f4,wk0)	
  
1114 format(/1x,'    orbital        norm      ')
1115 format(1x,i3,1x,a8,a1,e20.12)
  
end subroutine initGauss
      
