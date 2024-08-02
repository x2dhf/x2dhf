! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module initOrbitalsPotentials
  implicit none
contains
  ! ### initOrbPot ###
  !
  !     Initializes orbitals and potentials.
  !
  subroutine initOrbPot 
    use params
    use discrete
    use sharedMemory    
    use commons
    use prepGauss
    use diskInterface
    use coulombExchangePotentials
    use storedorbitals
    
    implicit none

    integer (KIND=IPREC) :: iorb
    character*20 caseinp,psiinp,coulinp,exchinp,caseout,psiout,coulout,exchout

    data caseinp,psiinp,coulinp,exchinp /'2dhf_input.dat', '2dhf_input.orb','2dhf_input.coul','2dhf_input.exch'/
    data caseout,psiout,coulout,exchout /'2dhf_output.dat', '2dhf_output.orb','2dhf_output.coul','2dhf_output.exch'/

    ! The program works on a set of standard input and output files containing orbitals,
    ! Coulomb potentials (with extentions orb and coul, respectively), exchange
    ! potentials (a file with extension exch or file fort.31, fort.32, ...) and a text
    ! file with data defining the case (*.dat).

    ! open separate files for reading and writing data defining a case
    open(iinp14,file=caseinp, status='unknown',form='formatted')
    open(iout24,file=caseout, status='unknown',form='formatted')

    ! open separate files for reading and writing orbitals and potentials

    ! default iform: in-one out-one
    if (inUnformatted) then 
       open(11,file=psiinp, status='unknown',form='unformatted')
       open(12,file=coulinp,status='unknown',form='unformatted')
       open(13,file=exchinp,status='unknown',form='unformatted')
       rewind(13)
    else
       open(11,file=psiinp, status='unknown',form='formatted')
       open(12,file=coulinp,status='unknown',form='formatted')
       open(13,file=exchinp,status='unknown',form='formatted')
       rewind(13)
    endif

    if (outUnformatted) then    
       open(21,file=psiout, status='unknown',form='unformatted')
       open(22,file=coulout,status='unknown',form='unformatted')
       open(23,file=exchout,status='unknown',form='unformatted')
       rewind(23)
    else
       open(21,file=psiout, status='unknown',form='formatted')
       open(22,file=coulout,status='unknown',form='formatted')
       open(23,file=exchout,status='unknown',form='formatted')
       rewind(23)
    endif
 
    rewind(11)
    rewind(12)
    rewind(13)
    rewind(21)
    rewind(22)
    rewind(23)
    
    initAddData=.false.
    if (linitFuncsHydrogen) then
       call initHyd
       initAddData=.true.
       do iorb=1,norb
          inhyd(iorb)=inhydlcao(iorb)
       enddo
    elseif (linitFuncsQRHF) then
       call initQRHF
       initAddData=.true.       
    elseif (linitFuncsLDA.or.linitFuncsHF) then
       if (abs(z1-z2).lt.homolevl.and. .not.lbreakCi) then 
          call initStoredOrbs4homo
       else
          call initStoredOrbs 
       endif
       initAddData=.true.       
       
    elseif (linitFuncsGauss) then
       ! initial values of orbitals are provided by the GAUSSIAN program
       write(*,*) 'Initializing orbitals using GAUSSIAN output'
       call prepareGaussian
       call initGauss
       initAddData=.true.       

    elseif (linitFuncsGaussC) then       
       ! gauss-c
       stop 'Support for customized GAUSSIAN format has been deprecated. Use the standard format.'

    elseif (initFuncsOED) then
       ! method oed
       !idump=0
       call readFromDisk 

       ! To generate a higher lying one-electron diatomic state all the lower
       ! ones have to be retrieved from disk file and kept for the sake of
       ! the orthogonalization.

       ! The current highest orbital is being initialized as a linear
       ! combination of hydrogenic atomic orbitals unless the previous
       ! case is being continued

       if (ini4.ne.0) then
          call initHyd
          initAddData=.true.       
       endif

    elseif (linitFuncsOld) then
       ! When the incomplete orbital disk file is encountered and the method is OED the
       ! values of ini4 is set to norb-norb_p
       call readFromDisk

       if (ini4.ne.0) then
          call initHyd
          initAddData=.true.       
       endif

    elseif (linitFuncsNoexch) then
       call readFromDisk 
       call initCoulomb
       call initExchange
       initAddData=.true.
       
    elseif (linitFuncsNodat) then
       ! old format of input data is requested (skipping reading 2dhf_input.dat file)
       ldatPresent=.false.
       call readFromDisk 
    endif

  end subroutine initOrbPot


  ! ### initGauss ###
  !
  !     Molecular orbitals are initialized as a linear combination of
  !     Gaussian orbitals. Parameters of the orbitals and the coeeficients
  !     are provided by the GAUSSIAN program.
  !
  !     Coulomb (HF) potentials are initialized as a linear combination of
  !     -ez1/r1 and -ez2/r2. For the initialization of exchange potentials
  !     see routine tfpot.
  !
  subroutine initGauss 
    use params
    use discrete
    use commons
    use sharedMemory
    use normOrtho
    use coulombExchangePotentials
    
    implicit none
    real (PREC), parameter :: normThreshold=1.0e-10_PREC
    integer (KIND=IPREC) :: i,igauss,igp,igrid,imu,in,inioff,iorb,ipb,ishift,m1, &
         ngrid,norbt

    real (PREC) :: c1,xnorm
    real (PREC), dimension(:,:), allocatable :: bf
    integer (KIND=IPREC) :: it, jt

    ! Values of contracted basis functions
    real (PREC), dimension(:,:), allocatable :: cbf
    ! m values of contracted basis functions
    integer (KIND=IPREC),dimension(:), allocatable :: cm
    ! Contracted basis function overlap
    real (PREC), dimension(:,:), allocatable :: cS
    real (PREC), dimension(:), pointer :: excp,excp1,f2,f4,psi,wgt2,wk0

    psi=>orbptr
    excp=>scratchptr
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr
    
    !     Initialization of molecular orbitals
    if (initFuncsOED) then
       ! If OED method is chosen only the highest orbital is
       ! initialized
       norbt=1
    else
       norbt=norb
    endif

    ! Evaluate basis functions
    allocate(bf(nni*mxnmu,npbasis))
    call evalGauss(bf)

    ! Evaluate overlap matrix over gaussian basis functions
    igauss=0
    if (idebug(562).ne.0) igauss=1
    if (igauss.ne.0) then
       write (*,*) 'Test finite basis primitive function overlaps on grid'
       do it=1,npbasis
          do jt=1,it
             ! Skip different m values
             if(mprim(it) .ne. mprim(jt)) cycle
             xnorm=0.0_PREC
             do i=1,nni*mxnmu
                xnorm=xnorm+bf(i,it)*bf(i,jt)*wgt2(i)*f4(i)
             end do
             !           if(abs(xnorm) .ge. 1e-6) then
             write (*,'(A,I3,A,I3,A,ES14.7)') 'S(',it,',',jt,') = ',xnorm
             !           end if
          end do
       end do

       write (*,*) 'Test finite basis contracted function overlaps on grid'
       ! Evaluate contracted basis function values on grid
       allocate(cbf(nni*mxnmu,ixref(npbasis)))
       allocate(cm(ixref(npbasis)))
       allocate(cS(ixref(npbasis),ixref(npbasis)))
       do it=1,ixref(npbasis)
          do i=1,nni*mxnmu
             cbf(i,it)=0.0_PREC
          end do
          do igp=1,npbasis
             if(ixref(igp).eq.it) then
                ! m value
                cm(it)=mprim(igp)
                do i=1,nni*mxnmu
                   cbf(i,it)=cbf(i,it) + coeff(igp)*bf(i,igp)
                end do
             end if
          end do
       end do
       do it=1,ixref(npbasis)
          do jt=1,it
             ! Skip different m blocks since those are ignored in code
             xnorm=0.0_PREC
             if(cm(it) .eq. cm(jt)) then
                ! Calculate norm
                do i=1,nni*mxnmu
                   xnorm=xnorm+cbf(i,it)*cbf(i,jt)*wgt2(i)*f4(i)
                end do
             end if
             if(abs(xnorm) <= normThreshold) xnorm=0.0_PREC
             cS(it,jt)=xnorm
             cS(jt,it)=xnorm
             write (*,'(1X,E12.6)',advance='no') cS(it,jt)
          end do
          write (*,*) ''
       end do
       deallocate(cbf)
       deallocate(cm)
       deallocate(cS)
    end if

    write(*,1114)
    do iorb=1,norbt
       ishift=i1b(iorb)-1
       ngrid= i1si(iorb)
       igp=ishift
       do i=1,ngrid
          psi(ishift+i)=0.0_PREC
       enddo

       do ipb=1,npbasis
          m1 = abs(mprim(ipb))
          ! Skip functions that have different m value
          if (m1.ne.mm(iorb)) cycle

          c1 =primcoef(iorb,ipb)

          if (abs(c1).gt.0.0_PREC) then
             ! Loop over grid
             do imu=1,mxnmu
                inioff=(imu-1)*nni
                do in=1,nni
                   igrid=inioff+in
                   igp=ishift+igrid
                   psi(igp)=psi(igp)+c1*bf(igrid,ipb)
                end do
             end do
          end if
       end do
       !     write (*,*) 'Orbital ',iorb,' overlap'
       call norm94 (iorb,xnorm)
       write (*,1115) iorn(iorb),bond(iorb),gusym(iorb),xnorm
    end do

    deallocate(bf)

    !  if (iprt.eq.1) write (*,1115) iorn(iorb),bond(iorb),gusym(iorb),xnorm

    write(*,*)

    if (idebug(560).ne.0) stop 'inigauss'

    ! initialize Coulomb and exchange potentials
    call initCoulomb
    call initExchange

1114 format(/1x,'    orbital        norm      ')
1115 format(1x,i3,1x,a8,a1,e20.12)

  end subroutine initGauss

  subroutine evalGauss(bf)
    use params
    use discrete
    use commons
    use utils
    implicit none

    integer (KIND=IPREC) :: ic1, igp, imu, in, inioff, ipb, l1, m1
    real (PREC) :: costh1, d1, expt, fnorm, r1, z, psi1
    real (PREC), dimension(:,:) :: bf

    ! Sanity check
    if(size(bf,1) .ne. nni*mxnmu .or. size(bf,2) .ne. npbasis) then
       write (*,*) 'bf matrix needs to be initialized, got size ',size(bf,1),' x ',size(bf,2)
    end if

    !     loop over basis set functions
    do ipb=1,npbasis
       ! Index of the center of the primitive
       !   ic1=1 basis function at centre Z1
       !   ic1=2 basis function at bond centre
       !   ic1=3 basis function at centre Z2
       ic1 =icgau(ipb)

       ! Primitive exponent
       d1  =primexp(ipb)
       ! (l,m) values
       l1  =lprim(ipb)
       m1  =abs(mprim(ipb))

       ! Normalization: primitive times spherical harmonic
       fnorm=fngau2(ipb)*shngau(ipb)

       ! Loop over grid
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=inioff+in

             ! for each grid point, i.e. for (vmu(imu),vni(ini))
             ! determine its distance |r_i| from the nucleus Z_i and
             ! cosine of the polar angle costh_i between the z axis
             ! and the vector r_i

             z=(r/2.0_PREC)*vxi(imu)*veta(in)

             if (ic1.eq.1) then
                ! Center on Z1
                r1=(r/2.0_PREC)*(vxi(imu)+veta(in))
                z=z+r/2.0_PREC

             elseif (ic1.eq.3) then
                ! Center on Z2
                r1=(r/2.0_PREC)*(vxi(imu)-veta(in))
                z=z-r/2.0_PREC

             elseif (ic1.eq.2) then
                ! Bond center
                r1=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)

             else
                write (*,*) 'invalid center ',ic1
                stop
             endif

             if (r1.lt.precis) then
                costh1=0.0_PREC
             else
                costh1=z/r1
                ! Sanity check
                if(costh1 < -1.0_PREC) costh1 = -1.0_PREC
                if(costh1 > 1.0_PREC) costh1 = 1.0_PREC
             endif

             ! Calculate exponential
             expt=exp(-d1*r1*r1)

             ! Value of basis function is
             if (r1.lt.precis) then
                if(l1.gt.0) then
                   psi1=0.0_PREC
                else
                   psi1=fnorm*expt*plegendg(l1,m1,costh1)
                end if
             else
                psi1=fnorm*r1**l1*expt*plegendg(l1,m1,costh1)
             endif

             ! Store value
             bf(igp,ipb)=psi1
          enddo
       enddo
    enddo
  end subroutine evalGauss

  ! ### initQRHF ###
  !
  !     This routine initializes molecular orbitals as linear combinations
  !     of Hartree-Fock functions (taken from the qrhf program) on centres
  !     A and B.  In the case of HF or HFS calculations Coulomb (exchange)
  !     potentials are approximated as a linear combination of
  !     Thomas-Fermi (1/r) potentials of the two centres. If method OED is
  !     chosen the potential functions are set to zero.
  !
  subroutine initQRHF 
    use params
    use discrete
    use sharedMemory    
    use commons
    use interpolate
    use normOrtho
    use utils
    use coulombExchangePotentials
    
    implicit none

    integer (KIND=IPREC) ::  igp,ihf1,ihf2,ilabel,imu,in,inioff,iorderl,ishift,mxmax1,mxmax2, &
         nhforb,ngridpts,nwf1,nwf2,ouf2dhf1,ouf2dhf2
    integer (KIND=IPREC) :: i,j,iorb,l1,m1,n1,l2,m2,n2
    real (PREC) :: costh1,costh2,psi1,psi2,psi1prev,psi2prev,shn1,shn2,r1t,r2t,rr,xnorm, &
         z,zhf1,zhf2

    parameter (nhforb=20,ngridpts=3000,iorderl=3,ouf2dhf1=8,ouf2dhf2=9)

    integer (KIND=IPREC),dimension(nhforb) :: nhf1,lhf1,nhf2,lhf2

    real (PREC), dimension(ngridpts) :: rhf1,rhf2,phf1t,phf2t
    real (PREC), dimension(nhforb) :: ehf1,qc1,ehf2,qc2
    !real (PREC), dimension(nhforb,ngridpts) :: phf1,phf2
    real (PREC), dimension(:,:), allocatable :: phf1,phf2

    character*8 :: atom1,term1,atom2,term2

    real (PREC), dimension(:), pointer :: excp,f4,psi,wgt2,wk0

    allocate (phf1(nhforb,ngridpts))
    allocate (phf2(nhforb,ngridpts))
    
    excp=>exchptr
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr
    
    ! Initialization of molecular orbitals

    print *,'... initializing molecular orbitals from QRHF functions ...'

    ilabel=0

    ! read HF functions for centre A
    if ( abs(z1)>epsilon(zero)) then
       open(ouf2dhf1,file='1dhf_centreA.orb',form='formatted',status='old')

       read(ouf2dhf1,*) atom1,term1,zhf1,nwf1,mxmax1

       if (mxmax1.gt.ngridpts) then
          write(*,*) "Error: too many grid points in the 1dHF data", &
               &           " for centre A -- increase the value of ngridpts parameter"
          stop 'initQRHF'
       endif

       read(ouf2dhf1,'(2x,20i5)') (nhf1(i),i=1,nwf1)
       read(ouf2dhf1,'(2x,20i5)') (lhf1(i),i=1,nwf1)
       read(ouf2dhf1,'(3x,20f5.0)') (qc1(i),i=1,nwf1)
       read(ouf2dhf1,'(24x,20d24.16)') (ehf1(i),i=1,nwf1)
       do j=1,mxmax1
          read(ouf2dhf1,'(20d24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
       enddo
       
#ifdef PRINT
! print=220: initQRHF: Orbitals on centre Z1: i,nhf1(i),lhf1(i), ehf1(i) 
       if (iprint(220).ne.0) then
          write(*,'(" Orbitals on centre Z1 (Z=",f4.0,"):")') zhf1
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf1
             write(*,'(4x,3i5,e16.6)') i,nhf1(i),lhf1(i), ehf1(i)
          enddo
       endif
#endif
 
#ifdef PRINT      
! print=222: initQRHF: Orbitals on centre Z1: rhf1(j),(phf1(i,j),i=1,nwf1) 
       if (iprint(222).ne.0) then
          do j=1,mxmax1
             write(*,'(20e24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
          enddo
       endif
#endif
       close(ouf2dhf1)
    endif
    
    ! read HF functions for centre B
    if ( abs(z2)>epsilon(zero)) then
       open(ouf2dhf2,file='1dhf_centreB.orb',form='formatted',status='old')

       read(ouf2dhf2,*) atom2,term2,zhf2,nwf2,mxmax2

       if (mxmax2.gt.ngridpts) then
          write(*,*) "Error: too many grid points in the 1dHF data", &
               " for centre B -- increase the value of ngridpts parameter"
          stop 'initQRHF'
       endif

       read(ouf2dhf2,'(2x,20i5)') (nhf2(i),i=1,nwf2)
       read(ouf2dhf2,'(2x,20i5)') (lhf2(i),i=1,nwf2)
       read(ouf2dhf2,'(3x,20f5.0)') (qc2(i),i=1,nwf2)
       read(ouf2dhf2,'(24x,20d24.16)') (ehf2(i),i=1,nwf2)
       do j=1,mxmax2
          read(ouf2dhf2,'(20d24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
       enddo

#ifdef PRINT
! print=220: initQRHF: Orbitals on centre Z2: i,nhf2(i),lhf2(i), ehf2(i) 
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z2 (Z=",f4.0,"):")') zhf2
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf2
             write(*,'(4x,3i5,e16.6)') i,nhf2(i),lhf2(i), ehf2(i)
          enddo
       endif
#endif

#ifdef PRINT       
! print=222: initQRHF: Orbitals on centre Z2: rhf2(j),(phf2(i,j),i=1,nwf2) o       
       if (iprint(222).ne.0) then
          do j=1,mxmax2
             write(*,'(20e24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
          enddo
       endif
#endif
       close(ouf2dhf2)
    endif
    
    ! loop over orbitals

    do iorb=1,norb

       ishift=i1b(iorb)-1
       n1=mgx(1,iorb)
       l1=mgx(2,iorb)
       m1=mgx(3,iorb)

       !print *, "iorb,n1,l1:",iorb,n1,l1
       
       ! normalization factor for spherical harmonics

       shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))
       if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))

       !if (co1(iorb).ne.0.0_PREC) then
       if ( abs(co1(iorb))>epsilon(zero)) then
       ihf1=0
          do i=1,nwf1
             if (nhf1(i).eq.n1.and.lhf1(i).eq.l1) ihf1=i
             !print *, "i,nhf1(i),lhf1(i)",i,nhf1(i),lhf1(i)
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

       ! normalization factor for spherical harmonics

       shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l2+1)*factor(l2-m2)/factor(l2+m2))

       !if (co2(iorb).ne.0.0_PREC) then
       if ( abs(co2(iorb))>epsilon(zero)) then
             
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

#ifdef PRINT
! print=221: initQRHF: n1,l1,m1,shn1,ihf1/n2,l2,m2,shn2,ihf2
       if (iprint(221).ne.0) then
          print *,'inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
          print *,'inihf: n2,l2,m2,shn2,ihf2',n2,l2,m2,shn2,ihf2
       endif
#endif
       ! loop over grid points

       psi1prev=0.0_PREC
       psi2prev=0.0_PREC
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=ishift+inioff+in

             ! for each grid point, e.i. for (vmu(imu),vni(ini))
             ! determine its distance |_r1| and |_r2| from the nuclei A
             ! and B and cosine of the polar angles costh1 and costh2
             ! between z axis and the vectors _r1 and _r2

             psi1=0.0_PREC
             psi2=0.0_PREC

             rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
             z=(r/2.0_PREC)*vxi(imu)*veta(in)
             r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))
             r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

             ! calculate radial part of the hydrogenic orbital centered
             ! on both the nuclei

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

             ! calculate radial part of the HF orbital taken from (R)HF
             ! program as P=R*r

             !if     (co1(iorb).ne.0.0_PREC) then
                if ( abs(co1(iorb))>epsilon(zero)) then
                if (r1t.ge.rhf1(mxmax1)) then
                   psi1=0.0_PREC
                elseif (r1t.ge.precis) then
                   psi1=flp(iorderl,mxmax1,rhf1,phf1t,r1t)/(r1t)*shn1*plegendg(l1,m1,costh1)
                   psi1prev=psi1
                else
                   psi1=psi1prev
                endif
             endif

             !if (co2(iorb).ne.zero) then
             if ( abs(co2(iorb))>epsilon(zero)) then
                if (r2t.ge.rhf2(mxmax2)) then
                   psi2=0.0_PREC
                elseif (r2t.ge.precis) then
                   psi2=flp(iorderl,mxmax2,rhf2,phf2t,r2t)/(r2t)*shn2*plegendg(l2,m2,costh2)
                   psi2prev=psi2
                else
                   psi2=psi2prev
                endif
             endif
             psi(igp)=co1(iorb)*psi1+co2(iorb)*psi2
          enddo
       enddo

#ifdef PRINT
! print=223: initQRHF: print out norms of LCAOs
       if (iprint(223).ne.0) then
          if (ilabel.eq.0) then
             ilabel=1
             write (*,*)
             write (*,*) 'Normalization of LCAOs'
          endif
          call norm94 (iorb,xnorm)
          write (*,1115) iorn(iorb),bond(iorb),gusym(iorb),xnorm
1115      format(i4,1x,a8,a1,3x,e22.16,2e16.2)
       endif
#endif
    enddo

    ! initialize Coulomb and exchange potentials
    call initCoulomb
    call initExchange

    deallocate (phf1)
    deallocate (phf2)

  end subroutine initQRHF

  ! ### initHyd ###
  ! 
  !     This routine initializes molecular orbitals as linear combinations
  !     of hydrogenic functions on centres A and B.  In the case of HF or
  !     HFS calculations Coulomb (exchange) potentials are approximated as
  !     a linear combination of Thomas-Fermi (1/r) potentials of the two
  !     centres; if method OED is chosen the potential functions are set
  !     to zero.
  !
  subroutine initHyd 
    use params
    use blas
    use commons
    use discrete
    use sharedMemory
    use coulombExchangePotentials
    use utils

    implicit none

    integer (KIND=IPREC) :: iorb,l1,l2,m1,m2,n1,n2,ngorb
    real (PREC) ez1,ez2
    real (PREC), dimension(:), pointer :: f4,psi

    f4=>supplptr(i4b(9):)

    ! Initialization of molecular orbitals
    
    print *,'... initializing orbitals from hydrogenic functions ...'

    ! loop over orbitals

    do iorb=1,norb
       if (inhyd(iorb).eq.0) then
          write(*,'(5x,"skipping  ",i4,1x,a8," (Is the LCAO entry correct?)")') iorn(iorb),bond(iorb)
          cycle
       endif
       !write(*,'(5x,i4,1x,a8)') iorn(iorb),bond(iorb)
       
       psi=>orbptr(i1b(iorb):)
       !call zeroArray(mxsize,psi(i1beg))
       call zeroArray(mxsize,psi)

       !if     (co1(iorb).ne.zero) then
          if ( abs(co1(iorb))>epsilon(zero)) then
          n1=mgx(1,iorb)
          l1=mgx(2,iorb)
          m1=mgx(3,iorb)
          ez1=eza1(iorb)
          call hydrogenOrbA(scratchptr,n1,l1,m1,ez1)
          !call daxpy(mxsize,co1(iorb),scratchptr,ione,psi(i1beg),ione)
          call daxpy(mxsize,co1(iorb),scratchptr,ione,psi,ione)
       endif

       
       !if (co2(iorb).ne.zero) then
          if ( abs(co2(iorb))>epsilon(zero)) then
          n2=mgx(4,iorb)
          l2=mgx(5,iorb)
          m2=mgx(6,iorb)
          ez2=eza2(iorb)
          call hydrogenOrbB(scratchptr,n2,l2,m2,ez2)
          !call daxpy (mxsize,co2(iorb),scratchptr,ione,psi(i1beg),ione)
          call daxpy (mxsize,co2(iorb),scratchptr,ione,psi,ione)
       endif
    enddo

    ! initialize Coulomb and exchange potentials
    call initCoulomb
    call initExchange
  end subroutine initHyd

  ! ### hydrogenOrbA ###
  !
  !     This routine returns an arrys containing a hydrogenic orbital
  !     centred on centre A for a given set of n,l,m and Z.
  !
  subroutine hydrogenOrbA (orbital,n1,l1,m1,ez1)
    use params
    use discrete
    use commons
    use utils

    implicit none
    integer (KIND=IPREC) :: igp,imu,in,inioff,n1,l1,m1

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
          !write(*,'(i5,5e12.4)') igp,fn1,plaguer(n1,l1,ez1,r1t),shn1,plegendg(l1,m1,costh1)
       enddo
    enddo
  end subroutine hydrogenOrbA
  
  subroutine hydrogenOrbB (orbital,n2,l2,m2,ez2)
    use params
    use discrete
    use commons
    use utils

    implicit none
    integer (KIND=IPREC) :: igp,imu,in,inioff,n2,l2,m2

    real (PREC) :: costh2,ez2,fn2,r2t,rr,shn2,z
    real (PREC), dimension(*) :: orbital

    ! normalization factor for Laguere polynomials
    fn2=laguerre_normalization(ez2,n2,l2)
    ! normalization factor for spherical harmonics
    shn2=sphharm_normalization(l2,m2)

    do in=1,nni
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          igp=inioff+in

          ! for each grid point, e.i. for (vmu(imu),vni(ini))
          ! determine its distance |_r1| and |_r2| from the nuclei A
          ! and B and cosine of the polar angles costh1 and costh2
          ! between z axis and the vectors _r1 and _r2
          rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
          z=(r/2.0_PREC)*vxi(imu)*veta(in)
          r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

          if (r2t.lt.precis) then
             !              jk 02/01
             costh2=-one
          else
             costh2=-(z-r/2.0_PREC)/r2t
             ! Sanity check
             if (costh2 < -1.0_PREC) costh2 = -1.0_PREC
             if (costh2 > 1.0_PREC) costh2 = 1.0_PREC
          endif

          ! calculate radial part of the hydrogenic orbital centered
          ! on the nucleus
          orbital(igp)=fn2*plaguer(n2,l2,ez2,r2t)*shn2*plegendg(l2,m2,costh2)
       enddo
    enddo
  end subroutine hydrogenOrbB

  ! ### laguerre_normalization ###
  !
  !     Normalization factor for Laguerre polynomials
  !
  function laguerre_normalization(z,n,l) result(fn)
    use params
    use commons
    use utils
    implicit none

    real(PREC), intent(in) :: z
    integer (KIND=IPREC),intent(in) :: n, l
    real(PREC) :: fn

    fn=(2.0_PREC*z/dble(n))**(3.0_PREC/2.0_PREC+dble(l))*sqrt(factor(n+l)/(2.0_PREC*dble(n)*factor(n-l-1)))/factor(2*l+1)
  end function laguerre_normalization

  ! ### sphharm_normalization ###
  !
  !     Normalization factor for spherical harmonics
  !
  function sphharm_normalization(l,m) result(shn1)
    use params
    use utils
    implicit none

    integer (KIND=IPREC),intent(in) :: l, m
    real(PREC) :: shn1

    if(m.eq.0) then
       shn1=1.0_PREC/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l+1.0_PREC)*factor(l-m)/factor(l+m))
    else
       shn1=(-1.0_PREC)**dble(m)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l+1.0_PREC)*factor(l-m)/factor(l+m))
    end if
  end function sphharm_normalization

  ! ### plaguer ###
  !
  !    Evaluates and returns a value of the Laguerre polynomial
  !
  function plaguer(n,l,cz,cr)
    use params
    use utils

    implicit none
    integer (KIND=IPREC) :: n,l,lp,np
    real (PREC) :: plaguer
    real (PREC) :: cz,cr,x

    np=n-l-1
    lp=2*l+2
    x=cz*cr/dble(n)
    plaguer=cr**dble(l)*exp(-x)*hypg1f1(-np,lp,x+x)

  end function plaguer
  
end module initOrbitalsPotentials
