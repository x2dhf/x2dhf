! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module scfUtils
  use params
  implicit none
  logical :: lmpoleInit,lorbitalLoop

!  logical :: lenergyThld,lnormThld
  logical, dimension(maxorb) :: lenergyThldFast,lnormThldFast

!  integer (KIND=IPREC) :: ibonus,iepoch,nepochs
!  parameter (ibonus=3,iepoch=10,nepochs=5000)
!  integer (KIND=IPREC), dimension(maxorb) :: iorbiter
  integer (KIND=IPREC) :: inde,indn
  real (PREC), dimension(maxorb) :: deltaee,deltaNorm,eeOld
  real (PREC) :: deltaeeMax,deltaNormMax
  real (PREC) :: orbEnergyThld,orbNormThld
  integer (KIND=IPREC), dimension(maxorb) :: energyThldHit,normThldHit
  
contains
  subroutine printOrbData(iorb)
    use commons
    integer (KIND=IPREC) :: iorb
    
    ! In every SCF iteration the following data are printed for every orbital:
    !   SCF iteration number
    !   orbital label
    !   orbital energy
    !   change of the orbital energy between the current and previous SCF iterations
    !   error of the orbital normalization constant
    !   the largest of overlap integrals calculated between the orbital and lower lying ones
 
    !   verboseLevel - level of output during scf process
    !       =1 at every nobackup SCF iteration only total energy is printed
    !       =2 at every SCF iteration orbital energy, orbital energy
    !          and normalization convergence of and non-orthogonality for
    !          the worst converged orbitals are printed
    !       =3 at every SCF iteration orbital energy, orbital energy
    !          and normalization convergence and non-orthogonality for
    !          each orbital are printed
    !       =4 additionaly, at every SCF iteration the orbital and normalization
    !          convergence rates for the worst converged orbital are printed
    
   if (mod(iscf,abs(saveScfData))==0) then    
       if (verboseLevel==2) then    
          if (iprint16.eq.0) then
             write(iout6,15110) iscf,iorn(inde),bond(inde),gusym(inde),&
                  ee(inde,inde),deltaeeMax,orbNorm(inde)-1.0_PREC
             write(iout6,15110) iscf,iorn(indn),bond(indn),gusym(indn),&
                  ee(indn,indn),deltaee(indn),orbNorm(indn)-1.0_PREC
          else
             write(iout6,25100) iscf,iorn(inde),bond(inde),gusym(inde),&
                  ee(inde,inde),deltaeeMax,orbNorm(inde)-1.0_PREC
             write(iout6,25100) iscf,iorn(indn),bond(indn),gusym(indn),&
                  ee(indn,indn),deltaee(indn),orbNorm(indn)-1.0_PREC
          endif
       endif
    endif
    
    if (lorbitalLoop.and.verboseLevel>=3) then
       if (iprint16.eq.0) then
          write(iout6,15100) iscf,iorn(iorb),bond(iorb),gusym(iorb),&
               ee(iorb,iorb),deltaee(iorb),orbNorm(iorb)-1.0_PREC,wstorthog(iorb)
       else
          write(iout6,25100) iscf,iorn(iorb),bond(iorb),gusym(iorb),&
               ee(iorb,iorb),deltaee(iorb),orbNorm(iorb)-1.0_PREC,wstorthog(iorb)
       endif
    endif

    if (.not.lorbitalLoop.and.verboseLevel==4) then
       if (iprint16.eq.0) then
          write(iout6,'(44x,1Pe16.2,1Pe16.2)') deltaeeMax,deltaNormMax
       else
          write(iout6,'(63x,1Pe44.2,1Pe44.2)') deltaeeMax,deltaNormMax
       endif
    endif
    
15100 format(i5,i4,1x,a8,a1,2x,1Pe23.16,1Pe16.2,1Pe16.2,1Pe16.2)
15110 format(i5,i4,1x,a8,a1,3x,1Pe23.16,1Pe16.2,1Pe16.2,i5)
25100 format(i5,i4,1x,a8,a1,2x,1Pe44.32,1Pe16.2,1Pe16.2,1Pe16.2)

  end subroutine printOrbData


  subroutine fastSCF(iorb)
    use commons
    integer (KIND=IPREC) :: iorb
    
    do iorb=1,norb
       if (ifixorb(iorb)==1) cycle
       
       if (.not.lenergyThldFast(iorb)) then
          if (abs(deltaee(iorb)).lt.orbEnergyThld) then
             lenergyThldFast(iorb)=.true.
             energyThldHit(iorb)=0
          endif
       endif
       
       if (.not.lnormThldFast(iorb)) then
          if (abs(deltaNormMax).lt.orbNormThld) then
             lnormThldFast(iorb)=.true.
             normThldHit(iorb)=0
          endif
       endif
       
       if (abs(deltaee(iorb)).lt.orbEnergyThld) then
          if (energyThldHit(iorb)==nscfextra) then
             ifixorb(iorb)=1
             !write(*,'(2x,"... orbital energy threshold reached for orbital",i4,1x,a8,a1)') &
             !     iorn(iorb),bond(iorb),gusym(iorb)
          else
             energyThldHit(iorb)=energyThldHit(iorb)+1
          endif
       endif
       
       if (abs(deltaNorm(iorb)).lt.orbNormThld) then
          if (normThldHit(iorb)==nscfextra) then
             ifixorb(iorb)=1
             !write(*,'(2x,"... orbital energy threshold reached for orbital",i4,1x,a8,a1)') &
             !     iorn(iorb),bond(iorb),gusym(iorb)
          else
             normThldHit(iorb)=normThldHit(iorb)+1
          endif
       endif
    enddo
 
  end subroutine fastSCF
  

  ! ### printTotalEnergy ####
  !
  !     Prints total energy and its components
  !
  subroutine printTotalEnergy
    use params
    use commons

    implicit none

    if (iprint16.eq.0) then
       write(iout6,6110)  etot
       write(iout6,6100)  evt

       write(iout6,6120)  virrat
       if (lpotGauss.and.iprint(590).ne.0) then
          write(iout6,6130)  etotFN
       endif
    else
       write(iout6,7110)  etot
       write(iout6,7100)  evt
       
       write(iout6,7120)  virrat
       if (lpotGauss.and.iprint(590).ne.0) then
          write(iout6,7130)  etotFN
       endif
    endif
    return

6100 format(5x,'total electronic energy: ',1Pe28.16)
6110 format(5x,'total energy:            ',1Pe28.16)
6120 format(5x,'virial ratio:            ',1Pe28.16)
6130 format(5x,'total energy + FN:       ',1Pe28.16)

7100 format(5x,'total electronic energy: ',1Pe44.32)
7110 format(5x,'total energy:            ',1Pe44.32)
7120 format(5x,'virial ratio:            ',1Pe44.32)
7130 format(5x,'total energy + FN:       ',1Pe44.32)

  end subroutine printTotalEnergy

  
  ! ### potAsympt ###
  !
  !     determines the asymptotic values of Coulomb and exchange
  !     potentials needed for a given Fock equation
  !
  subroutine potAsympt (iorb1,pot,excp)
    use params
    use commons
    use discrete    
    implicit none
    integer (KIND=IPREC) :: i3beg,ibeg,idel,iorb1,iorb2,k

    real (PREC), dimension(*) :: pot,excp

    !   determine values of Coulomb potentials in the asymptotic region
    !   from the multipole expansion

    if (itouch(iorb1).eq.1) then
       ibeg=i2b(iorb1)
       call coulAsympt(iorb1,pot(ibeg))
    endif

    if (DFT.or.HFS.or.SCMC) return    

    !   initialize amul array needed to evaluate exchange potential in
    !   the asymptotic region by calling zasyx

    do iorb2=1,norb
       if (itouch(iorb1).eq.0.and.itouch(iorb2).eq.0) goto 10
       k=k2(iorb1,iorb2)
       !      check if k=0
       if (k.eq.0) goto 10
       if (ilc(k).ne.0) then
          i3beg=i3b(k)
          idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
          if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)

          call exchAsympt (idel,k,excp(i3beg))

          if (ilc(k)==2) then
             idel=mgx(6,iorb2)+mgx(6,iorb1)
             k=k+norb*(norb+1)/2
             i3beg=i3beg+mxsize
             call exchAsympt (idel,k,excp(i3beg))
          endif
       endif
10     continue
    enddo
  end subroutine potAsympt

  subroutine potAsymptExp (iorb1,excp)
    use params
    use commons
    use discrete
    
    implicit none
    integer (KIND=IPREC) :: i3beg,ibeg,idel,iorb1,iorb2,k,ngrid

    real (PREC), dimension(*) :: excp

    !   determine values of Coulomb potentials in the asymptotic region
    !   from the multipole expansion

    if (itouch(iorb1).eq.1) then
       ibeg=i2b(iorb1)
       call coulAsympt(iorb1,excp(ibeg))
    endif

    if (DFT.or.HFS.or.SCMC) return        

    !   initialize amul array needed to evaluate exchange potential in
    !   the asymptotic region by calling zasyx

    do iorb2=1,norb
       if (itouch(iorb1).eq.0.and.itouch(iorb2).eq.0) cycle
       k=k2(iorb1,iorb2)
       ! check if k=0
       if (k.eq.0) cycle
       if (ilc(k).ne.0) then
          i3beg=i3b(k)
          idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
          if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)
          call exchAsympt (idel,k,excp(i3beg))
          if (ilc(k)==2) then
             idel=mgx(6,iorb2)+mgx(6,iorb1)
             k=k+norb*(norb+1)/2
             i3beg=i3beg+mxsize
             call exchAsympt (idel,k,excp(i3beg))
          endif
       endif
    enddo
  end subroutine potAsymptExp

  
  ! ### coulAsympt ###
  !
  !     mpole multipole moments (stored in cmulti) and associated Legendre
  !     functions are used to calculate for each i=1,nni 4 boundary values
  !     of Coulomb potential for a given orbital ($\tilde{V}^a_C$)),
  !     i.e. its values for j=mxnmu-3, mxnmu-2, mxnmu-1,mxnmu
  !
  subroutine coulAsympt(iorb,pot)
    use params
    use discrete
    use scfshr
    use commons
    
    implicit none
    integer (KIND=IPREC) :: i,iorb,itt,j,kk
    real (PREC) :: costh,rr,rr1
    real (PREC), dimension(*) :: pot

    do j=mxnmu-3,mxnmu
       itt=(j-1)*nni
       do i=1,nni
          kk=i+itt
          rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
          rr1=1.0_PREC/(rr*r2)
          costh=veta(i)*vxi(j)/rr
          pot(kk)=r2*vxi(j)*(rr1+vcoul(iorb,i,j,costh))
       enddo
    enddo
  end subroutine coulAsympt

  ! ### exchAsympt ###
  !
  !     Evaluates asymptotic values of a given exchange potential from the
  !     multipole expansion.
  !
  !     excptail array is used to provide asymptotic values for 'immersed'
  !     exchange potentials (see fill).
  !
  !     For odd values of q in P(k,q) (here q=idel) there should be no
  !     (-1)**q factor in eq. (19). That is why this factor is missing in
  !     this routine.
  !
  subroutine exchAsympt (idel,ipc,excp)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,idel,ipc,itt,j,kk

    real (PREC), dimension(*) :: excp
    real (PREC), dimension(maxmpole) :: pe

    do j=mxnmu-3,mxnmu
       itt=(j-1)*nni
       do i=1,nni
          kk=i+itt
          call vexch(i,j,idel,ipc,pe)
          excp(kk)= r2*vxi(j)*pe(mpole)
       enddo
    enddo
  end subroutine exchAsympt

  ! ### setCi ###
  !
  !     Determines Ci symmetry of a given orbital (or set of orbitals) and
  !     uses it to replace the values in the [pi/2,pi] region by these
  !     from the [0,pi/2] one.
  
  subroutine setCi (iorb,psi)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ihsym,ibeg,iorb,iorb1,iorb2,n2,nmut
    real (PREC), dimension(*) :: psi

    !   ihsym = 1 - symmetry g
    !   ihsym =-1 - symmetry u

    if (iorb.eq.0) then
       iorb1=1
       iorb2=norb
    else
       iorb1=iorb
       iorb2=iorb
    endif

    n2=nni/2+1

    do i=iorb1,iorb2
       ibeg = i1b(i)
       nmut = i1mu(i)
       ihsym=ihomo(i)
       call setCiOrb (n2,nmut,psi(ibeg),ihsym)
       
#ifdef PRINT
! print= 55: setCi: iorb ihomo(iorb) ihsym
       if (iprint(55).ne.0) then
          write(*,*) 'setCi: iorb ihomo(iorb) ihsym',iorn(i),bond(i),gusym(i),ihomo(i),ihsym
       endif
#endif
    enddo

  end subroutine setCi

  ! ### setCiOrb ###
  !
  !     Sets Ci symmetry of a given orbital: the values from the [0,pi/2] region replace
  !     those from [pi/2,pi] one
  !
  subroutine setCiOrb (n2,nmut,orb,ihsym)
    use params
    use discrete
    use commons

    implicit none
    logical, parameter :: lHomoSymInput=.true.
    integer (KIND=IPREC) :: ihsym,mi,mis,n2,ni,nis,nmut
    real (PREC), dimension(nni,nmut) :: orb

    if (lHomoSymInput) then
       ! set symmetry according to input data
       if (ihsym==1) then
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
       return
    endif

    ! check symmetry at a given point and impose the same symmetry at other ones
    mis=10
    nis=2

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
  end subroutine setCiOrb
  
  ! ### mpoleMoments ###
  !
  !     Recalculates multipole moment expansion coefficients every time demax(1),
  !     i.e. maximum error in orbital energy, is reduced by recalcMMfactor. Multipole
  !     moments used to calculate asymptotic values of Coulomb and exchange potentials are
  !     stored in cmulti and exc(di|qu|oc|he|5-8) arrays, respectively.
  !
  !     Coefficients and then asymptotic values are recalculated only for orbitals which
  !     undergo relaxation, i.e. those being touched (itouch=1).
  !
  subroutine mpoleMoments
    use params
    use commons
    use dateTime
    use sharedMemory
    
    implicit none
    integer (KIND=IPREC) :: iorb1,iorb2
    integer (KIND=IPREC) :: count1,count2,countRate
    real (PREC) :: time1,time2
    real (PREC), dimension(:), pointer :: cw_orb,cw_suppl,cw_sctch

    if (nel==1) return
    if (OED) return
    
    cw_orb=>orbptr
    cw_suppl=>supplptr
    cw_sctch=>scratchptr
    
    call getCpuTime(time1)
    call getRealTime (count1,countRate)

    call coulMom

    if (nel>1.and.(HF.or.TED)) then
       call exchMom
    endif

    call getCpuTime (time2)
    tmomen =tmomen + (time2-time1)

    call getRealTime (count2,countRate)
    tmomenReal=tmomenReal+dble(count2-count1)/dble(countRate)

    if (.not.lmpoleInit.and.verboseLevel>2) then
       write(iout6,'(" ... multipole moment expansion coefficients recalculated ...")')
    endif

  end subroutine mpoleMoments

  ! ### coulMom ###
  !
  !     Calculates multipole moments up to order k_max=7.
  !
  subroutine coulMom
    use params
    use discrete
    use scfshr
    use commons
    use utils
    use blas
    use sharedMemory
    
    implicit none
    integer (KIND=IPREC) :: i,ibeg,iorb,mu,n,ni
    real (PREC) :: costh,rr,xr,xw
    real (PREC), dimension(10) ::  dome
    real (PREC), dimension(:), pointer :: psi,f4,wgt2,wk1,wk2
    real (PREC), dimension(:), pointer :: dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8
    
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    psi=>orbptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)

    dd1 => legendreptr(         1:   mxsize)
    dd2 => legendreptr(  mxsize+1: 2*mxsize)
    dd3 => legendreptr(2*mxsize+1: 3*mxsize)
    dd4 => legendreptr(3*mxsize+1: 4*mxsize)
    dd5 => legendreptr(4*mxsize+1: 5*mxsize)
    dd6 => legendreptr(5*mxsize+1: 6*mxsize)
    dd7 => legendreptr(6*mxsize+1: 7*mxsize)
    dd8 => legendreptr(7*mxsize+1: 8*mxsize)

    ! dd1, dd2, dd3, ... are equal to P1*r, P2*r**2 and P3*r**2, ...,
    ! respectively (r==xr), where P1, P2, P3, ...  are Legendre polynomials
    ! of a specified order.
    
    ! Prepare integrands and calculate moments. Multiplication by f4 is
    ! necessary because the factor 2/(r*cosh(mu)) is incorporated in wgt2.

    do iorb=1,norb
       if (itouch(iorb).eq.0) cycle
       ibeg = i1b (iorb)
    
       call prod2 (mxsize,psi(ibeg:),psi(ibeg:),wk1)
       call prod  (mxsize,f4,wk1)

       if (lhomonucl .and. .not.lbreakCi) then
          cmulti (iorb) = 0.0_PREC
       else
          call prod2 (mxsize,dd1,wk1,wk2)
          cmulti (iorb)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif

       call prod2 (mxsize,dd2,wk1,wk2)
       cmulti (iorb+norb)=ddot(mxsize,wgt2,ione,wk2,ione)

       if (mpole.ge.3) then
          if (lhomonucl .and. .not.lbreakCi) then
             cmulti (iorb+2*norb)=0.0_PREC
          else
             call prod2 (mxsize,dd3,wk1,wk2)
             cmulti (iorb+2*norb)=ddot(mxsize,wgt2,ione,wk2,ione)
          endif
       endif

       if (mpole.ge.4) then
          call prod2 (mxsize,dd4,wk1,wk2)
          cmulti (iorb+3*norb) =ddot (mxsize,wgt2,ione,wk2,ione)
       endif

       if (mpole.ge.5) then
          if (lhomonucl .and. .not.lbreakCi) then
             cmulti (iorb+4*norb)=0.0_PREC
          else
             call prod2 (mxsize,dd5,wk1,wk2)
             cmulti (iorb+4*norb)=ddot(mxsize,wgt2,ione,wk2,ione)
          endif
       endif

       if (mpole.ge.6) then
          call prod2 (mxsize,dd6,wk1,wk2)
          cmulti (iorb+5*norb)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif

       if (mpole.ge.7) then
          if (lhomonucl .and. .not.lbreakCi) then
             cmulti (iorb+6*norb)=0.0_PREC
          else
             call prod2 (mxsize,dd7,wk1,wk2)
             cmulti (iorb+6*norb)=ddot(mxsize,wgt2,ione,wk2,ione)
          endif
       endif

       if (mpole.ge.8) then
          call prod2 (mxsize,dd8,wk1,wk2)
          cmulti (iorb+7*norb)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
    enddo
    
#ifdef PRINT
! print=161: multipole moments for orbital 
       if (iprint(161).ne.0) then
          write(*,*)
          write(*,*) 'momen0 - multipole moments for orbital ',iorb
          write(*,1111) (cmulti (iorb+i*norb),i=0,mpole-1)
1111      format(4e25.14)
       endif
#endif       

  end subroutine coulMom

  ! ### exchMom ###
  !
  !     Recalculates multipole moment expansion coefficients every time demax(1),
  !     i.e. maximum error in orbital energy, is reduced by recalcMMfactor. Multipole
  !     moments used to calculate asymptotic values of Coulomb and exchange potentials are
  !     stored in cmulti and exc(di|qu|oc|he|5-8) arrays, respectively.
  !
  !     Coefficients and then asymptotic values are recalculated only for orbitals which
  !     undergo relaxation, i.e. those being touched (itouch=1).
  !
  subroutine exchMom
    use blas
    use discrete
    use params
    use commons
    use dateTime
    use sharedMemory
    use scfshr
    use utils

#ifdef OPENMP    
    use omp_lib
#endif
    
    implicit none

    integer (KIND=IPREC) :: iorb1,iorb2
    integer (KIND=IPREC) :: i,ibeg1,ibeg2,idel,ipc,istart,istop,j,&
         mu,ni,maxThreads,nthread,nthreads
    real (PREC) ::  xrr,xw
    real (PREC), dimension(10) :: dome
!    real (PREC), dimension(nni,mxnmu) :: d1(nni,mxnmu),d2(nni,mxnmu),d3(nni,mxnmu),d4(nni,mxnmu),&
!         d5(nni,mxnmu),d6(nni,mxnmu),d7(nni,mxnmu),d8(nni,mxnmu)


    real (PREC), dimension(:,:), allocatable :: d1,d2,d3,d4,d5,d6,d7,d8

    real (PREC), dimension(:), pointer :: psi,f4,wgt2
    
#ifdef OPENMP
    real (PREC), dimension(:), allocatable :: wk1,wk2
#else
    real (PREC), dimension(:), pointer :: wk1,wk2
#endif

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    psi=>orbptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)

    allocate (d1(nni,mxnmu))
    allocate (d2(nni,mxnmu))    
    allocate (d3(nni,mxnmu))
    allocate (d4(nni,mxnmu))    
    allocate (d5(nni,mxnmu))
    allocate (d6(nni,mxnmu))    
    allocate (d7(nni,mxnmu))
    allocate (d8(nni,mxnmu))    
    
#ifdef OPENMP
    allocate (wk1(mxsize8))
    allocate (wk2(mxsize8))
    maxThreads=OMP_get_max_threads()
#else
    wk1 =>scratchptr(          1:   mxsize8)
    wk2 =>scratchptr(   mxsize8+1: 2*mxsize8)
#endif

    nthreads=min(nexchmm,maxThreads)
    !$OMP  PARALLEL NUM_THREADS(nthreads) DEFAULT(SHARED) &
    !$OMP& PRIVATE(i,ibeg1,ibeg2,idel,iorb1,iorb2,ipc,istart,istop,j,mu,&
    !$OMP&         nthread,ni,xrr,xw,dome,d1,d2,d3,d4,d5,d6,d7,d8,wk1,wk2)
    
#ifdef OPENMP

    nthread=OMP_get_thread_num()
    nthread=nthread+1
    !print *,"OPENMP-1",nexchmm,maxThreads,nthread
    if (nexchmm<=maxThreads) then
       istart=nthread
       istop=nthread
    else
       istart=(nthread-1)*(nexchmm/maxThreads+1)+1
       istop=nthread*(nexchmm/maxThreads+1)
       if (istop>nexchmm) istop=nexchmm
    endif

    if (istart>nexchmm) goto 9999
    
    ! do i=(nthread-1)*(nexchmm/maxThreads+1)+1,nthread*(nexchmm/maxThreads+1)
    do i=istart,istop
       if (i>nexchmm) exit
       
       ! if (i>nexchmm) then
       !    print *,"OPENMP3-",nexchmm,nthread,istart,istop,i
       !    exit
       ! endif
#else
    do i=1,nexchmm
#endif       
       idel=idelmm(i)
       iorb1=iorb1mm(i)
       iorb2=iorb2mm(i)    
       ipc=ipcmm(i)

       !print *,"OPENMP2-",nthread,i,ipc
       
       ibeg1=i1b (iorb1)
       ibeg2=i1b (iorb2)

       do mu=1,mxnmu
          do ni=1,nni
             xrr=r2*sqrt (vxisq(mu)+vetasq(ni)-1.0_PREC)
             call mulex(ni,mu,idel,dome)
             
             xw=xrr
             d1(ni,mu)=dome(1)*xw
             xw=xw*xrr
             d2(ni,mu)=dome(2)*xw
             
             if (mpole.ge.3) then
                xw=xw*xrr
                d3(ni,mu)=dome(3)*xw
             endif
             
             if (mpole.ge.4) then
                xw=xw*xrr
                d4(ni,mu)=dome(4)*xw
             endif
             
             if (mpole.ge.5) then
                xw=xw*xrr
                d5(ni,mu)=dome(5)*xw
             endif
             
             if (mpole.ge.6) then
                xw=xw*xrr
                d6(ni,mu)=dome(6)*xw
             endif
             
             if (mpole.ge.7) then
                xw=xw*xrr
                d7(ni,mu)=dome(7)*xw
             endif
             
             if (mpole.ge.8) then
                xw=xw*xrr
                d8(ni,mu)=dome(8)*xw
             endif
             
#ifdef PRINT
             ! print=165: xrr,xw,dome
             if (iprint(165).ne.0) then
                if ((iorb1.eq.1).and.(iorb2.eq.1).and.(ni.eq.nni/2).and.(mu.eq.mxnmu/2)) then
                   write(*,'(10e16.7)') xrr,xw
                   write(*,'(10e16.7)') (dome(j),j=1,8)
                   print *,'iorb1,iorb2,ni,mu,nnu,mxnmu',iorb1,iorb2,ni,mu,nni,mxnmu
                endif
             endif
#endif
          enddo
       enddo

       ! Now d1,d2,... are equal to N*P(1,q)*r, N*P(2,q)*r**2, ...,
       ! respectively. P(1,q), P(2,q),...  are associate Legendre polynomials
       ! of a specified order (k) and degree q (q==idel). N is equal to
       ! [(k-|q|)!/(k+|q|)!]^{1/2}.
       ! Note that N*N is present in eq. (19). One N is taken care of when
       ! evaluating Q_kq (Q_{k\Delta m}^{ab}) and the other when generating
       ! Pkm (P_{k|\Delta m|).
       
       ! Prepare integrands and calculate moments. Multiplication by f4 is
       ! necessary because of the factor 2/(r*cosh(mu)) incorporated in wgt2.
       
       call prod2 (mxsize,psi(ibeg1:),psi(ibeg2:),wk1)
       call prod  (mxsize,f4,wk1)
       
       call prod2 (mxsize,d1,wk1,wk2)
       excdi(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       
       call prod2 (mxsize,d2,wk1,wk2)
       excqu(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       
       if (mpole.ge.3) then
          call prod2 (mxsize,d3,wk1,wk2)
          excoc(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
       if (mpole.ge.4) then
          call prod2(mxsize,d4,wk1,wk2)
          exche(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
       if (mpole.ge.5) then
          call prod2 (mxsize,d5,wk1,wk2)
          exc5(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
       if (mpole.ge.6) then
          call prod2 (mxsize,d6,wk1,wk2)
          exc6(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
       if (mpole.ge.7) then
          call prod2 (mxsize,d7,wk1,wk2)
          exc7(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
       if (mpole.ge.8) then
          call prod2 (mxsize,d8,wk1,wk2)
          exc8(ipc)=ddot(mxsize,wgt2,ione,wk2,ione)
       endif
       
#ifdef PRINT
       ! print=166: excdi,excqu,excoc,exche,exc5,exc6,exc7,exc8
       if (iprint(166).ne.0) then
          !write(*,'("i,iorb1,iorb2,idel,ipc",4i5)') i,iorb1,iorb2,idel,ipc
          write(*,1000) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb2),bond(iorb2),gusym(iorb2),idel
1000      format(/,i4,1x,a8,a1,3x,i4,1x,a8,a1,3x,i5)
          write(*,1010) excdi(ipc),excqu(ipc),excoc(ipc),exche(ipc),exc5(ipc),exc6(ipc),exc7(ipc),exc8(ipc)
1010      format(4d25.14)
       endif
#endif

#ifdef OPENMP
    enddo
9999 continue    
    deallocate (wk1)
    deallocate (wk2)

    !$OMP END PARALLEL
#else
    enddo
#endif

  end subroutine exchMom

  ! ### mulex ###
  !
  !     This routine generates values of the associate Legendre functions
  !     P(k,q) multiplied by [(k-|q|)!/(k+|q|)!]^{1/2} for k=1,..,8 and
  !     q=mt.  See also routine which employs the same definition of the
  !     functions multi.
  !
  !     Note that the original routine (mulex) generates values scaled by
  !     (-1)**mt
  !
  subroutine mulex(i,j,mt,dome)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,mt,n
    real (PREC) :: costh,costh2,costh4,costh6,costh8,rr,sini,sini2,sini3
    real (PREC), dimension(10) :: dome

    dome(1)=0.0_PREC

    if (mt.gt.4) return

    rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
    costh=0.0_PREC
    if (abs(rr).gt.precis) then
       costh=veta(i)*vxi(j)/rr
    endif


    if (mt.eq.0) then
       !        m=0
       dome(1)=costh
       dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
       do n=2,mpole-1
          dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
       enddo
       return

    elseif (mt.eq.1) then
       !        m=1

       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)
       sini =sqrt(sini2)
       !        P_1^1
       dome(1)=-sini/sqrt(2.0_PREC)
       !        P_2^1
       dome(2)=-3.0_PREC*sini*costh/sqrt(6.0_PREC)

       if (mpole.lt.3) return
       !        P_3^1
       dome(3)=(3.0_PREC/2.0_PREC)*sini*(1.0_PREC-5.0_PREC*costh2)/sqrt(12.0_PREC)

       if (mpole.lt.4) return
       !        P_4^1
       dome(4)=(5.0_PREC/2.0_PREC)*sini*costh*(3.0_PREC-7.0_PREC*costh2)/sqrt(20.0_PREC)

       if (mpole.lt.5) return
       !        P_5^1
       costh4=costh2*costh2
       dome(5)=(15.0_PREC/8.0_PREC)*sini*(-1.0_PREC+14.0_PREC*costh2-21.0_PREC*costh4)/sqrt(30.0_PREC)
       if (mpole.lt.6) return
       !        P_6^1
       dome(6)=(21.0_PREC/8.0_PREC)*sini*costh*(-5.0_PREC+30.0_PREC*costh2-33.0_PREC*costh4)/sqrt(42.0_PREC)

       if (mpole.lt.7) return
       !        P_7^1
       costh6=costh4*costh2
       dome(7)=(7.0_PREC/16.0_PREC)*sini*(5.0_PREC-135.0_PREC*costh2+495.0_PREC*costh4-429.0_PREC*costh6)/sqrt(56.0_PREC)

       if (mpole.lt.8) return
       !        P_8^1
       dome(8)=(9.0_PREC/16.0_PREC)*sini*costh*(35.0_PREC-385.0_PREC*costh2+1001.0_PREC*costh4-715.0_PREC*costh6)/sqrt(72.0_PREC)
       return

    elseif (mt.eq.2) then
       !        m=2
       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)

       !        P_1^2
       dome(1)=0.0_PREC
       !        P_2^2
       dome(2)=3.0_PREC*sini2/sqrt(24.0_PREC)

       if (mpole.lt.3) return
       !        P_3^2
       dome(3)=15.0_PREC*costh*sini2/sqrt(120.0_PREC)

       if (mpole.lt.4) return
       !        P_4^2
       costh4=costh2*costh2
       dome(4)=(15.0_PREC/2.0_PREC)*(-1.0_PREC+8.0_PREC*costh2-7.0_PREC*costh4)/sqrt(360.0_PREC)

       if (mpole.lt.5) return
       !        P_5^2
       dome(5)=(105.0_PREC/2.0_PREC)*costh*(-1.0_PREC+4.0_PREC*costh2-3.0_PREC*costh4)/sqrt(840.0_PREC)

       if (mpole.lt.6) return
       !        P_6^2
       costh6=costh4*costh2
       dome(6)=(105.0_PREC/8.0_PREC)*(1.0_PREC-19.0_PREC*costh2+51.0_PREC*costh4-33.0_PREC*costh6)/sqrt(1680.0_PREC)

       if (mpole.lt.7) return
       !        P_7^2
       dome(7)=(63.0_PREC/8.0_PREC)*costh*(15.0_PREC-125.0_PREC*costh2+253.0_PREC*costh4-143.0_PREC*costh6)/sqrt(3024.0_PREC)

       if (mpole.lt.8) return
       !        P_8^2
       costh8=costh6*costh2
       dome(8)=(315.0_PREC/8.0_PREC)*(-1.0_PREC+34.0_PREC*costh2-176.0_PREC*costh4+ &
            286.0_PREC*costh6-143.0_PREC*costh8)/sqrt(5040.0_PREC)
       return

    elseif (mt.eq.3) then
       !        m=3

       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)
       sini =sqrt(sini2)
       sini3=sini2*sini

       dome(1)=0.0_PREC
       dome(2)=0.0_PREC

       if (mpole.lt.3) return
       !        P_3^3
       dome(3)=-15.0_PREC*sini3/sqrt(720.0_PREC)

       if (mpole.lt.4) return
       !        P_4^3
       dome(4)=-105.0_PREC*sini3*costh/sqrt(5040.0_PREC)

       if (mpole.lt.5) return
       !        P_5^3
       dome(5)=-(105.0_PREC/2.0_PREC)*sini3*(-1.0_PREC+9.0_PREC*costh2)/sqrt(20160.0_PREC)

       if (mpole.lt.6) return
       !        P_6^3
       dome(6)=-(315.0_PREC/2.0_PREC)*sini3*costh*(-3.0_PREC+11.0_PREC*costh2)/sqrt(60480.0_PREC)

       if (mpole.lt.7) return
       !        P_7^3
       costh4=costh2*costh2
       dome(7)=-(315.0_PREC/8.0_PREC)*sini3*(3.0_PREC-66.0_PREC*costh2+143.0_PREC*costh4)/sqrt(151200.0_PREC)

       if (mpole.lt.8) return
       dome(8)=-(3465.0_PREC/8.0_PREC)*sini3*costh*(3.0_PREC-26.0_PREC*costh2+39.0_PREC*costh4)/sqrt(332640.0_PREC)
       return

    elseif (mt.eq.4) then
       !        m=4

       dome(1)=0.0_PREC
       dome(2)=0.0_PREC
       dome(3)=0.0_PREC

       if (mpole.lt.4) return
       !        P_4^4
       costh2=costh*costh
       costh4=costh2*costh2
       dome(4)=105.0_PREC*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(40320.0_PREC)

       if (mpole.lt.5) return
       !        P_5^4
       dome(5)=945.0_PREC*costh*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(362880.0_PREC)

       if (mpole.lt.6) return
       !        P_6^4
       costh6=costh4*costh2
       dome(6)=(945.0_PREC/2.0_PREC)*(-1.0_PREC+13*costh2-23.0_PREC*costh4+11.0_PREC*costh6)/sqrt(1814400.0_PREC)

       if (mpole.lt.7) return
       !        P_7^4
       dome(7)=(3465.0_PREC/2.0_PREC)*costh*(-3.0_PREC+19*costh2-29.0_PREC*costh4+13.0_PREC*costh6)/sqrt(6652800.0_PREC)

       if (mpole.lt.8) return
       !        P_8^4
       costh8=costh6*costh2
       dome(8)=(10395.0_PREC/8.0_PREC)*(1.0_PREC-28.0_PREC*costh2+118.0_PREC*costh4- &
            156.0_PREC*costh6+65.0_PREC*costh8)/sqrt(19958400.0_PREC)
       return
    endif

  end subroutine mulex

  ! ### vcoul ###
  !
  !     Evaluates and returns the boundary value of the Coulomb potential
  !     for a given orbital at a given point ($\tilde{V}^a_C$))
  !
  function vcoul(iorb,i,j,costh)
    use params
    use discrete
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: i,iorb,j,kxk,m,n
    real (PREC) :: vcoul
    real (PREC) :: costh,rr,rr1,rr2,pe
    real (PREC), dimension(maxmpole) :: dome

    rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
    rr1=1.0_PREC/(rr*r2)

    dome(1)=costh
    dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.5_PREC
    do n=2,mpole-1
       dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
    enddo

    pe=0.0_PREC
    rr2=rr1
    do m=1,mpole
       kxk=iorb+(m-1)*norb
       pe=pe+cmulti(kxk)*dome(m)*(rr1**dble(m+1))
    enddo
    vcoul=pe

  end function vcoul

  ! ### vexch ###
  !
  !     Calculates the value of exchange potential from the multipole
  !     expansion.
  !
  subroutine vexch(i,j,mt,ipc,pe)
    use params
    use discrete
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ipc,j,mt,n
    real (PREC) :: costh,costh2,costh4,costh6,costh8,rr,sini,sini2,sini3,xrn
    real (PREC), dimension(maxmpole) :: dome,pe

    pe(1)=0.0_PREC

    if (mt.gt.4) return

    rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
    costh=0.0_PREC
    if (abs(rr).gt.precis) then
       costh=veta(i)*vxi(j)/rr
    endif

    rr=rr*r2

    if (mt.eq.0) then
       ! m=0
       dome(1)=costh
       dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
       do n=2,mpole-1
          dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
       enddo

       ! for m not equal zero the associated Legendre polynomials
       ! are multiplied by sqrt((k-|q|)!/(k+|q|)!) factor

    elseif (mt.eq.1) then
       ! m=1
       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)
       sini =sqrt(sini2)

       dome(1)=-sini/sqrt(2.0_PREC)
       dome(2)=-3.0_PREC*sini*costh/sqrt(6.0_PREC)

       if (mpole.lt.3) goto 100
       dome(3)=(3.0_PREC/2.0_PREC)*sini*(1.0_PREC-5.0_PREC*costh2)/sqrt(12.0_PREC)

       if (mpole.lt.4) goto 100
       dome(4)=(5.0_PREC/2.0_PREC)*sini*costh*(3.0_PREC-7.0_PREC*costh2)/sqrt(20.0_PREC)

       if (mpole.lt.5) goto 100
       costh4=costh2*costh2
       dome(5)=(15.0_PREC/8.0_PREC)*sini*(-1.0_PREC+14.0_PREC*costh2-21.0_PREC*costh4)/sqrt(30.0_PREC)

       if (mpole.lt.6) goto 100
       dome(6)=(21.0_PREC/8.0_PREC)*sini*costh*(-5.0_PREC+30.0_PREC*costh2-33.0_PREC*costh4)/sqrt(42.0_PREC)

       if (mpole.lt.7) goto 100
       costh6=costh4*costh2
       dome(7)=(7.0_PREC/16.0_PREC)*sini&
            *(5.0_PREC-135.0_PREC*costh2+495.0_PREC*costh4-429.0_PREC*costh6)/sqrt(56.0_PREC)
       if (mpole.lt.8) goto 100
       dome(8)=(9.0_PREC/16.0_PREC)*sini*costh&
            *(35.0_PREC-385.0_PREC*costh2+1001.0_PREC*costh4-715.0_PREC*costh6)/sqrt(72.0_PREC)
       !      minus sign for odd values of mt
       !       pe=-pe

    elseif (mt.eq.2) then
       !      m=2
       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)

       dome(1)=0.0_PREC
       dome(2)=3.0_PREC*sini2/sqrt(24.0_PREC)

       if (mpole.lt.3) goto 100
       dome(3)=15.0_PREC*costh*sini2/sqrt(120.0_PREC)

       if (mpole.lt.4) goto 100
       costh4=costh2*costh2
       dome(4)=(15.0_PREC/2.0_PREC)*(-1.0_PREC+8.0_PREC*costh2-7.0_PREC*costh4)/sqrt(360.0_PREC)

       if (mpole.lt.5) goto 100
       dome(5)=(105.0_PREC/2.0_PREC)*costh*(-1.0_PREC+4.0_PREC*costh2-3.0_PREC*costh4)/sqrt(840.0_PREC)

       if (mpole.lt.6) goto 100
       costh6=costh4*costh2
       dome(6)=(105.0_PREC/8.0_PREC)*&
            (1.0_PREC-19.0_PREC*costh2+51.0_PREC*costh4-33.0_PREC*costh6)/sqrt(1680.0_PREC)

       if (mpole.lt.7) goto 100
       dome(7)=(63.0_PREC/8.0_PREC)*costh&
            *(15.0_PREC-125.0_PREC*costh2+253.0_PREC*costh4-143.0_PREC*costh6)/sqrt(3024.0_PREC)

       if (mpole.lt.8) goto 100
       costh8=costh6*costh2
       dome(8)=(315.0_PREC/8.0_PREC)*(-1.0_PREC+34.0_PREC*costh2-176.0_PREC*costh4+&
            286.0_PREC*costh6-143.0_PREC*costh8)/sqrt(5040.0_PREC)

    elseif (mt.eq.3) then
       ! m=3
       costh2=costh*costh
       sini2=abs(1.0_PREC-costh2)
       ! FIXME sini =(sini2) or sini =sqrt(sini2)
       sini =sini2
       sini3=sini2*sini

       dome(1)=0.0_PREC
       dome(2)=0.0_PREC

       if (mpole.lt.3) goto 100
       dome(3)=-15.0_PREC*sini3/sqrt(720.0_PREC)

       if (mpole.lt.4) goto 100
       dome(4)=-105.0_PREC*sini3*costh/sqrt(5040.0_PREC)

       if (mpole.lt.5) goto 100
       dome(5)=-(105.0_PREC/2.0_PREC)*sini3*(-1.0_PREC+9.0_PREC*costh2)/sqrt(20160.0_PREC)

       if (mpole.lt.6) goto 100
       dome(6)=-(315.0_PREC/2.0_PREC)*sini3*costh*(-3.0_PREC+11.0_PREC*costh2)/sqrt(60480.0_PREC)

       if (mpole.lt.7) goto 100
       costh4=costh2*costh2
       dome(7)=-(315.0_PREC/8.0_PREC)*&
            sini3*(3.0_PREC-66.0_PREC*costh2+143.0_PREC*costh4)/sqrt(151200.0_PREC)

       if (mpole.lt.8) goto 100
       dome(8)=-(3465.0_PREC/8.0_PREC)*&
            sini3*costh*(3.0_PREC-26.0_PREC*costh2+39.0_PREC*costh4)/sqrt(332640.0_PREC)

    elseif (mt.eq.4) then
       ! m=4
       dome(1)=0.0_PREC
       dome(2)=0.0_PREC
       dome(3)=0.0_PREC

       if (mpole.lt.4) goto 100
       costh2=costh*costh
       costh4=costh2*costh2
       dome(4)=105.0_PREC*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(40320.0_PREC)

       if (mpole.lt.5) goto 100
       dome(5)=945.0_PREC*costh*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(362880.0_PREC)

       if (mpole.lt.6) goto 100
       costh6=costh4*costh2
       dome(6)=(945.0_PREC/2.0_PREC)*&
            (-1.0_PREC+13*costh2-23.0_PREC*costh4+11.0_PREC*costh6)/sqrt(1814400.0_PREC)

       if (mpole.lt.7) goto 100
       dome(7)=(3465.0_PREC/2.0_PREC)*&
            costh*(-3.0_PREC+19*costh2-29.0_PREC*costh4+13.0_PREC*costh6)/sqrt(6652800.0_PREC)

       if (mpole.lt.8) goto 100
       costh8=costh6*costh2
       dome(8)=(10395.0_PREC/8.0_PREC)*(1.0_PREC-28.0_PREC*costh2+118.0_PREC*costh4 &
            -156.0_PREC*costh6+65.0_PREC*costh8)/sqrt(19958400.0_PREC)
    endif

100 continue

    xrn=rr*rr
    pe(1)=        dome(1)*excdi(ipc)/xrn
    xrn=xrn*rr
    pe(2)=pe(1) + dome(2)*excqu(ipc)/xrn

    if (mpole.ge.3) then
       xrn=xrn*rr
       pe(3)=pe(2) + dome(3)*excoc(ipc)/xrn
    endif

    if (mpole.ge.4) then
       xrn=xrn*rr
       pe(4)=pe(3) + dome(4)*exche(ipc)/xrn
    endif

    if (mpole.ge.5) then
       xrn=xrn*rr
       pe(5)=pe(4) + dome(5)*exc5(ipc)/xrn
    endif

    if (mpole.ge.6) then
       xrn=xrn*rr
       pe(6)=pe(5) + dome(6)*exc6(ipc)/xrn
    endif

    if (mpole.ge.7) then
       xrn=xrn*rr
       pe(7)=pe(6) + dome(7)*exc7(ipc)/xrn
    endif

    if (mpole.ge.8) then
       xrn=xrn*rr
       pe(8)=pe(7) + dome(8)*exc8(ipc)/xrn
    endif
  end subroutine vexch

  ! ### alphaSCMC ###
  !
  !     Calculates self-consistent multiplicative constant; see Eq. 3 in
  !     Karasiev, Ludenia, eq.3 PR 64 (2002) 062510
  !
  subroutine alphaSCMC 
    use blas
    use commons
    use discrete
    use dftexc
    use diskInterface
    use exchContribs
    use params
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) ::  ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex1,&
         iorb,iorb1,iorb2,ipe1,ipe2,&
         isiorb,isiorb1,isiorb2,isipot,isipot1,isipot2,nmut,nmut1,nmut2

    real (PREC) :: ehfex,eps,oc,oc1,oc2,ocx1,ocx2,w,wdcoul,wex1,wex2
    real (PREC), dimension(:), pointer :: psi,excp,e,f0,wgt1,wgt2,&
              wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

    data eps /1.e-7_PREC/
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    wk6 =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    wk7 =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    wk8 =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    wk9 =>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10=>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11=>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12=>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13=>scratchptr(13*mxsize8+1:14*mxsize8)

    if (nel.lt.2) return

    ! contribution from coulomb interaction within the same shell

    wdcoul=0.0_PREC

    do iorb=1,norb
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)

       if (oc.gt.1.0_PREC) then
          call dcopy (isipot,excp(ibpot:),ione,wk2,ione)
          call prod  (isiorb,psi(iborb:),wk2)
          call prod  (isiorb,psi(iborb:),wk2)

          w=ddot(isiorb,wgt2,ione,wk2,ione)
          wdcoul=wdcoul+oc/2.0_PREC*w
       endif
    enddo

    !   contribution from coulomb and exchange interaction between shells
    wex1 =0.0_PREC
    wex2 =0.0_PREC

    do iorb1=1,norb
       iborb1=i1b(iorb1)
       isiorb1=i1si(iorb1)
       nmut1=i1mu(iorb1)
       ibpot1=i2b(iorb1)
       isipot1=i2si(iorb1)
       oc1=occ(iorb1)

       ipe1=mgx(6,iorb1)
       iex1 =k2(iorb1,iorb2)

       call dcopy (isiorb1,psi(iborb1:),ione,wk0,ione)
       call prod  (isiorb1,psi(iborb1:),wk0)

       ! calculate exchange interaction within pi, delta, etc. open shell
       if (ipe1.gt.0.and.abs(oc1-1.0_PREC).gt.eps) then
          call exint (iorb1,ocx1)
          ibex=i3b(iex1)
          call prod2  (mxsize,wk0,excp(ibex:),wk1)
          w=ddot (mxsize,wgt2,ione,wk1,ione)
          wex1=wex1+ocx1*w

#ifdef PRINT
! print= 82: orb v(x-within shell) 
          if (iprint(82).ne.0) then
             write(*,6020) iorb1,w
6020         format(' ','orb v(x-within shell) ',i4,2x,4e15.5)
          endif
#endif

#ifdef PRINT
! print= 83: orb v(x-within shell) 
          if (iprint(83).ne.0) then
             write(*,7020) iorn(iorb1),bond(iorb1),gusym(iorb1),ipe1,ocx1
7020         format(' ','v(x-within shell)  ',i4,1x,a8,a1,i4,f6.2)
          endif
#endif          
       endif

       do iorb2=iorb1+1,norb
          iborb2=i1b(iorb2)
          isiorb2=i1si(iorb2)
          nmut2=i1mu(iorb2)
          ibpot2=i2b(iorb2)
          isipot2=i2si(iorb2)
          oc2=occ(iorb2)

          ipe2=mgx(6,iorb2)
          ibex=i3b(k2(iorb1,iorb2))

          ! exchange interaction between shells (same lambda)

          call excont (iorb1,iorb2,ocx1,ocx2)
          call prod2 (mxsize,excp(ibex:),psi(iborb1:),wk1)
          call prod  (mxsize,psi(iborb2:),wk1)
          w=ddot(mxsize,wgt2,ione,wk1,ione)
          wex1=wex1+ocx1*w
#ifdef PRINT
! print= 82: exchange interaction between shells (same lambda)
          if (iprint(82).ne.0) then
             print *,'ibex,wex1,ocx1,w',ibex,wex1,ocx1,w
          endif
#endif          
          ! exchange interaction between shells (different lambda)

          if (ilc2(iorb1,iorb2)==2) then
             call prod2 (mxsize,excp(ibex+mxsize:),psi(iborb1:),wk1)
             call prod  (mxsize,psi(iborb2:),wk1)
             w=ddot (mxsize,wgt2,ione,wk1,ione)
             wex2=wex2+ocx2*w
#ifdef PRINT
! print= 83: exchange interaction between shells (different lambda)
             if (iprint(83).ne.0) then
                print *,'wex2,ocx2,w',wex2,ocx2,w
             endif
#endif
          endif
       enddo
    enddo

    ehfex=wdcoul+wex1+wex2

    alphaf=two/three
    !   exchange energy from DFT functionals
    call dftex (psi,excp,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)

    alphaf=abs(ehfex/edftex)*two/three

#ifdef PRINT
! print= 84: scmc: ehfex,edftex,alphaf
    if (iprint(84).ne.0) then
       write(*,'("scmc: ehfex,edftex,alphaf",3e15.6)') ehfex,edftex,alphaf
    endif
#endif

#ifdef PRINT
! print= 85: scmc: wdcoul,wex1,wex2
    if (iprint(85).ne.0) then
       write(*,'("scmc: wdcoul,wex1,wex2",3e15.6)') wdcoul,wex1,wex2
    endif
#endif
    
  end subroutine alphaSCMC

  ! ### dftex ###
  !
  !     Calculates exchange contribution due to various DFT functionals
  !
  subroutine dftex (psi,pot,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13)
    use params
    use discrete
    use commons
    use utils

    use blas
    use dftexc

    implicit none
    integer (KIND=IPREC) :: i,iborb,ibpot,iorb,isiorb,isipot
    real (PREC) :: oc,w,wex,wndc
    real (PREC), dimension(*) :: psi,pot,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown, &
         grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    wndc=0.0_PREC
    wex=0.0_PREC

    ! calculate the coulomb potential contribution from all orbitals
    ! (include 1/2 factor )

    do i=1,mxsize
       wk2(i)=0.0_PREC
    enddo

    do iorb=1,norb
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)/two
       call daxpy (isipot,oc,pot(ibpot),ione,wk2,ione)
    enddo

    ! contribution from the Coulomb interaction
    do iorb=1,norb
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       call prod2 (isiorb,psi(iborb),psi(iborb),wk0)
       call prod (isiorb,wk2,wk0)
       call dscal (isiorb,occ(iorb),wk0,ione)
       w=ddot(isiorb,wgt2,ione,wk0,ione)
       wndc=wndc+w
    enddo

    ! DFT exchange energy corrections
    edftex=0.0_PREC
    if     (idftex.eq.1) then
       edftex=exxalpha(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.2) then
       edftex=exbe88(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,   &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.3) then
       edftex=expw86(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,   &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    endif

  end subroutine dftex

end module scfUtils
