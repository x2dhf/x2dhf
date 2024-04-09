! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module doSCF
  use params
  implicit none
  logical :: lmpoleInit,lorbitalLoop

  logical :: lenergyThld,lnormThld
  logical, dimension(maxorb) :: lenergyThldFast,lnormThldFast

  integer (KIND=IPREC) :: ibonus,iepoch,nepochs
  parameter (ibonus=3,iepoch=10,nepochs=5000)
  integer (KIND=IPREC), dimension(maxorb) :: iorbiter,iscforder
  integer (KIND=IPREC) :: inde,indn
  real (PREC), dimension(maxorb) :: deltaee,deltaNorm,eeOld
  real (PREC) :: deltaeeMax,deltaNormMax,iscf
  real (PREC) :: orbEnergyThld,orbNormThld
  integer (KIND=IPREC), dimension(maxorb) :: energyThldHit,normThldHit
contains

  subroutine prepSCF 
    use blas
    use params
    use discrete
    use commons
    use fock
    use lagrangeMultipliers
    use scfUtils
    use totalEnergy
#ifdef LIBXC                           
    use totalEnergyLXC
#endif    
    use normOrtho
    use diskInterface
    use scfshr
    implicit none
    integer (KIND=IPREC) :: iorb,ipass,jorb

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    do iorb=1,norb
       if (HF) call fockHF (iorb)
       if (TED) call fockHF (iorb)
       
       ! use self-coded exchange-correlation functionals
       if (DFT.or.HFS.or.SCMC) call fockDFT(iorb)
       
#ifdef LIBXC
       ! use exchange-correlation functionals from libxc
       if (LXC) then
          if (lxcPolar) then
             call fockLXCpol(iorb)
          else
             call fockLXCunpol(iorb)
          endif
          !!!!call EaLXC(iorb)
          !call EaDFT(iorb)
       endif
#endif    
    enddo

    ! orbital energies and Lagrange multiplies are initialized
    ! or have been retrieved from the disk file

    ! asymptotic values of Coulomb and exchange potentials are
    ! calculated or have been retrieved from the disk file

    ! see routine initOrbPot for the explanation of the following command

    ! check orhogonality of orbitals

    ! if (lcoulexch) then
    !    call dcopy (norb*mxsize,cw_coul,ione,exchptr,ione)
    ! endif

#ifdef PRINT
! print= 30: prepSCF: checking orthogonality of orbitals
    if (iprint(30).ne.0) then
       write(*,'(1x,/"orthogonality of orbitals:")')
       do iorb=norb,1,-1
          call checkOrtho (.true.,iorb)
       enddo
    endif

#endif

    if (OED.or.initFuncsOED) then
       nel=1
    endif
    
    if (initAddData) then
       do iorb=1,norb
          jorb=norb+1-iorb
          call norm (jorb)
          call ortho (jorb)
       enddo

       do iorb=1,norb
          itouch(iorb)=1
       enddo

       ! initialize multipole moments
       lmpoleInit=.true.
       call mpoleMoments
       write(*,*) '... initializing multipole moment coefficients ...'

       do iorb=norb,1,-1
          if (HF.or.OED.or.TED) then
             call EaHF (iorb)
             eeOld(iorb)=ee(iorb,iorb)
             call EabHF (iorb)
          endif

          if (LXC) then
             call EaLXC(iorb)
             eeOld(iorb)=ee(iorb,iorb)
             call EabLXC(iorb)
          endif

          
          if (DFT.or.HFS.or.SCMC) then
             call EaDFT (iorb)
             eeOld(iorb)=ee(iorb,iorb)
             call EabDFT (iorb)
          endif
       enddo
       write(*,*) '... initializing Lagrange multipliers ...'
    endif

    ipass=0
    do iorb=1,norb
       if (ifixorb(iorb)==1.and.ipass==0) then
          ipass=1
          write(*,*)
          if (lfixorb4init) then
             write(*,'(5x,"orbitals kept frozen until potentials get converged:")')
          else
             write(*,'(5x,"orbitals kept frozen:")')
          endif
          write(*,'(6x,i2,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
       elseif (ifixorb(iorb)==1.and.ipass==1) then
          write(*,'(6x,i2,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
       endif
    enddo
    if (ipass==0.and.lfixorb) write(*,'(6x,"orbitals not relaxed")')     

    if (HF.or.OED) then
       call EtotalHF
    endif

    if (TED) then
       call EtotalTED
    endif

    if (DFT.or.HFS.or.SCMC) then
       call EtotalDFT 
    endif
    iorb=1
    
#ifdef LIBXC
    if (LXC) then
       call EtotalLXC
    endif
#endif

    call printTotalEnergy
    
  end subroutine prepSCF
 
  ! ### SCF ####
  !
  ! Controls SCF process.
  !
  ! In each SCF iteration the orbitals are relaxed in the reversed order
  ! (norb, norb-1,...,1). For each orbital the necessary Coulomb and
  ! exchange potentials are relaxed first. Not all exchange potentials
  ! required by the Fock equation are relaxed but only those which depend
  ! on the orbitals already relaxed in a given SCF iteration. The
  ! relaxation of Coulomb and exchange potentials can be done independently
  ! so this part of the SCF process can be parallelised (via the OpenMP or
  ! pthreads). Subsequently the given orbital is relaxed, normalized and
  ! orthogonalised (the Gram-Schmidt method is used).

  ! When the orbital loop is over the maximum change in the orbital energies (deltaeeMax)
  ! and norms (deltaNormMax) is determined.
  
  ! The multipole moment expansion coefficients are recalculated every time
  ! deltaEE, i.e. maximum error in orbital energy, changes by the
  ! recalcMMfactor factor. In case recalcMMfactor is set to a negative
  ! number these coefficients are not calculated.

  ! If deltaeeMax is less than the orbital energy threshold (orbEnergyThld) or
  ! deltaNormMax is less than the orbital norm threshold (orbNormThld) for the consecutive
  ! nscfExtra iterations, the SCF process terminates (the convergence
  ! criteria are applied in this order).

  ! The SCF proces also terminates when the maximum number of SCF
  ! iterations is reached or the maximum difference between orbital
  ! energies in two consecutive SCF iterations are too large (the SCF
  ! process is considered divergent).
  
  ! Every saveData SCF iteration the orbitals, potentials and some additional data are
  ! save to disk so that the SCF process could be resumed.
  !
  subroutine scf
    use params
    use discrete
    use commons
    use fock
    use lagrangeMultipliers
    use totalEnergy
#ifdef LIBXC    
    use totalEnergyLXC
#endif    
    use dateTime
    use diskInterface
    use normOrtho
    use scfUtils
    use sharedMemory
    use solver
    use diskInterface
    use diskInterfaceMisc
    implicit none

    integer (KIND=IPREC) :: i,iasympt,ibeg,ic,iend,iorb,istat,jorb,&
         modv,nei,next,noenergydec,nonormdec,nresets,nresets4no,nOrbEnergyThld,nOrbNormThld
    real (PREC) :: eeconv,eeconvprev,dnmaxprev,thren4init,time1,time2,tscf

    integer (KIND=IPREC),dimension(:), pointer :: cw_sor
    real (PREC), dimension(:), pointer ::  cw_coul,cw_exch,cw_suppl,cw_sctch,cw_scratch4lxc
    real (PREC), dimension(1000) ::  gDeltaEE
    
    integer (KIND=IPREC) :: mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c

    common /c_interface_17/ mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c

    real (PREC), dimension(:), pointer :: cw_orb
    common /c_orbptr/ cw_orb

    logical :: ldemax0,stop_x2dhf

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    cw_sor=>sorptr
    cw_coul=>exchptr
    cw_exch=>exchptr
    cw_suppl=>supplptr
    cw_sctch=>scratchptr
    cw_scratch4lxc=>scratch4lxcptr
    
    istat=1
    stop_x2dhf=.false.
    ldemax0=.true.
    ic=0
    orbEnergyThld=10.0_PREC**dble(-ienterm)
    orbNormThld=10.0_PREC**dble(-inoterm)
    thren4init=10.0_PREC**dble(-ienterm4init)
    nresets=0
    nresets4no=0
    nOrbEnergyThld=0

    do iorb=1,norb
       lenergyThldFast(iorb)=.false.
       lnormThldFast(iorb)=.false.
    enddo
    
    ! the SCF procedure starts here: set some parameters and go

    eeconv=10.0_PREC
    do i=1,norb
       itouch(i)=0
       demax(i)=100.0_PREC
       deltaeeMax=100.0_PREC
       deltaee(i)=0.0_PREC
       eeOld(i)=ee(i,i)
       deltaNorm(i)=0.0_PREC
       energyThldHit(i)=0
       normThldHit(i)=0
    enddo
    deltaNormMax=100.0_PREC
    iorb=1
    tscf=0.0_PREC
    tortho=0.0_PREC
    trayl =0.0_PREC
    tlagra=0.0_PREC
    tmomen=0.0_PREC
    ttoten=0.0_PREC

    trelaxCpu=0.0_PREC    
    trelaxCPUOrb=0.0_PREC
    trelaxCPUCoul=0.0_PREC
    trelaxCPUExch=0.0_PREC
    
    trelaxRealOrb=0.0_PREC
    trelaxRealCoul=0.0_PREC
    trelaxRealExch=0.0_PREC
    tmomenReal=0.0_PREC    
    lenergyThld=.false.
    lnormThld=.false.
    lmpoleInit=.false.    

    ! SCF/SOR process is stopped if orbital energy differences are
    ! either too large or too small

    eeconvprev=1.e35_PREC
    dnmaxprev=1.e35_PREC
    noenergydec=0
    nonormdec=0

  
    call getCpuTime (time1)
    tscf=time1

    iscf=0

    time1=0.0_PREC
    time2=0.0_PREC

    ! inforce C_i symmetry or orbitals
    if (lhomonucl .and. .not.lbreakCi) then
       call setci(izero,cw_orb)
    endif
    ! prepear an array needed for writing exchange potentials involving a given
    ! orbital

    !   ----- begin of scf loop -----

    
    if (maxscf.le.0) write(iout6,*) '... skipping scf iterations ...'

    do iscf=1,maxscf
       if (iscf.eq.1) then
          if (verboseLevel==1) then
             if (iprint16.eq.0) then
                write(iout6,15050)
             else
                write(iout6,25050)
             endif
          else
             if (iprint16.eq.0) then
                write(iout6,15000)
             else
                write(iout6,25000)
             endif
          endif
       endif

       ! ----- loop over orbitals (in reverse order) -----

       ! relaxation begins from the innermost orbital, i.e. the last one on
       ! the input data list unless iscforder specifies alternative
       ! ordering
       ibeg=1
       iend=norb
       iasympt=norb
       
       do i=ibeg,iend
          lorbitalLoop=.true.
          if (iscforder(i).ne.0) then
             iorb=iscforder(i)
          else
             iorb=norb-i+1
          endif
          
          itouch(iorb)=1
          if (OED.and.ifixorb(iorb).ne.0) cycle
          if (TED.and.ifixorb(iorb).ne.0) cycle
          if (lfastscf.and.ifixorb(iorb).ne.0) cycle          

          call getCpuTime (time1)

          ! determines Coulomb and exchange potentials in the asymptotic region
          ! from the multipole expansion

          ! if (iasympt.ne.0) then
          if (iasympt.ne.0.and.ifixorb(iorb).eq.0) then
             call potAsympt(iorb,cw_coul,cw_exch)
             iasympt=iasympt-1
          endif

          ! in case of near-degenerate orbitals, i.e. when performing
          ! calculations for homonuclear molecules without inforced
          ! symmetry (no homo label), especially when external electric
          ! field is applied, one can improve convergence by changing the
          ! direction of the SOR sweeps (forward/backward sweeps for
          ! even/odd SCF)

          isstart= 1
          isstop = ngrid1
          isstep = 1
          
          if (altSweeps) then
             if (mod(iorb,itwo).eq.0) then
                isstart= 1
                isstop = ngrid1
                isstep = 1
             else
                isstart= ngrid1
                isstop = 1
                isstep =-1
             endif
          endif
          ! these are used in sorpt.c
          isstartc=isstart
          isstopc=isstop
          
          ! perform relaxation of Coulomb and exchange potentials
          ! maxsorpot(iorb) SOR sweeps of each potential and
          ! maxsororb(iorb) sweeps for each orbital

          call relaxDriver (iorb)

          ! normalize the current orbital
          call getCpuTime (time1)
          if (ifixorb(iorb)==0) call norm (iorb)
          call getCpuTime (time2)

          tortho=tortho+(time2-time1)

          ! enforce C_i symmetry of iorb orbital
          if (lhomonucl .and. .not.lbreakCi) then
             call setci(iorb,cw_orb)
          endif

          ! orthogonalize the current orbital
          call getCpuTime (time1)
          if (norb.gt.1.and.ifixorb(iorb)==0) then
             call ortho (iorb)
          endif
          call getCpuTime(time2)
          tortho=tortho+(time2-time1)


          ! calculate diagonal energy value
          call getCpuTime(time1)

          ! if orbitals are kept fixed relaxation of potentials do not
          ! change exchangeptr needed in EaHF
          if (HF.or.OED) then
             if (ifixorb(iorb)/=0) call fockHF(iorb)
             call EaHF (iorb)
          endif

          if (TED) then
             if (ifixorb(iorb)/=0) call fockTED(iorb)
             call EaHF (iorb)
          endif

#ifdef LIBXC          
          if (LXC) then
             if (ifixorb(iorb)/=0) then
                if (lxcPolar) then
                   call fockLXCpol(iorb)
                else
                   call fockLXCunpol(iorb)
                endif
             endif
             if (lxcHyb) then
                call EaDFT (iorb)
                call EaLXC (iorb)
             else
                call EaLXC (iorb)
             endif
          endif
#endif             
          if (DFT.or.HFS.or.SCMC) then
             if (ifixorb(iorb)/=0) call fockDFT(iorb)
             call EaDFT (iorb)
          endif

          call getCpuTime(time2)
          trayl =trayl + (time2-time1)

          deltaee(iorb)=ee(iorb,iorb)-eeOld(iorb)
          
          ! Although it should not normally happen sometimes deltaee(iorb)
          ! becomes zero. In such cases its value is explicitly set to
          ! 2*orbEnergyThld.  (This procedure is especially important for OED
          ! calculations when converged orbitals get 'fixed' status.)
          if (abs(deltaee(iorb))<epsilon(zero)) then
             deltaee(iorb)=two*orbEnergyThld
             if (ldemax0) then
                write(*,'("  Warning: energy diff. set to ",1Pe7.1,"."$)') deltaee(iorb)
                write(*,'(" This message is shown only once!" )')
                ldemax0=.false.
             endif
          endif

          eeOld(iorb)=ee(iorb,iorb)

          ! sumNew=zero
          ! do iorb=1,norb
          !    !sumNew=sumNew+abs(deltaEE(iorb))                                                                    
          !    sumNew=sumNew+abs(deltaEE(iorb)/ee(iorb,iorb))
          ! enddo
          ! sumNew=sumNew/norb
          ! write(*,'("tuneSOR:",i5,2e14.4,i9)') iscf,sumNew,abs(log10(sumNew)),nsor4orb+nsor4pot
          
          ! calculate off-diagonal Lagrange multipliers
          call getCpuTime(time1)
          !if (HF) call EabHF (iorb)
          if (HF.or.OED.or.TED) call EabHF (iorb)
          if (DFT.or.HFS.or.SCMC) call EabDFT(iorb)
          if (LXC) call EabLXC (iorb)

          call getCpuTime(time2)
          tlagra =tlagra + (time2-time1)

          if (verboseLevel>=3) call printOrbData(iorb)

          ! Sometimes a lengthy calculation needs to be stopped earlier.  kill <pid> will
          ! do the job but you risk that the data being written will be lost. To avoid
          ! this and force the program to end gracefully make
          !   touch stop_x2dhf
          ! or simply
          !   ./xhf stop
          ! If this file is detected the program ends.

          inquire (file ='stop_x2dhf',exist=stop_x2dhf)
          if (stop_x2dhf) goto 110

       enddo
       ! ----- end of orbital loop -----

00110  continue

       lorbitalLoop=.false.
       
       inde=1
       deltaeeMax=0.0_PREC
       
       indn=1
       deltaNormMax=0.0_PREC
       
       if (norb.gt.0) then
          do nei=1,norb
             if (abs(deltaee(nei)).ge.abs(deltaeeMax)) then
                deltaeeMax=deltaee(nei)
                inde=nei
             endif
          enddo

          do nei=1,norb
             deltaNorm(nei)=abs(orbNorm(nei)-1.0_PREC)
             if(ifixorb(nei).eq.0) then
                if (deltaNorm(nei).ge.deltaNormMax) then
                   deltaNormMax=deltaNorm(nei)
                   if (abs(deltaNormMax)<epsilon(zero)) deltaNormMax=precis
                   indn=nei
                endif
             endif
          enddo
       endif
       
       ! inde and indn point to the worst converged orbital as far as
       ! orbital energies and norms are concerned, respectively
       
       if (abs(eeconv)<epsilon(zero)) eeconv=deltaeeMax

       if (verboseLevel>=2) call printOrbData(iorb)

       if (SCMC) call alphaSCMC

       ! Sometimes the convergence threshold (for energy or normalization)
       ! is set too low and cannot be satisfied on a given grid and as a
       ! result the SCF/SOR process continues in vain. To avoid these
       ! excessive iterations the convergence process is stopped if either
       ! deltaeeMax or deltaNormMax has not decreased over a given number
       ! (nlast) most recent iterations.

       ! skip nscf2skip iterations before starting the search for
       ! saturation; examine the last nenlast and nnolast iterations in case
       ! of orbital energy and normalization, respectively

       if (iscf.gt.nscf2skip) then
          if (abs(deltaeeMax).ge.abs(eeconvprev)) then
             noenergydec=noenergydec+1
          else
             eeconvprev=deltaeeMax
             noenergydec=0
          endif

          if (abs(deltaNormMax).ge.abs(dnmaxprev)) then
             nonormdec=nonormdec+1
          else
             dnmaxprev=deltaNormMax
             nonormdec=0
          endif

          if (noenergydec.gt.nenlast.or.nonormdec.gt.nnolast) then
             write(*,*) '          '
             write(*,*) '... solution cannot be further improved ...'
             goto 9999
          endif
       else
          eeconvprev=demax(1)
          dnmaxprev=deltaNormMax
       endif

       ! check if the scf process is to be continued
       if (abs(demax(1)).gt.orbEnergyIncThld) then
          write(*,*) '           '
          write(*,*) '... maximal orbital energy change is too large ...'
          goto 9999
       endif

       if (iscf==maxscf) then
          write(*,*) '          '
          write(*,*) '... scf iteration limit reached ... '
          goto 9999
       endif

       ! If fastSCF is on skip relaxing orbitals for which orbital energy
       ! relaxation threshold has been achieved
       if (lfastSCF) call fastSCF(iorb)

       if (lfixorb4init .and. abs(deltaeeMax).lt.thren4init) then
          ! stop SCF process if the orbital energy differences are less
          ! than the threshold for nscfextra consecutive iterations
          if (nOrbEnergyThld==nscfextra) then
             write(*,*) '          '
             write(*,*) '... initial relaxation of potentials with fixed orbitals reached the threshold ...'
             write(*,*) '... SCF continues with both orbitals and potentials being relaxed ...'
             do i=1,norb
                ifixorb(i)=0
             enddo
             lfixorb=.false.
             lfixorb4init=.false.
             eeconv=10.0_PREC
             deltaeeMax=100.0_PREC
             continue
          else
             nOrbEnergyThld=nOrbEnergyThld+1
          endif
       endif

       if (.not.lenergyThld) then
          if (abs(deltaeeMax).lt.orbEnergyThld) then
             lenergyThld=.true.
             eeconv=deltaeeMax
             iasympt=norb
             do iorb=1,norb
                itouch(iorb)=0
             enddo
             nOrbEnergyThld=0
          endif
       else
          if (abs(deltaeeMax).lt.orbEnergyThld) then
             ! stop SCF process if the orbital energy differences are less
             ! than the threshold for nscfExtra consecutive iterations
             if (nOrbEnergyThld==nscfExtra) then
                write(*,*) '          '
                write(*,*) '... orbital energy threshold reached ...'
                goto 9999
             else
                nOrbEnergyThld=nOrbEnergyThld+1
             endif
          elseif (abs(deltaNormMax).lt.orbNormThld) then
             ! stop SCF process if the norm differences are less
             ! than the threshold for nscfextra consecutive iterations
             if (nOrbNormThld==nscfextra) then
                write(*,*) '          '
                write(*,*) '... orbital normalization threshold reached ...'
                goto 9999
             else
                nOrbNormThld=nOrbNormThld+1
             endif
          else
             lenergyThld=.false.
             if (nresets>maxresets) then
                write(*,*) '          '
                write(*,*) '... orbital energies cannot be further improved ...'
                goto 9999
             else
                nresets=nresets+1
             endif
          endif
       endif

       if (.not.lfixorb) then
          if (.not.lnormThld) then
             if (abs(deltaNormMax).lt.orbNormThld) then
                lnormThld=.true.
                iasympt=norb
                do iorb=1,norb
                   itouch(iorb)=0
                enddo
                nOrbNormThld=0
             endif
          else
             if(abs(deltaNormMax).lt.orbNormThld) then
                ! stop SCF process if the orbital norms are less than the
                ! threshold for 3 consecutive iterations
                if (nOrbNormThld==nscfextra) then
                   write(*,*) '          '
                   write(*,*) '... orbital normalization threshold reached ...'
                   goto 9999
                else
                   nOrbNormThld=nOrbNormThld+1
                endif
             else
                lnormThld=.false.
                if (nresets4no>maxresets) then
                   write(*,*) '          '
                   write(*,*) '... orbital norms cannot be further improved ...'
                   goto 9999
                else
                   nresets4no=nresets4no+1
                endif
             endif
          endif
       endif

       if (stop_x2dhf) then
          write(*,*) '          '
          write(*,*) '... stop_x2dhf detected ... program is terminating ...'
          next=0
          goto 9999
       endif

       ! recalculate multipole moments expansion coefficients
       if (abs(eeconv/deltaeeMax)>recalcMMfactor) then
          eeconv=deltaeeMax
          call mpoleMoments
       endif

       iasympt=norb
       
       do iorb=1,norb
          itouch(iorb)=0
       enddo

       ! recalculate multipole coefficients and write functions to a disk
       ! file if saveScfData > 0
       if (saveScfData.gt.0 .and. mod(iscf,saveScfData)==0) then
          ! calculate 'intermediate' total energy
             !if (.not.DFT.and.lxcFuncs==0) call Etotal
          if (HF.or.OED) call EtotalHF
          if (TED) call eTotalTED             
          
#ifdef LIBXC                
          if (LXC) call EtotalLXC
#endif
          if (DFT.or.HFS.or.SCMC) then
             call EtotalDFT
          endif
          
          call getCpuTime(time2)
          ttoten=ttoten+(time2-time1)

          if (verboseLevel>1) then
             if (iprint16.eq.0) then
                write(iout6,'(/7x,"total energy: ",1Pe23.16/)') etot
             else
                write(iout6,'(/7x,"total energy: ",1Pe43.32/)') etot
             endif
          else
             if (iprint16.eq.0) then
                write(iout6,'(i5,16x,1Pe23.16)') iscf,etot
             else
                write(iout6,'(i5,16x,1Pe38.28)') iscf,etot
             endif

          endif
          !$OMP SECTIONS 
          call writeToDisk
          !$OMP END SECTIONS NOWAIT
       endif
       
       ! let's separate pintouts for consecutive SCF iterations  
       if (verboseLevel>=3.and.norb>1) write(iout6,'(1x)')
       
       ! In order to run a job in a batch mode (or in the background) with
       ! the output redirected into a file and to be able to monitor every
       ! single iteration one has to flush the output buffer to disk
       call flush(iout6)
    enddo
    !   ----- end of scf loop	-----

09999 continue

                    
    ! As long as the SCF process is not converged and is carried out in
    ! separate runs one can notice the discrepances between the orbital
    ! energies at the end of one run and the corresponding values
    ! calculated at the onset of the subsequent one (especially when label
    ! fixorb is used). As a remedy one has to recalculate all of the
    ! orbital energies (and off-diagonal Lagrange multipliers if necessary)
    ! with the final set of orbitals and potentials.
    if (norb.gt.1) then
       ! let's keep the norms for the final report generated in printResults 
       deltaee=orbNorm
       do iorb=1,norb
          jorb=norb+1-iorb
          call norm (jorb)
          call ortho (jorb)
       enddo
       orbNorm=deltaee
    endif

    do iorb=norb,1,-1
       if (HF) then
          call EaHF(iorb)
          call EabHF (iorb)
       endif
       
       if (DFT.or.HFS.or.SCMC) then
          call EaDFT(iorb) 
          call EabDFT(iorb)
       endif

       if (LXC) then
          call EaLXC(iorb) 
          call EabLXC(iorb)
       endif
    enddo

  ! write functions to a disk file 
    call mpoleMoments
    call writeToDisk
     
    call getCpuTime(time1)
    if (HF.OR.OED) then
       call EtotalHF
       write(*,*)
       call printTotalEnergy
    endif
    
    if (TED) then
       call eTotalTED
    endif

    
#ifdef LIBXC                
    if (LXC) then
       if (LXC) call EtotalLXC
       write(*,*)
       call printTotalEnergy
    endif
#endif
    if (DFT.or.HFS.or.SCMC) then
       call EtotalDFT
       write(*,*)
       call printTotalEnergy
    endif

    !write(*,'(5x,i6," (MC)SOR iterations")') nsor4orb+nsor4pot
    write(*,'(/5x,"(MC)SOR iterations:",12x,i8)') nsor4orb+nsor4pot    
    
    call getCpuTime(time2)
    ttoten=ttoten+(time2-time1)
 
    if (iout4dd==1) call writeDisk4dd
  
    if (iout4dft==1) call writeDisk4dft
  
    ! Calculate and write to disk the Pauli and von Weizsaecker kinetic potentials
    if (iout4kinpot==1) call writeDisk4kinpot
    return
    
15000 format(/'   scf  orbital',14x,'energy',15x,                 'energy diff.',6x,'1-norm',10x,'overlap',&
           /'   ---  -------',6x,'-----------------------',6x,'------------',5x,'---------',8x,'--------')
  
15050 format(/'   scf         ',10x,'total energy',&
           /'   ---         ',6x,'-----------------------')
  
15100 format(i5,i4,1x,a8,a1,2x,1Pe23.16,1Pe16.2,1Pe16.2,1Pe16.2)
15110 format(i5,i4,1x,a8,a1,3x,1Pe23.16,1Pe16.2,1Pe16.2,i5)
25000 format(/'  scf   orbital',27x,'energy',23x,'energy diff.',6x,'1-norm',10x,'overlap',&
           /'   ---  -------',11x,'---------------------------------------',&
           6x,'------------',5x,'---------',8x,'--------')
25050 format(/'  scf          ',20x,'total energy',&
           /'   ---',18x'-----------------------------------')
25100 format(i5,i4,1x,a8,a1,2x,1Pe44.32,1Pe16.2,1Pe16.2,1Pe16.2)
  end subroutine scf

  ! ### relaxDriver ###
  !
  !     Controls relaxation of Coulomb potentials, exchange potentials
  !     and orbitals.
  !
  !     By default the SOR method is used to relax both orbitals and
  !     potentials. Coulomb and exchange potentials associated with a given
  !     orbital are relaxed in parallel if OpenMP or pthread support is
  !     switched on (via OPENMP and PTHREAD directives, respectively).

  !     When the MCSOR method is selected (via MCSOR|MCSOR-O label) then
  !     the single-threaded version of the multi-colour SOR method is
  !     chosen to relax orbitals. When OPENMP/PTHREAD directive is present
  !     the parallelized version of MCSOR is employed when orbitals are
  !     relaxed. Coulomb and exchange potentials are relaxed in parallel
  !     using single-threaded SOR method. However, when MCSOR-CE label is
  !     used also the Coulomb/exchange potentials are relaxed by means of
  !     the parallelized version of the MCSOR routine. Labels MCSOR-O and
  !     MCSOR-CE were added to force the the usage of the MCSOR selectively
  !     for orbitals and Coulomb potentials, respectively.

  !     Relaxation of orbitals 
  !
  !     SOR: 
  !                       -->  orbSOR     -->  sor
  !        #ifdef OPENMP  -->  orbSOR     -->  sor
  !        #ifdef PTHREAD -->  orbSOR     -->  sor
  !        #ifdef TPOOL   -->  orbSOR     -->  sor  
  !
  !     MCSOR|MCSOR-O: 
  !                       -->  orbMCSOR   -->  mcsor (single-threaded)
  !        #ifdef OPENMP  -->  orbMCSOR   -->  mcsor (multi-threaded via OpenMP directives)
  !        #ifdef PTHREAD -->  orbMCSORPT -->  mcsor_pthread (multi-threaded via pthread)
  !        #ifdef TPOOL   -->  orbMCSORPT -->  mcsor_tpool (multi-threaded via pthread)
  !
  !
  !     Relaxation of Coulomb/exchange potentials 
  !
  !     SOR|MCSOR|MCSOR-O:  
  !                       -->  coulExchSOR    -->  sor
  !        #ifdef OPENMP  -->  coulExchSOR    -->  sor (each potential in a separate OpenMP thread)
  !        #ifdef PTHREAD -->  coulExchSORPT  -->  coulExch_pthread (each potential in a separate pthread)
  !        #ifdef TPOOL   -->  coulExchSORPT  -->  coulExch_ptool (each potential in a separate pthread)
  !
  !     MCSOR-CE: 
  !                       -->  exchMCSOR      -->  mcsor (single-threaded)
  !        #ifdef OPENMP  -->  coulExchMCSOR  -->  mcsor (multi-threaded via OpenMP directives; each 
  !                                                       potential relaxed in a separate OpenMP thread)
  !        #ifdef PTHREAD -->  coulExchSORPT  -->  coulExch_pthread (each potential relaxed in a separate
  !                                                                  pthread via SOR)
  !        #ifdef TPOOL   -->  coulExchSORPT  -->  coulExch_ptool (each potential relaxed in a separate
  !                                                                pthread via SOR)
  !
  !

  subroutine relaxDriver(iorb)
    use params
    use commons
    use dateTime
    use discrete
    use relaxOrbs
    use relaxPots
    
    implicit none
    integer (KIND=IPREC) :: i,iorb
    integer (KIND=IPREC) :: count1,count2,countRate
    real (PREC) time1,time2
    
    ! process the exchange potentials
    call getCpuTime (time1)   
    call getRealTime(count1,countRate)

    if (.not.OED.and.nexchpots(iorb)/=0) then
       if (.not.lfixexch) then
          ! mcsor-ce label present
          if (lpotmcsor) then
#ifdef OPENMP
             ! MCSOR can be used by the performance varies
             call coulExchMCSOR (iorb)
#elif ( defined PTHREAD || defined TPOOL )
             ! SOR must be used to avoid nested parallel regions
             call coulExchSORPT (iorb)
#else
             ! use MCSOR (single-threaded)
             call coulExchMCSOR (iorb)
#endif             
          else
#ifdef OPENMP
             ! every exchange potential is relaxed in separated thread via OpenMP
             call coulExchSOR (iorb)
#elif ( defined PTHREAD || defined TPOOL )
             call coulExchSORPT (iorb)
#else
             call coulExchSOR (iorb)
#endif
          endif
       endif
    endif
    
    call getRealTime (count2,countRate)
    call getCpuTime (time2)   
    trelaxCpuExch=trelaxCpuExch+time2-time1
    trelaxRealExch=trelaxRealExch+dble(count2-count1)/dble(countRate)
  
    call getCpuTime (time1)       
    call getRealTime (count1,countRate)
    if (.not.lfixorb) then
       if (lorbmcsor) then
#if ( defined PTHREAD || defined TPOOL )
          call orbMCSORPT(iorb)
#else
          call orbMCSOR(iorb)
#endif          
       else
          call orbSOR(iorb)
       endif
    endif

    call getRealTime (count2,countRate)
    call getCpuTime (time2)   
    trelaxCpuOrb=trelaxCpuOrb+time2-time1
    trelaxRealOrb=trelaxRealOrb+dble(count2-count1)/dble(countRate)

  end subroutine relaxDriver

end module doSCF
