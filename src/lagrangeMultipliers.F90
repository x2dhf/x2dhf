! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module lagrangeMultipliers
  implicit none

  interface EabDFT
     module procedure EabLXC
  end interface

contains
  ! ### Ea ###
  !
  !     A shortcut evaluation of the eigenvalue of h operator as
  !     <orb_a|h|orb_a> (orbitals are assumed to be normalized). It is
  !     assumed the the Coulomb and exchange contributions have already
  !     been calculated in the fock routine and are stored in coulombptr and
  !     exchangeptr arrays, respectively.
  subroutine EaHF(iorb)
    use params
    use discrete
    use scfshr
    use commons
    use utils
    use blas
    use integrals
    use sharedMemory
    implicit none
    integer (KIND=IPREC) :: iorb,ipc
         
    real (PREC) :: eshift,woneel,wtwoel,wtwoelCoul,wtwoelExch

    real (PREC), dimension(:), pointer :: e,excp,f0,psi,wgt1,wgt2,&
         wk0,wk1,wk2,wk3

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

    woneel=oneelij(iorb,iorb,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)

    if (idebug(807)==1) then
       print *,"Ea:      woneel=",woneel 
       woneel=oneelii4plot(iorb,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)
       print *,"Ea: woneel4plot=",woneel
       stop "Ea"
    endif

    call prod2 (mxsize,psi(i1b(iorb):),coulombptr(i2b(iorb):),wk0)
    call prod (mxsize,psi(i1b(iorb):),wk0)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk0,ione)

    call prod2 (mxsize,psi(i1b(iorb):),exchangeptr(i2b(iorb):),wk0)
    wtwoelExch=-ddot(mxsize,wgt2,ione,wk0,ione)
    wtwoel=wtwoelCoul+wtwoelExch

#ifdef PRINT    
! print= 64: EaHF: Coulomb and exchange energy contributions 
    if (iprint(64).ne.0) then
       write(*,'(4x,"EaHF: energy contributions for",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
       write(*,'("    1-electron-Coul         =",1Pe23.16)') wtwoelCoul
       write(*,'("    2-electron-Exch         =",1Pe23.16)') wtwoelExch
       write(*,'("  1+2-electron              =",1Pe23.16)') wtwoel
    endif
#endif

#ifdef PRINT    
! print= 65: EaHF: Coulomb and exchange energy contributions 
    if(iprint(65).ne.0) then
       !write(*,'(18x,"DFT: 2-electron energy contributions:")')
       write(*,'(4x,"<",i2,1x,a5,a1,"| J |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),wtwoelCoul
       write(*,'(4x,"<",i2,1x,a5,a1,"| K |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),wtwoelExch
       write(*,'(4x,"<",i2,1x,a5,a1,"|J+K|",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),wtwoel
    endif
#endif
    
    ipc=iorb+(iorb-1)*norb
    ee(iorb,iorb)=woneel+wtwoel

#ifdef PRINT        
    if(iprint(64).ne.0.or.iprint(65).ne.0) then
       write(*,'(4x,"<",i2,1x,a5,a1,"|1+2|",i2,1x,a5,a1,"> =",1Pe23.16)') iorn(iorb),bond(iorb),gusym(iorb),&
            iorn(iorb),bond(iorb),gusym(iorb),woneel+wtwoel
    endif
#endif
  end subroutine EaHF

  ! ### EaDFT ###
  !
  !     Evaluates an eigenvalue of the Fock equation for a given orbital
  !     using a DFT/LXC exchange potential
  !
  subroutine EaDFT(iorb)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use discrete
    use inout
    use integrals
    use sharedMemory
    use utils
    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,i1beg,i1beg1,i2beg,i2beg1,ipc,isym,nmut
    real (PREC) :: w,woneel,wtwoel,wtwoelCoul,wtwoelHyb1,wtwoelHyb2
    real (PREC), dimension(:), pointer :: e,excp,f0,psi,psi1,wgt1,wgt2,&
         wk0,wk1,wk2,wk3

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
    
    call zeroArray(mxsize,wk0)
    call zeroArray(mxsize,wk1)

    i1beg=i1b(iorb)
    i2beg=i2b(iorb)
    nmut=i1mu(iorb)

    psi1=>orbptr(i1beg:)
    ipc=iorb+(iorb-1)*norb
    !engo(ipc)=zero
    !ee(iorb,iorb)=zero
    wtwoelCoul=zero
    wtwoelHyb1=zero
    wtwoelHyb2=zero

    woneel=oneelij(iorb,iorb,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)

    call prod2 (mxsize,orbptr(i1b(iorb):),coulombptr(i2b(iorb):),wk0)
    call prod (mxsize,orbptr(i1b(iorb):),wk0)
    wtwoel=ddot(mxsize,wgt2,ione,wk0,ione)

    ee(iorb,iorb)=woneel+wtwoel    
    !engo(ipc) = woneel+wtwoel
    !eng(iorb) = engo(ipc)

#ifdef PRINT    
! print= 64: EaDFT: Coulomb and exchange energy contributions 
    if (iprint(64).ne.0) then
       write(*,'(4x,"EaDFT: energy contributions for",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
       write(*,'("    1-electron              =",1Pe23.16)') woneel
       write(*,'("    2-electron              =",1Pe23.16)') wtwoel
       write(*,'("  1+2-electron              =",1Pe23.16)') woneel+wtwoel
    endif
#endif
    
  end subroutine EaDFT

  ! ### EaDFT ###
  !
  !     Evaluates an eigenvalue of the Fock equation for a given orbital
  !     using a DFT/LXC exchange potential
  !
  subroutine EaLXC(iorb)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use discrete
    use inout
    use integrals
    use sharedMemory
    use utils
    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,i1beg,i1beg1,i2beg,i2beg1,ipc,isym,nmut
    real (PREC) :: w,woneel,wtwoel,wtwoelCoul,wtwoelHyb1,wtwoelHyb2
    real (PREC), dimension(:), pointer :: e,excp,f0,psi,psi1,wgt1,wgt2,&
         wk0,wk1,wk2,wk3
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
    
    call zeroArray(mxsize,wk0)
    call zeroArray(mxsize,wk1)

    i1beg=i1b(iorb)
    i2beg=i2b(iorb)
    nmut=i1mu(iorb)

    psi1=>orbptr(i1beg:)
    ipc=iorb+(iorb-1)*norb
    !engo(ipc)=zero
    ee(iorb,iorb)=zero
    wtwoelCoul=zero
    wtwoelHyb1=zero
    wtwoelHyb2=zero

    woneel=oneelij(iorb,iorb,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)
    wtwoel=twoelijLXC (iorb,iorb,psi,excp,excp,wgt2,wk0,wk1,wk2)

    !engo(ipc) = woneel+wtwoel
    !eng(iorb) = engo(ipc)
    ee(iorb,iorb)=woneel+wtwoel
    
  end subroutine EaLXC



  ! ### EaLxcDirect ###
  !
  !     Evaluates an eigenvalue of the Fock equation for a given orbital
  !     using a DFT/LXC exchange potential via coulombptr and exchangeptr
  !     arrays (see EaDirect).
  !
  subroutine EaLxcDirect(iorb)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use discrete
    use inout
    use integrals
    use sharedMemory
    use utils
    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,i1beg,i1beg1,i2beg,i2beg1,ipc,isym,nmut
    real (PREC) :: w,woneel,wtwoel,wtwoelCoul,wtwoelExch
    real (PREC), dimension(:), pointer :: e,excp,f0,psi,psi1,wgt1,wgt2,&
         wk0,wk1,wk2,wk3

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
    
    call zeroArray(mxsize,wk0)
    call zeroArray(mxsize,wk1)

    i1beg=i1b(iorb)
    i2beg=i2b(iorb)
    nmut=i1mu(iorb)

    psi1=>orbptr(i1beg:)
    ipc=iorb+(iorb-1)*norb
    !engo(ipc)=zero
    ee(iorb,iorb)=zero
    
    wtwoelCoul=zero
    wtwoelExch=zero

    woneel=oneelij(iorb,iorb,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)

    call prod2 (mxsize,psi(i1b(iorb):),coulombptr(i2b(iorb):),wk0)
    call prod (mxsize,psi(i1b(iorb):),wk0)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk0,ione)
    
    call prod2 (mxsize,psi(i1b(iorb):),exchangeptr(i2b(iorb):),wk0)
    wtwoelExch=-ddot(mxsize,wgt2,ione,wk0,ione)
    wtwoel=wtwoelCoul+wtwoelExch

    !print *,"EaLXCDirect: ",wtwoelCoul,wtwoelExch

    
#ifdef PRINT    
! print= 64: EaLxcDirect: Coulomb and exchange energy contributions 
    if (iprint(64).ne.0) then
       write(*,'(4x,"EaLxcDirect: energy contributions for",i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
       write(*,'("    1-electron              =",1Pe23.16)') woneel
       write(*,'("    2-electron              =",1Pe23.16)') wtwoel
       write(*,'("  1+2-electron              =",1Pe23.16)') woneel+wtwoel
    endif
#endif    
    
    ee(iorb,iorb)=woneel+wtwoel
    
  end subroutine EaLxcDirect

  !  ### EabHF ###
  !
  !      Evaluates off-diagonal Lagrange multipliers
  !
  !      ipc12=iorb1+(iorb2-1)*norb
  !      ipc21=iorb2+(iorb1-1)*norb
  !      e(iorb2,iorb1) = e(2,1) = engo(ipc21)=<1|T_k +V_n+V_C-V_x|2>
  !      e(iorb1,iorb2) = e(1,2) = engo(ipc12)=<2|T_k +V_n+V_C-V_x|1>
  !
  subroutine EabHF (iorb)
    use params
    use discrete
    use scfshr
    use commons
    use sharedMemory
    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,ipc12,ipc21
    real (PREC) :: engo12,engo21,engoprv12,engoprv21,ent,wocc
    real (PREC), external :: dot

    real (PREC), dimension(:), pointer :: e,excp,excp1,f0,psi,wgt1,wgt2,&
         wk0,wk1,wk2,wk3
    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(           1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)            

    if (ifixorb(iorb)==1) return
    if (iorb==norb.or.norb==1.or.nel==1) return

    do iorb1=iorb+1,norb
       ! if break is on and the two orbitals have different symmetry
       ! off-diagonal lm is not calculated
       if (.not.offDiagLM(iorb1,iorb)) cycle
       if (.not.loffDiagLM.and.(ihomo(iorb1)*ihomo(iorb)<0)) cycle

       ! store previous values in order to make damping of LM possible
       engoprv12=ee(iorb1,iorb)
       engoprv21=ee(iorb,iorb1)

       if (lmtype==0) then
          ! best convergence (100 SCF iterations for Li, grid=91/35)
          wocc=(occ(iorb1)+occ(iorb))/(occ(iorb1)*occ(iorb))
          call Eab1HF (iorb1,iorb)
          call Eab1HF (iorb,iorb1)
          ent=(ee(iorb1,iorb)+ee(iorb,iorb1))/wocc
          engo12=ent/occ(iorb)
          engo21=ent/occ(iorb1)
       elseif (lmtype==1) then
          ! second best convergence (123 SCF iterations for Li)
          call Eab1HF (iorb1,iorb)
          engo12=ee(iorb1,iorb)
          engo21=ee(iorb1,iorb)/occ(iorb1)
       elseif (lmtype.eq.2) then
          ! worst convergence (1000 SCF iterations for Li)
          call Eab1HF (iorb,iorb1)
          engo21=ee(iorb,iorb1)
          engo12=ee(iorb,iorb1)*occ(iorb1)
       else
          write(iout6,1000)
1000      format(/1x,'... off-diagonal Lagrange multiplies cannot be calculated ...'//)
          stop 'Eab'
       endif

       ee(iorb1,iorb)=engo12
       ee(iorb,iorb1)=engo21

#ifdef PRINT    
! print= 46: EabHF: off-diagonal Lagrange multipliers
       if (iprint(46).ne.0) then
          write(*,'(a8,i4,2e16.6,i4,a8,i4,a8,2e16.8)') 'Eab: ',lmtype,sflagra,dflagra, &
               iorn(iorb),bond(iorb1),iorn(iorb1),bond(iorb1),&
               ee(iorb1,iorb),ee(iorb,iorb1)
       endif
#endif
    enddo
  end subroutine EabHF

  
  ! ### Eab1HF ###
  ! <2|h|1>
  subroutine Eab1HF(iorb1,iorb2)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use discrete
    use integrals
    use inout
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,ihc,iorb,iorbRef,iorb1,iorb2,ipc12,ipc21,isym,kex

    real (PREC) :: w,woneel,wtwoel,wtwoelCoul,wtwoelExch
    real (PREC), dimension(:), pointer :: e,excp,f0,psi,wgt1,wgt2,&
         wk0,wk1,wk2,wk3

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

    wk0 =>scratchptr(           1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)            

    !ipc12=iorb1+(iorb2-1)*norb
    !engo(ipc12)=zero
    ee(iorb1,iorb2)=zero
    
    if (mm(iorb1).ne.mm(iorb2)) return
    if (ige(iorb1).ne.ige(iorb2)) return
    woneel=oneelij(iorb2,iorb1,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)    

    call prod2 (mxsize,psi(i1b(iorb1):),coulombptr(i2b(iorb1):),wk0)
    call prod (mxsize,psi(i1b(iorb2):),wk0)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk0,ione)

    call prod2 (mxsize,psi(i1b(iorb2):),exchangeptr(i2b(iorb1):),wk0)
    wtwoelExch=-ddot(mxsize,wgt2,ione,wk0,ione)
    wtwoel=wtwoelCoul+wtwoelExch
    ee(iorb1,iorb2)=woneel+wtwoel

#ifdef PRINT    
! print= 47: Eab1HF: 1- and 2-electron contributions to off-diagonal Lagrange multipliers
    if (iprint(47).ne.0) then    
       write(*,'(a10,i2,a5,a2,i2,a5,a1,3e12.4)') &
            " Eab1HF: <",iorn(iorb1),bond(iorb1)," |",iorn(iorb2),bond(iorb2),">",&
            woneel,wtwoel,woneel+wtwoel
    endif
#endif
  end subroutine Eab1HF

  ! ### EabLXC ###
  !
  !     Calculates the off-diagonal Lagrange multipliers in case of a
  !     local exchange potential.
  !
  subroutine EabLXC (iorb)
    use discrete
    use params
    use scfshr
    use commons
    use sharedMemory

    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,ipc12,ipc21
    real (PREC) :: engo12,engo21,ent,wocc,engoprv12,engoprv21

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    real (PREC), dimension(:), pointer :: e,excp,f0,psi,wgt1,wgt2,&
         wk0,wk1,wk2,wk3

    ! FIXME
    if (idebug(459).ne.0) return
    if (ifixorb(iorb)==1) return
    if (iorb==norb.or.norb==1.or.nel==1) return    
   
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

    do iorb1=iorb+1,norb
       ! if break is on and the two orbitals have different symmetry off-diagonal lm is
       ! not calculated if (ihomon.eq.2.and.ihomo(iorb1)*ihomo(iorb).lt.0) goto 10
       if (.not.offDiagLM(iorb1,iorb)) cycle
       if (.not.loffDiagLM.and.(ihomo(iorb1)*ihomo(iorb)<0)) cycle

       if (lmtype==0) then
          ee(iorb1,iorb)=zero
          ee(iorb,iorb1)=zero
          
          wocc=(occ(iorb1)+occ(iorb))/(occ(iorb1)*occ(iorb))
          call Eab1DFT (iorb1,iorb)
          call Eab1DFT (iorb,iorb1)
          ent=(ee(iorb1,iorb)+ee(iorb,iorb1))/wocc
          engo12=ent/occ(iorb)
          engo21=ent/occ(iorb1)
       elseif (lmtype==1) then
          call Eab1DFT (iorb1,iorb)
          engo12=ee(iorb1,iorb)
          engo21=ee(iorb1,iorb)/occ(iorb1)
       elseif (lmtype==2) then
          ! works for Li
          !call Eab2DFT (iorb,iorb1)
          call Eab1DFT (iorb1,iorb)
          engo21=ee(iorb,iorb1)
          engo12=ee(iorb,iorb1)*occ(iorb1)
       else
          write(iout6,'(/1x,"... off-diagonal Lagrange multiplies cannot be calculated ..."//)')
          stop 'EabDHF'
       endif

       ee(iorb1,iorb)=engo12
       ee(iorb1,iorb)=engo21

#ifdef PRINT    
! print= 46: EabLXC: off-diagonal Lagrange multipliers
       if (iprint(46).ne.0) then
          write(*,'(a8,i4,e16.6,i4,a8,i4,a8,2e16.8)') 'EabDFT: ',&
               lmtype,sflagra,iorn(iorb),bond(iorb1),iorn(iorb1),bond(iorb1),&
               ee(iorb1,iorb),ee(iorb,iorb1)
       endif
#endif
    enddo

  end subroutine EabLXC

  ! ### Eab1DFT ###
  !
  !     Evaluates off-diagonal Lagrange multipliers in cases of a local exchange
  !     approximation
  !
  !      ipc12=iorb1+(iorb2-1)*norb
  !      ipc21=iorb2+(iorb1-1)*norb
  !      e(iorb2,iorb1) = e(2,1) = engo(ipc21)=<1|T_k +V_n+V_C-V_x|2>
  !      e(iorb1,iorb2) = e(1,2) = engo(ipc12)=<2|T_k +V_n+V_C-V_x|1>
  !
  !      T_k +V_n+V_C-V_x|1> = E_1 + E_12 |2> 
  !      E_12=e(1,2)=<2|T_k +V_n+V_C-V_x|1>
  !
  subroutine Eab1DFT(iorb1,iorb2)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use integrals
    use inout
    use sharedMemory
    use utils

    implicit none
    integer (KIND=IPREC) :: iorbRef,iorb,iorb1,iorb2,ipc1,ipc2,&
         ipc12,ipc21,isym,nmut,kex,i3beg,i1beg,i2beg,ihc,iorbt
    real (PREC), dimension(:), pointer :: e,excp,excp1,excp2,f0,&
         psi,psi1,psi2,psi3,wgt1,wgt2,wk0,wk1,wk2,wk3
    real (PREC) :: w,woneel,wtwoel,oc,wtwoelCoul,wtwoelExch
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    if (mm(iorb1).ne.mm(iorb2)) return
    if (ige(iorb1).ne.ige(iorb2)) return

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    psi=>orbptr    
    psi1=>orbptr(i1b(iorb1):)
    psi2=>orbptr(i1b(iorb2):)

    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)            

    ee(iorb1,iorb2)=zero
    ee(iorb2,iorb1)=zero

    call zeroArray(mxsize,wk0)
    call zeroArray(mxsize,wk1)    

    woneel=oneelij(iorb2,iorb1,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)    

    call prod2 (mxsize,psi(i1b(iorb1):),coulombptr(i2b(iorb1):),wk0)
    call prod (mxsize,psi(i1b(iorb2):),wk0)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk0,ione)

    call prod2 (mxsize,psi(i1b(iorb2):),exchangeptr(i2b(iorb1):),wk0)
    wtwoelExch=-ddot(mxsize,wgt2,ione,wk0,ione)
    wtwoel=wtwoelCoul+wtwoelExch
    ee(iorb1,iorb2)=woneel+wtwoel
    ee(iorb2,iorb1)=ee(iorb1,iorb2)

#ifdef PRINT    
! print= 46: Eab1DFT: off-diagonal Lagrange multipliers
    if (iprint(46).ne.0) then
       write(*,'(a8,i4,e16.6,i4,a8,i4,a8,2e16.8)') 'Eab2DFT: ', &
            lmtype,sflagra,iorn(iorb2),bond(iorb2),iorn(iorb1),bond(iorb1),&
            ee(iorb2,iorb1),ee(iorb1,iorb2)
    endif
#endif    

  end subroutine Eab1DFT

  ! ### Eab2DFT ###
  !
  !     Evaluates the off-diagonal Lagrange multipliers in case of a local
  !     exchange approximation
  !
  subroutine Eab2DFT(iorb1,iorb2)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use sharedMemory
    use utils
    
    implicit none

    integer (KIND=IPREC) :: length
    integer (KIND=IPREC) :: i,iborb,ibpot,iorb1,iorb2,iborb1,iborb2,iorb,ipc12,ipc21,nmut
    real (PREC) :: engo1,engo2,engoprv1,engoprv2,oc,wtwoel
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    real (PREC), dimension(:), pointer :: e,excp,excp1,psi,psi1,psi2,wgt2,&
         wk0,wk1

    if (mgx(6,iorb1).ne.mgx(6,iorb2)) return
    if (ige(iorb1).ne.ige(iorb2)) return

    e=>supplptr(i4b(4):)
    excp=>exchptr
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)

    ee(iorb1,iorb2)=zero
    ee(iorb2,iorb1)=zero

    iborb1=i1b(iorb1)
    iborb2=i1b(iorb2)

    psi1=>orbptr(iborb1:)
    psi2=>orbptr(iborb2:)
    
    ! Add contributions from the Coulomb and local exchange potentials.
    ! In the local exchange approximtion the Coulomb potential includes
    ! also the contribution from a given orbital

    ! Calculate the Coulomb potential contribution from all the orbitals

    call zeroArray (mxsize,wk0)

    do iorb=1,norb
       iborb=i1b(iorb)
       ibpot=i2b(iorb)
       nmut=i1mu(iorb)
       oc=occ(iorb)
       excp1=>exchptr(ibpot:)
       call daxpy (mxsize,oc,excp1,ione,wk0,ione)
    enddo

    call prod  (mxsize,psi2,wk0)

    ! Multiply the local exchange potential by psi(iborb2)
    ! and add the result to the Coulomb potential

    excp1=>exchptr(length2-2*mxsize+1:)
    call prod2 (mxsize,psi2,excp1,wk1)
    call add (mxsize,wk0,wk1)
    
    if (lxcHyb) then
       ! add the (port) of exact exchange potential 
       call add (mxsize,excp1,wk1)
    endif

    ! To complete the integrand wk1 has to be multiplied by psi(iborb1)

    call prod (mxsize,psi1,wk1)
    wtwoel=ddot(mxsize,wgt2,ione,wk1,ione)

    ee(iorb1,iorb2)=wtwoel
    ee(iorb2,iorb1)=wtwoel
    
#ifdef PRINT    
! print= 46: Eab2DFT: off-diagonal Lagrange multipliers
    if (iprint(46).ne.0) then
       write(*,'(a8,i4,e16.6,i4,a8,i4,a8,2e16.8)') 'Eab2DFT: ', &
            lmtype,sflagra,iorn(iorb),bond(iorb1),iorn(iorb1),bond(iorb1),&
            ee(iorb1,iorb2),ee(iorb2,iorb1)
    endif
#endif
  end subroutine Eab2DFT
  
end module lagrangeMultipliers
