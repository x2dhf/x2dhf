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
! ### fockDFT ###

! Calculates the exchange potentials as the local Slater approximation.

module fockDFT_m
  implicit none
contains
  subroutine fockDFT(iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,&
       wk11,wk12,wk13)
    use params
    use discret
    use memory
    use scf
    use commons8
    use util

    use blas_m
    use exocc_m
    use fdftpot_m
    use fbe88_m
    use flypcs_m
    use fvwncs_m
    use fpw86_m
    use multf4_m
    use zeroArray_m

    
    implicit none
    integer :: i,iborb,iborb1,ibpot,ibpot1,iorb,iorb1,ipc,ngorb,ngorb1,ngpot,ngpot1,norb2,isiorb1

    real (PREC) :: const13,ocdown,ocup,tmpf,w
    real (PREC),dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,&
         wk11,wk12,wk13

    parameter (const13=1.0_PREC/3.0_PREC)

    ! contributions from one-electron terms
    iborb=i1b(iorb)
    ngorb=i1si(iorb)
    ibpot=i2b(iorb)
    ngpot=i2si(iorb)
    norb2=norb*norb

    call zeroArray(mxsize,fock1)
    call zeroArray(mxsize,fock2)
    ! (jk)
    call zeroArray(mxsize,excp(length3-mxsize))
    call zeroArray(mxsize,wk2)

    !   store wk0 in fock2 for a time being
    !    call copy (mxsize,wk0,ione,fock2,ione)
    ! use fock2 as wk0

    !   contributions from the local exchange approximation

    if (nel.gt.1) then
       
       ! store wk0 in fock2 for a time being
       ! call copy (mxsize,wk0,ione,fock2,ione)
       ! use fock2 as wk0
       
       ! contributions from the local exchange approximation
       select case (idftex) 
       case (1) ! LDA
          do iorb1=1,norb
             if (inhyd(iorb1).eq.1) cycle
             iborb1 =i1b (iorb1)
             isiorb1=i1si(iorb1)
             call exocc (iorb1,ocup,ocdown)
             call prodas (isiorb1,ocup , psi(iborb1),psi(iborb1),fock2)
             call prodas (isiorb1,ocdown,psi(iborb1),psi(iborb1),wk2)
          enddo

          do i=1,mxsize
             fock2(i)=(fock2(i))**const13
             wk2(i)=(wk2(i))**const13
          enddo
          call add (mxsize,fock2,wk2)

          ! multiply the local exchange potential by f4 to make it commensurate with the
          ! Coulomb potential

          call multf4(wk2)
          tmpf=fdftpot(alphaf)
          call  scal(mxsize,tmpf,wk2,ione)
       case (2) ! B88
          call fbe88(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       case (3) ! PW86
          write (*,*) 'FIXME'
          call abort
          ! call fpw86(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       case (4) ! PW91
          write(*,'(/"Error: PW91 exchange potential is not supported")')
          stop 'fockDFT'
       case default 
          write(*,'(/"Warning: unsupported exchange potential")')
          stop 'fockDFT'
       end select

       !      store the local exchange potential in a separate array
       call copy (mxsize,wk2,ione,fock1,ione)

       call zeroArray(mxsize,wk2)

       ! add contributions from correlations potentials
       if (idftcorr>0) then
          select case (idftcorr)
          case (1) ! LYP
             call flypcs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
             ! do i=1,mxsize
             !    print *,i,wk2(i)
             ! enddo
             
          case (2) ! VWN
             call fvwncs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
             
          case default 
             write(*,'(/"Warning: unsupported correlation potential")')
             stop 'fockDFT'
          end select
          
          ! add correlation contributions to the local exchange ones
          call add (mxsize,wk2,fock1)
       endif
    endif
    

    ! add local exchange and correlations contributions ones to one-electron contribution
    ! call copy (mxsize,wk2,ione,fock1,ione)
    ! call copy (mxsize,fock2,ione,wk0,ione)
    
    call zeroArray(mxsize,wk2)
    call zeroArray(mxsize,fock2)
    
    if (nel.gt.1) then
       ! add contributions from Coulomb and off-diagonal Lagrange multipliers.
       
       do iorb1=1,norb
          if (inhyd(iorb1).eq.1) cycle
          iborb1=i1b(iorb1)
          ngorb1=i1si(iorb1)
          ibpot1=i2b(iorb1)
          ngpot1=i2si(iorb1)
          
          ipc=iorb1+norb*(iorb-1)
          
          ! in the local exchange approximation the Coulomb potential also includes the
          ! contribution from the orbital
          call axpy (ngpot1,occ(iorb1),pot(ibpot1),ione,wk2,ione)
          
          if (iorb.ne.iorb1.and.engo(ipc).ne.0.0_PREC) then
             do i=1,ngorb1
                fock2(i)=fock2(i)+engo(ipc)*f4(i)*psi(iborb1+i-1)
             enddo
          endif
       enddo
       
       ! store the local exchange potential in exch array as its not used in HFS/DFT (to
       ! be used by EaDFT and EabDFT) at the end of excp array (important when scmc is on)
       call copy (mxsize,fock1,ione,excp(length3-mxsize),ione)
       
       ! add the coulomb potential to the local exchange one
       call add (mxsize,wk2,fock1)
       
       ! multiply coulomb/exchange potentials and off-diagonal Lagrange
       ! multipliers by f2
       call prod (mxsize,f2,fock1)
       call prod (mxsize,f2,fock2)
       ! nel.gt.1
    endif

    call copy (mxsize,f0,ione,wk2,ione)
    call axpy (mxsize,eng(iorb),f1,ione,wk2,ione)
    
    if (mm(iorb).ne.0) then
       !      e enter the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call axpy (mxsize,w,e,ione,wk2,ione)
    endif
    
    ! add Coulomb contributions to one-electron one (containing local
    ! exchange and correlation contributions)
    
    call add (mxsize,wk2,fock1)
    
  end subroutine fockDFT
end module fockDFT_m
