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
! ### fockLXC ###

! Calculates exchange and correlation potentials via calls to 
! libxcf90 and libxc libraries.

module fockLXC_m
  implicit none
contains
  subroutine fockLXC(iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,rhot,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,&
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
    use n2f_m
    use nfng_m
    use zeroArray_m

    use xc_f90_types_m
    use xc_f90_lib_m
    
    implicit none
    integer :: i,iborb,iborb1,ibpot,ibpot1,iorb,iorb1,ipc,ngorb,ngorb1,ngpot,ngpot1,norb2,isiorb1

    real (PREC) :: const13,ocdown,ocup,tmpf,w
    real (PREC),dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,rhot,&
         wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info
    integer :: func_id = 1
    
    parameter (const13=1.0_PREC/3.0_PREC)

    ! contributions from one-electron terms
    iborb=i1b(iorb)
    ngorb=i1si(iorb)
    ibpot=i2b(iorb)
    ngpot=i2si(iorb)
    norb2=norb*norb

    call zeroArray(mxsize,fock1)
    call zeroArray(mxsize,fock2)
    call zeroArray(mxsize,rhot)
    call zeroArray(mxsize,wk12)    
    ! (jk)
    call zeroArray(mxsize,excp(length3-mxsize))

    if (nel>1) then
       do iorb1=1,norb
          if (inhyd(iorb1).eq.1) cycle
          iborb1 =i1b (iorb1)
          call prodas (mxsize,occ(iorb1),psi(iborb1),psi(iborb1),rhot)
       enddo

       ! choose a given functional from libxc library     
       do func_id=1,lxcFuncs
          call xc_f90_func_init(xc_func, xc_info, lxcFuncs2use(func_id), XC_UNPOLARIZED)

          select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
             call xc_f90_lda_vxc(xc_func, mxsize, rhot(1), wk12(1))
          case(XC_FAMILY_GGA)
             ! calculate nabla rho nabla rho 
             call nfng (rhot,rhot,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
             call xc_f90_gga_vxc(xc_func, mxsize, rhot(1), wk10(1), wk12(1), wk3(1) )

             ! v_xc=vrho - 2 (nabla vsignma  nabla n + vsigma nabla^2 rho)
             
             ! vrho==wk12
             ! vsigma==wk3 
             ! wk11==(nabla vsignma  nabla n)
             call nfng (wk3,rhot,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
             ! wk8==nabla^2 rho
             call n2f (rhot,wk5,wk6,wk7,wk8)
             ! wk8==vsigma nabla^2 rho
             call prod (mxsize,wk3,wk8)
             ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 rho
             call add (mxsize,wk11,wk8)
             ! wk12==vrho - 2 (nabla vsignma  nabla n) + vsigma nabla^2 rho
             call axpy (mxsize,-two,wk8,ione,wk12,ione)
          case(XC_FAMILY_HYB_GGA)
             call nfng (rhot,rhot,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
             call xc_f90_gga_vxc(xc_func, mxsize, rhot(1), wk10(1), wk12(1), wk3(1) )
          case default
             write(*,'("Error! Undefined libxc functional.")')
             stop 'fockLXC'
          end select
          call xc_f90_func_end(xc_func)
          
          call add (mxsize,wk12,fock1)
       enddo
       ! take care of F4 factor
       call multf4(fock1)
    endif

    ! add local exchange and correlations contributions ones to one-electron contribution

    call zeroArray(mxsize,wk3)
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
          call axpy (ngpot1,occ(iorb1),pot(ibpot1),ione,wk3,ione)
          
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
       call add (mxsize,wk3,fock1)
       
       
       ! multiply coulomb/exchange potentials and off-diagonal Lagrange
       ! multipliers by f2
       call prod (mxsize,f2,fock1)
       call prod (mxsize,f2,fock2)
       ! nel.gt.1
    endif

    call copy (mxsize,f0,ione,wk3,ione)
    call axpy (mxsize,eng(iorb),f1,ione,wk3,ione)
    
    if (mm(iorb).ne.0) then
       !      e enter the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call axpy (mxsize,w,e,ione,wk3,ione)
    endif
    
    ! add Coulomb contributions to one-electron one (containing local
    ! exchange and correlation contributions)
    
    call add (mxsize,wk3,fock1)
    
  end subroutine fockLXC


end module fockLXC_m
