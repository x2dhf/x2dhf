! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### ehemsic ###
! ### SOURCE: exbe88.f and exbe88sup.f
!     Calculates SIC according to Samal and Harbola
!             (J.Phys.B:At.Mol.Opt.Phys. 38 (2005) 3765-3777
!					-- eqn. 19

module ehemsic_m
  implicit none
contains
  subroutine ehemsic (psi,pot,wgt2,wk1,wk0)
    use params
    use discret
    use commons8
    use blas_m

    implicit none
    integer :: i,iborb,ibpot,iorb,isiorb,iksic,isipot,nmut

    real (PREC) ::  const13,const43,cxlda,cxlsd,esic1,esic2,w

    real (PREC), dimension(*) ::  psi,pot,wgt2,wk0,wk1

    parameter (const13=1.0_PREC/3.0_PREC, const43=4.0_PREC/3.0_PREC)

    ! **** E^{SIC}_i[\phi_i] for each orbital
    do i=1,mxsize
       wk0(i) = 0.0_PREC
       wk1(i) = 0.0_PREC
    enddo
    write(*,'(" ** WARNING ** [SIC-(1)] MULTIPLY BY ocup = ")')
    ! contrib.f
    !     Coulomb potential contributions within the same shell
    esic1 = zero
    do iorb=1,norb
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       !         oc=occ(iorb)
       !        contribution from coulomb interaction within the same shell
       !        if (oc.gt.one) then
       if (isipot.ne.isiorb) then
          write(iout6,*) 'coulomb potentials and orbitals have to defined on the same number of subgrids'
          stop 'contrib'
       endif
       call copy (isipot,pot(ibpot),ione,wk0,ione)
       call prod  (isiorb,psi(iborb),wk0)
       call prod  (isiorb,psi(iborb),wk0)
       w=dot(isiorb,wgt2,ione,wk0,ione)
       !            esic1=w*oc/two
       esic1=w/two
       !            wdcoul=wdcoul+oc*(oc-one)/2.0_PREC*w
       !         endif
       write(*,'("     SIC-1            = ",i3,f20.12)') iorb,esic1

    enddo

    !-----------------
    ! E^{LSD}_x[\rho(\phi_i)] for each orbital
    CXLDA   = -(three/four)*(three/pii)**const13
    CXLSD   = CXLDA*two**const13

    do iksic=1,norb
       esic2 = zero
       do i=1,mxsize
          wk1(i)=0.0_PREC
       enddo

       do iorb=1,norb
          iborb=i1b(iorb)
          isiorb=i1si(iorb)
          nmut=i1mu(iorb)

          if(iorb.EQ.iksic)THEN
             call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
          endif
       enddo

       do i=1,mxsize
          wk1(i)=wk1(i)**const43
       enddo

       !     take care of f4 factor
       call multf4(wk1)
       esic2 = CXLSD*dot(mxsize,wgt2,ione,wk1,ione)

       write(*,'("     SIC-2            = ",i3,f20.12)') iksic,esic2
    enddo

    !-- ALTERNATE METHOD FOR \iint \frac{|\phi_i(r_1)|^2 * |\phi_i(r_1)|^2}
    !--                                 {|r_1 - r_2|} dr_1 dr_2
    !-- Same result as "wsic1"
    !--      do iorb=1,norb
    !--         ibpot=i2b(iorb)
    !--         isipot=i2si(iorb)
    !--         IF(iorb.EQ.iksic)THEN
    !--           oc=1.0_PREC
    !--         ELSE
    !--           oc=0.0_PREC
    !--        ENDIF
    !--         call axpy (isipot,oc,pot(ibpot),ione,wk1,ione)
    !--      enddo

    !     contribution from the Coulomb interaction
    !--      wndc = 0.0_PREC
    !--      do iorb=1,norb
    !--         iborb=i1b(iorb)
    !--         isiorb=i1si(iorb)
    !--         IF(iorb.EQ.iksic)THEN
    !--         call prod2 (isiorb,psi(iborb),psi(iborb),wk0)
    !--         call prod (isiorb,wk1,wk0)
    !--         w=dot(isiorb,wgt2,ione,wk0,ione)
    !--         wndc=wndc+w
    !--         ENDIF
    !--      enddo
    !--      esic1 = wndc
    !--      write(*,'("     Exchange [MLSD-SIC] (1)            = ",
    !--     &                                     f20.12)') esic1

  end subroutine ehemsic
end module ehemsic_m
