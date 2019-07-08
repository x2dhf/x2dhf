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
! ### orbMCSOR ###

!     Evaluates the Fock potential for a given orbital, sets up the
!     right-hand side of the Poisson equation for that orbital and
!     performs a few MCSOR iterations.

module orbMCSOR_m
  implicit none
contains
  subroutine orbMCSOR (iorb,cw_sor,psi,pot,excp,b,d,e,f0,f1,f2,f4,wgt2,  &
       lhs,rhs,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,fock1,fock2)
    use params
    use discret
    use scf
    use solver
    use commons8

    use fock_m
    use fockDFT_m
    use fockLXC_m
    use mcsor_m
    use orbAsymptDet_m
    use orbAsymptSet_m
    use putin_m
    use putin2_m
    use putin3_m
    use putin4_m
    use putout_m
    use putout34_m

    implicit none
    integer :: i,iborb,iborb1,ig,ii,ij,ioffs1,ioffst,iorb,itr1,itr2,j,nmut

    integer, dimension(*) :: cw_sor
    real (PREC), dimension(*) ::  psi,pot,excp,b,d,e,f0,f1,f2,f4,wgt2,lhs,rhs, &
         wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,fock1,fock2

    !     prepare right and left-hand side of the Poisson equation
    !     include the diagonal part of the differentiation operator in lhs

    if (islat.eq.0) then
       ! call fock (iorb,psi,pot,excp,e,f0,f1,f2,f4,wk0,wk1,wk2,wk3)
       call fock (iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,rhs,wk1,wk2)
    elseif (lxcFuncs>0) then
       ! use exchange-correlation functionals from libxc
       call fockLXC(iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,rhs, &
            wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,fock2,lhs)
    else
       ! use self-coded exchange-correlation functionals
       call fockDFT(iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,rhs, &
            wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,fock2,lhs)
    endif

    
    !     start the sor iterations
    !     fill wk2 array

    isym=isymOrb(iorb)
    iborb=i1b(iorb)

    do itr1=1,maxsor1
       do ig=1,i1ng(iorb)
          muoffs=0
          omega=ovforb(ig)
          omega1=1.0_PREC-omega
          ioffst=ioffs(ig)
          ioffs1=ioffst+1
          iborb1=iborb+ioffst

          do i=1,4
             dmu2t(i)=dmu2(i,ig)
             dmu1t(i)=dmu1(i,ig)
          enddo

          do i=1,ngsize(ig)
             lhs(ioffst+i)=fock1(ioffst+i)+diag(ig)
          enddo

          if (i1ng(iorb).eq.1) then
             !              icase=1
             ifill=1
             ngrd1 =ingr1(1,ig)
             ngrd6a=ingr1(2,ig)
             ngrd6b=ingr1(3,ig)
             ngrd7 =ingr1(4,ig)

             !              initialize aorb array needed to evaluate orbital values in
             !              the asymptotic region

             nmut=i1mu(iorb)
             do i=nmut-4,nmut
                ii=(i-1)*nni
                do j=1,nni
                   ij=ii+j
                   if (abs(f1(ij)).le.precis) then
                      wk1(ij)=0.0_PREC
                   else
                      wk1(ij)=lhs(ij)/f1(ij)
                   endif
                enddo
             enddo

             call orbAsymptDet (nmut,iorb,wk3,wk1)
             call putin (nni,nmu(ig),isym,psi(iborb),wk2)

             !               maxsor=maxsororb(iorb)
             do itr2=1,maxsororb(iorb)
                call mcsor(                                               &
                     wk2,lhs(ioffs1),rhs(ioffs1),b(ioffs1),d(ioffs1),     &
                     cw_sor(iadext(ig)),cw_sor(iadnor(ig)),               &
                     cw_sor(iadex1(ig)),cw_sor(iadex2(ig)),               &
                     cw_sor(iadex3(ig)))
             enddo
             call putout (nni,nmu(ig),psi(iborb),wk2)
             call orbAsymptSet (nmut,psi(iborb),wk3)
          else
             if (ig.eq.1) then
                !                 icase=2
                ifill=2
                ngrd1 =ingr2(1,ig)
                ngrd6a=ingr2(2,ig)
                ngrd6b=ingr2(3,ig)
                ngrd7 =ingr2(4,ig)

                call putin2 (nni,nmu(ig),psi(iborb),wk2)
                do itr2=1,maxsororb(iorb)
                   call mcsor                                          &
                        (wk2,lhs(ioffs1),rhs(ioffs1),b(ioffs1),        &
                        d(ioffs1),                                     &
                        cw_sor(iadext(ig+ngrids)),                     &
                        cw_sor(iadnor(ig+ngrids)),                     &
                        cw_sor(iadex1(ig+ngrids)),                     &
                        cw_sor(iadex2(ig+ngrids)),                     &
                        cw_sor(iadex3(ig+ngrids)))
                enddo
                call putout (nni,nmu(ig),psi(iborb),wk2)
             elseif (ig.ne.1.and.ig.ne.i1ng(iorb)) then
                !                icase=2
                ifill=3
                ngrd1 =ingr2(1,ig)
                ngrd6a=ingr2(2,ig)
                ngrd6b=ingr2(3,ig)
                ngrd7 =ingr2(4,ig)

                muoffs=iemu(ig-1)-1
                call putin3 (nni,nmu(ig),psi(iborb),wk2)
                do itr2=1,maxsororb(iorb)
                   call mcsor                                          &
                        (wk2,lhs(ioffs1),rhs(ioffs1),b(ioffs1),        &
                        d(ioffs1),                                     &
                        cw_sor(iadext(ig+ngrids)),                     &
                        cw_sor(iadnor(ig+ngrids)),                     &
                        cw_sor(iadex1(ig+ngrids)),                     &
                        cw_sor(iadex2(ig+ngrids)),                     &
                        cw_sor(iadex3(ig+ngrids)))
                enddo
                call putout34 (nni,nmu(ig),psi(iborb),wk2)
                muoffs=0
             elseif (ig.eq.i1ng(iorb)) then

                !                 initialize aorb array needed to evaluate orbital values in
                !                 the asymptotic region

                nmut=i1mu(iorb)
                do i=nmut-4,nmut
                   ii=(i-1)*nni
                   do j=1,nni
                      ij=ii+j
                      if (abs(f1(ij)).le.precis) then
                         wk1(ij)=0.0_PREC
                      else
                         wk1(ij)=lhs(ij)/f1(ij)
                      endif
                   enddo
                enddo
                !                  call iniaorb (nmut,wk0,wk1)

                ifill=4
                ngrd1 =ingr1(1,ig)
                ngrd6a=ingr1(2,ig)
                ngrd6b=ingr1(3,ig)
                ngrd7 =ingr1(4,ig)

                muoffs=iemu(ig-1)-1

                call putin4 (nni,nmu(ig),psi(iborb),wk2)
                do itr2=1,maxsororb(iorb)
                   call mcsor                                      &
                        (wk2,lhs(ioffs1),rhs(ioffs1),b(ioffs1),     &
                        d(ioffs1),                                  &
                        cw_sor(iadext(ig)),cw_sor(iadnor(ig)),      &
                        cw_sor(iadex1(ig)),cw_sor(iadex2(ig)),      &
                        cw_sor(iadex3(ig)))
                enddo
                call putout34 (nni,nmu(ig),psi(iborb),wk2)
                !                  call asymorb (nmut,psi(iborb),wk0)
                muoffs=0
             endif
          endif
       enddo
    enddo

  end subroutine orbMCSOR
end module orbMCSOR_m
