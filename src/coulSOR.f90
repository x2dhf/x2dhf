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
! ### coulSOR ###

!     Evaluates the Coulomb potential for a given orbital

subroutine coulSOR (iorb,cw_sor,psi,pot,bpot,d,f3,g,lhs,rhs,wk2)
  use params
  use discret
  use scf
  use solver
  use commons8

  implicit none
  integer :: i,iborb,iorb,ibpot,ibpot1,ig,ioffs1,ioffst,itr1,itr2,ngrid

  integer, dimension(*) :: cw_sor
  real (PREC), dimension(*) :: psi,pot,rhs,bpot,d,f3,g,lhs,wk2

  if (ifix(iorb).eq.1) return
  
  if (nel.eq.1) return
  
!   prepare right-hand side of Poisson's equation

  iborb=i1b(iorb)
  ngrid=i1si(iorb)
  call prod2 (ngrid,psi(iborb),psi(iborb),wk2)
  call prod2 (ngrid,wk2,g,rhs)

  !   even symmetry for coulomb potential
  isym=1
  ibpot=i2b(iorb)
  
  do itr1=1,maxsor1
     do ig=1,i1ng(iorb)
        omega=ovfcoul(ig)
        omega1=1.0_PREC-omega
        ioffst=ioffs(ig)
        ioffs1=ioffst+1
        ibpot1=ibpot+ioffst

        do i=1,4
           dmu2t(i)=dmu2(i,ig)
           dmu1t(i)=dmu1(i,ig)
        enddo
        
        !         prepare left-hand side of the poisson equation
        !         include the diagonal part of the differentiation 
        !         operator in lhs
        
        do i=1,ngsize(ig)
           lhs(ioffst+i)=f3(ioffst+i)+diag(ig)
        enddo
        
        if (i1ng(iorb).eq.1) then 
           !            icase=1
           ifill=1
           ngrd1 =ingr1(1,ig)
           ngrd6a=ingr1(2,ig)
           ngrd6b=ingr1(3,ig)
           ngrd7 =ingr1(4,ig)
           
           call putin (nni,nmu(ig),isym,pot(iborb),wk2)
           do itr2=1,maxsorpot(iorb)
              call sor (wk2,lhs(ioffs1),rhs(ioffs1),bpot(ioffs1),d(ioffs1),cw_sor(iadext(ig)), &
                   cw_sor(iadnor(ig)),cw_sor(iadex1(ig)),cw_sor(iadex2(ig)),cw_sor(iadex3(ig)))
           enddo
           call putout (nni,nmu(ig),pot(iborb),wk2)
        else
           if (ig.eq.1) then
              !               icase=2
              ifill=2
              ngrd1 =ingr2(1,ig)
              ngrd6a=ingr2(2,ig)
              ngrd6b=ingr2(3,ig)
              ngrd7 =ingr2(4,ig)
              
              call putin2 (nni,nmu(ig),pot(iborb),wk2)
              
              do itr2=1,maxsorpot(iorb)
                 call sor (wk2,lhs(ioffs1),rhs(ioffs1),bpot(ioffs1),d(ioffs1),cw_sor(iadext(ig+ngrids)),cw_sor(iadnor(ig+ngrids)),&
                      cw_sor(iadex1(ig+ngrids)),cw_sor(iadex2(ig+ngrids)),cw_sor(iadex3(ig+ngrids)))
              enddo
              call putout (nni,nmu(ig),pot(iborb),wk2)     
           elseif (ig.ne.1.and.ig.ne.i1ng(iorb)) then 
              !               icase=2
              ifill=3
              ngrd1 =ingr2(1,ig)
              ngrd6a=ingr2(2,ig)
              ngrd6b=ingr2(3,ig)
              ngrd7 =ingr2(4,ig)

              muoffs=iemu(ig-1)-1
              call putin3 (nni,nmu(ig),pot(iborb),wk2)
              
              do itr2=1,maxsorpot(iorb)
                 call sor (wk2,lhs(ioffs1),rhs(ioffs1),bpot(ioffs1),d(ioffs1), &
                      cw_sor(iadext(ig+ngrids)),cw_sor(iadnor(ig+ngrids)),cw_sor(iadex1(ig+ngrids)), &
                      cw_sor(iadex2(ig+ngrids)),cw_sor(iadex3(ig+ngrids)))
              enddo
              call putout34 (nni,nmu(ig),pot(iborb),wk2)     
              muoffs=0
           elseif (ig.eq.i1ng(iorb)) then 
              ifill=4
              ngrd1 =ingr1(1,ig)
              ngrd6a=ingr1(2,ig)
              ngrd6b=ingr1(3,ig)
              ngrd7 =ingr1(4,ig)
              
              muoffs=iemu(ig-1)-1
              
              call putin4 (nni,nmu(ig),pot(iborb),wk2)
              
              do itr2=1,maxsorpot(iorb)
                 call sor (wk2,lhs(ioffs1),rhs(ioffs1),bpot(ioffs1),d(ioffs1), &
                      cw_sor(iadext(ig)),cw_sor(iadnor(ig)),cw_sor(iadex1(ig)),cw_sor(iadex2(ig)),cw_sor(iadex3(ig)))
              enddo

              call putout34 (nni,nmu(ig),pot(iborb),wk2)
              muoffs=0
           endif
        endif
     enddo
  enddo
  
end subroutine coulSOR
