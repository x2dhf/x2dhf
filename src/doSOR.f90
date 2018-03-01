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
! ### doSOR ###

!     Controls relaxation of Coulomb potentials, exchange potentials
!     and orbitals

subroutine  doSOR (iorb,cw_sor,psi,pot,excp,borb,bpot,d,e,f0,f1,f2,f3,f4,g,wjac1,wjac2,wgt2,&
     lhs,rhs,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
  use params
  use commons8

  implicit none
  integer :: iorb
  integer, dimension(*) :: cw_sor
  real (PREC), dimension(*) ::  psi,pot,excp,wjac1,wjac2,borb,bpot,d,e, &
       f0,f1,f2,f3,f4,g,wgt2,lhs,rhs,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11

!   when ipoiss=3 then mcsor method is used to solve poisson's eqs for
!   potentials

  if (exlcoul.eq.0.0_PREC) then
     if (ipoiss.eq.1) then
        call coulSOR (iorb,cw_sor,psi,pot,excp,bpot,d,f3,g,lhs,rhs,wk3)
     elseif (ipoiss.eq.2.or.ipoiss.eq.3) then
        call coulMCSOR (iorb,cw_sor,psi,pot,bpot,d,f3,g,lhs,rhs,wk3)
     endif
  endif
  
  !   process the exchange potentials
  
  if (islat.eq.0) then
     if (exlexp.eq.0.0_PREC.or.exlexp.eq.2.0_PREC) then
        if (ipoiss.eq.1) then
           call exchSOR (iorb,cw_sor,psi,pot,excp,bpot,d,e,f3,g,lhs,rhs,wk3)
        elseif (ipoiss.eq.2.or.ipoiss.eq.3) then
           call exchMCSOR (iorb,cw_sor,psi,pot,excp,bpot,d,e,f3,g,lhs,rhs,wk3)
        endif
     endif
  endif
  
  if (islat.eq.1.and.iscmc.eq.1) then
     if (ipoiss.eq.1) then
        call exchSOR (iorb,cw_sor,psi,pot,excp,bpot,d,e,f3,g,lhs,rhs,wk3)
     elseif (ipoiss.eq.2.or.ipoiss.eq.3) then
        call exchMCSOR (iorb,cw_sor,psi,pot,excp,bpot,d,e,f3,g,lhs,rhs,wk3)
     endif
  endif
  
  if (exlorb.eq.0.0_PREC) then
     if (iprint(100).ne.0) call tail(psi)
     
     !       if (ifix(iorb).eq.1) return
     
     if (ipoiss.eq.1.or.ipoiss.eq.3) then
        if     (ioo.eq.1) then
           call orbSOR(iorb,cw_sor,psi,pot,excp,borb,d,e,f0,f1,f2,f4,wgt2,lhs,rhs,wk0,wk1,wk2,wk3,wk4, &
                wk5,wk6,wk7,wk8,wk9,wk10,wk11)
        ! elseif (ioo.eq.-1) then
        !    call orbSORrev(iorb,cw_sor,psi,pot,excp,borb,d,e,f0,f1,f2,f4,wgt2,lhs,rhs,wk0,wk1,wk2,wk3,wk4, &
        !         wk5,wk6,wk7,wk8,wk9,wk10,wk11)
        else
           write(*,*) 'Error: wrong ordering of mesh points'
           stop 'doSOR'
        endif
        
     elseif (ipoiss.eq.2) then
        call orbMCSOR(iorb,cw_sor,psi,pot,excp,borb,d,e,f0,f1,f2,f4,wgt2,lhs,rhs,wk0,wk1,wk2,wk3,wk4, &
             wk5,wk6,wk7,wk8,wk9,wk10,wk11)
     endif
  endif
  
end subroutine doSOR
