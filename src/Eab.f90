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
! c ### Eab ###
! c
! c     Calculates the off-diagonal Lagrange multipliers.

module Eab_m
  implicit none
contains
  subroutine Eab (iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
    use params
    use scf
    use commons8
    use Eab1HF_m

    implicit none

    integer :: iorb,iorb1,iorbbeg,ipc12,ipc21
    real (PREC) :: engo12,engo21,engoprv12,engoprv21,ent,wocc

    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3
    real (PREC), external :: dot

    if (nel.eq.1)  return
    if (norb.le.1) return
    if (iorb.eq.norb) return

    iorbbeg=iorb+1

    do iorb1=iorbbeg,norb
       !        if break is on and the two orbitals have different symmetry
       !        off-diagonal lm is not calculated
       !        if (ihomon.eq.2.and.ihomo(iorb1)*ihomo(iorb).lt.0) goto 10

       if ( (ilagra==0) .and. (ihomo(iorb1)*ihomo(iorb)<0 ) ) goto 10

       if (ilagra==0.and.nlm(iorb1,iorb)==0) goto 10

       if (ilagra==1.and.nlmf(iorb1,iorb)==0) goto 10

       if (nonorthog(iorb,iorb1)==1) go to 10

       ipc12=iorb1+(iorb-1)*norb
       ipc21=iorb+(iorb1-1)*norb

       !        store previous values in order to make damping of LM possible
       engoprv12=engo(ipc12)
       engoprv21=engo(ipc21)

       if (lmtype==0) then
          !           best convergence (100 SCF iterations for Li, grid=91/35)
          wocc=(occ(iorb1)+occ(iorb))/(occ(iorb1)*occ(iorb))
          call Eab1HF (iorb1,iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
          call Eab1HF (iorb,iorb1,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
          ent=(engo(ipc12)+engo(ipc21))/wocc
          engo12=ent/occ(iorb)
          engo21=ent/occ(iorb1)
       elseif (lmtype==1) then
          !           second best convergence (123 SCF iterations for Li)
          call Eab1HF (iorb1,iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
          engo12=engo(ipc12)
          engo21=engo(ipc12)/occ(iorb1)
       elseif (lmtype.eq.2) then
          !           worst convergence (1000 SCF iterations for Li)
          call Eab1HF (iorb,iorb1,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
          engo21=engo(ipc21)
          engo12=engo(ipc21)*occ(iorb1)
       else
          write(iout6,1000)
1000      format(/1x,'... off-diagonal Lagrange multiplies cannot be calculated ...'//)
          stop 'Eab'
       endif

       !        dflagra=0, sflagra=1 (obsolete)
       !        engo(ipc12)=sflagra*((1.0_PREC-dflagra)*engo12+dflagra*engoprv12)
       !        engo(ipc21)=sflagra*((1.0_PREC-dflagra)*engo21+dflagra*engoprv21)

       engo(ipc12)=engo12
       engo(ipc21)=engo21

       if (iprint(46).ne.0) then
          write(*,'(a8,i4,2e16.6,i4,a8,i4,a8,2e16.8)') 'Eab: ',lmtype,sflagra,dflagra, &
               iorn(iorb),bond(iorb1),iorn(iorb1),bond(iorb1),engo(ipc12),engo(ipc21)
       endif
10     continue
    enddo

  end subroutine Eab
end module Eab_m
