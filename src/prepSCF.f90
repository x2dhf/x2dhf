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
! ### prepSCF ####

!     Prepares an SCF process (orthogonalizes orbitals, calculates
!     orbital energies, etc)

module prepSCF_m
  implicit none
contains
  subroutine prepSCF (cw_sor,cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch)
    use params
    use discret
    use commons8

    use checkortho_m
    use ea_m
    use eaDFT_m
    use eab_m
    use eabDFT_m
    use etotal_m
    use etotalDFT_m
    use norm_m
    use ortho_m
    use prttoten_m
    use rfdexch_m
    use mpolemom_m

    implicit none
    integer :: i,iorb,iprtlev0,j,jorb

    integer, dimension(*) :: cw_sor
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch

    !     orbital energies and Lagrange multiplies are initialized
    !     or have been retrieved from the disk file

    !     asymptotic values of Coulomb and exchange potentials are
    !     calculated or have been retrieved from the disk file

    !     see routine initOrbPot for the explanation of the following command

       ! When flipping of orbitals is on then flip (i,j) pair if nfliporb(i,j)=1 or
       ! abs(denmax) is less then fliporbthresh

       if (ifliporb.eq.1) then
          do i=norb,1,-1
             do j=i-1,1,-1
                if (mgx(6,i).ne.mgx(6,j)) cycle
                if (eng(i).lt.eng(j).and.nfliporb(i,j).eq.0) cycle
                if (ige(i).ne.ige(j)) cycle
                write(iout6,16110) iorn(i),bond(i),gut(i),iorn(j),bond(j),gut(j)
                call fliporb(i,j,cw_orb,cw_coul,cw_sctch(i5b(1)))
             enddo
          enddo
       endif

    !     check orhogonality of orbitals
    if (iprint(30).ne.0) then
       print *,''
       print *,'orthogonality of orbitals:'
       iprint(20)=1
       do iorb=1,norb
          jorb=norb+1-iorb
          call checkOrtho (jorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
       enddo
       iprint(20)=0
    endif


    if (imethod.eq.3.or.imethod.eq.4.or.imethod.eq.5) islat=1
    if (imethod.eq.2.or.ini.eq.4) nel=1

    if (idump.ne.0) then
       do iorb=1,norb
          jorb=norb+1-iorb
          call norm (jorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
          call ortho (jorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
       enddo

       do i=1,norb
          itouch(i)=1
       enddo

       !        initialize multipole moments
       if (imethod.ne.2) then
          iprtlev0=iprtlev
          iprtlev=3
          call mpoleMom(cw_orb,cw_suppl,cw_sctch)
          write(*,*) '... initializing multipole moment coefficients ...'
          iprtlev=iprtlev0
       endif

       !        calculate diagonal and off-diagonal Lagrange multipliers
       do iorb=norb,1,-1
          if (iform.eq.0.or.iform.eq.2) call rfdexch(iorb,cw_exch)
          if (islat.eq.0) then
             call Ea (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),     &
                  cw_suppl(i4b(13)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)), &
                  cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
             engi(iorb)=eng(iorb)

             call Eab (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),    &
                  cw_suppl(i4b(13)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)), &
                  cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
          else
             call EaDFT (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),  &
                  cw_suppl(i4b(13)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)), &
                  cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
             engi(iorb)=eng(iorb)

             call EabDFT (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)), &
                  cw_suppl(i4b(13)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)), &
                  cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
          endif
       enddo

       write(*,*) '... initializing Lagrange multipliers ...'
    endif

    write(iout6,*)

    if (islat.eq.0) then
       call etotal (cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),      &
            cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)), &
            cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)), &
            cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),cw_sctch(i5b(13)),cw_sctch(i5b(14)))
    else
       call etotalDFT (cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),   &
            cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)), &
            cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)), &
            cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),cw_sctch(i5b(13)),cw_sctch(i5b(14)))
    endif
    call prttoten

16110 format(" ... flipping orbitals: ",i4,1x,a8,a1,3x,i4,1x,a8,a1)
    
  end subroutine prepSCF
end module prepSCF_m
