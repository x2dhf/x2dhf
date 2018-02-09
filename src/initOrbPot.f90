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
! ### initOrbPot ###
!
!     Initializes orbitals and potentials.

subroutine initOrbPot (cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch)
  use params
  use discret
  use commons8
  use prepGauss
  implicit none

  integer :: iorb

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch

  character*20 caseinp,psiinp,coulinp,exchinp,caseout,psiout,coulout,exchout

  data caseinp,psiinp,coulinp,exchinp /'2dhf_input.dat', '2dhf_input.orb','2dhf_input.coul','2dhf_input.exch'/
  data caseout,psiout,coulout,exchout /'2dhf_output.dat', '2dhf_output.orb','2dhf_output.coul','2dhf_output.exch'/

  !        the program works on a set of standard input and output files
  !        containing orbitals, Coulomb potentials (with extentions orb
  !        and coul, respectively), exchange potentials (a file with
  !        extension exch or file fort.31, fort.32, ...) and a text file
  !        with data defining the case (*.dat)

  !        open separate files for reading and writing data defining a case
  open(iinp14,file=caseinp, status='unknown',form='formatted')
  open(iout24,file=caseout, status='unknown',form='formatted')

  !        open separate files for reading and writing orbitals and potentials

  if (iform.eq.0) then
     if (inpform.eq.0) then
        open(11,file=psiinp, status='unknown',form='unformatted')
        open(12,file=coulinp,status='unknown',form='unformatted')
        iinp13=iout23
     else
        open(11,file=psiinp, status='unknown',form='formatted')
        open(12,file=coulinp,status='unknown',form='formatted')
        iinp13=iout23
     endif

     if (ioutform.eq.0) then
        open(21,file=psiout, status='unknown',form='unformatted')
        open(22,file=coulout,status='unknown',form='unformatted')
     else
        open(21,file=psiout, status='unknown',form='formatted')
        open(22,file=coulout,status='unknown',form='formatted')
     endif
  elseif (iform.eq.1) then
     if (inpform.eq.0) then
        open(11,file=psiinp, status='unknown',form='unformatted')
        open(12,file=coulinp,status='unknown',form='unformatted')
        open(13,file=exchinp,status='unknown',form='unformatted')
        rewind(13)
     else
        open(11,file=psiinp, status='unknown',form='formatted')
        open(12,file=coulinp,status='unknown',form='formatted')
        open(13,file=exchinp,status='unknown',form='formatted')
        rewind(13)
     endif

     if (ioutform.eq.0) then
        open(21,file=psiout, status='unknown',form='unformatted')
        open(22,file=coulout,status='unknown',form='unformatted')
     else
        open(21,file=psiout, status='unknown',form='formatted')
        open(22,file=coulout,status='unknown',form='formatted')
     endif

  elseif (iform.eq.2) then
     if (inpform.eq.0) then
        open(11,file=psiinp, status='unknown',form='unformatted')
        open(12,file=coulinp,status='unknown',form='unformatted')
     else
        open(11,file=psiinp, status='unknown',form='formatted')
        open(12,file=coulinp,status='unknown',form='formatted')
     endif

     if (ioutform.eq.0) then
        open(21,file=psiout, status='unknown',form='unformatted')
        open(22,file=coulout,status='unknown',form='unformatted')
        open(23,file=exchout,status='unknown',form='unformatted')
        rewind(23)
     else
        open(21,file=psiout, status='unknown',form='formatted')
        open(22,file=coulout,status='unknown',form='formatted')
        open(23,file=exchout,status='unknown',form='formatted')
        rewind(23)
     endif

  elseif (iform.eq.3) then
     if (inpform.eq.0) then
        open(11,file=psiinp, status='unknown',form='unformatted')
        open(12,file=coulinp,status='unknown',form='unformatted')
        open(13,file=exchinp,status='unknown',form='unformatted')
        rewind(13)
     else
        open(11,file=psiinp, status='unknown',form='formatted')
        open(12,file=coulinp,status='unknown',form='formatted')
        open(13,file=exchinp,status='unknown',form='formatted')
        rewind(13)
     endif

     if (ioutform.eq.0) then
        open(21,file=psiout, status='unknown',form='unformatted')
        open(22,file=coulout,status='unknown',form='unformatted')
        open(23,file=exchout,status='unknown',form='unformatted')
        rewind(23)
     else
        open(21,file=psiout, status='unknown',form='formatted')
        open(22,file=coulout,status='unknown',form='formatted')
        open(23,file=exchout,status='unknown',form='formatted')
        rewind(23)
     endif
  endif

  rewind(11)
  rewind(12)
  rewind(21)
  rewind(22)


  if (ini.eq.1) then
     !        'hydrogen'

     !        to initialize properly functions and orbital energies in the case of
     !        the HFS calculations islat (which is nonzero) parameter has to be kept
     !        zero before the program reaches the doSCF routine.

     !         if (islat.ne.0) islat=0

     call initHyd (cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
     idump=1
     do iorb=1,norb
        inhyd(iorb)=inhydlcao(iorb)
     enddo

     !        islat=1
  elseif (ini.eq.11) then

     !        'qrhf
     call initHF (cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
     idump=1

  elseif (ini.eq.2) then

     !        'gauss'

     !  initial values of orbitals are provided by the GAUSSIAN program

     write(*,*) 'Initializing orbitals using GAUSSIAN output'
     call prepare_Gaussian
     call initGauss(cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
     idump=1

  elseif (ini.eq.3) then

     !        'gauss-c'

     stop 'Support for customized GAUSSIAN format has been deprecated. Use the standard format.'

  elseif (ini.eq.4) then

     !        'method: OED'

     idump=0
     call initDisk (cw_orb,cw_coul,cw_exch,cw_sctch(i5b(1)))

     !        To generate a higher lying one-electron diatomic state all the lower
     !        ones have to be retrieved from disk file and kept for the sake of
     !        the orthogonalization.

     !        The current highest orbital is being initialized as a linear
     !        combination of hydrogenic atomic orbitals unless the previous
     !        case is being continued

     if (ini4.ne.0) then
        call initHyd (cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
        idump=1
     endif

  elseif (ini.eq.5) then

     !        'old'

     !        When the incomplete orbital disk file is encountered and the method is OED the
     !        values of ini4 is set to norb-norb_p
     idump=0

     call initDisk (cw_orb,cw_coul,cw_exch,cw_sctch(i5b(1)))

     if (ini4.ne.0) then
        call initHyd (cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
        idump=1
     endif

  elseif (ini.eq.6) then

     !        'noexch'

     idump=0
     call initDisk (cw_orb,cw_coul,cw_exch,cw_sctch(i5b(1)))
     print *,'... initializing exchange potentials ...'
     call initPot(cw_orb,cw_coul,cw_exch,cw_suppl(i4b(7)),cw_suppl(i4b(9)),cw_sctch(i5b(1)))

  elseif (ini.eq.55) then

     !        'nodat'

     !        old format of input data is requested (skipping reading 2dhf_input.dat file)
     idat=1
     idump=0
     call initDisk (cw_orb,cw_coul,cw_exch,cw_sctch(i5b(1)))

  endif

end subroutine initOrbPot
