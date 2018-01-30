! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### tail ####
!
!     examines the tails of all orbitals and determines the location of
!     the last maximum in mu variable, i.e.
!     f(max)=max{vmax(1),vmax(2),...,vmax(nni)}
!     and the maximum value of psi at the last mi=nmi grid line

 subroutine tail (psi)
  use params
  use discret
  use commons8

  implicit none
  integer :: ibeg,imaxg,inimax,ijk,iorb
  real (PREC) :: vpmax,vp
  real (PREC), dimension(*) :: psi

  do iorb=1,norb
     ibeg=i1b(iorb)-1
     vpmax=-1.e30_PREC
     imaxg =0
     
     vpmax=-1.e30_PREC
     do ini=1,nni
        ijk=(mxnmu-1)*nni+ini
        vp=abs(psi(ibeg+ijk))
        if (vp.gt.vpmax) then
           vpmax=vp
           inimax=ini
        endif
     enddo
     write(*,5110) 'tail: ',iorb,inimax,vpmax
  enddo

5000 continue
  
5100 format('   orbital     ni    max{psi(ni,mxnmui)}    ')
5110 format(1x,a10,2i5,2x,e12.2)

end subroutine tail
