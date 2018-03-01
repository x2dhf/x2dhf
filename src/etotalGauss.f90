! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1992 F.A. Parpia and I.P. Grant                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### etotalGauss ###
!
!     Evaluate the interaction energy between two nuclei with Gaussian
!     charge distributions Z1*Z2/r lowerGamma(half,eta12*r*r)/sqrt(pi),
!     1/eta12=1/eta1+1/eta2

subroutine etotalGauss
  use params
  use discret
  use commons8

  implicit none

  real (PREC) :: at3,atw,eta1,eta12,eta2,fmtoau,gammaLower,rrms,rrmsfm
  real (PREC), external :: dgamit,dgamma

  fmtoau=1.0e-13_PREC/ainfcm

  if (z1.eq.0.0_PREC) then
     eta1=1.e35_PREC
  endif

  if (z2.eq.0.0_PREC) then
     eta2=1.e35_PREC
  endif

  atw=z1atmass
  at3=atw**(one/three)
  rrmsfm = 0.8360_PREC*at3+0.5700_PREC
  rrms = rrmsfm*fmtoau
  eta1=3.0_PREC/(2.0_PREC*rrms*rrms)

  atw=z2atmass
  at3=atw**(one/three)
  rrmsfm = 0.8360_PREC*at3+0.5700_PREC
  rrms = rrmsfm*fmtoau
  eta2=3.0_PREC/(2.0_PREC*rrms*rrms)

  eta12=eta1*eta2/(eta1+eta2)
  !     lower incomplete gamma function
  gammaLower=(eta12*r*r)**(half)*dgamma(half)*dgamit(half,eta12*r*r)
  etotFN=etot-z1*z2/r+z1*z2*gammaLower/r/sqrt(pii)

end subroutine etotalGauss
