! ### zz1g ###
!
!     Evaluate the nuclear potential for Gaussian models.
!
!     This routine is based on nucpot routine from the GRASP2 package (ver. 1.00, 1992)

function zz1g(i,j)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j
  real (PREC) :: zz1g
  real (PREC) :: at3,atw,eta1,fmtoau,ri,rrms,rrmsfm

  fmtoau=1.0e-13_PREC/ainfcm

  if (z1.eq.0.0_PREC) then
     zz1g=0.0_PREC
     eta1=1.d35
     return
  endif

  !     set atomic weight

  atw=z1atmass
  at3=atw**(one/three)
  rrmsfm = 0.8360_PREC*at3+0.5700_PREC

  !     change units from fm into bohr

  rrms = rrmsfm*fmtoau

  !     exponent of the Gaussian distribution $\rho_0 exp(-\eta r^2)$
  !     rrms=(3/(2\eta)^{1/2)

  eta1=3.0_PREC/(2.0_PREC*rrms*rrms)

  ri = r*(vxi(i)+veta(j))/2.0_PREC

  !     set the normalization constant $Z=rho_0 (pi/eta)^{3/2}$

  !     rho0=z1/(pii/eta1)**1.50_PREC
  !     x=eta1*ri*ri
  !     dg15= dgamit(1.50_PREC,x)*x**1.50_PREC*dgamma(1.50_PREC)
  !     zz1g = 2.00_PREC*pii*rho0*(dg15/eta1**1.50_PREC+ri*exp(-x)/eta1)

  !     making use of F77 intrinsic function
  !     erf(x)=2/(Pi)^(1/2)*int_)^x exp(-t^2)dt
  !     if available

  zz1g=z1*erf(sqrt(eta1)*ri)

end function zz1g
