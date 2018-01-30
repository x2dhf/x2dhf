
!     optimal omega for 2nd-order stencil
function setOmegaPot2()
  use params
  use discret
  use commons8

  implicit none
  real (PREC) :: setOmegaPot2
  real (PREC) :: a,b,rho

  parameter (a=1.0_PREC,b=0.0_PREC)

  rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
  setOmegaPot2=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b

end function setOmegaPot2

