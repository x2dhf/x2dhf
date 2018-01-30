
!     gnuplot fit to omega_pot  (FH)
!        a = 0.568667         +/- 0.007234     (1.272%)
!        b = 0.863499         +/- 0.01428      (1.654%)

!     gnuplot fit to omega_pot  (LiH)
!        a = 0.52864          +/- 0.005446     (1.03%)
!        b = 0.942668         +/- 0.01075      (1.14%)

function setOmegaPot()
  use params
  use discret
  use commons8

  implicit none
  real (PREC) :: setOmegaPot
  real (PREC) :: a,b,rho

  parameter (a=0.5687_PREC,b=0.8635_PREC)
  !      parameter (a=0.5286,b=0.9427)

  rho=0.5*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
  setOmegaPot=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b
  
end function setOmegaPot

