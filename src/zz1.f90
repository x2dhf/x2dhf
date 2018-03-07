! ### zz1 ###
!
!     Evaluate the nuclear potential for Fermi models.
!
!     This routine is based on nucpot routine from the GRASP2 package
!     (ver. 1.00, 1992)
!
function zz1(i,j)
  use params
  use discret
  use commons8
  use es_m

  implicit none

  integer :: i,j
  real (PREC) :: zz1
  real (PREC) :: a,abc,abc2,abc3,at3,atw,c,cba,dmsas,en,facto1,fmtoau,h3,h3php,hpiac2,pi2,&
       rbc,ri,rmc,rmcba,rrms,rrmsfm,s2mcba,s3mcba,s2rcba,s3rcba,sabc3,t,tabc,tfm,thabc2,zbn

  if (z1.eq.0.0_PREC) then
     zz1=0.0_PREC
     return
  endif

  fmtoau=1.0e-13_PREC/ainfcm

  !    potential energy is V(i,j)=-zz1(i,j)/r1 -z2/r2
  !     set atomic weight

  atw=z1atmass
  at3=atw**(one/three)
  rrmsfm = 0.8360_PREC*at3+0.5700_PREC
  tfm = 2.300_PREC

  !     change units from fm into bohr
  rrms = rrmsfm*fmtoau
  t = tfm*fmtoau
  a = t/(4.00_PREC*log(3.00_PREC))
  facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
  c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)

  abc = a/c
  tabc = two*abc
  abc2 = abc*abc
  thabc2 = three*abc2
  abc3 = abc2*abc
  cba = c/a
  pi2 = pii*pii
  hpiac2 = half*pi2*abc2
  six = two*three
  h3 = half*three
  h3php = h3+hpiac2
  call es (-cba,s2mcba,s3mcba)
  sabc3 = six*abc3
  dmsas = -sabc3*s3mcba
  en = one + abc2*pi2 + dmsas
  zbn = z1/en

  ri = r*(vxi(i)+veta(j))/2.0_PREC
  rmc = ri-c
  rmcba = rmc/a
  rbc = ri/c
  if (rbc .le. one) then
     call es (rmcba,s2rcba,s3rcba)
     zz1 = zbn*( dmsas + sabc3*s3rcba+rbc*( h3php-thabc2*s2rcba-half*rbc*rbc) )
  else
     call es (-rmcba,s2rcba,s3rcba)
     zz1 = z1 * ( one+thabc2 * ( rbc *s2rcba+tabc*s3rcba ) / en )
  endif

  if (iprint(131).ne.0) then
     if (abs(zz1-z1).gt.1.0e-8_PREC) then
        write(*,*) i,j,ri,zz1
     endif
  endif
end function zz1
