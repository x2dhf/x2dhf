! ### zz1 ###
!
!     Evaluate the nuclear potential for Fermi models.  
!
!     This routine is based on nucpot routine from the GRASP2 package
!     (ver. 1.00, 1992)

function zz2(i,j)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j
  real (PREC) :: zz2
  real (PREC) :: a,abc,abc2,abc3,at3,atw,c,cba,dmsas,en,facto1,fmtoau,h3,h3php,hpiac2,pi2,&
       rbc,ri,rmc,rmcba,rrms,rrmsfm,s2mcba,s3mcba,s2rcba,s3rcba,sabc3,t,tabc,tfm,thabc2,zbn

  if (z2.eq.0.0_PREC) then
     zz2=0.0_PREC
     return
  endif

  fmtoau=1.0e-13_PREC/ainfcm
  
  !     Fermi distribution
  
  !     potential energy is V(i,j)=-zz1(i,j)/r1 -z2/r2
  !     set atomic weight for 
  
  atw=z2atmass
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
  zbn = z2/en
  
  ri = r*(vxi(i)-veta(j))/2.0_PREC
  rmc = ri-c
  rmcba = rmc/a
  rbc = ri/c
  if (rbc .le. one) then
     call es (rmcba,s2rcba,s3rcba)
     zz2 = zbn*( dmsas + sabc3*s3rcba+rbc*( h3php-thabc2*s2rcba-half*rbc*rbc) )
  else
     call es (-rmcba,s2rcba,s3rcba)
     zz2 = z2 * ( one+thabc2 * ( rbc *s2rcba+tabc*s3rcba ) / en )
  endif
  
  if (iprint(132).ne.0) then
     if (abs(zz2-z2).gt.1.0e-5_PREC) then
        write(*,*) i,j,ri,zz2
     endif
  endif
  
end function zz2
