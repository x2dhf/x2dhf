! ### nmucalc ###

!     Calculates nmu for given nmu and practical infinity values (a
!     single grid case is assumed) and adjusts nmu to be of the form
!     6k+1 and 5k+5+1 (if mcsor is chosen)

integer function nmucalc(n)
  use params
  use discret
  use scf
  use commons8


  implicit none

  integer :: ig,k5,k6,n,nk,ntemp


  real (PREC) :: hnu,xmi0
  parameter (k5=5, k6=6)
  
  !     Since in this version of the program the 7 point integration
  !     formula is used the number of mesh points in the nu and mu
  !     variables must be of the form 6*n+1
  
  !     Another restriction has to be imposed on the number of points in the
  !     mu variable if the mcsor scheme is to be used: nmu must be
  !     of the form kN+k+1, where k=iorder/2+1 for each grid (in this program
  !     iorder, i.e. the order of descritization is 8) 
  
  !     mxnmu --  maximum no. of grid points in mu variable
  !     hni - step in ni variable
  !     hmu - step in mu variable
  
  ig=1
  
  if (n.eq.0) then
     
     !        determine step size in nmu variable so as hmu is approximately
     !        equal to hni
     
     hnu=pii/dble(nni-1)
     
     xmi0=2.0_PREC*rgrid(ig)/r
     xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
     !        hmu(1)=xmi0/dble(mxnmu-1)
     mxnmu=xmi0/hnu+1
     nmu(ig)=mxnmu
     n=nmu(1)
  endif
  
  !     adjust, if necessary, the number of grid points in mu variable
  
12 nk=(n-1)/k6
  if (k6*nk.ne.n-1) then
     n=n-1
     goto 12
  endif
  
  ntemp=n
  if (ipoiss.eq.2.and.iorder(ngrids).eq.3) then
22   nk=(n-6)/k5
     if (k5*nk.ne.n-k5-1) then
        n=n-1
        goto 22
     endif
  endif
  
  if (n.gt.1.and.n.ne.ntemp) goto 12
  
  nmu(ig)=n
  mxnmu=nmu(ig)
  
  !     calculate hmu(1) for the adjusted value of mxnmu
  xmi0=2.0_PREC*rgrid(ig)/r
  xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
  hmu(1)=xmi0/dble(mxnmu-1)
  
  nmucalc=mxnmu
  
end function nmucalc
