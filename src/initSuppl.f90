! ### initSuppl ###

!     Initialises various supplementary arrays of case-dependent lengths
!     supported by cw_suppl (one-electron potentials, Jacobians,
!     integration and differentiation weights, etc)


subroutine initSuppl (borb,bpot,d,e,f0,f1,f2,f3,f4,g,wjac1,wjac2,wgt1,wgt2)
  use params
  use discret
  use scf
  use commons8
  use kh, only: poth3, potkh

  use plegendg_m
  
  implicit none

  integer :: i,ib,ie,ig,ii,imu,in,izz1,izz2,j,k

  real (PREC) :: atw1,atw2,costh,precis1,r1,rr,rr2,w1,w2,w4,w5,w6,w7,w8,w9,w10,wj1,wj2,wmup,xrr,wktmp,&
       xxplusyy,xxplusyy2,xy2,z,z1t,z2t,zcm

  real (PREC), dimension(nni,mxnmu) :: borb,bpot,d,e,f0,f1,f2,f3,f4,g,wjac1,wjac2
  real (PREC), dimension(mxsize) :: wgt1,wgt2

  real (PREC), dimension(9) :: aa1,aa2,a1,a2
  real (PREC), dimension(9,7,9) :: dc1,dc2

  real (PREC), external :: zz1,zz2,zz1g,zz2g,zgsz1,zgsz2,zgsz1g,zgsz2g

!     Coefficients of the first and second derivatives taken from from
!     the 8th-order Stirling interpolation formula

  data aa1/ 3.0_PREC, -32.0_PREC, 168.0_PREC, -672.0_PREC, 0.0_PREC, 672.0_PREC, -168.0_PREC, 32.0_PREC, -3.0_PREC /
  data aa2/ -9.0_PREC,  128.0_PREC, -1008.0_PREC, 8064.0_PREC, -14350.0_PREC, 8064.0_PREC,-1008.0_PREC, 128.0_PREC, -9.0_PREC /

  precis1=precis

!     iaddr4(1)=borb
!     iaddr4(2)=bpot
!     iaddr4(3)=d
!     iaddr4(4)=e
!     iaddr4(5)=f0
!     iaddr4(6)=f1  R^2/2 (cosh^2(mu) - cos^2(nu))=R^2/2 (xi^2 -eta^2)
!     iaddr4(7)=f2
!     iaddr4(8)=f3
!     iaddr4(9)=f4  R/2 cosh(mu)
!     iaddr4(10)=g
!     iaddr4(11)=wjac1
!     iaddr4(12)=wjac2
!     iaddr4(13)=wgt1
!     iaddr4(14)=wgt2

!     wjac1 - Jacobian for Rayleigh quotient
!     wjac2 - Jacobian for two-electron integrals

!     initialize arrays needed for the poisson equations

!     f1 - all potentials must be multiplied by this factor times -1

!     f4 - array of factors multiplying Coulomb and exchange
!	   potentials (to make the tilda counterparts), see eq.9 and 10

  if (ifefield.ne.0) then
     izz1=nint(z1)
     izz2=nint(z2)
     atw1=atweight(izz1)
     atw2=atweight(izz2)
     zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC
  else
     zcm=0.0_PREC
  endif

  wj1 =-r*pii/2.0_PREC
  wj2 = r*r*pii/2.0_PREC
  xrr=half*r
  ib=ibmu(1)

  do ig=1,ngrids
     ie=iemu(ig)
     diag(ig)=-14350.0_PREC/(5040.0_PREC)*( 1.0_PREC/(hmu(ig)*hmu(ig))+1.0_PREC/(hni*hni) )


     do i=ib,ie
        w1=vxi1(i)
        w2=vxi(i)*vxi(i)
        w9=xrr*vxi(i)

        ! ssf is included in differentiation array dmu1
        ! w4=vxi2(i)/ssf

        w4=vxi2(i)
        if(i.eq.1) then
           w5=0.0_PREC
           w6=0.0_PREC
        else
           w5= vxi(i)/vxi1(i)-2.0_PREC*vxi1(i)/vxi(i)
           w6= 1.0_PREC/(vxi1(i)*vxi1(i))
        endif

        w8=-r*r*r*pii*vxi(i)/2.0_PREC

        do j=1,nni
        w10=veta1(j)
        borb(j,i)=w4
        bpot(j,i)=w5
           d   (j,i)=veta2(j)
           if (w10.lt.precis) then
              e(j,i)= 0.0_PREC
           else
              e(j,i)=-w6 - 1.0_PREC/(w10*w10)
           endif
           w7=w2-veta(j)*veta(j)
           f1(j,i)= r*r*w7/2.0_PREC

           if (ipot.eq.0) then
              !                 point nuclei - Coulomb potential
              f0(j,i)=r*(vxi(i)*(z1+z2)+veta(j)*(z2-z1))
           elseif (ipot.eq.1) then
              !                 finite nuclei - Fermi model
              z1t=zz1(i,j)
              z2t=zz2(i,j)
              f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
           elseif (ipot.eq.2) then
              !                 finite nuclei - Gauss model
              z1t=zz1g(i,j)
              z2t=zz2g(i,j)
!              print *,i,j,z1t,z2t
              f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
           elseif (ipot.eq.3) then
              !                 OED - 3D Coulomb potential reduced to 2D external
              !                 potential must be multiplied by -f1 and added to f0
              z=(r/2.0_PREC)*vxi(i)*veta(j)
              xxplusyy2=(r/2.0_PREC)*vxi1(i)*veta1(j)
              f0(j,i)=-f1(j,i)*poth3(z,xxplusyy2,mpot,v0pot,apot,precis1)
              !                  print *,i,j,z,xxplusyy2,apot,mpot,f0(j,i)
           elseif (ipot.eq.4) then
              ! OED - Kramers-Henneberger potential
              z=(r/2.0_PREC)*vxi(i)*veta(j)
              xxplusyy2=(r/2.0_PREC)*vxi1(i)*veta1(j)
              f0(j,i)=-f1(j,i)*potkh(z,xxplusyy2,mpot,epspot,ompot,v0pot,apot,nsimp,precis1)
           elseif (ipot.eq.5) then
              !                 OED - Green, Sellin, Zachor model potential
              z1t=zgsz1(i,j)
              z2t=zgsz2(i,j)
              f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
           elseif (ipot.eq.55) then
              !                 OED - Green, Sellin, Zachor model potential + finite Gauss nucleus model
              z1t=zgsz1g(i,j)
              z2t=zgsz2g(i,j)
              f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))

           elseif (ipot.eq.9) then
              !                 harmonic potential (Hook's atom, harmonium) located at A centre
              r1=(r/2.0_PREC)*(vxi(i)+veta(j))
              f0(j,i) =-hook/two*(r1*r1)*f1(j,i)

           elseif (ipot.eq.10) then
              !                 harmonic potential (x^2+y^2+z^2)
              z=(r/2.0_PREC)*vxi(i)*veta(j)
              xy2=(r*r/4.0_PREC)*(vxisq(i)-one)*(one-vetasq(j))
              f0(j,i) =  -(xy2+z*z)*r*r*half*(vxisq(i)-vetasq(j))
           endif


           if (magfield.eq.1) then
              !                 external static magnetic field along z-axis for the
              !                 first orbital on the input list

              rr=(r/two)*sqrt(vxisq(i)+vetasq(j)-one)
              rr2=rr*rr
              xxplusyy=((r/2.0_PREC)*vxi1(i)*veta1(j))**2

              !                  z=(r/two)*vxi(i)*veta(j)
              !                  sintheta=(abs((one-z*z/rr2)))
              !                  if (rr.eq.zero) then
              !                     sintheta=one
              !                  endif


              f0(j,i)=f0(j,i)+(half*gammaf*dble(mgx(3,1))+one/eight*gammaf*gammaf*xxplusyy)*(-one)*f1(j,i)

              !     &                 *(-one)*r*r*half*(vxi(i)*vxi(i)-veta(j)*veta(j))
              !                  write(*,'(2i5,6d16.4)') j,i,z,rr,rr2*sintheta,xxplusyy
           endif

           if (ifefield.eq.1) then

              !                 when external electric field is present (z-zcm)*E
              !                 has to be added (times -f1)

              z=(r/2.0_PREC)*vxi(i)*veta(j)
              if (abs(z).lt.zcutoff) then
                 f0(j,i)=f0(j,i)-ffield*(z-zcm)*f1(j,i)
              else
                 f0(j,i)=f0(j,i)-ffield*(z-zcm)*exp(zcutoff-abs(z))*f1(j,i)
              endif

           elseif (ifefield.eq.2) then
              ! FIXME
              !                 if electric field gradient is nonzero the interaction
              !                 energy is -(Q_z/4)*(dE/dz)

              !                 if electric field gradient is nonzero the interaction
              !                 energy is  -R^2/3 P(2,0)*(dE/dz) (Greiner, Electrodynamics, p.91)

              imu=i
              in=j
              z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
              xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
              rr=sqrt(xxplusyy+z*z)
              if (abs(rr).lt.precis) then
                 costh=0.0_PREC
              else
                 costh=z/rr
              endif

              wktmp=rr**2*plegendg(itwo,izero,costh)
              z=(r/2.0_PREC)*vxi(i)*veta(j)
              if (abs(z).lt.zcutoff) then
                 f0(j,i)=f0(j,i)-fgrad*wktmp*f1(j,i)/3.0_PREC
              else
                 f0(j,i)=f0(j,i)-fgrad*wktmp*f1(j,i)/3.0_PREC*exp(zcutoff-abs(z))
              endif
           endif

           if (iharm2xy.eq.1) then
              !                 external three-dimensional harmonic oscilator in x-y-z coordinates
              !                 potential must be multiplied by -f1 and added to f0
              imu=i
              in=j
              f0(j,i)=f0(j,i)-f1(in,imu)*harm2xy*harm2xy/2.0_PREC*xrr*xrr*(vxi(i)*vxi(i)+veta(j)*veta(j)-1.0_PREC)

           elseif (iharm2xy.eq.2) then
              !                 external two-dimensional harmonic oscilator in x-y coordinates
              !                 potential must be multiplied by -f1 and added to f0
              imu=i
              in=j
              f0(j,i)=f0(j,i)-harm2xy*harm2xy/2.0_PREC*xrr*xrr*f1(in,imu)*vxi1(imu)*vxi1(imu)*veta1(in)*veta1(in)
           endif

           if(i.eq.1) then
              f2(j,i)=0.0_PREC
           else
              f2(j,i)=-r*w7/vxi(i)
           endif
           f3(j,i)=-2.0_PREC/w2
           g (j,i)= w8*w7
           f4(j,i)= w9

           !              Integration needs hmu*hni but the integration routine works
           !              for hni*hni. That is why ssf must be inserted here.

           !              wjac1(j,i)=wj1*w1*w5/ssf(ig)
           !              wjac2(j,i)=wj2*w1*w5*(w3+w6)/(ssf(ig)*w0)

           wjac1(j,i)=wj1*w1*w10
           wjac2(j,i)=wj2*w1*w10*(w1*w1+w10*w10)/vxi(i)
        enddo
     enddo
     ib=ie+1
     !        enddo over ig
  enddo

  !     initialize vector of integration weights
  !     (7-point formula	Abramovic 25.4.16.)

  wmup=0.0_PREC
  do ig=1,ngrids
     ib=ibmu(ig)
     ie=iemu(ig)
     w1 =hmu(ig)/140.0_PREC
     wmu(ib)=wmup+w1*41.0_PREC
     do i=ib+1,ie,6
        wmu(i  )=w1*216.0_PREC
        wmu(i+1)=w1* 27.0_PREC
        wmu(i+2)=w1*272.0_PREC
        wmu(i+3)=w1* 27.0_PREC
        wmu(i+4)=w1*216.0_PREC
        if (i+5.eq.ie) then
           wmu(i+5)=w1*41.0_PREC
        else
           wmu(i+5)=2.0_PREC*w1*41.0_PREC
        endif
     enddo
     wmup=wmu(ie)
  enddo

  w1    =hni/140.0_PREC
  wni(1)=w1*41.0_PREC
  do i=2,nni,6
     wni(i  )=w1*216.0_PREC
     wni(i+1)=w1* 27.0_PREC
     wni(i+2)=w1*272.0_PREC
     wni(i+3)=w1* 27.0_PREC
     wni(i+4)=w1*216.0_PREC
     if (i+5.eq.nni) then
        wni(i+5)=w1*41.0_PREC
     else
        wni(i+5)=2.0_PREC*w1*41.0_PREC
     endif
     !        enddo over i
  enddo

  !     multiply Jacobians for one and two-electron integrals by elements
  !     of wmu array

  do i=1,mxnmu
     w1=wmu(i)
     do j=1,nni
        wjac1(j,i)=w1*wjac1(j,i)
        wjac2(j,i)=w1*wjac2(j,i)
     enddo
  enddo

  !     prepare superarrays with the integration weights over ni variable
  !     and with jacobians included

  do i=1,mxnmu
     ii=(i-1)*nni
     do j=1,nni
        wgt1(ii+j)=wni(j)*wjac1(j,i)
        wgt2(ii+j)=wni(j)*wjac2(j,i)
     enddo
  enddo

  !     Initialize arrays needed for differentiation based on the 9-point
  !     interpolation formula

  !     Initialize dni1, dni2, dmu1 and dmu2

  w1=1.0_PREC/(840.0_PREC)
  w2=1.0_PREC/(5040.0_PREC)

  do ig=1,ngrids
     do  i=1,4
        dmu2(i,ig)=w2*aa2(i)/(hmu(ig)*hmu(ig))
        dmu1(i,ig)=w1*aa1(i)/hmu(ig)
     enddo
  enddo

  do i=1,4
     dni2(i)=w2*aa2(i)/(hni*hni)
     dni1(i)=w1*aa1(i)/hni
  enddo

  !     Initialize arrays needed for fast (matrix times matrix) evaluation
  !     of the first and second dervivative terms in the Laplasian (dmu, dnu) and
  !     the first derivatives for gradients (d1mu, d1nu)

  !     derivatives over mu variable
  ib=1
  do ig=1,ngrids
     do k=1,9
        a2(k)=w2*aa2(k)/(hmu(ig)*hmu(ig))
        a1(k)=w1*aa1(k)/hmu(ig)
     enddo

     ie=iemu(ig)
     do i=ib,ie
        do k=1,9
           dmu(k,i)=a2(k)+borb(1,i)*a1(k)
           d1mu(k,i)=a1(k)
        enddo
     enddo
     ib=ie+1
  enddo

  !     derivatives over nu variable
  do k=1,9
     a2(k)=w2*aa2(k)/(hni*hni)
     a1(k)=w1*aa1(k)/hni
  enddo

  do i=1,nni
     do k=1,9
        dni(k,i)=a2(k)+d(i,1)*a1(k)
        d1ni(k,i)=a1(k)
     enddo
  enddo

  ! FIXME obsolete in version 2.0
  if (ngrids.ne.1) then

     !        initialize arrays dc1 and dc2

     !        in order to check lagrpol use two ajasent grids with equal step size
     !        and compare coefficients with those of the Stirling formulae

     call lagrpolq (dc1,dc2)

     !        modify dmu array at the grid boundary region
     !        derivative coeff. from 8th-order Lagrange interpolation formula
     !        are stored in dc2(ngbound,imu,k) and dc1(ngbound,imu,k),
     !        where ngbound is the number of
     !        the grid boundaries (1 for 1-2, 2 for 2-3 etc), imu is one of the 7
     !        intergrid mu values and k points to a derivative coefficient.

     do ig=1,ngrids-1
        ib=iemu(ig)-3
        ie=iemu(ig)+3
        imu=0
        do i=ib,ie
           imu=imu+1
           do k=1,9
              !                 a2(k)=dc2(ig,imu,k)/(hmu(ig)*hmu(ig))
              !                 a1(k)=dc1(ig,imu,k)/hmu(ig)
              a2(k)=dc2(ig,imu,k)
              a1(k)=dc1(ig,imu,k)
           enddo

           do k=1,9
              dmu(k,i)=a2(k)+borb(1,i)*a1(k)
           enddo
        enddo
     enddo
  endif

end subroutine initSuppl
