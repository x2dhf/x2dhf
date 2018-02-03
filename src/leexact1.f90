! ### leexact1 ###
!
subroutine leexact1 (vt,orb)
  use params
  use discret
  use commons8

  implicit none
  integer :: idev2,imu,iorb,ipb,ipoints
  real (PREC) :: c1,d1,ddnom,dev1,dev2,dnum,e,r1,r1max,s1,s2,s3,t,w
  real (PREC), dimension(nni,mxnmu) :: orb,vt

  iorb=1
  w=eng(iorb)

  dev1=zero
  dev2=zero
  idev2=0
  r1max=0.20_PREC
  ipoints=0
  izero=0

  do ini=1,nni
     do imu=1,mxnmu
        r1=(r/2.0_PREC)*(vxi(imu)+veta(ini))
        dnum=zero
        ddnom=zero
        s1=zero
        s2=zero
        s3=zero
        if (abs(r1).gt.precis) then
           do ipb=npbasis,1,-1
              d1=primexp(ipb)
              c1=primcoef(iorb,ipb)
              e=(2.0_PREC*d1/pii)**(3.0_PREC/4.0_PREC)*exp(-d1*r1*r1)
              s1=s1+c1*(-2.0_PREC*d1*d1*r1*r1)*e
              s2=s2+c1*(3.0_PREC*d1)*e
              if (abs(r1).gt.precis) s3=s3+c1*(-z1/r1)*e
              ddnom=ddnom+c1*e
           enddo
           dnum=s1+s2+s3
           if (idbg(497).ne.0.and.abs(ddnom).gt.precis) then
              t=(vt(ini,imu)-dnum/ddnom*orb(ini,imu))**2
              dev1=dev1+(vt(ini,imu)-dnum/ddnom*orb(ini,imu))**2
           else
              t=0.0_PREC
           endif
           if (abs(r1).gt.r1max) then
              dev2=dev2+t
              idev2=idev2+1
           endif
        endif
     enddo
  enddo
  write(*,*) 'deviation from exact local energy',sqrt(dev1)
  write(*,*) 'deviation from exact local energy outside radius'
  write(*,*) idev2,r1max,sqrt(dev2)

end subroutine leexact1
