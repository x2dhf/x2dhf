! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### propet2 ####

!     Calculates dipole, quadrupole, octopole and hexadecapole moments
!     relative to the geometrical centre and the centre of mass.

subroutine propet2 (cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i1beg,igp,imu,in,inioff,iorb,izz1,izz2,k,mu,ni,ngorb

  real (PREC) :: atw1,atw2,cm2zz,costh,qq,qtot,qktot,qktot1,qkz1,qkz2,r1,rr,&
       sum1,sum2,sum3,xxplusyy,w,z,zcm

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5
  real (PREC), external :: dot,plegendg

  !     1 bohr electron = 2.541765 Debye -- Gaussian94 User's Reference
  !     data au2Debye /2.541765d0/
  !     data bohr2ang /0.529177249d0/


  do k=1,maxmpole
     elect(k)=0.0_PREC
     electDA(k)=0.0_PREC
     total(k)=0.0_PREC
     totalDA(k)=0.0_PREC
  enddo

  zcm=0.0_PREC
  qtot=0.0_PREC

  do iorb=1,norb
     qq=-occ(iorb)
     qtot=qtot+qq
     elect(1)=elect(1)+qq*cmulti(iorb)
     elect(2)=elect(2)+qq*cmulti(iorb+norb)
     if (mpole.ge.3) elect(3)=elect(3)+qq*cmulti(iorb+2*norb)
     if (mpole.ge.4) elect(4)=elect(4)+qq*cmulti(iorb+3*norb)
  enddo

  do k=1,mpole
     qkz1=z1*abs(r/2.0_PREC)**k*plegendg(k,izero,-1.0_PREC)
     qkz2=z2*abs(r/2.0_PREC)**k*plegendg(k,izero, 1.0_PREC)
     qktot1=elect(k)+qkz1+qkz2
     total(k)=qktot1
     totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
     electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)
  enddo

  if (iprint16.eq.0) then
     write(*,'(/," multipole moments relative to geometrical centre:"/)')

     write(*,'(36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif
  else
     write(*,'(/," multipole moments relative to geometrical centre:"/)')

     write(*,'(44x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif
  endif

  write(*,'(//," multipole moments relative to centre of mass:")')

  izz1=nint(z1)
  izz2=nint(z2)
  atw1=atweight(izz1)
  atw2=atweight(izz2)
  zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC
  write(iout6,1000) element(izz1),atweight(izz1),-r/2.0_PREC,element(izz2), atweight(izz2),r/2.0_PREC,zcm

01000 format(/,6x,'centre     atomic weight       z',/8x,a2,4x,f14.8,3x,f10.4,/, &
           8x,a2,4x,f14.8,3x,f10.4,/,2x,'centre-of-mass',13x,f18.14/)

  ! multipole moments relative to the centre of mass

  if (iprint(214).ne.0) then

     !        estimate round-off errors present in orbital contributions to
     !        multipole moments of a given order k

     write(iout6,'(4x,"k",5x,"orb",17x,"contribution",9x,"error" )')
  endif

  do k=1,mpole

!  loop over grid points

     do imu=1,mxnmu
        inioff=(imu-1)*nni
        do in=1,nni
           igp=inioff+in

           !              rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
           !              x=(r/2.0_PREC)*vxi1(imu)*veta1(in)

           z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
           xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
           rr=sqrt(xxplusyy+z*z)

           if (abs(rr).lt.precis) then
              costh=0.0_PREC
           else
              costh=z/rr
           endif
           wk0(igp)=rr**k*plegendg(k,izero,costh)
        enddo
     enddo

     qtot=0.0_PREC
     cm2zz=0.0_PREC
     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        qq=-occ(iorb)
        qtot=qtot+qq
        call prod2 (ngorb,cw_orb(i1beg),cw_orb(i1beg),wk3)

        nni2=nni/2
        sum1=0.0_PREC
        sum2=0.0_PREC
        do imu=mxnmu,1,-1
           inioff=(imu-1)*nni
           do in=1,nni2
              igp=inioff+in
              sum1=sum1+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo

           do in=nni2+1,nni
              igp=inioff+in
              sum2=sum2+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo
        enddo

        sum3=0.0_PREC
        do imu=1,mxnmu
           inioff=(imu-1)*nni
           do in=1,nni
              igp=inioff+in
              wk5(igp)=wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
              sum3=sum3+wk5(igp)
           enddo
        enddo

        if (iprint(214).ne.0) then
           write(iout6,'(i5,i4,1x,a8,a1,e28.16,e10.2)') k, iorn(iorb),bond(iorb),gut(iorb),sum3,sum1+sum2-sum3
        endif

        cm2zz=cm2zz+qq*(sum1+sum2)
     enddo

     elect(k)=cm2zz
     electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)

     qkz1=z1*abs(r/2.0_PREC+zcm)**k*plegendg(k,izero,-1.0_PREC)
     qkz2=z2*abs(r/2.0_PREC-zcm)**k*plegendg(k,izero,1.0_PREC)
     qktot=cm2zz+qkz1+qkz2
     total(k)=qktot
     totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
  enddo


  if (iprint16.eq.0) then
     write(*,'(/,36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif
  else
     write(*,'(/,44x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif
  endif

  if (iprint(216).ne.0) then
     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r_1^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              r1=(r/2.0_PREC)*(vxi(mu)+veta(ni))
              if (r1.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/r1
              endif
              wk1(imu+ni)=r1
              wk2(imu+ni)=r1*r1
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w
     enddo

     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)*r/2.0_PREC
              if (rr.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/rr
              endif
              wk1(imu+ni)=rr
              wk2(imu+ni)=rr*rr
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w

        call prod2 (ngorb,cw_orb(i1beg),cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <   > =',w
     enddo

     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r_2^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              r1=(r/2.0_PREC)*(vxi(mu)-veta(ni))
              if (r1.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/r1
              endif
              wk1(imu+ni)=r1
              wk2(imu+ni)=r1*r1
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w
     enddo
  endif

  return
end subroutine propet2


! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### propet2 ####

!     Calculates dipole, quadrupole, octopole and hexadecapole moments
!     relative to the geometrical centre and the centre of mass.

subroutine propet4 (cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i1beg,igp,imu,in,inioff,iorb,izz1,izz2,k,mu,ni,ngorb

  real (PREC) :: atw1,atw2,cm2zz,costh,qq,qtot,qktot,qktot1,qkz1,qkz2,r1,rr,&
       sum1,sum2,sum3,xxplusyy,w,z,zcm

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5
  real (PREC), external :: dot,plegendg

  !     1 bohr electron = 2.541765 Debye -- Gaussian94 User's Reference
  !     data au2Debye /2.541765d0/
  !     data bohr2ang /0.529177249d0/


  do k=1,maxmpole
     elect(k)=0.0_PREC
     electDA(k)=0.0_PREC
     total(k)=0.0_PREC
     totalDA(k)=0.0_PREC
  enddo

  zcm=0.0_PREC
  qtot=0.0_PREC

  do iorb=1,norb
     qq=-occ(iorb)
     qtot=qtot+qq
     elect(1)=elect(1)+qq*cmulti(iorb)
     elect(2)=elect(2)+qq*cmulti(iorb+norb)
     if (mpole.ge.3) elect(3)=elect(3)+qq*cmulti(iorb+2*norb)
     if (mpole.ge.4) elect(4)=elect(4)+qq*cmulti(iorb+3*norb)
  enddo

  do k=1,mpole
     qkz1=z1*abs(r/2.0_PREC)**k*plegendg(k,izero,-1.0_PREC)
     qkz2=z2*abs(r/2.0_PREC)**k*plegendg(k,izero, 1.0_PREC)
     qktot1=elect(k)+qkz1+qkz2
     total(k)=qktot1
     totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
     electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)
  enddo

  if (iprint16.eq.0) then
     write(*,'(/," multipole moments relative to geometrical centre:"/)')

     write(*,'(36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif
  else
     write(*,'(/," multipole moments relative to geometrical centre:"/)')

     write(*,'(44x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif
  endif

  write(*,'(//," multipole moments relative to centre of mass:")')

  izz1=nint(z1)
  izz2=nint(z2)
  atw1=atweight(izz1)
  atw2=atweight(izz2)
  zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC
  write(iout6,1000) element(izz1),atweight(izz1),-r/2.0_PREC,element(izz2), atweight(izz2),r/2.0_PREC,zcm

01000 format(/,6x,'centre     atomic weight       z',/8x,a2,4x,f14.8,3x,f10.4,/, &
           8x,a2,4x,f14.8,3x,f10.4,/,2x,'centre-of-mass',13x,f18.14/)

  ! multipole moments relative to the centre of mass

  if (iprint(214).ne.0) then

     !        estimate round-off errors present in orbital contributions to
     !        multipole moments of a given order k

     write(iout6,'(4x,"k",5x,"orb",17x,"contribution",9x,"error" )')
  endif

  do k=1,mpole

!  loop over grid points

     do imu=1,mxnmu
        inioff=(imu-1)*nni
        do in=1,nni
           igp=inioff+in

           !              rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
           !              x=(r/2.0_PREC)*vxi1(imu)*veta1(in)

           z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
           xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
           rr=sqrt(xxplusyy+z*z)

           if (abs(rr).lt.precis) then
              costh=0.0_PREC
           else
              costh=z/rr
           endif
           wk0(igp)=rr**k*plegendg(k,izero,costh)
        enddo
     enddo

     qtot=0.0_PREC
     cm2zz=0.0_PREC
     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        qq=-occ(iorb)
        qtot=qtot+qq
        call prod2 (ngorb,cw_orb(i1beg),cw_orb(i1beg),wk3)

        nni2=nni/2
        sum1=0.0_PREC
        sum2=0.0_PREC
        do imu=mxnmu,1,-1
           inioff=(imu-1)*nni
           do in=1,nni2
              igp=inioff+in
              sum1=sum1+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo

           do in=nni2+1,nni
              igp=inioff+in
              sum2=sum2+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo
        enddo

        sum3=0.0_PREC
        do imu=1,mxnmu
           inioff=(imu-1)*nni
           do in=1,nni
              igp=inioff+in
              wk5(igp)=wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
              sum3=sum3+wk5(igp)
           enddo
        enddo

        if (iprint(214).ne.0) then
           write(iout6,'(i5,i4,1x,a8,a1,e28.16,e10.2)') k, iorn(iorb),bond(iorb),gut(iorb),sum3,sum1+sum2-sum3
        endif

        cm2zz=cm2zz+qq*(sum1+sum2)
     enddo

     elect(k)=cm2zz
     electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)

     qkz1=z1*abs(r/2.0_PREC+zcm)**k*plegendg(k,izero,-1.0_PREC)
     qkz2=z2*abs(r/2.0_PREC-zcm)**k*plegendg(k,izero,1.0_PREC)
     qktot=cm2zz+qkz1+qkz2
     total(k)=qktot
     totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
  enddo


  if (iprint16.eq.0) then
     write(*,'(/,36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
     write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
        write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
     endif
  else
     write(*,'(/,44x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

     k=1
     write(*,'(5x,"dipole (Mu_z, k=1)          ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     k=2
     write(*,*)
     write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e44.32)') elect(k),total(k)
     write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)

     if (mpole.ge.3) then
        k=3
        write(*,*)
        write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif

     if (mpole.ge.4) then
        k=4
        write(*,*)
        write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e44.32)') elect(k),total(k)
        write(*,'(5x,"                            ",2e44.32)') electDA(k),totalDA(k)
     endif
  endif

  if (iprint(216).ne.0) then
     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r_1^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              r1=(r/2.0_PREC)*(vxi(mu)+veta(ni))
              if (r1.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/r1
              endif
              wk1(imu+ni)=r1
              wk2(imu+ni)=r1*r1
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w
     enddo

     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)*r/2.0_PREC
              if (rr.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/rr
              endif
              wk1(imu+ni)=rr
              wk2(imu+ni)=rr*rr
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w

        call prod2 (ngorb,cw_orb(i1beg),cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <   > =',w
     enddo

     write(*,*)
     write(*,*)
     write(*,*) 'Expectation values of r_2^k'

     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        write(*,*)
        write(*,*) 'orbital #',iorb
        do mu=1,mxnmu
           imu=(mu-1)*nni
           do ni=1,nni
              r1=(r/2.0_PREC)*(vxi(mu)-veta(ni))
              if (r1.lt.precis) then
                 wk0(imu+ni)=0.0_PREC
              else
                 wk0(imu+ni)=1.0_PREC/r1
              endif
              wk1(imu+ni)=r1
              wk2(imu+ni)=r1*r1
           enddo
        enddo
        call prod2 (ngorb,cw_orb(i1beg),wk0,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <1/r> =',w

        call prod2 (ngorb,cw_orb(i1beg),wk1,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              < r > =',w

        call prod2 (ngorb,cw_orb(i1beg),wk2,wk3)
        call prod  (ngorb,cw_orb(i1beg),wk3)
        call prod  (ngorb,f4,wk3)
        w=dot(ngorb,wgt2,ione,wk3,ione)
        write(*,*) '              <r^2> =',w
     enddo
  endif

! calculating total electronic density at nuclei

  elect(1)=0.0_PREC
  elect(2)=0.0_PREC
  do iorb=1,norb
     qq=-occ(iorb)
     i1beg=i1b(iorb)
     ngorb=i1si(iorb)
     elect(1)=elect(1)+qq*cw_orb(i1beg)*cw_orb(i1beg)
     elect(2)=elect(2)+qq*cw_orb(i1beg+nni-1)*cw_orb(i1beg+nni-1)
  enddo

  if (iprint16.eq.0) then
     write(*,'(/," total charge density at (0,0,-R/2):",e25.16 )') elect(2)
     write(*,'(  " total charge density at (0,0,+R/2):",e25.16 )') elect(1)


  else
     write(*,'(/," total charge density at (0,0,-R/2):",e25.16 )') elect(2)
     write(*,'(  " total charge density at (0,0,+R/2):",e25.16 )') elect(1)


  endif
  return
end subroutine propet4

