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
! ### initExWeights ###

!     Calculates the weights of exchange contributions to the total
!     energy expression for an open/closed shell configuration

subroutine initExWeights
  use params
  use scf
  use commons8

  character*8 :: alpha,beta
  data alpha/'+'/, beta/'-'/

  do i=1,2*maxorb*maxorb
     gec(i)=0.0_PREC
  enddo

  do i=1,norb
     do j=1,norb
        nlm(i,j)=0
     enddo
  enddo

  do i=1,4*norb
     spn(i)=spin(i)
  enddo

  do i=1,norb
     isumi=4*(i-1)
     !        gec(i)=occ(i)-1.0
     iloop(i)=2
     iqm(1+isumi)=1
     iqm(2+isumi)=1
     if (mgx(3,i).gt.0) then
        iloop(i)=4
        iqm(3+isumi)=-1
        iqm(4+isumi)=-1
     endif
  enddo

!     initialize spn and spin arrays
!     throw out the electrons + = alpha and - = beta
!     the orbital is locked = 1 (open shell case)
!                           = 0 (close shell case)

  do i=1,norb
     if (lock(i).eq.0) then
        jj=4*(i-1)
        spn(jj+1)=alpha
        spn(jj+2)=beta
        if (mgx(3,i).gt.0) then
           spn(jj+3)=alpha
           spn(jj+4)=beta
        endif
     endif
  enddo

  do i=1,4*norb
     spin(i)=spn(i)
  enddo

  !     loop now through all electrons

  !     gec-a=gec(kk) gec-b=gec(kk+isu2)
  if (iprint(170).ne.0) then
     write(6,1101)
01101 format(3x,'orb1    ',6x,'orb2       gec-a  gec-b  div')
  endif

  isu2=norb*norb
  do i=1,norb
     ixx=4*(i-1)
     div(i)=zero
     do ii=1,iloop(i)
        ix=ii+ixx
        if ((spn(ix).ne.alpha).and.(spn(ix).ne.beta)) goto 12
        div(i)=div(i)+one

        do j=1,norb
           jxx=4*(j-1)
           kk=i+norb*(j-1)

           mpl=abs(mgx(3,i)-mgx(3,j))
           mmi=abs(mgx(3,i)+mgx(3,j))

           do jj=1,iloop(j)
              jx=jxx+jj
              if ((spn(jx).ne.alpha).and.(spn(jx).ne.beta)) goto 15

!                 div(i)=div(i)+1.00_PREC

              if ((mpl.eq.mmi).and.(spn(ix).eq.spn(jx)).and.(i.ne.j)) then
                 gec(kk)=gec(kk)+1.0_PREC
              endif

              if ((mpl.ne.mmi).and.(spn(ix).eq.spn(jx))) then
                 if ((i.eq.j).and.(iqm(ix).ne.iqm(jx))) then
                    gec(kk)=gec(kk)+1.0_PREC
                 endif

                 if ((i.ne.j).and.(iqm(ix).eq.iqm(jx))) then
                    gec(kk)=gec(kk)+1.0_PREC
                 endif

                 if ((i.ne.j).and.(iqm(ix).ne.iqm(jx))) then
                    gec(kk+isu2)=gec(kk+isu2)+1.0_PREC
                 endif
              endif
15            continue
           enddo
        enddo
12      continue
     enddo

     do idiv=1,norb
        kk=i+norb*(idiv-1)
        if (abs(div(i)).gt.precis) then
           gec(kk)=gec(kk)/div(i)
           gec(kk+isu2)=gec(kk+isu2)/div(i)
        else
           gec(kk)=zero
           gec(kk+isu2)=zero
        endif
        if (iprint(172).ne.0) then
           write(6,1102)  iorn(i),bond(i),gut(i),iorn(idiv),bond(idiv),gut(idiv),gec(kk),gec(kk+isu2),div(i)
01102      format(i4,1x,a8,a1,i4,1x,a8,a1,3f7.3)
        endif
     enddo
  enddo

  ial=0
  ibe=0
  imag=0
  do i=1,norb
     isumi=4*(i-1)
     if (lock(i).eq.0) goto 100

     do ii=1,iloop(i)
        imu=ii+isumi
        if (spn(imu).eq.alpha) ial=ial+1
        if (spn(imu).eq.beta) ibe=ibe+1
        if (spn(imu).eq.alpha.or.spn(imu).eq.beta) then
           imag=imag+iqm(imu)
        endif
     enddo
100  continue
  enddo

!     Non-zero entries of nlm array indicate pairs of orbitals having
!     non-zero off-diagonal Lagrange multipliers. Each pair of orbitals
!     of the same symmetry is examined to see if an orthogonal
!     transformation can be applied and off-diagonal Lagrange
!     multipliers set to zero. One can make the multipliers to appear in
!     Fock equations by specifying the pair of orbitals in question via
!     largra label.

  do j=1,(norb-1)
     jxx=4*(j-1)
     do i=(j+1),norb
        nlm(j,i)=0
        nlm(i,j)=0

        if (mgx(6,j).ne.mgx(6,i)) goto 400
        if (gut(j).ne.gut(i)) goto 400

        if (iloop(j).ne.iloop(i)) goto 400
        ixx=4*(i-1)
        match=0
        do ii=1,iloop(i)
           ix=ii+ixx
           jx=ii+jxx
           if (spn(ix).eq.spn(jx)) match=match+1
        enddo
        if (match.ne.iloop(i)) then
           nlm(j,i)=1
           nlm(i,j)=1
        endif
400     continue
     enddo
  enddo

  do j=1,(norb-1)
     do i=(j+1),norb
        do k=1,2*norb,2
           if (lagraon(k).eq.j.and.lagraon(k+1).eq.i) then
              nlm(j,i)=1
              nlm(i,j)=1
           endif
        enddo
     enddo
  enddo

!     Testing mode only
!     Select a way in which off-diagonal Lagrange multipliers are
!     calculated between closed shell orbitals (i.e. when lagraon label
!     is present).

  lmtype=0
  if (idbg(336).eq.1) lmtype=1
  if (idbg(337).eq.1) lmtype=2

  return

  write(iout6,1001)
1001 format(/1x,'... error in spin orbitals ...'//)
  stop 'initExWeights'
end subroutine initExWeights
