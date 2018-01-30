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
! ### exchMCSOR ###

!     Evaluates exchange potentials involving a given orbital.

subroutine exchMCSOR (iorb,cw_sor,psi,pot,excp,bpot,d,e,f3,g,lhs,rhs,wk2)
  use params
  use discret
  use scf
  use solver
  use commons8

  implicit none
  integer :: i,ib1,ib2,ibexp,ibexp1,idel,idel1,ifirst,ig,in1,in2,ioffs1, &
       iorb,iorb1,ioffst,ipc,iswtch,itr1,itr2,ng1,ng2,ngexp,ngorb,ngrid
  real (PREC) :: xdel
  integer, dimension(*) :: cw_sor
  real (PREC), dimension(*) :: psi,pot,excp,rhs,bpot,d,e,f3,g,lhs,wk2

  if (nel.lt.2 ) return
  
  if     (exlexp.eq.0) then
     ifirst=1
  elseif (exlexp.eq.2) then
     !        recalculate only those exchange potentials which depend on
     !        orbitals modified so far (note the reverse order of relaxation
     !        in scf)
     ifirst=iorb
  endif
  
  ngorb=i1si(iorb)	  
  do iorb1=ifirst,norb
     iswtch=0
     
     if (iorb.eq.iorb1.and.mgx(6,iorb).eq.0 ) goto 10
     if ((iorb.eq.iorb1).and.(ilc(iorb*(iorb+1)/2).lt.1)) goto 10
     
     !        orbitals in increasing order
     
     if (iorb.lt.iorb1) then
        in1=iorb
        in2=iorb1	  
     else	    
        in1=iorb1
        in2=iorb
     endif

     ib1=i1b(in1)
     ng1=i1si(in1)
     ib2=i1b(in2)
     ng2=i1si(in2)	  
     
     ipc=in1+in2*(in2-1)/2
     iwexch(ipc)=-1
     
     ibexp=i3b(ipc)
     ngexp=i3si(ipc)
     
     idel=abs(mm(in1)-mm(in2))
     if (in1.eq.in2) idel=2*mm(in1)
     idel1=mm(in1)+mm(in2)
     
     !        prepare right-hand side of the Poisson equation
     
     do i=1,ngorb
        wk2(i)=zero
     enddo
     
     ngrid=min(ng1,ng2)
     call prod2 (ngrid,psi(ib1),psi(ib2),wk2)
     call prod2 (ngorb,wk2,g,rhs)
     
     if (mod(idel,itwo).eq.0) then
        isym= 1
     else
        isym=-1
     endif
     
1000 continue
     
     do itr1=1,maxsor1	  	
        do ig=1,i1ng(iorb)
           do i=1,i1si(iorb)
              wk2(i)=0.0_PREC
           enddo
           
           omega=ovfexch(ig)
           omega1=one-omega
           ioffst=ioffs(ig)	    
           ioffs1=ioffst+1
           ibexp1=ibexp+ioffst	  	  
           
           do i=1,4
              dmu2t(i)=dmu2(i,ig)
              dmu1t(i)=dmu1(i,ig)
           enddo
               
           !              prepare left-hand side of the poisson equation
           !              include the diagonal part of the differentiation 
           !              operator in lhs
           
           if (idel.eq.0) then
              do i=1,ngsize(ig)
                 lhs(ioffst+i)=f3(ioffst+i)+diag(ig)
              enddo
           else
              xdel=dble(idel*idel)
              do i=1,ngsize(ig)
                 lhs(ioffst+i)=f3(ioffst+i)+xdel*e(ioffst+i)+diag(ig)
              enddo
           endif
           
           if (i1ng(iorb).eq.1) then 
              !                 icase=1
              ifill=1
              ngrd1 =ingr1(1,ig)
              ngrd6a=ingr1(2,ig)
              ngrd6b=ingr1(3,ig)
              ngrd7 =ingr1(4,ig)	  	  
              call putin (nni,nmu(ig),isym,excp(ibexp),wk2)
              do itr2=1,maxsorpot(iorb)
                 call mcsor (wk2,lhs(ioffs1),rhs(ioffs1),bpot(ioffs1),d(ioffs1),                     &
                      cw_sor(iadext(ig)),cw_sor(iadnor(ig)),cw_sor(iadex1(ig)),cw_sor(iadex2(ig)),   &
                      cw_sor(iadex3(ig)))
              enddo
              call putout (nni,nmu(ig),excp(ibexp),wk2)     
           else
              if (ig.eq.1) then 
                 !                    icase=2
                 ifill=2
                 ngrd1 =ingr2(1,ig)
                 ngrd6a=ingr2(2,ig)
                 ngrd6b=ingr2(3,ig)
                 ngrd7 =ingr2(4,ig)	  	  	    
              elseif (ig.ne.1.and.ig.ne.i1ng(iorb)) then 
                 !                    icase=2
                 ifill=3
                 ngrd1 =ingr2(1,ig)
                 ngrd6a=ingr2(2,ig)
                 ngrd6b=ingr2(3,ig)
                 ngrd7 =ingr2(4,ig)	  	  	    
              elseif (ig.eq.i1ng(iorb)) then 
                 ifill=4
                 ngrd1 =ingr1(1,ig)
                 ngrd6a=ingr1(2,ig)
                 ngrd6b=ingr1(3,ig)
                 ngrd7 =ingr1(4,ig)	  	  	    
                 muoffs=iemu(ig-1)-1
              endif
           endif
        enddo
     enddo
     if (iswtch.eq.1) goto 10
     if (ilc(ipc).lt.2) goto 10
     iwexch(ipc)=-2
     idel=idel1
     ibexp=ibexp+ngexp
     ipc=ipc+norb*(norb+1)/2
     iswtch=1
     goto 1000
10   continue
  enddo
1021 format(/,i5,10e12.3,/(5x,10e12.3))

end subroutine exchMCSOR
