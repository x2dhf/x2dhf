! ### initAddr ### 

!     Determines dimensions of arrays and addresses of particular
!     orbitals and coulomb potentials in cw_orb and cw_coul arrays 
!     addresses of exchange potentials v(i,j) in the cw_exch array.

!     With each cw_sth array is associated corresponding address array
!     defining its division into subarrays holding individual
!     orbitals/potentials:

!     cw_orb   - i1xx(iorb)
!     cw_coul  - i2xx(iorb)
!     cw_exch  - i3xx(k)

!     cw_suppl - iaddr4(isuppl)
!     cw_sctch - iaddr5(iscratch)


!     i1b(iorb)  - the address of the first element of orbital iorb within
!                  cw_orb array (starting address of the iorb-th orbital)

!     i1e(iorb)  - the address of the last element of orbital iorb within
!                  cw_orb array (ending address of orbital iorb)

!     i1si(iorb) - the size of the cw_orb subarray holding values of the
!                  iorb-th orbital

!     i1ng(iorb) - number of grids used by the iorb-th orbital (must be
!                  1 in x2dhf ver. 2.0)
      


!     i2b(iorb)  - the address of the first element of Coulomb potential
!                  iorb within cw_coul array (starting address of the iorb-th
!                  potential)

!     i2e(iorb)  - the address of the last element of Coulomb potential
!                  iorb within cw_coul array (ending address of the iorb-th
!                  potential)

!     i2si(iorb) - the size of the cw_coul subarray holding values of
!                  the iorb-th Coulomb pptential

!     i2ng(iorb) - number of grids used by the iorb-th potential (must
!                  be 1 in x2dhf ver. 2.0)

!     i3b(k)     - the address of the first element of the exchange potential for
!                  the k-th pair of orbitals (iorb1,iorb2): k=iorb1+iorb2*(iorb2-1)/2

!     i3e(k)     - the address of the last element of the exchange potential for
!                  the k-th pair of orbitals (iorb1,iorb2)

!     i3si(k)    - the size of the cw_exch subarray holding values of 
!                  exchange potential for the k-th pair of orbitals

!     i3ng(k)    - number of grids needed to define the k-th exchange potential

!     ilc(k)     - 0, 1 or 2 indicating the number of exchange potentials for
!                  a given pair of orbitals

!     i3xk(iorb1,iorb2) 
!                - the index k of exchange potential corresponding to the pair of 
!                  orbitals iorb1 and iorb2

!     Subarrays of cw_suppl and cw_sctch are defined accordingly by i4xx
!     and i5xx arrays.

!     We assume no more than maxorb orbitals and Coulomb potentials and
!     consequently no more than (maxorb*(maxorb+1)/2) exchange
!     potentials.

subroutine initAddr 

  use params
  use discret
  use memory
  use commons8

  implicit none

  integer :: i,i3beg,i3bfin,ibeg,i3boffset,iend,iorb,iorb1,iorb2,irec,k,l,l1cur,l2cur,l3cur,l5cur,mxsize8,ngrid,ngds,nons

  i1b(1)=1
  ngrid=nni*iemu(i1ng(1))
  i1e(1)=ngrid
  i1si(1)=ngrid
  i1mu(1)=iemu(i1ng(1))
  
  !     current value of cw_orb length is stored in l1cur
  l1cur=ngrid
  do iorb=2,norb
     i1b(iorb)=i1e(iorb-1)+1
     ngrid=nni*iemu(i1ng(iorb))
     i1e(iorb)=i1b(iorb)+ngrid-1
     i1si(iorb)=ngrid
     i1mu(iorb)=iemu(i1ng(iorb))
     l1cur=l1cur+ngrid
  enddo
  
  !     it is assumed that potential functions are defined on the same number
  !     of grids as the orbitals
  
  i2b(1)=1
  ngrid=nni*iemu(i2ng(1))
  i2e(1)=ngrid
  i2si(1)=ngrid
  i2mu(1)=iemu(i1ng(1))
  
  !     current value of cw_coul length is stored in l2cur	
  l2cur=ngrid
  do iorb=2,norb
     i2b(iorb)=i2e(iorb-1)+1
     ngrid=nni*iemu(i2ng(iorb))
     i2e(iorb)=i2b(iorb)+ngrid-1
         i2si(iorb)=ngrid
         i2mu(iorb)=iemu(i2ng(iorb))
         l2cur=l2cur+ngrid
      enddo

!     dimension of i3xx arrays is norb*(norb+1)/2 to allow for two
!     exchange potentials between non-sigma orbitals
      
!     (obsolete) exchange potential corresponding to orbitals iorb1 and iorb2
!     defined on 1ing(iorb1) and i1ng(iorb2) grids respectively, is
!     defined on smaller of these grids, i.e
!     i3ng=min(i1ng(iorb1),i1ng(iorb2)) 

      iend=0
      l=0
      l3cur=0
      i3boffset=0
      do  iorb1=1,norb
         do  iorb2=iorb1,norb
            k=iorb1+iorb2*(iorb2-1)/2
            ilc(k)=0
            i3xk(iorb1,iorb2)=0
            i3xk(iorb2,iorb1)=0
            if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50
            ilc(k)=1
            l=l+1
            ngds=min(i1ng(iorb1),i1ng(iorb2))
            i3ng(k)=ngds
            ngrid=nni*iemu(ngds)
            iend=iend+ngrid
            i3e(k)=iend
            i3b(k)=i3e(k)-ngrid+1
            i3si(k)=ngrid
            i3mu(k)=iemu(i3ng(k))
            i3bfin=i3e(k)

            i3xk(iorb1,iorb2)=k
            i3xk(iorb2,iorb1)=k
            i3orb1(k)=iorb1
            i3orb2(k)=iorb2

            l3cur=l3cur+ngrid
            if (iorb1.eq.iorb2) goto 50
            if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50
            ilc(k)=2
            l=l+1
            iend=iend+ngrid
            l3cur=l3cur+ngrid
! FIXME 1830 should be replaced by a variable
            if (k.gt.1830) then
               write(*,*) 'initAddr:'
               write(*,*) 'address array for exchange potentials is too short'
               stop 'initAddr'	
            endif
00050       continue
         enddo
         i3boffset=i3bfin
      enddo
      
      nexch=l

!     When exchange potential functions are stored as separate files
!     each containing a single potential then their labels are generated
!     in the fashon defined below.

!     prepare a two dimension array containing numerical labels of
!     exchange potential functions being needed for every orbital
!     
!     all exchange potentials have to be defined on the same grid
!     
!     let k=i3xk(iorb1,iorb2) be an index associated uniquely with a given
!     pair of orbitals

!     if k=0  - no exchange potentials (no record)      

!     if k.ne.0 and ilc(k)=1  - exchange potential record label is stored in i3xrec1
!     if k.ne.0 and ilc(k)=2  - exchange potential record label is stored in i3xrec2
      
      ngrid=nni*mxnmu
      l=0
      do iorb1=1,norb
         ibeg=1
         do iorb2=iorb1,norb
            k=i3xk(iorb1,iorb2)
            if (k.ne.0) then
               if (ilc(k).ne.0) then
                  l=l+1
                  i3xrec1(k)=l
                  if (ilc(k).eq.2) then
                     l=l+1
                     i3xrec2(k)=l
                  endif
               endif
            endif
         enddo
      enddo
      
!     determine exchange potentials associated with a particular orbital
!     and store the record lables in i3rec; the number of exchange
!     potentials for each orbital is kept in i3nexp
      
      do iorb1=1,norb
         l=0
         ibeg=1
         do iorb2=1,maxorb
            i3breck(iorb1,iorb2)=0
         enddo
         do iorb2=1,norb
            k=i3xk(iorb1,iorb2)
            if (k.ne.0) then
               if (ilc(k).ne.0) then
                  l=l+1        
                  
                  if (l.gt.maxorb) then
                     print *,'initAddr: the second dimension of arrays i3xpair,i3brec can not exceed',maxorb
                     stop 'initAddr'
                  endif
                  
                  i3xpair(iorb1,l)=i3xrec1(k)
                  i3brec (iorb1,l)=ibeg
                  i3breck(iorb1,l)=k
                  
                  i3btv(iorb1,iorb2)=ibeg
                  i3btv(iorb2,iorb1)=ibeg
                  
!                 before a new portion of exchange potentials is read
!                 i3b has to be modified accordingly (k is determined by
!                 i3breck and must be non zero)
                  
                  i3nexcp(iorb1)=l
                  ibeg=ibeg+ngrid
                  if (ilc(k).eq.2) then
                     l=l+1      
                     i3xpair(iorb1,l)=i3xrec2(k)
                     i3brec(iorb1,l)=ibeg
                     i3breck(iorb1,l)=k
                     i3nexcp(iorb1)=l
                     ibeg=ibeg+ngrid
                  endif
               endif
            endif
         enddo
      enddo

!     is working array cw_suppl large enough?
      
      if (nsuppl*mxsize.gt.length4) then
         write(iout6,*) 'cw_suppl is too short!'
         write(iout6,*) 'declared length is:',length4
         write(iout6,*) 'needed length is nsuppl*mxsize'
         stop 'initAddr'
      endif

      do i=1,nsuppl
         i4b(i)=(i-1)*mxsize+1
         i4e(i)=iorb*mxsize
         i4si(i)=mxsize
         i4ng(i)=ngrids
      enddo
      
!     current value of cw_sctch length is stored in l5cur
      mxsize8=(nni+8)*(mxnmu+8)
      do i=1,nsctch
         i5b(i)=(i-1)*mxsize8+1
         i5e(i)=i*mxsize8
         i5si(i)=mxsize8
         i5ng(i)=ngrids
      enddo
      l5cur=i5e(20)

!     are working arrays large enough?

      if (l1cur.gt.length1) then
         write(iout6,*) 'cw_orb is too short!'
         write(iout6,*) 'declared length is:',length1      
         write(iout6,*) 'needed length is:',l1cur
         if (idbg(550).eq.0) stop 'initAddr'
      endif

      if (l2cur.gt.length2) then
         write(iout6,*) 'cw_coul is too short!'
         write(iout6,*) 'declared length is:',length2
         write(iout6,*) 'needed length is:',l2cur
         if (idbg(550).eq.0) stop 'initAddr'
      endif

!     l3cur is calculated as if all exchange potential functions were
!     to be kept in core. This is the case if iform=1 or iform=3.
      
      if (imethod.eq.1) then
         if (iform.eq.1.or.iform.eq.3) then
            if (l3cur.gt.length3) then
               write(iout6,*) 'cw_exch is too short!'
               write(iout6,*) 'declared length is:',length3
               write(iout6,*) 'needed length is:',l3cur
               if (idbg(550).eq.0) stop 'initAddr'
            endif
         else
            nons=0
            do iorb=1,norb
               if (ll(iorb).ne.0) nons=nons+1
            enddo
!           l3cur=nni*maxmu*(norb+nons)
            l3cur=mxsize*(norb+nons)
            if (l3cur.gt.length3) then
               write(iout6,*) 'cw_exch is too short!'
               write(iout6,*) 'declared length is:',length3
               write(iout6,*) 'needed length is:',l3cur
               stop 'initAddr'
            endif
         endif
      endif

      if (l5cur.gt.length5) then
         write(iout6,*) 'cw_sctch is too short!'
         write(iout6,*) 'declared length is:',length5
         write(iout6,*) 'needed length is',nsctch*mxsize8
      endif


!     checking i3bxx entries (potentials as separate files)

! FIXME improve layout
      if (iprint(10).ne.0) then
         print *,'initAddr: ngrid ',nni*mxnmu
         print *,'initAddr: nexch ',nexch

         print *, 'initAddr: iorb1,iorb2,k,i3xk,ilc(k),i3b'      
         do  iorb1=1,norb
            do  iorb2=iorb1,norb
               k=iorb1+iorb2*(iorb2-1)/2
               if (iorb1.eq.iorb2) print *,'==='
               write(*,1200) iorb1,iorb2,k,i3xk(iorb1,iorb2),ilc(k),i3b(k)
            enddo
            print *, '   '
         enddo

         print *,'iaddr: iorb,i3nexcp'
         do iorb1=1,norb
            print *,iorb1,i3nexcp(iorb1)
         enddo
         print *,' '
         
         print *,'iaddr: iorb // l,i3xpair,i3brec,i3breck'
         do iorb1=1,norb
            print *,'iorb1: ',iorb1
            do l=1,i3nexcp(iorb1)
               k=i3breck(iorb1,l)
               write(*,1210) l, i3xpair(iorb1,l),i3brec(iorb1,l),k
            enddo
         enddo
         
01200    format(2i4,i6,3i12)
01210    format(i4,5i12)
         
         print *,'iaddr: iorb1,iorb2,k,i3xrec1,i3xrec2'
         do  iorb1=1,norb
            do  iorb2=iorb1,norb
               k=i3xk(iorb1,iorb2)
               if (ilc(k).ne.0) then
                  print *,iorb1,iorb2,k,i3xrec1(k),i3xrec2(k)
               endif
            enddo
         enddo

         print *, '   '
         print *,'iorb1 // iorb2,k,irec,i3beg'
         do iorb1=1,norb
            print *,'    '
            print *, 'iorb1: ',iorb1
            do iorb2=1,i3nexcp(iorb1)
               k=i3breck(iorb1,iorb2)
               i3beg=i3brec(iorb1,iorb2)
               irec=i3xpair(iorb1,iorb2)
               ngrid=i3si(k)
               print *,iorb2,k,irec,i3beg
            enddo
         enddo
      endif

end subroutine initAddr







