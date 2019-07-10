module readLables_m
  implicit none
contains
  subroutine read_break
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none

    ! label: break 
    ! switch on symmetry breaking
    
    ibreak=1
    return
  end subroutine read_break

  subroutine read_config
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m

    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    use nmucalc_m
    use nnucalc_m
    
    implicit none
    !  character*12 char8,sigma,pi,delta,phi
    !  integer clo, cloe,i,isum,isum0,iput
    if (icompLEnc.ne.2) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'
       stop 'read_config'
    endif
    
    mt=0
    icompLEnc=icompLEnc+1

    call inInt(itotq)
    totq=z1+z2-dble(itotq)
    ! to allow for noninteger nucleus charge needed to force convergence
    ! in some open-shell cases read in the nuclei charges again

    z1t=0.0_PREC
    z2t=0.0_PREC
    call inFloat(z1t)
    call inFloat(z2t)
    if (z1t.gt.precis) z1=z1t
    if (z2t.gt.precis) z2=z2t
    
    ! read in symmetry and occupation and test if the orbital is locked = 1 (open shell
    ! case) = 0 (closed shell case)
    isum=0
100 call inCard
    call inInt(nbsym)
    call inStr(char8)
    if(char8.ne.sigma.and.char8.ne.pi.and.char8.ne.delta.and.char8.ne.phi) then
       write(iout6,*) 'Error: wrong symmetry or symmetry higher than phi'
       stop 'inputData'
    endif
    
    isum0=isum
    do i=1,nbsym
       isum=isum+1
       orbsym(isum)=char8
       bond(isum)=orbsym(isum)
       gut(isum)=space
    enddo
    
    if (isum.gt.maxorb) then
       write(iout6,*) 'too many orbitals (max. see maxorb)'
       stop 'inputData'
    endif
    ! in homonuclear case read in the u/g parity of each orbital
    ! if BREAK is on (ibreak=1) u/g labels are ignored
    ! if |z1-z2|<homolevl then the case is treated as a homonuclear one
    ! (see setDefaults)
    
    if (abs(z1-z2).lt.homolevl.and.ibreak.eq.0) then
       call inStr(gut(isum))
       isum1=isum0
       do i=1,nbsym
          isum1=isum1+1
          gut(isum1)=gut(isum)
       enddo
       
       cloe=4.0_PREC
       if (orbsym(isum).eq.sigma) cloe=2.0_PREC
       iput=4*(isum-1)
       call inStr(char8)
       if (char8.ne.endl) then
          spin(1+iput)=char8
          call inStr(char8)
          if (char8.ne.endl) then
             spin(2+iput)=char8
             call inStr(char8)
             if (char8.ne.endl) then
                spin(3+iput)=char8
                call inStr(char8)
                if (char8.ne.endl) then
                   spin(4+iput)=char8
                   call inStr(char8)
                endif
             endif
          endif
       endif
    else
       call inStr(guttmp)
       char8=guttmp
       if (guttmp.eq.'u'.or.guttmp.eq.'g') then
          gut(isum)=guttmp
          isum1=isum0
          do i=1,nbsym
             isum1=isum1+1
             gut(isum1)=gut(isum)
          enddo
          call inStr(char8)
       endif
       
       cloe=4.0_PREC
       if (orbsym(isum).eq.sigma) cloe=2.0_PREC
       iput=4*(isum-1)
       ! call inStr(char8)
       if (char8.ne.endl) then
          spin(1+iput)=char8
          call inStr(char8)
          if (char8.ne.endl) then
             spin(2+iput)=char8
             call inStr(char8)
             if (char8.ne.endl) then
                spin(3+iput)=char8
                call inStr(char8)
                if (char8.ne.endl) then
                   spin(4+iput)=char8
                   call inStr(char8)
                endif
             endif
          endif
       endif
    endif

    if (spin(1+iput).ne.'+'.and.spin(1+iput).ne.'-'.and.spin(1+iput).ne.'.') then
       do i=1,nbsym
          occ(isum0+i)= cloe
       enddo
    else
       do i=1,4
          if (spin(i+iput).eq.'+'.or.spin(i+iput).eq.'-') then
             do j=1,nbsym
                iput1=4*(isum0+j-1)
                occ(isum0+j)=occ(isum0+j)+one
                spin(i+iput1)=spin(i+iput)
             enddo
          endif
       enddo
    endif
    norb=isum
    no=norb
    
    ! initialize scaling factors for off-diagonal Lagrange multipliers
    
    do i=1,norb
       do j=i+1,norb
          sflagrat(i,j)=sflagra
       enddo
    enddo
    
    ! initialize SOR parameters
    
    do i=1,norb
       ifix(i)=0
       inhyd(i)=0
       maxsororb(i)=maxsor2
       maxsorpot(i)=maxsor3
    enddo
    
    if (char8.eq.endl) then
       totchar=0.0_PREC
       do iorb=1,norb
          totchar=totchar+occ(iorb)
       enddo
       
       ! Warning!
       ! on Cray Y-MP (F90) the calculated charge was 29. and the one
       ! read in -- 29.00000000000011.
       ! The statement if (totchar.eq.totq) had to be replaced by the following
       ! one
       
       ! FIXME r128
       if (abs(totchar-totq).gt.precis*1000.0_PREC) then
          write(iout6,*) 'Warning: mismatch in given andcalculated total charge:',totchar,totq
          ! FIXME !!!! OED is effected
          ! stop 'inputData'
       endif
       
       ! set some additional parameters
       
       do iorb=1,norb
          maxsororb(iorb)=maxsor2
          maxsorpot(iorb)=maxsor3
          if (z1.gt.precis) then
             co1(iorb)=1.0_PREC
          else
             co1(iorb)=0.0_PREC
          endif
          if (z2.gt.precis) then
             co2(iorb)=1.0_PREC
          else
             co2(iorb)=0.0_PREC
          endif
       enddo
       
       do iorb=1,norb
          if (i1ng(iorb).eq.0) then
             i1ng(iorb)=ngrids
             i2ng(iorb)=ngrids
          endif
       enddo
       
       do iorb=1,norb
          if (orbsym(iorb).eq.sigma) mt=0
          if (orbsym(iorb).eq.pi)    mt=1
          if (orbsym(iorb).eq.delta) mt=2
          if (orbsym(iorb).eq.phi)   mt=3
          
          ! asymuthal quantum number of a hydrogen orbital is taken
          ! in accordance with the symmetry of a molecular orbital
          
          mgx(3,iorb)=mt
          mgx(6,iorb)=mt
          lock(iorb)=0
          iopenshell=0
          clo=4.0_PREC
          if (orbsym(iorb).eq.sigma) clo=2.0_PREC
          clo=occ(iorb)/clo
          if (abs(clo-1.0_PREC).gt.1.d-06) then
             lock(iorb)=1
          endif
          mgi=(-1)**mgx(3,iorb)
       enddo
       return
    endif
    goto 100
    return
  end subroutine read_config
  
  subroutine read_conv
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    ! label: conv [nscf2skip [nnenlast [nnolast] ] ] 
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       nscf2skip=itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          nenlast=itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             nnolast=itmp
          endif
       endif
    endif
    return
  end subroutine read_conv

  subroutine read_debug
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: debug
    ! idbg  -  additional printout or action for testing purposes
    !          set debug flags. if an integer i is encountered debuf
    !          flag i is set, i.e. idbg(i)=1. At most 999 integers can be specified.
    
    ! FIXME
    ! List of used idbg flags:
    ! doSCF: 77 - assign maxsor for each orbital every iepoch iterations (call schedSOR)
    ! initExWeights: 335, 336 - select a way in which off-diagonal
    !                           Lagrange multipliers are calculated between closed shell orbitals
    
    !        locenergy: 496, 497
    !        printCase: 550 - mm3=3*mxsize
    !        initAddr: 550  - stop
    !        initGauss: 562 - call gauss_ovlap
    !        initGauss: 560, 561, 562, 565, 566
    !        initGauss: 560 - stop
    !        label: debug    
    !        idbg  -  additional printout or action for testing purposes
    !        set debug flags. if an integer i is encountered debuf flag i is
    !        set, i.e. idbg(i)=1. At most 999 integers can be specified.
    !                                  orbitals
    
    inzero=1
    do i=1,maxflags
       id(i)=0
       call inInt(id(i))
       if (id(i).gt.0) then
          inzero=inzero+1
          idbg(id(i))=1
       endif
    enddo
    return
  end subroutine read_debug


  subroutine read_dft
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m

    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m

    implicit none
    
    ! label: dft [lda|b88] [lyp] [vwn]
    ! if the DFT method is used one can choose exchange and correlation
    ! energy functionals
    ! call inInt(itmp) 
    
    ! if only dft label is present calculate exchange (LDA, B88, PW86,
    ! PW91) and correlation (LYP) contributions to total energy at the
    ! end of SCF process

    idft=1
    call inStr4lxc(char30)

    ! first try to find old functional labels 
    char8=trim(char30)
    idftex=0
    idftcorr=0
    do i=1,10
       if (char8.eq.cdftex(i)) idftex=i
       if (char8.eq.cdftcorr(i)) idftcorr=i
    enddo
    ! if (idftex.eq.0.and.idftcorr.eq.0) then
    !    write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
    !    stop 'read_dft'
    ! endif
    if (idftex.ne.0.or.idftcorr.ne.0) then
       call inStr(char8)
       if (char8.ne.endl) then
          do i=1,10
             if (char8.eq.cdftex(i)) idftex=i
             if (char8.eq.cdftcorr(i)) idftcorr=i
          enddo
          if (idftex.eq.0.and.idftcorr.eq.0) then
             write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
             stop 'read_dft'
          endif
       endif
       if ((idftex.gt.4).or.(idftcorr.gt.2)) then
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_dft'
       endif
       if (idftex.ne.0) then 
          islat =1
          imethod=4
       endif
    
       if (idftcorr.ne.0) then 
          islat=1
          imethod=4
          if (idftex.eq.0) idftex=1
       endif

       return
    endif
    
    !if (lxcFuncs.ne.0.or.idftex.ne.0.or.idftcorr.ne.0) idft=0
    
          ! some functionas are only valid for closed shell configurations
    if (trim(char30).ne."") then
       lxcFuncs=0
       do j1=1,2
          lxcFuncs=lxcFuncs+1
          lxcFuncs2use(lxcFuncs)=xc_f90_functional_get_number(trim(char30))
          if (lxcFuncs2use(lxcFuncs)<=0) then
             write(*,*) "Error! ",trim(char30),": no such libxc functional found"
             stop "inputData"
          endif
          
          call xc_f90_func_init(xc_func, xc_info, lxcFuncs2use(lxcFuncs), XC_UNPOLARIZED)
          select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             call inStr4lxc(char30)
             if (trim(char30).eq."") exit
          case default
             write(*,*) "Error! ", trim(char30),": unsupported libxc functional"
             stop 'inputData'
          end select
          call xc_f90_func_end(xc_func)
       enddo
       
       if (lxcFuncs>0) then
          islat=1
          imethod=4
          idft=0
          return
       endif
    endif
    return
  end subroutine read_dft

  subroutine read_exchio
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inStr_m
    
    implicit none

    call inStr(clabel1)
    call inStr(clabel2)
    if     (clabel1.eq.'in-many'.and.clabel2.eq.'out-many') then
       iform=0
    elseif (clabel1.eq.'in-one' .and.clabel2.eq.'out-many') then
       iform=1
    elseif (clabel1.eq.'in-many'.and.clabel2.eq.'out-one') then
       iform=2
    elseif (clabel1.eq.'in-one' .and.clabel2.eq.'out-one') then
       iform=3
    else
       write(iout6,*) 'Error: incorrect format of exchange potential file.'
       write(iout6,*) 'Try one of the folowing combinations: '
       write(iout6,*) 'IN-ONE OUT-ONE, IN-ONE OUT-MANY, IN-MANY OUT-ONE IN-MANY OUT-MANY'
       stop 'read_exchio' 
    endif
    if ((iform.eq.0.or.iform.eq.2).and.imethod.ne.1) then
       write(iout6,*) 'Error: incorrect format of exchange potential file.'
       write(iout6,*) 'Try one of the folowing combinations:'
       write(iout6,*) 'IN-ONE OUT-ONE, IN-ONE OUT-MANY, IN-MANY OUT-ONE IN-MANY OUT-MANY'
       stop 'read_exchio'
    endif
    ! extra consistancy checks 
    if (imethod.eq.3.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as HFS the second parameter must be set to 3'
       stop 'inputData'
    endif
    
    if (imethod.eq.4.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as DFT the second parameter must be set to 3'
       stop 'inputData'
    endif
    
    if (imethod.eq.5.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as SCMC  the second parameter must be set to 3'
       stop 'inputData'
    endif
    
    return
  end subroutine read_exchio

  subroutine read_fefield
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
        
    implicit none
    ! label: fefield
    
    ! finite field calculations
    ! 1 field strength
    ! 2 cutoff
    
    ifefield=1
    call inFloat(ffield)
    ! issue error message if homo is on
    if (ihomon.eq.2) then
       write(iout6,*) 'Error: FEFIELD and HOMO labels are mutually exclusive'
       stop 'inputData'
    endif
    return
  end subroutine read_fefield

  subroutine read_fermi
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    
    implicit none
    ! label: fermi
    
    call inFloat(z1atmass)
    if (z1atmass.le.0.0_PREC) then
       izz1=nint(z1)
       z1atmass=atweight(izz1)
    endif
    call inFloat(z2atmass)
    if (z2atmass.le.0.0_PREC) then
       izz2=nint(z2)
       z2atmass=atweight(izz2)
    endif
    if (z1atmass.le.0.0_PREC.and.z2atmass.le.0.0_PREC) then
       write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
       stop 'read_fermi'
    endif
    ifermi=1
    ipot=1
    return
  end subroutine read_fermi

  subroutine read_fix
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    ! exlwf exlpot- if nonzero wave functions, coulomb potentials or
    ! exlexp        exchange potentials are not relaxed; exlpot and 
    !               exlexp cannot be non zero if hydrogen orbitals 
    !               are used to start the scf process
    
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       iexlorb=itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          iexlcoul=itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             iexlexp=itmp
          endif
       endif
    endif
    
    exlorb=dble(iexlorb)
    exlcoul=dble(iexlcoul)
    exlexp=dble(iexlexp)
    return
  end subroutine read_fix
  
  subroutine read_fixcoul
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m 
    
    implicit none
    
    iexlcoul=1
    exlcoul=dble(iexlcoul)
    return
  end subroutine read_fixcoul

  subroutine read_fixexch
    use params
    use discret
    use scf
    use solver 
    use commons8
    use input_m
    
    implicit none

    iexlexp=1
    exlexp=dble(iexlexp)
    return
  end subroutine read_fixexch


  subroutine read_fixorb
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none

    inzero=0
    do i=1,maxflags
       id(i)=0
       call inInt(id(i))
       if (id(i).gt.0) then
          inzero=inzero+1
          if (inzero.gt.norb) then
             write(iout6,*) 'Error: too many items - see User''s guide'
             stop 'read_fixorb'
          endif
          ifix(norb-id(i)+1)=1
       endif
    enddo
    
    ! if no orbitals are give fix all of them
    if (inzero.eq.0.or.inzero==norb) then
       iexlorb=1
       exlorb=dble(iexlorb)
       do i=1,norb
          ifix(i)=1
       enddo
    endif
    return
  end subroutine read_fixorb

  subroutine read_gauss
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    implicit none
    
    ! label: gauss
    
    call inFloat(z1atmass)
    if (z1atmass.le.0.0_PREC) then
       izz1=nint(z1)
       z1atmass=atweight(izz1)
    endif
    call inFloat(z2atmass)
    if (z2atmass.le.0.0_PREC) then
       izz2=nint(z2)
       z2atmass=atweight(izz2)
    endif
    if (z1atmass.le.0.0_PREC.and.z2atmass.le.0.0_PREC) then
       write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
       stop 'read_gauss'
    endif
    ifermi=2
    ipot=2
    return
  end subroutine read_gauss

  subroutine read_grid
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inInt_m
    use nmucalc_m
    use nnucalc_m
    
    implicit none
    
    if (icompLEnc.ne.3) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'
       stop 'read_grid'
    endif
    
    icompLEnc=icompLEnc+1

    tmp1=0.0_PREC
    tmp2=0.0_PREC
    call inInt(nni)
    call inFloat(tmp1)
    call inFloat(tmp2)
    if (tmp2.ne.0.0_PREC) then
       ! nnu nmu R_infty
       n=nint(tmp1)
       nmu(ngrids)=nmucalc(n)
       rgrid(ngrids)=tmp2
    else
       ! nni R_infty
       if (ngrids.ne.1) then
          write(*,*)'inputData: missing item in GRID card'
          stop 'inputData'
       endif
       rgrid(ngrids)=tmp1
       n=0
       nmu(1)=nmucalc(n)
    endif
    
    rinf=rgrid(ngrids)
    nmutot=nmu(ngrids)

    ! TODO get rid of multiple grids altogether
    do iorb=1,norb
       i1ng(iorb)=ngrids
       i2ng(iorb)=ngrids
    enddo
    
    nni=nnucalc(nni)
    if (nni.lt.7) then
       write(*,*) 'Error: too few grid points in nu variable: ',nni
       stop 'read_grid'
    endif
    if (nni.gt.maxnu) then
       write(*,*) 'Error: too many grid points in nu variable: ',nni,' is greater than ',maxnu
       stop 'read_grid'
    endif
    
    if (nmutot.lt.7) then
       write(*,*) 'Error: too few grid points in mu variable: ',nmutot
       stop 'read_grid'
    endif
    if (nmutot.gt.maxmu) then
       write(*,*) 'Error: too many grid points in mu variable:',nmutot,' is greater than ',maxmu
       stop 'read_grid'
    endif
    
    ! do iorb=1,norb
    !    nmumax(iorb)=nmutot
    ! enddo
    return
  end subroutine read_grid

  subroutine read_header
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
   
    implicit none
    
    !     label: title
    call inCardh(header)
    call inStr(clabel)
    
    if (clabel.ne.'title'.and.clabel.ne.'TITLE') then
       write(iout6,'(/,"Error: label TITLE is missing. ")')
       write(iout6,*) 
       stop 'read_header'
    endif
    
    icompLEnc=icompLEnc+1
    return
  end subroutine read_header
  
  subroutine read_homo
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m   
    
    implicit none
    ! label: homo --  when this label is encountered the g/u symmetry of
    !                 orbitals is strictly imposed; in all other respects
    !                 a homonuclear case is treated as a heteronuclear one
    
    ! ihomon=0 --  heteronuclear case
    ! ihomon=1 --  homonuclear case |z1-z2|<homolevl (set in common block data)
    ! ihomon=2 --  homonuclear case with forced symmetry
    
    ihomon=2
    return
  end subroutine read_homo

  subroutine read_inout
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inStr_m
    
    implicit none
    ! label: inout i32|i64|r128
    ! force format of input and output data
    lengthintin=lengthint
    lengthfpin=lengthfp
    inpform=0
    ioutform=0
    inout32=0
    inout64=0
    inout128=0
    call inStr(clabel)
    if     (clabel.eq.'i32') then 
       lengthintin=4
       lengthfpin=8
       inpform=0
    elseif (clabel.eq.'i32f') then 
       lengthintin=4
       lengthfpin=8
       inpform=1
       formfp=formfp64
    elseif (clabel.eq.'i64') then
       lengthintin=8
       lengthfpin=8
       inpform=0
    elseif (clabel.eq.'i64f') then
       lengthintin=8
       lengthfpin=8
       inpform=1
       formfp=formfp64
    elseif (clabel.eq.'r128') then
       lengthintin=8
       lengthfpin=16
       inpform=0
    elseif (clabel.eq.'r128f') then
       lengthintin=8
       lengthfpin=16
       inpform=1
       formfp=formfp128
    else
       write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
       stop 'read_inout'
    endif
    
    call inStr(clabel)
    if     (clabel.eq.'i32') then 
       inout32=1
       ioutform=0
    elseif (clabel.eq.'i32f') then
       inout32=2
       ioutform=1
    elseif (clabel.eq.'i64') then
       inout64=1
       ioutform=0
    elseif (clabel.eq.'i64f') then
       inout64=2
       ioutform=1
       formfp=formfp64
    elseif (clabel.eq.'r128') then
       inout128=1
       ioutform=0
    elseif (clabel.eq.'r128f') then
       inout128=2
       formfp=formfp128
       ioutform=1
       if (lengthfp.ne.16) then
          write(*,*) 'Error: present build of the program does not support quadruple precision'
          stop 'read_inout'
       endif
    else
       write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
       stop 'read_inout'
    endif
  end subroutine read_inout

  subroutine read_interp
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       iord_nu_orb=itmp+1
       call inInt(iord_mu_orb)
       iord_mu_orb=iord_mu_orb+1
       
       call inInt(iord_nu_coul)
       iord_nu_coul=iord_nu_coul+1
       call inInt(iord_mu_coul)
       iord_mu_coul=iord_mu_coul+1
       
       call inInt(iord_nu_exch)
       iord_nu_exch=iord_nu_exch+1
       call inInt(iord_mu_exch)
       iord_mu_exch=iord_mu_exch+1
    endif
    if (((iord_nu_orb.ne.3).and.(iord_nu_orb.ne.5).and.(iord_nu_orb.ne.7).and.&
         (iord_nu_orb.ne.9)).or.&
         ((iord_mu_orb.ne.3).and.(iord_mu_orb.ne.5).and.(iord_mu_orb.ne.7).and.&
         (iord_mu_orb.ne.9)).or.&
         ((iord_nu_coul.ne.3).and.(iord_nu_coul.ne.5).and.(iord_nu_coul.ne.7).and.&
         (iord_nu_coul.ne.9)).or.&
         ((iord_mu_coul.ne.3).and.(iord_mu_coul.ne.5).and.(iord_mu_coul.ne.7).and.&
         (iord_mu_coul.ne.9)).or.&
         ((iord_nu_exch.ne.3).and.(iord_nu_exch.ne.5).and.(iord_nu_exch.ne.7).and.&
         (iord_nu_exch.ne.9)).or.&
         ((iord_mu_exch.ne.3).and.(iord_mu_exch.ne.5).and.(iord_mu_exch.ne.7).and.&
         (iord_mu_exch.ne.9))) &
         then
       write(iout6,*) 'Error: incorrect order of interpolation polynomial - allowed values are 2, 4, 6 or 8'
       stop 'read_interp'
    endif
    iinterp=1
    return
  end subroutine read_interp

  subroutine read_lagra
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
   
    implicit none
    ! label: lagra
    ! scaling and damping factors for off-diagonal Lagrange multipliers
    
    ilagra=1
    do i=1,maxorb*(maxorb-1)/2
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          id1=itmp1
       else
          return
       endif
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          id2=itmp1
       else
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_lagra'
       endif
       nlmf(norb-id1+1,norb-id2+1)=1
       nlmf(norb-id2+1,norb-id1+1)=1
    enddo
    return
  end subroutine read_lagra
  
  subroutine read_lcao
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use nmucalc_m
    use nnucalc_m
   
    implicit none

    if (icompLEnc.ne.5) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'
       stop 'read_lcao'
    endif
    !         icompLEnc=icompLEnc+1
    icompLAdd=icompLAdd+1

    inclorb=1
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       inclorb=itmp
    endif
    lcaoIncl=.true.

    ! skip lcao info when 'orbpot old'
    if (ini==5) then
       do iorb=1,norb
          call inCard
       enddo
       return
    endif
    ! FIXME3 
    if (ini.eq.1.and.(inclorb.ne.1.and.inclorb.ne.2)) then
       write(iout6,*) 'Error: incompatible parameters - no LCAO data present'
       stop 'read_lcao'
    endif
    
    if (imethod.eq.2) then
       iexlcoul=1
       exlcoul=dble(iexlcoul)
       iexlexp=1
       exlexp=dble(iexlexp)
       if (ini.ne.1) ini=4 
       facmul=10.0_PREC**(15)
    endif
    
    if (inclorb.eq.2) then
       
       ! 'hydrogen' initialization with screening
       
       do iorb=1,norb
          call inCard
          call inFloat(co1(iorb))
          call inInt(mgx(1,iorb))
          call inInt(mgx(2,iorb))
          call inFloat(eza1(iorb))
          eza1(iorb)=z1-eza1(iorb)
          call inFloat(co2(iorb))
          call inInt(mgx(4,iorb))
          call inInt(mgx(5,iorb))
          call inFloat(eza2(iorb))
          eza2(iorb)=z2-eza2(iorb)
          call inInt(inhyd(iorb))
          call inInt(maxsororb(iorb))
          if (maxsororb(iorb).eq.inpiexit) then
             maxsororb(iorb)=maxsor2
             maxsorpot(iorb)=maxsor3
          else
             call inInt(maxsorpot(iorb))
             if (maxsorpot(iorb).eq.inpiexit) then
                maxsorpot(iorb)=maxsor3  
             endif
          endif
       enddo
       
    elseif (inclorb.eq.1) then

       !  'hydrogen' initialization without screening
       do iorb=1,norb
          call inCard
          call inFloat(co1(iorb))
          call inInt(mgx(1,iorb))
          call inInt(mgx(2,iorb))
          call inFloat(eza1(iorb))
          eza1(iorb)=eza1(iorb)
          call inFloat(co2(iorb))
          call inInt(mgx(4,iorb))
          call inInt(mgx(5,iorb))
          call inFloat(eza2(iorb))
          eza2(iorb)=eza2(iorb)
          call inInt(inhyd(iorb))
          call inInt(maxsororb(iorb))
          if (maxsororb(iorb).eq.inpiexit) then
             maxsororb(iorb)=maxsor2
             maxsorpot(iorb)=maxsor3
          else
             call inInt(maxsorpot(iorb))
             if (maxsorpot(iorb).eq.inpiexit) then
                maxsorpot(iorb)=maxsor3
             endif
          endif
       enddo
    endif

    ! let's make sure that orbitals are to be initialized in case of inconsistent lcao data
    do i=1,norb
       inhyd(i)=1
    enddo
    
    do i=1,norb
       inhydlcao(i)=inhyd(i)
    enddo

    ! normalize mixing coefficients
    do iorb=1,norb
       co12=abs(co1(iorb))+abs(co2(iorb))
       co1(iorb)=co1(iorb)/co12
       co2(iorb)=co2(iorb)/co12
    enddo

    ! ini=5 enables to retrieve from complete/uncomplete dump file
    ! if only one (the highest= the first) orbital is missing
    ! the rest is retreived and the first is initialized as 
    ! a hydrogenic one (ini is changed from 5 to 4)

    return
  end subroutine read_lcao

  subroutine read_maxsor
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: maxsor
    
    do iorb=1,norb
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          call inInt(itmp2)
          if (itmp2.ne.inpiexit) then
             maxsororb(itmp1)=itmp2
          else
             write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
             stop 'read_maxsor'
          endif
       else
          goto 100
       endif
    enddo
100 continue
    return
  end subroutine read_maxsor

  subroutine read_mcsor
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m 
    use inInt_m
    use nmucalc_m
    use nnucalc_m
    
    implicit none
    ! label: mcsor
    !     maxsor2 - maximal number of (mc)sor iterations during relaxation
    !               of every orbital in an SCF cycle
    !   maxsorpot - maximal number of (mc)sor iterations during relaxation
    !               of every potential in an SCF cycle
    !     ipoiss  - choice of Poisson's equation solving method
    !        = 1    sor (modified LPS's code)
    !        = 2    mcsor - variables are relaxed in 5 sweeps each relaxing
    !               every 5th variable 
    
    ipoiss=2
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsor2 = itmp
       do i=1,norb
          maxsororb(i)=maxsor2
       enddo
       
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsor3 = itmp 
          do i=1,norb
             maxsorpot(i)=maxsor3
          enddo
          
       endif
    endif
    nni=nnucalc(nni)
    nmu(1)=nmucalc(nmu(1))
    nmutot=nmu(1)
    
    return
  end subroutine read_mcsor
  
  subroutine read_method
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inStr_m
    implicit none

    ! imethod=1 -- hf
    ! imethod=2 -- oed
    ! imethod=3 -- hfs
    ! imethod=4 -- dft
    ! imethod=5 -- scmc
    
    call inStr(char8)
    
    ihit=0
    do i=1,nmethods
       if (char8.eq.cmethod(i)) then
          ihit=ihit+1
          imethod=i
       endif
    enddo
    
    if (ihit.eq.0) then
       write(iout6,*) 'Error: incorrect method - allowed values are HF, HFS, DFT, OED or SCMS'
       stop 'read_method'
    endif
    ! HFS method is equivalent to DFT one with the LDA echange potential and the optimum
    ! value of alpha
    
    ! for HFS method set alpha to a (hopefully) reasonable value based
    ! on atomic optimum values due to Schwarz (see blk-data.f90)
    
    ! see winding up section of this routine
    
    if (imethod.eq.3) then
       islat =1
       idftex=1
    endif
    
    ! for DFT method the exchange potential is set to LDA and alpha is set (by default) to
    ! 2/3
    
    if (imethod.eq.4) then
       ! imethod=3
       islat =1
       idftex=1
       alphaf=two/three
    endif
    
    ! alpha parameter of the density-functional theory is calculated according to the
    ! Self-Consistent Multiplicative Constant method (see V.V.Karasiev and E.V.Ludenia,
    ! Self-consistent multiplicative constant method for the exchange energy in density
    ! functional theory, Phys. Rev. A 65 (2002) 062510
    
    if (imethod.eq.5) then
       ! imethod=3
       islat =1
       idftex=1
       iscmc=1   
       alphaf=two/three
    endif

    ! extra consistancy checks 
    if (imethod.eq.3.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as HFS the second parameter must be set to 3'
       stop 'inputData'
    endif
    
    if (imethod.eq.4.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as DFT the second parameter must be set to 3'
       stop 'inputData'
    endif
    
    if (imethod.eq.5.and.iform.ne.3) then
       write(iout6,*) 'Error: if method is chosen as SCMC  the second parameter must be set to 3'
       stop 'inputData'
    endif
    return
  end subroutine read_method

  subroutine read_multipol
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inInt_m
    
    implicit none
    ! label: multipol 
    ! define number of multipole moments to be calculated
    !     facmul  -  if facmul>0  multipole moment expansion
    !		 coefficients are recalculated	every time denmax,
    !		 i.e. maximum error in orbital energy, changes by
    !		 the factor facmul
    !		 if facmul<0 these coefficients are not calculated
    !      mpole  - number of multipole expansion coefficients
    
    call inFloat(facmul)
    if (facmul.lt.0.0_PREC) facmul=ten**(15)
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       mpole=itmp
       if (mpole.lt.2.or.mpole.gt.8) then
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_multipol'
       endif
    endif
    return
  end subroutine read_multipol

  subroutine read_nonortho
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    nonortho=1
    do i=1,maxorb*(maxorb-1)/2
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          id1=itmp1
       else
          return
       endif
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          id2=itmp1
       else
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_nonortho'
       endif
       nonorthog(norb-id1+1,norb-id2+1)=1
       nonorthog(norb-id2+1,norb-id1+1)=1
    enddo
    return
  end subroutine read_nonortho

  subroutine read_nuclei
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inStr_m
    
    implicit none
    if (icompLEnc.ne.1) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'       
       stop 'read_nuclei'
    endif
    
    icompLEnc=icompLEnc+1
    ! inuclei=1
    call inFloat(z1)
    call inFloat(z2)
    call inFloat(r)
    call inStr(char8)
    izz1=nint(z1)
    izz2=nint(z2)
    ! conversion factor due to Cohen and Taylor (1986), The 1986 Adjustment of the
    ! Fundamental Physical Constants
    if (char8.eq.angstrom) r=r/bohr2ang
    r2=r/two
    if ((abs(z1)+abs(z2)).eq.0.0_PREC.or.r.eq.0.0_PREC) then
       write(iout6,*) 'Error: incomplete input'
       stop 'read_nuclei'
    endif
    return
  end subroutine read_nuclei

  subroutine read_omega
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    
    implicit none
    ! label: omega
    !  ovforb   -  overelaxation parameter for orbitals and grids
    !  ovfcoul  -  overelaxation parameter for coulomb pot. and grids
    !  ovfexch  -  overelaxation parameter for exchange pot. and grids

    omegaIncl=.true.    
    call inFloat(ftmp1)
    call inFloat(ftmp2)
    if (ftmp1.ne.0.0_PREC.or.ftmp2.ne.0.0_PREC) then
       !         if (ftmp1.gt.0.0_PREC.and.ftmp2.ne.0.0_PREC) then
       if (ftmp1.ne.0.0_PREC) ovforb(1)=ftmp1
       if (ftmp2.ne.0.0_PREC) ovfcoul(1)=ftmp2
    else
       !           old format
       call inCard
       call inFloat(ftmp)
       if (ftmp.gt.0.0_PREC) then
          ovforb(1)=ftmp
       else
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_omega'
       endif
       
       do ig=2,ngrids
          ovforb(ig) = ovforb(1)
       enddo
       
       call inFloat(ftmp)
       if (ftmp.gt.0.0_PREC) then
          ovforb(2)=ftmp
          call inFloat(ftmp)
          if (ftmp.gt.0.0_PREC) then
             ovforb(3)=ftmp
          endif
       endif
       
       call inCard
       call inFloat(ftmp)
       ovfcoul(1)=ftmp
       if (ftmp.gt.0.0_PREC) then
          do ig=2,ngrids
             ovfcoul(ig) = ovfcoul(1)
          enddo
          
          call inFloat(ftmp)
          if (ftmp.gt.0.0_PREC) then
             ovfcoul(2)=ftmp
             call inFloat(ftmp)
             if (ftmp.gt.0.0_PREC) then
                ovfcoul(3)=ftmp
             endif
          endif
       else
          ! omega values cannot be set automatically for subgrids
          if (ngrids.gt.1) then
             write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
             stop 'read_omega'
          endif
       endif
       do ig=1,ngrids
          ovfexch(ig) = ovfcoul(ig)
       enddo
    endif
    return
  end subroutine read_omega

  subroutine read_omegaopt
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inInt_m
    
    implicit none
    !  label: omegaopt
    !  1  -  optimal overelaxation parameter for orbitals 
    !  2  -  optimal overelaxation parameter for Coulomb and exchange potentials
    !  3  -  optimal overelaxation parameters for orbitals and potentials

    !  1 - rather conservative (safe) omega values (equivalent to
    !  automatic parameters selection available via omega card (default)
    !  2 - near-optimal omega values for cases when good initial approximations 
    !      are available or when fixed orbitals/potentials calculations are performed

    !  two optional parameters for scaling the chosen omega values for orbitals and 
    !      potentials
    
    call inInt(iomega)
    if (iomega.eq.0) then
       iomega=1
    endif
    
    if (iomega.eq.1) then
       ovforb(1)=-1.0
       ovfcoul(1)=-1.0
    endif
    
    call inFloat(ftmp)
    if (ftmp.gt.0.0_PREC) then
       omegasfOrb=ftmp
       call inFloat(ftmp)
       if (ftmp.gt.0.0_PREC) then
          omegasfPot=ftmp
       endif
    endif
    return
  end subroutine read_omegaopt

  subroutine read_orbpot
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    
    if (icompLEnc.ne.4) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'
       stop 'read_orbpot'
    endif
    icompLEnc=icompLEnc+1
    
    call inStr(clabel)
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       inclorb=itmp
    endif
    if     (clabel.eq.'hydrogen') then
       ini=1
       icompLExp=1
    elseif (clabel.eq.'gauss') then
       ini=2
    elseif (clabel.eq.'gauss-c') then
       ini=3
    elseif (clabel.eq.'old') then
       ini=5
    elseif (clabel.eq.'noexch') then
       ini=6
    elseif (clabel.eq.'qrhf') then
       ini=11
    elseif (clabel.eq.'lda') then
       ini=12
       ldaIncl=.true.
    elseif (clabel.eq.'molcas') then
       ini=22
    elseif (clabel.eq.'nodat') then
       ini=55
    else
       write(iout6,*) 'Error: incorrect source of orbitals and potentials - try HYDROGEN, GAUSS, GAUSS-C, MOLCAS, OLD, NOEXCH, NODAT'
       stop 'read_orbpot'
    endif
    !        initial and orbpot labels share some code
    return
  end subroutine read_orbpot

  subroutine read_order
    
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    ! label: order
    
    ! iorder  - type of (mc)sor sweeps (ordering) for each grid
    !      = 1   - natural column-wise ordering for the mesh points
    !	   = 2   - default 'middle' type of sweep 
    !	   = 3   - natural row-wise
    !	   = 4   - reversed natural (column-wise) ordering (see mesh for details)
    
    !      = 10+iorder - in case of near-degenerate orbitals, i.e. when
    !          performing calculations for homonuclear molecules without
    !          inforced symmetry (no homo label), especially when external
    !          electric field is applied, one can improve convergence by
    !          changing the direction of the SOR sweeps (forward/backward
    !          sweeps for even/odd SCF)
    
    call inInt(iorder(ngrids))
    
    if (iorder(ngrids).gt.10) then
       iorder(ngrids)=iorder(ngrids)-10
       ialtsweeps=1
    endif
    
    if (iorder(ngrids).gt.4) then
       write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
       stop 'read_order'
    endif
    return
  end subroutine read_order

  ! subroutine read_out4dd
  !   use params
  !   use discret
  !   use scf
  !   use solver
  !   use commons8
  !   use input_m
    
  !   implicit none
  !   iout4dd=1
  !   return
  ! end subroutine read_out4dd

  subroutine read_potgsz
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    ! Green, Sellin, Zachor model potential 
    ipot=5
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_potgsz'
    endif
    return
  end subroutine read_potgsz

  subroutine read_potgszg
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    ! Green, Sellin, Zachor model potential + finite Gauss nucleus model 
    ipot=55
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_potgszg'
    endif

    return
  end subroutine read_potgszg

  subroutine read_potsap
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    ! Susi Lehtola: Superposition of Atomic Potentials
    ipot=66
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_potsap'
    endif
    
    return
  end subroutine read_potsap

  subroutine read_poth3
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inInt_m
    
    implicit none
    
    call inInt(mpot)
    call inFloat(apot)
    call inFloat(v0pot)
    ifermi=3
    ipot=3
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_poth3'
    endif
    return
  end subroutine read_poth3

  subroutine read_potharm
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    ! harmonic potential 
    ipot=10
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_potharm'
    endif
    return
  end subroutine read_potharm

  subroutine read_pothook
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    
    implicit none
    ! harmonic potential (Hook's atom, harmonium) 
    call inFloat(hook)
    ipot=9
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_pothook'
    endif
    return
  end subroutine read_pothook

  subroutine read_potkh
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    use inInt_m
    
    implicit none
    ! label: potkh
    
    apot=1.0_PREC
    v0pot=1.0_PREC
    nsimp=1000
    call inInt(mpot)
    call inFloat(epspot)
    call inFloat(ompot)
    call inFloat(tmp1)
    if (tmp1.ne.0.0_PREC) then
       apot=tmp1
       call inFloat(tmp2)
       if (tmp2.ne.0.0_PREC) then
          v0pot=tmp2
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             nsimp=itmp
          endif
       endif
    endif
    ifermi=4
    ipot=4
    if (imethod/=2) then
       write(iout6,'("Error: this potential cannot be used with ",a3," method! Try OED instead.")') cmethod(imethod)
       stop 'read_potkh'
    endif
    return
  end subroutine read_potkh

  subroutine read_print
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    
    implicit none
    ! label: print
    ! label followed by up to maxflags integers (1..999) specifying additional printouts
    ! FIXME
    ! list of used iprint flags:
    !   prepArray: 10 - check i3bxx arrays (separate files case)
    !   ortho: 20 - print overlap integrals
    !   ortho: 21 - print overlap integrals only for pairs of orbitals
    !               having non-zero off-diagonal Lagrange multipliers
    !   preSCF: 30 - check orthogonality
    !   ortho: 35 - force normalization of orbitals upon their orthogonalization         
    !   norm:  36 - check normalization
    !   mesh: 40,41     
    !   EabHF:  46  - print off-diagonal Lagrange multipliers 
    !   EabDFT: 48  - print off-diagonal Lagrange multipliers 
    !   doSCF: 50 - print orbital contributions to total energy
    !   setci: 55 - symmetry of orbitals
    !   EaHF: 60 - kinetic energy, nuclear energy  one-electron energy
    !   EaHF: 62,64,65,66 - two-electron contributions
    !   EaDFT: 67,68 one- and two-electron contributions
    !   EaDFT: 69 orbitals
    !   EaDFT: 71 wgt1, wgt2
    !   EaDFT: 71 f0
    !   EaDFT: 72 wk1, wk2
    !   EtotalDFT: 75 - test n2f and nfnf routines
    !   EtotalDFT: 77 - woneel,wndc,wex,etotal, woneel,wndc,wex,evt
    !   EtotalHF: 78, 79 - woneel,wndc,wex,etotal, woneel,wndc,wex,evt
    !   EtotalHForb: 80, 81 - woneel,wndc,wex,etotal, woneel,wndc,wex,evt
    !   scmc: 82,83 individual exchange contributions
    !   scmc: 84 - ehfex,edftex,alphaf
    !   scmc: 85 - wdcoul,wex1,wex2
    !   fock:  91 - one-electron contribution
    !   fock:  92 - coo0,coo1,coo2,engo(ideng)
    !   fock:  93 - Coulomb potential'
    !   fock:  95 - checking (E-V) values for a given orbital
    !   printResults: 100 - find last maximum in mu variable
    !   radialden: 110 - total radial density along the internuclear
    !                    axis from the centre A to minus infinity (-rinf)
    !   radialden: 111 - total radial density along the internuclear
    !                    axis from the centre B to plus infinity (+rinf)
    !   rfun: 114 - r values at grid points (to be implemented)
    !   rfun: 115 - input orbitals 
    !   rfun: 116 - input Coulomb potentials
    !   rfun: 117 - input exchange potentials
    !   wtdisknat: 121 - output orbitals 
    !   wtdisknat: 122 - output Coulomb potentials
    !   wtdisknat: 123 - output exchange potentials
    !   interpolq: 125: iadint2, cint2, iadint4, cint4
    !   interpolq: 126 - iadint2, cint2, iadint4, cint4, etc
    !   zz1: 131 - i,j,ri,zz1
    !   zz2: 132 - i,j,ri,zz2
    !   printResults: 140 - calculate the total electronic density at A
    !   printResults: 150 - calculate the derivative of the orbital
    !                       density with respect to z at A
    !   printResults: 152 - calculate multipole moments errors due to
    !                       orbital norms not being equal to 1
    !   coulMom: 160 - d1, d3, d5
    !            161 - multipole moments for a given orbital      
    !   EaDFT: 459
    !   exchMom: 165 - dome
    !            166 - excdi(ipc),excqu(ipc),excoc(ipc),exche(ipc) 
    !   initEXWeights: 170, 172 - ,idiv,gec(kk),gec(kk+isu2),div(i), etc
    !   initCBlocks: 180 - subgrids
    !   printResults:  191 - check orthogonality of orbitals
    !   printResults:  192 - calculate Euclidean norms of (T+V(n)+V-E)|i>
    !   printResults:  193 - check multipole expansion
    !   propet2: 214,216
    !   initHF: 220,221,222
    !   fermi: 230,231
    !   prepGauss: 553, 554, 555
    !   prepGaussCust: 563, 564, 565
    !   prttot: 590 print etotFN
    
    inzero=1
    do i=1,maxflags
       id(i)=0
       call inInt(id(i))
       if (id(i).gt.0) then
          inzero=inzero+1
          iprint(id(i))=1
       endif
    enddo
    return
  end subroutine read_print

  subroutine read_prtevery
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: prtevery
    ! customize printouts of two-dimentional arrays
    ! incrni  -  parameters of pmtx routine used to print two-dimensional arrays 
    ! incrmi	 row-wise. Every incrni row and every incrmi column is printed
    
    call inInt(incrni)
    call inInt(incrmu)
    return
  end subroutine read_prtevery

  ! subroutine read_reinit
  !   use params
  !   use discret
  !   use commons8
  !   use input_m
  !   use inInt_m
    
  !   implicit none

  !   ireinit=1
  !   inzero=0
  !   do i=1,norb
  !      id(i)=0
  !      call inInt(id(i))
  !      if (id(i).gt.0) then
  !         inzero=inzero+1
  !         if (inzero.gt.norb) then
  !            write(iout6,*) 'Error: too many items - see User''s guide'
  !            stop 'read_reinit'
  !         endif
  !         reinit(norb-id(i)+1)=1
  !      endif
  !   enddo
    
  !   return
  ! end subroutine read_reinit


  ! subroutine read_rydberg
  !   use params
  !   use discret
  !   use scf
  !   use solver
  !   use commons8
  !   use input_m
    
  !   implicit none
  !   irydberg=1
  !   return
  ! end subroutine read_rydberg
  
  subroutine read_scfexch
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    iexlexp=0
    exlexp=dble(iexlexp)
    return
  end subroutine read_scfexch

  subroutine read_scf
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: scf
    !        
    !   maxscf  - max. no. of scf iterations
    !   nobckup - no. of iterations between backups 
    !   ienterm - if max. error in orbital energy is less than 
    !             1/10**ienterm than scf iterations are terminated
    !   inoterm - if max. error in orbital norm is less than 
    !             1/10**inoterm than scf iterations are terminated
    !   iprtlev - level of output during scf process
    !             conv. rate, orbital energy, normalization of the worst
    !             converged orbital is printed every nobckup iterations
    !             if iprtlev=3 or every scf iteration if iprtlev=2 (default).
    !             if iprtlev=1 then conv. rate, orbital energy and normalization 
    !             of every orbital is printed in every scf iteration;
    !             in any case total energy is printed every nobckup iterations;
    
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxscf = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          nobckup = itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             ienterm = itmp
             call inInt(itmp)
             if (itmp.ne.inpiexit) then
                inoterm = itmp
                call inInt(itmp)
                if (itmp.ne.inpiexit) then
                   iprtlev=itmp
                endif
             endif
          endif
       endif
    endif
    return
  end subroutine read_scf

  subroutine read_scforder
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: scforder
    do iorb=1,norb
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          iscforder(iorb)=itmp
       else
          write(iout6,*) 'Error: missing or incorrect entry - see User''s guide'
          stop 'read_scforder'
       endif
    enddo
    return
  end subroutine read_scforder

  ! subroutine read_sgga
  !   use params
  !   use discret
  !   use scf
  !   use solver
  !   use commons8
  !   use input_m
    
  !   implicit none
  !   sgga=.true.
  !   return
  ! end subroutine read_sgga
  
  subroutine read_sor
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    
    implicit none
    ! label: sor
    !     maxsor2 - maximal number of (mc)sor iterations during relaxation
    !               of every orbital in an SCF cycle
    !   maxsorpot - maximal number of (mc)sor iterations during relaxation
    !               of every potential in an SCF cycle
    !     ipoiss  - choice of Poisson's equation solving method
    !        = 1    sor (modified LPS's code)
    !        = 2    mcsor - variables are relaxed in 5 sweeps each relaxing
    !               every 5th variable 
    
    ipoiss=1
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsor2 = itmp
       do i=1,norb
          maxsororb(i)=maxsor2
       enddo
       
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsor3 = itmp 
          do i=1,norb
             maxsorpot(i)=maxsor3
          enddo
       endif
    endif
    return
  end subroutine read_sor

  subroutine read_sormcsor
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inInt_m
    use nmucalc_m
    use nnucalc_m
    
    implicit none
    !        label: sormcsor
    !     maxsor2 - maximal number of (mc)sor iterations during relaxation
    !               of every orbital in an SCF cycle
    !   maxsorpot - maximal number of (mc)sor iterations during relaxation
    !               of every potential in an SCF cycle
    !     ipoiss  - choice of Poisson's equation solving method
    !        = 1    sor (modified LPS's code)
    !        = 2    mcsor - variables are relaxed in 5 sweeps each relaxing
    !               every 5th variable 
    !        = 3    SOR for orbitals and MCSOR for potentials 
    ipoiss=3
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsor2 = itmp
       do i=1,norb
          maxsororb(i)=maxsor2
       enddo
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsor3 = itmp 
          do i=1,norb
             maxsorpot(i)=maxsor3
          enddo
       endif
    endif
    nni=nnucalc(nni)
    nmu(1)=nmucalc(nmu(1))
    nmutot=nmu(1)
    return
  end subroutine read_sormcsor

  subroutine read_stop(ni_t,mu_t,no_t,nons_t)
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    
    implicit none
    integer :: ni_t,mu_t,no_t,nons_t
    if (icompLExp.ne.0.and.icompLAdd.eq.0) then
       write(iout6,*) 'Error: incorrect order of labels: try TITLE .. NUCLEI .. CONFIG .. GRID .. ORBPOT ...'
       stop 'read_stop'
    endif
    ! normalize mixing coefficients
    do iorb=1,norb
       co12=abs(co1(iorb))+abs(co2(iorb))
       co1(iorb)=co1(iorb)/co12
       co2(iorb)=co2(iorb)/co12
    enddo
    
    iform_t = iform
    ni_t    = nni
    mu_t    = nmutot
    no_t=norb
    nons_t=0
    do iorb=1,norb
       if (mgx(3,iorb).ne.0) nons_t=nons_t+1
    enddo
    
    write(*,*) '  ... end of input data  ...'
    write(*,*)'        '
    return
  end subroutine read_stop

  subroutine read_xalpha
    use params
    use discret
    use scf
    use solver
    use commons8
    use input_m
    use inFloat_m
    
    implicit none
    ! alphaf - if alphaf is not zero a DFT energy functional is used
    
    ftmp=0.0_PREC
    call inFloat(ftmp)
    if (ftmp.ne.0.0_PREC) then
       alphaf = ftmp
    endif
    return
  end subroutine read_xalpha

  
end module readLables_m

