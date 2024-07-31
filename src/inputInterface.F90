! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2023  Jacek Kobus 

module inputInterface
  use params
  implicit none
  integer (KIND=IPREC) :: i40,i80  
  parameter (i40=40, i80=80)
  integer (KIND=IPREC) :: jrec,jump
  integer (KIND=IPREC),dimension(i40) :: istrt,inumb
  character*1, dimension(i80) :: ia
  character*80 :: header

contains

  ! ### inCardTitle ###
  !
  !     This routine reads (and echoes) a single data card and
  !
  !       - changes upper case letters into a lower case ones
  !       - changes tabs (if any) into spaces
  !       - removes leading spaces
  !
  !     Then it scans the line for nonspace fields. The number of fields
  !     is stored in jump, the starting point of a given field in istrt(i)
  !     and the number of characters in that field in inumb(i). If the
  !     routine finds an exclamation mark or a hash then anything on the
  !     same line that follows these characters is treated as a comment.
  !
  subroutine inCardTitle(title)
    use params
    implicit none

    logical, parameter :: data_exist=.true.
    character*1 :: blank,exmark,hash
    character*1, dimension(80) :: iatmp
    character*1, dimension(80) :: title
    
    integer (KIND=IPREC) :: i,ido,ike,ile,isw,it
    
    data blank/' '/,exmark/'!'/,hash/'#'/

    do i=1,i80
       ia(i)=' '
    enddo

    jump = 0
    jrec = 0
    isw = 0
    do while (data_exist)
       read(iinp5,100,err=40,end=50) ia
       do i=1,i80
          iatmp(i)=ia(i)
       enddo

       ! check for blanks at left and remove them, fill with blanks at the end
       it=i80
       call lftpos(iatmp,it)

       ! check for tabs, and replace them by blanks
       call tabchk(iatmp,it)

       ! only the label is set into lower case 
       it=6
       call lowcase(iatmp,it)

       do i=1,i80
          ia(i)=iatmp(i)
       enddo
       
       ! check for an exclamation mark or a hash
       do ike=1,i80
          ile=ike
          if (ia(ile).eq.exmark.or.ia(ile).eq.hash) exit
       enddo
       
       if(ile.eq.1) cycle
       !  echo the slightly modified input data if not empty
       if (ia(1).ne.blank) write(iout6,101) ia
       
       if(ia(ile).eq.exmark.or.ia(ile).eq.hash) ile=ile-1
       
       do i = 1,ile
          if (ia(i).eq.blank) then
             isw=0
             cycle
          else
             if (isw.le.0) then
                jump = jump +1
                istrt(jump) = i
                inumb(jump) = 0
                isw=1
                inumb(jump) = inumb(jump) + 1
             else
                inumb(jump) = inumb(jump) + 1
             endif
          endif
       enddo
       
       if (jump.eq.0) then
          cycle
       else
          exit
       endif
    enddo

    ! store the title message in a separate array
    do i=7,i80
       title(i-6)=iatmp(i)
    enddo
    do i=75,i80
       title(i)=blank
    enddo
    return

40  write(iout6,140)
140 format(//1x,'Error encountered when reading input data.'//)
    stop 'inCardTitle'

50  write(iout6,150)
150 format(//1x,'Unexpected end of input data encountered. Is STOP label missing?'//)
    stop 'inCardTitle'
    
100 format(80a1)
101 format(2x,80a1)
  end subroutine inCardTitle
  
  ! ### inStr ###
  !
  !    Examines the contents of IA and extracts a character string of up
  !    to 8 characters. This string is stored in IBUF. The remaining
  !    non-blank characters (if any) are ignored.
  !
  subroutine inStr(guf)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n,nstrt
    character*1 iblnk
    character*8 :: ibufr,guf

    character*1, dimension(8) ::  ibuf
    character*1, dimension(60) :: iall

    equivalence (ibuf(1),ibufr)

    data iblnk/' '/
    data iall/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',&
         'u', 'v', 'w', 'x', 'y', 'z', '-', '-', '-', ' ',&
         'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J','K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',&
         'U', 'V', 'W', 'X', 'Y', 'Z', '-', '-', '-', ' '/

    do i=1,8
       ibuf(i)=iblnk
    enddo
    jrec = jrec + 1
    if(jrec .gt. jump) goto 11
    n = inumb(jrec)
    nstrt = istrt(jrec)
    if (n.gt.8) n=8
    do i = 1,n
       ibuf(i) = ia(nstrt)
       nstrt = nstrt + 1
    enddo

11  guf=ibufr

  end subroutine inStr

  ! ### inStr4lxc ###
  !
  subroutine inStr4lxc(guf)
    use params
    implicit none

    integer (KIND=IPREC) :: i,maxibuf,n,nstrt
    parameter (maxibuf=30)
    character*1 iblnk
    character*30 :: ibufr,guf

    character*1, dimension(maxibuf) :: ibuf
    character*1, dimension(60) :: iall

    equivalence (ibuf(1),ibufr)

    data iblnk/' '/
    data iall/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',&
         'u', 'v', 'w', 'x', 'y', 'z', '-', '-', '-', ' ',&
         'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J','K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',&
         'U', 'V', 'W', 'X', 'Y', 'Z', '-', '-', '-', ' '/

    do i=1,maxibuf
       ibuf(i)=iblnk
    enddo
    jrec = jrec + 1
    if(jrec .gt. jump) goto 11
    n = inumb(jrec)
    nstrt = istrt(jrec)
    if (n.gt.maxibuf) n=maxibuf
    do i = 1,n
       ibuf(i) = ia(nstrt)
       nstrt = nstrt + 1
    enddo

11  guf=ibufr

  end subroutine inStr4lxc

  ! ### inCard ###
  !
  !     This routine reads (and echoes) a single data card and
  !
  !       - changes upper case letters into a lower case ones
  !       - changes tabs (if any) into spaces
  !       - removes leading spaces
  !
  !     Then it scans the line for nonspace fields. The number of fields
  !     is stored in jump, the starting point of a given field in istrt(i)
  !     and the number of characters in that field in inumb(i). If the
  !     routine finds an exclamation mark or a hash then anything on the
  !     same line that follows these characters is treated as a comment.
  !
  subroutine inCard
    use params
    implicit none

    logical, parameter :: data_exist=.true.
    integer, parameter :: maxEmptyLines=50
    character*1 :: blank,exmark,hash
    character*1, dimension(80) :: iatmp

    integer (KIND=IPREC) :: i,ido,ike,ile,isw,it

    data blank/' '/,exmark/'!'/,hash/'#'/

    do i=1,i80
       ia(i)=' '
    enddo
    !     read a line
    jump = 0
    jrec = 0
    isw = 0
    do while (data_exist)
       !do ido=0,maxEmptyLines
       read(iinp5,100,err=40,end=50) ia
       
       do i=1,i80
          iatmp(i)=ia(i)
       enddo
       
       !     set into lower case
       !     FC3, ifort: error while passing i80 to a subroutine
       it=i80
       call lowcase(iatmp,it)
       
       !     check for tabs, and replace them by blanks
       call tabchk(iatmp,it)
       
       !     check for blanks at left and remove them, fill with blanks at the end
       call lftpos(iatmp,it)
       
       do i=1,i80
          ia(i)=iatmp(i)
       enddo
       
       !     check for an exclamation mark or a hash
       do ike=1,i80
          ile=ike
          if (ia(ile).eq.exmark.or.ia(ile).eq.hash) exit
       enddo
       
       if(ile.eq.1) cycle
       !  echo the slightly modified input data if not empty
       if (ia(1).ne.blank) write(iout6,101) ia
       
       if(ia(ile).eq.exmark.or.ia(ile).eq.hash) ile=ile-1
       
       do i = 1,ile
          if (ia(i).eq.blank) then
             isw=0
             cycle
          else
             if (isw.le.0) then
                jump = jump +1
                istrt(jump) = i
                inumb(jump) = 0
                isw=1
                inumb(jump) = inumb(jump) + 1
             else
                inumb(jump) = inumb(jump) + 1
             endif
          endif
       enddo
       
       if (jump.eq.0) then
          cycle
       else
          exit
       endif
    enddo
    
    return
       
40  write(iout6,140)
140 format(//1x,'Error encountered when reading input data.'//)
    stop 'inCard'

50  write(iout6,150)
150 format(//1x,'Unexpected end of input data encountered. Is STOP label missing?'//)
    stop 'inCard'

    
100 format(80a1)
101 format(2x,80a1)
    
  end subroutine inCard

  ! ### lftpos ###
  !
  !     Eliminates blanks to the left and left position chararcter string card.
  !
  subroutine lftpos(line,length)
    use params
    implicit none

    character*1, dimension(length) :: line

    integer (KIND=IPREC) :: ieff,ipos,length,ntest

    ieff = 0
    do ipos = 1, length
       if(ieff.gt.0) then
          ieff=ieff+1
          line(ieff) = line(ipos)
       end if
       if(line(ipos).ne.' '.and.ieff.eq.0) then
          ieff=1
          line(ieff) = line(ipos)
       end if
    end do

    !     fill end with trailing blanks

    do ipos = ieff+1,length
       line(ipos) = ' '
    end do

    ntest = 0
    if(ntest.ne.0) then
       write(6,*) ' Left adjusted character string '
       write(6,'(1H ,A)') line
    end if

    return
  end subroutine lftpos


  
  ! ### lowcase ###
  !
  !     Converts letters in a character string line to the lower case.
  !
  subroutine lowcase(line,length)
    use params
    implicit none
    character*1, dimension(26) :: lower,upper
    character*1, dimension(length) :: line

    integer (KIND=IPREC) :: i,icha,length

    data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/
    data upper/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

    do icha = 1, length
       do i = 1,26
          if ( line(icha).eq.upper(i) ) line(icha) = lower(i)
       end do
    end do

  end subroutine lowcase


  ! ### tabchk ###
  !
  !     Searches for a tab in the string line and replace it by a space.
  !
  subroutine tabchk(line,length)
    use params
    implicit none

    character*1 :: itab
    character*1, dimension(length) :: line

    integer (KIND=IPREC) :: i,length

    itab=char(9)
    do i=1,length
       if(line(i).eq.itab) then
          line(i)=' '
       end if
    end do

  end subroutine tabchk

  ! ### inCardG ###
  !
  !     This routine reads  a single data card and
  !
  !       - changes upper case letters into a lower case ones
  !       - changes tabs (if any) into spaces
  !       - removes leading spaces
  !
  !     Then it scans the line for nonspace fields. The number of fields
  !     is stored in jump, the starting point of a given field in istrt(i)
  !     and the number of characters in that field in inumb(i). If the
  !     routine finds an exclamation mark or a hash then anything on the
  !     same line that follows these characters is treated as a comment.
  !
  !     If iecho is nonzero the (modified) card is written to the standard output.
  !
  subroutine inCardG(iunit,iecho)
    use params
    implicit none

    character*1 :: blank,exmark,hash
    character*1, dimension(80) :: iatmp

    integer (KIND=IPREC) :: i,iecho,ike,ile,isw,it,iunit

    data blank/' '/,exmark/'!'/,hash/'#'/

    do i=1,i80
       ia(i)=' '
    enddo

    ! read a line

    jump = 0
    jrec = 0
    isw = 0
1   read(iunit,100,err=110,end=40) ia
100 format(80a1)
101 format(2x,80a1)
110 continue

    do i=1,i80
       iatmp(i)=ia(i)
    enddo

    ! set into lower case
    ! FC3, ifort: error while passing i80 to a subroutine
    it=i80
    call lowcase(iatmp,it)

    ! check for tabs, and replace them by blanks
    call tabchk(iatmp,it)

    ! check for blanks at left and remove them, fill with blanks at the end
    call lftpos(iatmp,it)

    do i=1,i80
       ia(i)=iatmp(i)
    enddo


    ! check for an exclamation mark or a hash
    do ike=1,i80
       ile=ike
       if (ia(ile).eq.exmark.or.ia(ile).eq.hash) go to 910
    enddo

910 if(ile.eq.1) go to 1
    ! echo the slightly modified input data
    if (iecho.ne.0) write(iout6,101) ia

    if(ia(ile).eq.exmark.or.ia(ile).eq.hash) ile=ile-1

    do i = 1,ile
       if (ia(i).eq.blank) then
          isw=0
          goto 30
       else
          if (isw.le.0) then
             jump = jump +1
             istrt(jump) = i
             inumb(jump) = 0
             isw=1
             inumb(jump) = inumb(jump) + 1
          else
             inumb(jump) = inumb(jump) + 1
          endif
       endif
30     continue
    enddo

    if (jump.eq.0) goto 1
    return

40  write(iout6,45)
45  format(//1x,'Error: unexpected end of input data. Is STOP label missing?'//)
    stop 'inCard'
  end subroutine inCardG

  ! ### inFloat ###
  !
  !     Extracts a floating point number from an input card.
  !
  subroutine inFloat(buf)
    use params
    implicit none

    integer (KIND=IPREC) :: i,ist,j,n,nstrt
    integer (KIND=IPREC) :: iexp,iexpdig,limit

    real (PREC) :: buf,exponent,fact,fact2
    character*1 :: itemp
    character*1, dimension(17) :: ichar

    parameter (limit=17)
    data ichar/'0','1','2','3','4','5','6','7','8','9','+','&','^','-','.','e','d'/

    jrec = jrec + 1
    buf=0.0_PREC
    if(jrec .gt. jump) return
    buf = 0.0_PREC
    fact2 = 0.0_PREC
    iexp=0
    n = inumb(jrec)
    fact = 1.0_PREC
    ist=istrt(jrec)
    nstrt = ist + n - 1
    iexpdig=0
    do i=1,n
       itemp = ia(nstrt)
       do j = 1,limit
          if(ichar(j).eq.itemp) goto 5
       enddo

       write(iout6,*) 'Error detected in inFloat'
       stop 'inFloat'

5      continue
       if (j.eq.16.or.j.eq.17) goto 12
       if (j.lt.11) goto 7
       if (j.le.14) goto 6
       if (iexp.eq.1) then
          fact2 = dble(i-1-iexpdig)
       else
          fact2 = dble(i-1)
       endif
       go to 9

12     continue
       if (iexp.eq.0) then
          buf=(0.10_PREC**fact2)*buf
          exponent=buf
          iexp=1
          buf=0
          fact=1.0_PREC
          fact2=0.0_PREC
          iexpdig=i
       endif
       goto 9

       buf = buf + dble(j-1) * fact
       fact=fact*10.0_PREC
       goto 9

6      continue
       !         if(nstrt.ne.ist.and.iexp.eq.0) go to 4
       if(j.eq.14) buf=-buf
       goto 9

7      buf = buf + dble(j-1) * fact
       fact=fact*10.0_PREC
9      continue
       nstrt = nstrt - 1
    enddo

    if (iexp.eq.1) then
       buf=(0.10_PREC**fact2)*buf
       buf=buf*10.0_PREC**exponent
       return
    endif
    buf=(0.10_PREC**fact2)*buf

  end subroutine inFloat
  
  ! ### inInt ###
  !
  !    Reads an integer from the array IA, starting at IA(istrt(jrec)) and
  !    continuing for inumb(jrec)) elements. Plus signs are ignored, the
  !    answer is accumulated in JBUF.
  subroutine inInt(jbuf)
    use params
    implicit none

    integer (KIND=IPREC) :: i,ifact,inpiexit,ist,j,jbuf,n,nstrt

    character*1 :: itemp
    character*1, dimension(14) :: ichar(14)


    data ichar/'0','1','2','3','4','5','6','7','8','9','+','&','^','-'/,inpiexit/-99999/

    jbuf = inpiexit
    jrec = jrec + 1
    if(jrec.gt.jump) goto 430
    jbuf=0
    n = inumb(jrec)
    ifact = 1
    ist=istrt(jrec)
    nstrt = ist + n - 1
    do i = 1,n
       itemp = ia(nstrt)
       do j=1,14
          if(ichar(j).eq.itemp) goto 45
       enddo
44     write(iout6,*) 'Error detected in inpi'
       stop 'inInt'

45     if(j.lt.11) goto 47
       if(nstrt.ne.ist) goto 44
       if(j.ge.14)jbuf=-jbuf
       goto 430
47     jbuf=jbuf+(j-1)*ifact
       ifact = ifact * 10
       nstrt=nstrt-1
    enddo

430 return

  end subroutine inInt

  subroutine inIntG(jbuf)
    use params
    implicit none

    integer (KIND=IPREC) :: i,ifact,inpiexit,ist,j,jbuf,n,nstrt

    character*1 :: itemp
    character*1, dimension(14) :: ichar(14)


    data ichar/'0','1','2','3','4','5','6','7','8','9','+','&','^','-'/,inpiexit/-99999/

    jbuf = inpiexit
    jrec = jrec + 1
    if(jrec.gt.jump)goto 430
    jbuf=0
    n = inumb(jrec)
    ifact = 1
    ist=istrt(jrec)
    nstrt = ist + n - 1
    do i = 1,n
       itemp = ia(nstrt)
       do j=1,14
          if(ichar(j).eq.itemp) goto 45
       enddo
44     continue
       !write(iout6,*) 'Error detected in inpi'
       !     stop 'inInt'
       return
45     if(j.lt.11)goto 47
       if(nstrt.ne.ist)goto 44
       if(j.ge.14)jbuf=-jbuf
       goto 430
47     jbuf=jbuf+(j-1)*ifact
       ifact = ifact * 10
       nstrt=nstrt-1
    enddo

430 return

  end subroutine inIntG

  ! ### inputData ###
  !
  !     Handles the input to x2DHF.  
  !
  subroutine inputData(ni_t,mu_t,no_t,nons_t)
    use params
    use discrete
    use doSCF
    use solver
    use commons
    use data4II
    
    implicit none

    integer (KIND=IPREC) :: ni_t,mu_t,no_t,nons_t
  
    character*8 readLabel
    character*20 label4m
    
    do i=1,maxdfts
       cdftex(i)=cdftext(i)
       cdftcorr(i)=cdftcorrt(i)
    enddo
    
    ! icompLEnc variables is used to quarantee the proper order of compulsory labels:
    icompLEnc=0

    ! if an extra specific card is needed (e.q. ORBPOT card with the
    ! hydrogen parameter requires LCAO card) it can be signalled by
    ! nonzero value of icompLExp; if the required label is spotted the
    ! value of icompLAdd must be increased
    icompLExp=0
    icompLAdd=0

    lcaoIncl=.false.
    ldaIncl=.false.
    omegaIncl=.false.
    
    write(*,*) '... start of input data ...'

    ! read the title of a current case (it must be the first input card)
    ! label: title

    !call inCardh(header)
    call inCardTitle(title)
    call inStr(clabel)
    if (clabel.ne.'title'.and.clabel.ne.'TITLE') then
       write(iout6,'(/2x,"Error: incorrect order of labels!")')
       write(iout6,'( 2x,"Order expected: TITLE METHOD NUCLEI CONFIG GRID ORBPOT"/)')       
       stop 'inputData'
    endif

    icompLEnc=icompLEnc+1

    clabel=""
    do while (clabel.ne.'stop')
       call inCard
       call inStr(clabel)
       call checkLabel
       if (clabel.eq.'altsweep') altSweeps=.true.
       if (clabel.eq.'break') then
          lbreakCi=.true.
       endif
       if (clabel.eq.'chktoten') call read_chktoten
       if (clabel.eq.'config') call read_config
       if (clabel.eq.'coulexch') call read_coulexch
       if (clabel.eq.'conv') call read_conv
       if (clabel.eq.'debug') call read_debug
       if (clabel.eq.'densthld') call read_densthld
       if (clabel.eq.'detnan') call read_detnan       
       if (clabel.eq.'dft') call read_dft
       if (clabel.eq.'fastscf') call read_fastscf       
       if (clabel.eq.'fefield') call read_fefield
       if (clabel.eq.'fermi') call read_fermi
       !if (clabel.eq.'fix') call read_fix
       if (clabel.eq.'fixnan') call read_fixnan
       if (clabel.eq.'fixpot') call read_fixpot
       if (clabel.eq.'fixorb') call read_fixorb       
       if (clabel.eq.'gauss') call read_gauss
       if (clabel.eq.'grid') call read_grid
       if (clabel.eq.'homo') call read_homo
       if (clabel.eq.'inout') call read_inout
       if (clabel.eq.'interp') call read_interp
       if (clabel.eq.'kinpot') call read_kinpot           
       if (clabel.eq.'lm'.or.clabel.eq.'lm0') then
          lmtype=0
          call read_offDiagLM
       endif

       if (clabel.eq.'lm1') then
          lmtype=1
          call read_offDiagLM
       endif

       if (clabel.eq.'lm2') then
          lmtype=2
          call read_offDiagLM
       endif


       if (clabel.eq.'extracul') call read_extracule
       if (clabel.eq.'intracul') call read_intracule
       if (clabel.eq.'lcao') call read_lcao
       if (clabel.eq.'lxcpolar') lxcPolar=.true.
       if (clabel.eq.'mcsor') call read_mcsor
       if (clabel.eq.'mcsor-o') call read_mcsor_o
       if (clabel.eq.'mcsor-ce') call read_mcsor_ce
       if (clabel.eq.'method') call read_method
       if (clabel.eq.'mmoments') call read_mmoments       
       if (clabel.eq.'multipol') call read_multipol
       if (clabel.eq.'nuclei') call read_nuclei    
       if (clabel.eq.'omega') call read_omega
       if (clabel.eq.'omegaz') lomegaz=.true.
       if (clabel.eq.'omegaopt') call read_omegaopt       
       if (clabel.eq.'orbpot') call read_orbpot
       if (clabel.eq.'order') call read_order
       if (clabel.eq.'out4pair') call read_out4pair       
       if (clabel.eq.'plot') call read_plot       
       if (clabel.eq.'potgsz') call read_potgsz
       if (clabel.eq.'potgszg') call read_potgszg
       if (clabel.eq.'potharm2') call read_potharm2
       if (clabel.eq.'potharm3') call read_potharm3       
       if (clabel.eq.'pothooke') call read_pothooke
       if (clabel.eq.'potcoul2') call read_potCoul2           
       if (clabel.eq.'potcoul3') call read_potCoul3    
       if (clabel.eq.'potkh') call read_potkh
       if (clabel.eq.'potsap') call read_potSAP       
       if (clabel.eq.'print') call read_print
       if (clabel.eq.'prtevery') call read_prtevery
       if (clabel.eq.'scf') call read_scf
       if (clabel.eq.'slowexch') call read_slowexch       
       if (clabel.eq.'sor') call read_sor
       if (clabel.eq.'sor4orb') call read_sor4orb
       if (clabel.eq.'sor4pot') call read_sor4pot
       if (clabel.eq.'stop') call read_stop(ni_t,mu_t,no_t,nons_t)
       if (clabel.eq.'tail') call read_tail
       if (clabel.eq.'xalpha') call read_xalpha
    enddo

    return

  end subroutine inputData

  ! ### checkLabel ###
  !
  !     Checks a label.  
  !
  subroutine checkLabel
    use data4II
    implicit none
    
    do i=1, nlabels
       if (clabel.eq.labellc(i)) return
    end do
    
    write(iout6,100) clabel,(labellc(i),i=1,nlabels)
100 format(/2x,'Error: label "',a8,'" is not supported. ',/,/2x,&
         'Try one of the following:',/,10(8x,a8,2x,a8,2x,a8,2x,a8,2x,a8/))
    write(iout6,*)
    stop 'checkLabel'
  end subroutine checkLabel

  subroutine checkLabelsOrder (labelPos)
    use params
    use data4II
    implicit none
    integer :: labelPos
    if (icompLEnc.ne.labelPos) then
       write(iout6,'(/2x,"Error: incorrect order of labels.")')
       write(iout6,'( 2x,"Order expected: TITLE, METHOD, NUCLEI, CONFIG, GRID, ORBPOT"/)')
       stop 'checkLabelsOrder'
    endif
  end subroutine checkLabelsOrder

  subroutine read_chktoten
    use params
    use discrete
    use commons
    
    implicit none
    lcheckTotEnergy=.true.
    return
  end subroutine read_chktoten

  
  subroutine read_config
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none

    call checkLabelsOrder(3)
    
    mt=0
    icompLEnc=icompLEnc+1

    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       itotq=itmp
    else
       itotq=0
    endif
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
       write(iout6,'(/2x,"Error: wrong symmetry or symmetry higher than phi")')
       stop 'read_config'
    endif
    
    isum0=isum
    do i=1,nbsym
       isum=isum+1
       orbsym(isum)=char8
       bond(isum)=orbsym(isum)
       gusym(isum)=space
    enddo
    
    if (isum.gt.maxorb) then
       write(iout6,*) 'too many orbitals (max. see maxorb)'
       stop 'read_config'
    endif
    ! in homonuclear case read in the u/g parity of each orbital
    ! if BREAK is on (lbreakCi=.true.) u/g labels are ignored
    ! if |z1-z2|<homolevl then the case is treated as a homonuclear one
    ! (see setDefaults)

    if (abs(z1-z2).lt.homolevl.and.lbreakCi) then
       call inStr(gusym(isum))
       isum1=isum0
       do i=1,nbsym
          isum1=isum1+1
          gusym(isum1)=gusym(isum)
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
       ! if u/g symmetry is assigned to orbitals their symmetry must be forced
       !ihomon=2
       !homoIncl=.true.
       call inStr(guttmp)
       char8=guttmp
       if (guttmp.eq.'u'.or.guttmp.eq.'g') then
       lhomonucl=.true.
          gusym(isum)=guttmp
          isum1=isum0
          do i=1,nbsym
             isum1=isum1+1
             gusym(isum1)=gusym(isum)
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
       ifixorb(i)=0
       inhyd(i)=1
       inDFT(i)=1       
    enddo
    
    if (char8.eq.endl) then
       totchar=0.0_PREC
       do iorb=1,norb
          totchar=totchar+occ(iorb)
       enddo
       
       if (abs(totchar-totq).gt.precis*1000.0_PREC) then
          write(iout6,'(/,2x,"Warning: mismatch in given and calculated total charge:",2f6.1/)') totchar,totq
       endif
       
       ! set some additional parameters
       
       do iorb=1,norb
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
       enddo
       return
    endif
    goto 100
    return
  end subroutine read_config
  
  subroutine read_conv
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
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


  subroutine read_coulexch
    use params
    use discrete
    use commons
    use data4II
    use solver    
    
    implicit none
    ! label: coulexch
    lcoulexch=.true.

    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       nthreads4coulexch = itmp
       ! nthreads4coulexch can be only set to 1 o maxpots
       if ( nthreads4coulexch.ne.1 ) then
          write(iout6,'(/2x,"Error: incorrect entry: must be 1 if specified.")')
          stop 'read_dft'
       endif
    endif

    return
  end subroutine read_coulexch

  
  subroutine read_debug
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none

    ! label: debug
    ! idbg  -  additional printout or action for testing purposes
    !          set debug flags. if an integer i is encountered debuf
    !          flag i is set, i.e. idebug(i)=1. At most 999 integers can be specified.
    
    ! List of used debug flags: xhf -D list
    inzero=1
    do i=1,maxflags
       id(i)=0
       call inInt(id(i))
       if (id(i).gt.0) then
          inzero=inzero+1
          idebug(id(i))=1
       endif
    enddo

    return
  end subroutine read_debug


  subroutine read_detnan
    use params
    use commons

    implicit none
    
    ! label: detNaN
    ! if the DFT method is used detect NaN values in Libxc potentials

    ldetectNaN=.true.
    return
  end subroutine read_detNaN


  subroutine read_densthld
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    call inFloat(ftmp)
    dens_threshold=ftmp
    return
  end subroutine read_densthld


  
  subroutine read_fixnan
    use params
    use commons

    implicit none
    
    ! label: fixNaN
    ! if the DFT method is used fix NaN values found in some Libxc potentials

    lfixNaN=.true.
    return
  end subroutine read_fixNaN
  
  subroutine read_dft
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
#ifdef LIBXC
    use, intrinsic :: iso_c_binding
    use xc_f90_lib_m
    use libxc_funcs_m    
#endif
 
    implicit none

#ifdef LIBXC    
    type(xc_f90_func_info_t) :: xc_info
    type(xc_f90_func_t) :: xc_func
    type(xc_f90_func_info_t) :: info
    integer(c_int) :: vmajor, vminor, vmicro, family_id, func_id, kind_id, err
    character(len=120) :: name, kind1, family, ref
#endif
    
    ! label: dft [lda|b88] [lyp] [vwn]
    ! if the DFT method is used one can choose exchange and correlation
    ! energy functionals
    ! call inInt(itmp) 
    
    ! if only dft label is present calculate exchange (LDA, B88, PW86,
    ! PW91) and correlation (LYP) contributions to total energy at the
    ! end of SCF process

    if (.not.DFT.and.HFinput) then
       write(iout6,'(/2x,"Warning: assuming DFT method (see User''s guide)."/)')
       HFinput=.false.
       !stop 'read_dft'
    endif
    call inStr4lxc(char30)

    ! first try to find the old functional labels 
    char8=trim(char30)
    idftex=0
    idftcorr=0
    do i=1,maxdfts
       if (char8.eq.cdftex(i)) idftex=i
       if (char8.eq.cdftcorr(i)) idftcorr=i
    enddo

    if (idftex/=0.or.idftcorr/=0) then
       if (TED) then
          write(iout6,'(/2x,"Error: no support for DFT potential(s) within TED metod."/)')
          stop 'read_dft'
       endif

       call inStr(char8)
       if (char8.ne.endl) then
          do i=1,maxdfts
             if (char8.eq.cdftex(i)) idftex=i
             if (char8.eq.cdftcorr(i)) idftcorr=i
          enddo
          if (idftex.eq.0.and.idftcorr.eq.0) then
             write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide.")')
             stop 'read_dft'
          endif
       endif
       
       if ((idftex.gt.2).or.(idftcorr.gt.2)) then
          write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide.")')
          stop 'read_dft'
       endif
       DFT=.true.
       HF=.false.
       return
    endif
    if (.not.lxcSupp) then
       write(iout6,'(/2x,"Error! No libxc support detected. Try ./x2dhfctl -bl.")') 
       stop "read_dft"
    endif

#ifdef LIBXC
    if (char30(1:3)=="hyb" .or. char30(1:6)=="xc_hyb") then
       lxcHyb=.true.
       HFinput=.true.
    endif
    if (trim(char30).ne."") then
       lxcFuncs=0
       do j1=1,2
          lxcFuncs=lxcFuncs+1
          lxcFuncs2use(lxcFuncs)=xc_f90_functional_get_number(trim(char30))
          if (lxcFuncs2use(lxcFuncs)<=0) then
             write(iout6,'(/2x,"Error! ",a20,": no such libxc functional found.")') trim(char30)
             stop "read_dft"
          endif
          if (lxcPolar) then
             call xc_f90_func_init(xc_func, lxcFuncs2use(lxcFuncs), XC_POLARIZED, err)
          else
             call xc_f90_func_init(xc_func, lxcFuncs2use(lxcFuncs), XC_UNPOLARIZED, err)
          endif
          
          info=xc_f90_func_get_info(xc_func)

          select case (xc_f90_func_info_get_family(info))
          case(XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_GGA)
             call inStr4lxc(char30)
             if (trim(char30).eq."") exit
          case default
             write(iout6,'(/2x,"Error! ",a8,": unsupported libxc functional."/)') trim(char30)
             stop 'read_dft'
          end select
       enddo
       call xc_f90_func_end(xc_func)
       if (lxcFuncs>0) then
          LXC=.true.
          DFT=.false.
          HF=.false.
          return
       else
          write(iout6,'(/2x,"Error! A libxc functional requested but no libxc support vailable."/)') 
          stop "read_dft"
       endif
    endif
#endif
  end subroutine read_dft
  
  subroutine read_fastscf
    use params
    use discrete
    use commons
    
    implicit none
    
    lfastscf=.true.
    return
  end subroutine read_fastscf
  
  subroutine read_fefield
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
        
    implicit none
    ! label: fefield
    
    ! finite field calculations
    ! 1 field strength
    ! 2 cutoff
    
    lfefield=.true.
    lmmoments=.true.
    
    call inFloat(ffield)
    ! issue error message if homo is on
    !if (ihomon.eq.2) then
    if (lhomonucl) then
       write(iout6,'(/2x,"Error: FEFIELD and g/u symmetry labels are mutually exclusive.")')
       write(iout6,'(/2x,"Try adding BREAK label."/)')
       stop 'read_fefield'
    endif
    return
  end subroutine read_fefield

  subroutine read_fermi
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
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
       write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
       stop 'read_fermi'
    endif
    lpotFermi=.true.
    lpotCoulomb=.false.

    if (z1<3.0_PREC.and.z2<3.0_PREC) then
       write(iout6,'(/2x,"Error: if Z1/Z2 is non-zero, then it must be greater than 2."/)')
       stop "read_fermi"
    endif
    return
  end subroutine read_fermi

  subroutine read_slowexch
    use params
    use discrete
    use commons
    
    implicit none
    
    lfastexch=.false.
    return
  end subroutine read_slowexch

  subroutine read_fixorb
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none

    inzero=0
    do i=1,maxflags
       id(i)=0
       call inInt(id(i))
       if (id(i).gt.0) then
          inzero=inzero+1
          if (inzero.gt.norb) then
             write(iout6,'(/2x,"Error: too many items - see User''s guide."/)')
             stop 'read_fixorb'
          endif
          ifixorb(norb-id(i)+1)=1
       endif
    enddo
    
    ! if no orbitals are given fix them all
    if (inzero.eq.0.or.inzero==norb) then
       lfixorb=.true.
       do i=1,norb
          ifixorb(i)=1
       enddo
    endif
    return
  end subroutine read_fixorb

  subroutine read_fixpot
    use params
    use discrete
    use scfshr
    use solver 
    use commons
    use data4II
    
    implicit none
    lfixcoul=.true.
    lfixexch=.true.
    return
  end subroutine read_fixpot
  
  subroutine read_gauss
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
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
       write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
       stop 'read_gauss'
    endif
    lpotGauss=.true.
    lpotCoulomb=.false.
    return
  end subroutine read_gauss

  subroutine read_grid
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    
    call checkLabelsOrder(4)
    icompLEnc=icompLEnc+1

    tmp1=0.0_PREC
    tmp2=0.0_PREC
    call inInt(nni)
    ! adjust nni first to provide correct value of hnu for nmucalc
    nni=nnucalc(nni)
    call inFloat(tmp1)
    call inFloat(tmp2)
    if (abs(tmp2)>epsilon(zero)) then
       ! nnu nmu R_infty
       n=nint(tmp1)
       nmu(ngrids)=nmucalc(n)
       rgrid(ngrids)=tmp2
    else
       ! nni R_infty
       if (ngrids.ne.1) then
          write(iout6,*)'read_grid: missing item on GRID card'
          stop 'read_grid'
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
    
    if (nni.lt.7) then
       write(iout6,'(/2x,"Error: too few grid points in nu variable."/)') 
       stop 'read_grid'
    endif
    if (nni.gt.maxnu) then
       write(iout6,'(/2x,"Error: too many grid points in nu variable: ",i4," is greater than ",i4/)') nni,maxnu
       stop 'read_grid'
    endif
    
    if (nmutot.lt.7) then
       write(iout6,'(/2x,"Error: too few grid points in mu variable."/)') 
       stop 'read_grid'
    endif
    if (nmutot.gt.maxmu) then
       write(iout6,'(/2x,"Error: too many grid points in mu variable: ",i4," is greater than ",i4/)') &
            nmutot,maxmu
       stop 'read_grid'
    endif

    return
  end subroutine read_grid

  subroutine read_homo
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II   
    
    implicit none
    ! label: homo --  when this label is encountered the g/u symmetry of
    !                 orbitals is strictly imposed; in all other respects
    !                 a homonuclear case is treated as a heteronuclear one
    
    ! lhomonucl=.false --  heteronuclear case
    ! lhomonucl=.true. --  homonuclear case (g/u symmetry is forced)
    
    lhomonucl=.true.
    return
  end subroutine read_homo

  subroutine read_inout
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! label: inout i32|i64|r128
    ! force format of input and output data
    lengthintin=lengthint
    lengthfpin=lengthfp
    call inStr(clabel)
    if     (clabel.eq.'i32') then 
       lengthintin=4
       lengthfpin=8
       inUnformatted=.true.
       inFormatI32=.true.
    elseif (clabel.eq.'i32f') then 
       lengthintin=4
       lengthfpin=8
       inUnformatted=.false.
       inFormatI32=.true.
       formfp=formfp64
    elseif (clabel.eq.'i64') then
       lengthintin=8
       lengthfpin=8
       inUnformatted=.true.
       inFormatI64=.true.
    elseif (clabel.eq.'i64f') then
       lengthintin=8
       lengthfpin=8
       inUnFormatted=.false.
       inFormatI64=.true.
       formfp=formfp64
    elseif (clabel.eq.'r128') then
       lengthintin=8
       lengthfpin=16
       inUnformatted=.true.
       inFormatR128=.true.
    elseif (clabel.eq.'r128f') then
       lengthintin=8
       lengthfpin=16
       inUnformatted=.false.
       inFormatR128=.true.
       formfp=formfp128
    else
       write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide"/)')
       stop 'read_inout'
    endif
    
    call inStr(clabel)
    if     (clabel.eq.'i32') then 
       outUnformatted=.true.
       outFormatI32=.true.
    elseif (clabel.eq.'i32f') then
       outFormatI32=.true.
       outUnformatted=.false.
    elseif (clabel.eq.'i64') then
       outUnformatted=.true.
       outFormatI64=.true.
    elseif (clabel.eq.'i64f') then
       outUnformatted=.false.
       outFormatI64=.true.
       formfp=formfp64
    elseif (clabel.eq.'r128') then
       outUnformatted=.true.
       outFormatR128=.true.
    elseif (clabel.eq.'r128f') then
       formfp=formfp128
       outUnformatted=.false.
       outFormatR128=.true.
       if (lengthfp.ne.16) then
          write(iout6,'(/2x,"Warning! The present build of the program does not support quadruple precision.")')
          write(iout6,'(2x,"Try running ./x2dhfctl -r 16."/)')
          stop 'read_inout'
       endif
    else
       write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide"/)')
       stop 'read_inout'
    endif
  end subroutine read_inout

  subroutine read_interp
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
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
       write(iout6,'(/2x,"Error: incorrect order of interpolation polynomial - allowed values are 2, 4, 6 or 8"/)')
       stop 'read_interp'
    endif
    iinterp=1
    return
  end subroutine read_interp


  subroutine read_kinpot
    use params
    use commons
    
    implicit none
    
    iout4kinpot=1

    return
  end subroutine read_kinpot

  
  subroutine read_offDiagLM
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
   
    implicit none
    ! label: lagra
    ! scaling and damping factors for off-diagonal Lagrange multipliers
    loffDiagLM=.true.
    do i=1,maxorb*(maxorb-1)/2,2
       call inInt(itmp1)
       call inInt(itmp2)
       if (itmp1.ne.inpiexit.and.itmp2.ne.inpiexit) then
          nlm(norb-itmp1+1,norb-itmp2+1)=1
          nlm(norb-itmp2+1,norb-itmp1+1)=1
          offDiagLM(norb-itmp2+1,norb-itmp1+1)=.true.
          offDiagLM(norb-itmp1+1,norb-itmp2+1)=.true.
       elseif (itmp1.eq.inpiexit.and.itmp2.eq.inpiexit) then
          return
       elseif (itmp1.ne.inpiexit.or.itmp2.ne.inpiexit) then
          write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
          stop "read_offDiagLM"
       endif
    enddo
  end subroutine read_offDiagLM

  subroutine read_lcao
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
   
    implicit none
    logical :: inhydStatus, inDFTStatus

    if (ldaIncl.or.hfIncl) then
       !write(iout6,'(/2x,"Error: label LCAO is incompatible with ORBPOT HF|LDA card ..."/)')
       write(iout6,'(6x,"Warning! Label LCAO is incompatible with ORBPOT HF|LDA card.")')
       write(iout6,'(6x,"This and the following",i3," LACO cards are read but ignored.")') norb
       do iorb=1,norb
          call inCard
       enddo
       return
    endif
    
    call checkLabelsOrder(6)

    ! icompLEnc=icompLEnc+1
    icompLAdd=icompLAdd+1

    inclorb=1
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       inclorb=itmp
    endif
    lcaoIncl=.true.

    if (linitFuncsHydrogen.and.(inclorb.ne.1.and.inclorb.ne.2)) then
       write(iout6,'(/2x,"Error: incompatible parameters - no LCAO data present."/)')
       stop 'read_lcao'
    endif

    ! method OED
    if (OED) then
       lfixexch=.true.
       if (.not.linitFuncsHydrogen) initFuncsOED=.true.
       recalcMMfactor=10.0_PREC**(15)
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
          if (inhyd(iorb).eq.inpiexit) then
             inhyd(iorb)=1
             inDFT(iorb)=1             
          else
             call inInt(inDFT(iorb))
             if (inDFT(iorb).eq.inpiexit) then
                inDFT(iorb)=1
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
          if (inhyd(iorb).eq.inpiexit) then
             inhyd(iorb)=1
             inDFT(iorb)=1             
          else
             call inInt(inDFT(iorb))
             if (inDFT(iorb).eq.inpiexit) then
                inDFT(iorb)=1
             endif
          endif
       enddo
    endif

    ! let's check that lcao data are consistent
    inhydStatus=.false.
    inDFTStatus=.false.
    do iorb=1,norb
       inhydlcao(iorb)=inhyd(iorb)
       if (inhyd(iorb)/=0) inhydStatus=.true.
       if (inDFT(iorb)/=0) inDFTStatus=.true.       
    enddo

    if (inhydStatus.eqv. .false. ) then
       write(iout6,'(/2x,"Warning: check consistency of 9th parameters on LCAO cards"/)')
       !stop "read_lcao"
    endif

    if (inDFTStatus.eqv. .false. ) then
       write(iout6,'(/2x,"Warning: check consistency of 10th parameters on LCAO cards"/)')
       !stop "read_lcao"       
    endif
    
    ! normalize mixing coefficients
    do iorb=1,norb
       co12=abs(co1(iorb))+abs(co2(iorb))
       if (abs(co12)>epsilon(zero)) then
          co1(iorb)=co1(iorb)/co12
          co2(iorb)=co2(iorb)/co12
       endif
    enddo

    ! ini=5 enables to retrieve from complete/uncomplete dump file if only
    ! one (the highest= the first) orbital is missing the rest is retreived
    ! and the first is initialized as a hydrogenic one (ini is changed from
    ! 5 to 4)

    return
  end subroutine read_lcao

  subroutine read_maxsor
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    
    do iorb=1,norb
       call inInt(itmp1)
       if (itmp1.ne.inpiexit) then
          call inInt(itmp2)
          if (itmp2.ne.inpiexit) then
             maxsororb(itmp1)=itmp2
          else
             write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
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
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
  
    implicit none
    ! label: mcsor
    !
    !       nthreads - number of threads used to relax grid points of a given colour 
    !   maxsororb(2) - maximal number of micro (mc)sor iterations during relaxation
    !                  of every orbital in an SCF cycle
    !   maxsorpot(2) - maximal number of micro (mc)sor iterations during relaxation
    !                  of every potential in an SCF cycle

    lorbmcsor=.true.
    lpotmcsor=.false.
    if (mcsorpt) then
       lmcsorpt=.true.
    endif

    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       nthreads = itmp
       if (nthreads>maxthreads4mcsor) then
          write(iout6,'(/2x,"Error: nthreads cannot exceed",i3," See params.f90 and sorpt.h."/)') maxthreads4mcsor
          stop 'readLabels'
       endif
    endif
    
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsororb(2)=itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsorpot(2)=itmp
       endif
    endif

    return
  end subroutine read_mcsor


  subroutine read_mcsor_o
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: mcsor-o
    !       nthreads - number of threads used to relax grid points of a given colour 
    !   maxsororb(2) - maximal number of micro (mc)sor iterations during relaxation
    !                  of every orbital in an SCF cycle
    !   maxsororb(1) - maximal number of macro (mc)sor iterations during relaxation
    !                  of every potential in an SCF cycle
    lorbmcsor=.true.
    if (mcsorpt) then
       lmcsorpt=.true.
    endif

    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       nthreads = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsor2 = itmp
          maxsororb(2)=itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             maxsor3 = itmp 
             maxsororb(1)=itmp
          endif
       endif
    endif

    return
  end subroutine read_mcsor_o

  subroutine read_mcsor_ce
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! label: mcsor-ce
    !    nthreads - number of threads used to relax grid points of a given colour 
    !   maxsorpot(2) - maximal number of micro mcsor iterations during relaxation
    !                  of every potential in an SCF cycle
    !   maxsorpot(1) - maximal number of macro mcsor iterations during relaxation
    !                  of every potential in an SCF cycle
    lpotmcsor=.true.
    if (mcsorpt) then
       lmcsorpt=.true.
    endif
    
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       nthreads = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsor2 = itmp
          maxsorpot(2)=itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             maxsor3 = itmp 
             maxsorpot(1)=itmp
          endif
       endif
    endif
       
    return
  end subroutine read_mcsor_ce

  subroutine read_method
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none

    ! imethod=1 -- hf
    ! imethod=2 -- oed
    ! imethod=3 -- hfs
    ! imethod=4 -- dft
    ! imethod=5 -- scmc

    call checkLabelsOrder(1)
    icompLEnc=icompLEnc+1
    
    call inStr(char8)
    
    ihit=0
    do i=1,nmethods
       if (trim(char8).eq.trim(cmethod(i))) then
          imethod=i
          if (i==1) HF=.true.
          if (i==2) OED=.true.
          if (i==3) HFS=.true.
          if (i==4) DFT=.true.
          if (i==5) SCMC=.true.
          if (i==6) TED=.true.
          ihit=1
          exit
       endif
    enddo
    if (ihit==0) then
       write(iout6,1000)
1000   format (/2x,"Error: incorrect method. Allowed values: ",&
            "DFT, HF, HFS, LXC, OED, SCMS or TED."/)
       stop 'read_method'
    endif
    
    ! HFS method is equivalent to DFT one with the LDA echange potential and the optimum
    ! value of alpha
    
    ! for HFS method set alpha to a (hopefully) reasonable value based
    ! on atomic optimum values due to Schwarz (see blk-data.f90)
    
    ! see winding up section of this routine

    if (HF) HFinput=.true.
    
    if (HFS) then
       DFT=.true.
       HF=.false.
       idftex=1
       if ( nint(z1).ge.nint(z2) )  then
          alphaf=alphaopt(nint(z1))
       else
          alphaf=alphaopt(nint(z2))
       endif
    endif

    ! for DFT method the exchange potential is set to LDA and alpha is set (by default) to
    ! 2/3
    if (DFT) then
       HF=.false.
       idftex=1
       alphaf=two/three
    endif
    
    ! alpha parameter of the density-functional theory is calculated according to the
    ! Self-Consistent Multiplicative Constant method (see V.V.Karasiev and E.V.Ludenia,
    ! Self-consistent multiplicative constant method for the exchange energy in density
    ! functional theory, Phys. Rev. A 65 (2002) 062510
    
    if (SCMC) then
       !DFT=.true.
       HF=.false.
       idftex=1
       alphaf=two/three
    endif

    return
  end subroutine read_method

  subroutine read_mmoments
    use params
    use commons
    
    implicit none
    lmmoments=.true.
  end subroutine read_mmoments

  
  subroutine read_multipol
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! label: multipol 
    ! define number of multipole moments to be calculated
    !     recalcMMfactor  -  if recalcMMfactor>0  multipole moment expansion
    !		 coefficients are recalculated	every time deltaEE,
    !		 i.e. maximum error in orbital energy, changes by
    !		 the recalcMMfactor factor
    !		 if recalcMMfactor<0 these coefficients are not calculated
    !      mpole  - number of multipole expansion coefficients
    
    call inFloat(recalcMMfactor)
    if (recalcMMfactor<epsilon(zero)) recalcMMfactor=ten**(15)
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       mpole=itmp
       if (mpole.lt.2.or.mpole.gt.maxmpole) then
          write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
          stop 'read_multipol'
       endif
    endif
    return
  end subroutine read_multipol

  subroutine read_nuclei
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none

    call checkLabelsOrder(2)
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
    if ((abs(z1)+abs(z2))<epsilon(zero).or.abs(r)<epsilon(zero)) then
       write(iout6,'(/2x,"Error: incomplete input."/)')
       stop 'read_nuclei'
    endif
    return
  end subroutine read_nuclei

  subroutine read_omega
    use commons
    use params
    use discrete
    use data4II
    use scfshr
    use solver
    
    implicit none
    ! label: omega
    !  ovforb   -  overelaxation parameter for orbitals 
    !  ovfcoul  -  overelaxation parameter for potentials

    omegaIncl=.true.    
    call inFloat(ftmp1)
    call inFloat(ftmp2)
    if (abs(ftmp1)>epsilon(zero).or.abs(ftmp2)>epsilon(zero)) then
       if (abs(ftmp1)>epsilon(zero)) ovforb=ftmp1
       if (abs(ftmp2)>epsilon(zero)) ovfcoul=ftmp2
    else
       ! old format
       call inCard
       call inFloat(ftmp)
       if (ftmp.gt.0.0_PREC) then
          ovforb=ftmp
       else
          write(iout6,'("Error: missing or incorrect entry - see User''s guide."/)')
          stop 'read_omega'
       endif
       
       call inCard
       call inFloat(ftmp)
       ovfcoul=ftmp
    endif
    ovfexch = ovfcoul
    return
  end subroutine read_omega

  subroutine read_omegaopt
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    !  label: omegaopt
    
    !  1 - near-optimal omega values (based on FH tests) for cases when
    !      good initial approximations are available or when fixed
    !      orbitals/potentials calculations are performed

    !  two optional parameters for scaling the chosen omega values for orbitals and 
    !  potentials
    call inInt(iomega)
    if (iomega.eq.0) then
       iomega=1
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
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none

    call checkLabelsOrder(5)
    icompLEnc=icompLEnc+1
    
    call inStr(clabel)
    if (clabel.ne.'old') then
       call inInt(itmp)
       ! if ienterm4init > 0 then initial SCF iterations are performed with orbitals being fixed
       ! until energy differences drop below 10**(-ienterm4init)
       ! by default ienterm4init=1
       if (itmp.ne.inpiexit) then
          ienterm4init=itmp
       endif
       !if (OED.or.TED) then
       if (OED) then
          ienterm4init=0
          lfixorb4init=.false.
       endif

       if (ienterm4init>0) then
          lfixorb4init=.true.
          lfixorb=.true.       
          do i=1,norb
             ifixorb(i)=1
          enddo
       endif
    endif
    
    if     (clabel.eq.'hydrogen') then
       !ini=1
       icompLExp=1
       linitFuncsHydrogen=.true.
       if (ienterm4init==0) then
          lfixorb4init=.false.
          lfixorb=.false.       
          do i=1,norb
             ifixorb(i)=0
          enddo
       endif
    elseif (clabel.eq.'gauss') then
       !ini=2
       linitFuncsGauss=.true.
    elseif (clabel.eq.'gauss-c') then
       !ini=3
       linitFuncsGaussC=.true.
    elseif (clabel.eq.'old') then
       !ini=5
       linitFuncsOld=.true.
    elseif (clabel.eq.'noexch') then
       !ini=6
       linitFuncsNoexch=.true.
    elseif (clabel.eq.'qrhf') then
       !ini=11
       linitFuncsQRHF=.true.
    elseif (clabel.eq.'lda') then
       !ini=12
       ldaIncl=.true.
       linitFuncsLDA=.true.
    elseif (clabel.eq.'ldasap') then
       !ini=12
       ldaIncl=.true.
       ldaSAPIncl=.true.
       linitFuncsLDA=.true.
       ienterm4init=1
       lfixorb4init=.true.
       lfixorb=.true.       
       do i=1,norb
          ifixorb(i)=1
       enddo
    elseif (clabel.eq.'hf') then
       !ini=13
       hfIncl=.true.
       linitFuncsHF=.true.
    elseif (clabel.eq.'molcas') then
       !ini=22
       linitFuncsMolcas=.true.
    elseif (clabel.eq.'nodat') then
       !ini=55
       linitFuncsNodat=.true.
    else
       write(iout6,1000)
1000   format(/2x,"Error: incorrect source of orbitals and potentials. Try: ",&
                 /9x,"HYDROGEN, GAUSS, GAUSS-C, HF, LDA, MOLCAS, NOEXCH, NODAT, OLD, QRHF. "/)
       stop 'read_orbpot'
    endif

    if (clabel.eq.'qrhf'.or.clabel.eq.'hf'.or.clabel.eq.'lda') then
       if (abs(z2)>epsilon(zero).and.(z2-z1)>epsilon(zero)) then
          write(iout6,1002)
1002      format (/2x,"Warning: Z1 < Z2. Swap atoms in the nuclei card to run the case."/)
          stop "read_orbpot"
       endif
    endif

  end subroutine read_orbpot

  subroutine read_order
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! label: order

    ! arg  - type of (mc)sor sweeps (ordering) for each grid
    !      = col-wise  - the natural column-wise ordering for mesh points
    !	   = middle    - the 'middle' type of sweep (default) 
    !	   = row-wise  - the natural row-wise ordering for mesh points
    !	   = rcol-wise - the reversed natural column-wise ordering (see mesh for details)

    
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
    character*8 :: arg
    call inStr(arg)
    meshOrdering=trim(arg)
    return
  end subroutine read_order

  subroutine read_out4pair
    use params
    use commons
    implicit none
    iout4pair=1
    return
  end subroutine read_out4pair

  subroutine read_potgsz
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! Green, Sellin, Zachor model potential 
    if (.not.OED) then        
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)') cmethod(imethod)
       stop 'read_potgsz'
    endif
    lpotGSZ=.true.
    lpotCoulomb=.false.
    return
  end subroutine read_potgsz

  subroutine read_potgszg
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! Green, Sellin, Zachor model potential + finite Gauss nucleus model 
    lpotGSZG=.true.
    lpotCoulomb=.false.
    if (.not.OED) then            
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)')&
            cmethod(imethod)
       stop 'read_potgszg'
    endif

    return
  end subroutine read_potgszg

  subroutine read_potsap
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! Susi Lehtola: Superposition of Atomic Potentials
    lpotSAP=.true.
    lpotCoulomb=.false.
    if (.not.OED) then                
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)')&
            cmethod(imethod)
       stop 'read_potsap'
    endif
    
    return
  end subroutine read_potsap

  subroutine read_potCoul2
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    
    call inInt(mpot)
    call inFloat(apot)
    call inFloat(v0pot)
    lpotCoul2=.true.
    lpotCoulomb=.false.

    if (.not.OED) then
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)') cmethod(imethod)
       stop 'read_potCoul2'
    endif
    return
  end subroutine read_potCoul2
  
  subroutine read_potCoul3
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    
    call inInt(mpot)
    call inFloat(apot)
    call inFloat(v0pot)
    lpotCoul3=.true.
    lpotCoulomb=.false.

    if (.not.OED) then    
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)') cmethod(imethod)
       stop 'read_potCoul3'
    endif
    return
  end subroutine read_potCoul3

  subroutine read_potharm2
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    real (PREC) tmp
    ! harmonic potential 
    lpotHarm2=.true.
    lpotCoulomb=.false.
    call inFloat(tmp)
    if (abs(tmp)>epsilon(zero)) then
       hooke=tmp
    endif
    
    if (.not.OED) then    
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)') &
            cmethod(imethod)
       stop 'read_potharm2'
    endif
    return
  end subroutine read_potharm2
  
  subroutine read_potharm3
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    real (PREC) tmp
    ! harmonic potential 
    lpotHarm3=.true.
    lpotCoulomb=.false.
    if (.not.OED) then        
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)') &
            cmethod(imethod)
       stop 'read_potharm3'
    endif

    call inFloat(tmp)
    if (abs(tmp)>epsilon(zero)) then
       hooke=tmp
    endif

    return
  end subroutine read_potharm3

  subroutine read_pothooke
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    real (PREC) tmp
    ! harmonic potential (Hooke's atom, harmonium) 
    call inFloat(tmp)
    if (abs(tmp)>epsilon(zero)) then
       hooke=tmp
    endif
    lpotHooke=.true.
    lpotCoulomb=.false.
    return
  end subroutine read_pothooke

  subroutine read_intracule
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    real (PREC) tmp
    ! intracular Hamiltonian
    call inFloat(tmp)
    if (abs(tmp)>epsilon(zero)) then
       hooke=tmp
    endif
    lintracule=.true.
    lpotCoulomb=.false.
    OED=.true.
    return
  end subroutine read_intracule


  subroutine read_extracule
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    real (PREC) tmp
    ! extracular Hamiltonian (spherical harmonic oscillator)
    call inFloat(tmp)
    if (abs(tmp)>epsilon(zero)) then
       hooke=tmp
    endif
    lextracule=.true.
    lpotCoulomb=.false.
    OED=.true.
    return
  end subroutine read_extracule

  subroutine read_plot
    use params
    use commons
    use data4II
    
    implicit none
    ! alphaf - if alphaf is not zero a DFT energy functional is used
    
    lplot=.true.
    iplot=1
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       iplot=itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          istep = itmp
       endif
    endif
    
  end subroutine read_plot

  
  subroutine read_potkh
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: potkh
    
    apot=1.0_PREC
    v0pot=1.0_PREC
    nsimp=1000
    call inInt(mpot)
    call inFloat(epspot)
    call inFloat(ompot)
    call inFloat(tmp1)
    if (abs(tmp1)>epsilon(zero)) then
       apot=tmp1
       call inFloat(tmp2)
       if (abs(tmp2)>epsilon(zero)) then
          v0pot=tmp2
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             nsimp=itmp
          endif
       endif
    endif
    lpotKH=.true.
    lpotCoulomb=.false.
    if (.not.OED) then        
       write(iout6,'(/2x,"Error: this potential cannot be used with ",a3," method! Try OED instead."/)')&
            cmethod(imethod)
       stop 'read_potkh'
    endif
    return
  end subroutine read_potkh

  subroutine read_print
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    ! label: print
    ! label followed by up to maxflags integers (1..999) specifying additional printouts

    ! list of used iprint flags: xh -P list
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
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: prtevery
    ! customize printouts of two-dimentional arrays
    ! incrni  -  parameters of pmtx routine used to print two-dimensional arrays 
    ! incrmi	 row-wise. Every incrni row and every incrmi column is printed
    
    call inInt(incrni)
    call inInt(incrmu)
    return
  end subroutine read_prtevery

  subroutine read_scf
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: scf
    !        
    !   maxscf  - max. no. of scf iterations
    !   saveFuncs - no. of iterations between backups 
    !   ienterm - if max. error in orbital energy is less than 
    !             1/10**ienterm than scf iterations are terminated
    !   inoterm - if max. error in orbital norm is less than 
    !             1/10**inoterm than scf iterations are terminated
    !   verboseLevel - level of output during scf process
    !       =1 at every nobackup SCF iteration only total energy is printed
    !       =2 at every SCF iteration orbital energy, orbital energy
    !          and normalization convergence of and non-orthogonality for
    !          the worst converged orbitals are printed
    !       =3 at every SCF iteration orbital energy, orbital energy
    !          and normalization convergence and non-orthogonality for
    !          each orbital are printed
    !       =4 additionaly, at every SCF iteration the orbital and normalization
    !          convergence rates for the worst converged orbital are printed
    
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxscf = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          saveScfData = itmp
          call inInt(itmp)
          if (itmp.ne.inpiexit) then
             ienterm = itmp
             call inInt(itmp)
             if (itmp.ne.inpiexit) then
                inoterm = itmp
                call inInt(itmp)
                if (itmp.ne.inpiexit) then
                   verboseLevel=itmp
                endif
             endif
          endif
       endif
    endif
    return
  end subroutine read_scf

  subroutine read_scforder
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: scforder
    do iorb=1,norb
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          iscforder(iorb)=itmp
       else
          write(iout6,'(/2x,"Error: missing or incorrect entry - see User''s guide."/)')
          stop 'read_scforder'
       endif
    enddo
    return
  end subroutine read_scforder

  subroutine read_sor
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: sor
    !     maxsororb(2) - maximal number of micro SOR iterations during relaxation
    !                    of an orbital in a single SCF cycle
    !     maxsorpot(2) - maximal number of micro SOR iterations during relaxation
    !                    of a potential in a single SCF cycle
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsororb(2) = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsorpot(2) = itmp
       endif
    endif
    return
  end subroutine read_sor


  subroutine read_sor4orb
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: sor4orb
    !     maxsororb(2) - maximal number of micro SOR iterations during relaxation
    !                    of an orbital in a single SCF cycle
    !     maxsororb(1) - maximal number of macro SOR iterations during relaxation
    !                    of an orbital in a single SCF cycle, i.e. every orbital
    !                    undergoes maxsor1*maxsor2 SOR iterations in a single SCF 
    !                    iteration
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsororb(2) = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsororb(1) = itmp 
       endif
    endif
    return
  end subroutine read_sor4orb

  subroutine read_sor4pot
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! label: sor4pot
    !     maxsorpot(2) - maximal number of micro SOR iterations during relaxation
    !                    of every potential in a single SCF cycle
    !     maxsorpot(1) - maximal number of macro SOR iterations during relaxation
    !                    of every potential in a single SCF cycle, i.e. every Coulomb
    !                    and exchange potential undergoes maxsor1*maxsor2 SOR iterations
    !                    in a single SCF iteration
    !                   
    call inInt(itmp)
    if (itmp.ne.inpiexit) then
       maxsorpot(2) = itmp
       call inInt(itmp)
       if (itmp.ne.inpiexit) then
          maxsorpot(1) = itmp 
       endif
    endif
    return
  end subroutine read_sor4pot
  
  subroutine read_stop(ni_t,mu_t,no_t,nons_t)
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    
    implicit none
    integer (KIND=IPREC) :: ni_t,mu_t,no_t,nons_t
    
    if (icompLExp.ne.0.and.icompLAdd.eq.0) then
       write(iout6,'(/2x,"Error: label LCAO is missing ...")')
       write(iout6,'( 2x,"It is needed by _ORBPOT hydrogen_ card ..."/)')
       stop 'read_stop'
    endif
    
    ! normalize mixing coefficients
    do iorb=1,norb
       co12=abs(co1(iorb))+abs(co2(iorb))
       if (abs(co12)>epsilon(zero)) then
          co1(iorb)=co1(iorb)/co12
          co2(iorb)=co2(iorb)/co12
       endif
    enddo
    
    iform_t = iform
    ni_t    = nni
    mu_t    = nmutot
    no_t=norb
    nons_t=0
    do iorb=1,norb
       if (mgx(3,iorb).ne.0) nons_t=nons_t+1
    enddo
    
    write(iout6,*) '... end of input data  ...'
    write(iout6,*)'        '
    return
  end subroutine read_stop

  subroutine read_tail
    use params
    use commons
    
    implicit none
    ! alphaf - if alphaf is not zero a DFT energy functional is used
    
    ltail=.true.
  end subroutine read_tail


  
  
  subroutine read_xalpha
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use data4II
    implicit none
    ! alphaf - if alphaf is not zero a DFT energy functional is used
    
    ftmp=0.0_PREC
    call inFloat(ftmp)
    if (abs(ftmp)>epsilon(zero)) then
       alphaf = ftmp
    endif
    return
  end subroutine read_xalpha

  ! ### nmucalc ###
  !
  !     Calculates nmu for given nmu and practical infinity values (a
  !     single grid case is assumed) and adjusts nmu to be of the form
  !     6k+1 and 5k+5+1 (if mcsor is chosen)
  !
  integer function nmucalc(n)
    use params
    use discrete
    use doSCF
    use commons
    implicit none

    integer (KIND=IPREC) :: ig,k5,k6,k10,n,nk,ntemp

    real (PREC) :: hnu,xmi0
    parameter (k5=5, k6=6, k10=10)

    ! Since in this version of the program the 7 point integration
    ! formula is used the number of mesh points in the nu and mu
    ! variables must be of the form 6*n+1.

    ! Another restriction has to be imposed on the number of points in the
    ! mu variable if the mcsor scheme is to be used: nmu must be
    ! of the form kN+k+1, where k=iorder/2+1 for each grid (in this program
    ! iorder, i.e. the order of discretization is 8).

    ! mxnmu --  maximum no. of grid points in mu variable
    ! hni - step in ni variable
    ! hmu - step in mu variable

    ig=1

    if (n.eq.0) then

       ! determine step size in nmu variable so as hmu is approximately
       ! equal to hni

       hnu=pii/dble(nni-1)

       xmi0=2.0_PREC*rgrid(ig)/r
       xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
       mxnmu=NINT(xmi0/hnu+1.0_PREC)
       nmu(ig)=mxnmu
       n=nmu(ig)
    endif

    ! adjust, if necessary, the number of grid points in mu variable
12  nk=(n-1)/k6
    if (k6*nk.ne.n-1) then
       n=n-1
       goto 12
    endif

    ntemp=n

22  continue
    if (meshOrdering/="middle".or.meshOrdering/="middlemt") then
       nk=(n-6)/k5
       if (k5*nk.ne.n-k5-1) then
          n=n-1
          goto 22
       endif
    else
       nk=(n-3)/k10
       if (k10*nk.ne.n-3) then
          n=n-1
          goto 22
       endif
    endif
    if (n.gt.1.and.n.ne.ntemp) goto 12

    nmu(ig)=n
    mxnmu=nmu(ig)

    ! calculate hmu(1) for the adjusted value of mxnmu
    xmi0=2.0_PREC*rgrid(ig)/r
    xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
    hmu(ig)=xmi0/dble(mxnmu-1)

    nmucalc=mxnmu

  end function nmucalc

  ! ### nnucalc ###
  !
  !     Adjust nnu to be of the form 6k+1 and 5k+5+1 (if mcsor is chosen)
  !
  integer function nnucalc(n)
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: k5,k6,k10,n,nk,ntemp

    parameter (k5=5, k6=6, k10=10)

    ! Since in this version of the program the 7 point integration formula
    ! is used the number of mesh points in the nu and mu variables must be
    ! of the form 6*n+1.

    ! Another restriction has to be imposed on the number of points in the
    ! mu variable if the mcsor scheme is to be used: nmu must be of the
    ! form kN+k+1, where k=iorder/2+1 for each grid (in this program
    ! iorder, i.e. the order of discretization is 8).

    ! Adjust, if necessary, the number of grid points in nu variable.

12  nk=(n-1)/k6
    if (k6*nk.ne.n-1) then
       n=n-1
       goto 12
    endif

    ntemp=n
22  continue
    if (meshOrdering/="middle".or.meshOrdering/="middlemt") then
       nk=(n-6)/k5
       if (k5*nk.ne.n-k5-1) then
          n=n-1
          goto 22
       endif
    else
       nk=(n-3)/k10
       if (k10*nk.ne.n-3) then
          n=n-1
          goto 22
       endif
    endif
    
    if (n.gt.1.and.n.ne.ntemp) goto 12
    if (n.eq.1) then
       write(*,'("Error: cannot find commensurate grid size in nu."/)')
       stop 'nnucalc'
    endif

    nnucalc=n
    
  end function nnucalc

end module inputInterface
