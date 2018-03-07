! ### inInt ###

!    Reads an integer from the array IA, starting at IA(istrt(jrec)) and
!    continuing for inumb(jrec)) elements. Plus signs are ignored, the
!    answer is accumulated in JBUF.

module inInt_m
  implicit none
contains
  subroutine inInt(jbuf)
    use params
    use card

    implicit none

    integer :: i,ifact,inpiexit,ist,j,jbuf,n,nstrt

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
    use card

    implicit none

    integer :: i,ifact,inpiexit,ist,j,jbuf,n,nstrt

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
end module inInt_m
