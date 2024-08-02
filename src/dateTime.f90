! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2024  Jacek Kobus 

module dateTime
  use params, only : IPREC, PREC
contains
  ! ### getDateTime ###
  !
  !     Returns current date and time as a string.
  !     This routine works for gfortran, g77 and ifort compilers
  !     Comment the lines
  !        stime=time()
  !        call ctime(stime,datetime)
  !     or replace them with appropriate equivalents
  !
  subroutine getDateTime(curDateTime)
    implicit none
    integer (KIND=IPREC) :: i
    character*80 :: curDateTime,datetimex
    character*8  :: date
    character*10 :: time
    character*1, dimension(80) :: str,times
    character*1, dimension(8) :: dates

    equivalence(datetimex,str(1))

    equivalence(date,dates(1))
    equivalence(time,times(1))

    do i=1,80
       str(i)=' '
    enddo

    ! stime=time()
    ! call ctime(stime,datetime)
    ! date_time works with g77, gfotran and ifort
    call date_and_time (date,time)

    ! copy data from date and time character variables into datetime one
    ! inserting appropriate separators

    str( 1)=dates(1)
    str( 2)=dates(2)
    str( 3)=dates(3)
    str( 4)=dates(4)
    str( 5)='/'
    str( 6)=dates(5)
    str( 7)=dates(6)
    str( 8)='/'
    str( 9)=dates(7)
    str(10)=dates(8)

    str(11)=' '
    str(12)=' '

    str(13)=times(1)
    str(14)=times(2)
    str(15)=':'

    str(16)=times(3)
    str(17)=times(4)
    str(18)=':'
    str(19)=times(5)
    str(20)=times(6)
    str(21)=times(7)
    str(22)=times(8)
    str(23)=times(9)


    curDateTime=datetimex

  end subroutine getdatetime

  ! ### getCpuTime ###
  !
  !     Returns cpu time in seconds with microseconds accuracy
  !
  subroutine getCpuTime(t)
    use params

    implicit none
    real (PREC) t
    call cpu_time (t)

  end subroutine getCpuTime

  ! ### getRealTime ###
  !
  !     Returns the counts of system clock in 1/countRate seconds
  !
  subroutine getRealTime(count,countRate)
    use params
    integer (KIND=IPREC) :: count,countRate,countMax
    call system_clock(count, countRate, countMax)
  end subroutine getRealTime

end module dateTime
