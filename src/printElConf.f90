module printElConf_m
  implicit none
contains
  subroutine printElConf
    use params
    use discret
    use commons8

    implicit none

    integer :: i,iorb,iorbr,iput,izz1,izz2,nsigma,nsigmag,nsigmau,npi,ndelta,nphi,nspins


    izz1=nint(z1)
    izz2=nint(z2)

    write(iout6,*)
    write(iout6,*) '  electronic configuration:'
    write(iout6,*)

    nsigma=0
    nsigmag=0
    nsigmau=0
    npi   =0
    ndelta=0
    nphi  =0

    do iorb=1,norb
       if (orbsym(iorb).eq.'sigma') nsigma=nsigma+1
       !        if (orbsym(iorb).eq.'sigma g') nsigmag=nsigmag+1
       !        if (orbsym(iorb).eq.'sigma u') nsigmau=nsigmau+1
       if (orbsym(iorb).eq.'pi'   ) npi   =npi   +1
       if (orbsym(iorb).eq.'delta') ndelta=ndelta+1
       if (orbsym(iorb).eq.'phi'  ) nphi  =nphi  +1
    enddo

    iorbr=0
    do iorb=1,norb
       if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma') iorbr=nsigma
       !        if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma g') iorbr=nsigmag
       !        if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma u') iorbr=nsigmau
       if (iorbr.eq.0.and.orbsym(iorb).eq.'pi') iorbr=npi
       if (iorbr.eq.0.and.orbsym(iorb).eq.'delta') iorbr=ndelta
       if (iorbr.eq.0.and.orbsym(iorb).eq.'phi') iorbr=nphi
       nspins=4
       if (orbsym(iorb).eq.'sigma') nspins=2
       iput=4*(iorb-1)
       write(iout6,1020) iorbr, orbsym(iorb),gut(iorb),(spin(i+iput),i=1,nspins)
       iorbr=iorbr-1
    enddo

    !     number of electrons and total charge

    write(iout6,*)
    write(iout6,1010) izz1+izz2-nel,nel,norb,norb,nexch

    if (ini.eq.1.or.ini.eq.4) then
       write(iout6,*)
       write(iout6,*) '  LCAO via hydrogenic functions:'
       write(iout6,*)
       write(iout6,1030)

       iorbr=0
       do iorb=1,norb
          if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma') iorbr=nsigma
          if (iorbr.eq.0.and.orbsym(iorb).eq.'pi') iorbr=npi
          if (iorbr.eq.0.and.orbsym(iorb).eq.'delta') iorbr=ndelta
          if (iorbr.eq.0.and.orbsym(iorb).eq.'phi') iorbr=nphi

          iput=4*(iorb-1)
          write(iout6,1032) iorbr, orbsym(iorb),gut(iorb),mgx(1,iorb),mgx(2,iorb),eza1(iorb),co1(iorb),&
               mgx(4,iorb),mgx(5,iorb),eza2(iorb),co2(iorb)
          iorbr=iorbr-1
       enddo
    endif

1010 format(/10x,'total charge            =',i4 &
         /10x,'number of',                    &
         /14x,'electrons           =',i4      &
         /14x,'orbitals            =',i4      &
         /14x,'Coulomb potentials  =',i4      &
         /14x,'exchange potentials =',i4)
    !  1010 format(10x,'number of electrons =',i5
    !      &     /10x,'total charge        =',i5
    !      &     /10x,'number of orbitals, Coulomb and exchange potentials',
    !      &     3i4)

1020 format(10x,i2,2x,a8,a1,2x,4a2)
1030 format(10x,' orbital           n1 l1   Z1    c1','       n2 l2   Z2    c2 '/)
1032 format(10x,i2,2x,a8,a1,5x,2i3,2f6.2,5x,2i3,2f6.2)

  end subroutine printElConf
end module printElConf_m
