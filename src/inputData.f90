! ### inputData ###
! Handles the input to x2DHF.  
module inputData_m
  implicit none
contains
  subroutine inputData(ni_t,mu_t,no_t,nons_t)
    use params
    use discret
    use scf
    use solver
    use commons8

    use input_m
    use readLables_m
    use inCard_m
    use inFloat_m
    use inInt_m
    use inStr_m
    use nmucalc_m
    use nnucalc_m
    use xc_f90_types_m
    use xc_f90_lib_m

    implicit none

    integer :: ni_t,mu_t,no_t,nons_t
  
    character*8 readLabel
    character*20 label4m
    
    do i=1,10
       cdftex(i)=cdftext(i)
       cdftcorr(i)=cdftcorrt(i)
    enddo
    
    do i=1,10
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

    call inCardh(header)
    call inStr(clabel)

    if (clabel.ne.'title'.and.clabel.ne.'TITLE') then
       write(iout6,'(/,"Error: label TITLE is missing. "/)')
       stop 'inputData'
    endif

    icompLEnc=icompLEnc+1

    clabel=""
    do while (clabel.ne.'stop')
       call inCard
       call inStr(clabel)
       call checkLabel
       if (clabel.eq.'break') call read_break
       if (clabel.eq.'config') call read_config
       if (clabel.eq.'conv') call read_conv
       if (clabel.eq.'debug') call read_debug
       if (clabel.eq.'dft') call read_dft
       if (clabel.eq.'exchio') call read_exchio
       if (clabel.eq.'fefield') call read_fefield
       if (clabel.eq.'fermi') call read_fermi
       if (clabel.eq.'fix') call read_fix
       if (clabel.eq.'fixorb') call read_fixorb
       if (clabel.eq.'fixcoul') call read_fixcoul
       if (clabel.eq.'fixexch') call read_fixexch
       if (clabel.eq.'gauss') call read_gauss
       if (clabel.eq.'grid') call read_grid
       if (clabel.eq.'homo') call read_homo
       if (clabel.eq.'inout') call read_inout
       if (clabel.eq.'interp') call read_interp    
       if (clabel.eq.'lagra') call read_lagra
       if (clabel.eq.'lcao') call read_lcao
       if (clabel.eq.'mcsor') call read_mcsor
       if (clabel.eq.'method') call read_method
       if (clabel.eq.'multipol') call read_multipol
       if (clabel.eq.'nuclei') call read_nuclei    
       if (clabel.eq.'omega') call read_omega
       if (clabel.eq.'orbpot') call read_orbpot
       if (clabel.eq.'order') call read_order
       if (clabel.eq.'potgsz') call read_potgsz
       if (clabel.eq.'potgszg') call read_potgszg
       if (clabel.eq.'potharm') call read_potharm
       if (clabel.eq.'pothook') call read_pothook
       if (clabel.eq.'poth3') call read_poth3    
       if (clabel.eq.'potkh') call read_potkh
       if (clabel.eq.'print') call read_print
       if (clabel.eq.'prtevery') call read_prtevery
       if (clabel.eq.'scf') call read_scf
       if (clabel.eq.'scfexch') call read_scfexch
       if (clabel.eq.'sor') call read_sor
       if (clabel.eq.'sormcsor') call read_sormcsor            
       if (clabel.eq.'stop') call read_stop(ni_t,mu_t,no_t,nons_t)
       if (clabel.eq.'xalpha') call read_xalpha
    enddo

    if (imethod.eq.3) then
       if ( nint(z1).ge.nint(z2) )  then
          alphaf=alphaopt(nint(z1))
       else
          alphaf=alphaopt(nint(z2))
       endif
    endif

    return

  end subroutine inputData
end module inputData_m
