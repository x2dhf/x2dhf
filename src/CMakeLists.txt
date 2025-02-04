set(X2DHF_F_SOURCES
  dgamit.F
  card.f90
  checkOrbSym.f90
  coulombExchangePotentials.f90
  data4II.f90
  dateTime.f90
  dftvxc.f90
  discret.f90
  dointerp.f90
  elocc.f90
  exchContribs.f90
  libxc_funcs.f90
  memory.f90
  nabla.f90
  orbAsympt.f90
  params.f90
  printUtils.f90
  scfshr.f90
  solver.f90
  sort.f90
  utils.f90
  blas.F90
  commons.F90
  detect.F90
  dftexc.F90
  diskInterface.F90
  diskInterfaceMisc.F90
  doSCF.F90
  dumpDataToDisk.F90
  fock.F90
  initCBAllocArrays.F90
  initOrbitalsPotentials.F90
  initVariables.F90
  inout.F90
  inputInterface.F90
  integrals.F90
  lagrangeMultipliers.F90
  libxc_master_f90.F90
  machine.F90
  normOrtho.F90
  prepGauss.F90
  printCase.F90
  relaxOrbs.F90
  relaxPots.F90
  sapPotential.F90
  scfUtils.F90
  sharedMemory.F90
  sormcsor.F90
  storedOrbitals.F90
  summary.F90
  totalEnergy.F90
  totalEnergyLXC.F90
  x2dhf.F90
)

set(X2DHF_PTHREADS_SOURCES
  c/sorpt.c
  c/mcsor_pthread.c
  c/coulExch_pthread.c
  c/relax_single_pot_pthread.c
  c/mcsor_single_colour_pthread.c
)

set(X2DHF_TPOOL_SOURCES
  c/sorpt.c
  c/mcsor_tpool.c
  c/coulExch_tpool.c
  c/relax_single_pot_tpool.c
  c/mcsor_single_colour_tpool.c
)

if(ENABLE_PTHREADS)
  if(USE_PTHREAD)
    set(X2DHF_SOURCES ${X2DHF_F_SOURCES} ${X2DHF_PTHREADS_SOURCES})
  else()
    set(X2DHF_SOURCES ${X2DHF_F_SOURCES} ${X2DHF_TPOOL_SOURCES})
  endif()
else()
  set(X2DHF_SOURCES ${X2DHF_F_SOURCES})
endif()

# Main program
add_executable(x2dhf ${X2DHF_SOURCES})
target_link_libraries(x2dhf -lxc)

# Install
install (TARGETS x2dhf DESTINATION bin)
