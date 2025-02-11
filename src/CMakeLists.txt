set(X2DHF_F_SOURCES
  dgamit.F
  blas.F90
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
  memory.f90
  nabla.f90
  orbAsympt.f90
  params.F90
  printUtils.f90
  scfshr.f90
  solver.f90
  sort.f90
  utils.f90
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
  machine.F90
  normOrtho.F90
  prepGauss.F90
  printCase.F90
  relaxOrbs.F90
  relaxPots.F90
  sapPotential.f90
  scfUtils.F90
  sharedMemory.F90
  sormcsor.F90
  storedOrbitals.F90
  summary.F90
  totalEnergy.F90
  totalEnergyLXC.F90
  x2dhf.F90
)

# Add libxc files to list
if(LIBXC_FOUND)
  set(X2DHF_F_SOURCES ${X2DHF_F_SOURCES}
    libxc_funcs.F90
    libxc_master.F90)
  include_directories(${LIBXC_INCLUDE_DIR})
  add_compile_definitions(-DLIBXC)
endif()

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
  # Include the present directory in the search path
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  if(USE_PTHREAD)
    set(X2DHF_SOURCES ${X2DHF_F_SOURCES} ${X2DHF_PTHREADS_SOURCES})
    add_compile_definitions(-DPTHREAD)
  else()
    set(X2DHF_SOURCES ${X2DHF_F_SOURCES} ${X2DHF_TPOOL_SOURCES})
    add_compile_definitions(-DTPOOL)
  endif()
else()
  set(X2DHF_SOURCES ${X2DHF_F_SOURCES})
endif()

if(BLAS_FOUND)
  # BLAS_FOUND signals external BLAS library found
  add_compile_definitions(-DBLAS)
  message("Using external BLAS")
else()
  # Use the internal library
  message("Using internal BLAS")
endif()

# Main program
add_executable(x2dhf ${X2DHF_SOURCES})
if(LIBXC_FOUND)
  target_link_libraries(x2dhf ${LIBXC_LIBRARY})
endif()
if(BLAS_FOUND)
  target_link_libraries(x2dhf ${BLAS_LIBRARIES})
endif()

# Install
install (TARGETS x2dhf DESTINATION bin)
