AC_INIT([x2dhf], [1.0])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE



AC_LANG_FORTRAN77
AC_PROG_F77(gfortran)
AC_F77_LIBRARY_LDFLAGS

AC_PROG_CC(cc)

#AC_PROG_FC([ifort ifc])

#AC_OUTPUT(Makefile src/Makefile)

AC_CONFIG_SRCDIR([src/x2dhf.f90])
AC_CONFIG_FILES([Makefile
 	         docs/Makefile
	         examples/Makefile
	         src/Makefile
		 utils/Makefile
		 tests/Makefile])

AC_OUTPUT
