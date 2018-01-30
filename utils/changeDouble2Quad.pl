#!/usr/bin/perl -w
################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################

# This script modifies all source files as follows:

# integer*4                  --> integer*8
# implicit integer*4         --> implicit integer*8
# real*8                     --> real*16
# implicit real*8            --> implicit real*16
# dble                       --> qext
# double precision constants --> quadrupole precision constants
#  
# In some files (see below) a specific treatment is needed.

# It may happen that routine getCpuTime will need an individual treatment.

open (ALLFILES, "ls *.f *.c *.raw *.inc|");

@allfiles=<ALLFILES>;
$allfiles_nb=@allfiles-1;

for ($i=0; $i<=$#allfiles; $i++) {
    open(FILE,"cat $allfiles[$i] |") || die "cannot open ", $allfiles[$i];
    @file=<FILE>;
    $file_nb=@file-1;
    close(FILE);
    open(FILE,">$allfiles[$i]") || die "cannot open ", $allfiles[$i];
    for ($j=0; $j<=$file_nb; $j++) {
	$t=$file[$j];

	if ($allfiles[$i] =~ /rheader/) {

	    if ( $t =~ /integer\*4\s+i4tmp1,i4tmp2|real\*8\s+r8tmp1,r8tmp2,r8tmp\(10\)/) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /rfun/) {
	    if ( $t =~ /integer\*4\s+i4tmp\(9750\)/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /integer\*8\s+i8tmp\(9750\)/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /real\*8\s+wk8\(\*\)/ ) {
		print (FILE $t);
		next;
	    }
	}


	if ($allfiles[$i] =~ /reada8/) {
	    if ( $t =~ /real\*8\s+a\(ndim\)/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /real\*8 wk8\(\*\)/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /wtdisk/) {
#	    if ( $t =~ /integer\*4\s+i4tmp/ ) {
#		print (FILE $t);
#		next;
#	    }
#	    if ( $t =~ /real\*8\s+r8tmp/ ) {
#		print (FILE $t);
#		next;
#	    }
	}


	if ($allfiles[$i] =~ /wtdisk32/) {
	    if ( $t =~ /integer\*4\s+i4tmp1/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /real\*8\s+r8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /wtdisk64/) {
	    if ( $t =~ /integer\*8\s+i8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /real\*8\s+r8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /wtdisk128/) {
	    if ( $t =~ /integer\*8\s+i8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	    if ( $t =~ /real\*16\s+r8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	}

        if ($allfiles[$i] =~ /igetrealtype/) {
            if ( $t =~ /integer\*4\s+i4tmp1/ ) {
                print (FILE $t);
                next;
            }

            if ( $t =~ /integer\*8\s+i8tmp1/ ) {
                print (FILE $t);
                next;
            }

            if ( $t =~ /real\*8\s+r8tmp1/ ) {
                print (FILE $t);
                next;
            }
        }

	$t =~ s/_REAL_/real\*16/;

	$t =~ s/      real\*8/      real*16/;
	$t =~ s/implicit real\*8/implicit real\*16/;

	$t =~ s/      integer\*4/      integer*8/;
	$t =~ s/implicit integer\*4/implicit integer\*8/;

	$t =~ s/(\d?\.0?\d*?)d(\-?\+?\d+)/$1q$2/g;

	$t =~ s/dble\(/qext\(/g;

# comment out some functions if they do not have their 64-bit
# alternatives (compiler dependent)

	$t =~ s/(\s+call flush)/c$1/;
        $t =~ s/(\s+stime)/c$1/;
        $t =~ s/(\s+call ctime)/c$1/;

	print (FILE $t);
    }
    close(FILE);
};
#


