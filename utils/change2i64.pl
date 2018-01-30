#!/usr/bin/perl -w
################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################

# This script modifies (with some exceptions) all source files as follows:

#    integer*4           --> integer*8
#    implicit integer*4  --> implicit integer*8


open (ALLFILES, "ls *.f *.c *.raw *.inc|");

@allfiles=<ALLFILES>;
$allfiles_nb=@allfiles-1;

for ($i=0; $i<=$allfiles_nb; $i++) {
    open(FILE,"cat $allfiles[$i] |") || die "cannot open ", $allfiles[$i];
    @file=<FILE>;
    $file_nb=@file-1;
    close(FILE);
    open(FILE,">$allfiles[$i]") || die "cannot open ", $allfiles[$i];
    for ($j=0; $j<=$file_nb; $j++) {
	$t=$file[$j];

	if ($allfiles[$i] =~ /rheader/) {
	    if ( $t =~ /integer\*4 i4tmp1,i4tmp2/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /rfun/) {
	    if ( $t =~ /integer\*4\s+i4tmp\(9750\)/ ) {
		print (FILE $t);
		next;
	    }
	}


	if ($allfiles[$i] =~ /wtdiski32/) {
	    if ( $t =~ /integer\*4\s+i4tmp1/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /wtdiski64/) {
	    if ( $t =~ /integer\*8\s+i8tmp1/ ) {
		print (FILE $t);
		next;
	    }
	}

	if ($allfiles[$i] =~ /wtdiski128/) {
	    if ( $t =~ /integer\*8\s+i8tmp1/ ) {
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
        }



	$t =~ s/      integer\*4/      integer*8/;
	$t =~ s/implicit integer\*4/implicit integer\*8/;

# comment out some functions if they do not have their 64-bit
# alternatives (compiler dependent)

# 	$t =~ s/(\s+call flush)/c$1/;
#         $t =~ s/(\s+stime)/c$1/;
#         $t =~ s/(\s+call ctime)/c$1/;

	print (FILE $t);
    }
    close(FILE);
};




