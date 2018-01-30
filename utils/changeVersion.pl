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


open (ALLFILES, "ls *.f *.c *.raw *.inc *.pl|");

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

	# $t =~ s/Copyright \(C\) 1996-2008 Jacek/Copyright \(C\) 1996-2009 Jacek/;
	# $t =~ s/Copyright \(C\) 1996-2006 Jacek/Copyright \(C\) 1996-2009 Jacek/;
	# $t =~ s/Copyright \(C\) 1997-2006 Jacek/Copyright \(C\) 1997-2009 Jacek/;
	# $t =~ s/Copyright \(C\) 1996,2008 Jacek/Copyright \(C\) 1996,2009 Jacek/;
	# $t =~ s/Copyright \(C\) 2008 Jacek/Copyright \(C\) 2009 Jacek/;
	# $t =~ s/Copyright \(C\) 2005-2008 Jacek/Copyright \(C\) 2005-2009 Jacek/;
	$t =~ s/version 2.0/version 2.1/;

	print (FILE $t);
    }
    close(FILE);
};




