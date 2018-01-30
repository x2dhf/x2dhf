#!/usr/bin/perl 

################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################

use Getopt::Long;

GetOptions ("help"      => \$help,
	    "homo"      => \$homo,
            "field=f"   => \$field,
            "deltae=f"  => \$deltaE,
            "deltaq=f"  => \$deltaQ,
            "deltaq1=f" => \$deltaQ1,
            "deltaq2=f" => \$deltaQ2);

if (  defined($help) ) {
    print <<'EOF';

This script facilitates the calculation of polarizabilities and hyperpolarizabilities from
a series of calculations performed by means of the 2DHF program. For heteronuclear systems
calculations must be performed for -2F, -F, 0, +F and +2F values of external static
electric field. For homonuclear molecules total energies and multipole moments for +F and
+2F field strengts may be obtained from symmetry considerations (see --homo opition). Edit
the script to specify the names of output files. In order to assess the accuracy of
calculated electric properties relative errors of total energies and dipole moments must
also be provided.

Usage:

   runelprop.pl [--help] [--homo] [--field=F] [--deltae=f] [--deltaq=f] [--deltaq1=f] [--deltaq2=f] < file


Options:

 --homo 
      assume a homonuclear case, i.e. use output data for only -2F, -F and 0 field
      strengths. If results for the +F and +2F values are also available (and
      heteronuclear case is assumed) they are used to assess the accuracy of total
      energies and multipole moments (see elprop.pl).

 --field 
      the value of external field strength must be given if it is not available in the
      second of the output file (relevant for the 2DHF programs prior to ver. September
      2002).

 --deltae 
      the relative error of total energies

 --deltaq
      the relative error of multipole moments 

 --deltaq1
      the relative error of dipole moments 

 --deltaq2
      the relative error of quadrupole moments 

   file
      a text file with names of the 2DHF output files with results for -2F, -F, 0, +f, +2F field strengths.
      Each name should appear in a separate line. Empty lines or lines beginning with a hash ('#') are skipped. 

EOF
exit;
} 


if ( defined($homo) ) {
    $homo="yes";
} else {
    $homo="no";
}

if ( defined($deltaE) ) {
    $de="--deltae=$deltaE";
} else {
    $de="";
}

if ( defined($deltaQ) ) {
    $dq="--deltaq=$deltaQ";
} else {
    $dq="";
}

if ( defined($deltaQ1) ) {
    $dq1="--deltaq1=$deltaQ1";
} else {
    $dq1="";
}

if ( defined($deltaQ2) ) {
    $dq2="--deltaq2=$deltaQ2";
} else {
    $dq2="";
}

if ( defined($field) ) {
    $f="--field=$field";
} else {
    $f="";
}

print $ARGV[0],"\n";

$pwd=`pwd`;
chomp($pwd);

$i=-1;
while (<>) {
    chomp;
    next if ( /^\#+/ );
    next if ( /^$/ );
    $i++;
    s/^\s+//;                          # filename cannot have leading spaces
    $fn[$i]=${pwd} . "/" . $_;
#    print $i,"  ",$fn[$i],"\n";
}

if ( $homo eq "yes" ) { 
    system ("elprop.pl $fn[0] $fn[1] $fn[2]        --homo $de $dq $dq1 $dq2 $f \n" ); 
} else {
    system ("elprop.pl $fn[0] $fn[1] $fn[2] $fn[3] $fn[4] $de $dq $dq1 $dq2 $f \n" ); 
}


exit;

=pod
# example of an input file 
# any comments
n2-645-250-2m3-final.lst
n2-645-250-1m3-final.lst
n2-645-250-0-final.lst
n2-645-250-1p3-final.lst
n2-645-250-2p3-final.lst

# empty lines are also skipped
#
=cut




