#!/usr/bin/perl  

# SPDX-License-Identifier: GPL-2.0+
# Copyright (C) 2002-2024 Jacek Kobus 

# This script facilitates the calculation of polarizabilities and
# hyperpolarizabilities from a series of calculations performed by means of
# the X2DHF program.

# Usage:

#   elprop.pl [--field=F] file_m2 file_m1 file_p0 file_p1 file_p2 [--deltae=f] [--deltaq=f]
#   elprop.pl [--field=F] --homo file_m2 file_m1 file_p0 [--deltae=f] [--deltaq=f]  
 
# where file_i are output files from the x2dhf program corresponding to
# calculations with the external finite electric field strengths -2F, -F,
# 0, +F and +2F, respectively. If possible, the F value is extracted from
# the second of these files. Otherwise it must be defined using the long
# option 'field'. For homonuclear systems calculations for the -2F, -F and
# 0 field values need only be performed. However, if available, results for
# the +F and +2F field strengths are used to assess the accuracy of total
# energies and multipole moments (see below).

# In order to assess the accuracy of calculated electric properties provide
# relative errors of total energies (deltae) and dipole moments (deltaq).

# The derivatives are evaluated by means of the 3 and 5-point central
# differences formulae in order to assess the accuracy of calculations.

# For further details see 

# J.Kobus, D.Moncrieff, and S.Wilson, A comparison of the electric moments
# obtained from finite basis set and finite difference {Hartree-Fock}
# calculations for diatomic molecules, Phys. Rev. A 62 (2000), 062503/1--9.

# J.Kobus, D.Moncrieff, and S.Wilson, Comparison of the polarizabilities and
# hyperpolarizabilities obtained from finite basis set and finite field, finite
# difference {Hartree-Fock} calculations for diatomic molecules, J. Phys. B:
# At. Mol. Opt. Phys. 34 (2001) 5127-5143.

#############################################################################################


use Getopt::Long;


#$debug1=0; # set to a nozero value for additional printouts
#$debug2=0; # set to a nozero value to monitor the dependence of alpha,
           # beta and gamma (h^4) values on the accuracy of the dipole moment

#$debug3=0; # set to a nonzero value to monitor accuracy of beta parameter
#$debug4=0; # set to a nonzero value to print total multipole moments
#$debug5=0; # set to a nonzero value to print electronic multipole moments

#$debug6=0; 

GetOptions ("homo"     => \$homo,
            "field=f"  => \$field,
            "deltae=f" => \$deltaE,
            "deltaq=f" => \$deltaQ,
            "deltaq1=f" => \$deltaQ1,
            "deltaq2=f" => \$deltaQ2,
            "debug1"    => \$debug1,
	    "debug2"    => \$debug2,
	    "debug3"    => \$debug3,
	    "debug4"    => \$debug4,
	    "debug5"    => \$debug5,
            "debug6"    => \$debug6);

if ( defined($homo) ) {
    $homo="yes";
    $narg=3;
} else {
    $homo="no";
    $narg=5;
}

if ( defined($debug1) ) {
    $debug1=1;
} else {
    $debug1=0;
}

if ( defined($debug2) ) {
    $debug2=1;
} else {
    $debug2=0;
}

if ( defined($debug3) ) {
    $debug3=1;
} else {
    $debug3=0;
}

if ( defined($debug4) ) {
    $debug4=1;
} else {
    $debug4=0;
}

if ( defined($debug5) ) {
    $debug5=1;
} else {
    $debug5=0;
}

if ( defined($debug6) ) {
    $debug6=1;
} else {
    $debug6=0;
}



sub printRelError {
    $text=shift;
    $value=shift;
    printf "%18s %-10.2e\n",$text, $value;
}

sub extract_data {

# examine the five output files (new format, version >= 1-2005) in turn

    for ($iarg=0; $iarg<$narg; $iarg++) {
	
	$x2dhf_out=$ARGV[$iarg]; 
	open(X2DHF_OUT,$x2dhf_out) || die " cannot open: $x2dhf_out";
	if ($debug1) {
	    printf "\n %6s %18s \n","File: ",$x2dhf_out;
	}

	while (<X2DHF_OUT>){

	    if ( /total energy uncertainty due to orbital norms not being equal/ ) {
		$_=<X2DHF_OUT>;


		
		$_ =  /absolute = (.*),\s+relative = (.*)/;
		$absoluteE[$iarg]=$1;
		$relativeE[$iarg]=$2;
		$absoluteE[$iarg]=~ s/\+\/\-(0\.\d\d.*)/\1/;
		$relativeE[$iarg]=~ s/%//;
		$relativeE[$iarg]=~ s/\+\/\-(0\.\d\d.*)/\1/;
		$relativeE[$iarg]=$relativeE[$iarg]/100.0;
#		print "xxxxxxxxxxxxxxxx ",$iarg," ",$absoluteE[$iarg],"\n";
#		print "xxxxxxxxxxxxxxxx ",$iarg," ",$relativeE[$iarg],"\n";
		next;
	    }

	    if ( $mmRelErrFromFiles eq "yes" ) {
		if ( /relative errors of moments due to orbital norms not being equal/ ) {
		    $_=<X2DHF_OUT>; 
		    $_ =  /\s+k=1\s+(.*)/;
		    $deltaQ1[$iarg]=$1;
		    
		    $_=<X2DHF_OUT>; 
		    $_ =  /\s+k=2\s+(.*)/;
		    $deltaQ2[$iarg]=$1;
		    if ($deltaQ1[$iarg] <= 1e-14) { $deltaQ1[$iarg]=1e-14; }
		    if ($deltaQ2[$iarg] <= 1e-14) { $deltaQ2[$iarg]=1e-14; }
#		    $mmRelErrFromFiles="yes";
		    next;
		}
	    }

	    if ( /multipole moments relative to geometrical centre:/ ) {
 		for ($i=1; $i<=2; $i++) {  # skip 2 lines 
 		    $_=<X2DHF_OUT>;
 		}

# extract total dipole moment (au)
		$_=<X2DHF_OUT>;
		@t = split (" ",$_);
		$Q_1tgeom[$iarg]=$t[4];
		if ($debug1) {
		    printf "%33s %25.16e\n","total dipole moment (au)         : ",$Q_1tgeom[$iarg];
		}
	    }

	    if ( $_ =~ /centre\s+atomic\s+weight\s+z/ ) {
# extract atomic masses and the centre of mass

		while (<X2DHF_OUT>) {  # skip some lines
		    last if (/dipole/);
		}

# extract dipole, quadrupole etc moments 
# electronic (au)           Q_1    Q_2  
# electronic (Debye-Ang**k) Q_1da  Q_2da 
# total (au)                Q_1t     Q_2t
# total (Debye-Ang**k)      Q_1tda   Q_2tda 
		
		if ($debug1) {
		    printf "%50s \n","multipole moments relative to centre of mass:";
		}

#		$_=<X2DHF_OUT>;
		@t = split (" ",$_);
		$Q_1[$iarg]=$t[3]; $Q_1t[$iarg]=$t[4];

		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_1da[$iarg]=$t[0]; $Q_1tda[$iarg]=$t[1];

		$_=<X2DHF_OUT>; 		
		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_2[$iarg]=$t[3]; $Q_2t[$iarg]=$t[4];

		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_2da[$iarg]=$t[0]; $Q_2tda[$iarg]=$t[1];

		$_=<X2DHF_OUT>; 		
		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_3[$iarg]=$t[3]; $Q_3t[$iarg]=$t[4];

		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_3da[$iarg]=$t[0]; $Q_3tda[$iarg]=$t[1];

		$_=<X2DHF_OUT>; 		
		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_4[$iarg]=$t[3]; $Q_4t[$iarg]=$t[4];

		$_=<X2DHF_OUT>; 		
		@t = split (" ",$_);
		$Q_4da[$iarg]=$t[0]; $Q_4tda[$iarg]=$t[1];

		if ($debug1) {
		    printf "%42s",  "dipole electronic total (au/Debye): ";
		    printf "%25.16e%25.16e\n", $Q_1[$iarg],$Q_1t[$iarg];
		    printf "%42s%25.16e%25.16e\n", ' ',$Q_1da[$iarg],$Q_1tda[$iarg];

		    printf "\n%42s","quadrupole electronic total (au/Debye): ";
		    printf "%25.16e%25.16e\n", $Q_2[$iarg],$Q_2t[$iarg];
		    printf "%42s%25.16e%25.16e\n", ' ', $Q_2da[$iarg],$Q_2tda[$iarg];

		    printf "\n%42s","octopole electronic total (au/Debye): ";
		    printf "%25.16e%25.16e\n", $Q_3[$iarg],$Q_3t[$iarg];
		    printf "%42s%25.16e%25.16e\n", ' ', $Q_3da[$iarg],$Q_3tda[$iarg];

		    printf "\n%42s","hexadecapole electronic total (au/Debye): ";
		    printf "%25.16e%25.16e\n", $Q_4[$iarg],$Q_4t[$iarg];
		    printf "%42s%25.16e%25.16e\n", ' ', $Q_4da[$iarg],$Q_4tda[$iarg];
		}
	    }
	    
	    if ( $_ =~ / ... saving data to disk .../) {
# extract total energy
		$_=<X2DHF_OUT>;
		$_=<X2DHF_OUT>;
		if ( $_ =~ /total energy/) {
		    ($t[0],$t[1],$t[2]) = split (" ",$_);
		    $t[2] =~ s/D/E/;
		    $E[$iarg]=$t[2];
		    if ($debug1) {
			printf "%18s %26.16e \n","total energy:",$E[$iarg];
		    }
		}
	    }
	    if ( $absoluteE[$iarg] >  $absoluteE ) { $absoluteE=$absoluteE[$iarg]; }
	    if ( $relativeE[$iarg] >  $relativeE ) { $relativeE=$relativeE[$iarg]; }

	}
	close(X2DHF_OUT);
    }
}



# check the number of command line arguments

foreach ($i=0; $i<=$narg-1;$i++) {
    if ( $ARGV[$i] eq '' ) {
	print ("\nError: missing argument(s) \n\n");
	print ("Usage: elprop.pl [--field=f] file_m2 file_m1 file_p0 file_p1 file_p2   \n"); 
	print ("       elprop.pl [--field=f] --homo file_m2 file_m1 file_p0         \n\n"); 
	exit;
    }
}


$Oh2="no"; # do not print alpha_zz and beta_zzz values of O(h^2) accuracy

if ( defined($deltaE) || defined($deltaQ) ) {$debug2=1};

$dbleForm="%8s %23.15e \n";
$quadForm="%8s %41.30e \n";

$dbleForm1="%8s %23.15e %20s %23.15e \n";
$dbleForm2="%8s %23.15e %12s %23.15e  \n";
#$dbleForm2="%8s %23.15e %12s %23.15e %12s %23.15e %12s %23.15e \n";
$quadForm1="%8s %41.30e %20s %41.30e \n";

$form=$dbleForm;
$form2=$dbleForm2;
$DeltaEdefault=1e-14;  
$DeltaQdefault=1e-12;  

if ( ! defined($deltaE) )  { $deltaE =$DeltaEdefault; }
if ( ! defined($deltaQ1) ) { $deltaQ1=$DeltaQdefault; }
if ( ! defined($deltaQ2) ) { $deltaQ2=$DeltaQdefault; }

$mmRelErrFromFiles="yes";


# extract the strength of a field and some other parameters from the
# second of the files

$x2dhf_out=$ARGV[1]; 

open(X2DHF_OUT,"<$x2dhf_out") || die " cannot open: $x2dhf_out";

$F="unset";

while (<X2DHF_OUT>){
    if ( $_ =~ /finite electric field:/) {
	@t=split(" ",$_);
	$F=abs($t[3]);
	last;
    }  
}
close(X2DHF_OUT);

# in case the F value is not included in an output file it can be
# given as the value of the long option field (e.q. elprop.pl --field=0.0001)



if ( $F eq "unset" ) {
    if (defined($field) ) {
	$F=$field;
    } else {
	print ("Error: finite field strength is zero; try using --field option \n");
	exit;
    }
}

# extract some other parameters from the file

$absoluteE=0.0;
$relativeE=0.0;


open(X2DHF_OUT,"<$x2dhf_out") || die " cannot open: $x2dhf_out";
while (<X2DHF_OUT>){
    if ( /nuclei/i) {
# extract Z1 Z2 
#	($t[1],$Z1,$Z2)=split(" ",$_);
        ($t[1],$Z1,$Z2,$R)=split(" ",$_);
        if ( /angstrom/i) {
            $R=$R/0.529177249;
        }
    }

    if ( $_ =~ /start of input data/) {
# # extract the description of a case and the internuclear distance
 	$_=<X2DHF_OUT>;
	$header=$_;
 	@t=split(" ",$_);
	$header=$t[1];
    }

    if ( /\s+total charge\s+=\s+(\-?\d+)/) {
# extract the grid size
	$totalCharge=$1;
#	print $totalCharge,"\n";
    }
    
    
    if ( $_ =~ /Grid:/) {
# extract the grid size
	$_=<X2DHF_OUT>;
	$_ =~ /.*\s+(\d{2,})\s+\(/;
	$nu=$1;
	$_=<X2DHF_OUT>;
	$_ =~ /\s+(\d{2,})\s+\(/;
	$mu=$1;
	$_=<X2DHF_OUT>;
	$_ =~ /\s+(\d+\.\d\d)/;
	$R_inf=$1;
    }

    if ( $_ =~ /.*\s+machine accuracy\s+=\s+(0\.\d{2})e(\-\d{2})/i) {
# extract accuracy of calculations
	$accuracy=$1*10**($2);

	if ( $accuracy < 1e-18 ) {
	    $form=$quadForm;
	    $form1=$quadForm1;
	    $DeltaEdefault=1e-30;  
            $DeltaQdefault=1e-30;  
	} else {
	    $form=$dbleForm;
	    $form1=$dbleForm1;
	    $DeltaEdefault=1e-12;  
            $DeltaQdefault=1e-12;  

	}
    }

    if ( $_ =~ /centre\s+atomic\s+weight\s+z/ ) {
# extract atomic masses and the centre of mass
	$_=<X2DHF_OUT>;
	@t = split (" ",$_);
	$ma=$t[1];
	$_=<X2DHF_OUT>;
	@t = split (" ",$_);
	$mb=$t[1];
	$_=<X2DHF_OUT>;
	@t = split (" ",$_);
	$z=$t[1];
	#if ($debug1) {print $ma,"  ",$mb, "  ",$z,"\n";}
    }
    
    if ( /electronic \(au\/Debye\-Ang\^k\)/ ) {	
	$ver=1; # output in new format
    }


}
close(X2DHF_OUT);

extract_data;


# printf "%12.2e %12.2e %12.2e\n\n", $DeltaE,$DeltaQ1,$DeltaQ2;


if ( $homo eq "yes" ) {
    $Q_2t[4] =  $Q_2t[0];
    $Q_2t[3] =  $Q_2t[1];
    $Q_1t[4] = -$Q_1t[0];
    $Q_1t[3] = -$Q_1t[1];
    $E[4]    =  $E[0];
    $E[3]    =  $E[1];
}


# $DeltaE is the estimated relative error of the total energy.

# $DeltaQ is the (relative) error of the worst converged orbital as indicated by
# its norm.

# For atoms and homonuclear molecules these values can be calculated from differences of
# the total energy and dipole moment values for +/-F.

# If DeltaEQ.data file is present then its content is used to set (or override) the
# default values of the total energy and the dipole and quadrupole moments. Its name must
# be given relative to the propet directory.

if ( -e "../deltaEQ.data" ) {
    open(DEFAULTS,"../deltaEQ.data") || die " cannot open: $!";
    $deltaE=<DEFAULTS>;
    $deltaQ1=<DEFAULTS>;
    $deltaQ2=<DEFAULTS>;
    close(DEFAULTS);
}

# Default values of $DelataE and $DeltaQ are accuracy dependent (see above)

$DeltaE  = $deltaE;
$DeltaE1 = $deltaE;
$DeltaE2 = $deltaE;


$DeltaQ1 = $deltaQ1;
$DeltaQ11= $DeltaQ1; 
$DeltaQ12= $DeltaQ1;

$DeltaQ2 = $deltaQ2;
$DeltaQ21= $DeltaQ2; 
$DeltaQ22= $DeltaQ2; 

#print $DeltaE,"\n";


# if ( $DeltaQ11 > $DeltaQ1 ) { $DeltaQ1=$DeltaQ11 }; 
# if ( $DeltaQ12 > $DeltaQ1 ) { $DeltaQ1=$DeltaQ12 }; 


# if ( $DeltaQ21 > $DeltaQ2 ) { $DeltaQ2=$DeltaQ21 }; 
# if ( $DeltaQ22 > $DeltaQ2 ) { $DeltaQ2=$DeltaQ22 }; 


for ($iarg=0; $iarg<$narg; $iarg++) {
    if ( $relativeE[$iarg] > $DeltaE ) {$DeltaE=$relativeE[$iarg];}

# In case of HeNe (very small Q1(F=0)) relatively large relative errors of Q1(F=0) value
# has to be singled out
 
    if ( $iarg != 2 && $mmRelErrFromFiles eq "yes" ) {
#    if ( $mmRelErrFromFiles eq "yes" ) {
	if ($deltaQ1[$iarg]   > $DeltaQ1 ) {$DeltaQ1=$deltaQ1[$iarg];}
	if ($deltaQ2[$iarg]   > $DeltaQ2 ) {$DeltaQ2=$deltaQ2[$iarg];}
    }
    if ($debug1) {
	printf "%-80s \n         %12.2e %12.2e %12.2e\n\n", $ARGV[$iarg],$relativeE[$iarg],$deltaQ1[$iarg],$deltaQ2[$iarg];
    }
}

if ( ($Z1 == $Z2) || ($Z1*$Z2 == 0) ) {

    $DeltaE1= abs(($E[1]-$E[3])/$E[0]);
    $DeltaE2= abs(($E[0]-$E[4])/$E[0]);

    $DeltaQ11= abs((abs($Q_1t[1])-abs($Q_1t[3]))/$Q_1t[1]);
    $DeltaQ12= abs((abs($Q_1t[0])-abs($Q_1t[4]))/$Q_1t[0]);

    $DeltaQ21= abs((abs($Q_2t[1])-abs($Q_2t[3]))/$Q_2t[1]);
    $DeltaQ22= abs((abs($Q_2t[0])-abs($Q_2t[4]))/$Q_2t[0]);

    if ( $DeltaE1 > $DeltaE ) { $DeltaE=$DeltaE1;}
    if ( $DeltaE2 > $DeltaE ) { $DeltaE=$DeltaE2;}

    if ( $deltaE > $DeltaE1 ) { $DeltaE1=$deltaE;}
    if ( $deltaE > $DeltaE2 ) { $DeltaE2=$deltaE;}


    if ( $DeltaQ11 > $DeltaQ1 ) { $DeltaQ1=$DeltaQ11;}
    if ( $DeltaQ12 > $DeltaQ1 ) { $DeltaQ1=$DeltaQ12;}

    if ( $DeltaQ21 > $DeltaQ2 ) { $DeltaQ2=$DeltaQ21;}
    if ( $DeltaQ22 > $DeltaQ2 ) { $DeltaQ2=$DeltaQ22;}

    $debug2=1;
} else {
     $DeltaQ11= $DeltaQ1; 
     $DeltaQ12= $DeltaQ1;

     $DeltaQ21= $DeltaQ2; 
     $DeltaQ22= $DeltaQ2; 
}

$changeQ1=$Q_1t[1]*(1+$DeltaQ1);
$changeQ2=$Q_2t[1]*(1+$DeltaQ2);

$changeQ11=$Q_1t[1]*(1+$DeltaQ11);
$changeQ21=$Q_2t[1]*(1+$DeltaQ21);

$changeQ12=$Q_1t[1]*(1+$DeltaQ12);
$changeQ22=$Q_2t[1]*(1+$DeltaQ22);

$changeE=$E[1]*(1+$DeltaE);
$changeE1=$E[1]*(1+$DeltaE1);
$changeE2=$E[1]*(1+$DeltaE2);

$changeQ13=$Q_1t[3]*(1+$DeltaQ1);
$changeQ14=$Q_1t[4]*(1+$DeltaQ1);

print "        ","\n";
print '=' x 90, "\n";
print " title ", $header,"\n\n";

print " grid  [",$nu,"x",$mu,"/",$R_inf,"]\n\n";

print " Z1    ",$Z1,"\n";
print " Z2    ",$Z2,"\n";
print " R     ",$R," (au)    ",$R*0.529177249, " (A)\n";
print " z     ",$z,"\n\n";
print " F     ",$F,"\n";

print "\n";


if ( $homo eq  "yes" ) { print " homonuclear case forced \n\n"};

printf $form,' E(-2F) ',$E[0];
printf $form,' E(-1F) ',$E[1];
printf $form,' E( 0F) ',$E[2];
printf $form,' E(+1F) ',$E[3];
printf $form,' E(+2F) ',$E[4];

print "\n";


printf $form2,' Q1(-2F)',$Q_1t[0], ' Q2(-2F)',$Q_2t[0];
printf $form2,' Q1(-1F)',$Q_1t[1], ' Q2(-1F)',$Q_2t[1];
printf $form2,' Q1( 0F)',$Q_1t[2], ' Q2( 0F)',$Q_2t[2];
printf $form2,' Q1(+1F)',$Q_1t[3], ' Q2(+1F)',$Q_2t[3];
printf $form2,' Q1(+2F)',$Q_1t[4], ' Q2(+2F)',$Q_2t[4];
print "\n\n";

# printf $form2,' Q1(-2F)',$Q_1t[0], ' Q2(-2F)',$Q_2t[0], ' Q3(-2F)',$Q_3t[0],   ' Q4(-2F)',$Q_4t[0];
# printf $form2,' Q1(-1F)',$Q_1t[1], ' Q2(-1F)',$Q_2t[1], ' Q3(-1F)',$Q_3t[1],   ' Q4(-1F)',$Q_4t[1];
# printf $form2,' Q1( 0F)',$Q_1t[2], ' Q2( 0F)',$Q_2t[2], ' Q3( 0F)',$Q_3t[2],   ' Q4( 0F)',$Q_4t[2];
# printf $form2,' Q1(+1F)',$Q_1t[3], ' Q2(+1F)',$Q_2t[3], ' Q3(+1F)',$Q_3t[3],   ' Q4(+1F)',$Q_4t[3];
# printf $form2,' Q1(+2F)',$Q_1t[4], ' Q2(+2F)',$Q_2t[4], ' Q3(+2F)',$Q_3t[4],   ' Q4(+2F)',$Q_4t[4];


print " Relative total energy uncertainty due to orbital norms not being equal 1:\n";
printRelError ( "deltaE(-2F)=  ",$relativeE[0] );
printRelError ( "deltaE(-1F)=  ",$relativeE[1] );
printRelError ( "deltaE( 0F)=  ",$relativeE[2] );
printRelError ( "deltaE(+1F)=  ",$relativeE[3] );
printRelError ( "deltaE(+2F)=  ",$relativeE[4] );
printf "       relative=   %8.2e\n",$relativeE;
print "\n\n";

#printf "     absolute =  %8.2e\n",$absoluteE;



print " Relative uncertainty of dipole and quadrupole moments due to orbital norms\n";
print " not being equal 1:\n";
if ( $mmRelErrFromFiles eq "yes" ) {
    printRelError ( "deltaQ1(-2F)= ",$deltaQ1[0] );
    printRelError ( "deltaQ1(-1F)= ",$deltaQ1[1] );
    printRelError ( "deltaQ1( 0F)= ",$deltaQ1[2] );
    printRelError ( "deltaQ1(+1F)= ",$deltaQ1[3] );
    printRelError ( "deltaQ1(+2F)= ",$deltaQ1[4] );
    print "\n";
    printRelError ( "deltaQ2(-2F)= ",$deltaQ2[0] );
    printRelError ( "deltaQ2(-1F)= ",$deltaQ2[1] );
    printRelError ( "deltaQ2( 0F)= ",$deltaQ2[2] );
    printRelError ( "deltaQ2(+1F)= ",$deltaQ2[3] );
    printRelError ( "deltaQ2(+2F)= ",$deltaQ2[4] );

}

printf "\n\n";
printf " Default relative errors of total energy, dipole and quadrupole moments:\n";
printRelError ( "deltaE  = ",$deltaE);
printRelError ( "deltaQ1 = ",$deltaQ1);
printRelError ( "deltaQ2 = ",$deltaQ2);
printf "\n\n";



if ( $Z1 == $Z2 || $Z1*$Z2 == 0 ) {
    printf " Relative errors of total energy, dipole and quadrupole moments due to +/-F symmetry:\n";

    printRelError ( "DeltaE(1F) = ",abs(($E[3]-$E[1])/$E[2]) );
    printRelError ( "DeltaE(2F) = ",abs(($E[4]-$E[0])/$E[2]) );
    print "\n";

    printRelError ( "DeltaQ1(1F) = ",abs(($Q_1t[3]+$Q_1t[1])/$Q_1t[1]) );
    printRelError ( "DeltaQ1(2F) = ",abs(($Q_1t[4]+$Q_1t[0])/$Q_1t[0]) );

    print "\n";
    printRelError ( "DeltaQ2(1F) = ",abs(($Q_2t[3]-$Q_2t[1])/$Q_2t[1]) );
    printRelError ( "DeltaQ2(2F) = ",abs(($Q_2t[4]-$Q_2t[0])/$Q_2t[0]) );

    # $DeltaQ12= abs((abs($Q_1t[0])-abs($Q_1t[4]))/$Q_1t[0]);
    # $DeltaQ22= abs((abs($Q_2t[0])-abs($Q_2t[4]))/$Q_2t[0]);


    # $DeltaQ11= abs((abs($Q_1t[1])-abs($Q_1t[3]))/$Q_1t[1]);
    # $DeltaQ21= abs((abs($Q_2t[1])-abs($Q_2t[3]))/$Q_2t[1]);
    printf "\n\n";
}



printf " Relative errors of total energy, dipole and quadrupole moments:\n";
printRelError ( "DeltaE  = ",$DeltaE);
printRelError ( "DeltaQ1 = ",$DeltaQ1);
printRelError ( "DeltaQ2 = ",$DeltaQ2);
printf "\n\n";

printRelError ( "DeltaE(1F) = ",abs(($E[3]-$E[1])/$E[2]) );
printRelError ( "DeltaE(2F) = ",abs(($E[4]-$E[0])/$E[2]) );
print "\n";

printRelError ( "DeltaQ1(1F) = ",abs(($Q_1t[3]+$Q_1t[1])/$Q_1t[1]) );
printRelError ( "DeltaQ1(2F) = ",abs(($Q_1t[4]+$Q_1t[0])/$Q_1t[0]) );

print "\n";
printRelError ( "DeltaQ2(1F) = ",abs(($Q_2t[3]-$Q_2t[1])/$Q_2t[1]) );
printRelError ( "DeltaQ2(2F) = ",abs(($Q_2t[4]-$Q_2t[0])/$Q_2t[0]) );

print "\n";
print "\n";

# z-component of the dipole moment calculated from the 3-point formula
$mu3=-($E[3]-$E[1])/(2.0*$F)+(-abs($R/2.0+$z)*$Z1+abs($R/2.0-$z)*$Z2);

# z-component of the dipole moment calculated from the 5-point formula
$mu5=($E[4]-8*$E[3]+8*$E[1]-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);

# print $DeltaE1,"\n";
# print $E[1],"  ",$DeltaE1,"  ",$changeE1,"\n";

if ( $Oh2 eq "yes" ) {
    printf "%5s  %32.16f %7s %24.16f %6s \n","mu_z",$mu3,"(E h^2)",$Q_1t[2],"(Q_1)    ";
}
printf "%5s  %32.16f %7s %24.16f %6s \n","mu_z",$mu5,"(E h^4)",$Q_1t[2],"(Q_1)    ";

if ($debug2 != 0 ) {
# z-component of the dipole moment calculated from the 5-point formula
    $mu5=($E[4]-8*$E[3]+8*$changeE-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);
    printf "%9s  %28.16f %7s %24.16f %6s \n","   +Delta",$mu5,"(E h^4)",$Q_1t[2]*(1+$DeltaQ1),"(Q_1)    ";

#    $mu5=($E[4]-8*$E[3]+8*$changeE1-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);
#    printf "%9s  %28.16f %7s %24.16f %6s \n","   +Delta",$mu5,"(E h^4)",$Q_1t[2]*(1+$DeltaQ11),"(Q_1)    ";

#    $mu5=($E[4]-8*$E[3]+8*$changeE2-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);
#    printf "%9s  %28.16f %7s %24.16f %6s \n","   +Delta",$mu5,"(E h^4)",$Q_1t[2]*(1+$DeltaQ12),"(Q_1)    ";

}



# alpha_zz calculated from the energy and dipole moment values by
# means of the 3-point formula (-1 0 1)

$alpha3e=-($E[3]+$E[1]-2.0*$E[2])/($F*$F);
$alpha3q= ($Q_1t[3]-$Q_1t[1])/(2.0*$F);

# alpha_zz calculated from the energy and dipole moment values by
# means of the 5-point formula (1 -8  8 -1)

$alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$E[1]+$E[0])/(12*$F*$F);
$alpha5q= (-$Q_1t[4]+8*$Q_1t[3]-8*$Q_1t[1]+$Q_1t[0])/(12.0*$F);
$alpha5qp=$alpha5q;

print "          ","\n";
if ( $Oh2 eq "yes" ) {
    printf "%9s  %28.16f %7s %24.16f %6s \n","alpha_zz",$alpha3e,"(E h^2)",$alpha3q,"(Q_1 h^2)";
}
printf "%9s  %28.16f %7s %24.16f %6s \n","alpha_zz",$alpha5e,"(E h^4)",$alpha5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$changeE+$E[0])/(12*$F*$F);
    $alpha5q= (-$Q_1t[4]+8*$Q_1t[3]-8*$changeQ11+$Q_1t[0])/(12.0*$F);
    printf "%7s  %28.16f  %31.16f     \n","   +Delta",$alpha5e,$alpha5q;

    # $alpha5q= (-$Q_1t[4]+8*$Q_1t[3]-8*$changeQ1+$Q_1t[0])/(12.0*$F);
    # printf "%7s  %28.16f  %31.16f     \n","   +Delta",$alpha5e,$alpha5q;

    $alpha5q= (-$changeQ14+8*$Q_1t[3]-8*$Q_1t[1]+$Q_1t[0])/(12.0*$F);
    printf "%12s  %23.16f  %31.16f     \n","   +Delta(Q14)",$alpha5e,$alpha5q;

    # $alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$changeE1+$E[0])/(12*$F*$F);
    # $alpha5q= (-$Q_1t[4]+8*$changeQ13-8*$Q_1t[1]+$Q_1t[0])/(12.0*$F);
    # printf "%7s  %28.16f  %31.16f     \n","   +Delta",$alpha5e,$alpha5q;            

    # $alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$changeE2+$E[0])/(12*$F*$F);
    # $alpha5q= (-$changeQ14+8*$Q_1t[3]-8*$Q_1t[1]+$Q_1t[0])/(12.0*$F);
    # printf "%7s  %28.16f  %31.16f     \n","   +Delta",$alpha5e,$alpha5q;            
}

# zzz-component of the first hyperpolarizability calculated from the
# dipole moment values by means of the 3 and 5-point formulae (1 -2  1
# and -1 16 -30 16 -1)

$beta5e=(-$E[4]+2*$E[3]-2*$E[1]+$E[0])/(2*$F*$F*$F);
$beta3q=($Q_1t[3]+$Q_1t[1]-2.0*$Q_1t[2])/($F*$F);
$beta5q=(-$Q_1t[4]+16*$Q_1t[3]-30*$Q_1t[2]+16*$Q_1t[1]-$Q_1t[0])/(12*$F*$F);
$beta5qp=$beta5q;

print ("\n");
if ( $Oh2 eq "yes" ) {
    printf "%9s  %61.16f %6s \n","beta_zzz",$beta3q,"(Q_1 h^2)";
}
printf "%9s  %28.16f %7s %24.16f %6s \n","beta_zzz",$beta5e,"(E h^4)",$beta5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $beta5e=(-$E[4]+2*$E[3]-2*$changeE+$E[0])/(2*$F*$F*$F);
    $beta5q=-($Q_1t[4]-16*$Q_1t[3]+30*$Q_1t[2]-16*$changeQ11+$Q_1t[0])/(12*$F*$F);
    printf "%14s  %23.16f %7s %24.16f     \n","  +Delta(Q11)",$beta5e,"       ",$beta5q;            

    $beta5q=-($changeQ14-16*$changeQ13+30*$Q_1t[2]-16*$Q_1t[1]+$Q_1t[0])/(12*$F*$F);
    printf "%14s  %23.16f %7s %24.16f     \n","  +Delta(Q13)",$beta5e,"       ",$beta5q;            

}

if ($debug3 != 0) {
    print ("\n");
    print ($Q_1t[2],"  ",$Q_1[2], "\n");
    print ($Q_1[4], "   ", -$Q_1[4]+$Q_1[0],"\n");
    
    print ($Q_1[3], "   ",$Q_1[3]-$Q_1[1]," \n");
    
    print ((-$Q_1[4]+$Q_1[0]), "  ",2*($Q_1[3]-$Q_1[1])," \n");
    
    print (((-$Q_1[4]+$Q_1[0])+2*($Q_1[3]-$Q_1[1]))/2," \n");
}

# zzzz-component of the second hyperpolarizability calculated from the
# dipole moment values by means of the 5-point formulae

$gamma5q=-(-$Q_1t[4]+2*$Q_1t[3]-2*$Q_1t[1]+$Q_1t[0])/(2*$F*$F*$F);
$gamma5qp=$gamma5q;

printf "\n";
printf "%10s  %59.16f %6s \n"," gamma_zzzz",$gamma5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $gamma5q=-(-$changeQ14+2*$Q_1t[3]-2*$Q_1t[1]+$Q_1t[0])/(2*$F*$F*$F);
    printf "%14s %57.16f     \n", " +Delta(Q14)",$gamma5q;            

    $gamma5q=-(-$Q_1t[4]+2*$Q_1t[3]-2*$changeQ11+$Q_1t[0])/(2*$F*$F*$F);
    printf "%14s %57.16f     \n", " +Delta(Q11)",$gamma5q;            

    # $gamma5q=-(-$Q_1t[4]+2*$changeQ13-2*$Q_1t[1]+$Q_1t[0])/(2*$F*$F*$F);
    # printf "%9s %62.16f     \n", "   +Delta",$gamma5q;            
}


# Az,zz and Bzz,zz can be calculated as field derivatives of quadrupole moment 
$Azzz3q= ($Q_2t[3]-$Q_2t[1])/(2.0*$F);
$Azzz5q= (-$Q_2t[4]+8*$Q_2t[3]-8*$Q_2t[1]+$Q_2t[0])/(12.0*$F);
$Azzz5qp=$Azzz5q;

print("\n");
if ( $Oh2 eq "yes" ) {
    printf "%9s  %28.16f %7s %24.16f %6s \n"," A_z,zz  ",$Azzz3q ,"(  h^2)",$Azzz5q,"(Q_2 h^4)";
} else {
    printf "%-44s    %24.16f %6s \n"," A_z,zz  ",$Azzz5q,"(Q_2 h^4)";
}


if ($debug2 != 0 ) {
    $Azzz5q= (-$Q_2t[4]+8*$Q_2t[3]-8*$changeQ2+$Q_2t[0])/(12.0*$F);
    printf "%9s %62.16f     \n", "   +Delta",$Azzz5q;            
}

$Bzzzz3q=  ($Q_2t[3]+$Q_2t[1]-2.0*$Q_2t[2])/($F*$F);
$Bzzzz5q=  (-$Q_2t[4]+16*$Q_2t[3]-30*$Q_2t[2]+16*$Q_2t[1]-$Q_2t[0])/(12*$F*$F);
$Bzzzz5qp=$Bzzzz5q;

print("\n");
if ( $Oh2 eq "yes" ) {
    printf "%9s  %28.16f %7s %24.16f %6s \n"," B_zz,zz ",$Bzzzz3q,"(  h^2)",$Bzzzz5q,"(Q_2 h^4)";
} else {
    printf "%-42s    %26.16f %6s \n"," B_zz,zz ",$Bzzzz5q,"(Q_2 h^4)";
}

if ($debug2 != 0 ) {
    $Bzzzz5q=  (-$Q_2t[4]+16*$Q_2t[3]-30*$Q_2t[2]+16*$changeQ2-$Q_2t[0])/(12*$F*$F);
    printf "%9s %62.16f     \n", "   +Delta",$Bzzzz5q;
}


print '=' x 90, "\n";



if ( $debug4 != 0 ) {

    print " Total multipole moments:\n";
    for ($i=0; $i<=4; $i++) {
	print '-' x 70,"\n";
	print " File ",$i+1,"\n";
	printf " Q_1\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_1t[$i], $Q_1tda[$i];
	printf " Q_2\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_2t[$i], $Q_2tda[$i];
	printf " Q_3\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_3t[$i], $Q_3tda[$i];
	printf " Q_4\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_4t[$i], $Q_4tda[$i];
    }

}

if ( $debug5 != 0 ) {
    print '-' x 70,"\n";
    print "\n Electronic multipole moments:\n";
    for ($i=0; $i<=4; $i++) {
	print '-' x 70,"\n";
	print " File ",$i+1,"\n";
	printf " Q_1\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_1[$i], $Q_1da[$i];
	printf " Q_2\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_2[$i], $Q_2da[$i];
	printf " Q_3\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_3[$i], $Q_3da[$i];
	printf " Q_4\t %20.16e au  %25.16e Debye-Ang^0 \n",  $Q_4[$i], $Q_4da[$i];
    }

	print '-' x 70,"\n";
}

if ( $debug6 != 0 ) {
# LiH+

#printf " %-6.4f %20.11e %20.11e %20.11e %20.11e %20.11e ", $R, $E[2],$Q_1t[2],$Q_2t[2],$Q_3t[2],$Q_4t[2];
    printf " %-6.4f %20.11e %20.11e %20.11e ", $R, $E[2],$Q_1t[2],$Q_2t[2];
    print   "  [",$nu,"x",$mu,";",$R_inf,"]\n";
    


    if ( $Z1 == $Z2 || $Z1*$Z2 == 0 ) {
	# atoms and homonuclear diatomics
	printf " %-6.5f %10.6f &%14.9f &%17.2e &%16.8f &%20.2e &%17.10f \\\\ %4s%10.2e \n", $R, $F, $alpha5qp,$beta5qp,$gamma5qp,$Azzz5qp,$Bzzzz5qp,"%",abs($F*$F*$alpha5qp/$E[2]);
	#    printf " F*F*alpha_zz/E %17.4e \n", abs($F*$F*$alpha5qp/$E[2]);
    } else {
	# heteronuclear diatomics
	printf " %-6.5f %10.6f &%12.9f &%17.7f &%15.2f &%21.8f &%17.7f \\\\ %4s%10.2e \n", $R, $F, $alpha5qp,$beta5qp,$gamma5qp,$Azzz5qp,$Bzzzz5qp,"%",abs($F*$F*$alpha5qp/$E[2]);
	printf " %-6.5f %10.6f &%12.9f &%17.7f &%15.2f &%21.8f &%17.7f \\\\ %4s%10.2e \n", $R*0.529177249, $F, $alpha5qp,$beta5qp,$gamma5qp,$Azzz5qp,$Bzzzz5qp,"%",abs($F*$F*$alpha5qp/$E[2]);
	#    printf " F*F*alpha_zz/E %17.4e \n", abs($F*$F*$alpha5qp/$E[2]);
    }

    print '=' x 90, "\n";
}



