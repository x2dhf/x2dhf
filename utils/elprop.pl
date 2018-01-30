#!/usr/bin/perl  
################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################
 

# This script facilitates the calculation of polarizabilities and hyperpolarizabilities
# from a series of calculations performed by means of the 2DHF program.

# Usage:

#   elprop.pl [--field=F] file_m2 file_m1 file_p0 file_p1 file_p2 [--deltae=f] [--deltaq=f]
#   elprop.pl [--field=F] --homo file_m2 file_m1 file_p0 [--deltae=f] [--deltaq=f]  
 
# where file_i are output files from the x2dhf program corresponding to
# calculations with the external finite electric field strengths -2F, -F, 0, +F
# and +2F, respectively. If possible, the F value is extracted from the second
# of these files. Otherwise it must be defined using the long option
# 'field'. For homonuclear systems calculations for the -2F, -F and 0 field
# values need only be performed. However, if available, results for the +F and
# +2F field strengths are used to assess the accuracy of total energies and
# multipole moments (see below).

# In order to assess the accuracy of calculated electric properties provide relative
# errors of total energies (deltae) and dipole moments (deltaq).

# The derivatives are evaluated by means of the 3 and 5-point central
# differences formulae in order to assess the accuracy of calculations.

# For further details see 

# J.~Kobus, D.~Moncrieff, and S.~Wilson, A comparison of the electric moments
# obtained from finite basis set and finite difference {Hartree-Fock}
# calculations for diatomic molecules, Phys. Rev. A 62 (2000), 062503/1--9.

# J.~Kobus, D.~Moncrieff, and S.~Wilson, Comparison of the polarizabilities and
# hyperpolarizabilities obtained from finite basis set and finite field, finite
# difference {Hartree-Fock} calculations for diatomic molecules, J. Phys. B:
# At. Mol. Opt. Phys. 34 (2001) 5127-5143.

#############################################################################################


use Getopt::Long;

GetOptions ("homo"     => \$homo,
            "field=f"  => \$field,
            "deltae=f" => \$deltaE,
            "deltaq=f" => \$deltaQ,
            "deltaq1=f" => \$deltaQ1,
            "deltaq2=f" => \$deltaQ2);


$ver=0;

if ( defined($homo) ) {
    $homo="yes";
    $narg=3;
} else {
    $homo="no";
    $narg=5;
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

$debug1=1; # set to nozero value for additional printouts
$debug2=0; # set to nozero value to monitor the dependence of alpha,
           # beta and gamma (h^4) values on the accuracy of the dipole moment

if ( defined($deltaE) || defined($deltaQ) ) {$debug2=1};

$dbleForm="%10s %25.15e \n";
$quadForm="%10s %41.30e \n";

$dbleForm1="%10s %25.15e %20s %25.15e \n";
$quadForm1="%10s %41.30e %20s %41.30e \n";

$form=$dbleForm;
$DeltaEdefault=1e-12;  
$DeltaQdefault=1e-12;  

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

open(X2DHF_OUT,"<$x2dhf_out") || die " cannot open: $x2dhf_out";
while (<X2DHF_OUT>){
    if ( /nuclei/i) {
# extract Z1 Z2 
#	($t[1],$Z1,$Z2)=split(" ",$_);
        ($t[1],$Z1,$Z2,$R)=split(" ",$_);
        if ( /angstroms/i) {
            $R=$R/0.529177249;
        }
    }

#     if ( $_ =~ /molecular system:/) {
# # extract the description of a case and the internuclear distance
# 	$_=<X2DHF_OUT>;
# 	$_=<X2DHF_OUT>;
# 	$header=$_;
# 	@t=split("=",$_);
# 	$t[2]=$t[1];
# 	@t=split(" ",$t[2]);
# 	$R=$t[0];
#     }

    if ( $_ =~ /grid:/) {
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
#	    if ($debug1) {print $ma,"  ",$mb, "  ",$z,"\n";}
    }

    if ( /electronic \(au\/Debye\-Ang\^k\)/ ) {	
	$ver=1; # output in new format
    }
}
close(X2DHF_OUT);

if ($ver == 0) {
    extract_data_0();
} else {
    extract_data_1();
}

sub extract_data_0 {

# examine the five output files (old format; versions < 1-2005) in
# turn

    for ($iarg=0; $iarg<$narg; $iarg++) {
	
	$x2dhf_out=$ARGV[$iarg]; 
	open(X2DHF_OUT,$x2dhf_out) || die " cannot open: $!";
	
	if ($debug1) {
	    printf "\n %6s %18s \n","File: ",$x2dhf_out;
	}
	while (<X2DHF_OUT>){
	    
	    if ( $_ =~ /centre\s+atomic\s+weight\s+z/ ) {
# extract atomic masses and the centre of mass
 		for ($i=1; $i<=7; $i++) {  # skip 7 lines 
 		    $_=<X2DHF_OUT>;
 		}


# extract dipole, quadrupole etc moments 
# electronic (au)           Q_1    Q_2  
# electronic (Debye-Ang**k) Q_1da  Q_2da 
# total (au)                Q_1t     Q_2t
# total (Debye-Ang**k)      Q_1tda   Q_2tda 
		
		if ($debug1) {
		    printf "%50s \n","multipole moments relative to centre of mass:";
		}
		
		@t = split (" ",$_);
		$Q_1t[$iarg]=$t[2]; $Q_2t[$iarg]=$t[3]; $Q_3t[$iarg]=$t[4]; $Q_4t[$iarg]=$t[5];
		if ($debug1) {
		    printf "%33s","electronic (au)         : ";
		    for ($i=2; $i<=5; $i++) {
			printf "%25.16e", $t[$i];
		    }
		    printf "\n";
	}
		
		$_=<X2DHF_OUT>;
		@t = split (" ",$_);
		$Q_1tda[$iarg]=$t[2]; $Q_2tda[$iarg]=$t[3]; $Q_3tda[$iarg]=$t[4]; $Q_4tda[$iarg]=$t[5];
		
		if ($debug1) {
		    printf "%33s","electronic (Debye-Ang^k): ";
		    for ($i=2; $i<=5; $i++) {
			printf "%25.16e", $t[$i];
		    }
		    printf "\n";
		}
		
		$_=<X2DHF_OUT>;
		@t = split (" ",$_);
		$Q_1t[$iarg]=$t[2]; $Q_2t[$iarg]=$t[3]; $Q_3t[$iarg]=$t[4]; $Q_4t[$iarg]=$t[5];
		
		if ($debug1) {
		    printf "%33s","total (au)              : ";
		    for ($i=2; $i<=5; $i++) {
			printf "%25.16e", $t[$i];
		    }
		    printf "\n";
		}
		
		$_=<X2DHF_OUT>;
		@t = split (" ",$_);
		$Q_1tda[$iarg]=$t[2]; $Q_2tda[$iarg]=$t[3]; $Q_3tda[$iarg]=$t[4]; $Q_4tda[$iarg]=$t[5];
		
		if ($debug1) {
		    printf "%33s","total (Debye-Ang^k)     : ";
		    for ($i=2; $i<=5; $i++) {
			printf "%25.16e", $t[$i];
		    }
		    printf "\n";
		}
	    }
	    
	    if ( $_ =~ / ... writing functions to disk .../) {
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
	}
	close(X2DHF_OUT);
    }
}


sub extract_data_1 {

# examine the five output files (new format, version >= 1-2005) in turn

    for ($iarg=0; $iarg<$narg; $iarg++) {
	
	$x2dhf_out=$ARGV[$iarg]; 
	open(X2DHF_OUT,$x2dhf_out) || die " cannot open: $x2dhf_out";
	
	if ($debug1) {
	    printf "\n %6s %18s \n","File: ",$x2dhf_out;
	}
	
	while (<X2DHF_OUT>){

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
	    
	    if ( $_ =~ / ... writing functions to disk .../) {
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
	}
	close(X2DHF_OUT);
    }
}

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

# For atoms and homonuclear molecules these values can be calculated
# from differences of the total energy and dipole moment values for +/-F.

# If $DelatE/$DeltaQ cannot be calculated the script tries to get these data from
# DeltaE.data/DeltaQ.data files. In the last resort the default values are used.

# Default values of $DelataE and $DeltaQ are accuracy dependent (see above)

 
if ( ($Z1 == $Z2) or ($Z1*$Z2 == 0) ) {
    $DeltaE1= abs(($E[1]-$E[3])/$E[0]);
    $DeltaE2= abs(($E[0]-$E[4])/$E[0]);
    $DeltaQ12= abs((abs($Q_1t[0])-abs($Q_1t[4]))/$Q_1t[0]);
    $DeltaQ22= abs((abs($Q_2t[0])-abs($Q_2t[4]))/$Q_2t[0]);


    $DeltaQ11= abs((abs($Q_1t[1])-abs($Q_1t[3]))/$Q_1t[1]);
    $DeltaQ21= abs((abs($Q_2t[1])-abs($Q_2t[3]))/$Q_2t[1]);

     if ( $DeltaE1 <= $DeltaE2 ) {
	$DeltaE=$DeltaE2; 
     } else {
 	$DeltaE=$DeltaE1;
     } 
 
     if ( $DeltaQ12 <= $DeltaQ11 ) {
 	$DeltaQ1=$DeltaQ11; 
     } else {
 	$DeltaQ1=$DeltaQ12;
     } 


     if ( $DeltaQ22 <= $DeltaQ21 ) {
 	$DeltaQ2=$DeltaQ21; 
     } else {
 	$DeltaQ2=$DeltaQ22;
     } 

    if ($debug1 == 1) {
	print "DeltaE1 = ",$DeltaE1,"\n";
	print "DeltaE2 = ",$DeltaE2,"\n";
	print "DeltaQ11= ",$DeltaQ11,"\n";
	print "DeltaQ12= ",$DeltaQ12,"\n";
	print "DeltaQ21= ",$DeltaQ21,"\n";
	print "DeltaQ22= ",$DeltaQ22,"\n";
    }
    $debug2=1;
} 

if ( defined($deltaE) ) {
    $DeltaE=$deltaE; 
}

if ( $DeltaE == 0 ) {$DeltaE= $DeltaEdefault;}


if ( defined($deltaQ) ) {
    $DeltaQ=$deltaQ; 
    $DeltaQ11=$deltaQ; 
    $DeltaQ12=$deltaQ; 
    $DeltaQ21=$deltaQ; 
}


if ( defined($deltaQ1) ) {
    $DeltaQ1=$deltaQ1; 
    $DeltaQ11=$deltaQ1; 
    $DeltaQ12=$deltaQ1; 
}

if ( defined($deltaQ2) ) {
    $DeltaQ21=$deltaQ2; 
}

if ( $DeltaQ1 == 0 ) {
    $DeltaQ1= $DeltaQdefault;
    $DeltaQ11=$DeltaQ1; 
    $DeltaQ12=$DeltaQ1; 
    $DeltaQ21=$DeltaQ1; 
}

$changeQ1=$Q_1[1]*(1+$DeltaQ1);
$changeQ2=$Q_2[1]*(1+$DeltaQ21);
$changeE=$E[1]*(1+$DeltaE);

print "          ","\n";
print '=' x 90, "\n";
print " molecular system:          ","\n";
print ($header," \n");
print " grid:    [",$nu,"x",$mu,";",$R_inf,"]\n\n";

print " Z1        ",$Z1,"\n";
print " Z2        ",$Z2,"\n";
print " R         ",$R,"\n";
print " z         ",$z,"\n\n";
print " F         ",$F,"\n";

print "\n";


if ( $homo eq  "yes" ) { print " homonuclear case forced \n\n"};

printf $form,' E(-2F) ',$E[0];
printf $form,' E(-1F) ',$E[1];
printf $form,' E( 0F) ',$E[2];
printf $form,' E(+1F) ',$E[3];
printf $form,' E(+2F) ',$E[4];

print "\n";

printf $form1,' Q1(-2F) ',$Q_1t[0],   ' Q2(-2F) ',$Q_2t[0];
printf $form1,' Q1(-1F) ',$Q_1t[1],   ' Q2(-1F) ',$Q_2t[1];
printf $form1,' Q1( 0F) ',$Q_1t[2],   ' Q2( 0F) ',$Q_2t[2];
printf $form1,' Q1(+1F) ',$Q_1t[3],   ' Q2(+1F) ',$Q_2t[3];
printf $form1,' Q1(+2F) ',$Q_1t[4],   ' Q2(+2F) ',$Q_2t[4];


printf "\n";
printf "Relative errors of total energy and dipole and quadrupole moments:\n";
printf "%7s %14.1e  \n",  " DeltaE",$DeltaE;
#printf "%7s %9.1e  \n",  " DeltaE*E(0)",abs($DeltaE*$E[2]);
printf "%7s %13.1e  \n",  " DeltaQ1",$DeltaQ1;
printf "%7s %13.1e  \n",  " DeltaQ2",$DeltaQ21;
printf "\n\n";

# z-component of the dipole moment calculated from the 3-point formula
$mu3=-($E[3]-$E[1])/(2.0*$F)+(-abs($R/2.0+$z)*$Z1+abs($R/2.0-$z)*$Z2);

# z-component of the dipole moment calculated from the 5-point formula
$mu5=($E[4]-8*$E[3]+8*$E[1]-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);

printf "%5s  %32.16f %7s %24.16f %6s \n","mu_z",$mu3,"(E h^2)",$Q_1t[2],"(Q_1)    ";
printf "%5s  %32.16f %7s %24.16f %6s \n","mu_z",$mu5,"(E h^4)",$Q_1t[2],"(Q_1)    ";


if ($debug2 != 0 ) {
# z-component of the dipole moment calculated from the 5-point formula

$mu5=($E[4]-8*$E[3]+8*$changeE-$E[0])/(12.0*$F)+(-($R/2.0+$z)*$Z1+($R/2.0-$z)*$Z2);
printf "%9s  %28.16f %7s %24.16f %6s \n","   +Delta",$mu5,"(E h^4)",$Q_1t[2]*(1+$DeltaQ1),"(Q_1)    ";
}


# alpha_zz calculated from the energy and dipole moment values by
# means of the 3-point formula (-1 0 1)

$alpha3e=-($E[3]+$E[1]-2.0*$E[2])/($F*$F);
$alpha3q= ($Q_1t[3]-$Q_1t[1])/(2.0*$F);

# alpha_zz calculated from the energy and dipole moment values by
# means of the 5-point formula (1 -8  8 -1)

$alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$E[1]+$E[0])/(12*$F*$F);
$alpha5q= (-$Q_1t[4]+8*$Q_1t[3]-8*$Q_1t[1]+$Q_1t[0])/(12.0*$F);

print "          ","\n";
printf "%9s  %28.16f %7s %24.16f %6s \n","alpha_zz",$alpha3e,"(E h^2)",$alpha3q,"(Q_1 h^2)";
printf "%9s  %28.16f %7s %24.16f %6s \n","alpha_zz",$alpha5e,"(E h^4)",$alpha5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $alpha5e= ($E[4]-16*$E[3]+30*$E[2]-16*$changeE+$E[0])/(12*$F*$F);
    $alpha5q= (-$Q_1[4]+8*$Q_1[3]-8*$changeQ1+$Q_1[0])/(12.0*$F);
    printf "%7s  %28.16f  %31.16f     \n","   +Delta",$alpha5e,$alpha5q;            
}

# zzz-component of the first hyperpolarizability calculated from the
# dipole moment values by means of the 3 and 5-point formulae (1 -2  1
# and -1 16 -30 16 -1)

$beta5e=(-$E[4]+2*$E[3]-2*$E[1]+$E[0])/(2*$F*$F*$F);
$beta3q=  ($Q_1t[3]+$Q_1t[1]-2.0*$Q_1t[2])/($F*$F);
$beta5q=  (-$Q_1t[4]+16*$Q_1t[3]-30*$Q_1t[2]+16*$Q_1t[1]-$Q_1t[0])/(12*$F*$F);

print ("\n");
printf "%9s  %61.16f %6s \n","beta_zzz",$beta3q,"(Q_1 h^2)";
printf "%9s  %28.16f %7s %24.16f %6s \n","beta_zzz",$beta5e,"(E h^4)",$beta5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $beta5e=(-$E[4]+2*$E[3]-2*$changeE+$E[0])/(2*$F*$F*$F);
    $beta5q=-($Q_1[4]-16*$Q_1[3]+30*$Q_1[2]-16*$changeQ1+$Q_1[0])/(12*$F*$F);
    printf "%9s  %28.16f %7s %24.16f     \n","  +Delta",$beta5e,"(E h^4)",$beta5q;            
}

$debug3=0;
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

$gamma5q=-(-$Q_1[4]+2*$Q_1[3]-2*$Q_1[1]+$Q_1[0])/(2*$F*$F*$F);


printf "\n";
printf "%10s  %59.16f %6s \n"," gamma_zzzz",$gamma5q,"(Q_1 h^4)";

if ($debug2 != 0 ) {
    $gamma5q=-(-$Q_1[4]+2*$Q_1[3]-2*$changeQ1+$Q_1[0])/(2*$F*$F*$F);
    printf "%9s %62.16f     \n", "   +Delta",$gamma5q;            
}


# Az,zz and Bzz,zz can be calculated as field derivatives of quadrupole moment 
$Azzz3q= ($Q_2[3]-$Q_2[1])/(2.0*$F);
$Azzz5q= (-$Q_2[4]+8*$Q_2[3]-8*$Q_2[1]+$Q_2[0])/(12.0*$F);

print("\n");
printf "%9s  %28.16f %7s %24.16f %6s \n","Az,zz   ",$Azzz3q ,"(  h^2)",$Azzz5q,"(Q_2 h^4)";

if ($debug2 != 0 ) {
    $Azzz5q= (-$Q_2[4]+8*$Q_2[3]-8*$changeQ2+$Q_2[0])/(12.0*$F);
    printf "%9s %62.16f     \n", "   +Delta",$Azzz5q;            
}

$Bzzzz3q=  ($Q_2[3]+$Q_2[1]-2.0*$Q_2[2])/($F*$F);
$Bzzzz5q=  (-$Q_2[4]+16*$Q_2[3]-30*$Q_2[2]+16*$Q_2[1]-$Q_2[0])/(12*$F*$F);

print("\n");
printf "%9s  %28.16f %7s %24.16f %6s \n","Bzz,zz  ",$Bzzzz3q,"(  h^2)",$Bzzzz5q,"(Q_2 h^4)";

if ($debug2 != 0 ) {
    $Bzzzz5q=  (-$Q_2[4]+16*$Q_2[3]-30*$Q_2[2]+16*$changeQ2-$Q_2[0])/(12*$F*$F);
    printf "%9s %62.16f     \n", "   +Delta",$Bzzzz5q;
}


print '=' x 90, "\n";

$debug=0;

if ( $debug != 0 ) {

    print " Total multipole moments:\n";
    for ($i=0; $i<=4; $i++) {
	print '-' x 70,"\n";
#	print " File ",$i+1,"\n";
	print " Q_1  ",  $Q_1t[$i], " au     ", $Q_1tda[$i],  " Debye-Ang^0 \n";
	print " Q_2  ",  $Q_2t[$i], " au     ", $Q_2tda[$i],  " Debye-Ang^1 \n";
	print " Q_3  ",  $Q_3t[$i], " au     ", $Q_3tda[$i],  " Debye-Ang^2 \n";
	print " Q_4  ",  $Q_4t[$i], " au     ", $Q_4tda[$i],  " Debye-Ang^3 \n";
    }

    print '-' x 70,"\n";
    print "\n Electronic multipole moments:\n";
    for ($i=0; $i<=4; $i++) {
	print '-' x 70,"\n";
	print " File ",$i+1,"\n";
	print " Q_1  ",  $Q_1t[$i], " au     ", $Q_1da[$i],  " Debye-Ang^0 \n";
	print " Q_2  ",  $Q_2t[$i], " au     ", $Q_2da[$i],  " Debye-Ang^1 \n";
	print " Q_3  ",  $Q_3t[$i], " au     ", $Q_3da[$i],  " Debye-Ang^2 \n";
	print " Q_4  ",  $Q_4t[$i], " au     ", $Q_4da[$i],  " Debye-Ang^3 \n";
    }
	print '-' x 70,"\n";
}
