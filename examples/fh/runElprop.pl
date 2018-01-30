#!/usr/bin/perl -w

# In case the electric field strength value F is not included in an output file (prior to
# the September 2002 versions of the 2DHF program) the long option '--field=F' should be
# used to supply the value.

$fncore="fh-";
$ext="lst";

$fn_m2=$fncore . "m2" . "." . $ext;
$fn_m1=$fncore . "m1" . "." . $ext;
$fn_0=$fncore  . "0"  . "." . $ext;
$fn_p1=$fncore . "p1" . "." . $ext;
$fn_p2=$fncore . "p2" . "." . $ext;

system ("elprop.pl  $fn_m2 $fn_m1 $fn_0 $fn_p1 $fn_p2 \n" );

# For homonuclear molecules three output files suffice
  
#system ("elprop.pl --homo $fn_m2 $fn_m1 $fn_0 \n" ); 
#system ("elprop.pl --homo --field=0.0001 $fn_m2 $fn_m1 $fn_0 \n" ); 

exit;

