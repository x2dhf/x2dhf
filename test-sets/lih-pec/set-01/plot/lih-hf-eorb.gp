set title 'Orbital energies'
set terminal png
set output "/home/jkob/github/x2dhf-v3/tests/lih-pec/set-01/plot/lih-hf-eorb.png"
set offsets 0.2,0.5,0.5,0.5
set auto fix
set xlabel "lih separation (au)"
set ylabel "HF orbital energies (absolute values in au)"
plot "/home/jkob/github/x2dhf-v3/tests/lih-pec/set-01/plot/lih-hf-eorb.data" u 1:2 w l t "1sigma" , "/home/jkob/github/x2dhf-v3/tests/lih-pec/set-01/plot/lih-hf-eorb.data" u 1:3 w l t "2sigma"
