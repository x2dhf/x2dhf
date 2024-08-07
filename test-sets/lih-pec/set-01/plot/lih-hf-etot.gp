set title 'Total energy'
set terminal png
set output "/home/jkob/github/x2dhf-v3/tests/lih-pec/set-01/plot/lih-hf-etot.png"
set offsets 0.2,0.5,0.5,0.5
set auto fix
set xlabel "lih separation (au)"
set ylabel "HF total energy (au)"
plot "/home/jkob/github/x2dhf-v3/tests/lih-pec/set-01/plot/lih-hf-etot.data" u 1:($2) w l t ""
