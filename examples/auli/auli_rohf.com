%Mem=64GB
%NProcShared=12
#P ROHF/UGBS1P units=au punch=mo gfinput integral(grid=450590)

Gold-lithium

0 1
Au
Li 1 R

R=4.303

