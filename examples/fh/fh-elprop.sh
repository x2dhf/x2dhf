#!/bin/bash

if [[ ! -e ./elprop.pl ]]; then
    ln -s ../../utils/elprop.pl .
fi

cat <<EOF | ../../utils/runelprop.pl --deltaq=1e-12 --deltae=1e-12
 fh-m2.lst
 fh-m1.lst
 fh-0.lst
 fh-p1.lst
 fh-p2.lst
EOF



