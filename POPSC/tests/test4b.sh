#! /bin/sh

echo "--------------------------------------------------------------"
echo " test4a.sh                                                    "
echo "--------------------------------------------------------------"

valgrind ../src/pops --mmcif 1F3R.cif.gz --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

