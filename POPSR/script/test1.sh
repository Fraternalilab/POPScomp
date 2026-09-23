#! /bin/sh
#-------------------------------------------------------------------------------
# test1.sh : run the POPScomp interface pipeline over the test structure
# The 'pops' program must be on the PATH.
#
# The expected count pins two fixes: the robust z-score is scaled by 'mad' alone
# ('mad' already applies the 1.4826 factor, so the former '1.4862 * mad' scaled
# twice and selected 6 residues), and every chain pair is analysed, not only the
# first one.
#-------------------------------------------------------------------------------

expected=8
out=tests/1F3R_A-1F3R_B_Qinterface.dat

rm -f "$out"
Rscript popscomp_interface.R --pdb tests/1F3R.pdb --workdir tests || exit 1

test -s "$out" || { echo "FAIL: $out was not written"; exit 1; }
## data lines, without the header
n=`tail -n+2 "$out" | wc -l`
echo "interface residues: $n (expected $expected)"
test "$n" -eq "$expected" || { echo "FAIL: $n interface residues, expected $expected"; exit 1; }
