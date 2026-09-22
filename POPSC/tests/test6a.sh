#! /bin/sh
#-------------------------------------------------------------------------------
# test6a : residues are counted
# Regression: the CA/N3 detection compared against space-padded atom names after
# the names had been stripped, so no residue was ever found ("residues = 0").
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6a : residue count                                       "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

../src/pops --pdb 5lff.pdb --outDirName "$d" --residueOut || exit 1

# 5lff is a 7-residue peptide
n=`sed -n '/=== RESIDUE SASAs ===/,$p' "$d"/pops.out | awk 'NF == 10 && $3 + 0 == $3' | wc -l`
echo "residues: $n"
test "$n" -eq 7 || { echo "FAIL: expected 7 residues, got $n"; exit 1; }
