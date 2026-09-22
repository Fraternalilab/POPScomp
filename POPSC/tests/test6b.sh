#! /bin/sh
#-------------------------------------------------------------------------------
# test6b : residue numbering does not change the SASA
# Regression: bonds were searched between residue numbers n and n+1, so a gap in
# the numbering (or descending numbering) lost the peptide bonds and the SASA
# dropped by ~15 %. Bonds are now found by residue position within a chain.
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6b : numbering gaps and residue order                    "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

total() { # total SASA = last line of three numbers
	test -s "$1" || return
	awk 'NF == 3 && $3 + 0 == $3 { t = $3 } END { if (t != "") printf "%.2f", t }' "$1"
}
# the extracted value must exist: an empty result would compare equal to another
# empty result and let the test pass without measuring anything
require() { # require <value> <description>
	test -n "$1" || { echo "FAIL: no $2"; exit 1; }
}


# doubled residue numbers: a gap after every residue
awk '/^ATOM/ { r = substr($0, 23, 4) * 2;
	$0 = substr($0, 1, 22) sprintf("%4d", r) substr($0, 27) } { print }' 1f3r.pdb > "$d"/gap.pdb
# chain B written before chain A
awk '/^ATOM/ && substr($0, 22, 1) == "B"' 1f3r.pdb > "$d"/swap.pdb
awk '/^ATOM/ && substr($0, 22, 1) == "A"' 1f3r.pdb >> "$d"/swap.pdb

../src/pops --pdb 1f3r.pdb --outDirName "$d" || exit 1
ref=`total "$d"/pops.out`
require "$ref" "reference SASA in $d/pops.out"

for v in gap swap; do
	../src/pops --pdb "$d"/$v.pdb --outDirName "$d" || exit 1
	t=`total "$d"/pops.out`
	require "$t" "SASA of $v in $d/pops.out"
	echo "$v: $t (reference $ref)"
	test "$t" = "$ref" || { echo "FAIL: $v SASA $t, expected $ref"; exit 1; }
done
