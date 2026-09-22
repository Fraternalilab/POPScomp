#! /bin/sh
#-------------------------------------------------------------------------------
# test6c : residues sharing a number across a chain boundary stay separate
# Regression: a new residue was started when the residue number or insertion code
# changed, but not when the chain changed, so the last residue of chain A and the
# first of chain B were merged into one residue (with Q(SASA) > 1).
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6c : residue numbers shared across chains                "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

residues() {
	test -s "$1" || { echo 0; return; }
	sed -n '/=== RESIDUE SASAs ===/,$p' "$1" | awk 'NF == 10 && $3 + 0 == $3' | wc -l
}

../src/pops --pdb 1f3r.pdb --outDirName "$d" --residueOut || exit 1
ref=`residues "$d"/pops.out`
test "$ref" -gt 0 || { echo "FAIL: no residues in $d/pops.out"; exit 1; }

# renumber chain B to start at the last residue number of chain A
last=`awk '/^ATOM/ && substr($0, 22, 1) == "A" { r = substr($0, 23, 4) + 0 } END { print r }' 1f3r.pdb`
first=`awk '/^ATOM/ && substr($0, 22, 1) == "B" { print substr($0, 23, 4) + 0; exit }' 1f3r.pdb`
off=`expr "$last" - "$first"`
awk -v off="$off" '/^ATOM/ { r = substr($0, 23, 4) + 0;
	if (substr($0, 22, 1) == "B") r += off;
	$0 = substr($0, 1, 22) sprintf("%4d", r) substr($0, 27) } { print }' 1f3r.pdb > "$d"/align.pdb

../src/pops --pdb "$d"/align.pdb --outDirName "$d" --residueOut || exit 1
n=`residues "$d"/pops.out`
echo "residues: $n (reference $ref)"
test "$n" -eq "$ref" || { echo "FAIL: $n residues, expected $ref"; exit 1; }
