#! /bin/sh
#-------------------------------------------------------------------------------
# test6g : results do not depend on the order of the chains in the file
# Regression: the buried SASA of an atom was accumulated from its SASA at the
# moment each contact was processed, and the overlap count was incremented twice
# for one of the two atoms, so both changed when the chains were reordered.
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6g : chain order independence                            "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

# atom records sorted by chain, residue number and atom name, without the atom
# serial number, so that the two chain orders can be compared line by line
table() {
	test -s "$1" || return
	sed -n "/=== ATOM $2 ===/,\$p" "$1" |
		awk 'NF >= 11 && $1 + 0 == $1 { $1 = ""; print }' | sort
}

awk '/^ATOM/ && substr($0, 22, 1) == "B"' 1f3r.pdb > "$d"/swap.pdb
awk '/^ATOM/ && substr($0, 22, 1) == "A"' 1f3r.pdb >> "$d"/swap.pdb

../src/pops --pdb 1f3r.pdb --outDirName "$d" --atomOut || exit 1
table "$d"/pops.out SASAs > "$d"/a.sasa
table "$d"/popsb.out bSASAs > "$d"/a.bsasa

../src/pops --pdb "$d"/swap.pdb --outDirName "$d" --atomOut || exit 1
table "$d"/pops.out SASAs > "$d"/b.sasa
table "$d"/popsb.out bSASAs > "$d"/b.bsasa

test -s "$d"/a.sasa || { echo "FAIL: no atom SASAs in $d/pops.out"; exit 1; }
test -s "$d"/a.bsasa || { echo "FAIL: no atom bSASAs in $d/popsb.out"; exit 1; }
cmp -s "$d"/a.sasa "$d"/b.sasa || { echo "FAIL: atom SASAs depend on chain order"; exit 1; }
cmp -s "$d"/a.bsasa "$d"/b.bsasa || { echo "FAIL: atom bSASAs depend on chain order"; exit 1; }
echo "atom SASAs and bSASAs: identical for both chain orders"
