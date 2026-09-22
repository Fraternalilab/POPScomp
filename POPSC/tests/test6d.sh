#! /bin/sh
#-------------------------------------------------------------------------------
# test6d : the input format does not change the SASA
# The PDB, PDBML and mmCIF readers share one residue and chain definition, so they
# report the same atoms. Note that the supplied 1f3r.pdb is the pre-remediation
# entry, in which the OD1/OD2 and CD/CE labels differ from 1F3R.cif, so its total
# differs from the other two by a fraction of an Angstrom squared: only the atom
# count is compared for the PDB file.
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6d : PDB, PDBML and mmCIF agree                          "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

atoms() { # atom count reported by the reader
	awk '/processed atoms/ { print $NF }' "$1"
}

../src/pops --pdb 1f3r.pdb --outDirName "$d" > "$d"/log.pdb || exit 1
a=`test -s "$d"/pops.out && awk 'NF == 3 && $3 + 0 == $3 { t = $3 } END { printf "%.2f", t }' "$d"/pops.out`
../src/pops --pdbml 1f3r.xml --outDirName "$d" > "$d"/log.xml || exit 1
b=`test -s "$d"/pops.out && awk 'NF == 3 && $3 + 0 == $3 { t = $3 } END { printf "%.2f", t }' "$d"/pops.out`
../src/pops --mmcif 1F3R.cif --outDirName "$d" > "$d"/log.cif || exit 1
c=`test -s "$d"/pops.out && awk 'NF == 3 && $3 + 0 == $3 { t = $3 } END { printf "%.2f", t }' "$d"/pops.out`

for v in "$a" "$b" "$c"; do
	test -n "$v" || { echo "FAIL: a reader produced no total SASA"; exit 1; }
done
echo "PDB $a / PDBML $b / mmCIF $c"
echo "atoms: PDB `atoms "$d"/log.pdb` / PDBML `atoms "$d"/log.xml` / mmCIF `atoms "$d"/log.cif`"
test "$b" = "$c" || { echo "FAIL: PDBML $b and mmCIF $c differ"; exit 1; }
test "`atoms "$d"/log.pdb`" = "`atoms "$d"/log.cif`" || { echo "FAIL: PDB and mmCIF atom counts differ"; exit 1; }
test "`atoms "$d"/log.xml`" = "`atoms "$d"/log.cif`" || { echo "FAIL: PDBML and mmCIF atom counts differ"; exit 1; }
