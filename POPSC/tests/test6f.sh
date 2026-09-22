#! /bin/sh
#-------------------------------------------------------------------------------
# test6f : --silent writes the output files
# Regression: the output file was opened inside 'if (! silent)', so with --silent
# the first write went to a NULL file pointer and POPS crashed.
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6f : --silent                                            "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0

../src/pops --pdb 5lff.pdb --outDirName "$d" --silent --atomOut --residueOut --chainOut || exit 1

for f in pops.out popsb.out sigma.out; do
	test -s "$d"/$f || { echo "FAIL: $f is missing or empty"; exit 1; }
done
echo "files written: pops.out popsb.out sigma.out"
