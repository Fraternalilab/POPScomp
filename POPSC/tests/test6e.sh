#! /bin/sh
#-------------------------------------------------------------------------------
# test6e : command line handling
# Regressions: an option that takes a file name swallowed the following option
# ('--sigmaOut --atomOut' wrote the SFE output to a file called '--atomOut');
# unknown options and missing arguments exited with status 0.
#-------------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " test6e : command line handling                               "
echo "--------------------------------------------------------------"

d=`mktemp -d` || exit 1
trap 'rm -rf "$d"' 0
p=`pwd`/../src/pops
t=`pwd`
cd "$d" || exit 1

# an option in place of a file name is refused, and no such file is written
"$p" --pdb "$t"/5lff.pdb --sigmaOut --atomOut > /dev/null 2>&1
test $? -ne 0 || { echo "FAIL: '--sigmaOut --atomOut' was accepted"; exit 1; }
test ! -f "./--atomOut" || { echo "FAIL: a file named '--atomOut' was written"; exit 1; }

# unknown option
"$p" --pdb "$t"/5lff.pdb --nosuchoption > /dev/null 2>&1
test $? -ne 0 || { echo "FAIL: unknown option exited with status 0"; exit 1; }

# a probe radius that is not a number
"$p" --pdb "$t"/5lff.pdb --rProbe abc > /dev/null 2>&1
test $? -ne 0 || { echo "FAIL: '--rProbe abc' was accepted"; exit 1; }

# --help succeeds
"$p" --help > /dev/null 2>&1 || { echo "FAIL: --help exited non-zero"; exit 1; }

# a file name is still accepted
"$p" --pdb "$t"/5lff.pdb --sigmaOut sfe.out > /dev/null 2>&1 || exit 1
test -s sfe.out || { echo "FAIL: --sigmaOut wrote no file"; exit 1; }
echo "command line handling: ok"
