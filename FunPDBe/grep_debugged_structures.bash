#! /bin/bash
## finds structures where error message suggests "'--coarse' switch"
cd /mnt/databases/JSON ; find ??/*.log -print0 | xargs --null -i bash -c "grep -l switch {}" > ~/POPScomp/FunPDBe/nullFiles_210106_debugged_coarse_210220.dat 

## finds structures where error message reports "Problem" with atom pair
cd /mnt/databases/JSON ; find ??/*.log -print0 | xargs --null -i bash -c "grep -l Problem {}" > ~/POPScomp/FunPDBe/nullFiles_210106_debugged_Problem_210220.dat 

## finds structures where error message reports "Problem" with too close atom pairs
cd /mnt/databases/JSON ; find ??/*.log -print0 | xargs --null -i bash -c "grep -l close {}" > ~/POPScomp/FunPDBe/nullFiles_210106_debugged_close_210220.dat 

