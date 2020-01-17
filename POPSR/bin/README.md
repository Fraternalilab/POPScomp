# Link the POPS executable
The R Shiny server expects the POPS executable to be
available under POPScomp/POPSR/bin/.

* Add a link to the POPS executable, for example:
```
cd POPScomp/POPSR/bin
ln -fs ../POPSC/src/pops
```

* Alternatively, if root or sudo permission is obtained,
add a link to a standard binary path:
```
cd /usr/local/bin
ln -s <path>/POPScomp/POPSC/src/bin
```
The latter adds the POPS executable to the computer's PATH variable
and 'pops' can be run from any path.

