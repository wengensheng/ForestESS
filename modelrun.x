#/bin/sh
gfortran src/esdvm.f90 src/main.F90 -o ess
./ess
gfortran src/BasalAreaAnalysis.f90 -o ana
./ana
