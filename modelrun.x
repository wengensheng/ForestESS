#/bin/sh
gfortran src/esdvm.F90 src/main.F90 -o ess
./ess
gfortran src/BasalAreaAnalysis.F90 -o ana
./ana

rm ess
rm ana
rm esdvm.mod
