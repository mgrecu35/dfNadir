gfortran -c -fPIC src/f90Types.Feb2015.f90
gfortran -c -fPIC src/f90DataTypes.Feb2015.f90
gfortran -c -fPIC src/random.f90
gfortran -c -fPIC src/nbin.f90
gfortran -c -fPIC src/cloudStub.f90
gfortran -c -fPIC -O3 src/b*f90
gfortran -c -fPIC -O3 src/gEnv.f90
gfortran -c -fPIC -O3 src/retTables2.f90
gfortran -c -fPIC src/readTables.Aug2015.nonsph.f90
gfortran -c -fPIC -O3 src/retTablesInt.half.f90
gfortran -c -fPIC -O3 src/allocateMem.Feb2015.f90
gcc -c -fPIC -O3 -I multiscatter-1.2.10/include src/multiscatter2_ascii.c

f2py -c -m pyHB2 src/multiscatter.f90  src/fhb1.py.f90 allocateMem.Feb2015.o beamConvP.o beamConvSet.o bisection2.o f90DataTypes.Feb2015.o f90Types.Feb2015.o  gEnv.o nbin.o random.o readTables.Aug2015.nonsph.o retTables2.o retTablesInt.half.o multiscatter2_ascii.o src/rosen.f src/gcloud.f src/emissivity.f band.o src/radtran.f multiscatter-1.2.10/lib/libmultiscatter.a

rm *.o
rm *.mod
