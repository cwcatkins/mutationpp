gfortran -c -J../../install/include/mutation++ -g energy_test.f
gfortran -o energy_test -J../../install/include/mutation++ -g energy_test.o ../../install/lib/libmutation++_fortran.so
