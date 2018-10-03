%for windows
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
mex -setup cpp
mex -output manip_dynamics -I'C:\MinGW\msys\1.0\local\include\' matlab_interface.c dynamics.c 'C:\MinGW\msys\1.0\local\lib\libgsl.a' 'C:\MinGW\msys\1.0\local\lib\libgslcblas.a'