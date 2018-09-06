# SoftManipulatorDynamics
This is the code for running the soft manipulator dynamics in Matlab.
The main dynamics code is written in C and is built to be compiled as a .mex file for Matlab.
The library depends on gsl and that must be installed independently.

To compile the mex file run:
```matlab
%cable driven
mex dynamicsStable.c -lgsl -lgslcblas -lm

%tca driven
mex dynamicsStableTCA.c -lgsl -lgslcblas -lm
```
