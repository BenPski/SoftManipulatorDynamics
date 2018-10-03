# SoftManipulatorDynamics
This is the code for running the soft manipulator dynamics in Matlab.
The main dynamics code is written in C and is built to be compiled as a .mex file for Matlab.
The library depends on gsl and that must be installed independently.

On linux getting the mex compilation working is straightforward.
Most distributions should have a C/C++ compiler available, only thing would be an issue is making sure it is compatible with your version of Matlab.
gsl should be available in the package manager or a simple installation from source.
However, when I compile it with gcc I get linker errors, but it works fine with g++. This may be due to how the header file defines the structs?
Then the compilation is:
```matlab
% the dynamics are unified into one function for now
mex GCC=/usr/bin/g++ -output manip_dynamics matlab_interface.c dynamics.C -lgsl -lgslcblas -lm
```

For Windows this is quite a bit more involved and easier to mess up. I only know the way I did it works and it took quite a while of mucking around to get it to work.

To start we need to get MinGW-w64 running properly as that will provide the C/C++ compiler. 
For this I start with the installation of [MinGW](http://mingw.org/) and make sure to include MSYS in the distribution (note this is different from MinGW-w64). The result of the this should be the appearance of C:\MinGW.

Here it is a good idea to add MinGW\bin and MinGW\msys\bin to your Path.

Next we need to install MinGW-w64. I was unable to compile from source, but I did get TDM-GCC running just fine and it should be equivalent. Download the installer from [TDM-GCC](http://tdm-gcc.tdragon.net/) and run it. The result should be the appearance of C:\TDM-GCC. Now I am not sure if you could just use TDM and avoid the previous steps, but I think you would still need to install MSYS.

Again add TDM-GCC\bin to your Path.

Now gsl needs to be installed. Download the latest release of [gsl](https://www.gnu.org/software/gsl/). Extract the gsl archive to be in the msys folder, this should be C:\MinGw\msys and it should be "next to" the 1.0 folder. Now we will be able to install it the typical linux way.

Open a command window (powershell, cmd, etc.). The steps to run are:
*enter bash
*append the linux utilities to the beginning of the PATH (order matters), these will be TDM, MinGw, and msys in that order
*configure, make, and install

So, in the command window enter:
```bash
bash
PATH=/c/TDM-GCC-64/bin:/c/MinGW/bin:/c/MinGW/msys/1.0/bin:$PATH
./configure
make
make install
```
This should install gsl into msys\1.0\local. You should be able to change the last three lines to ./configure && make && make install if you'd like.

With all of that setup you should be able to compile things like C programs on your Windows system like you would in linux, though it may need some finagling. However, to compile mex files from Matlab we still need some setup.

In order to compile with mex your Matlab environement needs to know the location of the C/C++ compiler. For me it did not find this automatically, so to let Matlab know where to look run:
```matlab
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64') %this has to be ran every time Matlab starts up if you are using mex, may be worth it to set it in your startup file
```
So, now Matlab should find the compiler, you can check by running:
```matlab
mex -setup
```
and if there is no error, it should be working fine.

Finally to compile the C files to get the mex files, Matlab does not know where the libraries and headers are for gsl and does not automatically recognize ".a" files, so we have to help it out. To properly compile the files run:
```matlab
%compile dynamics
mex -output manip_dynamics -I'C:\MinGW\msys\1.0\local\include\' matlab_interface.c dynamics.c 'C:\MinGW\msys\1.0\local\lib\libgsl.a' 'C:\MinGW\msys\1.0\local\lib\libgslcblas.a'

%or can execute something similar to the make_manip_dynamics.m script
make_manip_dynamics
```

After that it should have compiled and you'll see "manip_dynamics.mexw64", which will allow you to run the Matlab code. Then to run the dynamics, both simulations are the same just with different flags.

```matlab
%cable = 0
[g,xi,eta] = manip_dynamics(0,q,eta,xi,dt);

%tca = 1
[g,xi,eta] = manip_dynamics(0,q,eta,xi,dt);
```


# Example Usage
To use the code it generally goes, grab initial states and then run the dynamics iteratively for every step.

There are some definite and simple to fix issues currently, but they can be avoided. These are some of the matrices are transposed from what they should be, there are pretty much no error checks or sensible numbers checks, lots of dead code, and not really an issue, but right now the TCA dynamics includes a neural net for the temperature computations, this only works in the range of 0-3V as the input.

```matlab
%cable case
%g is configuration, xi is strain, and eta is velocity
[g,xi,eta] = initDynamics(10); %10 discretization points
%run one step of dynamics
[g,xi,eta] = dynamicsCable([0;0;0],eta,xi,0.01); %[0,0,0] are the cable tensions

%tca case, with the temperature neural net
q = [1;0;0]; %the actuation voltage
[g,xi,eta,tcaTemps] = initDynamics(10);
[g,xi,eta,tcaTemps] = dynamicsTCA(q,eta,xi,0.01,tcaTemps);
```


# Notes:
All the numerical simulations are done in C. So, can use the dynamics from C or matlab, matlab just has an interface for calling the dynamics stepping function. Working on a python interface and better interfaces overall. 

The matlab calling does not check for errors, so it is pretty easy to make it crash.
