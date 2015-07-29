Quick install.

1. Compile mex file:
For windows: 
%MATLAB\bin\mex trideigs.c %MATLAB\extern\lib\win32\lcc\libmwlapack.lib 

For linux: 
$MATLAB/bin/mex -llapack trideigs.c

You may don't need to compile if you using precompiled mex files:
For WindowsXP: trideigs.mexw32
For Linux x64: tridegs.mexa64 

2. Add the directory that contains compiled mex file to your path. 
3. Add tridiegs.m to the path also. 
4. In Matlab type: "help triegs" to read how to use this mex file.
5. See "Matlab version of qm1d" to see an example.
