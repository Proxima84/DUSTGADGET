# DUSTGADGET
Dust gas interaction is included in GADGET-2 code: https://wwwmpa.mpa-garching.mpg.de/gadget/
To work with the program code the package MPICH is need. The program code is compiled by a script: 
.\run.sh 
and run by the command: 
mpirun -np N ./main.exe cbd.param,
where N is number of processor cores for calculations.
The detailed documentation for the code is in the folder ./documentation
