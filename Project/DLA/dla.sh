mpicc dla.c -lm -o dla.exe
mpirun -H localhost:5 -np 5 dla.exe
