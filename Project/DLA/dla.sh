echo "Phrase 1: Compile"
mpicc dla.c -lm -o dla.exe
echo "Phrase 2: Run Simulation"
mpirun -H localhost:5 -np 5 dla.exe
echo "Phrase 3: Visualize"
python3 visualize.py
echo "Done!!!"