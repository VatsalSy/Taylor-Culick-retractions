qcc -fopenmp -Wall -O2 $file.c -o $file -lm
export OMP_NUM_THREADS=4
./$file 1e-2 1e-6 1e-3
