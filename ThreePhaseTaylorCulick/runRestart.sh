qcc -fopenmp -Wall -O2 $file.c -o $file -lm
export OMP_NUM_THREADS=4
./$file
