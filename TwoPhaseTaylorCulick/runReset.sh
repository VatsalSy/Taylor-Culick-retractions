if [ -f "log" ];
then
rm log
fi

if [ -f "dump" ];
then
rm dump
fi

if [ -f "$file" ];
then
rm $file
fi

if [ -d "intermediate" ]
then
rm -r intermediate
fi

qcc -fopenmp -Wall -O2 $file.c -o $file -lm
export OMP_NUM_THREADS=4
./$file 1e-2 1e-6 1e-3
