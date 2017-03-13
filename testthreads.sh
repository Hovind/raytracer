# !/bin/sh

bin="raytracer"
n="6"

for i in $(seq 1 8); do
	echo -en "threads\t$i"
	time mpirun -np $n $bin $i
	echo
done


