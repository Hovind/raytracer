# !/bin/sh

bin="raytracer"

for i in $(seq 2 10); do
	echo -en "hosts\t$i"
	time mpirun -np $i $bin
	echo
done


