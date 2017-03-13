# !/bin/sh

hostfile="hostfile.txt"
tmp="hosttmp.txt"

bin="raytracer"


for i in $(seq 2 10); do
	echo -en "hosts\t$i"
	head -n $i $hostfile > $tmp	
	time mpirun -np $i $bin --hostfile=$tmp
	echo
done


