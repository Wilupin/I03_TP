#! /bin/bash

fail=0
options=""

p=1
for i in $*
do
    var=`expr match "$i" '\([a-zA-Z]*.\)='`
    val=`expr match "$i" '.*=\([0-9.]*\)'`
  
    case $var in
	n)
	    n=$val
	    options="$options n=$n"
	    ;;
	m)
	    m=$val
	    options="$options m=$m"
	    ;;
	proc)
	    p=$val
	    ;;
	thread)
	    t=$val
	    options="$options threads=$t"
	    ;;
	dt)
	    dt=$val
	    options="$options dt=$dt"
	    ;;
	it)
	    it=$val
	    options="$options it=$it"
	    ;;
	*)
	    echo "unrecognised option : $i"
	    echo "    (not in n=, m=, proc=, thread=, dt=, it=)"
	    fail=1
	    ;;
    esac
done

if [ $fail != "0" ]
then
    echo "exit on error"
    exit 0
fi

echo $options
cat <<EOF > job_poisson.qsub

#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -pe mpi $p

echo "Poisson : $n x $m"
echo "Got \$NSLOTS processors."

mpirun --map-by node -report-bindings ./PoissonMPI $options

# mpirun -map-by core -bind-to core -report-bindings ./PoissonMPI $options

EOF

qsub job_poisson.qsub
