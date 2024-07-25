#! /bin/bash
#check the number of jobs you can submit now and submit these number of jobs
ninqueue=$(squeue|grep iciq7258|wc|awk '{ print $1 }')
ntosub=$((366-$ninqueue))
dir=$(pwd)
i=0
for x in $dir/*/unsub
do 
	if [ $ntosub -gt $i ]
	then
		cd ${x%unsub}
		echo ${x%unsub}
		if [ -f unsub ]
		then
			sbatch job
			rm unsub
		fi
		i=$(($i+1))
	fi
done
echo "submitted $i jobs"
