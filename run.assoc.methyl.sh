#!/bin/bash -l
<<<<<<< HEAD
#
=======
>>>>>>> 6e1117c... Initial commit
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

echo "sh run.assoc.methyl.sh $*"

cgid=$1
chr=$2

rcmd="/restricted/projectnb/gaw20/user/jiayiwu/assoc.geno/assoc.methyl.R"

cmd="R --vanilla --slave --args $cgid $chr < $rcmd"

echo $cmd
eval $cmd  || exit

chmod 777 /restricted/projectnb/gaw20/user/jiayiwu/assoc.geno/_qlog_/chr${chr}/$cgid.$chr.qlog
chmod 777 /restricted/projectnb/gaw20/user/jiayiwu/assoc.geno/chr${chr}/gaw20.${chr}.*.$cgid.assoc.out.gz

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="


