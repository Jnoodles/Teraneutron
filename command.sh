#!/bin/sh
pwd_=$PWD
jobDir=/local/storage/$SLURM_JOB_ID
mkdir jobs
cp tetraNeutron.C coalescence.h ${jobDir}/
cd ${jobDir}/
g++ tetraNeutron.C coalescence.h -o run1 -fopenmp `root-config --cflags --glibs` 
./run1
mv tetraneutron.root ${pwd_}/jobs/
rm -rf /local/storage/${SLURM_JOB_ID}
