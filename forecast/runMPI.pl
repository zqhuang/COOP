#!/usr/local/bin/perl
use Cwd;

#Use current directory as root
$curdir = cwd;


$ini = $ARGV[0];
$num = $ARGV[1];

open(Fout,">./scripts/script_MPI");
print Fout <<EMP;
#!/bin/csh -f
#PBS -N $ini
#PBS -l nodes=$num:ppn=8
#PBS -q workq
#PBS -r n
cd $curdir
mpirun -pernode ./DOCLIK $ini > ./scripts/$ini.log
EMP
close(Fout);

chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
chdir("../");
