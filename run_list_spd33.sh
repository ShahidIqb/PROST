#!/bin/bash

##....activate python 2 virtual environment for spd33...###
source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate py2

cd $2

#path to databases
NR90=$HOME/Downloads/uniref50db/uniref50
HHDB=$HOME/Downloads/uniprot20_2016_02/uniprot20_2016_02

#path to HHBLITS and psiblast

HHBLITS=$HOME/opt/anaconda3/envs/PROST_SEQvenv/bin/hhblits
psiblast=$HOME/opt/anaconda3/envs/PROST_SEQvenv/bin/psiblast
if [ ! -f $HHDB.cs219 ]; then HHDB=$HHDB_pub; fi
if [ ! -f $NR90.pal ]; then NR90=$NR90_pub; fi
#
#
PDIR=$(dirname $0)
xdir=$PDIR
xdir=".$xdir"
echo $xdir
ncpu=$OMP_NUM_THREADS
if [ "$ncpu" == "" ]; then ncpu=10; fi
if [ $# -lt 1 ]; then echo "usage: $0 *.seq"; exit 1; fi

for seq1 in $(shuf -e $1); do
	pro0=$(basename $seq1)
	pro1=$(basename $(basename $seq1 .seq) .pssm)
	[ -f $pro1.spd33 -o -f $pro0.spd33 ] && continue
	if [ ! -f $pro1.pssm -a ! -f $pro1.bla ]; then
		$psiblast -db $NR90 -num_iterations 3 -num_alignments 1 -num_threads $ncpu -query $seq1 -out  $pro1.bla -out_ascii_pssm ./$pro1.pssm #-out_pssm ./$pro1.chk
		[ ! -f $pro1.pssm ] && $xdir/script/seq2pssm.py $seq1 > $pro1.pssm  # using blosum when failed
		[ -f $pro1.pssm  ] && rm -f $pro1.bla
	fi
	if  [ ! -f $pro1.hhm ]; then
		$HHBLITS -i $seq1 -ohhm $pro1.hhm -d $HHDB -v0 -maxres 40000 -cpu $ncpu -Z 0 -o $pro1.hhr
		[ -f $pro1.hhm ] && rm -f $pro1.hhr
	fi
done
#
$xdir/script/spider3_pred.py $1