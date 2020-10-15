#!/bin/bash
Year=2020
Date=1010

srcdir=/eos/user/y/youzhiyong/RoutineCheck/eventCheck
iptlistdir=/eos/user/y/youzhiyong/RoutineCheck/eventCheck/iptlist/sample
outdir=/eos/user/y/youzhiyong/RoutineCheck/eventCheck/rootfile/sample
mkdir -p $outdir

#Tels="1 2 3 4 5 6"
Tels="1"
for Tel in ${Tels}
do
	Tel=`printf %02d $Tel`
	$srcdir/main $iptlistdir/WFCTA01.txt $outdir/$Year$Date.WFCTA${Tel}.check.root 2020 10 10
done
