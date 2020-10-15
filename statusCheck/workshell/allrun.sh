#!/bin/bash
source /cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/envb.sh
iptdir=/workfs/ybj/youzhiyong/RoutineCheck/statusCheck/iptfile
for ifile in $iptdir/*.txt
do
	ifile=`basename $ifile`
	tel=${ifile:14:2}
	tel=`printf %2d $tel`
	echo $ifile $tel
	./main $iptdir/$ifile $ifile.root $tel >ffffff
done
