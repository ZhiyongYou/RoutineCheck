#!/bin/bash
source /cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/envb.sh

year=`date +%Y -d "-16 day"`
date=`date +%m%d -d "-16 day"`
date2=`date +%m%d -d "-15 day"`
baseyear=2020
basedate=0126
confyear=2019
confdate=1010

checkfilefir=/scratchfs/ybj/youzhiyong/RoutineCheck/statusCheck
srcdir=/workfs/ybj/youzhiyong/RoutineCheck/statusCheck
basefiledir=/scratchfs/ybj/youzhiyong/telBase/createThresh/basefile/$baseyear/$basedate
conf=/workfs/ybj/youzhiyong/RoutineCheck/statusCheck/config/hotpixel/${confyear}${confdate}
iptfiledir=$checkfilefir/iptfile
outcheckdir=$checkfilefir/CheckResulte
outrootdir=$checkfilefir/rootfile

cd $srcdir
tels="1 2 3 4 5 6 7 10"
for itel in $tels
do
	tel=`printf %02d $itel`
	ls /eos/lhaaso/decode/wfcta/$year/$date/ES.*WFCTA${tel}*.${year}${date}[12]*.status.root >$iptfiledir/${year}${date}.WFCTA$tel.txt
	ls /eos/lhaaso/decode/wfcta/$year/$date2/ES.*WFCTA${tel}*.${year}${date2}0*.status.root >>$iptfiledir/${year}${date}.WFCTA$tel.txt
	basefile=wfcta${tel}.txt
	echo $itel $tel
	./main ${iptfiledir}/${year}${date}.WFCTA$tel.txt ${outrootdir}/${year}${date}.WFCTA$tel.root $outcheckdir ${basefiledir}/${basefile} $conf/hotpixel.txt $itel >>$checkfilefir/ffffff
done
