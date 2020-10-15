#!/bin/bash
Year=2020
Date=$1

srcdir=/workfs/ybj/youzhiyong/RoutineCheck/eventCheck
iptlistdir=/lhaasofs/user/youzhiyong/RoutineCheck/eventCheck/iptlist/$Year/$Date
outdir=/eos/user/y/youzhiyong/RoutineCheck/eventCheck/rootfile/$Year/$Date
jobdir=/lhaasofs/user/youzhiyong/RoutineCheck/eventCheck/job/$Year/$Date

mkdir -p $jobdir


echo '#!/bin/bash' >$jobdir/check_${Year}${Date}.sh
echo 'Year=2020' >>$jobdir/check_${Year}${Date}.sh
echo "Mon=${Date:0:2}" >>$jobdir/check_${Year}${Date}.sh
echo "Day=${Date:2:2}" >>$jobdir/check_${Year}${Date}.sh
echo 'Date=${Mon}${Day}' >>$jobdir/check_${Year}${Date}.sh
echo 'Tel=$1' >>$jobdir/check_${Year}${Date}.sh
echo 'echo hep_sub -g lhaaso $0 $Tel >&2' >>$jobdir/check_${Year}${Date}.sh
echo "" >>$jobdir/check_${Year}${Date}.sh
echo "srcdir=${srcdir}" >>$jobdir/check_${Year}${Date}.sh
echo "iptlistdir=${iptlistdir}" >>$jobdir/check_${Year}${Date}.sh
echo "outdir=${outdir}" >>$jobdir/check_${Year}${Date}.sh
echo "" >>$jobdir/check_${Year}${Date}.sh
echo 'mkdir -p $outdir' >>$jobdir/check_${Year}${Date}.sh
echo "" >>$jobdir/check_${Year}${Date}.sh
echo 'cd $srcdir' >>$jobdir/check_${Year}${Date}.sh
echo 'iTel=`printf %02d ${Tel}`' >>$jobdir/check_${Year}${Date}.sh
echo '$srcdir/main $iptlistdir/WFCTA${iTel}.txt $outdir/$Year$Date.WFCTA${iTel}.check.root ${Year} $Mon $Day' >>$jobdir/check_${Year}${Date}.sh
chmod 755 $jobdir/check_${Year}${Date}.sh
