#!/bin/bash
year=2020
date=0203
date2=`expr $date + 1`
date2=`printf %04d $date2`
tels="1 2 3 4 5 6 7 10"
for tel in $tels
do
#tel=$1
	tel=`printf %02d $tel`
	statusfiledir=/eos/lhaaso/decode/wfcta/$year/$date
	statusfiledir2=/eos/lhaaso/decode/wfcta/$year/$date2
	iptfiledir=/workfs/ybj/youzhiyong/RoutineCheck/statusCheck/iptfile
	ls $statusfiledir/*WFCTA${tel}*${year}${date}[12]*.status.root >${iptfiledir}/${year}${date}.WFCTA${tel}.txt
	ls $statusfiledir2/*WFCTA${tel}*${year}${date2}0*.status.root >>${iptfiledir}/${year}${date}.WFCTA${tel}.txt
done
