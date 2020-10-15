#!/bin/bash
year=2020
date=0203
date2=`expr $date + 1`
date2=`printf %04d $date2`
#tels="1 2 3 4 5 6 7 10"
tels="1"
for tel in $tels
do
#tel=$1
	tel=`printf %02d $tel`
	statusfiledir=/eos/lhaaso/decode/wfcta/$year/$date
	statusfiledir2=/eos/lhaaso/decode/wfcta/$year/$date2
	iptfiledir=/workfs/ybj/youzhiyong/RoutineCheck/statusCheck/iptfile/filebyfile
	for ifile in $statusfiledir/*WFCTA${tel}*${year}${date}[12]*.status.root
	do
		file=`basename $ifile`
		echo $ifile >${iptfiledir}/$file.txt
	done
	for ifile in $statusfiledir2/*WFCTA${tel}*${year}${date2}0*.status.root
	do
		file=`basename $ifile`
		echo $ifile >${iptfiledir}/$file.txt
	done

#ls $statusfiledir/*WFCTA${tel}*${year}${date}[12]*.status.root >${iptfiledir}/${year}${date}.WFCTA${tel}.txt
done
