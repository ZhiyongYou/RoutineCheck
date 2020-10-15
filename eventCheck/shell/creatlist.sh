#!/bin/bash
###############################################################
##creat check list
##usage: ./creatlist.sh Date (./creatlist.sh 1010)
###############################################################

Year=2020
Date=$1
Mon=${Date:0:2}
Day=${Date:2:2}
utc1=`date -d "$Year-$Mon-$Day 00:00:00" +%s`
utc2=`expr $utc1 + 86400`
Year2=`date -d @${utc2}  "+%Y"`
Date2=`date -d @${utc2}  "+%m%d"`
echo $Year/$Date---$Year2/$Date2

decodedir=/eos/lhaaso/decode/wfcta/$Year/$Date
decodedir2=/eos/lhaaso/decode/wfcta/$Year2/$Date2
srcdir=/workfs/ybj/youzhiyong/RoutineCheck/eventCheck
iptlisedir=/lhaasofs/user/youzhiyong/RoutineCheck/eventCheck/iptlist/$Year/$Date

mkdir -p $iptlisedir

#Tels="1 2 3 4 5 6"
#Tels="1 2 3 5 6 10"
Tels="1 2 3 5 6 10 11"
for Tel in ${Tels}
do
	Tel=`printf %02d $Tel`
	ls $decodedir/ES*WFCTA${Tel}*${Year}${Date}1[89]*.event.root >$iptlisedir/WFCTA${Tel}.txt
	ls $decodedir/ES*WFCTA${Tel}*${Year}${Date}2*.event.root >>$iptlisedir/WFCTA${Tel}.txt
	ls $decodedir2/ES*WFCTA${Tel}*${Year2}${Date2}0[0123456]*.event.root >>$iptlisedir/WFCTA${Tel}.txt
done
