#!/bin/bash
day=$1
indir=/eos/user/y/youzhiyong/RoutineCheck/eventCheck/rootfile/2020/$day
tels="1 2 3 4 5 6"
for tel in $tels
do
	./draw $indir/2020${day}.WFCTA0${tel}.check.root.root
done
