#!/bin/bash
datadir=../sacdata
stations=`ls -d ${datadir}/T*`
gausses=(1 2 3 5 8 10 12)
for gauss in ${gausses[@]};do
	#if [ $gauss != 3 ];then
	#	continue
	#fi
	for station in ${stations[@]};do
	        stationname=${station##*/}
		 #if [ $stationname != "T13" ];then
		 #	continue
		 #fi
		echo "current: "
		echo "    gauss=${gauss} stationname=${stationname}"  
	      ./rf_stack -S../sacdata/saclist_${stationname}_t9.txt -g30/90/5/4/0 -K3 -P2.5 -T0/-80/-10/60/-5 -G${gauss} -I../sacdata/${stationname} -O../rf_gcarcstack/${stationname}/gauss${gauss}	
	done
done
