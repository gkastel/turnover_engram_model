

function la_run {
	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	echo ./lamodel $LAPARAMS  
	./lamodel $LAPARAMS  
	#sleep 1
}


altparams=""

CmdParams="-P 4 -p 10 -R "
for toh in 0 ; do
	Suffix="REP{$toh}"
	for ws in 1440;  do
		for BT in 0 5 10 15 20;  do
			for run in  {0..9}; do
				
				LAPARAMS=" -o turnoverHotspots=${toh} -o nBranchesTurnover=${BT} ${CmdParams}  -T $ws -S 191$run -s ${Suffix}G_${BT}_${run} -G"
				la_run

				LAPARAMS=" -o turnoverHotspots=${toh}  -o nBranchesTurnover=${BT} ${CmdParams}  -T $ws -S 191$run -s ${Suffix}L_${BT}_${run} -L"
				la_run

			done
		done
	done
done



for npat in 10 ; do
	CmdParams="-P ${npat} -p 10 "
	for toh in 0 ; do
		Suffix="10CAP${npat}_{$toh}"
		for ws in 1440;  do
			for BT in 0 5 10 15 20;  do
				for run in  {0..9}; do
					
					LAPARAMS=" -o turnoverHotspots=${toh} -o nBranchesTurnover=${BT} ${CmdParams}  -T $ws -S 191$run -s ${Suffix}G_${BT}_${run} -G"
					la_run

					LAPARAMS=" -o turnoverHotspots=${toh}  -o nBranchesTurnover=${BT} ${CmdParams}  -T $ws -S 191$run -s ${Suffix}L_${BT}_${run} -L"
					la_run


				done
			done
		done
	done
done


