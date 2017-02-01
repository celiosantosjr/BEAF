while :
do 
	echo "***************************************************'
        System Optimzer -- updating in time (each 60s)
---------------------------------------------------"
	rm -rvf .smp* > .trash
	DIA=`date +"%d/%m/%Y"`
	HORA=`date +"%H:%M"`
	echo "$DIA -- $HORA"
	pidof cp > .smp1
	pidof gzip > .smp2
	pidof gunzip > .smp3
	pidof sed > .smp4
	pidof usearch > .smp5
	pidof blastx > .smp6
	pidof fastqc > .smp99
	pidof cutadapt > .smp09
	pidof blastn > .smp86
	pidof python > .smp7
	pidof perl > .smp12
	pidof mummer > .smp8
	pidof postnuc > .smp9
	pidof gedit > .smp10
	cat .smp* > .smpa
	for line in `cat .smpa`
	do
		sudo renice -20 "$line"
	done
	echo "To scape hit CTRL+C"
	sleep 60
done
rm -rf .trash .smp*

