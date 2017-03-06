while :
do 
	echo "***************************************************'
        System Optimzer -- updating in time (each 60s)
---------------------------------------------------"
	DIA=`date +"%d/%m/%Y"`
	HORA=`date +"%H:%M"`
	echo "Date: $DIA -- Time: $HORA"
	pidof cp > copy.id
	pidof gzip > gzip.id
	pidof gunzip > gunzip.id
	pidof sed > sed.id
	pidof usearch > usearch.id
	pidof blastx > blastx.id
	pidof blastn > blastn.id
	pidof blastp > blastp.id
	pidof makeblastdb > makeblastdb.id
	pidof cutadapt > cutadapt.id
	pidof fastqc > fastqc.id
	pidof python > python.id
	pidof perl > perl.id
	pidof mummer > mummer.id
	pidof postnuc > postnuc.id
	pidof gedit > gedit.id
	cat *.id > OptimalList.txt
	rm -rf *.id
	for ID in `cat OptimalList.txt`
	do
		sudo renice -20 "$ID"
	done
	rm -rf OptimalList
	echo "To scape hit CTRL+C"
	sleep 60
done
rm -rvf *.id; rm -rvf OptimalList.txt

