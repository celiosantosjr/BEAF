#!/usr/bin/env bash

address=$(cd "$(dirname "")" && pwd)/$(basename "")
echo "Starting program at folder ${address}"
ConfigFile=config.file

PathToSpades="/media/coppini/SDATA2/Bioinfo/MAGs/Lib/spades/bin/"
PathToQuast="/media/coppini/SDATA2/Bioinfo/MAGs/Lib/quast/"

References_Folder="/media/coppini/SDATA2/Bioinfo/MAGs/Reference_seqs"

threads="4"

CFLR="U"
ver="BEAF (full)"
ans="yes" #If you're sure you want to keep databases files from one loop to the next, change this to 'Y'. If you're sure you're not reusing .udb, change to 'N'


	# ======================================================================================================================================================================================== #
	# =====================================================================================BEAF MODULES======================================================================================= #
	# ======================================================================================================================================================================================== #



# BEAF works in a modular fashion, so that the main pipeline will call different functions (modules) as it runs. Each step works in a separate function so that the program can continue from a specific function if ended abruptly (see "Continue function" in README.md)


make_kp () # This function reorders the config.file in order to keep buckets in case they could be used more than once. This will prevent the program from copying, trimming and assessing the quality of trimmage more than once for the same data, reusing files.
{
if [[ "$CFLR" == "Y" ]]; then echo "******Reusing previous settings and files"; else
	rm -rf $address/*.kp; rm -rf $address/config.tmp; rm -rf $address/config.file1
	cat $ConfigFile | awk NF > $address/config.file1
	echo "# Checking if any buckets must be stored..."
	echo "890_abc.123_XYZ" > $address/LastR1.kp
	sort -S 50% --parallel=${threads} -k3,3 $address/config.file1 > $address/doconfig.kp
	while read T1 T2 R1 R2 Ref SubRef Out; do
		LastR1=`cat $address/LastR1.kp`	
		if [[ "$R1" == "$LastR1" ]]
		then
			echo "Y" >> $address/Keep_config.kp
		else
			echo "N" >> $address/Keep_config.kp
		fi
		echo "$R1" > $address/LastR1.kp
	done < $address/doconfig.kp
	echo "N" >> $address/Keep_config.kp
	sed -i -e 1,1d $address/Keep_config.kp
	paste $address/doconfig.kp $address/Keep_config.kp > $address/config.tmp
	rm -rf $address/*.kp; rm -rf $address/config.file1
	mv $address/config.tmp $address/config.kp
	echo "1" > $address/CR.step; CFLR="N"
fi
}

Check () # This function checks each line in config.file, checking for possible format errors.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "1" ]]; then echo "******Jumping autocheking system"; else
	rm -rf $address/*.check
	echo "# Checking your config file..."
	sort -S 50% --parallel=${threads} -k7,7 $address/config.kp > $address/config.check
	touch $address/LastOut.check
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		case $T1 in
			G|g|P|p|N|n|16S|16s|16)
			;;
			*)
				echo "# You're using a wrong config.file format. On the first column (T1), where you're currently using '$T1', use only G (for genome analysis), P (for protein analysis) or N (for protein nucleotide sequences analysis)"
				rm -rf $address/config.kp $address/config.check
				exit
			;;
		esac
		case $T2 in
			R|r)
				case $R2 in
					*.gz)
						if ! [ -s $R2 ];
						then
							echo "# Check your R2 file in '$R2'. The program either couldn't find the file or the file is empty."
							rm -rf $address/config.kp $address/config.check
							exit
						fi
					;;
					*)
						echo "# You're using a wrong config format. On the fourth column (R2), where you're currently using '$R2', you must use a gzipped file instead."
						rm -rf $address/config.kp $address/config.check
						exit
					;;
				esac
			;;
			I|i|F|f)
				touch $address/config.kp
			;;
			*)
				echo "# You're using a wrong config.file On the second column (T2), where you're currently using '$T2', use only R (for paired end fastq files), I (for interleaved fastq file) or F (for interleaved fasta file)"
				rm -rf $address/config.kp $address/config.check
				exit
			;;
		esac
		case $R1 in
			*.gz)
				if ! [ -s $R1 ];
				then
					echo "# Check your R1 file in '$R1'. The program either couldn't find the file or the file is empty."
					rm -rf $address/config.kp $address/config.check
					exit
				fi
			;;
			*)
				echo "# You're using a wrong config.file On the third column (R1), where you're currently using '$R1', you must use a gzipped file instead."
				rm -rf $address/config.kp $address/config.check
				exit
			;;
		esac
		LastOut=`cat $address/LastOut.check`
		if [[ "${Out}" == "$LastOut" ]]
			then
				echo "# You're using a wrong config.file On your seventh column (Out), you've used the same name for your output folder more than once, repeating '${Out}'"
				rm -rf $address/config.kp $address/config.check
				exit

		fi
		echo "${Out}" > $address/LastOut.check
	done < $address/config.check
	rm -rf $address/*.check
	echo "2" > $address/CR.step; CFLR="N"
fi
}

TimeHeader () # Generates the header for the Log.tsv file (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	echo "_________________________________________________________________________________________________________
Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs|AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntgSize|ORFs|AvgSizeORF|TotalSizeORF|StdDevORF|MaxORFSize" >> $address/Log.tsv
	echo "3" > $address/CR.step; CFLR="N"
fi
}

SoftTimeHeader () # Generates the header for the Log.tsv file (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	echo "_________________________________________________________________________________________________________
Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs|AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntgSize" >> $address/Log.tsv
	echo "3" > $address/CR.step; CFLR="N"
fi
}

Trim () # Trims adapters from Illumina data and merges sequences R1 and R2 into one file, when using pair-end (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "3" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping reads trimming opperations"; else
	rm -rf $address/Buckets
	mkdir $address/Buckets
	date -u +%s > $address/Buckets/datestarttrimming.tmp
	case $T2 in 
		R|r)
			case $T1 in
				16S|16s|16)
					echo "# We're now trimming your files. Only R1 will be used, as interleaved files may result in wrong abundance values..."
					cutadapt --cores=${threads} --quiet --minimum-length 80 --max-n 0.01 --quality-base 24 --trim-n -a AGATCGGAAGAGC -e 0.1 -O 5 -m 15 -o $address/Buckets/FastaQ-zcat.gz $R1 # > $address/Buckets/d.tmp # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
				;;
				*)
					echo "# Trimming and merging your files..."
					cutadapt --cores=${threads} --quiet --interleaved --minimum-length 80 --max-n 0.01 --quality-base 24 --trim-n -a AGATCGGAAGAGC -A AGATCGGAAGAGC -e 0.1 -O 5 -m 15 -o $address/Buckets/FastaQ-zcat.gz $R1 $R2 # > $address/Buckets/d.tmp # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
				;;
			esac
		;;
		I|i)
			echo "# Now we are trimming your files..."
			cutadapt --quiet --minimum-length 80 --max-n 0.01 --quality-base 24 --trim-n -a AGATCGGAAGAGC -e 0.1 -O 5 -m 15 -o $address/Buckets/FastaQ-zcat.gz $R1 > $address/Buckets/d.tmp # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
		;;
		F|f)
			gunzip -c < $R1 > $address/Buckets/FastaQ-zcat.fa
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/Buckets/datestarttrimming.tmp)" | bc -l > $address/Buckets/trimmingtime.nmb
	rm -rf $address/Buckets/d.tmp
	echo "4" > $address/CR.step; CFLR="N"
fi
}

QAnConversion () # Makes assessment of trimmage and converts file from fastq to fasta (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "4" ]]; then echo "******Skipping FASTQ handling / FASTQ conversion"; else
	date -u +%s > $address/Buckets/datestartfastqc.tmp
	case $T2 in
		R|r|I|i)
			echo "# Starting quality assessment of trimming..."
			rm -rf $address/Buckets/FASTQCresults; rm -rf $address/Buckets/FastaQ-zcat.fa
			mkdir $address/Buckets/FASTQCresults
			fastqc --quiet --threads $threads -f fastq -o $address/Buckets/FASTQCresults $address/Buckets/FastaQ-zcat.gz # FASTQ assessment is done here
			rm -rf $address/Buckets/FASTQCresults/*_fastqc
			echo "# We will convert merged file to fasta format."
			gunzip -c $address/Buckets/FastaQ-zcat.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/Buckets/FastaQ-zcat.fa
			echo "# Removing unwanted files..."
			rm -rf $address/Buckets/*.gz
		;;
		F|f)
			rm -rf $address/Buckets/FASTQCresults
			echo "# FASTA file type identified. Since FASTA does not have PHRED values, It wont be assessed by FASTQC algorithm."
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/Buckets/datestartfastqc.tmp)" | bc -l > $address/Buckets/fastqctime.nmb
	echo "5" > $address/CR.step; CFLR="N"
fi
}

CopyFile () # Soft - Copies Illumina data to a separate folder to work on it (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "3" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping files copying process"; else
	rm -rf $address/Buckets
	mkdir $address/Buckets
	case $T2 in 
		R|r)
			echo "# Copying file 1 from storage..."
			cp -r $R1 $address/Buckets
			echo "# Copying file 2 from storage..."
			cp -r $R2 $address/Buckets
		;;
		I|i|F|f)
			echo "# Copying file from storage..."
			cp -r $R1 $address/Buckets
		;;
	esac
	echo "4" > $address/CR.step; CFLR="N"
fi
}

SoftMergeRename () # Soft - Merges R1 and R2 into one file, and converts it from fastq to fasta, when needed (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "4" ]]; then echo "******Skipping files merging / conversion process"; else
	case $T2 in 
		R|r)
			rm -rf $address/Buckets/FastaQ-zcat.gz $address/Buckets/FastaQ-zcat.fa
			echo "Merging files"
			zcat $address/Buckets/*.gz | gzip -c > $address/Buckets/FastaQ-zcat.gz
			echo "# We will convert merged file to fasta format."
			gunzip -c $address/Buckets/FastaQ-zcat.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/Buckets/FastaQ-zcat.fa
			rm -rf $address/Buckets/*.gz
		;;
		I|i)
			rm -rf $address/Buckets/FastaQ-zcat.gz $address/Buckets/FastaQ-zcat.fa
			gunzip -c `ls $address/Buckets/*.gz` | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/Buckets/FastaQ-zcat.fa
		;;
		F|f)
			rm -rf $address/Buckets/FastaQ-zcat.fa
			gunzip -c <`ls $address/Buckets/*.gz`> $address/Buckets/FastaQ-zcat.fa
			rm -rf $address/Buckets/*.gz
		;;
	esac
	echo "5" > $address/CR.step; CFLR="N"
fi
}

BucketEngine () # Breaks the data into separate buckets to speed up the process 
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "5" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping Buckets System"; else
	if [[ -s $address/Buckets/bk_list.txt ]]
	then
		reads=`cat $address/Buckets/reads.nmb`
		rm -rf $address/Buckets/*.bk
	else
		if [[ -s $address/Buckets/buckets_list.txt ]];
		then
			buckets=`ls $address/Buckets/*.bk | wc -l`
			reads=`cat $address/Buckets/reads.nmb`
			echo "$buckets buckets were found from the previous run. We'll be using those for this run as well."
			echo "(The $buckets buckets had already been generated in a previous run. This is the time it took in that run for them to be generated, not in the current run)" > $address/Buckets/bucketpreviouslygeneratedmessage.tmp
		else
			if [[ -s $address/Buckets/reads.nmb ]]
			then
				touch $address/Buckets/reads.nmb
			else
				grep ">" $address/Buckets/FastaQ-zcat.fa | wc -l > $address/Buckets/reads.nmb
				reads=`cat $address/Buckets/reads.nmb`
			fi
			original_size=$(wc -c $address/Buckets/FastaQ-zcat.fa | sed 's/ .*//') # ; echo "original size $original_size"
			# sizeinKb=`expr $original_size / 1024`
			# sizeinMb=`expr $sizeinKb / 1024` # `expr $original_size / 1048576`
			buckets=`expr $original_size / 268435456 + 1` # ; echo "buckets $buckets" # (256Mb buckets)
			bucketsize=`expr $original_size / $buckets + 1024` # ; echo "bucketsize $bucketsize" # (+1Kb)
			date -u +%s > $address/Buckets/datestartbucketengine.tmp
			if [ "$buckets" -eq "1" ];
			then
				if [[ $original_size -ge 1048576 ]]
				then
					echo "# We identified $reads reads (file size is $original_size bytes, or $(expr $original_size / 1048576)Mb). It would take no buckets, avoiding this step."
				else
					if [[ $original_size -ge 1024 ]]
					then
						echo "# We identified $reads reads (file size is $original_size bytes, or $(expr $original_size / 1024)Kb). It would take no buckets, avoiding this step."
					else
						echo "# We identified $reads reads (file size is only $original_size bytes). It would take no buckets, avoiding this step."
					fi
				fi
				mv $address/Buckets/FastaQ-zcat.fa $address/Buckets/1.bk
			else
				echo "# We identified $reads reads (file size is $original_size bytes, or $(expr $original_size / 1048576)Mb). It will take $buckets bucket steps (the size of each bucket will be $bucketsize bytes, or $(expr $bucketsize / 1048576)Mb)."
				echo "### Starting operation of cutting and readapting..."
				echo "## Generating buckets..."
				cat $address/Buckets/FastaQ-zcat.fa | parallel --pipe --block $bucketsize --recstart ">" "cat >$address/Buckets/{#}.bk"
				ls $address/Buckets/*.bk | sort -S 50% --parallel=${threads} -V > $address/Buckets/bk_testlist.txt
				try=0
				while [[ "$(cat bk_testlist.txt | wc -l)" -gt "$buckets" && "$(echo $(wc -c $(tail -n 1 $address/Buckets/bk_testlist.txt) | sed 's/ .*//') + $(wc -c $(tail -n 2 $address/Buckets/bk_testlist.txt | head -n 1) | sed 's/ .*//') | bc)" -lt "$(wc -c $address/Buckets/1.bk | sed 's/ .*//')" ]] ; do
					cat $(tail -n 1 $address/Buckets/bk_testlist.txt) >> $(tail -n 2 $address/Buckets/bk_testlist.txt | head -n 1)
					rm -rf $(tail -n 1 $address/Buckets/bk_testlist.txt)
					ls $address/Buckets/*.bk | sort -S 50% --parallel=${threads} -V > $address/Buckets/bk_testlist.txt
				done
				while [[ $(cat $address/Buckets/bk_testlist.txt | wc -l) -gt $buckets && ${try} -le 5 ]] ; do
					if [[ -s "$address/Buckets/$(echo $buckets + 1 | bc).bk" ]]
					then
						if [[ "$(echo $(wc -c $address/Buckets/$buckets.bk | sed 's/ .*//') + $(wc -c $address/Buckets/$(echo "$buckets + 1" | bc).bk | sed 's/ .*//') | bc)" -lt "$(wc -c $address/Buckets/1.bk | sed 's/ .*//')" ]]
						then
							echo "Fixing buckets (try $try)"
							cat $address/Buckets/$(echo "$buckets + 1" | bc).bk >> $address/Buckets/$buckets.bk
							rm -rf $address/Buckets/$(echo "$buckets + 1" | bc).bk
						fi
					fi
					ls $address/Buckets/*.bk > $address/Buckets/bk_testlist.txt
					if [[ "$(cat bk_testlist.txt | wc -l)" -gt "$buckets" ]]
					then
						rm -rf $address/Buckets/*.bk
						rm -rf $address/Buckets/bk_list.txt $address/Buckets/bk_testlist.txt
						original_size=$(wc -c $address/Buckets/FastaQ-zcat.fa | sed 's/ .*//') # ; echo "original size of file: $original_size"
						buckets=`expr $original_size / 268435456 + 1` # ; echo "buckets: $buckets (256Mb buckets)"
						bucketsize=`expr (($original_size / $buckets) + (1048576 \* $try))` # ; echo "bucketsize: $bucketsize (+1Mb than previous try)"
						cat $address/Buckets/FastaQ-zcat.fa | parallel --pipe --block $bucketsize --recstart ">" "cat >$address/Buckets/{#}.bk"
						ls $address/Buckets/*.bk > $address/Buckets/bk_testlist.txt
					fi
					ls $address/Buckets/*.bk > $address/Buckets/bk_testlist.txt
					((try++))
				done
				ls $address/Buckets/*.bk | sort -S 50% --parallel=${threads} -V > $address/Buckets/bk_testlist.txt
				while [[ "$(cat bk_testlist.txt | wc -l)" -gt "$buckets" && "$(echo $(wc -c $(tail -n 1 $address/Buckets/bk_testlist.txt) | sed 's/ .*//') + $(wc -c $(tail -n 2 $address/Buckets/bk_testlist.txt | head -n 1) | sed 's/ .*//') | bc)" -lt "$(wc -c $address/Buckets/1.bk | sed 's/ .*//')" ]] ; do
					cat $(tail -n 1 $address/Buckets/bk_testlist.txt) >> $(tail -n 2 $address/Buckets/bk_testlist.txt | head -n 1)
					rm -rf $(tail -n 1 $address/Buckets/bk_testlist.txt)
					ls $address/Buckets/*.bk | sort -S 50% --parallel=${threads} -V > $address/Buckets/bk_testlist.txt
				done
				rm -rf $address/Buckets/bk_list.txt $address/Buckets/bk_testlist.txt
				mv $address/Buckets/bk_testlist.txt $address/Buckets/bk_list.txt
			fi
			echo "$(date -u +%s) - $(cat $address/Buckets/datestartbucketengine.tmp)" | bc -l > $address/Buckets/bucketenginetime.nmb
		fi
	fi
	echo "### Removing temporary files (stage 1)..."
	rm -rf $address/Buckets/*.txt; rm -rf $address/Buckets/*.gz # ; rm -rf $address/Buckets/*.fa
	ls $address/Buckets/*.bk > $address/Buckets/buckets_list.txt
	buckets="$(cat $address/Buckets/buckets_list.txt | wc -l)"
	if [[ $buckets != 1 ]]
	then
		echo "A total of ${buckets} were generated."
	fi
	cat $address/Buckets/buckets_list.txt | sed "s@$address/Buckets/@@" | sort -S 50% --parallel=${threads} -V > $address/Buckets/buckets_search.txt
	rm -rf $address/Buckets/bk_list.txt $address/Buckets/bk_testlist.txt
	buckets=`ls $address/Buckets/*.bk | wc -l`
	echo "6" > $address/CR.step; CFLR="N"
fi
}

Filter1 () # Creates udb files from reference if needed, and then align each bucket against reference file to filter for homology.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "6" ]]; then echo "******Skipping Filtering System 1"; else
	rm -rf $address/Buckets/FastaQ-zcat.fa
	date -u +%s > $address/OUTPUT/${Out}/datestartfilter1.tmp
	case $T1 in
		G|g)
			echo "# Starting searches..."
			echo "### Making database from reference genome... "
			cp -r $References_Folder/$Ref $address/Buckets/none
			cat $address/Buckets/none | awk NF > $address/Buckets/none1
			sed -i '/>/d' $address/Buckets/none1
			cat $address/Buckets/none1 | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $address/Buckets/md8
			rm -rf $address/Buckets/none $address/Buckets/none1
			usearch -makeudb_usearch $address/Buckets/md8 -output $address/Buckets/$Ref.udb -threads ${threads} > $address/Buckets/makeudblog.txt # udb from reference is created here for genomes
			rm -rf $address/Buckets/makeudblog.txt
			rm -rf $address/Buckets/md8
			rm -rf $address/Buckets/*.m7
			case $Keep in
				Y|y)
					for buck in `cat $address/Buckets/buckets_search.txt`; do
						if [[ -s $address/Buckets/$buck.m8 ]]
						then
							touch $address/Buckets/$buck.m8
						else
							echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
							usearch -usearch_global $address/Buckets/$buck -db $address/Buckets/$Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm for genome binning should be specified here
							rm -rf $address/Buckets/usearch.tmp
							mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
							sed -i -e 1,1d $address/Buckets/buckets_search.txt
						fi
					done
				;;
				*)
					for buck in `cat buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
						usearch -usearch_global $address/Buckets/$buck -db $address/Buckets/$Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm for genome binning should be specified here
						rm -rf $address/Buckets/usearch.tmp
						mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
						rm -rf $address/Buckets/$buck
						touch $address/Buckets/$buck
						sed -i -e 1,1d $address/Buckets/buckets_search.txt
					done
					rm -rf $address/Buckets/*.bk
				;;
			esac
			rm -rf $address/Buckets/*.udb; rm -rf $address/Buckets/*.m7; rm -rf $address/Buckets/buckets_search.txt
		;;
		P|p|N|n)
			echo "# Starting searches..."
			if [[ -s $address/Buckets/udblist ]]
			then
				touch $address/Buckets/udblist
			else
				rm -rf $address/Buckets/*.udb
				case $Ref in
					*.fa|*.fasta|*.fas|*.faa|*.fna|*.fsa|*.FA|*.FASTA|*.FAS|*.FAA|*.FNA|*.FSA) # Tests if reference is in fasta format
						echo "# Recognized $Ref file as fasta format. Making udb..."
						cat $References_Folder/$Ref | awk NF > $address/Buckets/none
						sed -i '/>/d' $address/Buckets/none
						cat $address/Buckets/none | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $address/Buckets/md8
						usearch -makeudb_usearch $address/Buckets/md8 -output $address/Buckets/$Ref.udb -threads ${threads} > $address/Buckets/makeudblog.txt
						rm -rf $address/Buckets/makeudblog.txt
						rm -rf $address/Buckets/md8
						ls $address/Buckets/*.udb > $address/Buckets/udblist
					;;
					*.udb) # In case reference is not in faste format, tests if it is an udb file
						echo "# Recognized $Ref file as .udb format."
						cp $Ref $address/Buckets
						ls $address/Buckets/*.udb > $address/Buckets/udblist
					;;
					*) 
						echo "Couldn't recognize $Ref file as neither fasta nor udb format. Will try to use it as udb regardless."
						cp $Ref $address/Buckets/$Ref.trying.udb
						ls $address/Buckets/$Ref.trying.udb > $address/Buckets/udblist
					;;
				esac
			fi
			dbinuse=`cat $address/Buckets/udblist` # Either a file originally in udb or a fasta format converted to udb will be used here
			rm -rf $address/Buckets/*.m7
			case $Keep in
				Y|y)
					for buck in `cat buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
						usearch -usearch_local $address/Buckets/$buck -db $dbinuse -strand both -id 0.25 -evalue 1e-5 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm for protein/gene binning should be specified here
						rm -rf $address/Buckets/usearch.tmp
						mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
						sed -i -e 1,1d $address/Buckets/buckets_search.txt
					done
				;;
				*)
					for buck in `cat buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
						usearch -usearch_local $address/Buckets/$buck -db $dbinuse -strand both -id 0.25 -evalue 1e-5 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm for protein/gene binning should be specified here
						rm -rf $address/Buckets/usearch.tmp
						mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
						rm -rf $address/Buckets/$buck
						touch $address/Buckets/$buck
						sed -i -e 1,1d $address/Buckets/buckets_search.txt
					done
					rm -rf $address/Buckets/*.bk
				;;
			esac
			rm -rf $address/Buckets/*.udb; rm -rf $address/Buckets/udblist; rm -rf $address/Buckets/hits; rm -rf $address/Buckets/*.m7; rm -rf $address/Buckets/buckets_search.txt
		;;
		16S|16s|16) 
			echo "# Starting searches..."
			rm -rf $address/Buckets/*.m7
			case $Keep in
				Y|y)
					for buck in `cat $address/Buckets/buckets_search.txt`; do
						if [[ -s $address/Buckets/$buck.m8 ]]
						then
							touch $address/Buckets/$buck.m8
						else
							echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
							usearch -usearch_global $address/Buckets/$buck -db $References_Folder/$Ref -strand both -id 0.95 -evalue 1e-20 --maxhits 1 --maxaccepts 1 --maxrejects 100 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm should be specified here
							rm -rf $address/Buckets/usearch.tmp
							mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
							sed -i -e 1,1d $address/Buckets/buckets_search.txt
						fi
					done
				;;
				*)
					for buck in `cat $address/Buckets/buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref and keeping buckets... ${buck%.bk}/$buckets"
						usearch -usearch_global $address/Buckets/$buck -db $References_Folder/$Ref -strand both -id 0.95 -evalue 1e-20 --maxhits 1 --maxaccepts 1 --maxrejects 100 -matched $address/Buckets/$buck.m7 -threads ${threads} > $address/Buckets/usearch.tmp # Parameters of reads search by Usearch algorithm should be specified here
						rm -rf $address/Buckets/usearch.tmp
						mv $address/Buckets/$buck.m7 $address/Buckets/$buck.m8
						rm -rf $address/Buckets/$buck
						touch $address/Buckets/$buck
						sed -i -e 1,1d $address/Buckets/buckets_search.txt
					done
					rm -rf $address/Buckets/*.bk
				;;
			esac
			rm -rf $address/Buckets/*.udb; rm -rf $address/Buckets/*.m7; rm -rf $address/Buckets/buckets_search.txt
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartfilter1.tmp)" | bc -l > $address/OUTPUT/${Out}/filter1time.nmb
	rm -rf $address/OUTPUT/${Out}/datestartfilter1.tmp
	cat $address/Buckets/*.m8 > $address/Buckets/hits
	echo "7" > $address/CR.step; CFLR="N"
fi
}

PreLogGen () # Generates Log.txt in Output/${Out} folder, with the results from the first homology search with usearch.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "7" ]]; then echo "******Skipping generation of pre-log"; else
	echo "## Calculating statistics..."
	hits=`grep ">" $address/Buckets/hits | wc -l`
	if [[ -s $address/Buckets/reads.nmb ]]
	then
		reads=`cat $address/Buckets/reads.nmb`
	else
		if [[ -s $address/OUTPUT/${Out}/reads.nmb ]]
		then
			reads=`cat $address/OUTPUT/${Out}/reads.nmb`
		fi
	fi
	ppm1=`expr 1000000 \* $hits / $reads`
	touch $address/Buckets/bucketpreviouslygeneratedmessage.tmp
	echo -e "RESULTS:
	File from: $R1\t$R2
	Results: ${Out}
	Reference: $Ref
	Reads: $reads
	Buckets: $buckets
	Hits: $hits
	Portion in ppm: $ppm1
" > $address/Buckets/Log.txt
	if [[ -s $address/Buckets/bucketpreviouslygeneratedmessage.tmp ]]
	then
		echo "
Trimming, Quality Analysis and Buckets generation were performed in a previous run. The time each of these process took can be seen below:" >> $address/Buckets/Log.txt
	else
		echo "The following steps were taken in this run, prior to initial filtering, with the respective processing time for each of them:" >> $address/Buckets/Log.txt
	fi
	echo "	Time for Trimming: $(cat $address/Buckets/trimmingtime.nmb)
	Time for Quality Analysis: $(cat $address/Buckets/fastqctime.nmb)
	Time to generate buckets: $(cat $address/Buckets/bucketenginetime.nmb)s $(cat $address/Buckets/bucketpreviouslygeneratedmessage.tmp)


Time for the first filter: $(cat $address/OUTPUT/${Out}/filter1time.nmb)s" >> $address/Buckets/Log.txt
	echo "$ppm1" > $address/OUTPUT/${Out}/ppm1.nmb
	cp -r reads.nmb $address/OUTPUT/${Out}
	mv Log.txt $address/OUTPUT/${Out}
	if [[ -d $address/Buckets/FASTQCresults ]]
	then
		cp -r $address/Buckets/FASTQCresults $address/OUTPUT/${Out}
	fi
	cntg="0"
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	rm -rf $address/Buckets/*.txt; rm -rf $address/Buckets/*.m8
	echo "8" > $address/CR.step; CFLR="N"
fi
}

G_Prepare_SPADES () # For genomes, prepares files in SPADES folder in order to start assemblage.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "8" ]]; then echo "******Skipping Assembly Preparation Module"; else
	echo "# Making contigs for ${Out}..."
	cp $address/Buckets/hits $address/Buckets/hits.fasta
	mv $address/Buckets/hits.fasta $address/OUTPUT/${Out}
	echo "9" > $address/CR.step; CFLR="N"
fi
}

G_SPADES1 () # For genomes, starts SPADES process using high kmers (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" ]]; then echo "******Skipping Assembly Module with High Kmers"; else
	date -u +%s > $address/OUTPUT/${Out}/datestartspades.tmp
	python $address/Lib/spades/bin/spades.py --threads ${threads} -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $address/OUTPUT/${Out}/hits.fasta -o $address/OUTPUT/${Out}/assembly_${Out} > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for genomes should be specified here, using high Kmers
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades.tmp)" | bc -l > $address/OUTPUT/${Out}/spadestime.nmb
	rm -rf $address/OUTPUT/${Out}/logspades.txt $address/OUTPUT/${Out}/datestartspades.tmp
	echo "10" > $address/CR.step; CFLR="N"
fi
}

G_SPADES2 () # For genomes, in case the first SPADES process with high kmers didn't work, it retries the assembly process using lower kmers. Soft version skips the first step so it will try only the low kmer assembly.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" && `cat $address/CR.step` != "10" ]]; then echo "******Skipping Assembly Mode with Lower Kmers"; else
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta ]]
	then
		echo "SPADES ran properly with high kmers"
	else
		echo "# SPADES couldn't find contigs for ${Out} with high kmers. Trying again with lower kmers."
		rm -rf $address/OUTPUT/${Out}/assembly_${Out}
		date -u +%s > $address/OUTPUT/${Out}/datestartspades2.tmp
		python $address/Lib/spades/bin/spades.py --threads ${threads} -k 11,15,21,25,31,35,41,45,51 --only-assembler -s $address/OUTPUT/${Out}/hits.fasta -o $address/OUTPUT/${Out}/assembly_${Out} > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for genomes should be specified here, using lower Kmers
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades2.tmp)" | bc -l > $address/OUTPUT/${Out}/spades2time.nmb
		rm -rf $address/OUTPUT/${Out}/logspades.txt $address/OUTPUT/${Out}/datestartspades2.tmp
	fi
	echo "11" > $address/CR.step; CFLR="N"
fi
}

GA () # For genomes, analyses the assembly doing an assessment of it with Quast, and then finds ORFs from the contigs generated on SPADES (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping Genome Analysis"; else
	if [[ -s $address/OUTPUT/${Out}/assessment.tar.gz ]]
	then
		echo "******Skipping QUAST service"
	else	
		if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta ]]
		then
			date -u +%s > $address/OUTPUT/${Out}/datestartquast.tmp
			echo "# Analyzing draft putative genome...\n"
			cp -r $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta $address/OUTPUT/${Out}/spades_contigs.fasta
			tar -zcvf $address/OUTPUT/${Out}/SPADES_results.tar.gz $address/OUTPUT/${Out}/assembly_${Out} --remove-files
			python $address/Lib/quast/metaquast.py --threads ${threads} --silent -R $References_Folder/$Ref -o $address/OUTPUT/${Out}/assessment $address/OUTPUT/${Out}/spades_contigs.fasta # Assessment of assemble is done in this step
			echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartquast.tmp)" | bc -l > $address/OUTPUT/${Out}/quasttime.nmb
			rm -rf $address/OUTPUT/${Out}/quastlog.txt $address/OUTPUT/${Out}/datestartquast.tmp
			echo "\n### Compressing results..."
			tar -zcvf assessment.tar.gz assessment --remove-files
			echo "# QUAST service is finished for the file going to OUTPUT/${Out}"
		else
			echo "# The proposed analysis of ${Out} could not continue due to problems in SPADES assembly."
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/spades_contigs.fasta ]]
	then
		echo "# Initiating ORF finding process for file going to ${Out}"
		rm -rf $address/OUTPUT/${Out}/ORFs.${Out}.fna
		date -u +%s > $address/OUTPUT/${Out}/datestartorffinder.tmp
		perl $address/Lib/bb.orffinder.pl --infile=$address/OUTPUT/${Out}/spades_contigs.fasta --outfile=$address/OUTPUT/${Out}/ORFs.${Out}.fna --minlen=200 --fasta > $address/OUTPUT/${Out}/perlog.txt # Parameters for genome ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffinder.tmp)" | bc -l > $address/OUTPUT/${Out}/orffindertime.nmb
		rm -rf $address/OUTPUT/${Out}/perlog.txt $address/OUTPUT/${Out}/datestartorffinder.tmp
		cntg=`grep ">" $address/OUTPUT/${Out}/spades_contigs.fasta | wc -l`
		if [[ -s $address/OUTPUT/${Out}/ORFs.${Out}.fna ]]
		then
			ORFs=`grep ">" $address/OUTPUT/${Out}/ORFs.${Out}.fna | wc -l`
			echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
			gzip $address/OUTPUT/${Out}/ORFs.${Out}.fna
			mv $address/OUTPUT/${Out}/ORFs.${Out}.fna.gz $address/OUTPUT/${Out}/ORFs.fa.gz
		else
			rm -rf $address/OUTPUT/${Out}/ORFs.${Out}.fna
			ORFs="0"
			echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
		fi
		echo "A total of $ORFs ORFs were found for reference $Ref, from $cntg contigs"
		gzip $address/OUTPUT/${Out}/spades_contigs.fasta
	else
		echo "# The proposed analysis of ${Out} could not continue due to problems in SPADES assembly."
		cntg="0"
	fi
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	echo "12" > $address/CR.step; CFLR="N"
fi
}

SoftGA () # Soft - For genomes, counts the number of contigs found after assemblage and organizes files (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping calculation of contigs and output compression"; else
	gzip $address/OUTPUT/${Out}/hits.fasta
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta ]]
	then
		cntg=`grep ">" $address/OUTPUT/${Out}/spades_contigs.fasta | wc -l`
		gzip $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta
		mv $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta.gz $address/OUTPUT/${Out}
	else
		cntg="0"
	fi
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	ORFs="[Warning: ORFs are not calculated in Soft version]"
	echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
	echo "12" > $address/CR.step; CFLR="N"
fi
}

BlastDBGen () # For proteins and genes, generates blast dbs from subreference files.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "8" ]]; then echo "******Skipping Preparation and generation of Blast Databases"; else
	cp $address/Buckets/hits $References_Folder/$SubRef
	rm -rf $References_Folder/$SubRef/SubRef_fasta.list; rm -rf $References_Folder/$SubRef/BlastDBlist
	ls $References_Folder/$SubRef/*.fa >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FA >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.fasta >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FASTA >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.fas >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FAS >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.faa >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FAA >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.fna >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FNA >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.fsa >> $References_Folder/$SubRef/SubRef_fasta.list
	ls $References_Folder/$SubRef/*.FSA >> $References_Folder/$SubRef/SubRef_fasta.list
	if [[ -s $References_Folder/$SubRef/SubRef_fasta.list ]]
	then
		date -u +%s > $address/OUTPUT/${Out}/datestartblastdbgen.tmp
		echo "# Recognized $SubRef files in fasta format. Making blast databases..."
		for sub in `cat $References_Folder/$SubRef/SubRef_fasta.list`; do
			case $T1 in
				P|p)
					makeblastdb -in $sub -dbtype prot -out $sub.db
				;;
				N|n)
					makeblastdb -in $sub -dbtype nucl -out $sub.db
				;;
			esac
		done
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartblastdbgen.tmp)" | bc -l > $address/OUTPUT/${Out}/blastdbgentime.nmb
	fi
	case $T1 in
		P|p)
			ls $References_Folder/$SubRef/*.psq | sed 's/.psq//' | sort -S 50% --parallel=${threads} -k1,1 > $References_Folder/$SubRef/BlastDBlist
		;;
		N|n)
			ls $References_Folder/$SubRef/*.nsq | sed 's/.nsq//' | sort -S 50% --parallel=${threads} -k1,1 > $References_Folder/$SubRef/BlastDBlist
		;;
	esac
	echo "9" > $address/CR.step; CFLR="N"
fi
}

Filter2 () # For proteins and genes, makes the second homology search, using blast, searching against each specific subreference.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" ]]; then echo "******Skipping Second Filtering System"; else
	date -u +%s > $address/OUTPUT/${Out}/datestartfilter2.tmp
	for sub in `cat $References_Folder/$SubRef/BlastDBlist`; do
		echo "# Searching against ${sub#"$References_Folder/$SubRef/"}..."
		case $T1 in
			P|p)
				dsp1=`date -u "+%s"`
				blastx -db $sub -query $address/Buckets/hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads $threads -outfmt 6 # Parameters of reads search by blast for proteins should be specified here
				dsp2=`date -u "+%s"`
				echo $dsp2 - $dsp1 |bc -l > $References_Folder/$SubRef/$sub.ft.time2
			;;
			N|n)
				dsp1=`date -u "+%s"`
				blastn -db $sub -query $address/Buckets/hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads $threads -outfmt 6 # Parameters of reads search by blast for genes should be specified here
				dsp2=`date -u "+%s"`
				echo $dsp2 - $dsp1 |bc -l > $References_Folder/$SubRef/$sub.ft.time2
			;;
		esac
		cat $sub.tmp | sort -S 50% --parallel=${threads} -k3,3 -k4,4 -n -r | awk '$3 > 90 && $4 > 25' | uniq > $References_Folder/$SubRef/$sub.ft
		touch $sub.ft
		rm -rf $sub.tmp
		sed -i -e 1,1d $References_Folder/$SubRef/BlastDBlist
	done
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartfilter2.tmp)" | bc -l > $address/OUTPUT/${Out}/filter2time.nmb
	rm -rf $References_Folder/$SubRef/BlastDBlist; rm -rf $References_Folder/$SubRef/*.tmp ; rm -rf $address/OUTPUT/${Out}/datestartfilter2.tmp
	echo "10" > $address/CR.step; CFLR="N"
fi
}

SaveDBs () # For proteins and genes, keeps subreference DBs in case ans (see line 3) is set to 'Y', moving fasta files to a folder used just to store them. If ans is 'N', removes blast DBs in this step.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "10" ]]; then echo "******Skipping the process of saving/removing databases"; else
	case $ans in
		Y|y|YES|Yes|yes)
			if [[ -s $References_Folder/$SubRef/SubRef_fasta.list ]]
			then
				if [[ -d $References_Folder/$SubRef/Fasta_files ]]
				then
					touch $References_Folder/$SubRef/Fasta_files
				else
					mkdir $References_Folder/$SubRef/Fasta_files
				fi
				echo "# Databases of subreference $SubRef now saved to References_seqs folder in blastdb format. Fasta files used to make the databases have been realocated to $References_Folder/$SubRef/Fasta_files."
				for file in `cat $References_Folder/$SubRef/SubRef_fasta.list`; do
					mv $file $References_Folder/$SubRef/Fasta_files
				done
			fi
		;;
		*)
			echo "# Removing blastdb databases generated using fasta files in the subreference folder (Reference_seqs/$SubRef)"
			if [[ -s $References_Folder/$SubRef/SubRef_fasta.list ]]
			then
				for file in `cat $References_Folder/$SubRef/SubRef_fasta.list`; do
					case $T1 in
						P|p)
							rm -rf $file.db.p*
						;;
						N|n)
							rm -rf $file.db.n*
						;;
					esac
				done
			fi
		;;
	esac
	rm -rf $References_Folder/$SubRef/SubRef_fasta.list
	echo "11" > $address/CR.step; CFLR="N"
fi
}

Extraction () # For proteins and genes, moves files to OUTPUT/${Out} to continue the process, and creates a python file to extract sequences.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping Extraction Module"; else
	echo "# Arranging data..."
	mv $References_Folder/$SubRef/*.ft $address/OUTPUT/${Out}
	mv $References_Folder/$SubRef/*.time2 $address/OUTPUT/${Out}
	mv $address/Buckets/hits $address/OUTPUT/${Out}
	mkdir $address/OUTPUT/${Out}/read_hits $address/OUTPUT/${Out}/contigs $address/OUTPUT/${Out}/blast_hits $address/OUTPUT/${Out}/ORFs
	mkdir $address/OUTPUT/${Out}/contigs/SPADES
	echo "12" > $address/CR.step; CFLR="N"
fi
}

PN_Prepare_SPADES () # For proteins and genes, prepares files for SPADES assembly, using python to extract blast matches from the hits, removing redundancy.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" ]]; then echo "******Skipping Assembly Preparation"; else
	rm -rf $address/OUTPUT/${Out}/*.rev; rm -rf $address/OUTPUT/${Out}/*.hits; rm -rf $address/OUTPUT/${Out}/*.rev-hits
	cut -f 1 $address/OUTPUT/${Out}/$File > $address/OUTPUT/${Out}/$File.rev
	python $address/Lib/ext.py $address/OUTPUT/${Out}/$File.rev $address/OUTPUT/${Out}/hits > $address/OUTPUT/${Out}/extpylog.txt
	rm -rf $address/OUTPUT/${Out}/extpylog.txt
	touch $File.rev-hits
	mv $File.rev-hits $address/OUTPUT/${Out}/$File.rev-hits
	if [[ -s $address/OUTPUT/${Out}/$File.rev-hits ]]
	then
		sq=`grep ">" $address/OUTPUT/${Out}/$File.rev-hits | wc -l`
		ppm2=`expr 1000000 \* $sq / $reads`
		mv $address/OUTPUT/${Out}/$File.rev-hits $address/OUTPUT/${Out}/$File.rev-hits.fasta
		cd-hit -i $address/OUTPUT/${Out}/$File.rev-hits.fasta -o $address/OUTPUT/${Out}/$File.rev-hits -c 1.00 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5 > $address/OUTPUT/${Out}/cdhitlog # Parameters for CD-Hit (removing redundancy) should be specified here
		rm -rf $address/OUTPUT/${Out}/$File.rev-hits.fasta $address/OUTPUT/${Out}/$File.rev-hits.clstr $address/OUTPUT/${Out}/cdhitlog
	else
		sq="0"
		ppm2="0"
	fi
	echo "$sq" > $address/OUTPUT/${Out}/sq.nmb
	echo "$ppm2" > $address/OUTPUT/${Out}/ppm2.nmb
	cp -r $address/OUTPUT/${Out}/$File.rev-hits $address/OUTPUT/${Out}/read_hits/$File.fasta
	mv $address/OUTPUT/${Out}/$File $address/OUTPUT/${Out}/blast_hits
	echo "13" > $address/CR.step; CFLR="N"
fi
}

PN_SPADES1 () # For proteins and genes, starts SPADES process using high kmers (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "13" ]]; then  echo "Skipping Assembly with High Kmers"; else
	$address/OUTPUT/${Out}/assembly_*
	date -u +%s > $address/OUTPUT/${Out}/datestartspades.tmp
	rm -rf $address/OUTPUT/${Out}/assembly_${Out}_$File
	python $address/Lib/spades/bin/spades.py --threads ${threads} -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $address/OUTPUT/${Out}/read_hits/$File.fasta -o $address/OUTPUT/${Out}/assembly_${Out}_$File > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for proteins and genes should be specified here, using high Kmers
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades.tmp)" | bc -l > $address/OUTPUT/${Out}/spadestime.nmb
	rm -rf $address/OUTPUT/${Out}/logspades.txt ; rm -rf $address/OUTPUT/${Out}/datestartspades.tmp
	echo "14" > $address/CR.step; CFLR="N"
fi
}

PN_SPADES2 () # For proteins and genes, in case the first SPADES process with high kmers didn't work, it retries the assembly process using lower kmers. Soft version skips the first step so it will try only the low kmer assembly.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "13" && `cat CR.step` != "14" ]]; then echo "******Skipping Assembly with Low Kmers"; else
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}_$File/spades_contigs.fasta ]]
	then
		echo "# No need to try with lower kmers "
	else
		echo "# Trying for $File for ${Out} with lower kmers"
		if [[ -s $address/OUTPUT/${Out}/read_hits/$File.fasta ]]
		then
			sq=`cat sq.nmb`
			ppm2=`cat ppm2.nmb`
		else
			sq="0"
			ppm2="0"
		fi
		rm -rf $address/OUTPUT/${Out}/assembly_${Out}_$File
		date -u +%s > $address/OUTPUT/${Out}/datestartspades2.tmp
		python $address/Lib/spades/bin/spades.py --threads ${threads} -k 9,11,13,15,17,19,21,31 --only-assembler -s $address/OUTPUT/${Out}/read_hits/$File.fasta -o $address/OUTPUT/${Out}/assembly_${Out}_$File > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for proteins and genes should be specified here, using lower Kmers
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades2.tmp)" | bc -l > $address/OUTPUT/${Out}/spades2time.nmb
		rm -rf $address/OUTPUT/${Out}/logspades.txt ; rm -rf $address/OUTPUT/${Out}/datestartspades2.tmp
	fi
	echo "15" > $address/CR.step; CFLR="N"
fi
}

PNA () # For proteins and genes, arranges data from SPADES assembly into OUTPUT/${Out}, finding the number of contigs for each subreference.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "15" ]]; then echo "******Skipping calculation of contigs"; else
	if [ -s $address/OUTPUT/${Out}/assembly_${Out}_$File/spades_contigs.fasta ]
	then
		cntg=`grep ">" spades_contigs.fasta | wc -l`
		echo "# SPADES worked on $File for ${Out}, finding $cntg contigs"
		mv $address/OUTPUT/${Out}/assembly_${Out}_$File/spades_contigs.fasta $address/OUTPUT/${Out}/contigs/contigs.$File.fasta
	else
		echo "# The proposed analysis could not continue due to problems in SPADES assembly."
		cntg="0"
		Warnings="WARNING: Did not run SPADES properly"
	fi
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	tar -zcvf $address/OUTPUT/${Out}/contigs/SPADES/SPADES_assembly_$File.tar.gz $address/OUTPUT/${Out}/assembly_${Out}_$File --remove-files
	rm -rf $address/OUTPUT/${Out}/assembly_*
	echo "16" > $address/CR.step; CFLR="N"
fi
}

PNORFs () # For proteins and genes, finds ORFs in each contig (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "16" ]]; then echo "******Skipping ORF finding process"; else
	if [[ -s $address/OUTPUT/${Out}/contigs/contigs.$File.fasta ]]
	then
		echo "# Initiating ORF finding process"
		rm -rf $address/OUTPUT/${Out}/ORFs.$File.fna
		date -u +%s > $address/OUTPUT/${Out}/datestartorffinder.tmp
		perl $address/Lib/bb.orffinder.pl --infile=$address/OUTPUT/${Out}/contigs/contigs.$File.fasta --outfile=$address/OUTPUT/${Out}/ORFs/ORFs.$File.fna --minlen=300 --fasta > $address/OUTPUT/${Out}/logorffinder.txt # Parameters for protein/gene ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffinder.tmp)" | bc -l > $address/OUTPUT/${Out}/orffindertime.nmb
		rm -rf $address/OUTPUT/${Out}/logorffinder ; rm -rf $address/OUTPUT/${Out}/datestartorffinder.tmp
		if [[ -s $address/OUTPUT/${Out}/ORFs/ORFs.$File.fna ]]
		then
			ORFs=`grep ">" $address/OUTPUT/${Out}/ORFs/ORFs.$File.fna | wc -l`
		else
			ORFs="0"
			rm -rf $address/OUTPUT/${Out}/ORFs/ORFs.$File.fna
		fi
		if [[ -s cntg.nmb ]]
		then 
			cntg=`cat cntg.nmb`
		fi
		echo "A total of $ORFs ORFs were found for file $File, from $cntg contigs"
	fi
	echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
	echo "17" > $address/CR.step; CFLR="N"
fi
}

PN_CalcStats () # For proteins and genes, calculates statistics of contigs and ORFs for each subreference (average sizes, total sizes of all ORFs/contigs, standard deviation of sizes and maximum size of contigs and ORFs).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "17" ]]; then echo "******Skipping Calculation of Statistics for Contigs and ORFs"; else
	AvgSizeCntg="0"
	TotalSizeCntg="0"
	StdDevCntg="0"
	AvgSizeORF="0"
	TotalSizeORF="0"
	StdDevORF="0"
	MaxCntg="0"
	MaxORF="0"
	echo "import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'rU')

for rec in SeqIO.parse(FastaFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print name, seqLen

FastaFile.close()" > $address/OUTPUT/countsize.py
	chmod +x $address/OUTPUT/countsize.py
	if [[ -s $address/OUTPUT/${Out}/contigs/contigs.$File.fasta ]]
	then
		rm -rf $address/OUTPUT/${Out}/contigs/*.count; rm -rf $address/OUTPUT/${Out}/contigs/*.cntg; rm -rf $address/OUTPUT/${Out}/contigs/$File.cntg.sizes
		cp $address/OUTPUT/${Out}/contigs/contigs.$File.fasta $address/OUTPUT/${Out}/contigs/contigs.$File.count
		sed -i 's/>.*/>size/g' $address/OUTPUT/${Out}/contigs/contigs.$File.count
		python $address/OUTPUT/countsize.py $address/OUTPUT/${Out}/contigs/contigs.$File.count > $address/OUTPUT/${Out}/contigs/$File.cntg
		rm -rf $address/OUTPUT/${Out}/contigs/contigs.$File.count
		grep "size" $address/OUTPUT/${Out}/contigs/$File.cntg | sed 's/size //g'> $address/OUTPUT/${Out}/contigs/$File.cntg.sizes
		rm -rf $address/OUTPUT/${Out}/contigs/$File.cntg
		AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		if [[ -s $address/OUTPUT/${Out}/ORFs/ORFs.$File.fna ]]
		then
			rm -rf $address/OUTPUT/${Out}/ORFs/*.count; rm -rf $address/OUTPUT/${Out}/ORFs/*.cntg; rm -rf $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes
			cp $address/OUTPUT/${Out}/ORFs/ORFs.$File.fna $address/OUTPUT/${Out}/ORFs/$File.count
			sed -i 's/>.*/>size/g' $address/OUTPUT/${Out}/ORFs/$File.count
			python $address/OUTPUT/countsize.py $address/OUTPUT/${Out}/ORFs/$File.count > $address/OUTPUT/${Out}/ORFs/$File.ORF
			rm -rf $address/OUTPUT/${Out}/ORFs/$File.count
			grep "size" $address/OUTPUT/${Out}/ORFs/$File.ORF | sed 's/size //g' > $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes
			rm -rf $address/OUTPUT/${Out}/ORFs/$File.ORF
			AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
		fi
	fi
	rm -rf $address/OUTPUT/countsize.py
	echo "18" > $address/CR.step; CFLR="N"
fi
}

rDNA_Chimeras ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping Dechimerization"; else
	if [ -s $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta ]
	then
		echo "# Analyzing 16S...\n"
		cp -r $address/OUTPUT/${Out}/assembly_${Out}/spades_contigs.fasta $address/OUTPUT/${Out}
	fi
	tar -zcvf $address/OUTPUT/${Out}/SPADES_results.tar.gz $address/OUTPUT/${Out}/assembly_${Out} --remove-files
	rm -rf $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/uparse.txt $address/OUTPUT/${Out}/otus.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16S.relabeled.fa
	date -u +%s > $address/OUTPUT/${Out}/datestartcluster.tmp
	awk 'BEGIN{RS=">";ORS=""}length($0)>75{print ">"$0}' $address/OUTPUT/${Out}/hits.fasta > $address/OUTPUT/${Out}/16S.fa
	# awk 'BEGIN{RS=">";ORS=""}length($0)>500{print ">"$0}' spades_contigs.fasta > 16Ssize.fasta
	usearch -fastx_relabel $address/OUTPUT/${Out}/16S.fa -prefix "${Out};" -fastaout $address/OUTPUT/${Out}/16S.relabeled.fa -keep_annots -sizeout -threads ${threads}
	usearch -cluster_fast $address/OUTPUT/${Out}/16S.relabeled.fa -id 0.97 --maxaccepts 0 --maxrejects 0 -sizeout -centroids $address/OUTPUT/${Out}/16Snr.fa -threads ${threads} # -uc 16S.clusters.uc
	usearch -cluster_otus $address/OUTPUT/${Out}/16Snr.fa -otus $address/OUTPUT/${Out}/otus.fa -uparseout $address/OUTPUT/${Out}/uparse.txt -relabel "${Out};" -minsize 2 -threads ${threads}
	grep -w "OTU" $address/OUTPUT/${Out}/uparse.txt | grep -vw "chimera" | awk '{ print $1 }' > $address/OUTPUT/${Out}/headers.txt
	python $address/Lib/ext.py $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/16Snr.fa
	mv headers.txt-16Snr.fa $address/OUTPUT/${Out}/16S.nonchimera
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartcluster.tmp)" | bc -l > $address/OUTPUT/${Out}/clustertime.nmb
	rm -rf $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/uparse.txt $address/OUTPUT/${Out}/otus.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16S.relabeled.fa $address/OUTPUT/${Out}/datestartcluster.tmp
	echo "12" > $address/CR.step; CFLR="N"
fi
}

rDNA_Filter2 ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" ]]; then echo "******Skipping GreenFiltering"; else
	echo ""
	if [ -s $address/OUTPUT/${Out}/otu_b6out.tsv ]
	then
		touch $address/OUTPUT/${Out}/otu_b6out.tsv
	else
		rm -rf $address/OUTPUT/${Out}/16S.tsv $address/OUTPUT/${Out}/16S.nm7 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.m7 $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.seq $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16Snr.fasta $address/OUTPUT/${Out}/otu_blast6out.tsv $address/OUTPUT/${Out}/otu_b6out.tsv $address/OUTPUT/${Out}/16S.table
		date -u +%s > $address/OUTPUT/${Out}/datestart16sfilter2.tmp
		if [[ -s $References_Folder/$SubRef ]]
		then
			echo "Using your SubReference $SubRef to search for OTUs."
		else
			echo "Using the full Green Genes database (clustered to 100% identity) to search for OTUs."
			SubRef=nr100_green.cur.udb
		fi
		usearch -usearch_global $address/OUTPUT/${Out}/16S.nonchimera -db $References_Folder/$SubRef -sizein -sizeout -strand plus -id 0.97 -evalue 1e-20 --maxaccepts 5000 --maxrejects 5000 -maxhits 1 -matched $address/OUTPUT/${Out}/16S.m7 -notmatched $address/OUTPUT/${Out}/16S.nm7 -otutabout $address/OUTPUT/${Out}/otu_table.seqidnum -blast6out $address/OUTPUT/${Out}/otu_blast6out.tsv -query_cov 0.9 -threads ${threads} # -uc $address/OUTPUT/${Out}/map.uc # Alterar variavel threads
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestart16sfilter2.tmp)" | bc -l > $address/OUTPUT/${Out}/16sfilter2time.nmb
		rm -rf 	$address/OUTPUT/${Out}/datestart16sfilter2.tmp
		mv $address/OUTPUT/${Out}/otu_blast6out.tsv $address/OUTPUT/${Out}/otu_b6out.tsv
	fi
	sed 's/;size=/	/g' $address/OUTPUT/${Out}/otu_b6out.tsv > $address/OUTPUT/${Out}/16S.tsv
	if [[ -s $address/OUTPUT/${Out}/16S.m8 ]]
	then
		touch $address/OUTPUT/${Out}/16S.m8
	else
		date -u +%s > $address/OUTPUT/${Out}/datestartrenaming.tmp
		if [[ -s $address/OUTPUT/${Out}/16S.table ]]
		then 
			touch $address/OUTPUT/${Out}/16S.table
		else
			sed -i 's/>/>16S_/' $address/OUTPUT/${Out}/16S.m7
			cat $address/OUTPUT/${Out}/16S.tsv > $address/OUTPUT/${Out}/16S.table
		fi
		counter="1"
		OTU=`cat $address/OUTPUT/${Out}/16S.table | wc -l`
		while read B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13; do
			name=`grep ^">$B3 " $References_Folder/16S/Green.list  | sed 's/.*.k__/k__/g'`
			new=">16S_$counter|ID_$B4|AligLength_$B5|eval_$B12|Bitscore_$B13|Size_${B2%;}|$name"
			echo -en "\r"; echo -en "Renaming OTUs ($counter/$OTU)   "
			sed -i "s@^>16S_$B1;size=.*@$new@" $address/OUTPUT/${Out}/16S.m7
			((counter+=1))
			sed -i 1,1d $address/OUTPUT/${Out}/16S.table
		done < $address/OUTPUT/${Out}/16S.table
		echo -e "Finished renaming OTUs"
	fi
	cat $address/OUTPUT/${Out}/16S.m7 > $address/OUTPUT/${Out}/16S.m8
	if [[ -s $address/OUTPUT/${Out}/16S.nm8 ]]
	then
		touch $address/OUTPUT/${Out}/16S.nm8
	else
		UnknownN=`grep ">" 16S.nm7 | wc -l`
		echo "Renaming $UnknownN Unknown sequences."
		awk -F ";" '/^>/{print ">16S_Unknown_" ++i ";" $3 ; next}{print}' < $address/OUTPUT/${Out}/16S.nm7 > $address/OUTPUT/${Out}/16S.nm8
		echo -e "Finished renaming Unknowns."
	fi
	mv $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.fasta
	mv $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/DiscardedSequences_possible16S.fasta
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartrenaming.tmp)" | bc -l > $address/OUTPUT/${Out}/renamingtime.nmb
	rm -rf $address/OUTPUT/${Out}/16S.table $address/OUTPUT/${Out}/16S.m7 $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.nm7 $address/OUTPUT/${Out}/16S.tsv $address/OUTPUT/${Out}/datestartrenaming.tmp
	echo "13" > $address/CR.step; CFLR="N"
fi
}

PCoA_Maker ()
{
	if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "13" ]]; then echo "******Skipping PCoA_Maker analysis"; else
	python $address/Lib/PcoA_maker.py $address/OUTPUT/${Out}/16S.fasta 
	rm -rf $address/OUTPUT/${Out}/16S.fasta
	echo "14" > $address/CR.step; CFLR="N"
fi
}

rDNA_Abundance ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "14" ]]; then echo "******Skipping Taxon Abundance Finding"; else
	head -1 $address/OUTPUT/${Out}/otu_table.seqidnum > $address/OUTPUT/${Out}/otu_table.header
	sed -i '1,1d' $address/OUTPUT/${Out}/otu_table.seqidnum
	sort -S 50% --parallel=${threads} -k1,1 $address/OUTPUT/${Out}/otu_table.seqidnum > $address/OUTPUT/${Out}/otu_table.ord
	join -1 1 -2 1 -o 1.2,2.2 $References_Folder/16S/ID_2_OTU.txt $address/OUTPUT/${Out}/otu_table.ord | sed 's/ /\t/g' | sort -S 50% --parallel=${threads} -V > $address/OUTPUT/${Out}/otu_table.otus
	awk '{seen[$1]+=$2}END{ for (id in seen) print id "\t" seen[id] }' $address/OUTPUT/${Out}/otu_table.otus | sort -S 50% --parallel=${threads} -V > $address/OUTPUT/${Out}/otu_table.sum
	cat $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/otu_table.sum > $address/OUTPUT/${Out}/otu_table.tsv
	rm -rf $address/OUTPUT/${Out}/otu_table.ord $address/OUTPUT/${Out}/otu_table.otus $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/otu_table.sum
	echo "16" > $address/CR.step; CFLR="N"
fi
}

rDNA_TaxonFinding ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "15" ]]; then echo "******Skipping Taxon Analysis"; else
	rm -rf $address/OUTPUT/${Out}/taxons
	mkdir $address/OUTPUT/${Out}/taxons
	while read otunum counts; do
		echo "$counts	$(grep -m 1 "$otunum" $References_Folder/16S/ID_2_Taxon.txt)" >> $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre
	done < $address/OUTPUT/${Out}/otu_table.tsv
	sed -i -e 1,1d $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre
	echo "Counts	SeqIDnumber	OTU	Kingdom	Phylum	Class	Order	Family	Genus	Species" > $address/OUTPUT/${Out}/taxons/associatedtaxons.header
	cat $address/OUTPUT/${Out}/taxons/associatedtaxons.header $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre > $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.tsv
	while read seqidnum counts; do
		echo "$counts	$(grep -m 1 ">${seqidnum} " $References_Folder/16S/Green.list)" >> $address/OUTPUT/${Out}/taxons/headers_counts.pre
	done < $address/OUTPUT/${Out}/otu_table.seqidnum
	echo "Counts	Full_Header" > $address/OUTPUT/${Out}/taxons/counts_fullheader.header
	cat $address/OUTPUT/${Out}/taxons/counts_fullheader.header $address/OUTPUT/${Out}/taxons/headers_counts.pre > $address/OUTPUT/${Out}/taxons/headers_counts.tsv
	rm -rf $address/OUTPUT/${Out}/taxons/*.pre $address/OUTPUT/${Out}/taxons/*.header # ; rm -rf 
	rm -rf $address/OUTPUT/${Out}/otu_table.seqidnum
	echo "15" > $address/CR.step; CFLR="N"
fi
}

CleaningTheMess () # Renames files properly, makes final calculations of statistics and assignment of variables, finishes printing the Log.txt
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" && `cat $address/CR.step` != "16"  && `cat $address/CR.step` != "18" && `cat $address/CR.step` != "16" ]]; then echo "******Skipping Self Organizing Module"; else
	rm -rf $address/Buckets/hits
	case $T1 in
		P|p|N|n)
			if [[ -s $address/OUTPUT/${Out}/log.tmps ]]
			then
				touch $address/OUTPUT/${Out}/log.tmps
			else
				cat $address/OUTPUT/${Out}/tmp1.* > $address/OUTPUT/${Out}/log.tmps
				rm -rf $address/OUTPUT/${Out}/tmp1.*
			fi
			mv $address/OUTPUT/${Out}/*.ft $address/OUTPUT/${Out}/blast_hits
			ls $address/OUTPUT/${Out}/blast_hits/*.ft > $address/OUTPUT/${Out}/blast_hits/list
			if [[ -s $address/OUTPUT/${Out}/blast_hits/list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/blast_hits/list`; do
					mv $file ${file%.f*}.blast_fmt6.tsv
				done
			fi
			rm -rf $address/OUTPUT/${Out}/blast_hits/list

			AvgSizeCntg="0"
			TotalSizeCntg="0"
			StdDevCntg="0"
			MaxCntg="0"
			touch $address/OUTPUT/${Out}/contigs/null000.sizes
			cat $address/OUTPUT/${Out}/contigs/*.sizes > $address/OUTPUT/${Out}/contigs/${Out}.allsizes
			if [[ -s $address/OUTPUT/${Out}/contigs/${Out}.allsizes ]]
			then
				AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
			fi
			echo "$AvgSizeCntg|$TotalSizeCntg|$StdDevCntg|$MaxCntg" > $address/OUTPUT/${Out}/CntgStats.nmb # Calculates all statistics of contigs for the reference as a whole.
			rm -rf $address/OUTPUT/${Out}/contigs/*sizes
			ls $address/OUTPUT/${Out}/contigs/*.fasta > $address/OUTPUT/${Out}/contigs/contigs_list
			if [[ -s contigs_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/contigs/contigs_list`; do
					mv $file ${file%.f*}.fasta
				done
			fi
			rm -rf $address/OUTPUT/${Out}/contigs/contigs_list

			AvgSizeORF="0"
			TotalSizeORF="0"
			StdDevORF="0"
			MaxORF="0"
			touch $address/OUTPUT/${Out}/ORFs/null000.sizes
			cat $address/OUTPUT/${Out}/ORFs/*.sizes > $address/OUTPUT/${Out}/ORFs/${Out}.allsizes
			if [[ -s $address/OUTPUT/${Out}/ORFs/${Out}.allsizes ]]
			then
				AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
			fi
			echo "$AvgSizeORF|$TotalSizeORF|$StdDevORF|$MaxORF" > $address/OUTPUT/${Out}/ORFStats.nmb # Calculates all statistics of ORFs for the reference as a whole.
			rm -rf $address/OUTPUT/${Out}/ORFs/*sizes
			ls $address/OUTPUT/${Out}/ORFs/*.fna > $address/OUTPUT/${Out}/ORFs/orfs_list
			if [[ -s $address/OUTPUT/${Out}/ORFs/orfs_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/ORFs/orfs_list`; do
					mv $file ${file%.f*}.fna
				done
			fi
			rm -rf $address/OUTPUT/${Out}/ORFs/orfs_list
			ls $address/OUTPUT/${Out}/read_hits/*.fasta > $address/OUTPUT/${Out}/read_hits/hits_list
			if [[ -s $address/OUTPUT/${Out}/read_hits/hits_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/read_hits/hits_list`; do
					mv $file ${file%.f*}.hits.fasta
				done
			fi
			rm -rf $address/OUTPUT/${Out}/read_hits/hits_list

			rm -rf $address/OUTPUT/${Out}/log.header $address/OUTPUT/${Out}/log.tmp1; rm -rf $address/OUTPUT/${Out}/cont_log*; rm -rf $address/OUTPUT/${Out}/ORF_log*
			sed -i 's/|/\t/g' $address/OUTPUT/${Out}/log.tmps
			rm -rf $address/OUTPUT/${Out}/cont_log*
			cut -f 4 $address/OUTPUT/${Out}/log.tmps > $address/OUTPUT/${Out}/cont_log$(( n %= 100001))
			if [[ -s $address/OUTPUT/${Out}/cont_log* ]]
			then
				cat $address/OUTPUT/${Out}/cont_log* > $address/OUTPUT/${Out}/c_int
				sed -i '/[a-z]/d' $address/OUTPUT/${Out}/c_int
				cntg=`awk '{s+=$1} END {print s}' $address/OUTPUT/${Out}/c_int`
				echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
			else
				if [[ -s $address/OUTPUT/${Out}/cntg.nmb ]]
				then 
					cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
					if [[ -s $address/OUTPUT/${Out}/CntgStats.nmb ]]
					then
						CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
					fi
				fi
			fi
			rm -rf $address/OUTPUT/${Out}/ORF_log*
			cut -f 9 $address/OUTPUT/${Out}/log.tmps > $address/OUTPUT/${Out}/ORF_log$(( n %= 100001))
			if [[ `cat $address/CR.mode` == "Soft" ]]; 
			then
				ORFs="[Warning: ORFs are not calculated in Soft version]"
				echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
				ORFStats="NA|NA|NA|NA"
				echo "$ORFStats" > $address/OUTPUT/${Out}/ORFStats.nmb
			else
				if [[ -s $address/OUTPUT/${Out}/ORF_log* ]]
				then
					cat $address/OUTPUT/${Out}/ORF_log* > $address/OUTPUT/${Out}/o_int
					sed -i '/[a-z]/d' $address/OUTPUT/${Out}/o_int
					ORFs=`awk '{s+=$1} END {print s}' $address/OUTPUT/${Out}/o_int`
					echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
					if [[ -s $address/OUTPUT/${Out}/ORFStats.nmb ]]
					then
						ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
					fi
				else
					if [[ -s $address/OUTPUT/${Out}/ORFs.nmb ]]
					then 
						ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
						if [[ -s $address/OUTPUT/${Out}/ORFStats.nmb ]]
						then
							ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
						fi
					fi
				fi
			fi
			echo "|Contigs: $cntg
|ORFs: $ORFs
-----------------------------------------------------------------------------
Subref_database|Hits_seq.|Ppm|Contigs|Total_Size_Contigs|Avg_Size_Contigs|Standard_Deviation_Contig_Sizes|Max_Contig_Size|ORFs|Total_Size_ORFs|Avg_ORF_Contigs|Standard_Deviation_ORF_Sizes|Max_ORF_Size|Blast_Time|SPADES_Time|Total_SubRef_Time|Status" > $address/OUTPUT/${Out}/log.header
			cat $address/OUTPUT/${Out}/log.header $address/OUTPUT/${Out}/log.tmps > $address/OUTPUT/${Out}/log.tmp1
			sed -i 's/|/\t/g' $address/OUTPUT/${Out}/log.tmp1
			rm -rf $address/OUTPUT/${Out}/log.header
		;;
		16S|16s|16)
			touch $address/OUTPUT/${Out}/Log.txt
		;;
		G|g)
			cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
			ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
			echo "	Contigs: $cntg
	ORFs: $ORFs" > $address/OUTPUT/${Out}/log.tmp1
		;;
	esac
	rm -rf $address/OUTPUT/${Out}/assembly_*
	rm -rf $address/OUTPUT/${Out}/fulltime.tmp
	if [[ -s $address/OUTPUT/${Out}/spadestime2.nmb ]]
	then
		echo "Time for first SPADES run (failed): $(cat $address/OUTPUT/${Out}/spadestime.nmb)s
Time for second SPADES run: $(cat $address/OUTPUT/${Out}/spades2time.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	else
		if [[ -s $address/OUTPUT/${Out}/spadestime.nmb ]]
		then
			echo "Time for the SPADES run: $(cat $address/OUTPUT/${Out}/spadestime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/quasttime.nmb ]]
	then
		echo "Time for Quast assessment: $(cat $address/OUTPUT/${Out}/quasttime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/orffindertime.nmb ]]
	then
		echo "Time for ORFs to be found by ORFfinder: $(cat $address/OUTPUT/${Out}/orffindertime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/blastdbgentime.nmb ]]
	then
		echo "Time for Blast databases to be created: $(cat $address/OUTPUT/${Out}/blastdbgentime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/filter2time.nmb ]]
	then
		echo "Time the second Usearch filter took to search subreferences: $(cat $address/OUTPUT/${Out}/filter2time.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/clustertime.nmb ]]
	then
		echo "Time for sequences from the first filter to be clusterized: $(cat $address/OUTPUT/${Out}/clustertime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/16sfilter2time.nmb ]]
	then
		echo "Time for the second Usearch filter to find results and generate the OTU table: $(cat $address/OUTPUT/${Out}/16sfilter2time.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/renamingtime.nmb ]]
	then
		echo "Time to rename sequences: $(cat $address/OUTPUT/${Out}/renamingtime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ $TimeLoss == "Y" ]]
	then
		d3="at least $(echo "$($date -u +%s) - $d1" |bc -l)"
	else
		d3=$(echo "$(date -u +%s) - $d1" |bc -l)
	fi
	if [[ -s $address/Buckets/bucketpreviouslygeneratedmessage.tmp ]]
	then
		echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds
	Total time adding time to trim, analyse sequence quality and generate buckets (which were all done in a previous run): : $(echo $d3 + $(cat $address/Buckets/trimmingtime.nmb) + $(cat $address/Buckets/fastqctime.nmb) + $(cat $address/Buckets/bucketenginetime.nmb) | bc -l) seconds" >> $address/OUTPUT/${Out}/fulltime.tmp
	else
		echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds
	Total time excluding time to trim, analyse sequence quality and generate buckets: $(echo $d3 - $(cat $address/Buckets/trimmingtime.nmb) - $(cat $address/Buckets/fastqctime.nmb) - $(cat $address/Buckets/bucketenginetime.nmb) | bc -l) seconds" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	touch $address/OUTPUT/${Out}/log.tmp1
	cat $address/OUTPUT/${Out}/log.tmp1 $address/OUTPUT/${Out}/fulltime.tmp >> $address/OUTPUT/${Out}/Log.txt
	# mv hits hits.fasta
	gzip $address/OUTPUT/${Out}/hits.fasta
	rm -rf $address/OUTPUT/${Out}/hits $address/OUTPUT/${Out}/hits.fasta
	if [[ -s $address/OUTPUT/${Out}/spades_contigs.fasta ]]
	then 
		gzip $address/OUTPUT/${Out}/spades_contigs.fasta
	fi
	case $T1 in 
		G|g)
			rm -rf $address/OUTPUT/${Out}/log.tmp1
		;;
		P|p|N|n)
			mv $address/OUTPUT/${Out}/log.tmp1 $address/OUTPUT/${Out}/SubRefs.tsv
		;;
		16S|16s|16)
			mv $address/OUTPUT/${Out}/spades_contigs.fasta.gz $address/OUTPUT/${Out}/assembled16S.fasta.gz
			rm -rf $address/OUTPUT/${Out}/log.tmp1
		;;
	esac
	rm -rf $address/OUTPUT/${Out}/spades_contigs.fasta $address/OUTPUT/${Out}/log.tmp $address/OUTPUT/${Out}/log.tmps $address/OUTPUT/${Out}/fulltime.tmp $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/*.time2; rm -rf $address/OUTPUT/${Out}/*.rev; rm -rf $address/OUTPUT/${Out}/*.hits ; rm -rf $address/OUTPUT/${Out}/*time.nmb ; rm -rf $address/OUTPUT/${Out}/*time.tmp
	echo "19" > $address/CR.step; CFLR="N"
fi
}

PrintResults () # Prints results of the whole output to the Log.tsv file
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "19" ]]; then echo "******Skipping Print Results"; else
	rm -rf $address/OUTPUT/${Out}/tmp1.*
	if [[ $TimeLoss == "Y" ]]
	then
		d3="at least $(echo "$(date -u +%s) - $d1" |bc -l)"
	else
		d3=$(echo "$(date -u +%s) - $d1" |bc -l)
	fi
	if [[ -s $address/OUTPUT/${Out}/cntg.nmb ]]
	then 
		cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
		if [[ -s $address/OUTPUT/${Out}/CntgStats.nmb ]]
		then
			CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
		fi
	fi
	if [[ `cat $address/CR.mode` == "Soft" ]]; 
	then
		ORFs="[Warning: ORFs are not calculated in Soft version]"
		echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
		ORFStats="NA|NA|NA|NA"
		echo "$ORFStats" > $address/OUTPUT/${Out}/ORFStats.nmb
	else
		if [[ -s $address/OUTPUT/${Out}/ORFs.nmb ]]
		then 
			ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
			if [[ -s $address/OUTPUT/${Out}/ORFStats.nmb ]]
			then
				ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
			fi
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/ppm1.nmb ]]
	then
		ppm1=`cat $address/OUTPUT/${Out}/ppm1.nmb`
	fi
	if [[ -s $address/Buckets/reads.nmb ]]
	then
		reads=`cat $address/Buckets/reads.nmb`
	else
		if [[ -s $address/OUTPUT/${Out}/reads.nmb ]]
		then
			reads=`cat $address/OUTPUT/${Out}/reads.nmb`
		fi
	fi
	echo "${Out}|$R1|$T1|$T2|$Ref|$SubRef|$d3|$reads|$buckets|$ppm1|$cntg|$CntgStats|$ORFs|$ORFStats" >> $address/Log.tsv
	sed -i 's/|/\t/g' $address/Log.tsv
	rm -rf $address/c_int $address/o_int; rm -rf $address/ORF_log*; rm -rf $address/cont_log*
	rm -rf $address/OUTPUT/${Out}/c_int $address/OUTPUT/${Out}/o_int; rm -rf $address/OUTPUT/${Out}/ORF_log*; rm -rf $address/OUTPUT/${Out}/cont_log*

	((d3days=${d3}/86400))
	((d3hours=(${d3}%86400)/3600))
	((d3minutes=((${d3}%86400)%3600)/60))
	((d3seconds=((${d3}%86400)%3600)%60))
	if [[ $d3days -ge 1 ]]
	then
		echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${totaldays}d${totalhours}h${totalminutes}m${totalseconds}s). ] \n"
	else
		if [[ $totalhours -ge 1 ]]
		then
			echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${totalhours}h${totalminutes}m${totalseconds}s). ] \n"
		else
			if [[ $totalminutes -ge 1 ]]
			then
				echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${totalminutes}m${totalseconds}s). ] \n"
			else
				echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds. ] \n"
			fi
		fi
	fi
	sed -i -e 1,1d $address/config.kp
	case $Keep in
		Y|y)
			touch Buckets
		;;
		*)
			rm -rf $address/Buckets
			rm -rf $address/OUTPUT/${Out}/reads.nmb
		;;
	esac
	echo "20" > $address/CR.step; CFLR="N"
fi
}

ErrorRevision () # Finds outputs for which SPADES couldn't find contigs and moves them to an 'Errors' folder.
{
	echo "21" > $address/CR.step; CFLR="N"
	echo "Starting Error Revision"
	while read T1 T2 R1 R2 Ref SubRef Out; do
		if [[ -d $address/OUTPUT ]]
		then
			touch $address/OUTPUT
		else
			mkdir $address/OUTPUT
		fi
		if [[ -d $address/OUTPUT/${Out} ]]
		then
			touch $address/OUTPUT/${Out}
		else
			mkdir $address/OUTPUT/${Out}
			echo "This line was skipped in the config.file" > $address/OUTPUT/${Out}/Error.msg
		fi
	done < $ConfigFile

	rm -rf $address/OUTPUT/list
	ls > $address/OUTPUT/list
	sed -i '/list/d' $address/OUTPUT/list
	sed -i '/Errors/d' $address/OUTPUT/list
	sed -i '/OUTPUT/d' $address/OUTPUT/list
	if [ -d $address/OUTPUT/Errors ]
	then
		mv $address/OUTPUT/Errors $address/OUTPUT/Errors.old
		mkdir $address/OUTPUT/Errors
		mv $address/OUTPUT/Errors.old $address/OUTPUT/Errors
	else
		mkdir $address/OUTPUT/Errors
	fi
	for folder in `cat $address/OUTPUT/list`; do
	        if [[ -s $address/OUTPUT/$folder/spades_contigs.fasta.gz ]]
		then
			touch $address/OUTPUT/$folder/spades_contigs.fasta.gz
	        else
	                if [[ -d $address/OUTPUT/$folder/contigs ]]
	                then
				grep "Contigs: " $address/OUTPUT/$folder/Log.txt | sed 's/Contigs: //' | sed '/[a-z]/d' sed '/[A-Z]/d' > $address/OUTPUT/$folder/test.txt
	                        if [[ $(cat $address/OUTPUT/$folder/test.txt) -gt 0 ]]
				then
					rm -rf $address/OUTPUT/$folder/test.txt
	                        else
					rm -rf $address/OUTPUT/$folder/test.txt
					mv $address/OUTPUT/$folder $address/OUTPUT/Errors
	                        fi
	                else
				if [[ -s $address/OUTPUT/16S.fasta && -s $address/OUTPUT/hits.fasta.gz && -s $address/OUTPUT/otu_b6out.tsv && -s $address/OUTPUT/otu_table.tsv && -s $address/OUTPUT/PCoA.png && -d $address/OUTPUT/taxons ]]
				then
					touch Log.txt
				else
					mv $address/OUTPUT/$folder $address/OUTPUT/Errors
				fi
	                fi
	        fi
		rm -rf $address/OUTPUT/Errors/$folder/test.txt
	done
	rm -rf $address/OUTPUT/list
	ls $address/OUTPUT/Errors/ > $address/OUTPUT/Errors/redolist
	sed -i '/redolist/d' $address/OUTPUT/Errors/redolist
	sed -i '/Errors.old/d' $address/OUTPUT/Errors/redolist
	for out in `cat $address/OUTPUT/Errors/redolist`; do
		grep "$out" $address/config.file >> $address/OUTPUT/Errors/config.redo
	done
	if [[ -s $address/OUTPUT/Errors/config.redo ]]
	then
		ErrorNumber=`cat $address/OUTPUT/Errors/config.redo | wc -l`
		echo "Error Revision complete: BEAF could not find contigs for $ErrorNumber references. These files can either have too few hits or may have suffered errors during the pipeline process, and are now stored in the folder 'Errors' inside the OUTPUT folder. You can find a config file in the same folder that can be used to rerun the process only for these specific files. It is recommended you check these files first ."
	else
		echo "Error Revision complete: BEAF found contigs for every reference it worked on."
	fi
}

	# ======================================================================================================================================================================================== #
	# =====================================================================================PIPELINE BEAF====================================================================================== #
	# ======================================================================================================================================================================================== #


BEAF () # The pipeline BEAF will run when using full version.
{
	echo "##### Running BEAF, full version #####"
	d0=`date -u "+%s"`
	make_kp
	Check
	TimeHeader
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		d1=`date -u "+%s"`
		echo -e "\n# Starting work in file $R1 with $Ref as reference, going to ${Out}\n"
		if [ -d $address/Buckets ]
		then
			ls $address/Buckets/*.bk > $address/Buckets/buckets_list.txt
			if [[ -s $address/Buckets/buckets_list.txt ]];
			then
				echo "# Using previous buckets"
			else
				Trim
				QAnConversion
			fi
		else
			mkdir $address/Buckets
			Trim
			QAnConversion
		fi
		if [[ -d $address/OUTPUT/${Out} ]]
		then
			echo "Folder ${Out} already exists, continuing work..."
		else
			if [[ -d $address/OUTPUT ]]
			then
				touch $address/OUTPUT
			else
				mkdir $address/OUTPUT
			fi
			mkdir $address/OUTPUT/${Out}
			echo "Folder ${Out} created in OUTPUT"
		fi
		BucketEngine
		buckets=`ls $address/Buckets/*.bk | wc -l`
		Filter1
		PreLogGen
		case $T1 in
			G|g)
				if [[ -s $address/Buckets/hits ]]
				then
					G_Prepare_SPADES
					G_SPADES1
					G_SPADES2
					GA
					if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
					then
						gzip $address/OUTPUT/${Out}/hits.fasta
					fi
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/Buckets/hits
					echo "12" > $address/CR.step
				fi
			;;
			P|p|N|n)
				if [[ -s $address/Buckets/hits ]]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					SaveDBs
					Extraction
					if [[ -s $address/Buckets/reads.nmb ]]
					then
						reads=`cat $address/Buckets/reads.nmb`
					else
						if [[ -s $address/OUTPUT/${Out}/reads.nmb ]]
						then
							reads=`cat $address/OUTPUT/${Out}/reads.nmb`
						fi
					fi	
					if [[ "$CFLR" == "Y" ]]
					then
						if [[ -s $address/OUTPUT/${Out}/list ]]
						then
							touch $address/OUTPUT/${Out}/list
						else
							ls $address/OUTPUT/${Out}/*.ft | sed "s@${address}/OUTPUT/${Out}/@@" > $address/OUTPUT/${Out}/list
						fi
					else
						ls $address/OUTPUT/${Out}/*.ft | sed "s@${address}/OUTPUT/${Out}/@@" > $address/OUTPUT/${Out}/list
					fi
					for File in `cat list`; do
						echo "# Working on file ${File} for ${Out}..."
						BTime="0"
						STime="0"
						TTime="0"
						d4=`date -u "+%s"`
						ppm2="0"
						cntg="0"
						sq="0"
						Warnings="[OK]"
						ORFs="0"
						if [ $CFLR == "Y" ]
						then
							if [[ -s $address/OUTPUT/${Out}/ppm2.nmb ]]
							then
								ppm2=`cat $address/OUTPUT/${Out}/ppm2.nmb`
							fi
							if [[ -s $address/OUTPUT/${Out}/cntg.nmb ]]
							then
								cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
							fi
							if [[ -s $address/OUTPUT/${Out}/sq.nmb ]]
							then
								sq=`cat $address/OUTPUT/${Out}/sq.nmb`
							fi
							if [[ -s $address/OUTPUT/${Out}/ORFs.nmb ]]
							then
								ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
							fi
						else
							rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb $address/OUTPUT/${Out}/ORFs.nmb
						fi
						if [[ -s $address/OUTPUT/${Out}/$File ]]
						then
							PN_Prepare_SPADES
							PN_SPADES1
							PN_SPADES2
							PNA
						else
							echo "# ${File} in ${Out} did not reach the minimum criteria to be considered homologus"
							Warnings="WARNING: Did not reach minimum criteria to be considered homologus"
							rm -rf $address/OUTPUT/${Out}/$File
						fi
						PNORFs
						d5=`date -u "+%s"`
						if [[ -s $address/OUTPUT/${Out}/$File.time2 ]]
						then
							BTime=`cat $address/OUTPUT/${Out}/$File.time2`
						fi
						if [[ $TimeLoss == "Y" && $CFLR == "Y" ]]
						then
							STime=$(echo $d5 - $d4 |bc -l)
							TTime="at least $(echo $BTime + $STime |bc -l)"
							STime="at least $(echo $d5 - $d4 |bc -l)"
						else
							STime=$(echo $d5 - $d4 |bc -l)
							TTime=$(echo $BTime + $STime |bc -l)
						fi
						rm -rf $address/OUTPUT/${Out}/*.time2
						if [[ -s $address/OUTPUT/${Out}/sq.nmb ]]
						then
							sq=`cat $address/OUTPUT/${Out}/sq.nmb`
						fi
						if [[ -s $address/OUTPUT/${Out}/ppm2.nmb ]]
						then
							ppm2=`cat $address/OUTPUT/${Out}/ppm2.nmb`
						fi
						if [[ -s $address/OUTPUT/${Out}/cntg.nmb ]]
						then
							cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
						fi
						if [[ -s $address/OUTPUT/${Out}/ORFs.nmb ]]
						then
							ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
						fi
						PN_CalcStats
						echo "${File%.f*}|$sq|$ppm2|$cntg|$TotalSizeCntg|$AvgSizeCntg|$StdDevCntg|$MaxCntg|$ORFs|$TotalSizeORF|$AvgSizeORF|$StdDevORF|$MaxORF|$BTime|$STime|$TTime|$Warnings" > $address/OUTPUT/${Out}/tmp1.$File
						sed -i -e 1,1d $address/OUTPUT/${Out}/list
					done
					rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb $address/OUTPUT/${Out}/ORFs.nmb
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference sequence."
					ppm2="0"
					echo "$ppm2" > $address/OUTPUT/${Out}/ppm2.nmb
					cntg="0"
					echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
					ORFs="0"
					echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
					echo "22" > $address/CR.step
				fi
			;;
			16S|16s|16)
				if [[ -s $address/Buckets/hits ]]
				then
					G_Prepare_SPADES
					G_SPADES1
					G_SPADES2
					rDNA_Chimeras
					if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
					then
						gzip $address/OUTPUT/${Out}/hits.fasta
					fi
					rDNA_Filter2
					PCoA_Maker
					rDNA_Abundance
					rDNA_TaxonFinding
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/Buckets/hits
					echo "16" > $address/CR.step
				fi
			;;
		esac
		CleaningTheMess
		PrintResults
		rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/*.nmb
		CFLR="N"; TimeLoss="N"
	done < $address/config.kp
	rm -rf $address/*.kp $address/*.nmb
	rm -rf $address/par.time
	
	ErrorRevision
	
	CFLR="N"
	rm -rf $address/CR.step $address/CR.mode

	d99=`date -u "+%s"`
	dtotal=$(echo "$d99 - $d0" |bc -l)
	((totaldays=${dtotal}/86400))
	((totalhours=(${dtotal}%86400)/3600))
	((totalminutes=((${dtotal}%86400)%3600)/60))
	((totalseconds=((${dtotal}%86400)%3600)%60))
	if [[ $totaldays -ge 1 ]]
	then
		echo "BEAF1011.65 worked for $dtotal seconds (${totaldays}d${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date "+%X")."
	else
		if [[ $totalhours -ge 1 ]]
		then
			echo "BEAF1011.65 worked for $dtotal seconds (${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date "+%X")."
		else
			if [[ $totalminutes -ge 1 ]]
			then
				echo "BEAF1011.65 worked for $dtotal seconds (${totalminutes}m${totalseconds}s), ending at $(date "+%X")."
			else
				echo "BEAF1011.65 worked for $dtotal seconds (${totalseconds}s), ending at $(date "+%X")."
			fi
		fi
	fi
}

SoftBEAF () # The pipeline BEAF will run when using soft version.
{
	echo "##### Running Soft BEAF, a faster but simpler version of the BEAF software #####"
	d0=`date -u "+%s"`
	make_kp
	Check
	SoftTimeHeader
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		d1=`date -u "+%s"`
		echo -e "\n# Starting work in file $R1 with $Ref as reference, going to ${Out}\n"
		if [[ -d $address/Buckets ]]
		then
			ls $address/Buckets/*.bk > $address/Buckets/buckets_list.txt
			if [[ -s $address/Buckets/buckets_list.txt ]]
			then
				echo "# Using previous buckets"
			else
				CopyFile
				SoftMergeRename
			fi
		else
			mkdir $address/Buckets
			CopyFile
			SoftMergeRename
		fi
		if [[ -d $address/OUTPUT/${Out} ]]
		then
			echo "Folder ${Out} already exists, continuing work..."
		else
			if [[ -d $address/OUTPUT ]]
			then
				touch $address/OUTPUT
			else
				mkdir $address/OUTPUT
			fi
			mkdir $address/OUTPUT/${Out}
			echo "Folder ${Out} created in OUTPUT"
		fi
		BucketEngine
		buckets=`ls $address/Buckets/*.bk | wc -l`
		Filter1
		PreLogGen
		case $T1 in
			G|g)
				if [[ -s $address/Buckets/hits ]]
				then
					G_Prepare_SPADES
					G_SPADES2
					SoftGA
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/Buckets/hits
					echo "12" > $address/CR.step
				fi
			;;
			P|p|N|n)
				if [[ -s $address/Buckets/hits ]]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					SaveDBs
					Extraction
					if [[ -s $address/Buckets/reads.nmb ]]
					then
						reads=`cat $address/Buckets/reads.nmb`
					else
						if [[ -s $address/OUTPUT/${Out}/reads.nmb ]]
						then
							reads=`cat $address/OUTPUT/${Out}/reads.nmb`
						fi
					fi					
					if [[ -s $address/OUTPUT/${Out}/list ]]
					then
						touch $address/OUTPUT/${Out}/list
					else
						ls $address/OUTPUT/${Out}/*.ft | sed "s@${address}/OUTPUT/${Out}/@@" > $address/OUTPUT/${Out}/list
					fi
					for File in `cat $address/OUTPUT/${Out}/list`; do
						echo "# Working on file $File for ${Out}..."
						BTime="0"
						STime="0"
						TTime="0"
						d4=`date -u "+%s"`
						ppm2="0"
						cntg="0"
						sq="0"
						Warnings="[OK]"
						if [ $CFLR == "Y" ]
						then
							if [ -s $address/OUTPUT/${Out}/ppm2.nmb ]
							then
								ppm2=`cat $address/OUTPUT/${Out}/ppm2.nmb`
							fi
							if [ -s $address/OUTPUT/${Out}/cntg.nmb ]
							then
								cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
							fi
							if [ -s $address/OUTPUT/${Out}/sq.nmb ]
							then
								sq=`cat $address/OUTPUT/${Out}/sq.nmb`
							fi
						else
							rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb $address/OUTPUT/${Out}/ORFs.nmb
						fi
						if [ -s $address/OUTPUT/${Out}/$File ]
						then
							PN_Prepare_SPADES
							PN_SPADES2
							PNA
						else
							echo "# $File in ${Out} did not reach the minimum criteria to be considered homologus"
							Warnings="WARNING: Did not reach minimum criteria to be considered homologus"
							rm -rf $address/OUTPUT/${Out}/$File
						fi
						d5=`date -u "+%s"`
						if [ -s $address/OUTPUT/${Out}/$File.time2 ]
						then
							BTime=`cat $address/OUTPUT/${Out}/$File.time2`
						fi
						if [[ $TimeLoss == "Y" && $CFLR == "Y" ]]
						then
							STime=$(echo $d5 - $d4 |bc -l)
							TTime="at least $(echo $BTime + $STime |bc -l)"
							STime="at least $(echo $d5 - $d4 |bc -l)"
						else
							STime=$(echo $d5 - $d4 |bc -l)
							TTime=$(echo $BTime + $STime |bc -l)
						fi
						rm -rf $address/OUTPUT/${Out}/*.time2
						if [[ -s $address/OUTPUT/${Out}/sq.nmb ]]
						then
							sq=`cat $address/OUTPUT/${Out}/sq.nmb`
						fi
						if [[ -s $address/OUTPUT/${Out}/ppm2.nmb ]]
						then
							ppm2=`cat $address/OUTPUT/${Out}/ppm2.nmb`
						fi
						if [[ -s $address/OUTPUT/${Out}/cntg.nmb ]]
						then
							cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
						fi
						echo "${File%.f*}|$sq|$ppm2|$cntg|$TotalSizeCntg|$AvgSizeCntg|$StdDevCntg|$MaxCntg|NA|NA|NA|NA|NA|$BTime|$STime|$TTime|$Warnings" > $address/OUTPUT/${Out}/tmp1.$File
						sed -i -e 1,1d $address/OUTPUT/${Out}/list
					done
					rm -rf $address/OUTPUT/${Out}/ORFs; rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference sequence."
					ppm2="0"
					echo "$ppm2" > $address/OUTPUT/${Out}/ppm2.nmb
					cntg="0"
					echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
					echo "21" > $address/CR.step
				fi
			;;
		esac
		CleaningTheMess
		PrintResults
		CFLR="N"; TimeLoss="N"
		rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/*.nmb
	done < $address/config.kp
	rm -rf $address/*.kp ; rm -rf $address/*.nmb 
	rm -rf $address/par.time
	
	ErrorRevision

	CFLR="N"
	rm -rf $address/CR.step $address/CR.mode

	d99=`date -u "+%s"`
	dtotal=$(echo "$d99 - $d0" |bc -l)
	echo "BEAF1011.65 worked for $dtotal seconds, ending at $(date "+%X")."
}

	# ======================================================================================================================================================================================== #
	# ========================================================================================MAIN============================================================================================ #
	# ======================================================================================================================================================================================== #

Main () # Interpreting user's input of whether he'll be using soft or full version.
{
	case $1 in
		S|s|-s|-S|Soft|soft|SOFT|-Soft|-soft|-SOFT|"SOFT")
			echo "Soft" > $address/CR.mode
			SoftBEAF
		;;
		*)
			echo "FullVersion" > $address/CR.mode
			BEAF
		;;
	esac
}

	# ======================================================================================================================================================================================== #
	# =====================================================================================PROGRAM START====================================================================================== #
	# ======================================================================================================================================================================================== #

# This is what will actually be running once you start BEAF, asking whether to use Soft or Full version, then redirecting to Main function, which will interpret the answer
# Once that is interpreted, it will start the proper pipeline, either BEAF or SoftBEAF, which will then run in a modular fashion, calling each function as it runs.

while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help|-help|\?|-\?)
			# show_help
			exit
		;;
		-s|-S|-Soft|-soft|-SOFT|--Soft|--soft|--SOFT|--Soft-Version|--soft-version|--Soft-version)
			ver="SOFT"
		;;
		--config|--Config|--CONFIG|-config|-Config|--CONFIG)
			if [[ -s ${2} ]]
			then
				ConfigFile=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
				echo "Using file $ConfigFile as config.file"
			else
				echo "Couldn't find the designated config.file '${2}'"
				exit
			fi
			shift
		;;
		--config=*|--Config=*|--CONFIG=*|-config=*|-Config=*|--CONFIG=*)
			if [[ -s ${1#*=} ]]
			then
				ConfigFile=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
				echo "Using file $ConfigFile as config.file"
			else
				echo "Couldn't find the designated config.file '${1#*=}'"
				exit
			fi
		;;
		--output|--Output|--OUTPUT|-output|-Output|-OUTPUT|-O|-o|-Out|-OUT|-out|--out|--Out|--OUT|--o|--O)
			if [[ -d "${2}" ]]
			then
				touch ${2}
			else
				mkdir ${2}
			fi
			address=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
			echo "Using folder $address as output folder"
			shift
		;;
		--output=*|--Output=*|--OUTPUT=*|-output=*|-Output=*|-OUTPUT=*|-O=*|-o=*|-Out=*|-OUT=*|-out=*|--out=*|--Out=*|--OUT=*|--o=*|--O=*)
			if [[ -d "${1#*=}" ]]
			then
				touch ${1#*=}
			else
				mkdir ${1#*=}
			fi
			address=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
			echo "Using folder $address as output folder"
		;;
		--force_continue|--Force_Continue|--Force_continue|--force_Continue|--FORCE_CONTINUE|--fc|--FC|--Fc|--fC|-fc|-FC|-Fc|-fC)
			CFLR="Y"
		;;
		--force_restart|--Force_Restart|--Force_restart|--force_Restart|--FORCE_RESTART|--forcerestart|--Forcerestart|--ForceRestart|--FORCERESTART|--fr|--FR|--Fr|--fR|-fr|-FR|-Fr|-fR)
			CFLR="N"
		;;
		-t|-T|--threads|--Threads|--THREADS|--t|--T)
			threads=${2}
			shift
		;;
		-t=*|-T=*|--threads=*|--Threads=*|--THREADS=*|--t=*|--T=*)
			threads=${1#*=}
		;;
	esac
	shift
done

if [[ $CFLR == "Y" ]]
then
	CFLR="Y"
	TimeLoss="Y"
	echo "BEAF will continue from last run."
	ver=`cat CR.mode`
	Main $ver
else
	if [[ $CFLR == "N" ]]
	then
		if [[ -s $address/CR.step ]]
		then
			echo "Deleting old files and starting a new run."
		else
			if [[ -s $address/config.kp ]]
			then
				echo "Deleting old files and starting a new run."
			fi
		fi
		rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/Buckets
		Main $ver
	else
		if [[ -s $address/CR.step && $CFLR == "U" ]]
		then
			echo "A file from your previous run was detected, indicating that your last run of BEAF was interrupted before the program finished. BEAF is capable of detecting at which step the program was interrupted and continuing interrupted runs, but note that if you removed files from the BEAF folder after the program was abruptly interrupted, the new run can be disrupted. Do you want to continue from where the program stopped? (Y/N)"
			read CFLR # continue from last run
			case $CFLR in
				Y|y|Yes|yes|continue|Continue)
					CFLR="Y"
					TimeLoss="Y"
					echo "BEAF will continue from last run."
					ver=`cat CR.mode`
					Main $ver
				;;
				*)
					CFLR="N"
					rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/Buckets
					# echo "BEAF will start a new run. 

		# How do you want to run BEAF?
		# -s: run Soft BEAF, a faster but simpler version of the BEAF software
		# -b: run full version, with all utilities"
		#			read ver
					Main $ver
				;;
			esac
		else
			if [[ -s $address/config.kp && $CFLR == "U" ]]
			then
				echo "A config.kp file from your previous run was detected, indicating that your last run of BEAF was interrupted before the program finished. Although BEAF couldn't detect at which point the program was interrupted (because the file 'CR.step' was removed), it was able to determine which OUTPUTs were already generated. BEAF is capable of continuing interrupted runs, but note that if you removed files from the BEAF folder after the program was abruptly interrupted, the new run can be disrupted. Do you want to continue from where the program stopped? (Y/N)"
				read CFLR # Should I Continue Last Run?
				case $CFLR in
					Y|y|Yes|yes|continue|Continue)
						CFLR="Y" # Continue From Last Run
						TimeLoss="Y"
						echo "1" > $address/CR.step
						echo "BEAF will continue from last run."
						ver=`cat CR.mode`
		#				read ver
						Main $ver
					;;
					*)
						CFLR="N" # Do not Continue From Last Run
						rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/Buckets
		#				echo "BEAF will start a new run. 
		#	
		#	How do you want to run BEAF?
		#	-s: run Soft BEAF, a faster but simpler version of the BEAF software
		#	-b: run full version, with all utilities"
		#		read ver
						Main $ver
					;;
				esac
			else
		#		echo "How do you want to run BEAF?
		# -s: run Soft BEAF, a faster but simpler version of the BEAF software
		# -b: run full version, with all utilities"
		#		read ver
				rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/Buckets
				Main $ver
			fi
		fi
	fi
fi
