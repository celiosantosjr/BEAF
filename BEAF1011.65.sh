#!/usr/bin/env bash

d0=`date -u "+%s"`
address=$(dirname $(readlink -f $0))
ans="no" #If you're sure you want to keep databases files from one loop to the next, change this to 'Y'. If you're sure you're not reusing .udb, change to 'N'

echo "# Your program is currently located in \"$address\""
echo "# Checking if any buckets must be stored..."
echo "X" > LastR1.kp
sort -k3,3 config.file > config.kp
while read T1 T2 R1 R2 Ref SubRef Out; do
	LastR1=`cat LastR1.kp`	
	if [[ "$R1" == "$LastR1" ]]
	then
		echo "Y" > LastKeep.kp
	else
		echo "N" > LastKeep.kp
	fi
	cat Keep_config.kp LastKeep.kp > Keep_config2.kp
	mv Keep_config2.kp Keep_config.kp
	echo "$R1" > LastR1.kp
done < config.kp
echo "N" > LastKeep.kp
cat Keep_config.kp LastKeep.kp > Keep_config2.kp
mv Keep_config2.kp Keep_config.kp
sed -i -e 1,1d Keep_config.kp
paste config.kp Keep_config.kp > config.tmp
rm -rf *.kp
mv config.tmp config.kp

echo "# Checking your config file..."

sort -k7,7 config.kp > config.check
while read T1 T2 R1 R2 Ref SubRef Out Keep; do
	cd $address
	case $T1 in
		G|g|P|p|N|n)
		;;
		*)
			echo "# You're using a wrong config.file format. On the first column (T1), where you're currently using '$T1', use only G (for genome analysis), P (for protein analysis) or N (for protein nucleotide sequences analysis)"
			rm -rf config.kp config.check
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
						rm -rf config.kp config.check
						exit
					fi
				;;
				*)
					echo "# You're using a wrong config format. On the fourth column (R2), where you're currently using '$R2', you must use a gzipped file instead."
					rm -rf config.kp config.check
					exit
				;;
			esac
		;;
		I|i|F|f)
		;;
		*)
			echo "# You're using a wrong config.file On the second column (T2), where you're currently using '$T2', use only R (for paired end fastq files), I (for interleaved fastq file) or F (for interleaved fasta file)"
			rm -rf config.kp config.check
			exit
		;;
	esac
	case $R1 in
		*.gz)
			if ! [ -s $R1 ];
			then
				echo "# Check your R1 file in '$R1'. The program either couldn't find the file or the file is empty."
				rm -rf config.kp config.check
				exit
			fi
		;;
		*)
			echo "# You're using a wrong config.file On the third column (R1), where you're currently using '$R1', you must use a gzipped file instead."
			rm -rf config.kp config.check
			exit
		;;
	esac
	LastOut=`cat LastOut.check`
	if [[ "$Out" == "$LastOut" ]]
		then
			echo "# You're using a wrong config.file On your seventh column (Out), you've used the same name for your output folder more than once, repeating '$Out'"
			exit
	fi
	echo "$Out" > LastOut.check
done < config.check
rm -rf *.check

echo "________________________________________________________________________________________________
Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs" > header
cat time.log header > time.log2
rm -rf header time.log
mv time.log2 time.log

while read T1 T2 R1 R2 Ref SubRef Out Keep; do
	d1=`date -u "+%s"`
	cd $address
	echo "
# Starting work in file $R1 with $Ref as reference, going to $Out
"
	case $T2 in 
		R|r)
			ls *.bk > buckets_list.txt
			if [ -s buckets_list.txt ];
			then
				echo "# Using previous buckets"
			else
				rm -rf *.gz
				echo "# Copying file 1 from storage..."
				cp $R1 $address
				echo "# Copying file 2 from storage..."
				cp $R2 $address
				echo "# Now we are trimming your files..."
				echo "# Starting file 1..."
				cutadapt --minimum-length 80 --quality-base 24 --trim-n -o R1.trimmed.fastq.gz $R1
				mv R1.trimmed.fastq.gz R1.trimmed.fastq.q
				echo "# Starting file 2..."
				cutadapt --minimum-length 80 --quality-base 24 --trim-n -o R2.trimmed.fastq.gz $R2
				mv R2.trimmed.fastq.gz R2.trimmed.fastq.q
				echo "# Removing intermediate files..."
				rm -rf *.gz
				mv R1.trimmed.fastq.q R1.trimmed.fastq.gz; mv R2.trimmed.fastq.q R2.trimmed.fastq.gz
				echo "# Now we are merging files, it could take several minutes away..."
				zcat *.gz | gzip -c > FastaQ-zcat.gz
				mv FastaQ-zcat.gz FastaQ-zcat.tmp
				rm -rf *.gz
				mv FastaQ-zcat.tmp FastaQ-zcat.gz
				echo "# Starting quality assessment of trimming..."
				mkdir FASTQCresults
				fastqc -f fastq -o FASTQCresults FastaQ-zcat.gz
				cd FASTQCresults
				rm -rf *_fastqc
				cd $address
				echo "# We will convert merged file to fasta format."
				gunzip -c FastaQ-zcat.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FastaQ-zcat.fa
				echo "# Removing unwanted files..."
				rm -rf *.gz
			fi
		;;
		I|i)
			ls *.bk > buckets_list.txt
			if [ -s buckets_list.txt ];
			then
				echo "# Using previous buckets"
			else
				rm -rf *.gz
				echo "# Copying file from storage..."
				cp $R1 $address
				echo "# Now we are trimming your files..."
				echo "# Starting file 1..."
				cutadapt --minimum-length 80 --quality-base 24 --trim-n -o I.trimmed.fastq.gz $R1
				mv I.trimmed.fastq.gz I.trimmed.fastq.q
				echo "# Removing intermediate files..."
				rm -rf *.gz
				mv I.trimmed.fastq.q I.trimmed.fastq.gz
				echo "# Starting quality assessment of trimming..."
				mkdir FASTQCresults
				fastqc -f fastq -o FASTQCresults I.trimmed.fastq.gz
				cd FASTQCresults
				rm -rf *_fastqc
				cd $address
				echo "# We will convert merged file to fasta format."
				gunzip -c I.trimmed.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FastaQ-zcat.fa
				echo "# Removing unwanted files..."
				rm -rf I.trimmed.fastq.gz
			fi
		;;
		F|f)
			ls *.bk > buckets_list.txt
			if [ -s buckets_list.txt ];
			then
				echo "# Using previous buckets"
			else
				echo "# Copying file from storage..."
				cp $R1 $address
				gunzip -c <`ls *.gz`> FastaQ-zcat.fa
				rm -rf *.gz
			fi
		;;
	esac
	if [ -s buckets_list.txt ];
	then
		buckets=`ls *.bk | wc -l`
		reads=`cat reads.nmb`
	else
		grep ">" FastaQ-zcat.fa | wc -l > reads.nmb
		reads=`cat reads.nmb`
		buckets=`expr $reads / 4500000`
		if [ "$buckets" -eq "0" ];
		then
			echo "# We identified $reads reads. It would take no buckets, avoiding this step."
			mv FastaQ-zcat.fa unique.bk
		else
			echo "# We identified $reads reads. It would take $buckets buckets step."
			echo "### Starting operation of cutting and readapting..."
			echo "## Generating buckets..."
			shuf -i 0-1000000 -n $buckets > bk_list.txt
			echo "# Sorting buckets and setting up configurations..."
			for n in `cat bk_list.txt`; do
				head -9000000 FastaQ-zcat.fa > $n.bk
				sed -i '1,9000000d' FastaQ-zcat.fa
			done
			if [ -s FastaQ-zcat.fa ];
			then
				mv FastaQ-zcat.fa last.bk
			fi
		fi
		echo "### Removing temporary files (stage 1)..."
		rm -rf *.txt; rm -rf *.fa; rm -rf *.gz
		ls *.bk > buckets_list.txt
		buckets=`ls *.bk | wc -l`
	fi
	case $T1 in
		G|g)	
			echo "# Starting searches..."
			echo "### Making database from reference genome... "
			cd $address/Reference_seqs
			cp $Ref none
			sed -i '/>/d' none
			cat none | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > md8
			rm -rf none
			usearch -makeudb_usearch md8 -output $Ref.udb
			rm -rf md8
			mv $Ref.udb $address
			cd $address
			case $Keep in
				Y|y)
					for buck in `cat buckets_list.txt`; do
						usearch -usearch_global $buck -db $Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $buck.m8
					done
				;;
				*)
					for buck in `cat buckets_list.txt`; do
						usearch -usearch_global $buck -db $Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $buck.m8
						rm -rf $buck 
					done
					rm -rf *.bk
				;;
			esac
			rm -rf *.udb
		;;
		P|p|N|n)
			cd $address/Reference_seqs
			mkdir TestExtension
			cp $Ref TestExtension
			cd $address/Reference_seqs/TestExtension
			ls *.fa > list.fa; ls *.fasta > list.fasta; ls *.fas > list.fas; ls *.faa > list.faa; ls *.fna > list.fna; ls *.fsa > list.fsa
			cat list.f* > list; rm -rf list.f*
			if [ -s list ]
			then
				echo "# Recognized $Ref file as fasta format. Making udb..."
				sed -i '/>/d' $Ref
				cat $Ref | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > md8
				usearch -makeudb_usearch md8 -output $Ref.udb
				rm -rf md8
				mv $Ref.udb $address
				cd $address/Reference_seqs
				rm -rf TestExtension
				cd $address
				case $Keep in
					Y|y)
						for buck in `cat buckets_list.txt`; do
							usearch -usearch_local $buck -db $Ref.udb -strand both -id 0.25 -evalue 1e-5 -matched $buck.m8
						done
					;;
					*)
						for buck in `cat buckets_list.txt`; do
							usearch -usearch_local $buck -db $Ref.udb -strand both -id 0.25 -evalue 1e-5 -matched $buck.m8
							rm -rf $buck
						done
						rm -rf *.bk
					;;
				esac
			else
				echo "# Recognized $Ref file as .udb format."
				mv $Ref $address
				cd $address/Reference_seqs
				rm -rf TestExtension
				cd $address
				case $Keep in
					Y|y)
						for buck in `cat buckets_list.txt`; do
							usearch -usearch_local $buck -db $Ref -strand both -id 0.25 -evalue 1e-5 -matched $buck.m8
						done
					;;
					*)
						for buck in `cat buckets_list.txt`; do
							usearch -usearch_local $buck -db $Ref -strand both -id 0.25 -evalue 1e-5 -matched $buck.m8
							rm -rf $buck
						done
						rm -rf *.bk
					;;
				esac
			fi
			rm -rf *.udb
		;;
	esac
	cat *.m8 > hits
	echo "### Removing temporary files (stage 2)..."
	rm -rf *.m8; rm -rf *.txt
	echo "## Calculating statistics..."
	hits=`grep ">" hits | wc -l`
	reads=`cat reads.nmb`
	ppm1=`expr 1000000 \* $hits / $reads`
	echo -e "RESULTS:
	File from: $R1\t$R2
	Results: $Out
	Reference: $Ref
	Reads: $reads
	Buckets: $buckets
	Hits: $hits
	Portion in ppm: $ppm1" > Log.txt
	mkdir $Out
	mkdir OUTPUT
	mv $Out OUTPUT
	mv Log.txt OUTPUT/$Out
	mv FASTQCresults OUTPUT/$Out
	case $T1 in
		G|g)
			if [ -s hits ]
			then
				echo "# Making contigs for $Out..."
				mv hits hits.fasta
				cp hits.fasta $address/spades/bin
				mv hits.fasta $address/OUTPUT/$Out
				cd $address/spades/bin
				rm -rf assembly_*
				python spades.py -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s hits.fasta -o assembly_$Out
				rm -rf hits.fasta
				cd $address/spades/bin/assembly_$Out
				if [ -s scaffolds.fasta ]
				then
					echo "SPADES ran properly with high kmers"
				else
					cd $address/OUTPUT/$Out
					echo "# Retrying for $Out with lower kmers"
					rm -rf $address/spades/bin/assembly_*
					cp hits.fasta $address/spades/bin
					cd $address/spades/bin
					rm -rf assembly_*
					python spades.py -k 11,15,21,25,31,35,41,45,51 --only-assembler -s hits.fasta -o assembly_$Out
					rm -rf hits.fasta
				fi
				cd $address/spades/bin/assembly_$Out
				if [ -s scaffolds.fasta ]
				then
					echo "# Analyzing draft putative genome..."
					cp scaffolds.fasta $address/quast
					cd $address/spades/bin
					rm -rf assembly_$Out
					cd $address/quast
					python metaquast.py -R $address/Reference_seqs/$Ref -o assessment scaffolds.fasta
					mv scaffolds.fasta $address/OUTPUT/$Out
					echo "### Compressing results..."
					tar -zcvf assessment.tar.gz assessment --remove-files
					mv assessment.tar.gz $address/OUTPUT/$Out/assessment.tar.gz
					echo "# QUAST service is finished for the file going to OUTPUT/$Out"
					cd $address/OUTPUT/$Out
					perl $address/bb.orffinder.pl --infile=scaffolds.fasta --outfile=ORFs.$Out.fna --minlen=600 --fasta
					cntg=`grep ">" scaffolds.fasta | wc -l`
					gzip ORFs.$Out.fna
					gzip scaffolds.fasta
				else
					echo "# The proposed analysis of $Out could not continue due to problems in SPADES assembly."
				fi
				cd $address/spades/bin
				rm -rf assembly_$Out
				cd $address/OUTPUT/$Out
				gzip hits.fasta
			else
				echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
			fi
		;;
		P|p|N|n)
			if [ -s hits ]
			then
				mkdir $address/OUTPUT/$Out
				echo "# Submitting to Blast per subreference family..."
				mv hits Reference_seqs/$SubRef
				cd $address/Reference_seqs/$SubRef
				ls *.fa > list.fa; ls *.fasta > list.fasta; ls *.fas > list.fas; ls *.faa > list.faa; ls *.fna > list.fna; ls *.fsa > list.fsa
				cat list.f* > glist; rm -rf list.f*
				if [ -s glist ]
				then
					echo "# Recognized $SubRef files in fasta format. Making blast databases..."
					for sub in `cat glist`; do
						case $T1 in
							P|p)
								makeblastdb -in $sub -dbtype prot -out $sub.db
							;;
							N|n)
								makeblastdb -in $sub -dbtype nucl -out $sub.db
							;;
						esac
					done
				fi
				case $T1 in
					P|p)
						ls *.psq | sed 's/.psq//' | sort -k1,1 > list
					;;
					N|n)
						ls *.nsq | sed 's/.nsq//' | sort -k1,1 > list
					;;
				esac
				for sub in `cat list`; do
					echo "# Searching against $sub..."
					case $T1 in
						P|p)
							dsp1=`date -u "+%s"`
							blastx -db $sub -query hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads 4 -outfmt 6
							dsp2=`date -u "+%s"`
							echo $dsp2 - $dsp1 |bc -l > $sub.ft.time2
						;;
						N|n)
							dsp1=`date -u "+%s"`
							blastn -db $sub -query hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads 4 -outfmt 6
							dsp2=`date -u "+%s"`
							echo $dsp2 - $dsp1 |bc -l > $sub.ft.time2
						;;
					esac
					cat $sub.tmp | sort -k3,3 -k4,4 -n -r | awk '$3 > 90 && $4 > 25' | uniq > $sub.ft
					touch $sub.ft
					rm -rf $sub.tmp
				done
				case $ans in
					Y|y|YES|Yes|yes)
						echo "# Databases of subreference $SubRef now saved to References_seqs folder in blastdb format. Fasta files used to make the databases have been realocated to Reference_seqs/$SubRef/Fasta_files."
						mkdir Fasta_files
						for file in `cat glist`; do
							mv $file Fasta_files
						done
						rm -rf glist
					;;
					*)
						"# Removing blastdb databases generated using fasta files in the subreference folder (Reference_seqs/$SubRef)"
						for file in `cat glist`; do
							case $T1 in
								P|p)
									rm -rf $file.db.p*
								;;
								N|n)
									rm -rf $file.db.n*
								;;
							esac
						done
						rm -rf glist
					;;
				esac
				rm -rf list
				echo "# Arranging data..."
				mv *.ft $address/OUTPUT/$Out
				mv hits $address/OUTPUT/$Out
				mv *.time2 $address/OUTPUT/$Out
				cd $address/OUTPUT/$Out
				echo "#exclamationmark/usr/bin/python

import string
import sys
ListOfIds = sys.argv[1]
fastafile = sys.argv[2]

try:
    ids=open(ListOfIds, 'r')
except IOError, e:
    print \"File error: \",ListOfIds
    pass


lignes = ids.readlines()
req=[]
for ligne in lignes:
    req.append(ligne.strip())

handle = open(fastafile)

bank={}
seqIDmap={}
seq_id = handle.next()
while (seq_id[0]!=\">\"):
    seq_id = handle.next()
while True:
    try:
        seq = handle.next()
        line = handle.next()
        while (line[0]!=\">\"):
            seq = seq+line
            line = handle.next()
        bank[seq_id]=seq
        IDclean=string.split(seq_id, \" \")[0][1:].strip()
        seqIDmap[IDclean]=seq_id
        seq_id = line # for the next
    except StopIteration:
        break
# last line
bank[seq_id]=seq
seqIDmap[string.split(seq_id, \" \")[0][1:].strip()]=seq_id

handle.close()

faName=fastafile.split(\"/\")[-1]
listName=ListOfIds.split(\"/\")[-1]
subsetName=listName+\"-\"+faName
subset = open(subsetName,\"w\")
nbNF=0
for i in req:
    try:
        subset.write(seqIDmap[i].strip()+\"\substn\")
        subset.write(bank[seqIDmap[i]].strip()+\"\substn\")
    except KeyError:
        print i, \"not found in fasta\"
        nbNF+=1

subset.close()

print
print nbNF, \"IDs (listed above) from\",listName, \"have not been found in\", faName
print
print \"the Subset fasta file\", subsetName, \"is now created\"" > ext.py
				sed -i 's/exclamationmark/!/' ext.py
				sed -i 's/substn/n/' ext.py
				mkdir read_hits
				mkdir contigs
				mkdir blast_hits
				mkdir ORFs
				ls *.ft > list
				for File in `cat list`; do
					cd $address/OUTPUT/$Out
					BTime="0"
					STime="0"
					TTime="0"
					echo "# Working on file $File for $Out..."
					d4=`date -u "+%s"`
					ppm2="0"
					cntg="0"
					sq="0"
					Warnings="[OK]"
					if [ -s $File ]
					then
						cd $address/OUTPUT/$Out
						cut -f 1 $File > $File.rev
						python ext.py $File.rev hits > $File.hits
						rm -rf $File.hits
						sq=`grep ">" $File.rev-hits | wc -l`
						ppm2=`expr 1000000 \* $sq / $reads`
						cp $File.rev-hits $address/spades/bin/$File.fasta
						mv $File.rev-hits read_hits/$File.fasta
						cd $address/spades/bin
						rm -rf $address/spades/bin/assembly_*
						python spades.py -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $File.fasta -o assembly_$Out_$File
						cd $address/spades/bin/assembly_$Out_$File
						if [ -s scaffolds.fasta ]
						then
							echo "# No need to try with lower kmers "
						else
							echo "# No contigs were found. Retrying for $File for $Out with lower kmers"
							cd $address/spades/bin
							sq=`grep ">" $File.fasta | wc -l`
							ppm2=`expr 1000000 \* $sq / $reads`
							python spades.py -k 9,11,13,15,17,19,21,31 --only-assembler -s $File.fasta -o assembly_$Out_$File
						fi
						cd $address/spades/bin/assembly_$Out_$File
						if [ -s scaffolds.fasta ]
						then
							cntg=`grep ">" scaffolds.fasta | wc -l`
							echo "# SPADES worked on $File for $Out, finding $cntg contigs"
							mv scaffolds.fasta contigs.$File.fasta
							mv contigs.$File.fasta $address/OUTPUT/$Out/contigs
							cd $address/spades/bin
							rm -rf assembly_$Out_$File; rm -rf $File.fasta
						else
							echo "# The proposed analysis could not continue due to problems in SPADES assembly."
							cd $address/spades/bin
							rm -rf assembly_$Out_$File; rm -rf $File.fasta
							cd $address/OUTPUT/$Out 
							Warnings="WARNING: Did not run SPADES properly"
						fi
					else
						echo "# $File in $Out did not reach the minimum criteria to be considered homologus"
						cd $address/OUTPUT/$Out
						Warnings="WARNING: Did not reach minimum criteria to be considered homologus"
						rm -rf $File
					fi
					cd $address/OUTPUT/$Out/contigs
					if [ -s contigs.$File.fasta ]
					then
						perl $address/bb.orffinder.pl --infile=contigs.$File.fasta --outfile=ORFs.$File.fna --minlen=600 --fasta
						ORFs=`grep ">" ORFs.$File.fna | wc -l`
						mv ORFs.$File.fna $address/OUTPUT/$Out/ORFs
					else
						ORFs="0"
					fi
					d5=`date -u "+%s"`
					cd $address/OUTPUT/$Out
					STime=$(echo $d5 - $d4 |bc -l)
					if [ -s $File.time2 ]
					then					
						BTime=`cat $File.time2`
					fi
					TTime=$(echo $BTime + $STime |bc -l)
					rm -rf *.time2
					cd $address/OUTPUT/$Out
					echo "${File/.f*/}|$sq|$ppm2|$cntg|$ORFs|$BTime|$STime|$TTime|$Warnings" > tmp1.$File
				done
				cd $address/OUTPUT/$Out
				rm -rf *.time2
				mv hits hits.fasta
				gzip hits.fasta
				rm -rf $address/spades/bin/assembly_*
				echo "-----------------------------------------------------------------------------
Subref_database|Hits_seq.|Ppm|Contigs|ORFs|Blast_Time|SPADES_Time|Total_SubRef_Time|Status" > log.header
				sed -i 's/|/\t/g' log.header
				cat log.header tmp1.* > log.tmp1
				sed -i 's/|/\t/g' log.tmp1
				cut -f 4 log.tmp1 > cont_log$(( n %= 100001))
				mv cont_log* $address
				rm -rf ext.py; rm -rf *.rev; rm -rf *.hits; rm -rf list; rm -rf tmp1.*; rm -rf log.header
				mv *.ft blast_hits
				cd $address/OUTPUT/$Out/blast_hits
				ls *.ft > list
				for file in `cat list`; do
					mv $file ${file/.f*/.blast.tsv}
				done
				rm -rf list
				cd $address/OUTPUT/$Out/contigs	
				ls *.fasta > contigs_list
				for file in `cat contigs_list`; do
					mv $file ${file/.f*/.fasta}
				done
				rm -rf contigs_list
				cd $address/OUTPUT/$Out/ORFs
				ls *.fna > orfs_list
				for file in `cat orfs_list`; do
					mv $file ${file/.f*/.fna}
				done
				rm -rf orfs_list
				cd $address/OUTPUT/$Out/read_hits
				ls *.fasta > hits_list
				for file in `cat hits_list`; do
					mv $file ${file/.f*/.hits.fasta}
				done
				rm -rf hits_list
				cd $address/OUTPUT/$Out
				d2=`date -u "+%s"`
				d3=$(echo "$d2 - $d1" |bc -l)
				echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds" > fulltime.tmp
				cat log.tmp1 fulltime.tmp > log.tmp
				cat Log.txt log.tmp > FL
				rm -rf Log.txt log.tmp log.tmp1 fulltime.tmp 
				mv FL Log.txt

			else
				echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference sequence."
				ppm2="0"
				cntg="0"
				ORFs="0"
			fi
			cd $address
		;;
	esac
	cd $address
	d2=`date -u "+%s"`
	d3=$(echo "$d2 - $d1" |bc -l)
	if [ -s cont_log* ]
	then
		cat cont_log* > int
		rm -rf cont_log*
		sed -i '/[a-z]/d' int
		cntg=`awk '{s+=$1} END {print s}' int`
		rm -rf int
	fi
	echo "$R1|$T1|$T2|$Ref|$SubRef|$d3|$reads|$buckets|$ppm1|$cntg" > par.time
	cat time.log par.time > time.log2; rm -rf time.log par.time; mv time.log2 time.log
	sed -i 's/|/\t/g' time.log
	echo -e "\n [BEAF12.01.17 worked in $R1 with reference as $Ref (output as $Out) for $d3 seconds] \n"
done < config.kp
rm -rf *.kp *.nmb 
rm -rf par.time
d99=`date -u "+%s"`
dtotal=$(echo "$d99 - $d0" |bc -l)
echo "BEAF1011.65 worked on files from config.file for $dtotal seconds."
