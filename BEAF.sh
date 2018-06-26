#!/usr/bin/env bash

address=$(cd "$(dirname "")" && pwd)/$(basename "")
ConfigFile=$(cd "$(dirname "config.file")" && pwd)/$(basename "config.file")

spades="/Path/To/SPADES/spades.py" ## spades="spades.py"
quast="/Path/To/QUAST/metaquast.py" ## "quast=metaquast.py"
usearch="/Path/To/usearch" ## usearch="/usr/bin/usearch"
fastqc="Path/To/fastqc" ## fastqc="fastqc"
cdhit="Path/To/cdhit" ## cdhit="cd-hit"
cutadapt="Path/To/cutadapt" ## cutadapt="cutadapt"
pyfasta="Path/To/pyfasta" ## pyfasta="pyfasta"
pigz="/Path/To/pigz" ## pigz="pigz"

LIB="/Path/To/Lib/" ## LIB="Lib"
orffinder="Path/To/bb.orffinder.pl" ## orffinder="$LIB/bb.orffinder.pl"
extpy="$LIB/ext.py"
SeqLength="$LIB/SeqLength.py"
PCAmaker="$LIB/PCA_maker.py"

BucketsFolder="Buckets"
ReferencesFolder="Reference_seqs" ## ReferencesFolder="Path/To/ReferencesFolder"

TrimFiles="True"
adapter="AGATCGGAAGAGC"

GenID=0.95
GenID1=0.85
GenID2=0.95
ProtID1=0.25
ProtID2=0.90
rDNA16SID1=0.90
rDNA16SID2=0.97

GENevalue="1e-20"
GENevalue1="1e-10"
GENevalue2="1e-20"
PROTevalue1="1e-5"
PROTevalue2="1e-5"
rDNAevalue1="1e-20"
rDNAevalue2="1e-20"

PROTqcovHits=50
PROTqcovORFs=50
rDNA16Sqcov=90

maxaccepts1=1
maxrejects1=32
maxaccepts2=5000
maxrejects2=5000
genmaxrejects2=64

rDNA16Scutoff=50
rDNAminCntgLength=500

strand="both"

GenomeSplitMethod="PyFasta"
GenomeCoverage=5
GenomeCoverage1=2.5
GenomeCoverage2=5
GenomeFragLength=100

threads=4
if [[ $(nproc) -gt 2 ]]
then
	threads=`echo "$(nproc) * 9 / 10" | bc` # If more than 2 threads in the computer, use 90% of threads
else
	if [[ $(nproc) -eq 2 ]] || [[ $(nproc) -eq 1 ]]
	then
		threads=1 # If only 1 or 2 threads in the computer, use only 1 thread
	else
		echo "Couldn't determine the number of threads in the computer. Please, specify threads manually." # If the system can't determine the number of threads as higher than 2, 2 nor 1, then ask for the threads to be determined manually, and, until it is, use the default of 4.
	fi
fi

CFLR="U"
ver="BEAF (full)"
KeepBlastDBs="no" #If you're sure you want to keep databases files from one loop to the next, change this to 'Y'. If you're sure you're not reusing blast databases, change to 'N'
KeepGenomeUDBs="no"
AlwaysKeepBuckets="no"
StopAfterMakingBuckets="no"

show_help ()
{
	echo "
	####	    BEAF - version 1.0.0 (Built on 22/05/2018)	     ####
	####	Avaiable at https://github.com/celiosantosjr/BEAF    ####

Usage: BEAF.sh [options]

Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

N=Integer (eg. 100)
n=Integer or float (eg. 90 or 0.90)
e=float or exponential (eg. 0.001 or 1e-3)

Tutorial options:
	-h, --help	Show this helpful help page
	--config_help	Show a help page in configuring your config.file

Basic options:
	--config [file]		Config file to be used [config.file]
	--output [folder]	Folder where output is sent [current folder]
	-t, --threads [N]	Number of threads [90% of avaiable threads]
	--reference [folder]	Folder where your reference files are located [Reference_seqs]

Advanced options for your run:
	--force_continue	BEAF continues a previous run
	--force_restart		BEAF deletes files from a previous run, starting anew.
	--disable_ref_folder	Use full paths for references in the config.file, disable the references folder (same as doing --reference "")
	--AlwaysKeepBuckets	Don't delete buckets after runs. Only use this if you are using the same sequencing file for all your runs.
	--DontTrim		Skip trimming of reads.
	--OnlyMakeBuckets	Stop after making buckets. Don't run the BEAF pipeline.
	--KeepBlastDBs		Don't delete blast databases after using them - for proteins and genes.
	--SedSplit		When making a genome database for genome searches, do not use PyFasta, simply cutting the genome at coverage 1x through UNIX commands

Homology search options
--Affecting all modules
	--maxaccepts [N]	Stop USEARCH search (2nd filter) after finding N results that match filter criteria [5000]
	--maxrejects [N]	Stop USEARCH search (apply to both filters) after finding N results that do not match filter criteria
	--maxrejects1 [N]	Stop USEARCH search (1st filter) after finding N results that do not match filter criteria [32]
	--maxrejects2 [N]	Stop USEARCH search (2nd filter) after finding N results that do not match filter criteria for protein/gene and taxonomy modes [5000]
	--genmaxrejects2 [N]	Stop USEARCH search (2nd filter) after finding N results that do not match filter criteria for genome mode [64]
--Genome
	--gen_id [n]		Minimum identity for genome analyses [0.95]
	--gen_id1 [n]		Minimum identity only for the 1st filter of genome analyses [0.85]
	--gen_id2 [n]		Minimum identity only for the 2nd filter of genome analyses (subrefs) [0.95]

	--gen_evalue1 [e]	Minimum e-value for the 1st filter of genome analyses
	--gen_evalue2 [e]	Minimum e-value for the 2nd filter of genome analyses

	--GenCoverage [n]	Genome coverage for the genome database used in the homology searches (if using subrefs, applies to both the 1st and 2nd filter)
	--GenCoverage1 [n]	Genome coverage for the genome database used in the 1st homology search (if using subrefs)
	--GenCoverage2 [n]	Genome coverage for the genome database used in the 2nd homology search (if using subrefs)
	--GenFragLength [N]	Length of fragments generated for the genome database used in the homology searches
	
--Protein/Gene
	--prot_id1 [n]		Minimum identity for the 1st filter of protein/gene analyses (general filter)
	--prot_id1 [n]		Minimum identity for the 2nd filter of protein/gene analyses (subref filter)

	--prot_evalue1 [e]	Minimum e-value for the 1st filter of protein/gene analyses
	--prot_evalue2 [e]	Minimum e-value for the 2nd filter of protein/gene analyses

	--prot_qcov1 [n]	Minimum query-coverage for the 1st filter of protein/gene analyses
	--prot_qcov2 [n]	Minimum query-coverage for the 2nd filter of protein/gene analyses

--Taxonomy
	--16S_id1 [N]		Minimum identity for the 1st filter of taxonomic analyses (general filter)
	--16S_id2 [N]		Minimum identity for the 2nd filter of taxonomic analyses (specific/OTU filter)

	--16S_evalue1 [e]	Minimum e-value for the 1st filter of taxonomic analyses
	--16S_evalue2 [e]	Minimum e-value for the 2nd filter of taxonomic analyses

	--16S_qcov [n]		Minimum query-coverage filter only for taxonomic analyses

	--16S_cutoff [N]	Minimum number of counts for an OTU to be accepted

	
"
}

show_config_help ()
{
	echo "STILL ON THE WORKS"
}

while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help|-help|\?|-\?)
			show_help
			exit
		;;
		--config_help|--help_config)
			show_config_help
			exit
		;;
		-s|-S|-Soft|-soft|-SOFT|--Soft|--soft|--SOFT|--Soft-Version|--soft-version|--Soft-version)
			ver="SOFT"
		;;
		--DontTrim|donttrim|Donttrim|dontTrim|DONTTRIM|DONTtrim|dontTRIM|noTrim|NoTrim|NOtrim|noTRIM|Notrim|NOTRIM)
			TrimFiles="False"
		;;
		--config|--Config|--CONFIG|-config|-Config|--CONFIG)
			if [[ -s ${2} ]]
			then
				ConfigFile=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
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
			else
				echo "Couldn't find the designated config.file '${1#*=}'"
				exit
			fi
		;;
		--output|--Output|--OUTPUT|-output|-Output|-OUTPUT|-O|-o|-Out|-OUT|-out|--out|--Out|--OUT|--o|--O|--address|--Address|--ADDRESS)
			if [[ -d "${2}" ]]
			then
				touch ${2}
			else
				mkdir ${2}
			fi
			address=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
			shift
		;;
		--output=*|--Output=*|--OUTPUT=*|-output=*|-Output=*|-OUTPUT=*|-O=*|-o=*|-Out=*|-OUT=*|-out=*|--out=*|--Out=*|--OUT=*|--o=*|--O=*|--address=*|--Address=*|--ADDRESS=*)
			if [[ -d "${1#*=}" ]]
			then
				touch ${1#*=}
			else
				mkdir ${1#*=}
			fi
			address=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
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
		--reference|--Reference|--REFERENCE|--references|--References|--REFERENCES|--ref|--Ref|--REF|--refs|--Refs|--REFS|-ref|-Ref|-REF|-reference|-references)
			if [[ -d "${2}" ]]
			then
				ReferencesFolder=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
			else
				ReferencesFolder="${2}"
			fi
			shift
		;;
		--reference=*|--Reference=*|--REFERENCE=*|--references=*|--References=*|--REFERENCES=*|--ref=*|--Ref=*|--REF=*|--refs=*|--Refs=*|--REFS=*|-ref=*|-Ref=*|-REF=*|-reference=*|-references=*)
			if [[ -d "${1#*=}" ]]
			then
				ReferencesFolder=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
			else
				ReferencesFolder="${1#*=}"
			fi
		;;
		--disable_ref_folder|--disable_ref)
			ReferencesFolder=""
		;;
		--gen_id)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							GenID1="1"
							GenID2="1"
						else
							GenID1="0.${2/./}"
							GenID2="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						GenID1="1"
						GenID2="1"
					else
						GenID1=0.${2#*.}
						GenID2=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--gen_id1)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							GenID1="1"
						else
							GenID1="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						GenID1="1"
					else
						GenID1=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--gen_id2)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							GenID2="1"
						else
							GenID2="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						GenID2="1"
					else
						GenID2=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--prot_id1)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							ProtID1="1"
						else
							ProtID1="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						ProtID1="1"
					else
						ProtID1=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--prot_id2)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							ProtID2="1"
						else
							ProtID2="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						ProtID2="1"
					else
						ProtID2=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--16S_id1)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							rDNA16SID1="1"
						else
							rDNA16SID1="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						rDNA16SID1="1"
					else
						rDNA16SID1=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--16S_id2)
			if [[ $(echo "${2} > 0" | bc -l) -eq 1 ]]
			then 
				if [[ $(echo "${2} > 1" | bc -l) -eq 1 ]]
				then
					if [[ $(echo "${2} > 100" | bc -l) -eq 1 ]]
					then
						echo "Invalid identity for genome filter. Input was over 100 (gen_id=${2}). Ignoring this input."
					else
						if [[ ${2} -eq 100 ]]
						then
							rDNA16SID2="1"
						else
							rDNA16SID2="0.${2/./}"
						fi
					fi
				else
					if [[ $(echo "${2} == 1" | bc -l) -eq 1 ]]
					then
						rDNA16SID2="1"
					else
						rDNA16SID2=0.${2#*.}
					fi
				fi
			else
				echo "Invalid identity for genome filter. Input was gen_id=${2}. Ignoring this input."
			fi
			shift
		;;
		--prot_evalue1)
			PROTevalue1=$2
			shift
		;;
		--prot_evalue2)
			PROTevalue2=$2
			shift
		;;
		--gen_evalue1)
			GENevalue1=$2
			shift
		;;
		--gen_evalue2)
			GENevalue2=$2
			shift
		;;
		--16S_evalue1)
			rDNAevalue1=$2
			shift
		;;
		--16S_evalue2)
			rDNAevalue2=$2
			shift
		;;
		--prot_qcov)
			if [[ "$2" -gt "1" ]]
			then
				PROTqcovORFs=$2
				PROTqcovHits=$2
			else
				if [[ "$2" -gt "0" ]]
				then
					PROTqcovORFs=$(echo "$2 * 100" | bc -l)
					PROTqcovHits=$(echo "$2 * 100" | bc -l)
				fi
			fi
			shift
		;;
		--prot_qcov1)
			if [[ "$2" -gt "1" ]]
			then
				PROTqcovHits=$2
			else
				if [[ "$2" -gt "0" ]]
				then
					PROTqcovHits=$(echo "$2 * 100" | bc -l)
				fi
			fi
			shift
		;;
		--prot_qcov2)
			if [[ "$2" -gt "1" ]]
			then
				PROTqcovORFs=$2
			else
				if [[ "$2" -gt "0" ]]
				then
					PROTqcovORFs=$(echo "$2 * 100" | bc -l)
				fi
			fi
			shift
		;;
		--16S_qcov)
			if [[ "$2" -gt "1" ]]
			then
				rDNA16Sqcov=$2
			else
				if [[ "$2" -gt "0" ]]
				then
					rDNA16Sqcov=$(echo "$2 * 100" | bc -l)
				fi
			fi
			shift
		;;
		--16S_cutoff)
			rDNA16Scutoff=$2
			shift
		;;
		--16S_min_cntg_length)
			rDNAminCntgLength=$2
			shift
		;;
#		--maxaccepts1)
#			maxaccepts1=$2
#			shift
#		;;
		--maxaccepts2|--maxaccepts)
			maxaccepts2=$2
			shift
		;;
		--maxrejects1)
			maxrejects1=$2
			shift
		;;
		--maxrejects2)
			maxrejects2=$2
			shift
		;;
		--allmaxrejects2)
			maxrejects2=$2
			genmaxrejects2=$2
			shift
		;;
		--genmaxrejects2)
			genmaxrejects2=$2
			shift
		;;
		--maxrejects)
			maxrejects1=$2
			maxrejects2=$2
			shift
		;;
		--allmaxrejects)
			genmaxrejects2=$2
			maxrejects1=$2
			maxrejects2=$2
			shift
		;;
		# --PyFastaSplit|--pyfastasplit|--pfsplit|--pfs)
		# 	GenomeSplitMethod="PyFasta"
		# ;;
		--SedSplit|--sedsplit|--SEDSPLIT|--Sedsplit|--SEDsplit|--SEDSplit|--sedSPLIT|--sedSplit)
			GenomeSplitMethod="SED"
			GenomeCoverage=1
		;;
		--GenomeSplitMethod)
			case $2 in
				SED|Sed|sed)
					GenomeSplitMethod="SED"
				;;
				PF|PyFasta|PYFASTA|Pyfasta|PyFASTA|PYfasta|PYFasta|pyFASTA|pyFasta|pf|Pf)
					GenomeSplitMethod="PyFasta"
				;;
				*)
					echo "Couldn't recognize genome split method. Using PyFasta."
					GenomeSplitMethod="PyFasta"
				;;
			esac
		;;
		--GenCoverage)
			GenomeCoverage=$2
			GenomeCoverage1=$2
			GenomeCoverage2=$2
			shift
		;;
		--GenCoverage1)
			GenomeCoverage1=$2
			shift
		;;
		--GenCoverage2)
			GenomeCoverage2=$2
			shift
		;;
		--GenFragLength)
			GenomeFragLength=$2
			shift
		;; 
		--KeepBlastDBs|--keepblastdbs|--keepbdb|--keepbdbs|--Keepbdb|--Keepbdbs)
			KeepBlastDBs="yes"
		;;
		--KeepGenomeUDBs|--keepgenomeudbs|--keepudb|--keepudbs|--Keepudb|--Keepudbs)
			KeepGenomeUDBs="yes"
		;;
		--OnlyMakeBuckets)
			StopAfterMakingBuckets="Y"
		;;
		--AlwaysKeepBuckets|--KeepBuckets|--AKB)
			AlwaysKeepBuckets="Y"
		;;
		*)
			echo "Couldn't recognize command '${1}'. Ignoring it."
			sleep 1
		;;
	esac
	shift
done

echo "Starting program in folder ${address}, at $(date "+%X, %e/%m/%Y")"


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
	case $AlwaysKeepBuckets in 
		Y|y|Yes|yes|YES)
			echo "Y" >> $address/Keep_config.kp
		;;
		*)
			echo "N" >> $address/Keep_config.kp
		;;
	esac
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
			G|g)
			;;
			P|p|N|n)
			;;
			16S|16s|16)
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
			*.fasta|*.fa|*.fna|*.fsa|*.FASTA|*.FA|*.FNA|*.FSA|*.Fasta|*.Fa|*.Fsa|*.Fna)
				case $T2 in
					F|f)
					;;
					*)
						echo "# You're using a wrong config.file. You're using a fasta file ($R1), which means you must use the fasta option (F) instead of $T2."
						rm -rf $address/config.kp $address/config.check
						exit
					;;
				esac
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
	if [[ ! -s $address/$(basename $ConfigFile) ]]
	then
		cp ${ConfigFile} ${address}/$(basename "${ConfigFile}")
		if [[ ! -s $address/$(basename $ConfigFile) ]]
		then
			ConfigFile="${address}/$(basename ${ConfigFile})"
		fi
	fi
	echo "2" > $address/CR.step; CFLR="N"
fi
}

TimeHeader () # Generates the header for the Log.tsv file (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	echo "______|________|_____|_____|_________|______|____|_____|_______|____|_______|___________|_____________|__________|___________|____|__________|____________|_________|__________
Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs|AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntgSize|ORFs|AvgSizeORF|TotalSizeORF|StdDevORF|MaxORFSize" >> $address/Log.tsv
	echo "3" > $address/CR.step; CFLR="N"
fi
}

SoftTimeHeader () # Generates the header for the Log.tsv file (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	echo "______|________|_____|_____|_________|______|____|_____|_______|____|_______|___________|_____________|__________|___________|____|__________|____________|_________|__________
Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs|AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntgSize" >> $address/Log.tsv
	echo "3" > $address/CR.step; CFLR="N"
fi
}

Trim () # Trims adapters from Illumina data and merges sequences R1 and R2 into one file, when using pair-end (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "3" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping reads trimming opperations"; else
	rm -rf $address/${BucketsFolder}
	BucketsFolder="Bucket_$(basename $R1 | sed 's/.gz//' | sed 's/\..*//')"
	mkdir $address/${BucketsFolder}
	date -u +%s > $address/${BucketsFolder}/datestarttrimming.tmp
	case $TrimFiles in
		False|false|FALSE)
			case $T2 in
				R|r|I|i)
					cp $R1 $address/${BucketsFolder}/
					cp $R2 $address/${BucketsFolder}/
					${pigz} -p ${threads} -d -c $address/${BucketsFolder}/*.gz | ${pigz} -p ${threads} -c > $address/${BucketsFolder}/FastaQ-full.gz
				;;
				F|f)
				;;
			esac
		;;
		*)
			case $T2 in 
				R|r)
					case $T1 in
						16S|16s|16)
							echo "# We're now trimming your files. Only R1 will be used, as interleaved files may result in wrong abundance values..."
							${cutadapt} --cores=${threads} --quiet --minimum-length 80 --max-n 0.008 --trim-n -a $adapter -e 0.1 -O 5 -q 24,24 -o $address/${BucketsFolder}/FastaQ-full.gz $R1 > $address/${BucketsFolder}/cutadaptlog.txt # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
						;;
						*)
							echo "# Trimming and merging your files..."
							${cutadapt} --cores=${threads} --quiet --interleaved --minimum-length 80 --max-n 0.008 --trim-n -a $adapter -A $adapter -e 0.1 -O 5 -q 24,24 -o $address/${BucketsFolder}/FastaQ-full.gz $R1 $R2 > $address/${BucketsFolder}/cutadaptlog.txt # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
						;;
					esac
				;;
				I|i)
					echo "# We are now trimming your files..."
					${cutadapt} --cores=${threads} --quiet --minimum-length 80 --max-n 0.008 --trim-n -a $adapter -e 0.1 -O 5 -q 24,24 -o $address/${BucketsFolder}/FastaQ-full.gz $R1 > $address/${BucketsFolder}/cutadaptlog.txt # Parameters of reads trimming should be specified here. Trimming universal Illumina adapter
				;;
				F|f)
					case $R1 in
						*.gz)
							${pigz} -p $threads -d -c $R1 > $address/${BucketsFolder}/FastaQ-full.fa
						;;
						*.fasta|*.fa|*.fna|*.fsa|*.FASTA|*.FA|*.FNA|*.FSA|*.Fasta|*.Fa|*.Fsa|*.Fna)
							cp $R1 $address/${BucketsFolder}/FastaQ-full.fa
						;;
					esac
				;;
			esac
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/${BucketsFolder}/datestarttrimming.tmp)" | bc -l > $address/${BucketsFolder}/trimmingtime.nmb
	rm -rf $address/${BucketsFolder}/cutadaptlog.txt
	echo "4" > $address/CR.step; CFLR="N"
fi
}

QAnConversion () # Makes assessment of trimmage and converts file from fastq to fasta (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "4" ]]; then echo "******Skipping FASTQ handling / FASTQ conversion"; else
	date -u +%s > $address/${BucketsFolder}/datestartfastqc.tmp
	case $T2 in
		R|r|I|i)
			echo "# Starting quality assessment of trimming..."
			rm -rf $address/${BucketsFolder}/FASTQCresults; rm -rf $address/${BucketsFolder}/FastaQ-full.fa
			mkdir $address/${BucketsFolder}/FASTQCresults
			${fastqc} --quiet --threads $threads -f fastq -o $address/${BucketsFolder}/FASTQCresults $address/${BucketsFolder}/FastaQ-full.gz > $address/${BucketsFolder}/fastqclog.txt # FASTQ assessment is done here
			rm -rf $address/${BucketsFolder}/FASTQCresults/*_fastqc ; rm -rf $address/${BucketsFolder}/fastqclog.txt
			echo "# We will convert merged file to fasta format."
			${pigz} -p ${threads} -d -c $address/${BucketsFolder}/FastaQ-full.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/${BucketsFolder}/FastaQ-full.fa
			echo "# Removing unwanted files..."
			rm -rf $address/${BucketsFolder}/*.gz
		;;
		F|f)
			rm -rf $address/${BucketsFolder}/FASTQCresults
			echo "# FASTA file type identified. Since FASTA does not have PHRED values, It wont be assessed by FASTQC algorithm."
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/${BucketsFolder}/datestartfastqc.tmp)" | bc -l > $address/${BucketsFolder}/fastqctime.nmb
	echo "5" > $address/CR.step; CFLR="N"
fi
}

CopyFile () # Soft - Copies Illumina data to a separate folder to work on it (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "3" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping files copying process"; else
	rm -rf $address/${BucketsFolder}
	BucketsFolder="Bucket_$(basename $R1 | sed 's/.gz//' | sed 's/\..*//')"
	mkdir $address/${BucketsFolder}
	case $T2 in 
		R|r)
			echo "# Copying file 1 from storage..."
			cp -r $R1 $address/${BucketsFolder}
			echo "# Copying file 2 from storage..."
			cp -r $R2 $address/${BucketsFolder}
		;;
		I|i|F|f)
			echo "# Copying file from storage..."
			cp -r $R1 $address/${BucketsFolder}
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
			rm -rf $address/${BucketsFolder}/FastaQ-full.gz $address/${BucketsFolder}/FastaQ-full.fa
			echo "Merging files"
			${pigz} -p ${threads} -d -c $address/${BucketsFolder}/*.gz | ${pigz} -p ${threads} -c > $address/${BucketsFolder}/FastaQ-full.gz
			echo "# We will convert merged file to fasta format."
			${pigz} -p ${threads} -d -c $address/${BucketsFolder}/FastaQ-full.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/${BucketsFolder}/FastaQ-full.fa
			rm -rf $address/${BucketsFolder}/*.gz
		;;
		I|i)
			rm -rf $address/${BucketsFolder}/FastaQ-full.gz $address/${BucketsFolder}/FastaQ-full.fa
			${pigz} -p ${threads} -d -c `ls $address/${BucketsFolder}/*.gz` | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $address/${BucketsFolder}/FastaQ-full.fa
		;;
		F|f)
			rm -rf $address/${BucketsFolder}/FastaQ-full.fa
			${pigz} -p ${threads} -d -c `ls $address/${BucketsFolder}/*.gz`> $address/${BucketsFolder}/FastaQ-full.fa
			rm -rf $address/${BucketsFolder}/*.gz
		;;
	esac
	echo "5" > $address/CR.step; CFLR="N"
fi
}

BucketEngine () # Breaks the data into separate buckets to speed up the process 
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "5" && `cat $address/CR.step` != "24" ]]; then echo "******Skipping Buckets System"; else
	if [[ -s $address/${BucketsFolder}/bk_list.txt ]]
	then
		reads=`cat $address/${BucketsFolder}/reads.nmb`
		# rm -rf $address/${BucketsFolder}/*.bk
	else
		if [[ -s $address/${BucketsFolder}/buckets_list.txt ]];
		then
			buckets=`ls $address/${BucketsFolder}/*.bk | wc -l`
			reads=`cat $address/${BucketsFolder}/reads.nmb`
			echo "$buckets buckets were found from the previous run. We'll be using those for this run as well."
			echo "(The $buckets buckets had already been generated in a previous run. This is the time it took in that run for them to be generated, not in the current run)" > $address/${BucketsFolder}/bucketpreviouslygeneratedmessage.tmp
		else
			if [[ -s $address/${BucketsFolder}/reads.nmb ]]
			then
				touch $address/${BucketsFolder}/reads.nmb
			else
				grep ">" $address/${BucketsFolder}/FastaQ-full.fa | wc -l > $address/${BucketsFolder}/reads.nmb
				reads=`cat $address/${BucketsFolder}/reads.nmb`
			fi
			original_size=$(wc -c $address/${BucketsFolder}/FastaQ-full.fa | sed 's/ .*//') # ; echo "original size $original_size"
			# sizeinKb=`expr $original_size / 1024`
			# sizeinMb=`expr $sizeinKb / 1024` # `expr $original_size / 1048576`
			buckets=`expr $original_size / 268435456 + 1` # ; echo "buckets $buckets" # (256Mb buckets)
			bucketsize=`expr $original_size / $buckets + 1024` # ; echo "bucketsize $bucketsize" # (+1Kb)
			date -u +%s > $address/${BucketsFolder}/datestartbucketengine.tmp
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
				mv $address/${BucketsFolder}/FastaQ-full.fa $address/${BucketsFolder}/1.bk
			else
				echo "# We identified $reads reads (file size is $original_size bytes, or $(expr $original_size / 1048576)Mb). It will take $buckets bucket steps (the size of each bucket will be $bucketsize bytes, or $(expr $bucketsize / 1048576)Mb)."
				echo "### Starting operation of cutting and readapting..."
				echo "## Generating buckets..."
				cat $address/${BucketsFolder}/FastaQ-full.fa | parallel --pipe --block $bucketsize --recstart ">" "cat >$address/${BucketsFolder}/{#}.bk"
				ls $address/${BucketsFolder}/*.bk | sort -S 50% --parallel=${threads} -V > $address/${BucketsFolder}/bk_testlist.txt
				try=0
				while [[ "$(cat $address/${BucketsFolder}/bk_testlist.txt | wc -l)" -gt "$buckets" && "$(echo $(wc -c $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt) | sed 's/ .*//') + $(wc -c $(tail -n 2 $address/${BucketsFolder}/bk_testlist.txt | head -n 1) | sed 's/ .*//') | bc)" -lt "$(wc -c $address/${BucketsFolder}/1.bk | sed 's/ .*//')" ]] ; do
					cat $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt) >> $(tail -n 2 $address/${BucketsFolder}/bk_testlist.txt | head -n 1)
					rm -rf $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt)
					ls $address/${BucketsFolder}/*.bk | sort -S 50% --parallel=${threads} -V > $address/${BucketsFolder}/bk_testlist.txt
				done
				while [[ $(cat $address/${BucketsFolder}/bk_testlist.txt | wc -l) -gt $buckets && ${try} -le 5 ]] ; do
					if [[ -s "$address/${BucketsFolder}/$(echo $buckets + 1 | bc).bk" ]]
					then
						if [[ "$(echo $(wc -c $address/${BucketsFolder}/$buckets.bk | sed 's/ .*//') + $(wc -c $address/${BucketsFolder}/$(echo "$buckets + 1" | bc).bk | sed 's/ .*//') | bc)" -lt "$(wc -c $address/${BucketsFolder}/1.bk | sed 's/ .*//')" ]]
						then
							echo "Fixing buckets (try $try)"
							cat $address/${BucketsFolder}/$(echo "$buckets + 1" | bc).bk >> $address/${BucketsFolder}/$buckets.bk
							rm -rf $address/${BucketsFolder}/$(echo "$buckets + 1" | bc).bk
						fi
					fi
					ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/bk_testlist.txt
					if [[ "$(cat $address/${BucketsFolder}/bk_testlist.txt | wc -l)" -gt "$buckets" ]]
					then
						rm -rf $address/${BucketsFolder}/*.bk
						rm -rf $address/${BucketsFolder}/bk_list.txt $address/${BucketsFolder}/bk_testlist.txt
						original_size=$(wc -c $address/${BucketsFolder}/FastaQ-full.fa | sed 's/ .*//') # ; echo "original size of file: $original_size"
						buckets=`expr $original_size / 268435456 + 1` # ; echo "buckets: $buckets (256Mb buckets)"
						bucketsize=`expr (($original_size / $buckets) + (1048576 \* $try))` # ; echo "bucketsize: $bucketsize (+1Mb than previous try)"
						cat $address/${BucketsFolder}/FastaQ-full.fa | parallel --pipe --block $bucketsize --recstart ">" "cat >$address/${BucketsFolder}/{#}.bk"
						ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/bk_testlist.txt
					fi
					ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/bk_testlist.txt
					((try++))
				done
				ls $address/${BucketsFolder}/*.bk | sort -S 50% --parallel=${threads} -V > $address/${BucketsFolder}/bk_testlist.txt
				while [[ "$(cat $address/${BucketsFolder}/bk_testlist.txt | wc -l)" -gt "$buckets" && "$(echo $(wc -c $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt) | sed 's/ .*//') + $(wc -c $(tail -n 2 $address/${BucketsFolder}/bk_testlist.txt | head -n 1) | sed 's/ .*//') | bc)" -lt "$(wc -c $address/${BucketsFolder}/1.bk | sed 's/ .*//')" ]] ; do
					cat $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt) >> $(tail -n 2 $address/${BucketsFolder}/bk_testlist.txt | head -n 1)
					rm -rf $(tail -n 1 $address/${BucketsFolder}/bk_testlist.txt)
					ls $address/${BucketsFolder}/*.bk | sort -S 50% --parallel=${threads} -V > $address/${BucketsFolder}/bk_testlist.txt
				done
				rm -rf $address/${BucketsFolder}/bk_list.txt
				mv $address/${BucketsFolder}/bk_testlist.txt $address/${BucketsFolder}/bk_list.txt
			fi
			echo "$(date -u +%s) - $(cat $address/${BucketsFolder}/datestartbucketengine.tmp)" | bc -l > $address/${BucketsFolder}/bucketenginetime.nmb
		fi
	fi
	echo "### Removing temporary files (stage 1)..."
	rm -rf $address/${BucketsFolder}/*.txt; rm -rf $address/${BucketsFolder}/*.gz # ; rm -rf $address/${BucketsFolder}/*.fa
	ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/buckets_list.txt
	buckets="$(cat $address/${BucketsFolder}/buckets_list.txt | wc -l)"
	cat $address/${BucketsFolder}/buckets_list.txt | sed "s@$address/${BucketsFolder}/@@" | sort -S 50% --parallel=${threads} -V > $address/${BucketsFolder}/buckets_search.txt
	rm -rf $address/${BucketsFolder}/bk_list.txt $address/${BucketsFolder}/bk_testlist.txt
	buckets=`ls $address/${BucketsFolder}/*.bk | wc -l`
	if [[ $buckets != 1 ]]
	then
		echo "A total of ${buckets} buckets were generated, with an average size of $(echo $(wc -c $address/${BucketsFolder}/*.bk | tail -n 1 | awk '{ print $1 }') / $buckets / 1048576 | bc)Mb."
	fi
	BucketTotalTime=$(echo "$(date -u +%s) - $d0" |bc -l)
	if [[ "${BucketTotalTime}" -ge "$(echo $(cat $address/${BucketsFolder}/trimmingtime.nmb) + $(cat $address/${BucketsFolder}/fastqctime.nmb) + $(cat $address/${BucketsFolder}/bucketenginetime.nmb) | bc -l)" 	]]
	then
		echo "$BucketTotalTime" > $address/${BucketsFolder}/BucketTotalTime.nmb
	else
		echo $(cat $address/${BucketsFolder}/trimmingtime.nmb) + $(cat $address/${BucketsFolder}/fastqctime.nmb) + $(cat $address/${BucketsFolder}/bucketenginetime.nmb) | bc -l > $address/${BucketsFolder}/BucketTotalTime.nmb
		BucketTotalTime=$(cat $address/${BucketsFolder}/BucketTotalTime.nmb)
	fi
	case $StopAfterMakingBuckets in
		Y|y|Yes|yes|YES)
			rm -rf $address/CR.step $address/CR.mode $address/config.kp
			d99=`date -u "+%s"`
			dtotal=$(echo "$d99 - $d0" |bc -l)
			((totaldays=${dtotal}/86400))
			((totalhours=(${dtotal}%86400)/3600))
			((totalminutes=((${dtotal}%86400)%3600)/60))
			((totalseconds=((${dtotal}%86400)%3600)%60))
			if [[ $totaldays -ge 1 ]]
			then
				echo "BEAF OnlyMakeBuckets worked for $dtotal seconds (${totaldays}d${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date +%X)."
			else
				if [[ $totalhours -ge 1 ]]
				then
					echo "BEAF1011.65 worked for $dtotal seconds (${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date +%X)."
				else
					if [[ $totalminutes -ge 1 ]]
					then
						echo "BEAF OnlyMakeBuckets worked for $dtotal seconds (${totalminutes}m${totalseconds}s), ending at $(date +%X)."
					else
						echo "BEAF OnlyMakeBuckets worked for $dtotal seconds (${totalseconds}s), ending at $(date +%X)."
					fi
				fi
			fi
			exit
		;;
		*)
		;;
	esac
	echo "6" > $address/CR.step; CFLR="N"
fi
}

Filter1 () # Creates udb files from reference if needed, and then align each bucket against reference file to filter for homology.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "6" ]]; then echo "******Skipping Filtering System 1"; else
	GenomeCoverage=${GenomeCoverage2}
	rm -rf $address/${BucketsFolder}/FastaQ-full.fa
	date -u +%s > $address/OUTPUT/${Out}/datestartfilter1.tmp
	case $T1 in
		G|g)
			if [[ -s $address/${BucketsFolder}/ReferenceGenome.udb ]]
			then
				echo "# Continuing searches..."
			else
				if [[ -s $Ref ]]
				then
					GenomeCoverage=${GenomeCoverage2}
					case $Ref in 
						*.udb|*.UDB|*.uDB|*.Udb)
							echo "### Using provided USEARCH database: $Ref"
							cp $Ref $address/${BucketsFolder}/ReferenceGenome.udb
						;;
						*)
							echo "### Making database from reference genome $Ref..."
							cat $Ref | awk NF > $address/${BucketsFolder}/RefGenome.fasta
							case $GenomeSplitMethod in
								PyFasta|Pyfasta|pyfasta|PYFASTA|pyFasta|PyFASTA|pyFASTA|PYfasta|PYFasta|PF|pf|Pf|pF)
									Overlap1=$(echo "$GenomeFragLength - $(echo "$GenomeFragLength / ${GenomeCoverage}" | bc)" | bc)
									echo "##### Splitting sequences with PyFasta Split using an overlap of ${Overlap1} bases in order to provide a coverage of ${GenomeCoverage}x."
									${pyfasta} split -n 1 -k ${GenomeFragLength} -o ${Overlap1} $address/${BucketsFolder}/RefGenome.fasta > $address/${BucketsFolder}/pyfastalog.txt
									mv $address/${BucketsFolder}/RefGenome.split.*.fasta $address/${BucketsFolder}/RefGen_compressed.fa
									rm -rf $address/${BucketsFolder}/RefGenome.fasta.flat $address/${BucketsFolder}/RefGenome.fasta.gdx $address/${BucketsFolder}/RefGenome.fasta $address/${BucketsFolder}/pyfastalog.txt
									echo "New udb" > $address/${BucketsFolder}/newudb.txt
								;;
								*)	
									echo "##### Splitting genome into $GenomeFragLength basepair sequences - no overlaps, 1x coverage"
									cat $address/${BucketsFolder}/RefGenome.fasta | sed '/>/d' | tr -d '\n' | sed "s/.\{${GenomeFragLength}\}/&\n>\n/g" | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $address/${BucketsFolder}/RefGen_compressed.fa
									rm -rf $address/${BucketsFolder}/RefGenome.fasta
									if [[ $(tail -n 1 $address/${BucketsFolder}/RefGen_compressed.fa | grep ">" | wc -l) -gt 0 ]]
									then
										sed -i '$ d' $address/${BucketsFolder}/RefGen_compressed.fa
									fi					
								;;
							esac
							rm -rf $address/${BucketsFolder}/none $address/${BucketsFolder}/RefGenome.fasta
						;;
					esac
				else
					if [[ -d $ReferencesFolder/$SubRef ]]
					then
						echo "BEAF couldn't find your Reference Genome '$Ref', so it will create a Reference Genome based on the sum of your SubReferences."
						GenomeCoverage=${GenomeCoverage1}
						rm -rf $ReferencesFolder/$SubRef/SubRef_fasta.list ; rm -rf $ReferencesFolder/$SubRef/Pooled_SubRefs.fasta
						ls $ReferencesFolder/$SubRef | grep -v "^SubRef_.*_cov.*.fasta$" | grep -vw "Pooled_SubRefs.fasta" | grep -Ei "(.fasta|.fa|.faa|.fas|.fna|.fsa|.ffn|.frn|.mpfa)$" | sed "s#$ReferencesFolder/$SubRef/##" >> $ReferencesFolder/$SubRef/SubRef_fasta.list
						case $GenomeSplitMethod in
							PyFasta|Pyfasta|pyfasta|PYFASTA|pyFasta|PyFASTA|pyFASTA|PYfasta|PYFasta|PF|pf|Pf|pF)
								Overlap1=$(echo "$GenomeFragLength - $(echo "$GenomeFragLength / ${GenomeCoverage}" | bc)" | bc)
								echo "##### Splitting sequences with PyFasta Split using an overlap of ${Overlap1} bases in order to provide a coverage of ${GenomeCoverage}x."
								for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
									cp $ReferencesFolder/$SubRef/$sub $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
									${pyfasta} split -n 1 -k ${GenomeFragLength} -o ${Overlap1} $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta > $address/${BucketsFolder}/pyfastalog.txt
									cat $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.split.*.fasta > $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}_compressed.fa
									if [[ "$(echo "0${GenomeCoverage}" == "0${GenomeCoverage}" | bc)" -eq "1" ]]
									then
										cp $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}_compressed.fa $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
									fi
									rm -rf $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.split.*.fasta ; rm -rf $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta.flat $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta.gdx $address/${BucketsFolder}/pyfastalog.txt $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
									cat $ReferencesFolder/$SubRef/$sub >> $ReferencesFolder/$SubRef/Pooled_SubRefs.fasta
								done
								cat $address/${BucketsFolder}/SubRef_*_cov${GenomeCoverage}_compressed.fa > $address/${BucketsFolder}/RefGen_compressed.fa
								rm -rf $address/${BucketsFolder}/SubRef_${sub%.f*}_cov${GenomeCoverage}_compressed.fa
								echo "New udb" > $address/${BucketsFolder}/newudb.txt
								Ref="$ReferencesFolder/$SubRef/Pooled_SubRefs.fasta"
								echo "Pooled your SubReferences into a new Reference file at $Ref"
							;;
							*)	
								echo "##### Splitting genome into $GenomeFragLength basepair sequences - no overlaps, 1x coverage"
								for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
									cat $ReferencesFolder/$SubRef/$sub >> $ReferencesFolder/$SubRef/Pooled_SubRefs.fasta
								done
								cat $ReferencesFolder/$SubRef/Pooled_SubRefs.fasta | awk NF | sed '/>/d' | tr -d '\n' | sed "s/.\{${GenomeFragLength}\}/&\n>\n/g" | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $address/${BucketsFolder}/RefGen_compressed.fa
								if [[ $(tail -n 1 $address/${BucketsFolder}/RefGen_compressed.fa | grep ">" | wc -l) -gt 0 ]]
								then
									sed -i '$ d' $address/${BucketsFolder}/RefGen_compressed.fa
								fi
								Ref="$SubRef/Pooled_SubRefs.fasta"
								echo "Pooled your SubReferences into a new Reference file at $Ref"
							;;
						esac
					else
						echo "Couldn't find your Reference Genome '$Ref'"
					fi	
				fi
				echo "##### Making a USEARCH database."
				# ${cdhit} -i $address/${BucketsFolder}/RefGen_compressed.fa -o $address/${BucketsFolder}/RefGen_nr.fa -c 0.98 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5 > $address/${BucketsFolder}/cdhitlog # Parameters for CD-Hit (removing redundancy) should be specified here
				# rm -rf $address/${BucketsFolder}/cdhitlog
				${usearch} -makeudb_usearch $address/${BucketsFolder}/RefGen_compressed.fa -output $address/${BucketsFolder}/ReferenceGenome.udb -quiet > $address/${BucketsFolder}/makeudblog.txt # udb from reference is created here for genomes
				rm -rf $address/${BucketsFolder}/makeudblog.txt
				rm -rf $address/${BucketsFolder}/RefGen_compressed.fa
				rm -rf $address/OUTPUT/${Out}/*.m7
				echo "# Starting searches..."
			fi
			if [[ -d $ReferencesFolder/$SubRef ]]
			then
				GenID=$GenID1
				GENevalue=$GENevalue1
			else
				GenID=$GenID2
				GENevalue=$GENevalue2
			fi
			case $Keep in
				Y|y|yes|Yes|YES)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						if [[ -s $address/OUTPUT/${Out}/$buck.m8 ]]
						then
							touch $address/OUTPUT/${Out}/$buck.m8
						else
							echo -en "\r"; echo -e "# Searching against reference $Ref (id $GenID, evalue ${GENevalue} and maxrejects ${maxrejects1}) and keeping buckets... ${buck%.bk}/$buckets"
							${usearch} -usearch_global $address/${BucketsFolder}/$buck -db $address/${BucketsFolder}/ReferenceGenome.udb -strand both -id $GenID -evalue ${GENevalue} --maxaccepts 1 --maxrejects ${maxrejects1} -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm for genome binning should be specified here
							rm -rf $address/${BucketsFolder}/usearch.tmp
							mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
							sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
						fi
					done
				;;
				*)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref (id $GenID, evalue ${GENevalue} and maxrejects ${maxrejects1}) and deleting buckets... ${buck%.bk}/$buckets"
						${usearch} -usearch_global $address/${BucketsFolder}/$buck -db $address/${BucketsFolder}/ReferenceGenome.udb -strand both -id $GenID -evalue ${GENevalue} --,axaccepts 1 --maxrejects ${maxrejects1} -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm for genome binning should be specified here
						rm -rf $address/${BucketsFolder}/usearch.tmp
						mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
						rm -rf $address/${BucketsFolder}/$buck
						touch $address/${BucketsFolder}/$buck
						sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
					done
					rm -rf $address/${BucketsFolder}/*.bk
				;;
			esac
			case $KeepGenomeUDBs in
				Y|y|YES|Yes|yes)
					if [[ -s $ReferencesFolder/BEAF_GenomeUDBs_Folder/$(basename $Ref).udb ]]
					then
						touch $address/${BucketsFolder}/ReferenceGenome.udb
					else
						if [[ -s $address/${BucketsFolder}/newudb.txt ]]
						then
							if [[ -d $ReferencesFolder/BEAF_GenomeUDBs_Folder ]]
							then
								touch $ReferencesFolder/BEAF_GenomeUDBs_Folder
							else
								mkdir $ReferencesFolder/BEAF_GenomeUDBs_Folder
							fi
							mv $address/${BucketsFolder}/ReferenceGenome.udb $ReferencesFolder/BEAF_GenomeUDBs_Folder/$(basename $Ref).udb
						fi
					fi
				;;
				*)
					touch $address/${BucketsFolder}/ReferenceGenome.udb
				;;
			esac
			rm -rf $address/${BucketsFolder}/*.udb; rm -rf $address/OUTPUT/${Out}/*.m7; rm -rf $address/${BucketsFolder}/buckets_search.txt $address/${BucketsFolder}/usearch.tmp $address/${BucketsFolder}/RefGenome.fasta $address/${BucketsFolder}/ReferenceGenome.udb
		;;
		P|p|N|n)
			echo "# Starting searches..."
			if [[ -s $address/${BucketsFolder}/udblist ]]
			then
				touch $address/${BucketsFolder}/udblist
			else
				rm -rf $address/${BucketsFolder}/*.udb
				case $Ref in
					*.fa|*.fasta|*.fas|*.faa|*.fna|*.fsa|*.FA|*.FASTA|*.FAS|*.FAA|*.FNA|*.FSA) # Tests if reference is in fasta format
						echo "# Recognized $Ref file as fasta format. Making udb..."
						cat $Ref | awk NF > $address/${BucketsFolder}/none
						sed -i '/>/d' $address/${BucketsFolder}/none
						cat $address/${BucketsFolder}/none | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $address/${BucketsFolder}/md8
						if [[ $(tail -n 1 $address/${BucketsFolder}/md8 | grep ">" | wc -l) -gt 0 ]]
						then
							sed -i '$ d' $address/${BucketsFolder}/md8
						fi
						${usearch} -makeudb_usearch $address/${BucketsFolder}/md8 -output $address/${BucketsFolder}/$Ref.udb -quiet > $address/${BucketsFolder}/makeudblog.txt
						rm -rf $address/${BucketsFolder}/makeudblog.txt
						rm -rf $address/${BucketsFolder}/md8
						ls $address/${BucketsFolder}/*.udb > $address/${BucketsFolder}/udblist
					;;
					*.udb) # In case reference is not in faste format, tests if it is an udb file
						echo "# Recognized $Ref file as .udb format."
						cp $Ref $address/${BucketsFolder}
						ls $address/${BucketsFolder}/*.udb > $address/${BucketsFolder}/udblist
					;;
					*) 
						echo "Couldn't recognize $Ref file as neither fasta nor udb format. Will try to use it as udb regardless."
						cp $Ref $address/${BucketsFolder}/$Ref.trying.udb
						ls $address/${BucketsFolder}/$Ref.trying.udb > $address/${BucketsFolder}/udblist
					;;
				esac
			fi
			dbinuse=`cat $address/${BucketsFolder}/udblist` # Either a file originally in udb or a fasta format converted to udb will be used here
			rm -rf $address/OUTPUT/${Out}/*.m7
			case $Keep in
				Y|y)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref (id $ProtID1, evalue ${PROTevalue1} and maxrejects ${maxrejects1}) and keeping buckets... ${buck%.bk}/$buckets"
						${usearch} -usearch_local $address/${BucketsFolder}/$buck -db $dbinuse -strand both -id $ProtID1 -evalue ${PROTevalue1} --maxrejects ${maxrejects1} -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm for protein/gene binning should be specified here
						rm -rf $address/${BucketsFolder}/usearch.tmp
						mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
						sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
					done
				;;
				*)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref (id $ProtID1, evalue ${PROTevalue1} and maxrejects ${maxrejects1}) and deleting buckets... ${buck%.bk}/$buckets"
						${usearch} -usearch_local $address/${BucketsFolder}/$buck -db $dbinuse -strand both -id $ProtID1 -evalue ${PROTevalue1} --maxrejects ${maxrejects1} -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm for protein/gene binning should be specified here
						rm -rf $address/${BucketsFolder}/usearch.tmp
						mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
						rm -rf $address/${BucketsFolder}/$buck
						touch $address/${BucketsFolder}/$buck
						sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
					done
					rm -rf $address/${BucketsFolder}/*.bk
				;;
			esac
			rm -rf $address/${BucketsFolder}/*.udb; rm -rf $address/${BucketsFolder}/udblist; rm -rf $address/OUTPUT/${Out}/*.m7; rm -rf $address/${BucketsFolder}/buckets_search.txt
		;;
		16S|16s|16) 
			echo "# Starting searches..."
			rm -rf $address/OUTPUT/${Out}/*.m7
			if [[ ! -s $Ref ]]
			then
				if [[ -s $ReferencesFolder/$Ref ]]
				then
					Ref="$ReferencesFolder/$Ref"
				else
					if [[ -s $ReferencesFolder/nr97_Greengenes.udb ]]
					then
						echo "Using a clustered Green Genes database (clustered to 97% identity) to search for OTUs."
						Ref=nr97_Greengenes.udb
					else
						echo "Couldn't find reference $Ref in References Folder $ReferencesFolder."
					fi
				fi
			fi
			case $Keep in
				Y|y)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						if [[ -s $address/OUTPUT/${Out}/$buck.m8 ]]
						then
							touch $address/OUTPUT/${Out}/$buck.m8
						else
							echo -en "\r"; echo -e "# Searching against reference $Ref (id $rDNA16SID1, evalue ${rDNAevalue1} and maxrejects ${maxrejects1}) and keeping buckets... ${buck%.bk}/$buckets"
							${usearch} -usearch_global $address/${BucketsFolder}/$buck -db $Ref -strand both -id $rDNA16SID1 -evalue ${rDNAevalue1} --maxrejects ${maxrejects1} --maxhits 1 -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm should be specified here
							rm -rf $address/${BucketsFolder}/usearch.tmp
							mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
							sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
						fi
					done
				;;
				*)
					for buck in `cat $address/${BucketsFolder}/buckets_search.txt`; do
						echo -en "\r"; echo -e "# Searching against reference $Ref (id $rDNA16SID1, evalue ${rDNAevalue1} and maxrejects ${maxrejects1}) and deleting buckets... ${buck%.bk}/$buckets"
						${usearch} -usearch_global $address/${BucketsFolder}/$buck -db $Ref -strand both -id $rDNA16SID1 -evalue ${rDNAevalue1} --maxrejects ${maxrejects1} --maxhits 1 -matched $address/OUTPUT/${Out}/$buck.m7 -threads ${threads} -quiet > $address/${BucketsFolder}/usearch.tmp # Parameters of reads search by Usearch algorithm should be specified here
						rm -rf $address/${BucketsFolder}/usearch.tmp
						mv $address/OUTPUT/${Out}/$buck.m7 $address/OUTPUT/${Out}/$buck.m8
						rm -rf $address/${BucketsFolder}/$buck
						touch $address/${BucketsFolder}/$buck
						sed -i -e 1,1d $address/${BucketsFolder}/buckets_search.txt
					done
					rm -rf $address/${BucketsFolder}/*.bk
				;;
			esac
			rm -rf $address/${BucketsFolder}/*.udb; rm -rf $address/OUTPUT/${Out}/*.m7; rm -rf $address/${BucketsFolder}/buckets_search.txt
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartfilter1.tmp)" | bc -l > $address/OUTPUT/${Out}/filter1time.nmb
	rm -rf $address/OUTPUT/${Out}/datestartfilter1.tmp
	cat $address/OUTPUT/${Out}/*.m8 > $address/OUTPUT/${Out}/hits.fasta
	echo "7" > $address/CR.step; CFLR="N"
fi
}

PreLogGen () # Generates Log.txt in Output/${Out} folder, with the results from the first homology search with usearch.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "7" ]]; then echo "******Skipping generation of pre-log"; else
	echo "## Calculating statistics..."
	hits=`grep ">" $address/OUTPUT/${Out}/hits.fasta | wc -l`
	## cp $address/OUTPUT/${Out}/hits.fasta $address/OUTPUT/${Out}/1st_filter_hits.fasta
	if [[ -s $address/${BucketsFolder}/reads.nmb ]]
	then
		reads=`cat $address/${BucketsFolder}/reads.nmb`
	else
		if [[ -s $address/OUTPUT/${Out}/reads.nmb ]]
		then
			reads=`cat $address/OUTPUT/${Out}/reads.nmb`
		fi
	fi
	ppm1=`expr 1000000 \* $hits / $reads`
	touch $address/${BucketsFolder}/bucketpreviouslygeneratedmessage.tmp
	echo -e "RESULTS:
	BEAF mode: $T1 / $T2
	File from:	$R1\t$R2
	Reference:	$Ref
	SubReference(s):	$SubRef
	Results:	${Out}

	Line:	$CommandLine

	Total Reads:	$reads
	Buckets:	$buckets
	Hits:	$hits
	Portion in ppm:	$ppm1
	Percentage of hits from total number of reads:	$(echo "scale=4; $ppm1 / 10000" | bc -l | awk '{print $0 * 1.000000}')%

	Threads:	${threads}
" > $address/${BucketsFolder}/Log.txt
	if [[ -s $address/${BucketsFolder}/bucketpreviouslygeneratedmessage.tmp ]]
	then
		echo "
Trimming, Quality Analysis and Buckets generation were performed in a previous run. The time each of these process took can be seen below:" >> $address/${BucketsFolder}/Log.txt
	else
		echo "The following steps were taken in this run, prior to initial filtering, with the respective processing time for each of them:" >> $address/${BucketsFolder}/Log.txt
	fi
	echo "	Time for Trimming: $(cat $address/${BucketsFolder}/trimmingtime.nmb)
	Time for Quality Analysis: $(cat $address/${BucketsFolder}/fastqctime.nmb)
	Time to generate buckets: $(cat $address/${BucketsFolder}/bucketenginetime.nmb)s $(cat $address/${BucketsFolder}/bucketpreviouslygeneratedmessage.tmp)


Time for the first filter: $(cat $address/OUTPUT/${Out}/filter1time.nmb)s" >> $address/${BucketsFolder}/Log.txt
	echo "Time so far:	$(echo "$(date -u +%s) - $d1" |bc -l)
Time without bucket: $(echo "$(date -u +%s) - $d1 - $(cat $address/${BucketsFolder}/BucketTotalTime.nmb)" | bc -l)" >> $address/${Out}/TotalTimeAfterFirstFilter
	echo "$ppm1" > $address/OUTPUT/${Out}/ppm1.nmb
	touch $address/${BucketsFolder}/reads.nmb
	cp -r $address/${BucketsFolder}/reads.nmb $address/OUTPUT/${Out}
	mv $address/${BucketsFolder}/Log.txt $address/OUTPUT/${Out}
	if [[ -d $address/${BucketsFolder}/FASTQCresults ]]
	then
		cp -r $address/${BucketsFolder}/FASTQCresults $address/OUTPUT/${Out}
	fi
	cntg="0"
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	rm -rf $address/${BucketsFolder}/*.txt; rm -rf $address/OUTPUT/${Out}/*.m8
	echo "8" > $address/CR.step; CFLR="N"
fi
}

SubGenomesMakeDB ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "8" ]]; then echo "******Skipping generation of SubReferences databases..."; else
	if [[ -d $ReferencesFolder/$SubRef ]]
	then
		GenomeCoverage=${GenomeCoverage2}
		echo "#Making SubReference databases..."
		rm -rf $ReferencesFolder/$SubRef/SubRef_fasta.list
		ls $ReferencesFolder/$SubRef | grep -v "^SubRef_.*_cov.*.fasta$" | grep -vw "Pooled_SubRefs.fasta" | grep -Ei "(.fasta|.fa|.faa|.fas|.fna|.fsa|.ffn|.frn|.mpfa)$" | sed "s#$ReferencesFolder/$SubRef/##" | sort | uniq >> $ReferencesFolder/$SubRef/SubRef_fasta.list
		if [[ -s $ReferencesFolder/$SubRef/SubRef_fasta.list ]]
		then
			Overlap2=$(echo "$GenomeFragLength - $(echo "$GenomeFragLength / ${GenomeCoverage}" | bc)" | bc)
			for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
				if [[ ! -s $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta ]]
				then
					cat $ReferencesFolder/$SubRef/$sub | awk NF > $ReferencesFolder/$SubRef/$sub.nf.fasta
					case $GenomeSplitMethod in
						PyFasta|Pyfasta|pyfasta|PYFASTA|pyFasta|PyFASTA|pyFASTA|PYfasta|PYFasta|PF|pf|Pf|pF)
							${pyfasta} split -n 1 -k ${GenomeFragLength} -o ${Overlap2} $ReferencesFolder/$SubRef/$sub.nf.fasta > $ReferencesFolder/$SubRef/pyfastalog.txt
							mv $ReferencesFolder/$SubRef/$sub.nf.split.*.fasta $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
							rm -rf $ReferencesFolder/$SubRef/$sub.nf.fasta.flat $ReferencesFolder/$SubRef/$sub.nf.fasta.gdx $ReferencesFolder/$SubRef/$sub.nf.fasta $ReferencesFolder/$SubRef/pyfastalog.txt
						;;
						*)	
							cat $ReferencesFolder/$SubRef/$sub.nf.fasta | sed '/>/d' | tr -d '\n' | sed "s/.\{${GenomeFragLength}\}/&\n>\n/g" | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
							if [[ $(tail -n 1 $ReferencesFolder/$SubRef/$sub.cov${GenomeCoverage}.fasta | grep ">" | wc -l) -gt 0 ]]
							then
								sed -i '$ d' $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta
							fi					
						;;
					esac
					rm -rf $ReferencesFolder/$SubRef/$sub.nf.fasta
				fi
				echo "Making UDBs"
				${usearch} -makeudb_usearch $ReferencesFolder/$SubRef/SubRef_${sub%.f*}_cov${GenomeCoverage}.fasta -output $ReferencesFolder/$SubRef/${sub%.f*}.cov${GenomeCoverage}.udb -quiet > $address/${BucketsFolder}/makeudblog.txt # udb from reference is created here for genomes
				rm -rf $address/${BucketsFolder}/makeudblog.txt
				rm -rf $ReferencesFolder/$SubRef/$sub.nf*
			done
		fi
		ls $ReferencesFolder/$SubRef/*.cov${GenomeCoverage}.udb | sed "s#$ReferencesFolder/$SubRef/##" | sed "s#.cov${GenomeCoverage}.udb##" | sort -S 50% --parallel=${threads} -k1,1 > $ReferencesFolder/$SubRef/udb_cov${GenomeCoverage}_List
	fi
	echo "9" > $address/CR.step; CFLR="N"
fi
}

SubGenomesFilter2 ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" ]]; then echo "******Skipping SubReferences Genomes filtering..."; else
	if [[ -d $ReferencesFolder/$SubRef ]]
	then
		mkdir $address/OUTPUT/${Out}/SubReferences
		for sub in `cat $ReferencesFolder/$SubRef/udb_cov${GenomeCoverage}_List`; do
			mkdir $address/OUTPUT/${Out}/SubReferences/${sub}
			echo -e "# Searching against reference $SubRef/${sub} (id $GenID2, evalue ${GENevalue2} and maxrejects ${genmaxrejects2})..."
			rm -rf $address/OUTPUT/${Out}/SubReferences/$sub.hits.fasta
			${usearch} -usearch_global $address/OUTPUT/${Out}/hits.fasta -db $ReferencesFolder/$SubRef/$sub.cov${GenomeCoverage}.udb -strand both -id $GenID2 -evalue ${GENevalue2} --maxaccepts 1 --maxrejects ${genmaxrejects2} -matched $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta -threads ${threads} -quiet
			if [[ ! -s $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta ]] 
			then
				echo "  Couldn't find matches for SubReference $sub."
			else
				echo "  A total of $(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta) reads were found for subreference $sub"
			fi
		done
	fi
	echo "9" > $address/CR.step; CFLR="N"
fi
}

G_SPADES1 () # For genomes, starts SPADES process using high kmers (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" ]]; then echo "******Skipping Assembly Module with High Kmers"; else
	date -u +%s > $address/OUTPUT/${Out}/datestartspades.tmp
	echo "Starting SPADES analyses..."
	python ${spades} --threads ${threads} -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $address/OUTPUT/${Out}/hits.fasta -o $address/OUTPUT/${Out}/assembly_${Out} > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for genomes should be specified here, using high Kmers
	case $T1 in 
		G|g)
			if [[ -d $ReferencesFolder/$SubRef ]]
			then
				for sub in `cat $ReferencesFolder/$SubRef/udb_cov${GenomeCoverage}_List`; do
					if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta ]]
					then
						python ${spades} --threads ${threads} -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta -o $address/OUTPUT/${Out}/SubReferences/${sub}/assembly_${sub} > $address/OUTPUT/${Out}/logspades.txt
					fi
				done
			fi
		;;
		*)
			touch $address/OUTPUT/${Out}/assembly_${Out}
		;;
	esac
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades.tmp)" | bc -l > $address/OUTPUT/${Out}/spadestime.nmb
	rm -rf $address/OUTPUT/${Out}/logspades.txt $address/OUTPUT/${Out}/datestartspades.tmp
	echo "10" > $address/CR.step; CFLR="N"
fi
}

G_SPADES2 () # For genomes, in case the first SPADES process with high kmers didn't work, it retries the assembly process using lower kmers. Soft version skips the first step so it will try only the low kmer assembly.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" && `cat $address/CR.step` != "10" ]]; then echo "******Skipping Assembly Mode with Lower Kmers"; else
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]]
	then
		cntg=`grep ">" $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta | wc -l`
		echo "SPADES ran properly with high kmers, finding ${cntg} contigs."
	else
		echo "# SPADES couldn't find contigs for ${Out} with high kmers. Trying again with lower kmers."
		rm -rf $address/OUTPUT/${Out}/assembly_${Out}
		date -u +%s > $address/OUTPUT/${Out}/datestartspades2.tmp
		python ${spades} --threads ${threads} -k 11,15,21,25,31,35,41,45,51 --only-assembler -s $address/OUTPUT/${Out}/hits.fasta -o $address/OUTPUT/${Out}/assembly_${Out} > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for genomes should be specified here, using lower Kmers
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades2.tmp)" | bc -l > $address/OUTPUT/${Out}/spades2time.nmb
		rm -rf $address/OUTPUT/${Out}/logspades.txt $address/OUTPUT/${Out}/datestartspades2.tmp
		if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]]
		then
			cntg=`grep ">" $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta | wc -l`
			echo "SPADES ran properly with lower kmers, finding ${cntg} contigs."
		else
			echo "No contigs were found for ${Out} even when using lower kmers. SPADES assemble failed.
"
		fi
	fi
	case $T1 in 
		G|g)
			if [[ -d $ReferencesFolder/$SubRef ]]
			then
				for sub in `cat $ReferencesFolder/$SubRef/udb_cov${GenomeCoverage}_List`; do
					if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub}/assembly_${sub}/contigs.fasta ]]
					then 
						touch $address/OUTPUT/${Out}/SubReferences/${sub}/assembly_${sub}/contigs.fasta
					else
						if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta ]]
						then
							python ${spades} --threads ${threads} -k 11,15,21,25,31,35,41,45,51 --only-assembler -s $address/OUTPUT/${Out}/SubReferences/${sub}/$sub.hits.fasta -o $address/OUTPUT/${Out}/SubReferences/assembly_${sub} > $address/OUTPUT/${Out}/logspades.txt
						fi
					fi
				done
			fi
		;;
		*)
			touch $address/OUTPUT/${Out}/assembly_${Out}
		;;
	esac
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
		date -u +%s > $address/OUTPUT/${Out}/datestartquast.tmp
		if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]]
		then
			echo -e "# Analyzing draft putative genome of full reference...\n"
			cp -r $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta $address/OUTPUT/${Out}/contigs.fasta
			# tar -zcf $address/OUTPUT/${Out}/SPADES_results.tar.gz -C $address/OUTPUT/${Out}/ assembly_${Out} --remove-files
			rm -rf $address/OUTPUT/${Out}/assembly_*
			python ${quast} --threads ${threads} --silent -r $Ref -o $address/OUTPUT/${Out}/assessment $address/OUTPUT/${Out}/contigs.fasta > $address/OUTPUT/${Out}/quastlog.txt # Assessment of assemble is done in this step
			rm -rf $address/OUTPUT/${Out}/quastlog.txt
			if [[ -s $address/OUTPUT/${Out}/assessment ]]
			then
				echo "QUAST Assessment complete:
	Duplication Ratio:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/Duplication_ratio.txt | cut -d " " -f3)
	Genome Fraction:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/Genome_fraction.txt | cut -d " " -f3)
	LGA50:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/LGA50.txt | cut -d " " -f3)
	NGA50:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/NGA50.txt | cut -d " " -f3)
	Number of misassemblies:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/num_misassemblies | cut -d " " -f3)
	Mismatches per 100kbp:	$(tail -n1 $address/OUTPUT/${Out}/assessment/summary/TXT/num_mismatches_per_100_kbp.txt | cut -d " " -f3)" > $address/OUTPUT/${Out}/Quastlog.tmp
				echo "${Out} $(head -n2 $address/OUTPUT/${Out}/Quastlog.tmp | tail -n1)"
				echo -e "### Compressing results..."
				tar -zcf $address/OUTPUT/${Out}/assessment.tar.gz -C $address/OUTPUT/${Out}/ assessment --remove-files
			fi
			echo "# QUAST service is finished for the file going to OUTPUT/${Out}"
		else
			echo "# The proposed analysis of ${Out} could not continue due to problems in SPADES assembly."
		fi
		case $T1 in 
			G|g)
				if [[ -d $ReferencesFolder/$SubRef ]]
				then
					for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
						if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assembly_${sub%.f*}/contigs.fasta ]]
						then
							echo "# Analyzing draft putative genome of subreference $SubRef/${sub%.f*} using QUAST..."
							cp -r $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assembly_${sub%.f*}/contigs.fasta $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta
							rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assembly_${sub%.f*}
							python ${quast} --threads ${threads} --silent -r $ReferencesFolder/$SubRef/$sub -o $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta > $address/OUTPUT/${Out}/quastlog.txt # Assessment of assemble is done in this step
							if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment ]]
							then
								echo "QUAST Assessment of ${sub%.f*} complete:
	Duplication Ratio:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/Duplication_ratio.txt | cut -d " " -f3)
	Genome Fraction:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/Genome_fraction.txt | cut -d " " -f3)
	LGA50:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/LGA50.txt | cut -d " " -f3)
	NGA50:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/NGA50.txt | cut -d " " -f3)
	Number of misassemblies:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/num_misassemblies.txt | cut -d " " -f3)
	Mismatches per 100kbp:	$(tail -n1 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment/summary/TXT/num_mismatches_per_100_kbp.txt | cut -d " " -f3)
	Total number of hits:	$(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/${sub%.f*}.hits.fasta)
	Percentage of hits from total number of reads:	$(echo "scale=3; $(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/${sub%.f*}.hits.fasta) / $(cat $address/${BucketsFolder}/reads.nmb)" | bc -l)" > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/Quastlog.tmp
								echo "${sub%.f*} $(head -n2 $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/Quastlog.tmp | tail -n1)"
								echo -e "### Compressing results..."
								tar -zcf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/assessment.tar.gz -C $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ assessment --remove-files
							else
								echo "QUAST Assessment of subreference $SubRef/${sub%.f*} failed."
							fi
								echo "# QUAST service is finished for the file going to OUTPUT/${Out}/SubReferences/${sub%.f*}/"

						else
							echo "# The proposed analysis of ${Out} subrefere $SubRef/${sub%.f*} could not continue due to problems in SPADES assembly."
						fi
					done
				fi
			;;
			*)
				touch $address/OUTPUT/${Out}/Log.txt
			;;
		esac

		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartquast.tmp)" | bc -l > $address/OUTPUT/${Out}/quasttime.nmb
		rm -rf $address/OUTPUT/${Out}/quastlog.txt ; rm -rf $address/OUTPUT/${Out}/datestartquast.tmp
	fi
	date -u +%s > $address/OUTPUT/${Out}/datestartorffinder.tmp
	if [[ -s $address/OUTPUT/${Out}/contigs.fasta ]]
	then
		echo "# Initiating ORF finding process for file going to ${Out}"
		rm -rf $address/OUTPUT/${Out}/ORFs.${Out}.fasta
		perl ${orffinder} --infile=$address/OUTPUT/${Out}/contigs.fasta --outfile=$address/OUTPUT/${Out}/ORFs.${Out}.fasta --minlen=150 --fasta > $address/OUTPUT/${Out}/perlog.txt # Parameters for genome ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
		rm -rf $address/OUTPUT/${Out}/perlog.txt $address/OUTPUT/${Out}/datestartorffinder.tmp
		cntg=`grep ">" $address/OUTPUT/${Out}/contigs.fasta | wc -l`
		if [[ -s $address/OUTPUT/${Out}/ORFs.${Out}.fasta ]]
		then
			ORFs=`grep ">" $address/OUTPUT/${Out}/ORFs.${Out}.fasta | wc -l`
			echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
			mv $address/OUTPUT/${Out}/ORFs.${Out}.fasta $address/OUTPUT/${Out}/ORFs.fa
		else
			rm -rf $address/OUTPUT/${Out}/ORFs.${Out}.fasta
			ORFs="0"
			echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
		fi
		echo "A total of $ORFs ORFs were found for reference $Ref, from $cntg contigs"
	else
		echo "# The proposed analysis of ${Out} could not continue due to problems in SPADES assembly."
		cntg="0"
	fi
	if [[ -d $ReferencesFolder/$SubRef ]]
	then
		for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do 
			if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta ]]
			then
				echo "# Initiating ORF finding process for file going to ${Out}"
				rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta
				perl ${orffinder} --infile=$address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta --outfile=$address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta --minlen=150 --fasta > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/perlog.txt # Parameters for genome ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
				rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/perlog.txt
				grep ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta | wc -l > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/cntg.nmb
				if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta ]]
				then
					grep ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta | wc -l > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFS.nmb
					mv $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.fa
				else
					rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.${Out}.fasta
					echo "0" > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.nmb
				fi
			else
				echo "0" > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/cntg.nmb
			fi
		done
	fi
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffinder.tmp)" | bc -l > $address/OUTPUT/${Out}/orffindertime.nmb
	AvgSizeCntg="0"
	TotalSizeCntg="0"
	StdDevCntg="0"
	MaxCntg="0"
	AvgSizeORF="0"
	TotalSizeORF="0"
	StdDevORF="0"
	MaxORF="0"
	if [[ -s $address/OUTPUT/${Out}/contigs.fasta ]]
	then
		rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/contigs.cntg.sizes
		python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/contigs.fasta > $address/OUTPUT/${Out}/contigs.cntg.sizes
		AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		echo "$AvgSizeCntg|$TotalSizeCntg|$StdDevCntg|$MaxCntg" > $address/OUTPUT/${Out}/CntgStats.nmb # Calculates all statistics of contigs for the reference as a whole.
		CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
		rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/contigs.cntg.sizes
		if [[ -s $address/OUTPUT/${Out}/ORFs.fa ]]
		then
			rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/ORFs.ORF.sizes
			python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/ORFs.fa > $address/OUTPUT/${Out}/ORFs.ORF.sizes
			AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/ORFs.ORF.sizes`
			TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/ORFs.ORF.sizes`
			StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/ORFs.ORF.sizes`
			MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/ORFs.ORF.sizes`
			echo "$AvgSizeORF|$TotalSizeORF|$StdDevORF|$MaxORF" > $address/OUTPUT/${Out}/ORFStats.nmb # Calculates all statistics of ORFs for the reference as a whole.
			ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
			rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/ORFs.ORF.sizes
		else
			echo "0|0|0|0" > $address/OUTPUT/${Out}/ORFStats.nmb
		fi
	else
		echo "0|0|0|0" > $address/OUTPUT/${Out}/CntgStats.nmb
		echo "0|0|0|0" > $address/OUTPUT/${Out}/ORFStats.nmb
	fi
	if [[ -d $ReferencesFolder/$SubRef ]]
	then
		for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do 
			AvgSizeCntg="0"
			TotalSizeCntg="0"
			StdDevCntg="0"
			MaxCntg="0"
			AvgSizeORF="0"
			TotalSizeORF="0"
			StdDevORF="0"
			MaxORF="0"
			if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta ]]
			then
				rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.count; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes
				python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes
				AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes`
				TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes`
				StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes`
				MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes`
				echo "$AvgSizeCntg|$TotalSizeCntg|$StdDevCntg|$MaxCntg" > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/CntgStats.nmb # Calculates all statistics of contigs for the reference as a whole.
				rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.count; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.cntg.sizes
				if [[ -s $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.fa ]]
				then
					rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.count; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes
					python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.fa > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes
					AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes`
					TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes`
					StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes`
					MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes`
					echo "$AvgSizeORF|$TotalSizeORF|$StdDevORF|$MaxORF" > $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFStats.nmb # Calculates all statistics of ORFs for the reference as a whole.
					rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.count; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.ORF.sizes
				fi
			fi
			echo "${sub%.f*}|$(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}.hits.fasta)|$(expr 1000000 \* $(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}.hits.fasta) / $reads)|$(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/contigs.fasta)|$TotalSizeCntg|$AvgSizeCntg|$StdDevCntg|$MaxCntg|$(grep -c ">" $address/OUTPUT/${Out}/SubReferences/${sub%.f*}/ORFs.fa)|$TotalSizeORF|$AvgSizeORF|$StdDevORF|$MaxORF" >> $address/OUTPUT/${Out}/tmp1.subs
		done
	fi
	ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
	CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	echo "12" > $address/CR.step; CFLR="N"
fi
}

SoftGA () # Soft - For genomes, counts the number of contigs found after assemblage and organizes files (soft version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping calculation of contigs and output compression"; else
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]]
	then
		cntg=`grep ">" $address/OUTPUT/${Out}/contigs.fasta | wc -l`
		${pigz} -p ${threads} $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta
		mv $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta.gz $address/OUTPUT/${Out}
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
	rm -rf $ReferencesFolder/$SubRef/SubRef_fasta.list; rm -rf $ReferencesFolder/$SubRef/BlastDBlist
	ls $ReferencesFolder/$SubRef | grep -Ei "(.fasta|.fa|.faa|.fas|.fna|.fsa|.ffn|.frn|.mpfa)$" | sed "s#$ReferencesFolder/$SubRef/##" >> $ReferencesFolder/$SubRef/SubRef_fasta.list
	if [[ -s $ReferencesFolder/$SubRef/SubRef_fasta.list ]]
	then
		date -u +%s > $address/OUTPUT/${Out}/datestartblastdbgen.tmp
		echo "# Recognized $SubRef files in fasta format. Making blast databases..."
		for sub in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
			case $T1 in
				P|p)
					makeblastdb -in $ReferencesFolder/$SubRef/$sub -dbtype prot -out $ReferencesFolder/$SubRef/${sub%.f*}.bdb > $ReferencesFolder/$SubRef/makeblastdblog.txt
					rm -rf $ReferencesFolder/$SubRef/makeblastdblog.txt
				;;
				N|n)
					makeblastdb -in $ReferencesFolder/$SubRef/$sub -dbtype nucl -out $ReferencesFolder/$SubRef/${sub%.f*}.bdb > $ReferencesFolder/$SubRef/makeblastdblog.txt
					rm -rf $ReferencesFolder/$SubRef/makeblastdblog.txt
				;;
			esac
		done
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartblastdbgen.tmp)" | bc -l > $address/OUTPUT/${Out}/blastdbgentime.nmb
	fi
	case $T1 in
		P|p)
			ls $ReferencesFolder/$SubRef/*.psq | sed 's/.psq//' | sed "s#$ReferencesFolder/$SubRef/##" | sort -S 50% --parallel=${threads} -k1,1 > $ReferencesFolder/$SubRef/BlastDBlist
		;;
		N|n)
			ls $ReferencesFolder/$SubRef/*.nsq | sed 's/.nsq//' | sed "s#$ReferencesFolder/$SubRef/##" | sort -S 50% --parallel=${threads} -k1,1 > $ReferencesFolder/$SubRef/BlastDBlist
		;;
	esac
	echo "9" > $address/CR.step; CFLR="N"
fi
}

Filter2 () # For proteins and genes, makes the second homology search, using blast, searching against each specific subreference.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "9" ]]; then echo "******Skipping Second Filtering System"; else
	date -u +%s > $address/OUTPUT/${Out}/datestartfilter2.tmp
	for sub in `cat $ReferencesFolder/$SubRef/BlastDBlist`; do
		echo "# Searching against ${sub#"$ReferencesFolder/$SubRef/"} at evalue ${PROTevalue2} and filtering for ${ProtID2} identity and at least ${PROTqcovHits}% query coverages..."
		date -u +%s > $address/OUTPUT/${Out}/datestarthitsfiltering.tmp
		case $T1 in
			P|p)
				blastx -db $ReferencesFolder/$SubRef/$sub -query $address/OUTPUT/${Out}/hits.fasta -out $ReferencesFolder/$SubRef/${sub%.bdb}.tmp -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for proteins should be specified here
			;;
			N|n)
				blastn -db $ReferencesFolder/$SubRef/$sub -query $address/OUTPUT/${Out}/hits.fasta -out $ReferencesFolder/$SubRef/${sub%.bdb}.tmp -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for genes should be specified here
			;;
		esac
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestarthitsfiltering.tmp)" | bc -l >> $address/OUTPUT/${Out}/${sub%.bdb}.ft.time2
		cat $ReferencesFolder/$SubRef/${sub%.bdb}.tmp | awk -v identity=${ProtID2#*.} -v qcov=${PROTqcovHits} -v evalue=${PROTevalue2} '{if ($3>=identity && $11<=evalue && $12>=qcov) print $0}' | sort -S 50% --parallel=${threads} -k3,4 -n -r | uniq > $ReferencesFolder/$SubRef/${sub%.bdb}.ft
		touch $ReferencesFolder/$SubRef/${sub%.bdb}.ft
		rm -rf $ReferencesFolder/$SubRef/${sub%.bdb}.tmp ; rm -rf $ReferencesFolder/$SubRef/hits
		sed -i -e 1,1d $ReferencesFolder/$SubRef/BlastDBlist
	done
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartfilter2.tmp)" | bc -l > $address/OUTPUT/${Out}/filter2time.nmb
	rm -rf $ReferencesFolder/$SubRef/BlastDBlist; rm -rf $ReferencesFolder/$SubRef/*.tmp ; rm -rf $address/OUTPUT/${Out}/datestartfilter2.tmp
	echo "10" > $address/CR.step; CFLR="N"
fi
}

Extraction () # For proteins and genes, moves files to OUTPUT/${Out} to continue the process.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "10" ]]; then echo "******Skipping Extraction Module"; else
	echo "# Arranging data..."
	mv $ReferencesFolder/$SubRef/*.ft $address/OUTPUT/${Out}
	mkdir $address/OUTPUT/${Out}/fasta_hits $address/OUTPUT/${Out}/blast_hits $address/OUTPUT/${Out}/contigs $address/OUTPUT/${Out}/blast_contigs $address/OUTPUT/${Out}/ORFs $address/OUTPUT/${Out}/blast_ORFs
	# mkdir $address/OUTPUT/${Out}/contigs/SPADES
	echo "11" > $address/CR.step; CFLR="N"
fi
}

PN_Prepare_SPADES () # For proteins and genes, prepares files for SPADES assembly, using python to extract blast matches from the hits, removing redundancy.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping Assembly Preparation"; else
	rm -rf $address/OUTPUT/${Out}/*.rev; rm -rf $address/OUTPUT/${Out}/*.hits; rm -rf $address/OUTPUT/${Out}/*.rev-hits.fasta
	cut -f 1 $address/OUTPUT/${Out}/$File > $address/OUTPUT/${Out}/$File.rev
	python ${extpy} $address/OUTPUT/${Out}/$File.rev $address/OUTPUT/${Out}/hits.fasta $address/OUTPUT/${Out}/$File.rev-hits.fasta > $address/OUTPUT/${Out}/extpylog.txt
	cat $address/OUTPUT/${Out}/tsv_b6.header.txt $address/OUTPUT/${Out}/$File > $address/OUTPUT/${Out}/blast_hits/$File
	rm -rf $address/OUTPUT/${Out}/extpylog.txt $address/OUTPUT/${Out}/$File.rev
	touch $address/OUTPUT/${Out}/$File.rev-hits.fasta
	if [[ -s $address/OUTPUT/${Out}/$File.rev-hits.fasta ]]
	then
		sq=`grep ">" $address/OUTPUT/${Out}/$File.rev-hits.fasta | wc -l`
		ppm2=`expr 1000000 \* $sq / $reads`
		${cdhit} -i $address/OUTPUT/${Out}/$File.rev-hits.fasta -o $address/OUTPUT/${Out}/$File.rev-hits -c 1.00 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5 > $address/OUTPUT/${Out}/cdhitlog # Parameters for CD-Hit (removing redundancy) should be specified here
		rm -rf $address/OUTPUT/${Out}/$File.rev-hits.fasta $address/OUTPUT/${Out}/$File.rev-hits.clstr $address/OUTPUT/${Out}/cdhitlog
	else
		rm -rf $address/OUTPUT/${Out}/$File.rev-hits.fasta
		sq="0"
		ppm2="0"
	fi
	echo "$sq" > $address/OUTPUT/${Out}/sq.nmb
	echo "$ppm2" > $address/OUTPUT/${Out}/ppm2.nmb
	cp -r $address/OUTPUT/${Out}/$File.rev-hits.fasta $address/OUTPUT/${Out}/fasta_hits/$File.fasta
	rm -rf $address/OUTPUT/${Out}/$File ; rm -rf $address/OUTPUT/${Out}/*.rev-hits.fasta
	echo "12" > $address/CR.step; CFLR="N"
fi
}

PN_SPADES1 () # For proteins and genes, starts SPADES process using high kmers (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" ]]; then  echo "Skipping Assembly with High Kmers"; else
	rm -rf $address/OUTPUT/${Out}/assembly_*
	date -u +%s > $address/OUTPUT/${Out}/datestartspades.tmp
	rm -rf $address/OUTPUT/${Out}/assembly_${Out}_$File
	python ${spades} --threads ${threads} -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $address/OUTPUT/${Out}/fasta_hits/$File.fasta -o $address/OUTPUT/${Out}/assembly_${Out}_$File > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for proteins and genes should be specified here, using high Kmers
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades.tmp)" | bc -l >> $address/OUTPUT/${Out}/spadestime.nmb
	rm -rf $address/OUTPUT/${Out}/logspades.txt ; rm -rf $address/OUTPUT/${Out}/datestartspades.tmp
	echo "13" > $address/CR.step; CFLR="N"
fi
}

PN_SPADES2 () # For proteins and genes, in case the first SPADES process with high kmers didn't work, it retries the assembly process using lower kmers. Soft version skips the first step so it will try only the low kmer assembly.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" && `cat $address/CR.step` != "13" ]]; then echo "******Skipping Assembly with Low Kmers"; else
	if [[ -s $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta ]]
	then
		echo "# No need to try with lower kmers "
		echo "0" >> $address/OUTPUT/${Out}/spades2time.nmb
	else
		echo "# Trying for ${File%.ft} for ${Out} with lower kmers"
		if [[ -s $address/OUTPUT/${Out}/fasta_hits/$File.fasta ]]
		then
			sq=`cat $address/OUTPUT/${Out}/sq.nmb`
			ppm2=`cat $address/OUTPUT/${Out}/ppm2.nmb`
		else
			sq="0"
			ppm2="0"
		fi
		rm -rf $address/OUTPUT/${Out}/assembly_${Out}_$File
		date -u +%s > $address/OUTPUT/${Out}/datestartspades2.tmp
		python ${spades} --threads ${threads} -k 9,11,13,15,17,19,21,31 --only-assembler -s $address/OUTPUT/${Out}/fasta_hits/$File.fasta -o $address/OUTPUT/${Out}/assembly_${Out}_$File > $address/OUTPUT/${Out}/logspades.txt # Parameters for SPADES assembly for proteins and genes should be specified here, using lower Kmers
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartspades2.tmp)" | bc -l >> $address/OUTPUT/${Out}/spades2time.nmb
		rm -rf $address/OUTPUT/${Out}/logspades.txt ; rm -rf $address/OUTPUT/${Out}/datestartspades2.tmp
	fi
	echo "14" > $address/CR.step; CFLR="N"
fi
}

PNA () # For proteins and genes, arranges data from SPADES assembly into OUTPUT/${Out}, finding the number of contigs for each subreference.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "14" ]]; then echo "******Skipping calculation of contigs"; else
	if [ -s $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta ]
	then
		cntg=`grep ">" $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta | wc -l`
		echo "# SPADES worked on ${File%.ft} for ${Out}, finding $cntg contigs. Blasting those contigs against your subreference file now."
		date -u +%s > $address/OUTPUT/${Out}/datestartcontigfiltering.tmp
		case $T1 in
			P|p)
				blastx -db $ReferencesFolder/$SubRef/${File%.ft}.bdb -query $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta -out $address/OUTPUT/${Out}/contigs/contigs.${File}_blast_contigs.tsv -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for proteins should be specified here
			;;
			N|n)
				blastn -db $ReferencesFolder/$SubRef/${File%.ft}.bdb -query $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta -out $address/OUTPUT/${Out}/contigs/contigs.${File}_blast_contigs.tsv -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for genes should be specified here
			;;
		esac
		cat $address/OUTPUT/${Out}/contigs/contigs.${File}_blast_contigs.tsv | awk -v identity=${ProtID2#*.} -v qcov=${PROTqcovHits} '{if ($3>=identity && $12>=qcov) print $0}' | sort -S 50% --parallel=${threads} -k3,4 -n -r | uniq > $address/OUTPUT/${Out}/contigs/contigs.${File}_qcovs_contigs.tsv
		rm -rf $address/OUTPUT/${Out}/contigs/contigs.${File}_blast_contigs.tsv
		cut -f 1 $address/OUTPUT/${Out}/contigs/contigs.${File}_qcovs_contigs.tsv > $address/OUTPUT/${Out}/contigs/contigs.${File}_headers.tsv
		python ${extpy} $address/OUTPUT/${Out}/contigs/contigs.${File}_headers.tsv $address/OUTPUT/${Out}/assembly_${Out}_$File/contigs.fasta $address/OUTPUT/${Out}/contigs/contigs.$File.fasta > $address/OUTPUT/${Out}/extpylog.txt
		rm -rf $address/OUTPUT/${Out}/contigs/contigs.${File}_headers.tsv $address/OUTPUT/${Out}/extpylog.txt
		cat $address/OUTPUT/${Out}/tsv_b6.header.txt $address/OUTPUT/${Out}/contigs/contigs.${File}_qcovs_contigs.tsv > $address/OUTPUT/${Out}/blast_contigs/contigs.${File}.b6out_qcovs.tsv
		rm -rf $address/OUTPUT/${Out}/contigs/contigs.${File}_qcovs_contigs.tsv
		cntg=`grep ">" $address/OUTPUT/${Out}/contigs/contigs.$File.fasta | wc -l`
		echo "# After blasting those contigs, it resulted on a total of $cntg contigs."
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartcontigfiltering.tmp)" | bc -l >> $address/OUTPUT/${Out}/contigfilteringtime.nmb
	else
		echo "# The proposed analysis could not continue due to problems in SPADES assembly."
		echo "0" >> $address/OUTPUT/${Out}/contigfilteringtime.nmb
		cntg="0"
		Warnings="WARNING: Did not run SPADES properly"
	fi
	echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
	# tar -zcf $address/OUTPUT/${Out}/contigs/SPADES/SPADES_assembly_$File.tar.gz -C $address/OUTPUT/${Out}/ assembly_${Out}_$File --remove-files
	rm -rf $address/OUTPUT/${Out}/assembly_*
	echo "15" > $address/CR.step; CFLR="N"
fi
}

PNORFs () # For proteins and genes, finds ORFs in each contig (full version).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "15" ]]; then echo "******Skipping ORF finding process"; else
	if [[ -s $address/OUTPUT/${Out}/contigs/contigs.$File.fasta ]]
	then
		echo "# Initiating ORF finding process for contigs of subreference ${File%.ft}"
		rm -rf $address/OUTPUT/${Out}/${File}_ORFs.fasta
		date -u +%s > $address/OUTPUT/${Out}/datestartorffinder.tmp
		perl ${orffinder} --infile=$address/OUTPUT/${Out}/contigs/contigs.$File.fasta --outfile=$address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta --minlen=150 --fasta > $address/OUTPUT/${Out}/logorffinder.txt # Parameters for protein/gene ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffinder.tmp)" | bc -l >> $address/OUTPUT/${Out}/orffindertime.nmb
		rm -rf $address/OUTPUT/${Out}/logorffinder.txt ; rm -rf $address/OUTPUT/${Out}/datestartorffinder.tmp
		ORFs=$(grep ">" $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta | wc -l)
	else
		if [[ -s $address/OUTPUT/${Out}/fasta_hits/$File.fasta ]]
		then
			echo "# Initiating ORF finding process for hits to subreference ${File%.ft}. As SPADES could not generate contigs, the sequences that had hits with the subreference will be used to search for ORFs."
			rm -rf $address/OUTPUT/${Out}/${File}_ORFs.fasta
			date -u +%s > $address/OUTPUT/${Out}/datestartorffinder.tmp
			perl ${orffinder} --infile=$address/OUTPUT/${Out}/fasta_hits/$File.fasta --outfile=$address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta --minlen=150 --fasta > $address/OUTPUT/${Out}/logorffinder.txt # Parameters for protein/gene ORF finding should be specified here. If user wants to find orfs bigger or smaller just change parameter "minlen" to the minimum length required
			echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffinder.tmp)" | bc -l >> $address/OUTPUT/${Out}/orffindertime.nmb
			rm -rf $address/OUTPUT/${Out}/logorffinder.txt ; rm -rf $address/OUTPUT/${Out}/datestartorffinder.tmp
			ORFs=$(grep ">" $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta | wc -l)
		else
			ORFs="0"
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta ]]
	then
		echo "Found $ORFs ORFs for ${File%.ft}. Blasting those ORFs against it's subreference at evalue ${PROTevalue2} and filtering for ${ProtID2} identity and at least ${PROTqcovORFs}% query coverages... "
		date -u +%s > $address/OUTPUT/${Out}/datestartorffiltering.tmp
		case $T1 in
			P|p)
				blastx -db $ReferencesFolder/$SubRef/${File%.ft}.bdb -query $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta -out $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.tsv -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for proteins should be specified here
			;;
			N|n)
				blastn -db $ReferencesFolder/$SubRef/${File%.ft}.bdb -query $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta -out $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.tsv -evalue ${PROTevalue2} -strand both -max_target_seqs 1 -num_threads $threads -outfmt '6 std qcovs' # Parameters of reads search by blast for genes should be specified here
			;;
		esac
		cat $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.tsv | awk -v identity=${ProtID2#*.} -v qcov=${PROTqcovORFs} '{if ($3>=identity && $12>=qcov) print $0}' | sort -S 50% --parallel=${threads} -k3,4 -n -r | uniq > $address/OUTPUT/${Out}/ORFs/${File}_qcovs_ORFs.tsv
		rm -rf $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.tsv
		cut -f 1 $address/OUTPUT/${Out}/ORFs/${File}_qcovs_ORFs.tsv > $address/OUTPUT/${Out}/ORFs/${File}_headers.tsv
		python ${extpy} $address/OUTPUT/${Out}/ORFs/${File}_headers.tsv $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta > $address/OUTPUT/${Out}/extpylog.txt
		rm -rf $address/OUTPUT/${Out}/ORFs/${File}_headers.tsv $address/OUTPUT/${Out}/extpylog.txt
		cat $address/OUTPUT/${Out}/tsv_b6.header.txt $address/OUTPUT/${Out}/ORFs/${File}_qcovs_ORFs.tsv > $address/OUTPUT/${Out}/blast_ORFs/ORFs.${File}.b6out_qcovs.tsv
		rm -rf $address/OUTPUT/${Out}/ORFs/${File}_qcovs_ORFs.tsv
		ORFs=`grep ">" $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta | wc -l`
		echo "# After blasting those ORFs, it resulted on a total of $ORFs ORFs."
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartorffiltering.tmp)" | bc -l >> $address/OUTPUT/${Out}/orffilteringtime.nmb
		rm -rf $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.fasta $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta $address/OUTPUT/${Out}/ORFs/${File}_qcovs_ORFs.tsv ; rm -rf $address/OUTPUT/${Out}/datestartorffiltering.tmp
	else
		echo "0" >> $address/OUTPUT/${Out}/orffilteringtime.nmb
		ORFs="0"
	fi
	if [[ -s $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta ]]
	then
		ORFs=`grep ">" $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta | wc -l`
	else
		ORFs="0"
	fi
	rm -rf $address/OUTPUT/${Out}/ORFs/${File}_blast_ORFs.fasta $address/OUTPUT/${Out}/ORFs/${File}_not_filtered_ORFs.fasta ; rm -rf $address/OUTPUT/${Out}/datestartorffiltering.tmp
	if [[ -s cntg.nmb ]]
	then 
		cntg=`cat cntg.nmb`
	fi
	echo "A total of $ORFs ORFs were found for file $File, from $cntg contigs"
	echo "$ORFs" > $address/OUTPUT/${Out}/ORFs.nmb
	echo "16" > $address/CR.step; CFLR="N"
fi
}

PN_CalcStats () # For proteins and genes, calculates statistics of contigs and ORFs for each subreference (average sizes, total sizes of all ORFs/contigs, standard deviation of sizes and maximum size of contigs and ORFs).
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "16" ]]; then echo "******Skipping Calculation of Statistics for Contigs and ORFs"; else
	AvgSizeCntg="0"
	TotalSizeCntg="0"
	StdDevCntg="0"
	AvgSizeORF="0"
	TotalSizeORF="0"
	StdDevORF="0"
	MaxCntg="0"
	MaxORF="0"
	if [[ -s $address/OUTPUT/${Out}/contigs/contigs.$File.fasta ]]
	then
		rm -rf $address/OUTPUT/${Out}/contigs/*.count; rm -rf $address/OUTPUT/${Out}/contigs/*.cntg; rm -rf $address/OUTPUT/${Out}/contigs/$File.cntg.sizes
		python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/contigs/contigs.$File.fasta > $address/OUTPUT/${Out}/contigs/$File.cntg.sizes
		AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs/$File.cntg.sizes`
		if [[ -s $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta ]]
		then
			rm -rf $address/OUTPUT/${Out}/ORFs/*.count; rm -rf $address/OUTPUT/${Out}/ORFs/*.cntg; rm -rf $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes
			python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/ORFs/${File}_ORFs.fasta > $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes
			AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
			MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/ORFs/$File.ORF.sizes`
		fi
	fi
	echo "17" > $address/CR.step; CFLR="N"
fi
}

SaveDBs () # For proteins and genes, keeps subreference DBs in case KeepBlastDBs is set to 'Y', moving fasta files to a folder used just to store them. If KeepBlastDBs is 'N', removes blast DBs in this step.
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "15" && `cat $address/CR.step` != "17" ]]; then echo "******Skipping the process of saving/removing databases"; else
	case $KeepBlastDBs in
		Y|y|YES|Yes|yes)
			if [[ -s $ReferencesFolder/$SubRef/SubRef_fasta.list ]]
			then
				if [[ -d $ReferencesFolder/$SubRef/Fasta_files ]]
				then
					touch $ReferencesFolder/$SubRef/Fasta_files
				else
					mkdir $ReferencesFolder/$SubRef/Fasta_files
				fi
				echo "# Databases of subreference $SubRef now saved to References_seqs folder in blastdb format. Fasta files used to make the databases have been realocated to $ReferencesFolder/$SubRef/Fasta_files."
				for file in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
					mv $ReferencesFolder/$SubRef/$file $ReferencesFolder/$SubRef/Fasta_files
				done
			fi
		;;
		*)
			echo "# Removing blastdb databases generated using fasta files in the subreference folder (Reference_seqs/$SubRef)"
			if [[ -s $ReferencesFolder/$SubRef/SubRef_fasta.list ]]
			then
				for file in `cat $ReferencesFolder/$SubRef/SubRef_fasta.list`; do
					case $T1 in
						P|p)
							rm -rf $ReferencesFolder/$SubRef/$file.bdb.p*
						;;
						N|n)
							rm -rf $ReferencesFolder/$SubRef/$file.bdb.n*
						;;
					esac
				done
			fi
		;;
	esac
	rm -rf $ReferencesFolder/$SubRef/SubRef_fasta.list
	echo "18" > $address/CR.step; CFLR="N"
fi
}

rDNA_Contig ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "11" ]]; then echo "******Skipping contig analysis"; else
	AvgSizeCntg="0"
	TotalSizeCntg="0"
	StdDevCntg="0"
	MaxCntg="0"
	if [ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]
	then
		cp -r $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta $address/OUTPUT/${Out}
	fi
	if [[ -s $address/OUTPUT/${Out}/contigs.fasta ]]
	then
		cntg=`grep "^>" $address/OUTPUT/${Out}/contigs.fasta | wc -l`
		echo "$cntg" > $address/OUTPUT/${Out}/cntg.nmb
		rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/contigs.cntg.sizes
		python ${SeqLength} --OnlyPrintLength $address/OUTPUT/${Out}/contigs.fasta > $address/OUTPUT/${Out}/contigs.cntg.sizes
		AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs.cntg.sizes`
		echo "$AvgSizeCntg|$TotalSizeCntg|$StdDevCntg|$MaxCntg" > $address/OUTPUT/${Out}/CntgStats.nmb # Calculates all statistics of contigs for the reference as a whole.
		CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
		rm -rf $address/OUTPUT/${Out}/*.count; rm -rf $address/OUTPUT/${Out}/*.cntg; rm -rf $address/OUTPUT/*.sizes ; rm -rf $address/OUTPUT/${Out}/contigs.cntg.sizes
	fi
	echo "12" > $address/CR.step; CFLR="N"
fi
}

rDNA_Chimeras ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" ]]; then echo "******Skipping Dechimerization"; else
	if [ -s $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta ]
	then
		cp -r $address/OUTPUT/${Out}/assembly_${Out}/contigs.fasta $address/OUTPUT/${Out}
	fi
	# tar -zcf $address/OUTPUT/${Out}/SPADES_results.tar.gz -C $address/OUTPUT/${Out}/ assembly_${Out} --remove-files
	rm -rf $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/uparse.txt $address/OUTPUT/${Out}/otus.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16S.relabeled.fa ; rm -rf $address/OUTPUT/${Out}/assembly_*
	echo "# Analyzing 16S...
Clustering sequences to 97% identity..."
	date -u +%s > $address/OUTPUT/${Out}/datestartcluster.tmp
	python ${SeqLength} --Cut-off-MinLength 75 $address/OUTPUT/${Out}/hits.fasta > $address/OUTPUT/${Out}/16S.fa #only selects sequences with a length of 75 or more
	mkdir $address/OUTPUT/${Out}/contigs
	python ${SeqLength} --Cut-off-MinLength ${rDNAminCntgLength} $address/OUTPUT/${Out}/contigs.fasta | sed 's/_length_/;length=/' | sed 's/_cov_.*/;size=&;/' | sed 's/_cov_//' | awk -v Out=${Out} -F "[;=>]" '/^>/{round=sprintf("%.0f",$6);print ">" Out ";" $2 ";" $3 "=" $4 ";" $5 "=" round ; next}{print}' > $address/OUTPUT/${Out}/contigs/contigs_over_${rDNAminCntgLength}.fa
	python ${SeqLength} --Cut-off-MaxLength ${rDNAminCntgLength} $address/OUTPUT/${Out}/contigs.fasta | sed 's/_length_/;length=/' | sed 's/_cov_.*/;size=&;/' | sed 's/_cov_//' | awk -v Out=${Out} -F "[;=>]" '/^>/{round=sprintf("%.0f",$6);print ">" Out ";" $2 ";" $3 "=" $4 ";" $5 "=" round ";Discarded_due_to_length" ; next}{print}' > $address/OUTPUT/${Out}/contigs/contigs_under_${rDNAminCntgLength}.fa
	${usearch} -fastx_relabel $address/OUTPUT/${Out}/16S.fa -prefix "${Out};" -fastaout $address/OUTPUT/${Out}/16S.relabeled.fa -keep_annots -quiet > $address/OUTPUT/${Out}/usearchlog.txt
	${usearch} -cluster_fast $address/OUTPUT/${Out}/16S.relabeled.fa -id ${rDNA16SID2} --maxaccepts 0 --maxrejects 0 -sizeout -centroids $address/OUTPUT/${Out}/16Snr.fa -threads ${threads}  -quiet  > $address/OUTPUT/${Out}/usearchlog.txt 
	${usearch} -cluster_otus $address/OUTPUT/${Out}/16Snr.fa -otus $address/OUTPUT/${Out}/otus.fa -uparseout $address/OUTPUT/${Out}/uparse.txt -relabel "${Out};" -minsize 2 -quiet  > $address/OUTPUT/${Out}/usearchlog.txt
	grep -w "OTU" $address/OUTPUT/${Out}/uparse.txt | grep -vw "chimera" | cut -f 1 > $address/OUTPUT/${Out}/headers.txt
	python ${extpy} $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16S.nonchimera > $address/OUTPUT/${Out}/extpylog.txt
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartcluster.tmp)" | bc -l > $address/OUTPUT/${Out}/clustertime.nmb
	rm -rf $address/OUTPUT/${Out}/headers.txt $address/OUTPUT/${Out}/uparse.txt $address/OUTPUT/${Out}/otus.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16S.relabeled.fa $address/OUTPUT/${Out}/datestartcluster.tmp $address/OUTPUT/${Out}/usearchlog.txt $address/OUTPUT/${Out}/extpylog.txt
	echo "13" > $address/CR.step; CFLR="N"
fi
}

rDNA_Filter2 ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "13" ]]; then echo "******Skipping GreenFiltering"; else
	echo ""
	if [ -s $address/OUTPUT/${Out}/otu_b6out.tsv ]
	then
		touch $address/OUTPUT/${Out}/otu_b6out.tsv
	else
		rm -rf $address/OUTPUT/${Out}/16S.tsv $address/OUTPUT/${Out}/16S.nm7 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.m7 $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.seq $address/OUTPUT/${Out}/16S.fa $address/OUTPUT/${Out}/16Snr.fa $address/OUTPUT/${Out}/16Snr.fasta $address/OUTPUT/${Out}/otu_blast6out.tsv $address/OUTPUT/${Out}/otu_b6out.tsv $address/OUTPUT/${Out}/16S.table
		date -u +%s > $address/OUTPUT/${Out}/datestart16sfilter2.tmp
		if [[ ! -s $ReferencesFolder/$SubRef ]]
		then
				if [[ -s $ReferencesFolder/nr100_Greengenes.udb ]]
				then
					echo "Using the full Green Genes database (clustered to 100% identity) to search for OTUs."
					SubRef=nr100_Greengenes.udb
				else
					echo "Couldn't find subreference $SubRef in References Folder $ReferencesFolder to search for OTUs."
				fi
		fi
		if [[ -s $address/OUTPUT/${Out}/16S.nonchimera ]]
		then
			echo "Starting second filter (searching for homology to $SubRef database) at ${rDNA16SID2} identity, evalue of ${rDNAevalue2}. Using maxaccepts ${maxaccepts2} and maxrejects ${maxrejects2}"
			${usearch} -usearch_global $address/OUTPUT/${Out}/16S.nonchimera -db $ReferencesFolder/$SubRef -sizein -sizeout -strand plus -id ${rDNA16SID2} -evalue ${rDNAevalue2} --maxaccepts ${maxaccepts2} --maxrejects ${maxrejects2} -maxhits 1 -matched $address/OUTPUT/${Out}/16S.m7 -notmatched $address/OUTPUT/${Out}/16S.nm7 -otutabout $address/OUTPUT/${Out}/otu_table.seqidnum -blast6out $address/OUTPUT/${Out}/otu_blast6out.tsv -query_cov 0.${rDNA16Sqcov} -threads ${threads} -quiet -uc $address/OUTPUT/${Out}/map_nobarcode.uc > $address/OUTPUT/${Out}/usearchlog.txt  
			mv $address/OUTPUT/${Out}/otu_blast6out.tsv $address/OUTPUT/${Out}/otu_b6out.tsv
		fi
		if [[ -s $address/OUTPUT/${Out}/contigs/contigs_over_${rDNAminCntgLength}.fa ]]
		then
			echo "Starting filter for contigs, with the same options."
			${usearch} -usearch_global $address/OUTPUT/${Out}/contigs/contigs_over_${rDNAminCntgLength}.fa -db $ReferencesFolder/$SubRef -sizein -sizeout -strand plus -id ${rDNA16SID2} -evalue ${rDNAevalue2} --maxaccepts ${maxaccepts2} --maxrejects ${maxrejects2} -maxhits 1 -matched $address/OUTPUT/${Out}/contigs/16S_contig.m7 -notmatched $address/OUTPUT/${Out}/contigs/16S_contig.nm7 -otutabout $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum -blast6out $address/OUTPUT/${Out}/contigs/otu_blast6out_contig.tsv -query_cov 0.${rDNA16Sqcov} -threads ${threads} -quiet -uc $address/OUTPUT/${Out}/contigs/map_nobarcode_contig.uc > $address/OUTPUT/${Out}/usearchlog.txt
			mv $address/OUTPUT/${Out}/contigs/otu_blast6out_contig.tsv $address/OUTPUT/${Out}/contigs/otu_b6out_contig.tsv
		fi
		echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestart16sfilter2.tmp)" | bc -l > $address/OUTPUT/${Out}/16sfilter2time.nmb
		rm -rf 	$address/OUTPUT/${Out}/datestart16sfilter2.tmp $address/OUTPUT/${Out}/usearchlog.txt

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
		totalmatched=`cat $address/OUTPUT/${Out}/16S.table | wc -l`
		while read B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13; do
			name=`grep ^">$B3 " $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.headers.txt  | sed 's/.*.k__/k__/g'`
			new=">16S_$counter|ID_$B4|AligLength_$B5|eval_$B12|Bitscore_$B13|Size_${B2%;}|$name"
			echo -en "\r"; echo -en "Renaming matched sequences ($counter/$totalmatched)   "
			sed -i "s@^>16S_$B1;size=.*@$new@" $address/OUTPUT/${Out}/16S.m7
			((counter+=1))
			sed -i 1,1d $address/OUTPUT/${Out}/16S.table
		done < $address/OUTPUT/${Out}/16S.table
		echo -e "
Finished renaming matched sequences"
	fi
	if [[ -s $address/OUTPUT/${Out}/contigs/16S_contig.m8 ]]
	then
		touch $address/OUTPUT/${Out}/contigs/16S_contig.m8
	else
		if [[ -s $address/OUTPUT/${Out}/contigs/16S_contig.table ]]
		then 
			touch $address/OUTPUT/${Out}/contigs/16S_contig.table
		else
			if [[ -s $address/OUTPUT/${Out}/contigs/16S_contig.m7 ]]
			then
				sed -i 's/>/>16S_/' $address/OUTPUT/${Out}/contigs/16S_contig.m7
				cat $address/OUTPUT/${Out}/contigs/otu_b6out_contig.tsv > $address/OUTPUT/${Out}/contigs/16S_contig.table
			fi
		fi
		counter="1"
		totalmatched=`cat $address/OUTPUT/${Out}/contigs/16S_contig.table | wc -l`
		if [[ -s $address/OUTPUT/${Out}/contigs/16S_contig.table ]]
		then
			while read B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B12; do
				name=`grep ^">$B2 " $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.headers.txt  | sed 's/.*.k__/k__/g'`
				echo -en "\r"; echo -en "Renaming matched contig sequences ($counter/$totalmatched)   "
				sed -i "s@^>$B1.*@&;$name;ID_$B3;AligLength_$B5;eval_$B12;Bitscore_$B13@" $address/OUTPUT/${Out}/contigs/16S_contig.m7
				((counter+=1))
				sed -i 1,1d $address/OUTPUT/${Out}/contigs/16S_contig.table
			done < $address/OUTPUT/${Out}/contigs/16S_contig.table
			echo -e "
	Finished renaming matched contig sequences"
		else
			echo "No contigs matched with reference $Ref"
		fi
	fi
	cat $address/OUTPUT/${Out}/16S.m7 > $address/OUTPUT/${Out}/16S.m8
	cat $address/OUTPUT/${Out}/contigs/16S_contig.m7 > $address/OUTPUT/${Out}/contigs/16S_contig.m8
	if [[ -s $address/OUTPUT/${Out}/16S.nm8 ]]
	then
		touch $address/OUTPUT/${Out}/16S.nm8
	else
		UnknownN=`grep ">" $address/OUTPUT/${Out}/16S.nm7 | wc -l`
		echo "Renaming $UnknownN Unknown sequences."
		awk -F ";" '/^>/{print ">16S_Unknown_" ++i ";" $3 ; next}{print}' < $address/OUTPUT/${Out}/16S.nm7 > $address/OUTPUT/${Out}/16S.nm8
		echo -e "Finished renaming Unknowns."
	fi
	mv $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.fasta
	cat $address/OUTPUT/${Out}/contigs/contigs_under_${rDNAminCntgLength}.fa $address/OUTPUT/${Out}/16S.nm8 > $address/OUTPUT/${Out}/DiscardedSequences_possible16S.fasta
	rm -rf $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/contigs/contigs_under_${rDNAminCntgLength}.fa
	mv $address/OUTPUT/${Out}/contigs/16S_contig.m8 $address/OUTPUT/${Out}/contigs/16S_contigs.fasta
#	grep ">" $address/OUTPUT/${Out}/contigs/16S_contigs.fasta | wc -l
	mv $address/OUTPUT/${Out}/contigs/16S_contig.nm7 $address/OUTPUT/${Out}/contigs/Discarded_Contig_Sequences_possible16S.fasta
	cat $address/OUTPUT/${Out}/map_nobarcode.uc | awk -F "\t" -v Out=${Out} '{if ($9 ~ /;$/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "barcodelabel=" Out ";\t" $10 ; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 ";barcodelabel=" Out ";\t" $10}' > $address/OUTPUT/${Out}/map.uc
	if [[ "$(cat $address/OUTPUT/${Out}/otu_table.seqidnum | wc -l)" -le 2 ]]
	then
		echo "#OTU ID	${Out}" > $address/OUTPUT/${Out}/otu_table.seqidnum
		cut -f 9,10 $address/OUTPUT/${Out}/map.uc | sed 's#.*;size=##' | sed 's#;barcodelabel=.*;\t#\t>#' | grep -v "\*" | sed 's# .*##' | awk '{seen[$2]+=$1}END{ for (id in seen) print id "\t" seen[id] }'  >> $address/OUTPUT/${Out}/otu_table.seqidnum
	fi
	cat $address/OUTPUT/${Out}/contigs/map_nobarcode_contig.uc | awk -F "\t" -v Out=${Out} '{if ($9 ~ /;$/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "barcodelabel=" Out ";\t" $10 ; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 ";barcodelabel=" Out ";\t" $10}' > $address/OUTPUT/${Out}/contigs/map_contig.uc
	if [[ "$(cat $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum | wc -l)" -le 2 ]]
	then
		echo "#OTU ID	${Out}" > $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum
		cut -f 9,10 $address/OUTPUT/${Out}/contigs/map_contig.uc | sed 's#.*;size=##' | sed 's#;barcodelabel=.*;\t#\t>#' | grep -v "\*" | sed 's# .*##' | awk '{seen[$2]+=$1}END{ for (id in seen) print id "\t" seen[id] }'  >> $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum
	fi
	echo "$(date -u +%s) - $(cat $address/OUTPUT/${Out}/datestartrenaming.tmp)" | bc -l > $address/OUTPUT/${Out}/renamingtime.nmb
	rm -rf $address/OUTPUT/${Out}/16S.table $address/OUTPUT/${Out}/16S.m7 $address/OUTPUT/${Out}/16S.m8 $address/OUTPUT/${Out}/16S.nm8 $address/OUTPUT/${Out}/16S.nm7 $address/OUTPUT/${Out}/16S.tsv $address/OUTPUT/${Out}/map_nobarcode.uc $address/OUTPUT/${Out}/contigs/contigs_over_${rDNAminCntgLength}.fa $address/OUTPUT/${Out}/contigs/contigs_under_${rDNAminCntgLength}.fa $address/OUTPUT/${Out}/contigs/16S_contig.table $address/OUTPUT/${Out}/contigs/16S_contig.m7 $address/OUTPUT/${Out}/contigs/16S_contig.m8 $address/OUTPUT/${Out}/contigs/16S_contig.nm8 $address/OUTPUT/${Out}/contigs/16S_contig.nm7 $address/OUTPUT/${Out}/contigs/16S_contig.tsv $address/OUTPUT/${Out}/contigs/map_nobarcode_contig.uc $address/OUTPUT/${Out}/datestartrenaming.tmp 
	echo "14" > $address/CR.step; CFLR="N"
fi
}

PCA_Maker ()
{
	if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "14" ]]; then echo "******Skipping PCA_Maker analysis"; else
	echo -n "Starting PCA analysis..."
	if [[ "$(grep -m2 -c '>' $address/OUTPUT/${Out}/16S.fasta)" -gt 1 ]]
	then
		python ${PCAmaker} $address/OUTPUT/${Out}/16S.fasta $address/OUTPUT/${Out}/PCA.png > $address/OUTPUT/${Out}/PCAlog.txt
	fi
	if [[ "$(grep -m2 -c '>' $address/OUTPUT/${Out}/contigs/16S_contigs.fasta)" -gt 1 ]]
	then
		python ${PCAmaker} $address/OUTPUT/${Out}/contigs/16S_contigs.fasta $address/OUTPUT/${Out}/contigs/PCA_contig.png > $address/OUTPUT/${Out}/PCAlog.txt
	fi
	rm -rf $address/OUTPUT/${Out}/16S.nonchimera $address/OUTPUT/${Out}/PCAlog.txt
	echo "    Done!"
	echo "15" > $address/CR.step; CFLR="N"
fi
}

rDNA_Abundance ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "15" ]]; then echo "******Skipping Taxon Abundance Finding"; else
	head -1 $address/OUTPUT/${Out}/otu_table.seqidnum > $address/OUTPUT/${Out}/otu_table.header
	sed -i '1,1d' $address/OUTPUT/${Out}/otu_table.seqidnum
	sed -i '1,1d' $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum
	sort -S 50% --parallel=${threads} -k1,1 $address/OUTPUT/${Out}/otu_table.seqidnum > $address/OUTPUT/${Out}/otu_table.ord
	sort -S 50% --parallel=${threads} -k1,1 $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum > $address/OUTPUT/${Out}/contigs/otu_table_contig.ord
	if [[ -s $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_OTU.tsv ]]
	then
		join -1 1 -2 1 -o 1.2,2.2 $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_OTU.tsv $address/OUTPUT/${Out}/otu_table.ord | sed 's/ /\t/g' | sort -S 50% --parallel=${threads} -V > $address/OUTPUT/${Out}/otu_table.otus
		awk '{seen[$1]+=$2}END{ for (id in seen) print id "\t" seen[id] }' $address/OUTPUT/${Out}/otu_table.otus | sort -S 50% --parallel=${threads} -V | awk -v cutoff=${rDNA16Scutoff} '{if ($2>cutoff) print}' > $address/OUTPUT/${Out}/otu_table.sum
		cat $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/otu_table.sum > $address/OUTPUT/${Out}/otu_table.tsv
		sq=`grep -c ">" $address/OUTPUT/${Out}/16S.fasta`
		echo "$sq reads matched to $(cat $address/OUTPUT/${Out}/otu_table.sum | wc -l) OTUs"
		if [[ -s $address/OUTPUT/${Out}/contigs/otu_table_contig.ord ]]
		then
			join -1 1 -2 1 -o 1.2,2.2 $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_OTU.tsv $address/OUTPUT/${Out}/contigs/otu_table_contig.ord | sed 's/ /\t/g' | sort -S 50% --parallel=${threads} -V > $address/OUTPUT/${Out}/contigs/otu_table_contig.otus
			awk '{seen[$1]+=$2}END{ for (id in seen) print id "\t" seen[id] }' $address/OUTPUT/${Out}/contigs/otu_table_contig.otus | sort -S 50% --parallel=${threads} -V | awk -v cutoff=${rDNA16Scutoff} '{if ($2>(cutoff / 10)) print}' > $address/OUTPUT/${Out}/contigs/otu_table_contig.sum
			cntg=`grep -c ">" $address/OUTPUT/${Out}/contigs/16S_contigs.fasta`
			echo "$cntg contigs matched to $(cat $address/OUTPUT/${Out}/contigs/otu_table_contig.sum | wc -l) OTUs"
			cat $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/contigs/otu_table_contig.sum > $address/OUTPUT/${Out}/contigs/otu_table_contig.tsv
		else
			rm -rf $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum $address/OUTPUT/${Out}/contigs/otu_table_contig.ord
		fi
	else
		if [[ -s $address/OUTPUT/${Out}/otu_table.ord ]]
		then
			sq=`grep -c ">" $address/OUTPUT/${Out}/16S.fasta`
			echo "$sq reads matched to $(cat $address/OUTPUT/${Out}/otu_table.ord | wc -l) OTUs"
			cat $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/otu_table.ord > $address/OUTPUT/${Out}/otu_table.tsv
		else
			rm -rf $address/OUTPUT/${Out}/otu_table.ord $address/OUTPUT/${Out}/otu_table.tsv
		fi
		if [[ -s $address/OUTPUT/${Out}/contigs/otu_table_contig.ord ]]
		then
			cntg=`grep -c ">" $address/OUTPUT/${Out}/contigs/16S_contigs.fasta`
			echo "$cntg contigs matched to $(cat $address/OUTPUT/${Out}/contigs/otu_table_contig.ord | wc -l) OTUs"
			cat $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/contigs/otu_table_contig.ord > $address/OUTPUT/${Out}/contigs/otu_table_contig.tsv
		else
			rm -rf $address/OUTPUT/${Out}/contigs/otu_table_contig.ord $address/OUTPUT/${Out}/contigs/otu_table_contig.tsv
		fi
	fi
	rm -rf $address/OUTPUT/${Out}/otu_table.ord $address/OUTPUT/${Out}/otu_table.otus $address/OUTPUT/${Out}/otu_table.header $address/OUTPUT/${Out}/otu_table.sum $address/OUTPUT/${Out}/contigs/otu_table_contig.ord $address/OUTPUT/${Out}/contigs/otu_table_contig.otus $address/OUTPUT/${Out}/contigs/otu_table_contig.sum
	echo "16" > $address/CR.step; CFLR="N"
fi
}

rDNA_TaxonFinding ()
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "16" ]]; then echo "******Skipping Taxon Analysis"; else
	rm -rf $address/OUTPUT/${Out}/taxons $address/OUTPUT/${Out}/contigs/taxons
	if [[ -s $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_Taxon.tsv ]]
	then
		if [[ -s $address/OUTPUT/${Out}/otu_table.tsv ]]
		then
			mkdir $address/OUTPUT/${Out}/taxons
			while read otunum counts; do
				echo "$counts	$(grep -m 1 -w "${otunum#>}" $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_Taxon.tsv)" >> $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre
			done < $address/OUTPUT/${Out}/otu_table.tsv
			sed -i -e 1,1d $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre
		fi
		if [[ -s $address/OUTPUT/${Out}/contigs/otu_table_contig.tsv ]]
		then
			mkdir $address/OUTPUT/${Out}/contigs/taxons
			while read otunum counts; do
				echo "$counts	$(grep -m 1 -w "${otunum#>}" $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_Taxon.tsv)" >> $address/OUTPUT/${Out}/contigs/taxons/otu_table_contig_with_associated_taxons.pre
			done < $address/OUTPUT/${Out}/contigs/otu_table_contig.tsv
			sed -i -e 1,1d $address/OUTPUT/${Out}/contigs/taxons/otu_table_contig_with_associated_taxons.pre
		fi
		echo "Counts	SeqIDnumber	OTU	Kingdom	Phylum	Class	Order	Family	Genus	Species" > $address/OUTPUT/${Out}/taxons/associatedtaxons.header
		if [[ -s $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre ]]
		then
			cat $address/OUTPUT/${Out}/taxons/associatedtaxons.header $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.pre > $address/OUTPUT/${Out}/taxons/otu_table_with_associated_taxons.tsv
		fi
		if [[ -s $address/OUTPUT/${Out}/contigs/taxons/otu_table_contig_with_associated_taxons.pre ]]
		then
			cat $address/OUTPUT/${Out}/taxons/associatedtaxons.header $address/OUTPUT/${Out}/contigs/taxons/otu_table_contig_with_associated_taxons.pre > $address/OUTPUT/${Out}/contigs/taxons/otu_table_contig_with_associated_taxons.tsv
		fi
		echo "Counts	Full_Header" > $address/OUTPUT/${Out}/taxons/counts_fullheader.header
		if [[ -s $address/OUTPUT/${Out}/otu_table.seqidnum ]]
		then
			while read seqidnum counts; do
				echo "$counts	$(grep -w -m 1 ">${seqidnum#>}" $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.headers.txt)" >> $address/OUTPUT/${Out}/taxons/headers_counts.pre
			done < $address/OUTPUT/${Out}/otu_table.seqidnum
			if [[ -s $address/OUTPUT/${Out}/taxons/headers_counts.pre ]]
			then
				cat $address/OUTPUT/${Out}/taxons/counts_fullheader.header $address/OUTPUT/${Out}/taxons/headers_counts.pre > $address/OUTPUT/${Out}/taxons/headers_counts.tsv
			fi
		fi
		if [[ -s $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum ]]
		then
			while read seqidnum counts; do
				echo "$counts	$(grep -w -m 1 ">${seqidnum#>}" $ReferencesFolder/Taxonomy_Index/${SubRef%.*}.headers.txt)" >> $address/OUTPUT/${Out}/contigs/taxons/headers_counts_contig.pre
			done < $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum
			if [[ -s $address/OUTPUT/${Out}/contigs/taxons/headers_counts_contig.pre ]]
			then
				cat $address/OUTPUT/${Out}/taxons/counts_fullheader.header $address/OUTPUT/${Out}/contigs/taxons/headers_counts_contig.pre > $address/OUTPUT/${Out}/contigs/taxons/headers_counts_contig.tsv
			fi
		fi
	else
		echo "Couldn't find the file with taxons on '$ReferencesFolder/Taxonomy_Index/${SubRef%.*}.ID_2_Taxon.tsv'"
	fi
	rm -rf $address/OUTPUT/${Out}/taxons/*.pre $address/OUTPUT/${Out}/taxons/*.header ; rm -rf $address/OUTPUT/${Out}/contigs/taxons/*.pre
	rm -rf $address/OUTPUT/${Out}/otu_table.seqidnum $address/OUTPUT/${Out}/contigs/otu_table_contig.seqidnum
	echo "17" > $address/CR.step; CFLR="N"
fi
}

CleaningTheMess () # Renames files properly, makes final calculations of statistics and assignment of variables, finishes printing the Log.txt
{
if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "12" && `cat $address/CR.step` != "16"  && `cat $address/CR.step` != "18" && `cat $address/CR.step` != "17" ]]; then echo "******Skipping Self Organizing Module"; else
	case $T1 in
		P|p|N|n)
			if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
			then
				${pigz} -p ${threads} $address/OUTPUT/${Out}/hits.fasta
			fi
			if [[ -s $address/OUTPUT/${Out}/hits.fasta.gz ]]
			then
				rm -rf $address/OUTPUT/${Out}/hits.fasta
			fi
			rm -rf $address/OUTPUT/${Out}/tsv_b6.header.txt
			if [[ -s $address/OUTPUT/${Out}/log.tmps ]]
			then
				touch $address/OUTPUT/${Out}/log.tmps
			else
				cat $address/OUTPUT/${Out}/tmp1.* > $address/OUTPUT/${Out}/log.tmps
				sed -i 's/|/\t/g' $address/OUTPUT/${Out}/log.tmps
				rm -rf $address/OUTPUT/${Out}/tmp1.*
			fi
			ls $address/OUTPUT/${Out}/blast_hits/*.ft > $address/OUTPUT/${Out}/blast_hits/list
			if [[ -s $address/OUTPUT/${Out}/blast_hits/list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/blast_hits/list`; do
					mv $file ${file%.ft}.blast_fmt6.tsv
				done
			fi
			rm -rf $address/OUTPUT/${Out}/blast_hits/list
			cntg="0"
			AvgSizeCntg="0"
			TotalSizeCntg="0"
			StdDevCntg="0"
			MaxCntg="0"
			touch $address/OUTPUT/${Out}/contigs/null000.sizes
			cat $address/OUTPUT/${Out}/contigs/*.sizes > $address/OUTPUT/${Out}/contigs/${Out}.allsizes
			if [[ -s $address/OUTPUT/${Out}/contigs/${Out}.allsizes ]]
			then
				cntg=`awk 'BEGIN{s=0;}{s=s+$4;}END{print s;}' $address/OUTPUT/${Out}/log.tmps`
				AvgSizeCntg=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				TotalSizeCntg=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				StdDevCntg=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
				MaxCntg=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/contigs/${Out}.allsizes`
			fi
			echo "$AvgSizeCntg|$TotalSizeCntg|$StdDevCntg|$MaxCntg" > $address/OUTPUT/${Out}/CntgStats.nmb # Calculates all statistics of contigs for the reference as a whole.
			CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
			# echo "Contigs: $cntg ; ContigStats: $CntgStats ; AvgSizeCntg: $AvgSizeCntg ; TotalSizeCntg: $TotalSizeCntg ; StdDevCntg: $StdDevCntg ;MaxCntg: $MaxCntg"
			rm -rf $address/OUTPUT/${Out}/contigs/*.sizes ; rm -rf $address/OUTPUT/${Out}/contigs/${Out}.allsizes
			ls $address/OUTPUT/${Out}/contigs/*.fasta > $address/OUTPUT/${Out}/contigs/contigs_list
			if [[ -s contigs_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/contigs/contigs_list`; do
					mv $file ${file%.ft.fasta}.fasta
				done
			fi
			rm -rf $address/OUTPUT/${Out}/contigs/contigs_list
			ORFs="0"
			AvgSizeORF="0"
			TotalSizeORF="0"
			StdDevORF="0"
			MaxORF="0"
			touch $address/OUTPUT/${Out}/ORFs/null000.sizes
			cat $address/OUTPUT/${Out}/ORFs/*.sizes > $address/OUTPUT/${Out}/ORFs/${Out}.allsizes
			if [[ -s $address/OUTPUT/${Out}/ORFs/${Out}.allsizes ]]
			then
				ORFs=`awk 'BEGIN{s=0;}{s=s+$9;}END{print s;}' $address/OUTPUT/${Out}/log.tmps`
				AvgSizeORF=`awk 'BEGIN{s=0;}{s=s + $1;}END{printf "%.5f", s/NR;}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				TotalSizeORF=`awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				StdDevORF=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { printf "%.5f", sqrt(mean2 / NR); }' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
				MaxORF=`awk 'BEGIN{x=0};$0>x{x=$0};END{print x}' $address/OUTPUT/${Out}/ORFs/${Out}.allsizes`
			fi
			echo "$AvgSizeORF|$TotalSizeORF|$StdDevORF|$MaxORF" > $address/OUTPUT/${Out}/ORFStats.nmb # Calculates all statistics of ORFs for the reference as a whole.
			ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
			# echo "Contigs: $ORFs ; ContigStats: $ORFStats ; AvgSizeORF: $AvgSizeORF ; TotalSizeORF: $TotalSizeORF ; StdDevORF: $StdDevORF ; MaxORF: $MaxORF"
			rm -rf $address/OUTPUT/${Out}/ORFs/*.sizes ; rm -rf $address/OUTPUT/${Out}/ORFs/${Out}.allsizes
			ls $address/OUTPUT/${Out}/ORFs/*.fasta > $address/OUTPUT/${Out}/ORFs/orfs_list
			if [[ -s $address/OUTPUT/${Out}/ORFs/orfs_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/ORFs/orfs_list`; do
					mv $file ${file%.ft.fasta}.fasta
				done
			fi
			rm -rf $address/OUTPUT/${Out}/ORFs/orfs_list
			ls $address/OUTPUT/${Out}/fasta_hits/*.fasta > $address/OUTPUT/${Out}/fasta_hits/hits_list
			if [[ -s $address/OUTPUT/${Out}/fasta_hits/hits_list ]]
			then
				for file in `cat $address/OUTPUT/${Out}/fasta_hits/hits_list`; do
					mv $file ${file%.ft.fasta}.hits.fasta
				done
			fi
			rm -rf $address/OUTPUT/${Out}/fasta_hits/hits_list

			rm -rf $address/OUTPUT/${Out}/log.header $address/OUTPUT/${Out}/log.tmp1; rm -rf $address/OUTPUT/${Out}/cont_log*; rm -rf $address/OUTPUT/${Out}/ORF_log*
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
					else
						CntgStats="0|0|0|0"
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
					else
						ORFStats="0|0|0|0"
					fi
				else
					if [[ -s $address/OUTPUT/${Out}/ORFs.nmb ]]
					then 
						ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
						if [[ -s $address/OUTPUT/${Out}/ORFStats.nmb ]]
						then
							ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
						else
							ORFStats="0|0|0|0"
						fi
					fi
				fi
			fi
			echo "|Contigs: $cntg
|ORFs: $ORFs
-----------------------------------------------------------------------------
Subref_database|Hits_seq.|Ppm|Contigs|Total_Size_Contigs|Avg_Size_Contigs|Standard_Deviation_Contig_Sizes|Max_Contig_Size|ORFs|Total_Size_ORFs|Avg_Size_ORF|Standard_Deviation_ORF_Sizes|Max_ORF_Size|Hits_Blast_Time|SPADES_Time|Contig_Blast_time|ORF_Blast_time|Total_SubRef_Time|Status" > $address/OUTPUT/${Out}/log.header
			cat $address/OUTPUT/${Out}/log.header $address/OUTPUT/${Out}/log.tmps | sed 's/|/\t/g' > $address/OUTPUT/${Out}/log.tmp1
			rm -rf $address/OUTPUT/${Out}/log.header
		;;
		16S|16s|16)
			touch $address/OUTPUT/${Out}/Log.txt
			if [[ -s $address/OUTPUT/${Out}/16S.fasta && -s $address/OUTPUT/${Out}/DiscardedSequences_possible16S.fasta && -s $address/OUTPUT/${Out}/contigs.fasta ]]
			then
				rm -rf $address/OUTPUT/${Out}/hits.fasta $address/OUTPUT/${Out}/hits.fasta.gz
			else
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					${pigz} -p ${threads} $address/OUTPUT/${Out}/hits.fasta
				fi
			fi
			cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
			if [[ -s $address/OUTPUT/${Out}/CntgStats.nmb ]]
			then
				CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
			else
				CntgStats="0|0|0|0"
			fi
			echo "Contigs:	$cntg" >> $address/OUTPUT/${Out}/log.tmp1
			rm -rf $address/OUTPUT/${Out}/hits.fasta
		;;
		G|g)
			if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
			then
				${pigz} -p ${threads} $address/OUTPUT/${Out}/hits.fasta
			fi
			rm -rf $address/OUTPUT/${Out}/hits.fasta

			cntg=`cat $address/OUTPUT/${Out}/cntg.nmb`
			if [[ -s $address/OUTPUT/${Out}/CntgStats.nmb ]]
			then
				CntgStats=`cat $address/OUTPUT/${Out}/CntgStats.nmb`
			else
				CntgStats="0|0|0|0"
			fi
			ORFs=`cat $address/OUTPUT/${Out}/ORFs.nmb`
			if [[ -s $address/OUTPUT/${Out}/ORFStats.nmb ]]
			then
				ORFStats=`cat $address/OUTPUT/${Out}/ORFStats.nmb`
			else
				ORFStats="0|0|0|0"
			fi
			if [[ -d $ReferencesFolder/$SubRef ]]
			then
				if [[ -s $address/OUTPUT/${Out}/log.tmps ]]
				then
					touch $address/OUTPUT/${Out}/log.tmps
				else
					echo "	Contigs: $cntg
	ORFs: $ORFs
AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntg
$(cat $address/OUTPUT/${Out}/CntgStats.nmb)
AvgSizeORF|TotalSizeORF|StdDevORF|MaxORF
$(cat $address/OUTPUT/${Out}/ORFStats.nmb)

-----------------------------------------------------------------------------
Subref_database|Hits_seq.|Ppm|Contigs|Total_Size_Contigs|Avg_Size_Contigs|Standard_Deviation_Contig_Sizes|Max_Contig_Size|ORFs|Total_Size_ORFs|Avg_Size_ORF|Standard_Deviation_ORF_Sizes|Max_ORF_Size
$(cat $address/OUTPUT/${Out}/tmp1.subs)" | sed 's/|/\t/g' > $address/OUTPUT/${Out}/log.tmp1
					rm -rf $address/OUTPUT/${Out}/log.header
				fi
			else
				echo "	Contigs: $cntg
	ORFs: $ORFs
AvgSizeCntg|TotalSizeCntg|StdDevCntg|MaxCntg
$(cat $address/OUTPUT/${Out}/CntgStats.nmb)
AvgSizeORF|TotalSizeORF|StdDevORF|MaxORF
$(cat $address/OUTPUT/${Out}/ORFStats.nmb)" | sed 's/|/\t/g' > $address/OUTPUT/${Out}/log.tmp1
			fi
			rm -rf $address/OUTPUT/${Out}/cntg.header $address/OUTPUT/${Out}/ORF.header
		;;
	esac
	rm -rf $address/OUTPUT/${Out}/assembly_*
	rm -rf $address/OUTPUT/${Out}/fulltime.tmp
	touch $address/OUTPUT/${Out}/spadestime.nmb
	touch $address/OUTPUT/${Out}/spadestime2.nmb
	if [[ $(cat $address/OUTPUT/${Out}/spadestime2.nmb | wc -l) -eq 1 ]]
	then
		echo "Time for first SPADES run (failed): $(cat $address/OUTPUT/${Out}/spadestime.nmb | awk '{ sum += $1 } END { print sum }')s
Time for second SPADES run: $(cat $address/OUTPUT/${Out}/spades2time.nmb | awk '{ sum += $1 } END { print sum }')s" >> $address/OUTPUT/${Out}/fulltime.tmp
	else
		if [[ -s $address/OUTPUT/${Out}/spadestime.nmb ]]
		then
			echo "Time for the SPADES run: $(cat $address/OUTPUT/${Out}/spadestime.nmb $address/OUTPUT/${Out}/spadestime2.nmb | awk '{ sum += $1 } END { print sum }')s" >> $address/OUTPUT/${Out}/fulltime.tmp
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/quasttime.nmb ]]
	then
		echo "Time for Quast assessment: $(cat $address/OUTPUT/${Out}/quasttime.nmb)s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/orffindertime.nmb ]]
	then
		echo "Time for ORFs to be found by ORFfinder: $(cat $address/OUTPUT/${Out}/orffindertime.nmb | awk '{ sum += $1 } END { print sum }')s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/contigfilteringtime.nmb ]]
	then
		echo "Time for Contigs to be filtered against databases: $(cat $address/OUTPUT/${Out}/contigfilteringtime.nmb | awk '{ sum += $1 } END { print sum }')s" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	if [[ -s $address/OUTPUT/${Out}/orffilteringtime.nmb ]]
	then
		echo "Time for ORFs to be filtered against databases: $(cat $address/OUTPUT/${Out}/orffilteringtime.nmb | awk '{ sum += $1 } END { print sum }')s" >> $address/OUTPUT/${Out}/fulltime.tmp
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
		d3="at least $(echo "$(date -u +%s) - $d1" |bc -l)"
	else
		d3=$(echo "$(date -u +%s) - $d1" |bc -l)
	fi
	BucketTotalTime=$(cat $address/${BucketsFolder}/BucketTotalTime.nmb)
	if [[ -s $address/${BucketsFolder}/bucketpreviouslygeneratedmessage.tmp ]]
	then
		echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds
	Total time adding time to trim, analyse sequence quality and generate buckets (which were all done in a previous run): : $(echo $d3 + $BucketTotalTime | bc -l) seconds" >> $address/OUTPUT/${Out}/fulltime.tmp
	else
		echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds
	Total time excluding time to trim, analyse sequence quality and generate buckets: $(echo $d3 - $BucketTotalTime | bc -l) seconds" >> $address/OUTPUT/${Out}/fulltime.tmp
	fi
	echo "

" >> $address/OUTPUT/${Out}/log.tmp1
	echo "

" >> $address/OUTPUT/${Out}/fulltime.tmp
	echo "

" >> $address/OUTPUT/${Out}/Quastlog.tmp
	cat $address/OUTPUT/${Out}/fulltime.tmp $address/OUTPUT/${Out}/log.tmp1 $address/OUTPUT/${Out}/Quastlog.tmp >> $address/OUTPUT/${Out}/Log.txt
	if [[ -s $address/OUTPUT/${Out}/contigs.fasta ]]
	then 
		${pigz} -p ${threads} $address/OUTPUT/${Out}/contigs.fasta
	fi
	if [[ -s $address/OUTPUT/${Out}/ORFs.fa ]]
	then
		${pigz} -p ${threads} $address/OUTPUT/${Out}/ORFs.fa
	else
		if [[ -s $address/OUTPUT/${Out}/ORFs.${Out}.fasta ]]
		then
			mv $address/OUTPUT/${Out}/ORFs.${Out}.fasta $address/OUTPUT/${Out}/ORFs.fa
			${pigz} -p ${threads} $address/OUTPUT/${Out}/ORFs.fa
		else
			rm -rf $address/OUTPUT/${Out}/ORFs.${Out}.fasta $address/OUTPUT/${Out}/ORFs.fa
		fi
	fi
	case $T1 in 
		G|g)
			rm -rf $address/OUTPUT/${Out}/log.tmp1
		;;
		P|p|N|n)
			mv $address/OUTPUT/${Out}/log.tmp1 $address/OUTPUT/${Out}/SubRefs.tsv
		;;
		16S|16s|16)
			if [[ -s $address/OUTPUT/${Out}/contigs.fasta.gz ]]
			then
				mv $address/OUTPUT/${Out}/contigs.fasta.gz $address/OUTPUT/${Out}/assembled16S.fasta.gz
			fi
			rm -rf $address/OUTPUT/${Out}/log.tmp1
		;;
	esac
	rm -rf $address/OUTPUT/${Out}/contigs.fasta $address/OUTPUT/${Out}/log.tmp $address/OUTPUT/${Out}/fulltime.tmp $address/OUTPUT/${Out}/Quastlog.tmp $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/*.time2; rm -rf $address/OUTPUT/${Out}/*.rev; rm -rf $address/OUTPUT/${Out}/*.hits ; rm -rf $address/OUTPUT/${Out}/*time.nmb ; rm -rf $address/OUTPUT/${Out}/*time.tmp $address/OUTPUT/${Out}/log.tmps
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
		else
			CntgStats="0|0|0|0"
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
			else
				ORFStats="0|0|0|0"
			fi
		fi
	fi
	if [[ -s $address/OUTPUT/${Out}/ppm1.nmb ]]
	then
		ppm1=`cat $address/OUTPUT/${Out}/ppm1.nmb`
	fi
	if [[ -s $address/${BucketsFolder}/reads.nmb ]]
	then
		reads=`cat $address/${BucketsFolder}/reads.nmb`
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
	if [[ ${d3days} -ge 1 ]]
	then
		echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${d3days}d${d3hours}h${d3minutes}m${d3seconds}s). ] \n"
	else
		if [[ ${d3hours} -ge 1 ]]
		then
			echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${d3hours}h${d3minutes}m${d3seconds}s). ] \n"
		else
			if [[ ${d3minutes} -ge 1 ]]
			then
				echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds (${d3minutes}m${d3seconds}s). ] \n"
			else
				echo -e "\n [ BEAF12.01.17 worked in $R1 with reference as $Ref (output as ${Out}) for $d3 seconds. ] \n"
			fi
		fi
	fi
	sed -i -e 1,1d $address/config.kp
	case $Keep in
		Y|y)
			touch $address/${BucketsFolder}
		;;
		*)
			rm -rf $address/${BucketsFolder}
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
	ls $address/OUTPUT > $address/OUTPUT/list
	sed -i '/list/d' $address/OUTPUT/list
	sed -i '/Errors/d' $address/OUTPUT/list
	sed -i '/Errors.old/d' $address/OUTPUT/list
	sed -i '/OUTPUT/d' $address/OUTPUT/list
	sed -i '/Previous_Runs/d' $address/OUTPUT/list
	sed -i '/redolist/d' $address/OUTPUT/list
	if [[ -d $address/OUTPUT/Errors ]]
	then
		if [[ -d $address/OUTPUT/Errors.old ]]
		then
			mv $address/OUTPUT/Errors.old $address/OUTPUT/Errors
		fi
		mv $address/OUTPUT/Errors $address/OUTPUT/Errors.old
		mkdir $address/OUTPUT/Errors
		mv $address/OUTPUT/Errors.old $address/OUTPUT/Errors
	else
		mkdir $address/OUTPUT/Errors
	fi
	for folder in `cat $address/OUTPUT/list`; do
	        if [[ -s $address/OUTPUT/$folder/contigs.fasta.gz ]]
		then
			touch $address/OUTPUT/$folder/contigs.fasta.gz
	        else
			if [[ -s $address/OUTPUT/$folder/16S.fasta && -s $address/OUTPUT/$folder/hits.fasta.gz && -s $address/OUTPUT/$folder/otu_b6out.tsv && -s $address/OUTPUT/$folder/otu_table.tsv && -s $address/OUTPUT/$folder/PCA.png && -s $address/OUTPUT/$folder/Log.txt && -d $address/OUTPUT/$folder/taxons ]]
			then
				touch $address/OUTPUT/${Out}/Log.txt
			else
			        if [[ -d $address/OUTPUT/$folder/contigs ]]
			        then
					if [[ "$(cat $address/OUTPUT/$folder/contigs/taxons/otu_table_contig_with_associated_taxons.tsv | wc -l)" -gt 1 ]]
					then
						touch $address/OUTPUT/$folder/contigs/taxons
					else
					        if [[ $(cat $address/OUTPUT/$folder/Log.txt | grep "Contigs: " | sed 's/Contigs: //' | sed 's/ //' | sed 's/\t//' | sed 's/	//' | sed '/[a-z]/d' | sed '/[A-Z]/d') -gt 0 ]]
						then
							touch $address/OUTPUT/$folder/contigs
					        else
							# rm -rf $address/OUTPUT/$folder/contigs
							subrefhits=0
							subrefhits=$(cat $address/OUTPUT/$folder/blast_hits/*.tsv | wc -l)
							if [[ ${subrefhits} -gt 0 ]]
							then
								echo "BEAF couldn't find contigs for ${folder}, although a total of ${subrefhits} hits were found for it's subreferences."
							else
								mv $address/OUTPUT/$folder $address/OUTPUT/Errors
							fi
					        fi
						if [[ -d $address/OUTPUT/$folder/ORFs ]]
						then
							if [[ $(cat $address/OUTPUT/$folder/Log.txt | grep "ORFs: " | sed 's/ORFs: //' | sed 's/ //' | sed 's/\t//' | sed 's/	//' | sed '/[a-z]/d' | sed '/[A-Z]/d') -gt 0 ]]
							then
								touch $address/OUTPUT/$folder/ORFs/
							else
								rm -rf $address/OUTPUT/$folder/ORFs/
							fi
						fi
					fi
	                	else
					mv $address/OUTPUT/$folder $address/OUTPUT/Errors
				fi
	                fi
	        fi
		rm -rf $address/OUTPUT/Errors/$folder/test.txt
	done
	rm -rf $address/OUTPUT/list
	ls $address/OUTPUT/Errors/ > $address/OUTPUT/Errors/redolist
	sed -i '/list/d' $address/OUTPUT/Errors/redolist
	sed -i '/Errors/d' $address/OUTPUT/Errors/redolist
	sed -i '/Errors.old/d' $address/OUTPUT/Errors/redolist
	sed -i '/OUTPUT/d' $address/OUTPUT/Errors/redolist
	sed -i '/Previous_Runs/d' $address/OUTPUT/Errors/redolist
	sed -i '/redolist/d' $address/OUTPUT/Errors/redolist
	for out in `cat $address/OUTPUT/Errors/redolist`; do
		grep "$out" $ConfigFile >> $address/OUTPUT/Errors/config.redo
	done
	if [[ -s $address/OUTPUT/Errors/redolist ]]
	then
		ErrorNumber=`cat $address/OUTPUT/Errors/redolist | wc -l`
		echo "Error Revision complete: BEAF could not find proper results for $ErrorNumber references. These files can either have too few hits or may have suffered errors during the pipeline process, and are now stored in the folder 'Errors' inside the OUTPUT folder. You can find a config file in the same folder that can be used to rerun the process only for these specific files. It is recommended you check these files first."
	else
		echo "Error Revision complete: BEAF found results for every reference it worked on."
		rm -rf $address/OUTPUT/Errors/redolist
		if [[ -d $address/OUTPUT/Errors/Errors.old ]]
		then
			mv $address/OUTPUT/Errors/Errors.old $address/OUTPUT/
			rm -rf $address/OUTPUT/Errors
			mv $address/OUTPUT/Errors/Errors.old $address/OUTPUT/Errors
		else
			if [[ "$(ls $address/OUTPUT/Errors/ | wc -l)" -gt 0 ]]
			then
				touch $address/OUTPUT/Errors
			else
				rm -rf $address/OUTPUT/Errors
			fi
		fi
	fi
	rm -rf $address/OUTPUT/list
}

	# ======================================================================================================================================================================================== #
	# =====================================================================================PIPELINE BEAF====================================================================================== #
	# ======================================================================================================================================================================================== #


BEAF () # The pipeline BEAF will run when using full version.
{
	echo "
##### Running BEAF, full version #####

Config file = ${ConfigFile}
Output folder = ${address}
References folder = ${ReferencesFolder}

Threads = ${threads}
"
	d0=`date -u "+%s"`
	make_kp
	Check
	TimeHeader
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		d1=`date -u "+%s"`
		case $T1 in
			G|g)
				if [[ -s $ReferencesFolder/$Ref ]]
				then 
					Ref="$ReferencesFolder/$Ref"
				else
					if [[ ! -s $Ref ]]
					then
						if [[ ! -d ${ReferencesFolder}/$Subref ]]
						then
							echo "Couldn't find your Reference $Ref (neither as full path nor in the References_Folder)"
							exit exit exit
						fi
					fi
				fi
				if [[ -d ${ReferencesFolder}/${SubRef} ]]
				then
					if [[ -s $Ref ]]
					then
						echo -e "\n# Starting work in file $R1 with $Ref as reference and genomes in $ReferencesFolder/$SubRef as subreferences, at $(date '+%X %e/%m/%Y'), going to ${Out}\n"
					else
						echo -e "\n# Starting work in file $R1 with a pooled reference made from genomes in $ReferencesFolder/$SubRef (subreferences) at $(date '+%X %e/%m/%Y'), going to ${Out}\n"
					fi
				else
					echo -e "\n# Starting work in file $R1 with $Ref as reference, at $(date '+%X %e/%m/%Y'), going to ${Out}\n"
				fi
 				CommandLine="$0 --config=$ConfigFile --output=$address --threads=$threads --references $ReferencesFolder --maxaccepts2=$maxaccepts2 --maxrejects1=$maxrejects1 --maxrejects2=$maxrejects2 --gen_id1=$GenID1 --gen_id2=$GenID2 --gen_evalue1=$GENevalue1 --gen_evalue2=$GENevalue2 --GenomeSplitMethod=$GenomeSplitMethod --GenCoverage1=$GenomeCoverage1 --GenCoverage2=$GenomeCoverage2 --GenFragLength=$GenomeFragLength" 
			;;
			P|p|N|n)
				if [[ -s $ReferencesFolder/$Ref ]]
				then 
					Ref="$ReferencesFolder/$Ref"
				else
					if [[ ! -s $Ref ]]
					then
						echo "Couldn't find your Reference $Ref (neither as full path nor in the References_Folder)"
						exit exit exit
					fi
				fi
				echo -e "\n# Starting work in file $R1 with $Ref as reference for the first filter and databases in folder ${SubRef} as subreferences, at $(date '+%X %e/%m/%Y'), going to ${Out}\n"
				CommandLine="$0 --config=$ConfigFile --output=$address --threads=$threads --references $ReferencesFolder --maxaccepts2=$maxaccepts2 --maxrejects1=$maxrejects1 --maxrejects2=$maxrejects2 --prot_id1=$ProtID1 --prot_id2=$ProtID2 --prot_evalue1=$PROTevalue1 --prot_evalue2=$PROTevalue2 --prot_qcov1=$PROTqcovHits --prot_qcov2=$PROTqcovORFs"
			;;
			16S|16s|16)
				if [[ ! -s $ReferencesFolder/$Ref ]]
				then
					if [[ ! -s $Ref && -s $ReferencesFolder/nr97_Greengenes.udb ]]
					then
						echo "Using a clustered Green Genes database (clustered to 97% identity) to search for OTUs."
						Ref=nr97_Greengenes.udb
					else
						echo "Couldn't find reference $Ref in References Folder $ReferencesFolder."
					fi
				fi

				if [[ ! -s $ReferencesFolder/$SubRef ]]
				then
						if [[ ! -s $SubRef && -s $ReferencesFolder/nr100_Greengenes.udb ]]
						then
							echo "Using the full Green Genes database (clustered to 100% identity) to search for OTUs."
							SubRef=nr100_Greengenes.udb
						else
							echo "Couldn't find subreference $SubRef in References Folder $ReferencesFolder to search for OTUs."
						fi
				fi
				echo -e "\n# Starting work in file $R1 with $Ref as reference for the first filter and ${SubRef} as reference for the second filter, at $(date '+%X %e/%m/%Y'), going to ${Out}\n"
 				CommandLine="$0 --config=$ConfigFile --output=$address --threads=$threads --references $ReferencesFolder --maxaccepts2=$maxaccepts2 --maxrejects1=$maxrejects1 --maxrejects2=$maxrejects2 --16S_id1=$rDNA16SID1 --16S_id2=$rDNA16SID2 --16S_evalue1=$rDNAevalue1 --16S_evalue2=$rDNAevalue2 --16S_qcov=$rDNA16Sqcov --16S_cutoff=$rDNA16Scutoff --16S_min_cntg_length=$rDNAminCntgLength"
			;;
		esac
		BucketsFolder="Bucket_$(basename $R1 | sed 's/.gz//' | sed 's/\..*//')"
		if [[ -d $address/${BucketsFolder} ]]
		then
			ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/buckets_list.txt
			if [[ -s $address/${BucketsFolder}/buckets_list.txt ]];
			then
				echo "# Using previous buckets"
			else
				Trim
				QAnConversion
			fi
		else
			mkdir $address/${BucketsFolder}
			Trim
			QAnConversion
		fi
		if [[ -d $address/OUTPUT/${Out} ]]
		then
			if [[ $CFLR == "N" ]]
			then
				echo "@ Folder ${Out} already exists. Moving it to folder 'Previous_Runs' inside the OUTPUT folder."
				if [[ -d $address/OUTPUT/Previous_Runs ]]
				then
					touch $address/OUTPUT/Previous_Runs
				else
					mkdir $address/OUTPUT/Previous_Runs
				fi
				if [[ -d $address/OUTPUT/Previous_Runs/Previous_${Out} ]]
				then
					mv $address/OUTPUT/Previous_Runs $address/OUTPUT/Previous_Runs.old
					mkdir $address/OUTPUT/Previous_Runs
					mv $address/OUTPUT/Previous_Runs.old $address/OUTPUT/Previous_Runs
				fi
				mv   $address/OUTPUT/${Out} $address/OUTPUT/Previous_Runs/Previous_${Out}
				mkdir $address/OUTPUT/${Out}
				echo "A new folder ${Out} has been created in OUTPUT"
			else
				echo "Folder ${Out} already exists, continuing work..."
			fi
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
		buckets=`ls $address/${BucketsFolder}/*.bk | wc -l`
		Filter1
		PreLogGen
		case $T1 in
			G|g)
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					SubGenomesMakeDB
					SubGenomesFilter2
					G_SPADES1
					G_SPADES2
					GA
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/OUTPUT/${Out}/hits.fasta
					echo "12" > $address/CR.step
				fi
			;;
			P|p|N|n)
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					Extraction
					if [[ -s $address/${BucketsFolder}/reads.nmb ]]
					then
						reads=`cat $address/${BucketsFolder}/reads.nmb`
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
					echo "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs" > $address/OUTPUT/${Out}/tsv_b6.header.txt
					for File in `cat $address/OUTPUT/${Out}/list`; do
						echo "# Working on file ${File%.ft} for ${Out}..."
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
							echo "# ${File%.ft} in ${Out} did not reach the minimum criteria to be considered homologus"
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
							TTime="at least $(echo $d5 - $d4 |bc -l)"
						else
							TTime=$(echo $d5 - $d4 |bc -l)
						fi
						STime=$(echo "$(tail -n 1 $address/OUTPUT/${Out}/spades2time.nmb) + $(tail -n 1 $address/OUTPUT/${Out}/spadestime.nmb)" | bc -l)
						rm -rf $address/OUTPUT/${Out}/$File.time2
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
						echo "${File%.f*}|$sq|$ppm2|$cntg|$TotalSizeCntg|$AvgSizeCntg|$StdDevCntg|$MaxCntg|$ORFs|$TotalSizeORF|$AvgSizeORF|$StdDevORF|$MaxORF|$BTime|$STime|$(cat $address/OUTPUT/${Out}/contigfilteringtime.nmb | awk '{ sum += $1 } END { print sum }')|$(cat $address/OUTPUT/${Out}/orffilteringtime.nmb | awk '{ sum += $1 } END { print sum }')|$TTime|$Warnings" > $address/OUTPUT/${Out}/tmp1.$File
						sed -i -e 1,1d $address/OUTPUT/${Out}/list
					done
					SaveDBs
					cntg="0"
					ORFs="0"
					rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb $address/OUTPUT/${Out}/ORFs.nmb
				else
					rm -rf $address/OUTPUT/${Out}/hits.fasta
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
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					if [[ "$CFLR" == "Y" && `cat $address/CR.step` != "8" ]]; then echo ""; else
						echo "9" > $address/CR.step; CFLR="N"
					fi
					G_SPADES1
					G_SPADES2
					rDNA_Contig
					rDNA_Chimeras
					rDNA_Filter2
					PCA_Maker
					rDNA_Abundance
					rDNA_TaxonFinding
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/OUTPUT/${Out}/hits.fasta
					echo "17" > $address/CR.step
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
		echo "BEAF1011.65 worked for $dtotal seconds (${totaldays}d${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date "+%X %e/%m/%Y")."
	else
		if [[ $totalhours -ge 1 ]]
		then
			echo "BEAF1011.65 worked for $dtotal seconds (${totalhours}h${totalminutes}m${totalseconds}s), ending at $(date "+%X %e/%m/%Y")."
		else
			if [[ $totalminutes -ge 1 ]]
			then
				echo "BEAF1011.65 worked for $dtotal seconds (${totalminutes}m${totalseconds}s), ending at $(date "+%X %e/%m/%Y")."
			else
				echo "BEAF1011.65 worked for $dtotal seconds (${totalseconds}s), ending at $(date "+%X %e/%m/%Y")."
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
		if [[ -d $ReferencesFolder/$SubRef ]]
		then
			echo -e "\n# Starting work in file $R1 with $Ref as reference and $SubRef as subreference, going to ${Out}\n"
		else
			echo -e "\n# Starting work in file $R1 with $Ref as reference, going to ${Out}\n"
		fi
		BucketsFolder="Bucket_$(basename $R1 | sed 's/.gz//' | sed 's/\..*//')"
		if [[ -d $address/${BucketsFolder} ]]
		then
			ls $address/${BucketsFolder}/*.bk > $address/${BucketsFolder}/buckets_list.txt
			if [[ -s $address/${BucketsFolder}/buckets_list.txt ]]
			then
				echo "# Using previous buckets"
			else
				CopyFile
				SoftMergeRename
			fi
		else
			mkdir $address/${BucketsFolder}
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
		buckets=`ls $address/${BucketsFolder}/*.bk | wc -l`
		Filter1
		PreLogGen
		case $T1 in
			G|g)
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					G_SPADES2
					SoftGA
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf $address/OUTPUT/${Out}/hits.fasta
					echo "12" > $address/CR.step
				fi
			;;
			P|p|N|n)
				if [[ -s $address/OUTPUT/${Out}/hits.fasta ]]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					Extraction
					if [[ -s $address/${BucketsFolder}/reads.nmb ]]
					then
						reads=`cat $address/${BucketsFolder}/reads.nmb`
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
							TTime="at least $(echo $d5 - $d4 |bc -l)"
						else
							TTime=$(echo $d5 - $d4 |bc -l)
						fi
						rm -rf $address/OUTPUT/${Out}/$File.time2
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
						echo "${File%.f*}|$sq|$ppm2|$cntg|$TotalSizeCntg|$AvgSizeCntg|$StdDevCntg|$MaxCntg|$(cat $address/OUTPUT/${Out}/contigfilteringtime.nmb | awk '{ sum += $1 } END { print sum }')|NA|NA|NA|NA|NA|NA|$BTime|$STime|$TTime|$Warnings" > $address/OUTPUT/${Out}/tmp1.$File
						sed -i -e 1,1d $address/OUTPUT/${Out}/list
					done
					SaveDBs
					rm -rf $address/OUTPUT/${Out}/ORFs; rm -rf $address/OUTPUT/${Out}/list; rm -rf $address/OUTPUT/${Out}/ppm2.nmb $address/OUTPUT/${Out}/cntg.nmb $address/OUTPUT/${Out}/sq.nmb
				else
					rm -rf $address/OUTPUT/${Out}/hits.fasta
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
	echo "
###########################################################
BEAF1011.65 worked for $dtotal seconds, ending at $(date "+%X %e/%m/%Y").
###########################################################
"
}

	# ======================================================================================================================================================================================== #
	# ========================================================================================MAIN============================================================================================ #
	# ======================================================================================================================================================================================== #

Main () # Interpreting user's input of whether he'll be using soft or full version.
{
	address=${address%/}
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
		if [[ -s $address/CR.step ]] || [[ -s $address/config.kp ]] || [[ -s $address/CR.mode ]] || [[ -d $address/${BucketsFolder} ]]
		then
			echo "Forced restart. Deleting old files and starting a new run."
		fi
		rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/${BucketsFolder}
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
					rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/${BucketsFolder}
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
				read CFLR # Should I Continue From Last Run?
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
						rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/${BucketsFolder}
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
				rm -rf $address/CR.step $address/CR.mode; rm -rf $address/*.kp; rm -rf $address/*.nmb; rm -rf $address/*.tmp; rm -rf $address/${BucketsFolder}
				Main $ver
			fi
		fi
	fi
fi
