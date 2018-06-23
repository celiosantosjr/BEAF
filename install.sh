address=$(cd "$(dirname "")" && pwd)/$(basename "")
LIB=Lib
ReferencesFolder=Reference_seqs
usearch="usearch" ## Path/To/USEARCH
cdhit="cd-hit" ## Path/To/CD-Hit

UpgradeBasics=0
InstallLib=0
InstallORFFinder=0
InstalCDHit=0
Installcutadapt=0
InstallFastQC=0
Installpigz=0
InstallQUAST=0
InstallSPADES=0
InstallPyFasta=0

DownloadGreengenes=0
MakeGreengenesTaxonomyIndex=0
MakeTaxonomyIndex_ggpattern=01


# spades="/Path/To/SPADES"
# quast="/Path/To/QUAST"
# usearch="/Path/To/USEARCH"
# fastqc="/Path/To/FastQC"
# cdhit="/Path/To/cd-hit"
# cutadapt="/Path/To/cutadapt"
# pyfasta="/Path/To/PyFasta"
# pigz="/Path/To/pigz"
# orffinder="/Path/To/bb.orffinder.pl"

# LIB="/Path/To/Lib/"
# extpy="$LIB/ext.py"
# SeqLength="$LIB/SeqLength.py"
# PCAmaker="$LIB/PCA_maker.py"


while [[ $# -gt 0 ]]
do
	case $1 in
		help|?|--help|-h)
			echo "Choose your install options (Default in brackets [default])

Basic path options:
--USEARCH=(Path/To/usearch)	Use the path to USEARCH in your system. [usearch]
--CD-Hit=(Path/To/cd-hit)	Use the path to CD-Hit in your system [cd-hit]
--address=(Path)	Choose where BEAF will be installed [current directory]
--LIB=(Path)	Choose where BEAF sub-scripts and third-party programs will be installed [Lib]
--Refs=(Path)	Path to your References Folder - or where to create one [Reference_seqs]
--Find	This script will search for each of the programs used in BEAF on your computer, and print their paths. This is useful if you're not sure whether you have the softwares installed or not, or if you don't remember their paths. This might take a while.
--FindGreengenes	Same as Find, except that instead of looking for software, it will look for files related to the Greengenes database. This might take a while.

Install options (choose what to install):
--Full	Installs everything, including programs and Greengenes database
--Software	Installs all BEAF sub-scripts and third party software
--UpgradeBasic	Upgrades basic libraries and softwares used by BEAF
--InstallLib	Installs only BEAF sub-scripts and bb.ORFFinder in LIB
--InstallORFFinder	Installs only bb.ORFFinder
--InstallCDHit	Installs CD-Hit
--Installcutadapt	Installs only cutadapt
--InstallFastQC	Installs only FastQC
--Installpigz	Installs only pigz
--InstallQUAST	Installs only QUAST
--InstallSPADES	Installs only SPADES
--InstallPyFasta	Installs only PyFasta

Database options (choose what to download as databases):
--DownloadGreenGenes	Downloads Greengenes database and clusters it using cd-hit (to 100% and 97% identity), then makes a USEARCH database (.udb) file for each clustered file and indexes them for BEAF Taxon Finding process
--MakeGreengenesTaxonomyIndex	Skips the download and clustering steps. This just looks for the files (ReferencesFolder)/Taxonomy/nr100_Greengenes.fasta and (ReferencesFolder)/Taxonomy/nr97_Greengenes.fasta - if you don't have these files already clustered in your computer, inside of a directory named Taxonomy in your references directory, it will fail to index.
--MakeTaxonomyIndex_ggpattern	Indexes the taxonomy of any fasta files inside (ReferencesFolder)/Taxonomy following Greengenes taxonomy. Notice that Greengenes taxonomy follows a pattern, where headers are the following:
	>[Sequence_ID] [Name_or_description_of_the_sequence] k__[Kingdom]; p__[Phylum]; c__[Class]; o__[Order]; f__[Family]; g__Genus; s__[Species]; otu_[OTU_ID]
	If you're using a file with other taxonomy pattern, consider manually converting it to Greengenes pattern of headers before using this tool; indexing it yourself by following the instructions in --help_taxonomy_index; or simply not using the Taxon finding tool of BEAF in your analyses."
		;;
		--help_taxonomy_index)
		echo "If you're using a database other than Greengenes, and with a pattern for headers different than Greengenes (below), you should consider making a Taxonomy index so you can use BEAF's taxon finding tool properly.
	>[Sequence_ID] [Name_or_description_of_the_sequence] k__[Kingdom]; p__[Phylum]; c__[Class]; o__[Order]; f__[Family]; g__[Genus]; s__[Species]; otu_[OTU_ID]

Many databases will use another pattern for headers, and although BEAF will be able to use headers lacking either the name/description or/and the OTU_ID, the automatic indexing won't be able to index other patterns, failing at generating files ID_2_OTU.tsv and ID_2_Taxon.tsv.
In order to make those files, all that is needed is:

ID_2_OTU.tsv
	[Sequence_ID] in the first column of the file ; OTU_ID in the second column.
		Generally, the OTU_ID can be created by using a specific String for each different genus/species/strain/OTU group in the database. For many databases, this could be done by using sed, awk and grep to find sequences and give them similar names. For example, in database where every sequence is labeled as:
			>[Sequence_ID];[Genus];[Species]
			>[Sequence_ID];[Kingdom];[Phylum];[Class];[Order];[Family];[Genus];[Species]
				One could simply use awk or sed to replace the first instance of /;/ for a tab delimiter (\\t or /	/), and the rest would already be considered a proper taxon, then remove the '>' at the beginning of the header.
		For a database where headers consisted of the following:
			>[Sequence_ID] [random_text] [Species] [random_text]
				One could delete everything between the first and second spaces, then after the last space, and finally be left only with the sequence ID and species, removing '>'.
		Of course, many databases lack proper formatting or may prove difficult to parse and properly index - one of the many reasons Greengenes was selected as our default database. In those cases, one can simply consider each sequence as a OTU of it's own, and not use the 'clustering/grouping' of sequences from same Taxon into a single 'Taxon_otu' that BEAF usually does.

ID_2_Taxon.tsv
	[Sequence_ID] in the first column of the file ; OTU_ID in the second column ; Taxonomic groups in each column that follows, as in
	[Sequence_ID]	[OTU_ID]	[Kingdom]	[Phylum]	[Class]	[Order]	[Family]	[Genus]	[Species]
		Again, there are many ways this could be done, but it will depend entirely on your database of choice."
		;;
		--USEARCH|--Usearch|--usearch)
			usearch="${2}"
		;;
		--USEARCH=*|--Usearch=*|--usearch=*)
			usearch="${1#*=}"
		;;
		--CDHIT|--CD-Hit|--cdhit|--cd-hit|--CD-hit)
			cdhit="$2"
		;;
		--CDHIT=*|--CD-Hit=*|--cdhit=*|--cd-hit=*|--CD-hit=*)
			cdhit="${1#*=}"
		;;
		--address|--Address|--ADDRESS|-address|-Address|-ADDRESS|-a|-A)
			if [[ -d "${2}" ]]
			then
				touch ${2}
			else
				mkdir ${2}
			fi
			address=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
		;;
		--address=*|--Address=*|--ADDRESS=*|-address=*|-Address=*|-ADDRESS=*|-a=*|-A=*)
			if [[ -d "${1#*=}" ]]
			then
				touch ${1#*=}
			else
				mkdir ${1#*=}
			fi
			address=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
		;;
		--lib|--Lib|--LIB|-lib|-Lib|-LIB|--library|--Library|--LIBRARY|-Library|-library|-LIBRARY|-L|-l)
			if [[ -d "${2}" ]]
			then
				touch ${2}
			else
				mkdir ${2}
			fi
			LIB=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
		;;
		--lib|--Lib|--LIB|-lib|-Lib|-LIB|--library|--Library|--LIBRARY|-Library|-library|-LIBRARY|-L|-l)
			if [[ -d "${1#*=}" ]]
			then
				touch ${1#*=}
			else
				mkdir ${1#*=}
			fi
			LIB=$(cd "$(dirname "${1#*=}")" && pwd)/$(basename "${1#*=}")
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
		--Find|--find|--FIND|-Find|-find|-FIND|--FindSoftware|--FindSoftwares|--FindProgram|--FindPrograms|--find|--findsoftware|--findsoftwares|--findprogram|--findprograms|--FIND|-Find|-find|-FIND)
			echo "
	## Searching for ORFFinder ##"
			find / -type f -name bb.orffinder.pl 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			find / -type f -name orffinder.pl 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for CD-Hit ##"
			find / -type f -name cd-hit 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			find / -type f -name cdhit 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for cutadapt ##"
			find / -type f -name cutadapt 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for FastQC ##"
			find / -type f -name fastqc 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			find / -type f -name FastQC 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			find / -type f -name FASTQC 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for pigz ##"
			find / -type f -name pigz 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for QUAST ##"
			find / -type f -name metaquast.py 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for SPADES ##"
			find / -type f -name spades.py 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			echo "
	## Searching for PyFasta ##"
			find / -type f -name pyfasta 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
		;;
		--FindGreengenes|--FindGreenGenes|--Findgreengenes|--findgreengenes|--FINDGREENGENES|-FindGreengenes|-Findgreengenes|-FindGreenGenes|-findgreengenes|-FINDGREENGENES)
			echo "
	## Searching for Greengenes ##"
			find / -type f -name *Greengenes* 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
			find / -type f -name *GreenGenes* 2>&1 | grep -v "Permission denied" | grep -v "No such file or directory"
		;;
		--UpgradeBasics)
			UpgradeBasics=1
		;;
		--InstallLib)
			InstallLib=1
		;;
		--InstallORFFinder)
			InstallORFFinder=1
		;;
		--InstallCDHit)
			InstalCDHit=1
		;;
		--Installcutadapt)
			Installcutadapt=1
		;;
		--InstallFastQC)
			InstallFastQC=1
		;;
		--Installpigz)
			Installpigz=1
		;;
		--InstallQUAST)
			InstallQUAST=1
		;;
		--InstallSPADES)
			InstallSPADES=1
		;;
		--InstallPyFasta)
			InstallPyFasta=1
		;;
		--InstallSoftware|--Software)
			UpgradeBasics=1
			InstallLib=1
			InstallORFFinder=1
			InstalCDHit=1
			Installcutadapt=1
			InstallFastQC=1
			Installpigz=1
			InstallQUAST=1
			InstallSPADES=1
			InstallPyFasta=1
		;;
		--FullInstall|--Full)
			UpgradeBasics=1
			InstallLib=1
			InstallORFFinder=1
			InstalCDHit=1
			Installcutadapt=1
			InstallFastQC=1
			Installpigz=1
			InstallQUAST=1
			InstallSPADES=1
			InstallPyFasta=1
			DownloadGreengenes=1
			MakeGreengenesTaxonomyIndex=1
			MakeTaxonomyIndex_ggpattern=1
		;;
		--DownloadGreengenes)
			DownloadGreengenes=1
			MakeGreengenesTaxonomyIndex=1
		;;
		--MakeGreengenesTaxonomyIndex)
			MakeGreengenesTaxonomyIndex=1
		;;
		--MakeTaxonomyIndex_ggpattern)
			MakeTaxonomyIndex_ggpattern=1
		;;
		*)
			echo "Couldn't recognize command '${1}'. Ignoring it."
		;;
	esac
	shift
done

cd $address
if [[ -d "$LIB" ]]
then
	touch $LIB
else
	mkdir $LIB
fi
cd $LIB
LIB="$(pwd)"

if [[ ! $(${usearch} --version) ]]
then
	echo "USEARCH may not be properly installed."
else
	echo "Your current version of USEARCH is $(${usearch} --version)"
fi

echo "If you're running the install program for BEAF ($0), at this point you should already have USEARCH installed in your computer (and preferably in your PATH).
It is also recommended that you update your avaiable packages and versions by using the command
	'sudo apt-get update'
Upgrading your current packages may also be recommended (although that is only good practice)
	'sudo apt-get upgrade'"



if [[ $UpgradeBasics -eq 1 ]]
then
	sudo apt-get install --upgrade python python3 ## Download Python
	sudo apt-get install --upgrade python3-pip python-pip cmake ncbi-blast+ python-biopython git parallel ## Download basic tools

	sudo pip install --upgrade pip
	sudo pip install --upgrade setuptools biopython scikit-learn scipy numpy ## Download BioPython
	sudo pip3 install --upgrade pip
	sudo apt-get install -f
fi

if [[ $InstallLib -eq 1 ]]
then
	wget https://github.com/celiosantosjr/BEAF/blob/master/Lib/SeqLength.py
	wget https://github.com/celiosantosjr/BEAF/blob/master/Lib/PCA_maker.py
	wget https://github.com/celiosantosjr/BEAF/blob/master/Lib/ext.py
	if [[ ! -s bb.orffinder.pl ]]
	then
		wget https://github.com/vikas0633/perl/blob/master/orffinder.pl
		mv orffinder.pl bb.orffinder.pl
	fi
fi

if [[ $InstallORFFinder -eq 1 ]]
then
	if [[ ! -s bb.orffinder.pl ]]
	then
		wget https://github.com/vikas0633/perl/blob/master/orffinder.pl
		mv orffinder.pl bb.orffinder.pl
	fi
fi

if [[ $InstalCDHit -eq 1 ]]
then
	sudo apt-get install --upgrade cd-hit
	sudo apt-get install -f
fi

if [[ $Installcutadapt -eq 1 ]]
then
	if [[ $UpgradeBasics -eq 0 ]]
	then
		sudo apt-get install --upgrade python python3 ## Download Python
		sudo apt-get install --upgrade python3-pip python-pip 
		sudo pip install --upgrade pip
		sudo pip3 install --upgrade pip
	fi
	sudo pip install --upgrade cutadapt ## Download cutadapt
	sudo pip3 install --upgrade cutadapt ## Upgrade cutadapt for multiple cores
	sudo apt-get install -f
fi

if [[ $InstallFastQC -eq 1 ]]
then
	cd $LIB
	git clone https://github.com/s-andrews/FastQC
	cd FastQC
	fastqc="$(pwd)/fastqc"
	cd $LIB
fi

if [[ $Installpigz -eq 1 ]]
then
	cd $LIB
	wget https://zlib.net/pigz/pigz-2.4.tar.gz ## Download pigz
	tar -xzf pigz-2.4.tar.gz
	cd pigz-2.4
	sudo make
	pigz="$(pwd)/pigz"
	cd $LIB
fi

if [[ $InstallQUAST -eq 1 ]]
then
	cd $LIB
	git clone https://github.com/ablab/quast ## Download QUAST
	cd quast
	chmod +x install_full.sh 
	sudo sh install_full.sh install_full
	quast="$(pwd)/metaquast.py"     
	cd $LIB
fi

if [[ $InstallSPADES -eq 1 ]]
then
	cd $LIB
	git clone https://github.com/ablab/spades ## Download SPADES
	cd spades
	chmod +x make-targz.sh
	sudo sh ./make-targz.sh
	gunzip SPAdes-$(cat assembler/VERSION).gz
	cd SPAdes-$(cat assembler/VERSION)
	chmod +x spades_compile.sh
	sudo sh ./spades_compile.sh
	spades="$(pwd)/spades.py"
	cd $LIB
fi

if [[ $InstallPyFasta -eq 1 ]]
then
	cd $LIB
	git clone https://github.com/brentp/pyfasta ## Download PyFasta
	cd pyfasta
	sudo python setup.py install
	cd $LIB
fi

cd $address ## Leaving Lib


if [[ $DownloadGreengenes -eq 1 ]]
then
	cd $address
	if [[ -d "$ReferencesFolder" ]]
	then
		touch $ReferencesFolder
	else
		mkdir $ReferencesFolder
	fi
	cd $ReferencesFolder
	ReferencesFolder="$(pwd)"
	if [[ -d "Taxonomy" ]]
	then
		touch Taxonomy
	else
		mkdir Taxonomy
	fi
	cd Taxonomy
	wget http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/current_GREENGENES_gg16S_unaligned.fasta.gz ## Download Greengenes
	tar -xzf current_GREENGENES_gg16S_unaligned.fasta.gz
	rm -rf current_GREENGENES_gg16S_unaligned.fasta.gz
	cd $address
	cd $ReferencesFolder
	${cdhit} -i Taxonomy/current_GREENGENES_gg16S_unaligned.fasta -o Taxonomy/nr100_Greengenes.fasta -c 1.00 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5
	${usearch} -makeudb_usearch Taxonomy/nr100_Greengenes.fasta -output nr100_Greengenes.udb
	${cdhit} -i Taxonomy/current_GREENGENES_gg16S_unaligned.fasta -o Taxonomy/nr97_Greengenes.fasta -c 0.97 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5
	${usearch} -makeudb_usearch Taxonomy/nr97_Greengenes.fasta -output nr97_Greengenes.udb
	cd ..
fi

if [[ $MakeGreengenesTaxonomyIndex -eq 1 ]]
then
	cd $address
	if [[ -d "$ReferencesFolder" ]]
	then
		touch $ReferencesFolder
	else
		mkdir $ReferencesFolder
	fi
	cd $ReferencesFolder
	ReferencesFolder="$(pwd)"
	if [[ -d "Taxonomy" ]]
	then
		touch Taxonomy
	else
		mkdir Taxonomy
	fi
	if [[ -d "Taxonomy_Index" ]]
	then
		touch Taxonomy_Index
	else
		mkdir Taxonomy_Index
	fi
	grep ">" $ReferencesFolder/Taxonomy/nr100_Greengenes.fasta > Taxonomy_Index/nr100_Greengenes.headers.txt
	grep ">" $ReferencesFolder/Taxonomy/nr97_Greengenes.fasta > Taxonomy_Index/nr97_Greengenes.headers.txt
	cd Taxonomy_Index
	cat nr100_Greengenes.headers.txt | sed 's# .*; otu_#\totu_#' | sed 's#>##' | sort > nr100_Greengenes.ID_2_OTU.tsv
	cat nr100_Greengenes.headers.txt | sed 's#>##' | sed 's# .*k__#\t#' | sed 's#; p__#\t#' | sed 's#; c__#\t#' | sed 's#; o__#\t#' | sed 's#; f__#\t#' | sed 's#; g__#\t#' | sed 's#; s__#\t#' | sed 's#; otu_#\totu_#' | sed 's#; Unclassified#\tUnclassified#' | awk -F "otu_" '{print "otu_" $2 "\t" $1}' | awk -F "\t" '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > nr100_Greengenes.ID_2_Taxon.tsv
	cat nr97_Greengenes.headers.txt | sed 's# .*; otu_#\totu_#' | sed 's#>##' | sort > nr97_Greengenes.ID_2_OTU.tsv
	cat nr97_Greengenes.headers.txt | sed 's#>##' | sed 's# .*k__#\t#' | sed 's#; p__#\t#' | sed 's#; c__#\t#' | sed 's#; o__#\t#' | sed 's#; f__#\t#' | sed 's#; g__#\t#' | sed 's#; s__#\t#' | sed 's#; otu_#\totu_#' | sed 's#; Unclassified#\tUnclassified#' | awk -F "otu_" '{print "otu_" $2 "\t" $1}' | awk -F "\t" '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > nr97_Greengenes.ID_2_Taxon.tsv
	cd $address
fi

if [[ $MakeTaxonomyIndex -eq 1 ]]
then
	cd $address
	if [[ -d "$ReferencesFolder/Taxonomy" ]]
	then
		if [[ -d "Taxonomy_Index" ]]
		then
			touch Taxonomy_Index
		else
			mkdir Taxonomy_Index
		fi
		ls $ReferencesFolder/Taxonomy | grep -Ei "(.fasta|.fa|.faa|.fas|.fna|.fsa|.ffn|.frn|.mpfa)$" >> $ReferencesFolder/Taxonomy/FastaLIST_MakeIndex.txt
		for FastaFile in $(cat $ReferencesFolder/Taxonomy/FastaLIST_MakeIndex.txt); do
			if [[ ! -s $ReferencesFolder/${FastaFile%.f*}.udb ]]
			then
				if [[ $(cat $ReferencesFolder/Taxonomy/$FastaFile | grep ">" | sed 's# .*##' | sort | uniq -d) -gt 0 ]]
				then
					cat $ReferencesFolder/Taxonomy/$FastaFile | awk '/^>/{print ">SeqID" NR " BEAFwillDelete<" $0; next} {print}' | sed 's# BEAFwillDelete<># #' > $ReferencesFolder/Taxonomy/${FastaFile}.UNIQ_SEQ_IDS
					echo "$FastaFile had multiple sequences/headers with the same Sequence_ID (from the header format below)
	>[Sequence_ID] [Name_or_description_of_the_sequence] k__[Kingdom]; p__[Phylum]; c__[Class]; o__[Order]; f__[Family]; g__[Genus]; s__[Species]; otu_[OTU_ID]
This file was rewritten to include a Sequence_ID, by generating a sequence ID to each header and then renaming each header to >[Sequence_ID] [old header]."
					rm -rf $ReferencesFolder/Taxonomy/${FastaFile}
					mv $ReferencesFolder/Taxonomy/${FastaFile}.UNIQ_SEQ_IDS $ReferencesFolder/Taxonomy/${FastaFile}
				fi
				${usearch} -makeudb_usearch $ReferencesFolder/Taxonomy/$FastaFile -output $ReferencesFolder/Taxonomy/${FastaFile%.f*}.udb
			fi
			grep ">" $ReferencesFolder/Taxonomy/$FastaFile > $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.headers.txt
			if [[ $(cat $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.headers.txt | wc -l) -eq $(cat $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.headers.txt | grep "otu_" | wc -l) ]]
			then
				cat $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.headers.txt | sed 's# .*; otu_#\totu_#' | sed 's#>##' | sort > $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.ID_2_OTU.tsv
				cat $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.headers.txt | sed 's#>##' | sed 's# .*k__#\t#' | sed 's#; p__#\t#' | sed 's#; c__#\t#' | sed 's#; o__#\t#' | sed 's#; f__#\t#' | sed 's#; g__#\t#' | sed 's#; s__#\t#' | sed 's#; otu_#\totu_#' | sed 's#; Unclassified#\tUnclassified#' | awk -F "otu_" '{print "otu_" $2 "\t" $1}' | awk -F "\t" '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > $ReferencesFolder/Taxonomy_Index/${FastaFile%.f*}.ID_2_Taxon.tsv
			else
				echo "$FastaFile is not in Greengenes pattern format of 
	>[Sequence_ID] [Name_or_description_of_the_sequence] k__[Kingdom]; p__[Phylum]; c__[Class]; o__[Order]; f__[Family]; g__Genus; s__[Species]; otu_[OTU_ID]
This file wasn't properly indexed for Taxonomy finding. You may try to index it yourself by creating files ID_2_OTU.tsv and ID_2_Taxon.tsv manually, or simply let BEAF skip the step where it finds specific taxons to match to each sequence...
More information can be found using install.sh --help_taxonomy_index"
			fi
		done
	fi
	cd $address
fi

cd $address
rm -rf .test.file
cd Reference_seqs
ls *tar.gz > list
for file in `cat list`; do
	tar -zxvf $file
done
rm -rf list
ls *gz > list
for file in `cat list`; do
	gunzip $file
done
cd ..
chmod u+x BEAF.sh #nome do BEAF aqui
echo "Testing complete installation..."
echo "G	R	Test_sample/Alistipes_putredinis_DSM_17216.fna.fastq.gz	NA	Alistipes_putredinis_DSM_17216.fna	NA	Test_genome1
N	I	Test_sample/Alistipes_putredinis_DSM_17216.fna.fastq.gz	NA	transposon.fasta	transposon	Test_nt1
P	I	Test_sample/Alistipes_putredinis_DSM_17216.fna.fastq.gz	NA	bdg.fa	BDG	Test_prot1" > config.file
sh ./BEAF.sh > Run.log
rm -rf Test_sample
cd Reference_seqs
rm -rf *.fa *.fasta *.fna BDG DNA_pol
cd ..
echo "See in OUTPUT folder if all files read good."
echo "To test again, please download Test_sample and Reference_seqs folders again from source."
echo "########### Finished ###########"
