address=$(cd "$(dirname "")" && pwd)/$(basename "")
LIB=Lib
ReferencesFolder=Reference_seqs

Installcutadapt=0
Installpigz=0
InstallQUAST=0
InstallSPADES=0
InstallPyFasta=0

DownloadGreengenes=0
MakeTaxonomyIndex=0

while [[ $# -gt 0 ]]
do
	case $1 in
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
		--Installcutadapt|--cutadapt)
			Installcutadapt=1
		;;
		--Installpigz|--pigz)
			Installpigz=1
		;;
		--InstallQUAST|--QUAST)
			InstallQUAST=1
		;;
		--InstallSPADES|--SPADES)
			InstallSPADES=1
		;;
		--InstallPyFasta|--PyFasta)
			InstallPyFasta=1
		;;
		--DownloadGreengenes|--Greengenes)
			DownloadGreengenes=1
			MakeTaxonomyIndex=1
		;;
		--MakeGreengenesTaxonomyIndex)
			MakeGreengenesTaxonomyIndex=1
		;;
		--InstallSoftware|--Software)
			Installcutadapt=1
			Installpigz=1
			InstallQUAST=1
			InstallSPADES=1
			InstallPyFasta=1
		;;
		--FullInstall|--Full)
			Installcutadapt=1
			Installpigz=1
			InstallQUAST=1
			InstallSPADES=1
			InstallPyFasta=1
			DownloadGreengenes=1
			MakeGreengenesTaxonomyIndex=1
		;;
		*)
			echo "Couldn't recognize command '${1}'. Ignoring it."
		;;
	esac
	shift
done



mkdir $LIB
cd $LIB
LIB="$(pwd)"

echo "At this point you should have installed Usearch in your computer to ensure complete and successful instalation"

sudo apt-get install --upgrade python python3 ## Download Python
sudo apt-get install --upgrade python3-pip python-pip cmake ncbi-blast+ python-biopython git cd-hit ## Download basic tools

sudo pip install --upgrade pip
sudo pip install --upgrade setuptools biopython scikit-learn scipy numpy ## Download BioPython and cutadapt


if [[ $Installcutadapt -eq 1 ]]
then
	sudo pip install --upgrade cutadapt
	sudo pip3 install --upgrade pip
	sudo pip3 install --upgrade cutadapt ## Upgrade cutadapt for multiple cores
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
	git clone https://github.com/brentp/pyfasta ## Download PyFasta
	cd pyfasta
	sudo python setup.py install
	cd $LIB
fi

cd $address ## Leaving Lib


if [[ $DownloadGreengenes -eq 1 ]]
then
	mkdir $ReferencesFolder
	cd $ReferencesFolder
	ReferencesFolder="$(pwd)"
	wget http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/current_GREENGENES_gg16S_unaligned.fasta.gz ## Download Greengenes
	tar -xzf current_GREENGENES_gg16S_unaligned.fasta.gz
	cd-hit -i current_GREENGENES_gg16S_unaligned.fasta -o nr100_Greengenes.udb -c 1.00 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5
	cd-hit -i current_GREENGENES_gg16S_unaligned.fasta -o nr97_Greengenes.udb -c 0.97 -aS 1.0 -g 1 -d 0 -M 0 -T 0 -n 5
	cd ..
fi

if [[ $MakeGreengenesTaxonomyIndex -eq 1 ]]
then
	cd $ReferencesFolder
	ReferencesFolder="$(pwd)"
	mkdir Taxonomy
	grep ">" $ReferencesFolder/nr100_Greengenes.udb > Taxonomy/nr100_Greengenes.headers.txt
	cat nr100_Greengenes.headers.txt | sed 's# .*; otu_#\totu_#' | sed 's#>##' | sort > nr100_Greengenes.ID_2_OTU.tsv
	cat nr100_Greengenes.headers.txt | sed 's#>##' | sed 's# .*k__#\t#' | sed 's#; p__#\t#' | sed 's#; c__#\t#' | sed 's#; o__#\t#' | sed 's#; f__#\t#' | sed 's#; g__#\t#' | sed 's#; s__#\t#' | sed 's#; otu_#\totu_#' | sed 's#; Unclassified#\tUnclassified#' | awk -F "otu_" '{print "otu_" $2 "\t" $1}' | awk -F "\t" '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > nr100_Greengenes.ID_2_Taxon.tsv
	cat nr97_Greengenes.headers.txt | sed 's# .*; otu_#\totu_#' | sed 's#>##' | sort > nr97_Greengenes.ID_2_OTU.tsv
	cat nr97_Greengenes.headers.txt | sed 's#>##' | sed 's# .*k__#\t#' | sed 's#; p__#\t#' | sed 's#; c__#\t#' | sed 's#; o__#\t#' | sed 's#; f__#\t#' | sed 's#; g__#\t#' | sed 's#; s__#\t#' | sed 's#; otu_#\totu_#' | sed 's#; Unclassified#\tUnclassified#' | awk -F "otu_" '{print "otu_" $2 "\t" $1}' | awk -F "\t" '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > nr97_Greengenes.ID_2_Taxon.tsv
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
