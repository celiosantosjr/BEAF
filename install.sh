
### ADD PYFASTA


echo "At this point you should have installed Usearch in your computer to ensure complete and successful instalation"
which usearch > .test.file
if [ -s .test.file ]
then
	sudo apt-get install --upgrade python python3
	sudo apt-get install --upgrade ncbi-blast+ python-biopython git python3-pip python-pip cmake

	sudo pip install -U pip

	sudo pip install -U biopython scikit-learn scipy numpy cutadapt

	sudo pip3 install --upgrade cutadapt

	git clone https://github.com/ablab/quast                        
	git clone https://github.com/ablab/spades
	sudo sh quast/install_full.sh install_full
	cd spades
	SpadesVersion="$(cat assembler/VERSION)"
	chmod +x make-targz.sh
	sh ./make-targz.sh
	gunzip SPAdes-${SpadesVersion}.gz
	cd SPAdes-${SpadesVersion}
	sudo sh ./spades_compile.sh
	cd ../..
else
	echo "You actually do not have installed usearch, please install usearch."
fi
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
