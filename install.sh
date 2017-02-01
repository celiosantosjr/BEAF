echo "At this point you should have installed Usearch in /usr/bin folder or in /bin folder to ensure complete and successful instalation"
which usearch > .test.file
if [ -s .test.file ]
then
	sudo apt-get install hmmer ncbi-blast+ python-biopython git
	sudo apt-get upgrade python-biopython
	pip install -U pip
	pip install --user --upgrade cutadapt
	git clone https://github.com/ablab/quast                        
	git clone https://github.com/ablab/spades
	git clone https://github.com/dsenalik/bb
	cd bb
	mv bb.orffinder ../bb.orffinder.pl
	cd ..
	chmod +x bb.orffinder.pl
	rm -rf bb
	cd quast
	sh ./install_full.sh
	cd ..
	cd spades
	chmod +x make-targz.sh
	sh ./make-targz.sh
	cd Reference_seqs
	ls *tar.gz > list
	for file in `cat list`
	do
		tar -zxvf $file
	done
	rm -rf list
	gunzip *.gz
	cd ..
	chmod u+x BEAF1011.65.sh
	chmod u+x soft_BEAF1011.65.sh
	echo "Testing complete installation..."
	sh ./BEAF1011.65.sh > Run.log
	sh ./soft_BEAF1011.65.sh > SoftRun.log
	rm -rf Test_sample
	cd Reference_seqs
	rm -rf *.fa *.fasta *.fna BDG DNA_pol list
	cd ..
	echo "See in OUTPUT folder if all files read good."
	echo "To test again, please download Test_sample and Reference_seqs folders again from source."
	echo "########### Finished ###########"
else
	echo "You actually do not have installed usearch, please install usearch in /usr/bin or /bin folders"
fi
rm -rf .test.file
