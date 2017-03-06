echo "#### At this point you should have installed Usearch in /usr/bin folder or in
 /bin folder to ensure complete and successful instalation ####"
which usearch > .test.file
if [ -s .test.file ]
then
	sudo apt-get install hmmer ncbi-blast+ python-biopython git
	sudo apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib
	sudo apt-get upgrade python-biopython
	pip install -U pip
	pip install --user --upgrade cutadapt
	git clone https://github.com/ablab/quast                        
	git clone https://github.com/dsenalik/bb
	cd bb
	mv bb.orffinder ../bb.orffinder.pl
	chmod +x bb.orffinder.pl
	cd ../quast
	rm -rf ../bb
	chmod +x install_full.sh
	sh ./install_full.sh
	cd ..
	wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz
	tar -xzf SPAdes-3.10.1-Linux.tar.gz
	mv SPAdes-3.10.1-Linux spades
	rm -rvf *.tar.gz
	cd Reference_seqs
	ls *tar.gz > list
	for file in `cat list`
	do
		tar -zxvf $file
	done
	rm -rf list
	gunzip *.gz
	cd ..
	chmod u+x BEAF.sh
	chmod u+x optmizer.sh
	echo "########### Finished ###########"
	echo "PLEASE USE THE FOLLOWING COMMANDS:

$ ./BEAF.sh

Option -b. Look in OUTPUT folder and observe if all results are consistent
After this, repeat the operation using -s. Repeat the verification."
else
	echo "##### You actually do not have installed usearch, 
please install usearch in /usr/bin or /bin folders ######"
fi
