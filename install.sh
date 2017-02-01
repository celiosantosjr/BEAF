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
	unzip quast-master.zip 
	cd quast-master
	sh ./install_full.sh
	cd ..
	unzip spades*
	ls spades* > .list_spades
	for folder in `cat .list_spades`;
	do
		cd $folder
		chmod +x make-targz.sh
		sh ./make-targz.sh
	done
	rm -rf .list_spades
	mv .list_test list
	./BEAF10.11.65.sh
	rm -rf list
else
	echo "You actually do not have installed usearch, please install usearch in /usr/bin or /bin folders"
fi
rm -rf .test.file
cd Reference_seqs
ls *tar.gz > list
for file in `cat list`
do
tar -zxvf $file
done
rm -rf list
ls *gz > list
for file in `cat list`
do
gunzip $file
done
cd ..
chmod u+x BEAF1011.65.sh
echo "Testing complete installation..."
sh ./BEAF1011.65.sh > Run.log
rm -rf Test_sample
cd Reference_seqs
rm -rf *.fa *.fasta *.fna BDG DNA_pol
cd ..
echo "See in OUTPUT folder if all files read good."
echo "To test again, please download Test_sample and Reference_seqs folders again from source."
echo "########### Finished ###########"
