#!/usr/bin/python

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

if sys.argv[1] == '--AddLength':
	for header, sequence in SimpleFastaParser(open(sys.argv[2], 'rU')):
		print '>' + header + ';length=' + str(len(sequence)) + '\n' + sequence    		# With this, simply add `;length=(sequence length)` to the header of each sequence (printing the fasta sequence)
else:
	if sys.argv[1] == '--AddLength_Cut-off-MinLength':
		for header, sequence in SimpleFastaParser(open(sys.argv[3], 'rU')):
			if ( int(len(sequence)) >= int(sys.argv[2]) ):
				print '>' + header + ';length=' + str(len(sequence)) + '\n' + sequence    		# same as above, but only prints sequences with length equal or higher to a cut-off (second argument given)
	else:
		if sys.argv[1] == '--AddLength_Cut-off-MaxLength':
			for header, sequence in SimpleFastaParser(open(sys.argv[3], 'rU')):
				if ( int(len(sequence)) <= int(sys.argv[2]) ):
					print '>' + header + ';length=' + str(len(sequence)) + '\n' + sequence    		# same as above, but only prints sequences with length equal or smaller to a cut-off (second argument given)
		else:
			if sys.argv[1] == '--Cut-off-MinLength':
				for header, sequence in SimpleFastaParser(open(sys.argv[3], 'rU')):
					if ( int(len(sequence)) >= int(sys.argv[2]) ):
						print '>' + header + '\n' + sequence    		# Only prints sequences with length equal or higher to a cut-off (second argument given), without adding anything
			else:
				if sys.argv[1] == '--Cut-off-MaxLength':
					for header, sequence in SimpleFastaParser(open(sys.argv[3], 'rU')):
						if ( int(len(sequence)) <= int(sys.argv[2]) ):
							print '>' + header + '\n' + sequence    		# Only prints sequences with length equal or smaller to a cut-off (second argument given), without adding anything
				else:
					if sys.argv[1] == '--OnlyPrintLength':
						for header, sequence in SimpleFastaParser(open(sys.argv[2], 'rU')):
							print len(sequence)								# With this, only print the length of each sequence, without their headers
					else:
						for header, sequence in SimpleFastaParser(open(sys.argv[1], 'rU')):
							print len(sequence)								# With this, only print the length of each sequence, without their headers



# from Bio import SeqIO
# for SeqInfo in SeqIO.parse(open(sys.argv[1], 'rU'), 'fasta'):
#	print (len(SeqInfo))            								# With this, only print the length of each sequence, without their headers
#	print SeqInfo.description + "; length=" + str(len(SeqInfo))     				# With this, prints the full name of the sequence, followed by `; length=(sequence length)` (prints only headers)
#	print ">" + SeqInfo.description + "; length=" + str(len(SeqInfo)) + "\n" + SeqInfo.seq    	# With this, simply add `; length=(sequence length)` to the header of each sequence (printing the fasta sequence)
#	print SeqInfo.name, len(SeqInfo)   								# With this, prints the name of the sequence (anything before a space or ;), followed by it's length (only headers)

