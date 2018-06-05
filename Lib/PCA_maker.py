#!/usr/bin/python

from Bio import SeqIO
import sys
from scipy.spatial import distance
import numpy as np
from sklearn.decomposition import PCA
import matplotlib
matplotlib.get_cachedir()
import matplotlib.pyplot as plt
import time

def calc_tetra(seq_record):

	kmers = {}
	for a in ['A', 'C', 'G', 'T']:
		for b in ['A', 'C', 'G', 'T']:
			for c in ['A', 'C', 'G', 'T']:
				for d in ['A', 'C', 'G', 'T']:
					kmers[a+b+c+d] = 0

	start = 0
	end = 4	
	for i in range(0,len(str(seq_record.seq))):
		if len(str(seq_record.seq[start:end])) == 4:
			try:
				kmers[str(seq_record.seq[start:end])] += 1
			except:
				pass	
		start += 1
		end += 1
	
	total = sum(kmers.values())
	for k in kmers.keys():
		kmers[k] = float(kmers[k]) / float(total)

	return kmers

def calc_parameters(input_file):

	print "Calculating 4mer frequencies..."
	kmers = {}
	sz = {}
	input_handle = open(input_file, "rU")
	for record in SeqIO.parse(input_handle, "fasta") :
		kmers[record.id] = calc_tetra(record)
		sz[record.id] = len(record.seq)

	return kmers, sz

def plot_pca(sz):

	print "Running PCA on 4mer distances..."
	array = []
	for t in sorted(kmers.keys()):
		temp = []
		for tet in sorted(kmers[t].keys()):
			temp.append(kmers[t][tet])
		array.append(temp)

	kmers_np = np.array(array)
	pca = PCA(n_components=2)
	fit = pca.fit(kmers_np).transform(kmers_np)

	for sr_c in sz:
		sz[sr_c] = float(sz[sr_c]) / float(max(sz.values())) * 250 + 35
	sz_list = []
	for sr_c in sorted(sz.keys()):
		sz_list.append(sz[sr_c])

	print "Plotting PCA graph..."
	fig, ax = plt.subplots()
	ax.scatter(fit[:,0], fit[:,1], s=sz_list, c='#1f78b4', alpha=0.9)
	positions = {}
	i = 0 
	avg_size = sum(sz.values()) / len(sz.values())
	for name in positions.keys():
		x = positions[name][0]
		y = positions[name][1]
	fig.savefig(output)
	fig.clf()
	ax.cla()
	plt.close(fig)

if __name__ == '__main__':

	start = time.time()
	input_file = sys.argv[1]
	output = sys.argv[2]
	kmers, sz = calc_parameters(input_file)
	try: plot_pca(sz)
	except:	pass
	end = time.time()
	print "Run time: "
	print(end - start)
