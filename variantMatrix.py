from __future__ import print_function
import collections
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.axes_grid1 import ImageGrid
from numpy.random import RandomState
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
import unirest
import requests 
import json
from multiprocessing import Pool
from finch import *
import operator
import random


class VariantMatrix:

	def __init__(self, dataset, client, vs_ids, cs_ids, num_people, chrom, start, end ):

		self.client = client
		self.vs_ids = vs_ids
		self.cs_ids = cs_ids
		self.num_people = num_people
		self.chrom = chrom
		self.start = start
		self.end = end
		self.dataset = dataset

		#for i in range(0,len(vs_ids)):
		#`	print ( str(vs_ids[i]) + " " + str(cs_ids[i]) )


		self.variantMatrix = self.variant_matrix()


	def __str__(self):

		return str(self.variantMatrix)

	def __repr__(self):
		return ('\n'.join(','.join(str(m) for j in i) for i in self.variantMatrix))


	def generate_dictionary(self):
		"""
		this function generates a dictionary of call set id keys and values for the index
		that the call set id appeared in throughout a list of variants within a given range on the human genome.
   
		num_people - the number of people to analyze
		num_chrom - the chromosome number
		begin - the start position on the genome
		stop - the end position on the genome
		"""


		dictionary = {}
		#count = 0
		#for vsid in self.vs_ids[0:num_people]:
		#	variants = list(self.client.search_variants(vsid, call_set_ids= self.cs_ids,
		#		start=begin, end=stop, reference_name = num_chrom))
		#	i = 0
		#	for v in variaints:
		#		for call in v.calls:
		#			i += 1
		#			dicitionary[call.call_set_id] = i
		#return (dictionary)
		#i = 0
		#for csid in self.cs_ids:
		#	i += 1
		#	dictionary[csid] = i

		#dictionary = {}

		#for variant_set_id in self.vs_ids:
		variants = self.client.search_variants(self.vs_ids[0], call_set_ids= self.cs_ids[0:self.num_people],
		   		start=self.start, end=self.end, reference_name = self.chrom )

		i = 0
		for v in variants:
			for call_inner in v.calls:
				i += 1
				dictionary[call_inner.call_set_id] = i

		return (dictionary)


	def normalize_indexes(self,dictionary):
		"""
		initializes the indexes of the variant dictionay to be {0,...,n} where n is the number of call set ids

		dictionary - call set id keys, call set id index values
		"""
		keys = []
		keys = dictionary.keys()

		new_dictionary = {}
		i = 0
		j = 0
		k = 0
		for i in range(len(dictionary)):
			l = keys[j]
			dictionary[l] = k
			j += 1
			k += 1

		return dictionary


	def variant_matrix(self):
		"""
		num_people - the number of people to analyze
		num_chrom - the chromosome number
		begin - the start position on the genome
		stop - the end position on the genome
		"""

		# get dictionary of call set ids for the number of people
		csid_dict = {}
		csid_dict = self.generate_dictionary()

		self.normalize_indexes(csid_dict)

		h = len(csid_dict)
		w = len(csid_dict)
		dimension = h*w

		# initialize callMatrix with zeroes
		callMatrix = {}
		callMatrix = [[0 for x in range(w)] for y in range(h)]

		#i = 0
		#for vsid1 in self.vs_ids[0:num_people]:
		#	cs1 = self.cs_ids[i]
		#	variants1 = list(self.client.search_variants(vsid1, call_set_ids= [cs1],
		#		start=begin, end=stop, reference_name = num_chrom))
		#	i += 1
		#	j = 0
		#	for vsid2 in self.vs_ids[0:num_people]:
		#		cs2 = self.cs_ids[j]
		#		variants2 = list(self.client.search_variants(vsid2, call_set_ids= [cs2],
		#			start=begin, end=stop, reference_name = num_chrom))
		#		j += 1
		#		for v1 in variants1:
		#			for v2 in variants2:
		#				if v1.names[0] == v2.names[0]:
		#					for call_out in v1.calls:
		#						for call_in in v2.calls:
		#							outer_index = csid_dict[call_out.call_set_id]
		#							inner_index = csid_dict[call_in.call_set_id]
		#							A1 = call_out.genotype[0]
		#							A2 = call_out.genotype[1]
		#							B1 = call_in.genotype[0]
		#							B2 = call_in.genotype[1]
									# True for the shared genotype
		#							if (((A1 + A2 ) > 0) and ((B1 + B2) > 0)):
		#								callMatrix[outer_index][inner_index] += 1

		#							elif (((A1 + A2 ) == 0) and ((B1 + B2) == 0)):
		#								callMatrix[outer_index][inner_index] += 1

		#for variant_set in self.vs_ids[0:self.num_people]:
			
		#	vs = list(self.client.search_variants(variant_set, call_set_ids= self.cs_ids,
		#		start=self.start, end=self.end, reference_name = self.chrom))

			# iterate through each persons call set id and increment cells in the matrix for shared genotypes
		#	for v in vs:
		#		for call_outer in v.calls:
		#			for call_inner in v.calls:
		#				outer_index = csid_dict[call_outer.call_set_id]
		#				inner_index = csid_dict[call_inner.call_set_id]
				 		# can be 0 or 1
		#				A1 = call_outer.genotype[0]
		#				A2 = call_outer.genotype[1]
		#				B1 = call_inner.genotype[0]
		#				B2 = call_inner.genotype[1]
						# True for the shared genotype
		#				if (((A1 + A2 ) > 0) and ((B1 + B2) > 0)):
		#					callMatrix[outer_index][inner_index] += 1

		#				elif (((A1 + A2 ) == 0) and ((B1 + B2) == 0)):
		#					callMatrix[outer_index][inner_index] += 1

		#for variant_set_id in self.vs_ids:
		vs = self.client.search_variants(self.vs_ids[0], call_set_ids= self.cs_ids[0:self.num_people],
		   	start=self.start, end=self.end, reference_name = self.chrom)

		# iterate through each persons call set id and increment cells in the matrix for shared genotypes
		for v in vs:
			for call_outer in v.calls:
				for call_inner in v.calls:
					outer_index = csid_dict[call_outer.call_set_id]
					inner_index = csid_dict[call_inner.call_set_id]
					# can be 0 or 1
					A1 = call_outer.genotype[0]
					A2 = call_outer.genotype[1]
					B1 = call_inner.genotype[0]
					B2 = call_inner.genotype[1]
					# True for the shared genotype
					if (((A1 + A2 ) > 0) and ((B1 + B2) > 0)):
						callMatrix[outer_index][inner_index] += 1

					elif (((A1 + A2 ) == 0) and ((B1 + B2) == 0)):
						callMatrix[outer_index][inner_index] += 1

		return (callMatrix)


	def visualize_matrix_color(self):

		"""
		visualizes the matrix in rgb color

		num_people - the number of people to analyze
		num_chrom - the chromosome number
		begin - the start position on the genome
		stop - the end position on the genome
		"""

		fig = plt.figure(1,(10.,10.))
		grid = ImageGrid(fig, 111,
					nrows_ncols=(1,1),
					axes_pad=0.1)

		ax = grid[0]
		ax.set_title('Color matrix comparing the occurence of shared variants obtained from call set ids\n', fontsize=14, fontweight='bold')
		ax.set_xlabel("Individual's call set ids", fontsize=12)
		ax.set_ylabel("Individual's call set ids", fontsize=12)
		ax.imshow(self.variantMatrix, origin = "lower", interpolation="nearest")

		plt.show()

	
	def k_means_clustering_on_variant_matrix(self,clusters):
		"""
		returns a list of labels given by the k-means clustering algorithm which indicate which
		individuals should be clustered together

		num_people - the number of people to analyze
		num_chrom - the chromosome number
		begin - the start position on the genome
		stop - the end position on the genome
		clusters - the number of clusters to form with this run on kmeans
		"""

		# number of individuals ()
		N = self.num_people
		# randomization
		rand = RandomState(1123581321)
		# matrix generation
		# initialize KMeans clustering object
		kmeans_obj = KMeans(n_clusters=clusters, n_init=10,
					init='k-means++', precompute_distances=True,
					tol=1e-4, random_state= rand)
		# Compute clustering and transform matrix to cluster-distance space.
		labels = kmeans_obj.fit_predict(self.variantMatrix)

		label_list = []
		for i in labels:
			label_list.append(i)

		return label_list

	def assign_labels_to_indexes(self,clusters):

		"""
		assigns the labels from the kmeans clustering to the individuals in a variant matrix

		num_people - the number of people to analyze
		num_chrom - the chromosome number
		begin - the start position on the genome
		stop - the end position on the genome
		clusters - the number of clusters to form with this run on kmeans

		returns a dictionary of call set id keys, (label, matrix index( individual )) pairs values
		"""

		label_list = self.k_means_clustering_on_variant_matrix(self.num_people,
													  self.num_chrom, self.begin, self.end, clusters)
		dictionary = self.generate_dictionary(self.num_people, self.num_chrom, self.begin, self.end)
		dictionary = self.normalize_indexes(dictionary)

		i = 0
		j = 0
		k = 0
		label_index_dict = {}
		label_index_pair_list = []
		for i in label_list:
			pair = (i,j)
			label_index_pair_list.append(pair)
			j += 1


		label_index_pair_list.sort()

		i=0
		j=0
		k=0
		for i,j in dictionary.items():
			l = label_index_pair_list[k]
			dictionary[i] =l
			k += 1

		return (dictionary)

	def map_call_set_ids_to_population_group(self):
		"""
		This allows for the individuals being analyzed to be identified based on 1000 genomes population group.
		"""
		population_map = {}

		for call_set in self.client.search_call_sets(dataset.id):

			bio_sample = self.client.get_bio_sample(call_set.bio_sample_id)
			population_map[call_set.id] = bio_sample.info['Population'].values[0].string_value


		return population_map


	def compare(self):

		"""
		compares populations with the labels that they have been assigned to

		returns population occurences
		"""

		matrix = {}
		label_dictionary = {}
		population_dictionary= {}

		label_dictionary = self.assign_labels_to_indexes(50, "1", 100000, 200000, 6)
		population_dictionary = self.map_call_set_ids_to_population_group()

		population_occurence = []

		for ld_key, ld_val in label_dictionary.iteritems():
			for pd_key, pd_val in population_dictionary.iteritems():
				if ( ld_key == pd_key):
					print ("Population: {} has cluster {}".format(pd_val, ld_val[0]))
					population_occurence.append((pd_val, ld_val[0]))

		return population_occurence

