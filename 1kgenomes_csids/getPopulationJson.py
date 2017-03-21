from __future__ import print_function
import ga4gh.client.client as client
client = client.HttpClient("http://1kgenomes.ga4gh.org")
import json
from pprint import pprint 
import random 

"""
******************** get datasets and variant_sets ******************* 
"""
datasets = list(client.search_datasets())

# get 1000 gneomes dataset 
dataset = client.get_dataset(datasets[0].id)

release = None
functional = None
for variant_set in client.search_variant_sets(dataset_id=dataset.id):
	if variant_set.name == "phase3-release":
		release = variant_set
	else:
		functional = variant_set

callsi = list(client.search_call_sets(release.id))

variant_sets = list(client.search_variant_sets(dataset.id))
variant_set_id = variant_sets[0].id

call_set_ids = []

for csi in callsi: call_set_ids.append(csi.id)

"""
*********************************************************************** 
"""
# for getting data on populations 
def map_call_set_ids_to_population_group(pop_group):
	"""
	This allows for the individuals being analyzed to be identified
	based on 1000 genomes population group. 
	"""
	population_map = {}

	for call_set in client.search_call_sets(release.id):

		bio_sample = client.get_biosample(call_set.biosample_id)
		pop = bio_sample.info['Population'].values[0].string_value
		if pop == pop_group: population_map[call_set.id] = pop_group 

	return population_map

def main():

	pop_group = input("Enter the name of your 1000 genomes dataset: ")
	population_map = map_call_set_ids_to_population_group(pop_group)
	with open(	pop_group + ".json",'w') as output_file:
		json.dump(population_map,output_file)

if __name__=="__main__":
	main()
