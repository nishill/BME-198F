from __future__ import print_function
import ga4gh.client.client as client
simons_client = client.HttpClient("http://10.50.100.241/")
import json
"""
******************** get datasets and variant_sets ******************* 
"""
datasets = list(simons_client.search_datasets())
dataset = simons_client.get_dataset(datasets[0].id)

simons_csid_dict = {}
for biosample in simons_client.search_biosamples(dataset_id=dataset.id):
    for variant_set in simons_client.search_variant_sets(dataset_id=dataset.id):
        if variant_set.name == biosample.info['Individual_id'].values[0].string_value:
            print(variant_set.name)
            print ( variant_set.id)
            callset = list(simons_client.search_call_sets(variant_set.id))[0]
            print(callset.id)
            print (biosample.info['Name'].values[0].string_value)
            simons_csid_dict[callset.id] = (biosample.info['Name'].values[0].string_value,variant_set.id)

with open(	"simons_csids.json",'w') as output_file:
	json.dump(simons_csid_dict,output_file)

