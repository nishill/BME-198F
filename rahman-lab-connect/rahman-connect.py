from __future__ import print_function
import ga4gh.client.client as client
rahman_client = client.HttpClient("http://52.160.96.216/ga4gh")

"""
******************** get datasets and variant_sets ******************* 
"""
datasets = list(rahman_client.search_datasets())

dataset = rahman_client.get_dataset(datasets[0].id)

release = None
functional = None
for variant_set in rahman_client.search_variant_sets(dataset_id=dataset.id):
    if variant_set.name == "phase3-release":
        release = variant_set
    else:
        functional = variant_set

"""
*********************************************************************** 
"""

def main():
	callsi = list(rahman_client.search_call_sets(functional.id))
	variant_sets = list(rahman_client.search_variant_sets(dataset.id))		
	variant_set_id = variant_sets[0].id
	print ( variant_set_id )	

	call_set_ids = []
	callsi = list(rahman_client.search_call_sets(functional.id))
	for csi in callsi:call_set_ids.append(csi.id)

	call_set_id = callsi[0]
	
	print (dataset)	
	print ( call_set_id)	

	exampleVariant = rahman_client.search_variants(variant_set_id=functional.id, start= 87512392 , end= 87512442, reference_name="1").next()
	if bool(exampleVariant.names) == False:
		print("Variant name: ?")
	else:
		print("Variant name: {}".format(exampleVariant.names[0]))
	print("Start: {}, End: {}".format(exampleVariant.start, exampleVariant.end))
	print("Reference bases: {}".format(exampleVariant.reference_bases))
	print("Alternate bases: {}".format(exampleVariant.alternate_bases))
	print("Number of calls: {}".format(len(exampleVariant.calls)))	


if __name__=="__main__":
	main()
