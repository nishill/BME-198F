# BME 198F - independent study course at UCSC

*To test and use GA4GH APIs*

```
wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.0.7.tar.gz
tar -xvf virtualenv-12.0.7.tar.gz
python virtualenv-12.0.7/virtualenv.py venv
source venv/bin/activate
pip install -r requirements.txt
```

## Files:
>
	populations.ipynb - ipython notebook for analysis and testing of the GA4GH APIS through
>							the 1000 genomes and SGDP datasets

	variantMatrix.py - variantMatrix module used as a replication tool for testing the use cases
							betweent two datasets
	
## Directories:
>
	1kgenomes_csids - premade callset id value, population key json dictionaries 
>						for speeding up the data analysis process in the 1000 Genomes dataset

	simons_csids - premade callset id value, population key json dictionaries
>						for speeding up the data analysis process in the SGDP dataset

	metadata_processing - metadata processing etl script for processing the metadata in the SGDP
							dataset

## Statistical Tests:
>
	Random Sampling genes - I randomly looked at a region on each chromosome in the genome and then randomly sampled 10 genes on each chromosom.
							Then, I displayed the variant allele frequencies for each variant in the randomly chosen gene. This process is a good way to
							summarize information about a genome as well as capability to compare population specific information between individuals genomes.
	![Settings Window](https://raw.github.com/nishill/BME-198F/master/raw-allele-frequencies.png)

>
	Bootstrapping - here I random sampled variant information in a subpopulation in order to gain insight into what the populations behavior is as a whole. I 
					graped the allele frequencies through a histogram and then fit the Gaussian distribution to the plot. 


## Cross Dataset Comparisons:

In this section I will explain the difficulty in comparing two variant sets when utilizing 
the GA4GH application programming interfaces. I hope that this will lead to improvements in the data model
so that the GA4GH APIs will become a helpful tool in comparing genomic data across different datasets. 

	I hosted the Simons Genome Diversity Project (SGDP) data through the GA4GH APIs. The code I have written
in the populations ipython notebook uses the functionality of the GA4GH APIs to test how researchers might use
two large-scale genomic datasets to gain meaningful insights through the use of the GA4GH APIs. With the GA4GH I 
had previously written code which creates a variant call adjacency matrix with the 1000 Genomes project APIs. To test
the SGDP analysis framework I attempted to replicate the variant matrix with the SGDP data. The result revealed 
differences in how the variant call data is structured as well as how the GA4GH APIs could be improved to 
accomodate these structural differences.

Specifically, all of the variant call information in the 1000 Genomes Project data is located in a single variant set.
In the SGDP dataset, there is a variant set for each individual. This creates an incompatibility in how these two datasets
could be analyzed together. Call set ids in one variant set cannot be compared to call set ids in another variant set. For example, 
in the 1000 genomes dataset, variant call data can be compared across individuals in the variant set because they all belong to a single dataset. 
In the SGDP dataset, individuals variant call data cannot be compared between other variant sets. Other challenges arose also
because some of the variants in the SGDP dataset did not have variant names associated with them so there was no way to compare one
variant between in one variant set to the same variant in another variant set.  

I tried to resolve this issue by altering my metadata processing script to include all of the individuals SGDP callset information 
in each variant set though could not proceed because of itegrity constraints in the sql database. 

I would like to accomodate in resolving this barrier in the future so it is easier for users of the API to compare different genomic datasets.  
