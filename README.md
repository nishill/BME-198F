# BME 198F - independent study course at UCSC

>
	This repository contains all of the code used for independent coursework with the GA4GH in
	the Genomics Institute at UCSC, winter quarter 2017. The first task I completed was to host
	Nanzeen Rahman's variant call data through the GA4GH APIs with a virtual machine. The second 
	task I completed was to host the Simons Foundation Genome Diversity (SGDP) on another 
	vritual machine, accessible through the GA4GH APIs. I obtained access to the SGDP dataset
	through the Open Science Grid (OSG). The final task I completed was to do statistical analysis
	of populations in both the 1000 Genomes and SGDP datasets. This allowed thorough testing of
	the APIs and demonstrations of their application to Bioinformatics.  



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
	populations.ipynb:
	ipython notebook for analysis and testing of the GA4GH APIS through the 1000 genomes and SGDP datasets

	variantMatrix.py: 
	variantMatrix module used as a replication tool for testing the use cases betweent two datasets
	
## Directories:
>
	1kgenomes_csids:
	premade callset id value, population key json dictionaries for speeding up the data analysis process in the 1000 Genomes dataset

	simons_csids: 
	premade callset id value, population key json dictionaries for speeding up the data analysis process in the SGDP dataset

	metadata_processing:
	metadata processing etl script for processing the metadata in the SGDP dataset

	rahman-connect:
	quick code demonstration of connecting to the GA4GH Nanzeen Rahman lab data

	pics:
	images of Bioinformatics models and error message

## Statistical Tests:
>
	Random Sampling genes:
	I randomly looked at a region on each chromosome in the genome and then randomly sampled 10 genes on each chromosome.
	I then displayed the variant allele frequencies for each variant in the randomly chosen gene. This process is a good 
	way to summarize information about a genome as well as capability to compare population specific information between 
	individuals genomes.

![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/raw-allele-frequencies.png)

>
	Bootstrapping:
	Bootstrapping allows you to measure the accuracy in your dataset. I used it in order to gain insight into whole 
	subpopulations. To bootstrap I random sampled variant information in a subpopulation for 100 iterations. Then I 
	graphed the resulting allele frequencies for each iteration and made a histogram of the data. I then constructed the 
	Gaussian Distribution with parameters mean(data) and standard_deviation(data). I did this for the Great Britain and 
	Japanese subpopulations in the 1000 Genomes dataset. These variants are in the gene SPTA1.  

![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/gbr_rv.png)
![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/jpt_rv.png)
![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/gbr_cv.png)
![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/jpt_cv.png)

## Cross Dataset Comparisons:

In this section I will explain some challenges I had when comparing two variant sets through the 
GA4GH application programming interfaces. I hope that this will lead to improvements in the data model
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
in the 1000 genomes dataset, variant call data can be compared across individuals in the variant set because they all belong 
to a single dataset. In the SGDP dataset, individuals variant call data cannot be compared between other variant sets. Other 
challenges arose when some of the variants in the SGDP dataset did not have variant names associated with them so there was no 
way to compare one variant in one variant set to the same variant in another variant set.  

I tried to resolve this issue by altering my metadata processing script to include all of the individuals SGDP callset information 
in each variant set though could not proceed because of itegrity constraints in the sql database. 

I would like to make call set information available through all vairant sets in the future so it is easier for 
users of the API to compare different genomic datasets. 

Output of Variant Call Matrix in 1kgenomes Chinese subpopulation
![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/CHSvcm.png) 

Error message from constructing Variant Call Matrix in SGDP dataset
![Settings Window](https://github.com/nishill/BME-198F/blob/master/pics/ga4gh_error.png) 
