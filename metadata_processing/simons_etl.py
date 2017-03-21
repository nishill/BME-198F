"""
Simons Diversity Project Metadata Processing: simons_etl.py
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import shutil
import json
import pysam
import utils
import sys
import generate_gff3_db
import tempfile
import zipfile
import csv
import datetime
import string 
import sqlite3
import urllib 
import glob

import ga4gh.server.datarepo as datarepo  # NOQA
import ga4gh.server.datamodel.references as references  # NOQA
import ga4gh.server.datamodel.datasets as datasets  # NOQA
import ga4gh.server.datamodel.variants as variants  # NOQA
import ga4gh.server.datamodel.reads as reads  # NOQA
import ga4gh.server.datamodel.ontologies as ontologies  # NOQA
import ga4gh.server.datamodel.sequence_annotations as sequenceAnnotations  # NOQA
import ga4gh.server.datamodel.bio_metadata as biodata  # NOQA
from ga4gh.server.datamodel.variants import CallSet
# save_files_locally()
# Requires wget
def save_files_locally(data):
    print("Gonna download {} indexes!".format(len(data)))
    for row in data:
        download_url = make_address(row['name'], os.path.basename(row['indexUrl']))
        os.system("wget {}".format(download_url))


# parses population metadata in csv file and returns
# biosample dictionaries
def parse_file_biosamples(filename):
  bio_samples = []
  individuals = []
  print("Loading bio data tsv")
  with open(filename, 'rU') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      for row in reader:
        description = "{}{}".format(
          row['ID'],
          row['Description'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        biosample = {
             "name": row['Name'],
             "description": description,
             "disease": None,  # Ontology term
             "created": datetime.datetime.now().isoformat(),
             "updated": datetime.datetime.now().isoformat(),
             "info": info
        }
        if row['Gender'] == 'M':
           sex = {
               "id": "PATO:0020001",
               "term": "male genotypic sex",
               "sourceName": "PATO",
               "sourceVersion": "2015-11-18"
        }
        elif row['Gender'] == 'F':
          sex = {
            "id": "PATO:0020002",
            "term": "female genotypic sex",
            "sourceName": "PATO",
            "sourceVersion": "2015-11-18"
          }
        else:
          sex = None
        bio_samples.append(biosample)
  
  return bio_samples

# parses population metadata in csv file and returns
# individual dictionaries
def parse_file_individuals(filename):
  individuals = []
  print("Loading individual data tsv")
  with open(filename, 'rU') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      for row in reader:
        description = "{}{}".format(
          row['ID'],
          row['Description'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        if row['Sex'] == 'M':
           sex = {
               "id": "PATO:0020001",
               "term": "male genotypic sex",
               "sourceName": "PATO",
               "sourceVersion": "2015-11-18"
        }
        elif row['Sex'] == 'F':
          sex = {
            "id": "PATO:0020002",
            "term": "female genotypic sex",
            "sourceName": "PATO",
            "sourceVersion": "2015-11-18"
          }
        else:
          sex = None
        individual = {
               "name": row['Name'],
               "description": description,
               "created": datetime.datetime.now().isoformat(),
               "updated": datetime.datetime.now().isoformat(),
               "species": {
                   "term": "Homo sapiens",
                   "id": "NCBITaxon:9606",
                   "sourceName": "http://purl.obolibrary.org/obo",
                   "sourceVersion": "2016-02-02"
               },
               "sex": sex,
               "info": info
        }
        individuals.append(individual)
  
  return individuals

# main():
# populates database relations with data from each person
# in both the individual and biosample directories 
#@utils.Timed()
def main():

    # Set for using hg38 rather than hg19
    # reference_set_path = '/mnt/ga4gh/repo_data/hg38.fa.gz'
    reference_set_path = '/mnt/ga4gh/repo_data/hs37d5.fa.gz'	
   
    bio_tsv_location = 'SGDP_metadata.279public.21signedLetter.samples.Biosample.tsv'
    ind_tsv_location = 'SGDP_metadata.279public.21signedLetter.samples.individual.tsv' 
  
    bio_samples = parse_file_biosamples( bio_tsv_location )
    individuals = parse_file_individuals( ind_tsv_location )
    repoPath = os.path.join("repo2.db")
    repo = datarepo.SqlDataRepository(repoPath)
    if ( os.path.isfile("repo2.db") == True ):
        os.system( "rm repo2.db" )
    repo.open("w")
    repo.initialise()
    
    dataset = datasets.Dataset("Simons")
    dataset.setDescription("Variants from the Simons Foundation Genome Diversity Project")
    repo.insertDataset(dataset)
    
    print("Inserting biosamples")
    new_bio_samples = []
    for bio_sample in bio_samples:
      new_bio_sample = biodata.Biosample(dataset, unicode(bio_sample['name'], errors='replace'))
      new_bio_sample.populateFromJson(json.dumps(bio_sample))
      repo.insertBiosample(new_bio_sample)
      new_bio_samples.append(new_bio_sample)
    
    print("Inserting individuals")
    new_individuals= []
    for individual in individuals:
      new_individual = biodata.Individual(dataset, unicode(individual['name'], errors='replace' ))
      new_individual.populateFromJson(json.dumps(individual))
      repo.insertIndividual(new_individual)
      new_individuals.append(new_individual)
    
    print("Adding reference set (takes a while)")
    reference_set = references.HtslibReferenceSet("NCBI37")
    reference_set.populateFromFile(reference_set_path)
    reference_set.setDescription("NCBI37 assembly of the human genome")
    reference_set.setNcbiTaxonId(9606)
    reference_set.setSourceUri("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz")
    for reference in reference_set.getReferences():
        reference.setNcbiTaxonId(9606)
    repo.insertReferenceSet(reference_set)
    
    seq_ontology = ontologies.Ontology("/mnt/ga4gh/repo_data/so-xp")	
    ontology_file_path = '/mnt/ga4gh/repo_data/so-xp-simple.obo'
    seq_ontology.populateFromFile(ontology_file_path)
    seq_ontology._id = "so-xp"
    repo.insertOntology(seq_ontology)
    repo.addOntology(seq_ontology)
 
    vcf_directory = os.path.dirname('/mnt/ga4gh/data/vcf/') 
    pattern = os.path.join(vcf_directory, "*.vcf.gz")
    for vcfFile in glob.glob(pattern):
        name = vcfFile.replace("/mnt/ga4gh/data/vcf/","")
        name = name.replace(".annotated.nh2.variants.vcf.gz","") 
        print (name)
        variant_set = variants.HtslibVariantSet(dataset, name)
        variant_set.setReferenceSet(reference_set)
        variant_set.populateFromFile([vcfFile], [vcfFile + ".tbi"])
        variant_set.checkConsistency()
        for call_set in variant_set.getCallSets():
            for bio_sample in new_bio_samples:
                if bio_sample.getLocalId() == call_set.getLocalId():
                    call_set.setBioSampleId(bio_sample.getId())
        
        repo.insertVariantSet(variant_set)

        name = name+"-annotated-nh2"
        print (name)
        variant_set2 = variants.HtslibVariantSet(dataset, name)
        variant_set2.setReferenceSet(reference_set)
        variant_set2.populateFromFile([vcfFile], [vcfFile + ".tbi"])
        variant_set2.checkConsistency()
        repo.insertVariantSet(variant_set2)
        for annotation_set in variant_set2.getVariantAnnotationSets():
            print ( str(annotation_set) + "found" )
            annotation_set.setOntology(seq_ontology)
            repo.insertVariantAnnotationSet(annotation_set)
		 
    repo.commit()
    print ( "database filled!")			

if __name__ == "__main__":
    main()
