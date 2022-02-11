# <i>Neotoma albigula</i> Parasite-Microbiome Study
This repository contains R code used to process parasitism and 16S microbiome data.
Repository is associated with publication: 
Doolin, ML, Weinstein, SB, Dearing, MD. <i>In review</i>. Pinworms are associate with taxonomic but not functional differences in the gut microbiome of white-throated woodrats (<i>Neotoma albigula</i>). Journal of Parasitology. 
Work associated with the University of Utah.

## Project Description
Parasites are known to change the host microbiome in controlled systems, and we set investigated whether natural parasite infections in a wild system would have similar impacts on gut microbial communities. For this work, we characterized gut microbial communities, parasite infection prevalence and egg shedding intensities, and general animal health data. We investigated animal health and microbial communities in light of parasite diagnoses. In this repository, we include R code and metadata necessary to reproduce our findings. Raw microbial 16S sequence data is deposited in the GenBank Sequence Read Archive under Bioproject BioProject PRJNA757000, BioSample SAMN20931663. 

## Files

AlbigulaPaper16SScript.R -- R script containing 16S sequence processing pipeline, from raw sequences to diversity analyses.

AlbigulaPaperPsiteFunctionScript.R -- Script containing code to perform parasite prevalence and model testing for association between parasitism and digestive functional metrics. 

AllSampsAvg_Metadat.csv -- Metadata associated with 16S sequence processing

OTUtabps3.csv -- OTU table for phyloseq object ps3, which is the basis of further subsetting and analysis

Taxatabps3.csv -- Taxa table for phyloseq object ps3

Psite Prev.csv -- Summary table of parasite prevalences by sampling period
