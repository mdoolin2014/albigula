# This is being run again my M. Doolin for Castle Valley microbiome data
# Adapted from my previous version that was used to run May 2019 CV data and from
# Whitewater 16S Sequence Analysis by Dylan Klure.
# Last updated 7/9/20
# this script is adapted from the following tutorials: 
# https://benjjneb.github.io/dada2/tutorial.html
# https://joey711.github.io/phyloseq/
# Make sure all packages and R is up-to-date first! 
###In this version, looking at all May and October sequences including spike and 
### regular diet animals together. 

#Updated R in July 2020, running 
##R version 4.0.2 (2020-06-22)
##Platform: x86_64-apple-darwin17.0 (64-bit)
##Running under: macOS Mojave 10.14.3   

#Other potential good tutorials on analyzing data:
#Whole process: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
#Statistical analysis: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/



##### Uninstalling and reinstalling all packages and BiocManager #####
##Having trouble with packages actually working after updating R to 4.0.2. Uninstalled and reinstalled
## R and RStudio. Didn't fucking work. Now going to have to uninstall all Bioconductor 
## packages.
#from tutorial here: https://www.r-bloggers.com/how-to-remove-all-user-installed-packages-in-r/
ip <- as.data.frame(installed.packages())
head(ip)
ip <- subset(ip, !grepl("MRO", ip$LibPath))
ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
path.lib <- unique(ip$LibPath)
pkgs.to.remove <- ip[,1]
head(pkgs.to.remove)
pkgs.to.remove
#remove the packages
sapply(pkgs.to.remove, remove.packages, lib = path.lib)


chooseCRANmirror()
install.packages("BiocManager", version="3.11")
BiocManager::install(version="3.11")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version="3.11")
BiocManager::install(version="3.11")
BiocManager::install("dada2", version="3.11")
BiocManager::install("Biostrings", version="3.11")
BiocManager::install("BiocParallel", version="3.11")
BiocManager::install("ShortRead", version="3.11")
BiocManager::install("SummarizedExperiment", version="3.11")
BiocManager::install("metagenomeSeq", version="3.11")


##### Coming back to the 16S script #####



#required packages: 
library(Rcpp)
library(ggplot2)
library(ape)
library(gridExtra)
library(SummarizedExperiment)
library(dada2) 
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(RSQLite)
library(ShortRead)
library(Biostrings)


# verify everything is running
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


#Currently using sequences that are from sequencing runs in June 2019 and January 2020
# from the Univ of Illinois at Chicago sequencing core.
# load fastq files into R
# the files need to be unzipped and all in the same folder 
# must have both the forward and reverse reads 
setwd("/Users/mdoolin/Desktop/Working_R_Docs/")
sequencepath <- "/Users/mdoolin/Desktop/Working_R_Docs/All_albigula_seqs/"


# Make sure the files are actually there: 
list.files(sequencepath)

# great, now: Sort to ensure fwd and reverse reads are in same order
fnFs <- sort(list.files(sequencepath, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(sequencepath, pattern="_R2_001.fastq.gz"))

# specify full path to the fnFs and fnRs
fnFs <- file.path(sequencepath, fnFs)
fnRs <- file.path(sequencepath, fnRs)
fnFs[42:45]
fnRs[1:3]

plotQualityProfile(fnFs[60:69])
plotQualityProfile(fnRs[1:9])




##### removing primers before further processing #################
# CutAdapt: removing primers before further processing in dada2
# tell R where to find cutadapt - must download in terminal before this step
# https://cutadapt.readthedocs.io/en/stable/installation.html'
# code from: https://benjjneb.github.io/dada2/ITS_workflow.html
cutadapt <- "/Users/mdoolin/miniconda3/envs/qiime2-2019.7/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R
## version 2.5

# New 16S V4 primer sequences: 
# http://www.earthmicrobiome.org/protocols-and-standards/16s/
FWD <- "GTGYCAGCMGCCGCGGTAA" #  fwd, 515F updated 
REV <- "GGACTACNVGGGTWTCTAAT" # rev, 806R updated 

# to verify these primers are actually in the dataset: 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# remove sequences with ambiguous nucleotides: 
fnFs.filtN <- file.path(sequencepath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(sequencepath, "filtN", basename(fnRs))


filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, compress=TRUE, maxN = 0, 
              multithread= TRUE, verbose = TRUE)


# Counting the # of primers that appear in the fwd and rev reads
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#May 2019 numbers
#                 Forward Complement Reverse RevComp
#FWD.ForwardReads   73849          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads   72792          0       0       0

#Again, this is weirdly happening that only the first part of the sequences in the list,
# in this case, the no oxalate spike data, are showing up. But hopefully, like last time,
# This will work.
#All October and May sequences.
##                 Forward Complement Reverse RevComp
##FWD.ForwardReads   73849          0       0       0
##FWD.ReverseReads       0          0       0       0
##REV.ForwardReads       0          0       0       0
##REV.ReverseReads   72792          0       0       0



#Think about upping the quality truncation score from 10 to 30. 

# cutadapt portion: Removal of primers 
path.cut <- file.path(sequencepath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 



# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "--discard-untrimmed", "--minimum-length", 8,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
save.image("~/Desktop/Working_R_Docs/working.RData")
#Maggie, PROBLEM August 2019: All sequences are being removed at the cutadapt step! Need to figure
#out what is going on here, and why none of them are retained. 
#
##Resolved the issue by just rezipping all files. Unzipped files did not process correctly.
##For future, always keep all files zipped. 

# check that the primers were actually removed. Should all be 0: 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
#Great, all 0!

#Make new files with all sequences with primers removed.
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz")) 
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz"))
head(cutFs)
head(cutRs)
#Ok, those look good.
View(cutFs)
#Put these files where they should be in the output folder.
cutFs <- file.path(path.cut, cutFs)
cutRs <- file.path(path.cut, cutRs)
cutFs[1:3]
cutRs[1:3]
#Ok, now there are those files there, in the right spot. Now, let's get samples names.


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
View(sample.names)

# check quality of reads after primer removal: 
plotQualityProfile(fnFs[1:15]) # before primer removal
plotQualityProfile(cutFs[1:15]) # after primer removal 
plotQualityProfile(fnRs[1:3])
plotQualityProfile(cutRs[40:48]) # reverse after primer removal 
#Things are looking good. Quality should be above 20, and are basically all above
# 30 for the whole length.


save.image("~/Desktop/Working_R_Docs/working.RData")
##### Read Filtering and Trimming ##########################
# Filtering and Trimming of Reads for Quality Control 
# define the filtered file path, these parameters do not trim only filter: 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

f_path <- file.path(sequencepath, "filtered") 
if(!file_test("-d", f_path)) dir.create(f_path)
filtFs <- file.path(f_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(f_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Defining filtering parameters:  
# max of 2 expected errors per reads (Edgar and Flyvbjerg 2015)
# filter out reads < 50 bp 
# minimum quality score = 10 s
# rm.phix: Illumina uses a small DNA virus as a positive sequencing control, this step
# filters this viral DNA out of the dataset if it is present. 

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN=0, maxEE=c(2,5), minLen = 50, rm.phix=TRUE, minQ = 20, truncQ = 0,  
                     compress=TRUE, verbose = TRUE, multithread=TRUE, matchIDs = TRUE) 

head(out)
View(out)
#This is where I could see average reads across samples before filtering. 

# plot quality after filtering 
# make sure quality drop improves  
plotQualityProfile(filtFs[1:5])
plotQualityProfile(filtRs[1:9])

# inferring sequence variants 
#Dereplication: combining identical sequence reads in "unique sequences" with 
#corresponding abundance
# sample inference will take longer amounts of time with greater 
# depth and increased sample #'s 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# call the derep-class objects by their sample names: 
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# estimate error rates (sequencing error vs biological variation)
# heads up this again may take some time: 
errF <- learnErrors(filtFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose = TRUE)

# plot errors: 
# black line should follow sequence data points
# indicating computer effectively learned to predict errors 
plotErrors(errF)
plotErrors(errR)

# I recommend saving the work space after each time-consuming step 
save.image("~/Desktop/Working_R_Docs/uptoerrors.RData")

# now use the built algorithm to infer true sequence variants (at 100% identity)
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# visualize the dada class object after denoising: 
dadaFs[(1)]

# merge the forward and reverse reads  
# default minimum overlap = 12 bp 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE) 

#inspect the merges: 
head(mergers[(1)])

# make a sequence table and inspect it: 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
##[1]   74 6058 for all Oct (reg and spike) and May data, all run together
table(nchar(getSequences(seqtab)))

#May 2019
## 189  252  253  254  255 
##   1 1788 1736   21   17 

#All October seqs, regular and spiked. 
##178  189  251  252  253  254  255 
##  1    1    1 5530 3359   48   14 

#All Oct and May
##178  252  253  254  255 
##  1 3308 2706   30   13




# remove artificial chimeras: 
# Chimeras are: artifact sequences generated during the sequencing process 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, 
                                    verbose=TRUE)
dim(seqtab.nochim)
##[1]   74 2041

sum(seqtab.nochim)/sum(seqtab)
## 0.9550072 MD on 7/10/20 for All Oct and May seqs


#Make a table where you're looking at all of the read depths across all steps. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)


##### Optional stuff that I'm choosing to skip bc not working yet.######
#This part isn't working. Skipping this. Not the right version of R
# Optional step: Generate a SamplexSequence Individual FASTA for SILVang
install.packages("amplicR")
library(amplicR)
#amplicR is only for more recent versions of R. Mine is still 3.5.3
install.packages("remotes")
library(remotes)
table2fasta2(seqtab.nochim, seq.list = NULL, dir.out = "~/Desktop/SILVang", verbose = TRUE)
##### Back to the rest of the prep steps #####



#assign taxonomy: 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim))

# optional: if single sample: remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)


# set up the matrix for taxonomy assignment: 
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



#taxa training step: 
#Trying to use v138, let's see what happens this time...it worked!
Taxa <- assignTaxonomy(seqtab.nochim, 
          "/Users/mdoolin/Desktop/Working_R_Docs/silva_nr_v138_train_set.fa", multithread=TRUE, 
          verbose=TRUE)



# Optional next step assign species at 100% alignment to ASVs: 
# if interested: follow instructions here: https://benjjneb.github.io/dada2/tutorial.html

Taxa <- addSpecies(Taxa, "~/Desktop/Working_R_Docs/silva_species_assignment_v138.fa") 




# Next mandatory step: Inspect the taxa table: 
taxa.print <- Taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

# looks great, lets save the workspace 
save.image("~/Desktop/Working_R_Docs/uptotaxatable.RData") 

##### Making a Microbial Phylogenetic Tree. Do this. ################
# using the Decipher package in R and FastTree in the terminal
# code source + instructions: http://www.microbesonline.org/fasttree/
# this can be downloaded on a linux or Mac OS system and run in the terminal: 
# FastTree will rapdily generate a Maximum Likelihood Tree from alignments 
# packages: 
BiocManager::install("DECIPHER", version="3.11")
BiocManager::install("phyloseq", version="3.11")
library(DECIPHER)
library(phyloseq)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
# Multiple Sequence Alignment via DECIPHER (Better for Large Datasets:)
mult <- AlignSeqs(DNAStringSet(seqs), anchor=NA) 
# write out alignment file in FASTA format for FastTree in Terminal 
# make sure to use writeFASTA (shortreads) and not write.FASTA (ape) 
writeFasta(mult, "~/Desktop/Working_R_Docs/msa2.fasta")

## Terminal Commands for Downloading FastTree and generating ML Tree: 
# Installation: 
#   conda install -c bioconda fasttree
#   conda install -c bioconda/label/cf201901 fasttree
# Change the working directory to the desktop where the msa file is: 
#   cd Desktop/Working_R_Docs
# Upload the msa file and create ML tree: 
#   FastTree -gtr -nt < msa2.fasta > MLtree2.nwk


# Upload ML tree made in FastTree back into R 
tree <- read_tree("~/Desktop/Working_R_Docs/MLtree2.nwk")
head(taxa_names(tree)) # make sure names of sequences were not altered 

sample.names
plot_tree(tree, ladderize="left", label.tips="OTU")
#Lol, yikes this worked but it a tree of all of the microbial taxa in the whole 
# phyloseq object so is a crazy, illegible plot. Really what I would care to do
# is to make a UPGMA tree to look at the relationship between samples.
# Potential source for this is https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/


##### Setting up the first Phyloseq Object ########



# uploading Metadata Sheet: Make sure this is in an R acceptable format 
#Also make sure that the working directory is already set. 
samdf <- read.csv("AllSampsAvg_Metadat.csv")
head(rownames(seqtab.nochim)) #row names will need to be in this format
# this will propagate the sequence IDs to match the metadata sheet 
samdf1<-cbind(rownames(seqtab.nochim), samdf)
samdf2 <- data.frame(samdf1, row.names = 1)
#View this to make sure that everything matches up
View(samdf2)

# now combine the metadata, Taxonomy, ASVs, and tree into a phyloseq object 
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(Taxa),
                sample_data(samdf2), phy_tree(tree))


# Optional: if processing a single sample: 
##ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(Taxa)) 

# look at phylo object: 
ps1

sample_data(ps1)  

#All May and All October sequences, including all gut regions from May, and both diet
# types from October. 
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 2041 taxa and 74 samples ]
##sample_data() Sample Data:       [ 74 samples by 18 sample variables ]
##tax_table()   Taxonomy Table:    [ 2041 taxa by 6 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 2041 tips and 2039 internal nodes ]



# look at what kingdoms are in the data set: 
get_taxa_unique(ps1, "Kingdom")

# remove reads that are non-bacteria, ie. mitochondria and chloroplasts from the dataset 
ps2 <- subset_taxa(ps1, Kingdom == "Bacteria" &  
                     Family  != "Mitochondria" & 
                     Class   != "Chloroplast" )

ps2

# Now down to 1572 taxa across the 74 samples in All May and Oct data..
 

#Making the final phyloseq object to be subsetted and/or rarefied
#remove singletons from the datasets: (Large chance these are artifacts)
ps3<-prune_taxa(taxa_sums(ps2) > 1, ps2)

ps3
# for all May and Oct Castle Valley data, now down to 1513 taxa in ps3.  
# but this still includes 



# make a DF with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps3))
smean <- mean(sample_sums(ps3))
smean 
##13934.36 for all May and October. I decided to up the quality score when I filtered to 
# 20 instead of 10, so it cut out many more seqs than in previous analyses of these seqs. 

sd(sample_sums(ps3))
##6229.949

#Cool, save it to the working directory.
saveRDS(ps3, "~/Desktop/Working_R_Docs/ps3RDS.rds")
ps3 <- readRDS("ps3RDS.rds")

#Exporting as Taxa and OTU table with all samples. 
OTU1 = as(otu_table(ps3), "matrix")
if(taxa_are_rows(ps3)){OTU1 <- t(OTU1)}
OTUdf.norare = as.data.frame(OTU1)
write.csv(OTUdf.norare,"~/Desktop/Working_R_Docs/OTUtabps3.csv")

Taxa1 = as(tax_table(ps3), "matrix")
if(taxa_are_rows(ps3)){Taxa1 <- t(Taxa1)}
Taxadf.norare = as.data.frame(Taxa1)
write.csv(Taxadf.norare,"~/Desktop/Working_R_Docs/Taxatabps3.csv")

#Going back through in a final pass, and the exported OTU table is correct as the high 
# quality-filtered ps3. This includes all May samples from all gut sections, and all 
# October samples from both regular and spike diet. 


##### Rarefying and Subsetting phyloseq objects, saving as csv's #####
#Making phyloseq objects by subsetting data by parasite and diet type

#RAREFYING WHOLE DATASET
psrare = rarefy_even_depth(ps3)
psrare <- prune_taxa(taxa_sums(psrare) >= 1, psrare)
psrare
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 1131 taxa and 74 samples ]
##sample_data() Sample Data:       [ 74 samples by 18 sample variables ]
##tax_table()   Taxonomy Table:    [ 1131 taxa by 6 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 1131 tips and 1129 internal nodes ]

mean(sample_sums(ps3))
##mean for ps3 is 13934.36
##30339 is max read number at ps3

mean(sample_sums(psrare))
##917
##All samples rarefied to 917 seqs. Yikes, that's pretty low. 
##There are still the samples from other gut sections in here, which is why 
## there are issues with this low, low count.


OTU1 = as(otu_table(psrare), "matrix")
if(taxa_are_rows(psrare)){OTU1 <- t(OTU1)}
OTUdf.psrare = as.data.frame(OTU1)
write.csv(OTUdf.psrare,"~/Desktop/Working_R_Docs/OTUtabAllrare.csv")

Taxa1 = as(tax_table(psrare), "matrix")
if(taxa_are_rows(psrare)){Taxa1 <- t(Taxa1)}
Taxadf.psrare = as.data.frame(Taxa1)
write.csv(Taxadf.psrare,"~/Desktop/Working_R_Docs/TaxatabAllrare.csv")
#Saving phyloseq objects separately as RDS files to more easily start work with them
# after closing R.
saveRDS(psrare, "~/Desktop/Working_R_Docs/AllrareRDS.rds")


#Going to make a ps4 physeq object - just October and May feces samples from the regular 
# diet trial, of only nematode and uninfected animals. This will be the main object I use
# for my analyses, because this is the main comparison I want to do. 
ps4 <- subset_samples(ps3, gut_region=="feces" & Parasites != "Coccidian" &
                        Parasites != "Coinfection" & Trial=="reg")
ps4 <- prune_taxa(taxa_sums(ps4) >= 1, ps4)
ps4
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 1033 taxa and 32 samples ]
##sample_data() Sample Data:       [ 32 samples by 18 sample variables ]
##tax_table()   Taxonomy Table:    [ 1033 taxa by 6 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 1033 tips and 1031 internal nodes ]

sample_variables(ps4)
sample_data(ps4)
mean(sample_sums(ps4))
## mean 14971.19
## sd 5613.269
saveRDS(ps4, "ps4RDS.rds")

OTU1 = as(otu_table(ps4), "matrix")
if(taxa_are_rows(ps4)){OTU1 <- t(OTU1)}
OTUdf.ps4 = as.data.frame(OTU1)
write.csv(OTUdf.ps4,"~/Desktop/Working_R_Docs/OTUtabps4.csv")

Taxa1 = as(tax_table(ps4), "matrix")
if(taxa_are_rows(ps4)){Taxa1 <- t(Taxa1)}
Taxadf.ps4 = as.data.frame(Taxa1)
write.csv(Taxadf.ps4,"~/Desktop/Working_R_Docs/Taxatabps4.csv")

setwd("~/Desktop/Working_R_Docs/")
ps4 <- readRDS("ps4RDS.rds")



#New object with October and May, feces only, both diet types, nematode and 
# uninfected only.
ps5 <- subset_samples(ps3, gut_region=="feces" & Parasites != "Coccidian" &
                        Parasites != "Coinfection")
ps5 <- prune_taxa(taxa_sums(ps5) >= 1, ps5)



######Other objects

#####Different objects that include October and May#####

#Take out all of the other gut regions so there's only poop left, but still has 
# all parasite infection statuses and spike trial.
pspoop <- subset_samples(ps3, gut_region=="feces")
pspoop <- prune_taxa(taxa_sums(pspoop) >= 1, pspoop)
saveRDS(pspoop, "pspoop.rds")


#Now making a phyloseq object with ONLY NEMATODE AND UNINFECTED INCLUDES SPIKE AND ALL
# GUT REGIONS.
#This is for analyses recommended, just looking at most prevalent psite.
psnococc <- subset_samples(ps3, Parasites != "Coccidian" &
                             Parasites != "Coinfection")
psnococc<-prune_taxa(taxa_sums(psnococc) >= 1, psnococc)
sampdat <- sample_data(psnococc)

#579 is the only one in the regular data and not in the spike data. She is in here, but 
# could take her out to have direct comparison if I wanted to do that.

OTU1 = as(otu_table(psspike), "matrix")
if(taxa_are_rows(psspike)){OTU1 <- t(OTU1)}
OTUdf.psspike = as.data.frame(OTU1)
write.csv(OTUdf.psspike,"~/Desktop/Working_R_Docs/OTUtabAllpsspike.csv")

Taxa1 = as(tax_table(psspike), "matrix")
if(taxa_are_rows(psspike)){Taxa1 <- t(Taxa1)}
Taxadf.psspike = as.data.frame(Taxa1)
write.csv(Taxadf.psspike,"~/Desktop/Working_R_Docs/Taxatabpsspike.csv")

saveRDS(psnococc, "~/Desktop/Working_R_Docs/Allnococc.rds")


#####Subsetting samples for May only#####

#Only May samples, from feces, with just nematode and no parasite
psmay <- subset_samples(ps4, month=="may")
psmay <- prune_taxa(taxa_sums(psmay) >= 1, psmay)

sd(sample_sums(psmay))
##mean for psmay is 18978.67
##sd is 2016.473

#But also want to be able to compare all of the gut sections looked at in May, 
# so make another object for that, from ps3.
psmay1 <- subset_samples(ps3, month=="may")
psmay1 <- prune_taxa(taxa_sums(psmay1) >= 1, psmay1)



#####Subsetting October samples#####

#October samples with nematode and uninfected, regular and spike diet
#This is what I'll use to compare October individuals to themselves, before and 
# after the dietary oxalate challenge. 
psalloct <- subset_samples(ps3, gut_region=="feces" & Parasites != "Coccidian" &
                             Parasites != "Coinfection" & month!="may")
psalloct <- prune_taxa(taxa_sums(psalloct) >= 1, psalloct)


#Only October samples, with just nematode and uninfected, without spike trial.
psoct <- subset_samples(ps4, month=="oct")
psoct <- prune_taxa(taxa_sums(psoct) >= 1, psoct)
sample_data(psoct)

sd(sample_sums(psoct))
## mean for psoct is 12566.7
## sd is 5729.401

#How about just making only the oxalate spike group a separate phyloseq object
#Not rarefied, still includes 
psspike <- subset_samples(ps3, Trial=='spike' & Parasites != "Coccidian" &
                            Parasites != "Coinfection")
sample_names(psspike)
psspike<-prune_taxa(taxa_sums(psspike) >= 1, psspike)

saveRDS(psspike, "~/Desktop/Working_R_Docs/OctSpikeRDS.rds")
max(sample_sums(psspike))
##mean read count is 10525.89
##sd is 4100.456
##range from 3100-17496 reads
#Wrote out OTU and taxa table to see what this super high SD is about...

#Then an October, with all parasitic infections, without spike trial.
sample_data(ps3)$Trial
psnorm <- subset_samples(ps3, Trial=='reg' & month=="oct")
psnorm <- prune_taxa(taxa_sums(psnorm) >= 1, psnorm)


#Saved all phyloseq objects as RDS files so that I can just load them instead
 #of loading the whole workspace every time.


##### Trying cumulative sum scaling (CSS) on my data using metagenomeSeq package#####
# http://www.metagenomics.wiki/tools/16s/norm/css#:~:text=Cumulative%20Sum%20Scaling%20(CSS)%20is,in%20total%20counts%20between%20samples.
#Starting with the line of code from above, pulling out the OTU table
OTU1 = as(otu_table(psnococc), "matrix")
OTU_read_count = as.data.frame(t(OTU1)) #transposed it and made into a dataframe
dim(OTU_read_count) 
## 1413   59   #  (this means 1413 taxa and 59 samples)
class(OTU_read_count)
## "data.frame" #Great, so now I can proceed, no error.

library(metagenomeSeq)  #if I haven't already done this.

SeqObj = newMRexperiment(OTU_read_count) 
# CSS normalization
SeqObj_CSS = cumNorm(SeqObj, p=cumNormStatFast(SeqObj))
# convert CSS normalized data into data.frame-formatted OTU table (log transformed data)
OTU_read_count_CSS = data.frame(MRcounts(SeqObj_CSS, norm=TRUE, log=TRUE))

#Ok this worked, so where to go from here...make this the OTU table for a new version of 
# psnococc
otu_table(psnococc1) <- matr

#Need to turn into a matrix to add to a phyloseq object.
matr <- data.matrix(OTU_read_count_CSS)

#And also need the taxa table and metadata.
Taxa1 = as(tax_table(psnococc), "matrix")
sampdat <- sample_data(psnococc)

#Merge the OTU and taxa tables into a phyloseq object
hmm <- phyloseq(matr, Taxa1, sampdat)

plot_tree(psnococc, color="Parasites", label.tips="", ladderize="left", plot.margin=0.3)



##### Generating BIOM Files from PS object #######################
# if you need BIOM files follow the below script 
# Exporting the new phyloseq object into biom files. 
install.packages("biomformat")
library(biomformat)
library(readxl)

setwd("/Users/mdoolin/Desktop/Working_R_Docs/")
Taxa <-as(tax_table(psrare),"matrix")
tax_cols <- colnames(Taxa)
Taxa <-as.data.frame(Taxa)
Taxa$taxonomy<-do.call(paste, c(Taxa[tax_cols], sep=";"))
for(co in tax_cols) Taxa[co]<-NULL
write.table(Taxa, "~/Desktop/TaxaTable.txt", quote=FALSE, col.names=FALSE, sep="\t")

OTUtable <-t(as(otu_table(psrare),"matrix")) # 't' to transform if taxa_are_rows=FALSE

# re-upload corrected file with concatenated taxa names 
correctedmerged <- read.csv("OTUtabOctMayrare.csv")

# create a biom file from OTU table 
otu_biomOctMayrare <-make_biom(data=correctedmerged)
write_biom(otu_biomOctMayrare,"OctMayrare.biom")
#Will need to do this for certain diff abundance analysis methods. 


##### Tax4Fun stuff, temporarily giving up on picrust2 and coming back to SILVA scripts#############
# KEGG Analysis via Tax4Fun in R
# tutorial: https://cran.r-project.org/web/packages/themetagenomics/vignettes/functional_prediction.html
# Last updated 3/22/2020 by Dylan Klure, then 3/25/20 by Maggie Doolin
# package source: https://cran.r-project.org/web/packages/themetagenomics/index.html
# session info: 
# R version 3.5.3 (2019-03-11)
library(themetagenomics)
library(dada2)
library(phyloseq)

############Tax4Fun tutorial stuff 
#####Tax4Fun can be used with SILVA taxonomy calls unlike PiCRUST2
# data needs to be in a list in this format: 
sample <- DAVID
sample$ABUND[1:2, 1:2] # this is the OTU table format 
# 00001 00002
# ERR531441  7232 44702
# ERR531442  1807     0

sample$TAX[1:2, 1:2] #taxa table format.
# Kingdom    Phylum         
# 00001 "Bacteria" "Bacteroidetes"
# 00002 "Bacteria" "Firmicutes"

sample$META[1:2, 1:2]
#               ID Day
# ERR531441 ERR531441 233
# ERR531442 ERR531442 169

#####Dylan's Tax4Fun script, and my alterations to his stuff.
# load in a phyloseq object and seperate out its contents 
##ps <- readRDS("C:/Users/Dylan/OneDrive/Desktop/WWprelim/RDS/ps.rds")

#My phyloseq object is already in the environment.
#Doing this for all data, and also looking at things with psnococc

OTU <- as.data.frame(psnococc@otu_table)
Taxa <- as.data.frame(psnococc@tax_table)
Meta <- as.data.frame(psnococc@sam_data)

# acquiring the KEGG SILVA reference database 
tmp <- tempdir()
download_ref(tmp,reference='silva_ko',overwrite=FALSE)

# running Tax4fun 
system.time(FUNCTIONS <- t4f(OTU,rows_are_taxa=FALSE,tax_table=Taxa,
                             reference_path=tmp,type='uproc',short=TRUE,
                             cn_normalize=FALSE,sample_normalize=FALSE,drop=FALSE,
                             verbose = TRUE))
#For psrare
##  user  system elapsed 
## 3.506   0.090   3.610 

#For psnococc All Oct seqs
##3.504   0.087   3.614 

# output is a data class object containing a matrix of gene counts, a list of 
# functional metadata that corresponds to the functional table, and a metadata for 
# the ran samples 

# explore the predicted functional outputs in Microbiome analyst: 
# need to write out data as txt files. 
T4Fprofile <- FUNCTIONS$fxn_table
T4Fprofile <- data.frame(t(FUNCTIONS$fxn_table))
write.csv(T4Fprofile, "Tax4FunProfileAllOctnococc.csv")
write.csv(Meta,"Tax4FunMETAAllOctnococc.csv")







