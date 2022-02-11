#Code to process data from M.L. Doolin, S.B. Weinstein, M.D. Dearing —— 
#Neotoma albigula parasite-microbiome study, Journal of Parasitology
#Submitted for publication February 2022. 
#16S sequence processing, using DADA2 and phyloseq pipeline.

##R version 4.0.2 (2020-06-22)
##Platform: x86_64-apple-darwin17.0 (64-bit)
##Running under: macOS Mojave 10.14.3   

############## Import, filter, trim, and make first phyloseq objects #####

#This portion of this script is adapted from the following tutorials: 
# https://benjjneb.github.io/dada2/tutorial.html
# https://joey711.github.io/phyloseq/
#Modified by M.L. Doolin 2029-2022

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



##### Locate sequence files #####

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
#Look at the forward and reverse read qualities.
plotQualityProfile(fnFs[60:69])
plotQualityProfile(fnRs[1:9])




##### Look at sequences, first filtering, and remove primers #####

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

# Make sequence path for cutadapt-trimmed files. 
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


# Run Cutadapt to remove the primers. 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "--discard-untrimmed", "--minimum-length", 8,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


##### Work with files to filter after primer removal #####


#Make new files with all sequences with primers removed.
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz")) 
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz"))
View(cutFs) #Make sure your files are there. 

#Put these files where they should be in the output folder.
cutFs <- file.path(path.cut, cutFs)
cutRs <- file.path(path.cut, cutRs)

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

# Filtering and Trimming of Reads for Quality Control 
# define the filepaths for the next round of filtering, 
# these parameters do not trim only filter: 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
f_path <- file.path(sequencepath, "filtered") 
if(!file_test("-d", f_path)) dir.create(f_path)
filtFs <- file.path(f_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(f_path, paste0(sample.names, "_R_filt.fastq.gz"))


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



##### Infer sequence variants with dereplication #####

#Dereplication: combining identical sequence reads in "unique sequences" with 
#corresponding abundance
# sample inference will take longer amounts of time with greater 
# depth and increased sample #'s 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# call the derep-class objects by their sample names: 
names(derepFs) <- sample.names
names(derepRs) <- sample.names


##### Learn errors #####

# Estimate error rates (sequencing error vs biological variation)
# This may take some time: 
errF <- learnErrors(filtFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose = TRUE)

# plot errors: 
# black line should follow sequence data points
# indicating computer effectively learned to predict errors 
plotErrors(errF)
plotErrors(errR)


##### Infer true seq variants and merge forward and reverse reads #####

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

table(nchar(getSequences(seqtab)))




##### Remove chimeras and look at read depth through filtering steps #####

# remove artificial chimeras: 
# Chimeras are: artifact sequences generated during the sequencing process 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, 
                                    verbose=TRUE)
dim(seqtab.nochim) #How many unique ASVs do you have after chimera removal?

sum(seqtab.nochim)/sum(seqtab) #Ratio of sequences retained


#Make a table where you're looking at all of the read depths across all steps. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)


##### Assign taxonomy #####

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

   ##### Optional next step assign species at 100% alignment to ASVs ##### 
# if interested: follow instructions here: https://benjjneb.github.io/dada2/tutorial.html

Taxa <- addSpecies(Taxa, "~/Desktop/Working_R_Docs/silva_species_assignment_v138.fa") 

# Next mandatory step: Inspect the taxa table: 
taxa.print <- Taxa 
rownames(taxa.print) <- NULL
head(taxa.print)



##### Making a Microbial Phylogenetic Tree for phylogeny-based analyses ################
# using the Decipher package in R and FastTree in terminal (bash script)
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
# Multiple Sequence Alignment via DECIPHER (Better for Large Datasets)
mult <- AlignSeqs(DNAStringSet(seqs), anchor=NA) 
# write out alignment file in FASTA format for FastTree in Terminal 
# make sure to use writeFASTA (shortreads) and not write.FASTA (ape) 
writeFasta(mult, "~/Desktop/Working_R_Docs/msa2.fasta")

   ##### Terminal ommands to download FastTree and generate tree ##### 
# Installation: 
#   conda install -c bioconda fasttree
#   conda install -c bioconda/label/cf201901 fasttree
# Change the working directory to the desktop where the msa file is: 
#   cd Desktop/Working_R_Docs
# Upload the msa file and create tree: 
#   FastTree -gtr -nt < msa2.fasta > MLtree2.nwk


# Read tree made in FastTree back into R 
tree <- read_tree("~/Desktop/Working_R_Docs/MLtree2.nwk")
head(taxa_names(tree)) # make sure names of sequences were not altered 

sample.names
plot_tree(tree, ladderize="left", label.tips="OTU") #Take a look at your tree.


##### Read in metadata and set up phyloseq objects ########

   ##### Metadata upload #####
# uploading Metadata sheet: Make sure this is in an R acceptable format 
samdf <- read.csv("AllSampsAvg_Metadat.csv")
head(rownames(seqtab.nochim)) #row names will need to be in this format
# this will propagate the sequence IDs to match the metadata sheet 
samdf1<-cbind(rownames(seqtab.nochim), samdf)
samdf2 <- data.frame(samdf1, row.names = 1)
#View this to make sure that everything matches up
View(samdf2)


   ##### Create phyloseq object #####
# now combine the metadata, Taxonomy, ASVs, and tree into a phyloseq object 
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(Taxa),
                sample_data(samdf2), phy_tree(tree))

sample_data(ps1) # call metadata to make sure you're in there correctly.
get_taxa_unique(ps1, "Kingdom") # look at what kingdoms are in the data set: 



   ##### Trimming and subsetting phyloseq objects #####

# remove reads that are non-bacteria, ie. mitochondria and chloroplasts from the dataset 
ps2 <- subset_taxa(ps1, Kingdom == "Bacteria" &  
                     Family  != "Mitochondria" & 
                     Class   != "Chloroplast" )
#Making the final phyloseq object to be subsetted and/or rarefied
#remove singletons from the datasets: (Large chance these are artifacts)
ps3<-prune_taxa(taxa_sums(ps2) > 1, ps2)

# make a DF with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps3))
smean <- mean(sample_sums(ps3))
smean 

#Save phyloseq object that will be subsetted for future quick access. 
saveRDS(ps3, "~/Desktop/Working_R_Docs/ps3RDS.rds")
ps3 <- readRDS("ps3RDS.rds")

#Subset physeq object to have only feces, pinworm and uninfected, regular
# diet (not OX spike diet)
ps4 <- subset_samples(ps3, gut_region=="feces" & Parasites != "Coccidian" &
                        Parasites != "Coinfection" & Trial=="reg")
ps4 <- prune_taxa(taxa_sums(ps4) >= 1, ps4)
#This is the physeq object I'll use for May/October regular diet analysis.


#Make a physeq object that has only May fecal samples. 
psmay <- subset_samples(ps4, month=="may")
psmay <- prune_taxa(taxa_sums(psmay) >= 1, psmay)
sd(sample_sums(psmay))
##mean for psmay is 18978.67
##sd is 2016.473


#Make a physeq object that only has October pinworm and uninfected, regular diet. 
psoct <- subset_samples(ps4, month=="oct")
psoct <- prune_taxa(taxa_sums(psoct) >= 1, psoct)
sd(sample_sums(psoct))
## mean for psoct is 12566.7
## sd is 5729.401

#Make a physeq object with just the OX spike diet samples, pinworm and uninfected.
psspike <- subset_samples(ps3, Trial=='spike' & Parasites != "Coccidian" &
                            Parasites != "Coinfection")
sample_names(psspike)
psspike<-prune_taxa(taxa_sums(psspike) >= 1, psspike)
max(sample_sums(psspike))
##mean read count is 10525.89
##sd is 4100.456
##range from 3100-17496 reads




   ##### Example exporting taxa and OTU tables as csv's #####

#Exporting as Taxa and OTU table with all samples. 
OTU1 = as(otu_table(ps3), "matrix")
if(taxa_are_rows(ps3)){OTU1 <- t(OTU1)}
OTUdf.norare = as.data.frame(OTU1)
write.csv(OTUdf.norare,"~/Desktop/working_R_Docs/OTUtabps3.csv")

Taxa1 = as(tax_table(ps3), "matrix")
if(taxa_are_rows(ps3)){Taxa1 <- t(Taxa1)}
Taxadf.norare = as.data.frame(Taxa1)
write.csv(Taxadf.norare,"~/Desktop/Working_R_Docs/Taxatabps3.csv")











############## Working with phyloseq objects for diversity metrics #####

library(tidyverse)
library(ggplot2)
library(ape)
library(RAM)
library(phyloseq)
library(dplyr)
library(vegan)

#General processing pipeline found here
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/composition-plots.html
#Modified by M.L. Doolin 2019-2022

##### Create sample data df outside of phyloseq #####

dftmp <- data.frame(estimate_richness(ps4, measures = "Observed"),
                    sum = sample_sums(ps4), estimate_richness(ps4, measures="Shannon"),
                    btools::estimate_pd(ps4))  #Last part adds phylo div and spec richness 
dfps4 <- cbind(sample_data(ps4),dftmp)


##### Richness and alpha diversity metrics and visualization #####


   ##### View richness, compare across parasitism groups via t-test #####
#Visualize observed and alpha diversity metrics
plot_richness(ps4, x='Parasites', color="Rat", shape="Parasites", measures=c('chao1', "Observed"))
?plot_richness


#Calculate t-test to compare alpha div across groups.
results = estimate_richness(ps4, measures="Shannon")
#"results" is now a table of all of the alpha div metrics for all samples. Or can specify:
#measures = c("Observed", 'Shannon') to include those metrics.
d = sample_data(ps4)

# calculate t-test
a = results[d[,'Parasites'] == 'Nematode',]
b = results[d[,'Parasites'] == 'Uninfected',]
t.test(b,a)

# Is there a linear relationship between observed richness and read count?
# Do this with unrarefied data to see what's up.
ggplot(dfps4, aes(sum, Observed,  color=Parasites))+
  geom_point(aes(shape=mp),size = 4, alpha = 0.75) +
  labs(x = "Filtered reads", y = "Observed ASVs")+ 
  theme_minimal() + theme(text=element_text(size=20)) +
  geom_text(mapping = aes(label = Rat), size = 3, vjust = 1.5) 



   ##### Stepwise model choice for modeling alpha diversity #####

#Choosing model fit to do more comprehensive/appropriate stats on ADiv
#From https://danstich.github.io/stich/classes/BIOL217/index.html
#Tutorial on model selection.
null=lm(Observed~1, data=dfps4)
full=lm(Observed~Parasites + month + gdna.conc, data=dfps4)
step(null, scope=list(lower=null, upper=full), direction ='backward')



   ##### Look at stacked bar plots with microbiome package #####
#Tutorial from: https://microbiome.github.io/tutorials/Composition.html
# using microbiome package.
library(microbiome)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(RColorBrewer)

#Also, for some of the stuff later on in this, need to use microbiomeutilities package
install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")# Make sure we use functions from correct package
transform <- microbiome::transform

#Make sure that whatever data you're going to average by is as.factor and not character.
sample_data(ps4)$P1 <- as.factor(sample_data(ps4)$Parasites)
sample_data(ps4)$m1 <- as.factor(sample_data(ps4)$month)
sample_data(ps4)$mp <- paste0(sample_data(ps4)$month, "", sample_data(ps4)$Parasites)

# Merge rare taxa to speed up examples
pseq <- transform(ps4, "compositional")
pseq1 <- aggregate_rare(pseq, level = "Family", detection = 0.01, prevalence = 0.1)
# This last line was making rare taxa merge together. If I do this at Genus level,
#  I will get almost all unassigned.
# Detection would mean the relative abundance in the sample, and prevalence would be
# the proportion of samples in which the taxon is found.


#To get every sample with a bar and sort the families in the order you want:
p <- microbiome::plot_composition(pseq1,
   taxonomic.level = "Family",
   x.label = "Parasites", 
   sample.sort = "mp",
   otu.sort=c("Muribaculaceae", "Lactobacillaceae", "Lachnospiraceae", "Erysipelotrichaceae", "Ruminococcaceae", 
                          "Monoglobaceae", "Oscillospiraceae", "Bifidobacteriaceae", "Bacteroidales_RF16_group", "Other")
   ) +
   guides(fill = guide_legend(ncol = 1)) +
   scale_y_percent() + 
   labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance ps4",
       subtitle = "subtitle", 
       caption = "Parasitism by month") + 
   theme_ipsum() + scale_fill_brewer("Family", palette = "RdBu") + ggpubr::rotate_x_text() 
print(p)


#To average bars by group
p <- plot_composition(pseq, otu.sort= "abundance",
                      average_by = "m1", transform = "compositional") +
  theme_ipsum() + scale_fill_brewer(palette = "PRGn") +
  ##theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("ps4 Stacked by month")
p +coord_flip()



##### Beta diversity metrics #####
library(vegan) #make sure this is loaded. 

   ##### Ordinate and visualize #####

#Use vegan to create ordinations.
#supported ord methods: ord_meths = "DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")

psmay1 <- subset_samples(psmay, Rat != "CV559") #to see May data without the one
# weird animal that resembles October animals.

# Transform data to proportions as appropriate for B-diversity analysis 
pstrans1 <- transform_sample_counts(psoct, function(otu) otu/sum(otu)) 
#And now to look at these with NO COCCIDIANS OR COINFECTIONS. 
ps.ordi <- ordinate(pstrans1, "PCoA", "bray")
ps.ordi2 <- ordinate(pstrans1, "NMDS", "bray")
ps.ordi3 <- ordinate(pstrans1, "PCoA", "unifrac", weighted = TRUE)
ps.ordi4 <- ordinate (pstrans1, "PCoA", "unifrac", weighted = FALSE)
ps.ordi5 <- ordinate(pstrans1, "NMDS", "unifrac", weighted = TRUE)
ps.ordi6 <- ordinate(pstrans1, "NMDS", "unifrac", weighted = FALSE)
plot_scree(ps.ordi4)


#Plotting the various visualizations
plot_ordination(pstrans1, ps.ordi, color="Parasites") +
  geom_point(aes(shape=Parasites), size = 10) + scale_shape_manual(values=c(16, 17))+
  theme(text=element_text(size=20)) + theme_bw() +
  scale_colour_manual(values=c("#000000", "#000000")) +
  #geom_text(mapping = aes(label = Rat), size = 4, vjust = 2) +
  ggtitle("psoct bray PCoA") + stat_ellipse(level=0.85, type="norm") 


#Trying something a little different with plotting, a split plot.
plot_ordination(pstrans1, ps.ordi, type="split", color = "Phylum", shape="Parasites", 
                label="Rat") +
  theme(text=element_text(size=20)) +
  ggtitle("All October Bray-Curtis PCoA")



   ##### Check out the dispersion of the data to see which stats to use #####

#Are the data suitable for parametric testing (i.e. PERMANOVA)? Look 
# at data dispersion.
?betadisper
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
x <- distance(ps4, "bray")
y <- as(sample_data(ps4), "data.frame")
groups <- y[["Parasites"]]
groups
mod <- betadisper(x, groups)
anova(mod)
plot(mod)


   ##### Run an ANOSIM (non-parametric) #####

#Anosim by Trial, for psnococc Bray-curtis
x = get_variable(pstrans1, "Parasites")
y = get_variable(pstrans1, "month")
ANOSIM = anosim(phyloseq::distance(pstrans1, "unifrac", ), x)
summary(ANOSIM) 



   ##### Run PERMANOVA (parametric) #####

metadata <- as(sample_data(ps4), "data.frame")
adonis2(phyloseq::distance(ps4, method="bray") ~  Parasites + month + Parasites*month,
        data = metadata)


##### Looking for specific taxa in OTU table######

# Example: subset to only Eggerthellaceae in dataset 
psegg = subset_taxa(psrare, Family == "Eggerthellaceae")
plot_bar(psegg, x="Date", fill="Genus") +xlab(" ") +
  ggtitle("Abundance of Eggerthellaceae across Samples")

p = plot_bar(psegg, x="Date", fill ="Genus")  
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle("Abundance of Eggerthellaceae")

# removing low abundance reads for easier plotting 
psegg1 <- prune_taxa(taxa_sums(psegg) > 50, psegg)

plot_tree(psegg1, color = "Date", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.3, ladderize = TRUE)

# create a sequence table for eggerthellaceae 
TaxaEgg = as(tax_table(psegg), "matrix")
if(taxa_are_rows(psegg)){TaxaEgg <- t(Taxaegg)}
Taxadf.psegg = as.data.frame(TaxaEgg)
write.csv(Taxadf.psegg,"~/Desktop/Eggerthellaceae.csv")


##### Differential abundance analysis with DESeq2 #####

library(DESeq2)
library(microbiome)

#Tutorial on DESeq2 from phyloseq.
#https://joey711.github.io/phyloseq-extensions/DESeq2.html


#Need to make sure you're doing pairwise comparisons. Make sure data structure
# and sample variables are logically set out to be able to do those comparisons.


#Using phyloseq_to_deseq2 to do the hard work and transform data into a form
# that is able to be analyzed by this program, which was made for RNA seq data
dsps = phyloseq_to_deseq2(ps4, ~ month)
dsps = DESeq(dsps, test="Wald", fitType="parametric")
res = results(dsps, cooksCutoff = FALSE)
summary(res)


#Shrink estimates to account for low read counts and other issues. 
#BiocManager::install("apeglm") #version 1.10.0
#install.packages("ahsr") #version 2.2.47
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="month_oct_vs_may", type="ashr")
# can use apeglm or ashr as type. Note that they give different LFC values
# that are several orders of magnitude different, even if total ASVs are the same.
summary(resLFC)

#Define your cut-offs for the adjusted p-value for output.
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps4)[rownames(sigtab), ], "matrix"))
dim(sigtab)


#Now to look at the differentially enriched taxa in ggplot
#library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

#Going to do Family instead.
#Genus isn't all that informative just because NA is the biggest group,
# likely from Muribaculaceae. 
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#Plot it
X <- ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ggtitle("DESeq2 ps4 ~ month w/lfcShrink")
X
View(sigtab)

#Write out results in spreadsheet of diff abundant ASVs so you can reference later
write.csv(sigtab, "~/Desktop/Working_R_Docs/ps4DESeq.csv")







