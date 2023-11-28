# ENVS-6452_Ass4_APP_Dada2
This is an amplicon processing pipeline for my assignment.

# Loading all the required libraries

library ("dada2")
library ("phyloseq")
library ("Biostrings")
library ("ggplot2")

# You can install these packages using the following commands: (if you don't have them installed)

install.packages("dada2")
install.packages("phyloseq")
install.packages("Biostrings")
install.packages("ggplot2")

# Setting working directory
setwd ("C:/Users/_Location for Directory_")

# Defining the working directory
path<- "C:/Users/_Location for Directory_"

# Checking files in the directory
list.files (path)

# If filenames are different, change them accordingly
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extracting sample names from filenames, assuming filenames have format: SAMPLENAME_XXX.fastq 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Visualizing the quality profile of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in the "filtered/" subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names,"_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=FALSE) 

head(out)

# Checking errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting errors
plotErrors(errF, nominalQ = TRUE)


# Looking for unique forward and reverse reads we have in the sample
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
 

# Merging pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Constructing a sequence table named seqtab
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspecting distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# removing chimeras (after looking at perfect sequences)
seqtab.nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE, verbose=TRUE) 
dim(seqtab.nochim)

# summarizing 
sum(seqtab.nochim)/sum(seqtab)

# tracking reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Comparing with dataset downloaded
taxa<- "C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/silva_nr99_v138.1_train_set.fa.gz"
       
# Assigning taxonomy
taxa <-assignTaxonomy(seqtab.nochim,C:Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/silva_nr99_v138.1_train_set.fa.gz, multithread=FALSE)

# Removing sequence rownames for display only
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
       
# Saving as csv file
       
write.csv (taxa, file="C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/taxa.csv")
write.csv (seqtab.nochim, file="C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/seqtab_nochim.csv")
       
# Making objects and opening files
taxa<- read.csv(file="C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/taxa.csv", sep=',',row.names=1)
       
seqtab_nochim <- read.csv(file="C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/seqtab_nochim.csv", sep=',',row.names=1)
       
# make as matrices for phyloseq
setab.nochim <- as.matrix(seqtab.nochim)
taxa<- as.matrix (taxa)

# transpose the file
flipped_seqtab.nochim<- as.data.frame (t (seqtab.nochim))

# Look for dimensions
dim (flipped_seqtab.nochim)

# Merging files
mergedfiles<- cbind (flipped_seqtab.nochim, taxa)

# Saving the file
write.csv (mergedfiles, file="C:/Users/BHUMIKA/Desktop/FALL23/BIOINFORMATICS/Assignment3/mergedfiles.csv")

# Making a custom dataframe
samples.out<- rownames (seqtab.nochim)
samdf <- data.frame (samples.out)
rownames (samdf)<- samples.out

# Handoff data to phyloseq

ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows = FALSE),
sample_data (samdf),
tax_table (taxa))


# Using biostrings
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


# Plotting graphs
plot_bar (ps)
plot_bar (ps , fill = "Phylum")
# making object p
p<- plot_bar (ps)
#add ggplot over it


# change to dataframe for plooting graph with ggplot2
ps.table <-psmelt(ps)

# use factor$ command-take column phylum 
ps.table$Phylum <-factor(ps.table$Phylum)

ggplot (data=ps.table, mapping=aes (x = Sample, y= Abundance))+ geom_bar(aes(fill=Phylum), stat= "identity",position="stack")

# Transforming data to relative abundance
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# Convert the phyloseq object ps.prop to a melted data frame
ps.proptable<- psmelt(ps.prop)

# Creating Object Phylum
ps.proptable$Phylum <-factor(ps.proptable$Phylum)

# ploting graph
ggplot (data=ps.proptable, mapping=aes (x = Sample, y= Abundance))+ geom_bar(aes(fill=Phylum), stat= "identity",position="stack")
