## RNA-seq Analysis with DESeq2
## Written by Jan Forster, modified by Natasha Picciani 
## Available at https://github.com/snakemake-workflows/rna-seq-star-deseq2

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)

# Load counts as matrix
mycounts<- as.matrix(read.table(snakemake@input[["counts"]],header=T,row.names=1,comment.char='#'))
mycounts <- mycounts[,!colnames(mycounts) %in% c("Chr","Start","End","Strand","Length")] #remove unwanted columns from FeatureCounts output file
colnames(mycounts) <- gsub("results.star.mapping.|.Aligned.out.sam","",colnames(mycounts)) #remove unwanted patterns from sample names
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
mode(mycounts) <- "numeric"

# Check if formatting is correct
all(rownames(samples) %in% colnames(mycounts)) # must be TRUE
all(rownames(samples) == colnames(mycounts)) # must be TRUE
#mycounts<-mycounts[,rownames(samples)] # run if previous check returns FALSE or a WARNING message

# Create a DESEqData Object
dds <- DESeqDataSetFromMatrix(countData=mycounts, colData=coldata, design=as.formula(snakemake@params[["model"]]))

# Run the DESeq pipeline
dds <- DESeq(dds)

# Save dds object as RDS
saveRDS(dds, file=snakemake@output[[1]])

# Export normalized counts to a csv file
norm.counts <- counts(dds, normalized=TRUE) 
write.csv(norm.counts, file=snakemake@output[[2]])