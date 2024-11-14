library("Rsubread")
library("dplyr")

# Take in scratch directory from sbatch script and make file names
args = commandArgs(trailingOnly=TRUE)
gtf_pr <- paste0(args[1],"/regions/hg38_refseq_pauseratio_pauseregion.gtf")
gtf_rog <- paste0(args[1],"/regions/hg38_refseq_pauseratio_restofgene.gtf")
bam_input <- paste0(args[1],"/",args[2])
outfile_pr <- paste0(args[1],"/results/pause_region_counts.txt")
outfile_rog <- paste0(args[1],"/results/restofgene_counts.txt")

# Put file and sample lists into variable
inputData <- read.table(bam_input, sep="\t")
fileList <- as.character(t(inputData["V1"]))
sampleList <- as.character(t(inputData["V2"]))

gtf_table <- read.table(gtf_pr)
colnames(gtf_table) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9",
                         "GeneID","V11","V12","TranscriptID","V14","V15","V16")

# Run featureCounts on these BAM files
coverage_pr <- featureCounts(files=fileList,
                    annot.ext=gtf_pr,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="gene_length",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=FALSE,
                    allowMultiOverlap=TRUE,
                    countMultiMappingReads=FALSE,
                    fracOverlap=1,
                    isPairedEnd=FALSE,
                    strandSpecific=1,
                    nthreads=32)

# Run featureCounts on these BAM files
coverage_rog <- featureCounts(files=fileList,
                    annot.ext=gtf_rog,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="gene_length",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=FALSE,
                    allowMultiOverlap=TRUE,
                    countMultiMappingReads=FALSE,
                    fracOverlap=1,
                    isPairedEnd=FALSE,
                    strandSpecific=1,
                    nthreads=32)


# Write out annotation and count data as tab delimited txt files
colnames(coverage_pr$counts) <- sampleList
colnames(coverage_rog$counts) <- sampleList

pr_out <- data.frame(GeneID=gtf_table["GeneID"],
                     TranscriptID=gtf_table["TranscriptID"],
                     Length=coverage_pr$annotation["Length"],
                     coverage_pr$counts,
                     stringsAsFactors=False)

rog_out <- data.frame(GeneID=gtf_table["GeneID"],
                      TranscriptID=gtf_table["TranscriptID"],
                      Length=coverage_rog$annotation["Length"],
                      coverage_rog$counts,
                      stringsAsFactors=False)

write.table(pr_out,file=outfile_pr,row.names=FALSE,sep='\t',quote=FALSE)
write.table(rog_out,file=outfile_rog,row.names=FALSE,sep='\t',quote=FALSE)
