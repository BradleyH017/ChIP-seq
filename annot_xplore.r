# Bradley May 2022
# Annotating and exploring the ChIP-seq data
# Make a conda environment for this
# conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/ChIP-seq_env r-base r-signac bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase bioconductor-rtracklayer bioconductor-fgsea
# source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/ChIP-seq_env

# Set up
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/data")
#myPaths <- .libPaths()
#myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
#myPaths <- c(myPaths[2], myPaths[1])
#.libPaths(myPaths)


# Read in the narrowPeak results for each antibody and cell line
ab <- c("OCAT", "POU2F3")
line <- c("NCIH211", "NCIH526")

# Load the data into a named list
conds <- rep("", length(ab)*length(line))
for(a in seq_along(ab)){
  for(b in seq_along(line)){
    if(a==1){
      conds[a*b] <- paste(ab[a], line[b], sep = "_")
    } else {
      conds[a+b] <- paste(ab[a], line[b], sep = "_")
    }
  }
}

# First read in the narrowPeak file
# From https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
# First 6 columns of narrowPeak file are: chromosome, start, end, name, score, strand
# Last 4 are: signalValue(Overall enrichment for region), pValue, qValue, peak (point source for the peak; o-based offset from the start position),
# Following https://charlesjb.github.io/How_to_import_narrowPeak/
require(rtracklayer)
cond=1
narrowPeak_file <- c("OCAT_NCIH211/GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_peaks.narrowPeak")
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
gr_narrowPeak <- import(narrowPeak_file, format="BED", extraCols=extraCols_narrowPeak, genome="hg38")

# Now read in the annotated version of this
anno <- read.delim(paste(narrowPeak_file, "annot", sep = "."), header=T)
trans <- c("POU2F3", "ITPRID1", "BMX", "SH2D6", "SH2D7", "CHAT", "HTR3E", "TRPM5", "AZGP1", "OGDHL", "ACTG1P22", "AVIL", "KLK13", "IRAG2", "PSTPIP2", "PIK3CG", "PLCG2", "TAS1R3")
anno[anno$Gene.Name %in% trans,]

# Now combine with a dataframe version of the GRanges object
# This had to be initially loaded as granges as is converted from 0 base 1 base, so incorrect otherwise
df_narrowPeak <- data.frame(chr=gr_narrowPeak@seqnames, cbind(as.data.frame(gr_narrowPeak@ranges), as.data.frame(gr_narrowPeak@elementMetadata)))
colnames(anno)[which(colnames(anno) == "Chr")] <- "chr"
colnames(anno)[which(colnames(anno) == "Start")] <- "start"
colnames(anno)[which(colnames(anno) == "End")] <- "end"
want <- c("chr", "start", "end", "Gene.Name")
anno_merge <- anno[,colnames(anno) %in% want]
df_narrowPeak$peak_pos <- paste(df_narrowPeak$chr, ":", df_narrowPeak$start, "-", df_narrowPeak$end, sep = "")
anno_merge$peak_pos <- paste(anno_merge$chr, ":", anno_merge$start, "-", anno_merge$end, sep = "")
anno_merge <- anno_merge[,-c(1:3)]
np <- merge(df_narrowPeak, anno_merge, by = "peak_pos")
np$gene.peak_pos <- paste(np$Gene.Name, np$peak_pos, sep = "_")

# Now we have the gene names nearest the peaks called, can see if this includes any of out trans-eQTLs
# 1. Subset for those that are significant
np_sig <- np[np$qValue > 10^-0.01,] # all are
library(fgsea)
gene_list <- np_sig$signalValue
names(gene_list) <- np_sig$gene.peak_pos
# Need to use any of the trans genes that have peaks
trans_peak_pos <- np$gene.peak_pos
trans_peak_pos <- trans_peak_pos[grep(paste(trans, collapse = "|"), trans_peak_pos)]
trans_list <- vector("list")
trans_list[['rs3087967_trans-eQTL_targets']] <- trans_peak_pos
fgseaRes <- fgseaMultilevel(trans_list, gene_list, minSize = 0, maxSize = 500, scoreType = "pos", eps = 0)
# Plot this
fgseaRes_df <- as.data.frame(fgseaRes[,-8])
enr_anno <- paste("NES=", signif(fgseaRes_df$NES,2), "\np=", signif(fgseaRes_df$pval,2), sep = "")
enr_plot <- plotEnrichment(trans_list[[1]], gene_list) +
  theme_bw() +
  ggtitle("Enrichment of trans-eQTL targets in peaks") +
  xlab("Rank") + ylab("Enrichment score") +
  theme(plot.title=element_text(face="bold", size = 14)) +
  annotate(geom="text", x=10000, y=0.4, label=enr_anno, color="black", size = 6, hjust=1)




# Volcano plot of the enrichment results with interesting genes annotated
library(ggplot2)
library(ggrepel)
np$trans <- ifelse(np$gene.peak_pos %in% trans_list[['rs3087967_trans-eQTL_targets']], np$Gene.Name, "")
np$trans_yes <- ifelse(np$gene.peak_pos %in% trans_list[['rs3087967_trans-eQTL_targets']], "Yes", "No")
volcano <- ggplot(np, aes(x=signalValue, y=qValue, col=trans_yes, label=trans)) +
  geom_point() +
  theme_bw() +
  geom_text_repel(size=3.5, max.overlaps=Inf, label.size=NA, fill=NA, seed=1234, nudge_x = .15, nudge_y = 0, box.padding = 0, xlim=c(1.2,NA)) +
  scale_color_manual(values=c("black", "red")) +
  ggtitle("Peaks near trans-eQTL targets") +
  theme(plot.title=element_text(face="bold", size = 14))


library(patchwork)
volcano + enr_plot

# Looking at the bigWig file - observing the peaks themselves
# Just as a test - have a look at "SH2D6_chr2:85427323-85427674"     "SH2D6_chr2:85429212-85429634", "SH2D6_chr2:85431199-85431535"
# SH2D6 gene is chr2:85,429,448-85,437,029
# Will plot chr2:85425000-85437029
# Following https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/presentations/slides/Viz_part_2.html#9
library(Gviz)
genomeAxis <- GenomeAxisTrack(name="SH2D6")
#plotTracks(genomeAxis,from=0,to=12100)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
customFromTxDb <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                  chromosome="chr2")

bw.file <- "OCAT_NCIH211/GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_tag.ucsc.bigWig"
allChromosomeCoverage <- import.bw(bw.file,
                                   as="GRanges")
accDT <- DataTrack(allChromosomeCoverage,chomosome="chr2")
plotTracks(c(accDT,customFromTxDb,genomeAxis),
           from=85425000,to=85437029,
           chromosome="chr2",type=c("histogram"))
