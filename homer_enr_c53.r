# Bradley May 2022
# Motif enrichment using homer
# Follow annot_explore_all.r
# source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/ChIP-seq_env

# Load packages and data
library('Biostrings')
library('BSgenome.Hsapiens.UCSC.hg38')
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/results")
trans_set="c11_cor_genes"
peaks <- read.csv(paste0(trans_set, "/peaks_near_eQTL_targets_c53.csv"), row.names=1, header=T)
peaks <- peaks[,1]

# Remove C11orf53, COLCA1 and COLCA2
peaks <- peaks[-grep("C11orf53|COLCA1|COLCA2", peaks)]

# Make list and fill with sequence
seq_list <- vector("list", length = length(peaks))
chr_search <- unlist(strsplit(peaks, "\\:"))[c(T,F)]
chr_search <- unlist(strsplit(chr_search, "\\_"))[c(F,T)]
start_search <- unlist(strsplit(peaks, "\\:"))[c(F,T)]
start_search <- as.numeric(unlist(strsplit(start_search, "\\-"))[c(T,F)])
end_search <- unlist(strsplit(peaks, "\\:"))[c(F,T)]
end_search <- as.numeric(unlist(strsplit(end_search, "\\-"))[c(F,T)])
names(seq_list) <- peaks
for(peak in seq_along(peaks)){
  print(peaks[peak])
  seq_list[[peak]] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, chr_search[peak], start_search[peak], end_search[peak]))
}
peak_seq <- matrix(nrow=length(peaks), ncol=2)
colnames(peak_seq) <- c("peak_name", "sequence")
peak_seq <- as.data.frame(peak_seq)
peak_seq$peak_name <- peaks
seqs <- unlist(seq_list)
names(seqs) <- NULL
peak_seq$sequence <- unlist(seqs)
write.table(peak_seq, paste0(trans_set, "/motif_enr/peak_sequences_c53.txt"), sep = "\t", row.names=F, quote=F)

### Save this as a bed file for the motfi discovery with HOMER
bed <- data.frame(chromsome=chr_search, start=start_search, end=end_search, id=peaks, not=rep("", length(peaks)), Strand=rep(".", length(peaks)))
write.table(bed, paste0(trans_set, "/motif_enr/peak_c53.bed"), sep = "\t", row.names=F, quote=F)

## then from the command line (assuming using the c11_cor_genes eQTLs - OTHERWISE USE THE RNASEQ ONES)
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/results/c11_cor_genes/motif_enr
# One time only, install genome into homer perl /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/ChIP-seq_env/share/homer/.//configureHomer.pl -install hg38
findMotifsGenome.pl peak_c53.bed hg38 ./ -size 200

# Finding which sequences have the motif we find here. Make a matrix of these and plot a heatmap
hits <- list(POU2F3="ATGCAAAT|TTTGCAT", EGR1="TGTGGGTGC|GCACCCACA", RUNX="GCAGTGC|GCACTGC")
hit_genes <- vector("list", length=length(hits))
names(hit_genes) <- names(hits)
tot <- nrow(peak_seq)
for(hit in seq_along(hits)){
  temp <- peak_seq[grep(hits[[hit]], peak_seq$sequence),]
  hit_genes[[hit]] <- unique(unlist(strsplit(temp$peak_name, "\\_"))[c(T,F)])
}
targs <- unique(unlist(strsplit(peaks,"\\_"))[c(T,F)])
res_mat <- matrix(nrow=length(hit_genes), ncol=length(targs))
dimnames(res_mat) <- list(names(hits), targs)
for(r in 1:nrow(res_mat)){
  for(c in 1:ncol(res_mat)){
    if(targs[c] %in% hit_genes[[r]]){
      res_mat[r,c] <- 1
    }
  }
}
# Plot
library(ggplot2)
library(reshape2)
library(dplyr)
res_mat <- as.data.frame(res_mat)
res_mat$motif <- rownames(res_mat)
df <- melt(res_mat)
colnames(df) <- c("x", "y", "value")
df$value <- ifelse(is.na(df$value), 0, df$value)
#df$y <- factor(df$y, levels=levels(df$y)[order(levels(df$y), decreasing=T)])
df$x <- factor(df$x, levels = c("RUNX", "EGR1", "POU2F3"))
ppi=300
png(file="../../Fig_plots/Thesis/Chapter2/C11orf53_ChIP_targets.png", width=9*ppi, height=3*ppi,res=ppi)
ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile(color="black") +
  coord_flip() +
  scale_fill_gradient2(low = "#FFFFFF", high = "#CD0000") +
  theme_minimal() +
  theme(panel.background = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x=element_text(size=12))
dev.off()




# make a matrix of this to plot



####### Now do the same for POU2F3
peaks <- read.csv("peaks_near_eQTL_targets_pou.csv", row.names=1, header=T)
peaks <- peaks[,1]

# Make list and fill with sequence
seq_list <- vector("list", length = length(peaks))
chr_search <- unlist(strsplit(peaks, "\\:"))[c(T,F)]
chr_search <- unlist(strsplit(chr_search, "\\_"))[c(F,T)]
start_search <- unlist(strsplit(peaks, "\\:"))[c(F,T)]
start_search <- as.numeric(unlist(strsplit(start_search, "\\-"))[c(T,F)])
end_search <- unlist(strsplit(peaks, "\\:"))[c(F,T)]
end_search <- as.numeric(unlist(strsplit(end_search, "\\-"))[c(F,T)])
names(seq_list) <- peaks
for(peak in seq_along(peaks)){
  print(peaks[peak])
  seq_list[[peak]] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, chr_search[peak], start_search[peak], end_search[peak]))
}
peak_seq <- matrix(nrow=length(peaks), ncol=2)
colnames(peak_seq) <- c("peak_name", "sequence")
peak_seq <- as.data.frame(peak_seq)
peak_seq$peak_name <- peaks
seqs <- unlist(seq_list)
names(seqs) <- NULL
peak_seq$sequence <- unlist(seqs)
write.table(peak_seq, "motif_enr_pou/peak_sequences_pou.txt", sep = "\t", row.names=F, quote=F)

### Save this as a bed file for the motfi discovery with HOMER
bed <- data.frame(chromsome=chr_search, start=start_search, end=end_search, id=peaks, not=rep("", length(peaks)), Strand=rep(".", length(peaks)))
write.table(bed, "motif_enr_pou/peak_pou.bed", sep = "\t", row.names=F, quote=F)

## then from the command line
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/results/motif_enr_pou
findMotifsGenome.pl peak_pou.bed hg38 ./ -size 200
