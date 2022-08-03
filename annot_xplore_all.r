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

# 15/07/22 edit. Decide which set of trans-eQTL targets to perform the enrichment on
trans_set <- "HT12_eQTLs"
if(trans_set == "HT12_eQTLs"){
  trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")
}
if(trans_set == "RNAseq_eQTLs"){
  trans <- c("POU2F3", "ITPRID1", "BMX", "SH2D6", "SH2D7", "CHAT", "HTR3E", "TRPM5", "AZGP1", "OGDHL", "ACTG1P22", "AVIL", "KLK13", "IRAG2", "PSTPIP2", "PIK3CG", "PLCG2", "TAS1R3")
}

require(rtracklayer)
library(ggplot2)
library(ggrepel)
np_list <- vector("list", length = length(conds))
plots <- vector("list", length = length(np_list))
titles <- vector("list", length = 4)
enr_plots <- vector("list", length=4)
for(s in seq_along(enr_plots)){
  enr_plots[[s]] <- vector("list", length=2)
}
for(cond in seq_along(conds)){
  # Set up
  files=list.files(conds[cond])
  np_files <- files[grep(".narrowPeak$", files)]
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
  np_list[[cond]] <- vector("list", length = length(np_files))
  plots[[cond]] <- vector("list", length=length(np_files))
  titles[[cond]] <- vector("list", length=length(np_files))
  # Load raw files
  for(rep in seq_along(np_files)){
    narrowPeak_file <- paste(conds[cond], np_files[rep], sep = "/")
    gr_narrowPeak <- import(narrowPeak_file, format="BED", extraCols=extraCols_narrowPeak, genome="hg38")
    # Now read in the annotated version of this
    anno <- read.delim(paste(narrowPeak_file, "annot", sep = "."), header=T)
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
    np_list[[cond]][[rep]] <- np
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
    if(conds[cond] == "OCAT_NCIH211"){
      enr_plot <- plotEnrichment(trans_list[[1]], gene_list) +
        theme_bw() +
        ggtitle("Enrichment of trans-eQTL targets in peaks") +
        xlab("Rank") + ylab("Enrichment score of 11q23.1 trans-eQTL targets") +
        theme(plot.title=element_text(face="bold", size = 14)) +
        annotate(geom="text", x=(0.95*length(gene_list)), y=0.35, label=enr_anno, color="black", size = 6, hjust=1) + ggtitle(paste(gsub("OCAT_", "", conds[cond]), rep, sep = "."))
    } else {
      enr_plot <- plotEnrichment(trans_list[[1]], gene_list) +
        theme_bw() +
        ggtitle("Enrichment of trans-eQTL targets in peaks") +
        xlab("Rank") + ylab("Enrichment score of 11q23.1 trans-eQTL targets") +
        theme(plot.title=element_text(face="bold", size = 14)) +
        annotate(geom="text", x=(0.95*length(gene_list)), y=0.25, label=enr_anno, color="black", size = 4, hjust=1) + ggtitle(paste(gsub("OCAT_", "", conds[cond]), rep, sep = "."))
    }
    enr_plots[[cond]][[rep]] <- enr_plot
      # Volcano plot of the enrichment results with interesting genes annotated
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
    title=np_files[rep]
    title <- paste(unlist(strsplit(title, "\\_"))[c(F,T,T,T,T,F,F)], collapse = "_")
    titles[[cond]][[rep]] <- title
    plots[[cond]][[rep]] <- volcano + enr_plot + plot_annotation(title = title,  theme = theme(plot.title = element_text(size = 18, face="bold")))
    ppi=300
    png(file=paste("../results/", trans_set, "/", title, "_volcano_fgsea.png", sep = ""), width=10*ppi, height=5*ppi, res=ppi)
    plot(plots[[cond]][[rep]])
    dev.off()
    print("Plotted")
  }
}

# Have a peak at the enr_plots
library(patchwork)
enr_saves <- enr_plots
enr_saves[[1]][[2]] <- enr_saves[[1]][[2]] + theme(axis.title.y = element_blank())
enr_saves[[2]][[1]] <- enr_saves[[2]][[1]] + theme(axis.title.y = element_blank())
enr_saves[[2]][[2]] <- enr_saves[[2]][[2]] + theme(axis.title.y = element_blank())
ppi=300
png(file="../../Fig_plots/11q23_transeQTLs_C53_ChIPseq_fig.png", width=15*ppi, height=4*ppi,res=ppi)
enr_saves[[1]][[1]] + enr_saves[[1]][[2]] + enr_saves[[2]][[1]] + enr_saves[[2]][[2]] + plot_layout(ncol = 4)
dev.off()

# Summarizing which eQTL targets are picked up in which cells
# Find the list of unique peaks in each experiment
cis_trans <- c("C11orf53", "COLCA1", "COLCA2", trans)
trans_peaks <- vector("list", length = length(conds))
in_trans <- vector("list", length = length(conds))
in_trans_broad <- vector("list", length = length(conds))
for(cond in seq_along(conds)){
  trans_peaks[[cond]] <- vector("list", length = length(plots[[cond]]))
  in_trans[[cond]] <- vector("list", length = length(plots[[cond]]))
  in_trans_broad[[cond]] <- vector("list", length = length(plots[[cond]]))
  for(rep in seq_along(trans_peaks[[cond]])){
    trans_peaks[[cond]][[rep]] <- np_list[[cond]][[rep]][np_list[[cond]][[rep]]$Gene.Name %in% cis_trans,]
    in_trans[[cond]][[rep]] <- levels(factor(trans_peaks[[cond]][[rep]]$gene.peak_pos))
    in_trans_broad[[cond]][[rep]] <- levels(factor(trans_peaks[[cond]][[rep]]$Gene.Name))
  }
}
all_in_trans <- unlist(in_trans)
all_in_trans <- all_in_trans[!duplicated(all_in_trans)]
all_in_trans <- all_in_trans[order(all_in_trans)]
all_in_trans_broad <- unlist(in_trans_broad)
all_in_trans_broad <- all_in_trans_broad[!duplicated(all_in_trans_broad)]
all_in_trans_broad <- all_in_trans_broad[order(all_in_trans_broad)]

# plot the broad enrichment per gene
res_mat <- matrix(nrow=length(all_in_trans_broad), ncol = 8)
rownames(res_mat) <- all_in_trans_broad
colnames(res_mat) <- unlist(titles)
colnames(res_mat) <- gsub(".narrowPeak", "", colnames(res_mat))
colnames(res_mat) <- gsub(".narrowPeak", "", colnames(res_mat))
colnames(res_mat) <- gsub("_bt2_snstv.sorted_trimmed_rd.bam", "", colnames(res_mat))
colnames(res_mat) <- gsub("snstv.sorted_trimmed_rd.bam_peaks", "", colnames(res_mat))
for(cond in seq_along(conds)){
  for(rep in seq_along(np_list[[cond]])){
    col_mat <- 2*cond-(2-rep)
    for(gene in seq_along(all_in_trans_broad)){
      if(all_in_trans_broad[gene] %in% np_list[[cond]][[rep]]$Gene.Name){
        res_mat[gene,col_mat] <- max(np_list[[cond]][[rep]][np_list[[cond]][[rep]]$Gene.Name == all_in_trans_broad[gene],]$signalValue)
      } else {
        res_mat[gene,col_mat] <- 0
      }
    }
  }
}


# Plot heatmap
library(ComplexHeatmap)
hm <- Heatmap(res_mat, name = "SignalValue", column_title = "Expt",
          row_title = "Gene", column_title_gp = gpar(fontsize=14, fontface="bold"),
          cluster_rows=T, cluster_columns=T, column_names_rot = 45,
          row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize=12))
png(file=paste("../results/", trans_set, "/Heatmap_eQTL_targets_per_expr.png", sep = ""), width=8*ppi, height=8*ppi, res=ppi)
draw(hm, padding = unit(c(2, 15, 2, 15), "mm"))
dev.off()

# Save results
write.csv(res_mat, paste("../results/", trans_set, "/QTL_targets_per_expr.csv", sep = ""))
write.csv(all_in_trans, paste("../results/", trans_set, "/peaks_near_eQTL_targets.csv", sep = ""))


# Save a list of peaks found in just the C11orf53 bound sites
in_trans_c53 <- in_trans[c(1,2)]
all_in_trans_53 <- unlist(in_trans_c53)
all_in_trans_53 <- all_in_trans_53[!duplicated(all_in_trans_53)]
all_in_trans_53 <- all_in_trans_53[order(all_in_trans_53)]
write.csv(all_in_trans_53, paste("../results/", trans_set, "/peaks_near_eQTL_targets_c53.csv", sep = ""))

# Save a list of peaks found in just the POU2F3 bound sites
in_trans_pou <- in_trans[c(3,4)]
all_in_tans_pou <- unlist(in_trans_pou)
all_in_tans_pou <- all_in_tans_pou[!duplicated(all_in_tans_pou)]
all_in_tans_pou <- all_in_tans_pou[order(all_in_tans_pou)]
write.csv(all_in_tans_pou, paste("../results/", trans_set, "/peaks_near_eQTL_targets_pou.csv", sep = ""))


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
