##### Bradley July 2022
##### Plotting the bigwig files in R to visualise C11orf53 binding
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/ChIP-seq_env
##### Following https://bioconductor.org/packages/release/bioc/vignettes/wiggleplotr/inst/doc/wiggleplotr.html

# Set up
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/data")
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
#BiocManager::install("wiggleplotr")
library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")


# Extracting transcript ids from ensemble (testing)
library("ensembldb")
library("EnsDb.Hsapiens.v86")
#plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "SH2D6",
#                             transcript_ids = c("ENST00000340326",
#                             "ENST00000389938", "ENST00000469800",
#                             "ENST00000477170", "ENST00000481395",
#                             "ENST00000481426", "ENST00000488219",
#                             "ENST00000488657", "ENST00000651736",
#                             "ENST00000651962"))



# Downloading ALL of the neccesary data for all interseting genes
trans_ens <- c("ENSG00000137709","ENSG00000180347","ENSG00000102010","ENSG00000152292","ENSG00000183476","ENSG00000070748","ENSG00000186038","ENSG00000070985","ENSG00000160862","ENSG00000197444","ENSG00000271615","ENSG00000135407","ENSG00000167759","ENSG00000118308","ENSG00000152229","ENSG00000105851","ENSG00000197943","ENSG00000169962")
names(trans_ens) <- c("POU2F3","ITPRID1","BMX","SH2D6","SH2D7","CHAT","HTR3E","TRPM5","AZGP1","OGDHL","ACTG1P22","AVIL","KLK13","IRAG2","PSTPIP2","PIK3CG","PLCG2","TAS1R3")



  ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "jan2020.archive.ensembl.org")
  ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)
  listFilters(ensembl_dataset)
  attributes = listAttributes(ensembl_dataset)
  selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                          "external_gene_name", "strand",
                          "gene_biotype", "transcript_biotype")
  data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
  data = dplyr::rename(data,
                       transcript_id = ensembl_transcript_id,
                       gene_id = ensembl_gene_id,
                       gene_name = external_gene_name)
  temporary_file = tempfile(pattern = "file", tmpdir = "../results/bigWig_plots/", fileext = ".rds")
  saveRDS(data, temporary_file)
  transcript_metadata = readRDS(temporary_file)
  my_filter <- list(ensembl_gene_id=trans_ens)
  txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host="jan2020.archive.ensembl.org",
                           filter=my_filter)
  txdb_file = tempfile(pattern = "file", tmpdir = "../results/bigWig_plots/", fileext = ".rds")
  saveDb(txdb, txdb_file)

seqlevelsStyle(txdb) <- "UCSC"
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

# Now have the neccessary files for all genes we want to plot
# Construct the sample data file
sample_data <- as.data.frame(matrix(nrow=6, ncol=4))
colnames(sample_data) <- c("sample_id", "replicate", "scaling_factor", "bigWig")
sample_data$sample_id <- c(rep("OCAT_NCIH211", 3), rep("OCAT_NCIH526", 3))
sample_data$replicate <- rep(c("1", "2", "Input"), 2)
sample_data$scaling_factor <- rep(1, nrow(sample_data))
for(i in 1:nrow(sample_data)){
  f <- list.files(sample_data$sample_id[i])
  if(sample_data$replicate[i] %in% c("1", "2")){
    f <- f[grep(paste0("Ab", sample_data$replicate[i]), f)]
  } else {
    f <- f[grep("Input", f)]
  }
  f <- f[grep(".bigWig", f)]
  sample_data$bigWig[i] <- paste(sample_data$sample_id[i], f, sep = "/")
}
all(file.exists(sample_data$bigWig))

# Make the trackdata for subsequent plotting
track_data = dplyr::mutate(sample_data, track_id = paste(sample_id, replicate, sep = "."), colour_group = sample_id)
track_data$track_id <- gsub("OCAT_", "", track_data$track_id)

# Plot SH2D6 for example
selected_transcripts = transcript_metadata %>%
  dplyr::filter(gene_name == "HTR3E", transcript_biotype == "protein_coding")
tx_ids = selected_transcripts$transcript_id
tx_ids = tx_ids[tx_ids %in% names(cdss)]
#plotTranscripts(exons[tx_ids],
#                cdss[tx_ids],
#                transcript_metadata, rescale_introns = TRUE)

# Plot the coverage#Plot only two transcripts of the gens
library(ggplot2)
library(patchwork)
raw_plots <- vector("list", length = length(nrow(track_data)))
for(s in 1:nrow(track_data)){
  if(track_data$sample_id[s] == "OCAT_NCIH211"){
    temp <- plotCoverage(exons[tx_ids], cdss[tx_ids],
                transcript_metadata, track_data[s,],
                fill_palette = "#E9181D", return_subplots_list=T)
  } else {
    temp <- plotCoverage(exons[tx_ids], cdss[tx_ids],
                transcript_metadata, track_data[s,],
                fill_palette = "#18354B", return_subplots_list=T)
  }
#  raw_plots[[s]] <- temp[[1]] + theme_bw() + ylab("Number of reads")
  raw_plots[[s]] <- temp[[1]] + ylab("") + ylim(0,28)
}
ppi=300
png(file="../../Fig_plots/Thesis/Chapter2/HTR3E_ChIP_C53.png", height=12*ppi, width=10*ppi, res=ppi)
raw_plots[[1]] + raw_plots[[2]] + raw_plots[[3]] + raw_plots[[4]] + raw_plots[[5]] + raw_plots[[6]] + temp[[2]] +  plot_layout(nrow=7, heights = c(1,1,1,1,1,1,1.5))
dev.off()

# Merging replicates within line (using the mean of each)
track_data = dplyr::mutate(sample_data, track_id = paste(sample_id, replicate, sep = "."), colour_group = sample_id)
track_data$track_id <- gsub("OCAT_", "", track_data$track_id)
track_data$track_id <- ifelse(track_data$replicate == "Input", track_data$track_id, gsub("OCAT_", "", track_data$sample_id))
track_data$sample_id <- track_data$track_id

trans_ens <- trans_ens[-which(names(trans_ens) == "ACTG1P22")]
trans_ens <- trans_ens[-which(names(trans_ens) == "IRAG2")]
plot_list <- vector("list", length = length(trans_ens))
for(gene in seq_along(plot_list)){
  selected_transcripts = transcript_metadata %>%
    dplyr::filter(gene_name == names(trans_ens)[gene], transcript_biotype == "protein_coding")
  tx_ids = selected_transcripts$transcript_id
  tx_ids = tx_ids[tx_ids %in% names(cdss)]
  p <- plotCoverage(exons[tx_ids], cdss[tx_ids],
              transcript_metadata, track_data,
              fill_palette = c("#E9181D", "#18354B"))
  plot_list[[gene]] <- p + ylab("")
}

# Plot overlaying the replicates
track_data$track_id <- ifelse(track_data$replicate == "Input", "Input", "C11orf53_ChIP")
plot_list <- vector("list", length = length(trans_ens))
for(gene in seq_along(plot_list)){
  selected_transcripts = transcript_metadata %>%
    dplyr::filter(gene_name == names(trans_ens)[gene], transcript_biotype == "protein_coding")
  tx_ids = selected_transcripts$transcript_id
  tx_ids = tx_ids[tx_ids %in% names(cdss)]
  p <- plotCoverage(exons[tx_ids], cdss[tx_ids],
              transcript_metadata, track_data,
              fill_palette = c("#E9181D", "#18354B"), alpha=0.5, heights=c(0.75, 0.25))
  plot_list[[gene]] <- p + ylab("")
}

#### try this but removing x and y lables, using ggpubr
library(patchwork)
ppi=300
png(file="../../Fig_plots/Thesis/Chapter2/wrapped_C53_ChIP_read_plots.png", width=10*ppi, height=15*ppi, res=ppi)
wrap_plots(plot_list, ncol=2)
dev.off()




########################## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### Summarising the results so that we have only two tracks per cell line (ChIP and input)
###### Requires merging of the bed bigwig files
# From command line
# conda install -c bioconda ucsc-bigwigtowig
# Following https://www.biostars.org/p/448502/
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/data/
mkdir test_OCAT_NCIH211
cp OCAT_NCIH211/*.bigWig test_OCAT_NCIH211/
cd test_OCAT_NCIH211
# convert to wig
for f in *.bigWig
do
echo $f
bigWigToWig $f $f.wig
done
# convert to bed
# conda install -c bioconda bedops
for f in *.wig
do
echo $f
wig2bed < $f > $f.bed
done

# Now merge the bed files for the antibody samples
bedmap --echo --echo-map-score --delim '\t' <(bedops --partition GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_tag.ucsc.bigWig.wig.bed GSM5657743_NCIH211_OCA-T_Ab2_S7_R1_001.fastq.gz_tag.ucsc.bigWig.wig.bed) <(bedops --everything GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_tag.ucsc.bigWig.wig.bed GSM5657743_NCIH211_OCA-T_Ab2_S7_R1_001.fastq.gz_tag.ucsc.bigWig.wig.bed) > GSM5657743_NCIH211_OCA-T_both.bed

# Convert this back into a bedgraph
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' GSM5657743_NCIH211_OCA-T_both.bed > GSM5657743_NCIH211_OCA-T_both.bedgraph

#In case the BED file was not sorted properly:
#sort -k1,1 -k2,2n GSM5657743_NCIH211_OCA-T_both.bedgraph > GSM5657743_NCIH211_OCA-T_both_sorted.bedgraph

# Get chrom sizes
mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg38 -e 'select chrom,size from chromInfo' > myChrom.sizes
sed '1d' myChrom.sizes > myChrom2.sizes
rm myChrom.sizes
mv myChrom2.sizes myChrom.sizes

#Finally, use the UCSC bedGraphToBigWig tool:
# conda install -c bioconda ucsc-bedgraphtobigwig
bedGraphToBigWig GSM5657743_NCIH211_OCA-T_both.bedgraph myChrom.sizes GSM5657743_NCIH211_OCA-T_both.bw


###### Now in R, repeat the analysis (hopefully having combined the bigwig files)
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/ChIP-seq/data")
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
#BiocManager::install("wiggleplotr")
library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")


# Extracting transcript ids from ensemble (testing)
library("ensembldb")
library("EnsDb.Hsapiens.v86") # Check this

trans_ens <- c("ENSG00000137709","ENSG00000180347","ENSG00000102010","ENSG00000152292","ENSG00000183476","ENSG00000070748","ENSG00000186038","ENSG00000070985","ENSG00000160862","ENSG00000197444","ENSG00000271615","ENSG00000135407","ENSG00000167759","ENSG00000118308","ENSG00000152229","ENSG00000105851","ENSG00000197943","ENSG00000169962", "ENSG00000176920")
names(trans_ens) <- c("POU2F3","ITPRID1","BMX","SH2D6","SH2D7","CHAT","HTR3E","TRPM5","AZGP1","OGDHL","ACTG1P22","AVIL","KLK13","IRAG2","PSTPIP2","PIK3CG","PLCG2","TAS1R3", "FUT2")
trans_txdb <- vector("list", length=)


  ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "jan2020.archive.ensembl.org")
  ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)
  listFilters(ensembl_dataset)
  attributes = listAttributes(ensembl_dataset)
  selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                          "external_gene_name", "strand",
                          "gene_biotype", "transcript_biotype")
  data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
  data = dplyr::rename(data,
                       transcript_id = ensembl_transcript_id,
                       gene_id = ensembl_gene_id,
                       gene_name = external_gene_name)
  temporary_file = tempfile(pattern = "file", tmpdir = "../results/bigWig_plots/", fileext = ".rds")
  saveRDS(data, temporary_file)
  transcript_metadata = readRDS(temporary_file)
  my_filter <- list(ensembl_gene_id=trans_ens)
  txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host="jan2020.archive.ensembl.org",
                           filter=my_filter)
  txdb_file = tempfile(pattern = "file", tmpdir = "../results/bigWig_plots/", fileext = ".rds")
  saveDb(txdb, txdb_file)

seqlevelsStyle(txdb) <- "UCSC"
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

# Now have the neccessary files for all genes we want to plot
# Construct the sample data file
sample_data <- as.data.frame(matrix(nrow=8, ncol=4))
colnames(sample_data) <- c("sample_id", "replicate", "scaling_factor", "bigWig")
sample_data$sample_id <- c(rep("OCAT_NCIH211", 4), rep("OCAT_NCIH526", 4))
sample_data$replicate <- rep(c("1", "2", "Input", "both"), 2)
sample_data$scaling_factor <- rep(1, nrow(sample_data))
for(i in 1:nrow(sample_data)){
  if(i != 4 & i !=8){
    f <- list.files(sample_data$sample_id[i])
  } else {
    f <- list.files(paste0("test_", sample_data$sample_id[i]))
  }
  if(sample_data$replicate[i] %in% c("1", "2")){
    f <- f[grep(paste0("Ab", sample_data$replicate[i]), f)]
    f <- f[grep(".bigWig", f)]
    sample_data$bigWig[i] <- paste(sample_data$sample_id[i], f, sep = "/")
  }
  if(sample_data$replicate[i] == "Input"){
    f <- f[grep("Input", f)]
    f <- f[grep(".bigWig", f)]
    sample_data$bigWig[i] <- paste(sample_data$sample_id[i], f, sep = "/")
  }
  if(sample_data$replicate[i] == "both"){
    f <- f[grep("_both.bw", f)]
    sample_data$bigWig[i] <- paste0("test_", sample_data$sample_id[i], "/", f)
  }
}
all(file.exists(sample_data$bigWig))


# Plot this for all but the last one
track_data = dplyr::mutate(sample_data, track_id = paste(sample_id, replicate, sep = "."), colour_group = sample_id)
track_data$track_id <- gsub("OCAT_", "", track_data$track_id)

# Plot SH2D6 for example
selected_transcripts = transcript_metadata %>%
  dplyr::filter(gene_name == "HTR3E", transcript_biotype == "protein_coding")
tx_ids = selected_transcripts$transcript_id
tx_ids = tx_ids[tx_ids %in% names(cdss)]
#plotTranscripts(exons[tx_ids],
#                cdss[tx_ids],
#                transcript_metadata, rescale_introns = TRUE)

# Plot the coverage#Plot only two transcripts of the gens
track_data <- track_data[-8,]
library(ggplot2)
library(patchwork)
raw_plots <- vector("list", length = length(nrow(track_data)))
for(s in 1:nrow(track_data)){
  if(track_data$sample_id[s] == "OCAT_NCIH211"){
    temp <- plotCoverage(exons[tx_ids], cdss[tx_ids],
                transcript_metadata, track_data[s,],
                fill_palette = "#E9181D", return_subplots_list=T)
  } else {
    temp <- plotCoverage(exons[tx_ids], cdss[tx_ids],
                transcript_metadata, track_data[s,],
                fill_palette = "#18354B", return_subplots_list=T)
  }
#  raw_plots[[s]] <- temp[[1]] + theme_bw() + ylab("Number of reads")
  raw_plots[[s]] <- temp[[1]] + ylab("") + ylim(0,28)
}

raw_plots[[1]] + raw_plots[[2]] + raw_plots[[3]] + raw_plots[[4]] + raw_plots[[5]] + raw_plots[[6]] + raw_plots[[7]] + temp[[2]] +  plot_layout(nrow=8, heights = c(1,1,1,1,1,1,1.5))



####### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############### BELOW this is scrap
