pdf(file=NULL)
#Import libraries
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(ggplot2)

#Set working directory
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Oncosolid_QC_plot.R work_dir")
}

wd <- args[1]

######## process sensitive control data ######################
#oncominedata<-read.csv("sc_filtered_variants.tsv",sep = "\t")
# now handles two controls
QC_file.list = list.files(path = wd, "*_filtered_variants.tsv", full.names = T)
run_id <- read.csv(QC_file.list[1],sep = "\t")[1,1]
process_SC_control <- function(x) {
  oncominedata<-read.csv(x,sep = "\t")
  if (grepl( "o2", oncominedata[1,2], fixed = TRUE)) {
    chip = 2
  } else { chip = 1 }
  oncominedata$locus<-paste0(oncominedata$Genes,":", oncominedata$Coding)
  
  DNAqcdata<-subset(oncominedata,  locus  %in% c (
    "NRAS:c.182A>G",
    "IDH1:c.394C>T",
    "PIK3CA:c.1633G>A",
    "KIT:c.2447A>T",
    "EGFR:c.2236_2250delGAATTAAGAGAAGCA",
    "EGFR:c.2573T>G",
    "EGFR:c.2369C>T",
    "BRAF:c.1799T>A",
    "GNAQ:c.626A>C",
    "KRAS:c.35G>A",
    "ERBB2:c.2313_2324dup",
    "GNA11:c.626A>T" )  )
  
  run_id <- DNAqcdata[1,1]
  DNAqcdata$AF <- as.numeric(as.character(DNAqcdata$X..Frequency))
  print("getting DNAQC")
  DNAQC <- ggplot(DNAqcdata, aes(x= locus, y=AF,color= locus)) +
    geom_point(size=6)+
    #ggtitle(paste("Oncomine Solid",run_id,"chip#", chip, "DNA QC: 12 Sites Expected",sep = " ")) +
    ggtitle(paste("Oncomine Solid",run_id,"DNA QC", paste("(chip#",chip,")",sep=""),":12 Sites Expected",sep = " ")) +
    xlab("Genomic Location") + ylab("Allele Frequency") +
    geom_hline(yintercept=3, linetype="dashed",  color = "red", linewidth=1)+
    geom_hline(yintercept=15, linetype="dashed",  color = "blue", linewidth=1)+
    theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
    scale_color_manual(values=c("goldenrod3", "brown2", "deepskyblue1", "plum1", "turquoise2", "orchid", "greenyellow","seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2" ))+
    guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))
  
  fusions_df <- subset(oncominedata, Type == "FUSION" | Type == "RNAExonVariant")
  fusions_df$Fusions <- paste(fusions_df$Genes,"(",fusions_df$Exon,")",sep="")
  RNAqcdata<-subset(fusions_df, Fusions  %in% c (
    "EML4|ALK(13|20)",
    "KIF5B|RET(24|11)",
    "NCOA4|RET(7|12)",
    "CD74|ROS1(6|34)",
    "SLC34A2|ROS1(4|34)",
    "TPM3|NTRK1(7|10)",
    "FGFR3|BAIAP2L1(17|2)",
    "PAX8|PPARG(9|2)",
    "FGFR3|TACC3(17|11)",
    "ETV6|NTRK3(5|15)",
    "LMNA|NTRK1(2|11)",
    "SLC45A3|BRAF(1|8)",
    "TMPRSS2|ERG(1|2)",
    "MET|MET(13|15)"  )  )
  
  RNAqcdata$Read.Counts[is.na(RNAqcdata$Read.Counts)] <- 1
  print("making RNAQC")
  RNAQC <- ggplot(RNAqcdata, aes(x=Fusions, y=log2(as.numeric(as.character(Read.Counts))),color= Fusions)) +
    geom_point(size=6)+ xlab("Fusions(Exome)") + ylab("log2 Read Counts") +
    #ggtitle(paste("Oncomine Solid SC-QC", run_id, "RNA QC: 14 Events Expected",sep = " ")) +
    ggtitle(paste("Oncomine Solid SC-QC", run_id, "RNA QC",paste("(chip#",chip,")",sep=""),": 14 Events Expected",sep = " ")) +
    geom_hline(yintercept=log2(100), linetype="dashed",  color = "red", linewidth=1)+
    theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text =  element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
    scale_color_manual(values=c("seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2", "brown", "deepskyblue3", "plum3", "turquoise", "orchid2", "greenyellow","salmon2" ))+
    guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))
  return(list(DNAQC,RNAQC))
}
#QCplot = list(DNAQC, RNAQC)
QCplot = lapply(QC_file.list, process_SC_control)
pdf(paste("Oncomine-Solid-SC-QC-",run_id,".pdf",sep=""))
QCplot
dev.off()

####### process SNP QC ##############
#read unfilted TSV files
files_to_read <- list.files(path=paste0(wd,"/Variants"), pattern = "_Non-Filtered.*-oncomine\\.tsv$", recursive = TRUE)

#Read all tsv files in directory
all_files <- lapply(files_to_read,function(x) {
  read.csv(file = paste(paste0(wd,"/Variants"),x,sep="/"), quote="", sep = '\t', header = TRUE, skip=4)
})


#Combine file content list and file name list
all_lists <- mapply(c, files_to_read,all_files, SIMPLIFY = FALSE)

#Unlist all lists and add file name as new column
all_result <- rbindlist(all_lists, fill = T)

#Rename new column
names(all_result)[1] <- "Sample"

#Shorten Sample name by removing all characters after the first underscore
all_result$Sample<-gsub("_.*","",all_result$Sample)

#List of locus
snps <- read.csv(paste0(wd,"/data/ONCsolid.SNP_locus.csv"))
snps <- snps$Locus

#add the 2 CNV number replacing the VAF
combined_final <- distinct(all_result)
combined_final <- as.data.frame(sapply(all_result, function(x) gsub("\"", "", x)))

data0 <- subset(combined_final, !(X.call. == "NEG") & !(X.call. == "NOCALL")& (X.rowtype. == "snp"))
data0<-data0[,c("Sample"  , "X.CHROM." ,"X.POS." , "X.INFO.A.AF.", "X.INFO...CI.")]
#data0<-combined_final[,c("Sample"  , "X.CHROM." ,"X.POS." , "X.INFO.A.AF.", "X.INFO...CI.")]
data0$X..Frequency<-paste0(data0$X.INFO.A.AF., data0$X.INFO...CI.)
data0$Locus <- paste(data0$X.CHROM.,data0$X.POS.,sep = ":")
data0$X..Frequency<-gsub("NA","",data0$X..Frequency)
data0 <- data0[,c("Sample","Locus","X..Frequency")]
#data0 = data0[!duplicated(data0$Locus),]
res.targeloci<-subset(distinct(data0), Locus  %in% snps)
data.e<-spread(res.targeloci, Sample, X..Frequency)
data.e[is.na(data.e)] <- 0

write.csv(data.e, file=paste("QC-SNPs-final-",run_id,".csv",sep=""), row.names=FALSE)
write.csv(combined_final, file=paste("combined_final_",run_id,".csv",sep=""), row.names=FALSE)

imbalancefusion<-subset(combined_final,  (grepl(glob2rx("5p3pAssays") , X.rowtype.) )
                        & (as.numeric(X.INFO.1.5P_3P_ASSAYS.) > 0)
                        & !(X.call. == "NEG") & !(X.call. == "NOCALL"))
write.csv(imbalancefusion, file=paste("imbalancefusion_final_",run_id,".csv",sep=""), row.names=FALSE)

# expcontrol QC
expcontrol<-subset(combined_final,  (grepl(glob2rx("ExprControl") , X.rowtype.) ))
write.csv(expcontrol, file=paste("expcontrol_",run_id,".csv",sep=""), row.names=FALSE)

#Read.Counts
expQC<- ggplot(expcontrol, aes(y=log2(as.numeric(as.character(X.INFO...READ_COUNT.))), x=Sample, color=Sample)) +
  geom_point(size=3)+ ylab("log2 Read Counts") +
  theme(axis.text.x = element_text(size = 6, angle = 65))+
  facet_wrap( ~X.INFO.1.GENE_NAME., ncol =2)+
  
  geom_hline(yintercept=log2(1800), linetype="dashed",  color = "red", linewidth=0.5)+
  geom_hline(yintercept=log2(700), linetype="dashed",  color = "blue", linewidth=0.5)+
  ggtitle("Case Expression Control ONC.Solid QC") +
  guides(col = guide_legend ( ncol=1 )) +
  theme(legend.text=element_text(size=7))

pdf(paste("Case_Expression_Control_ONC_Solid_QC_",run_id,".pdf",sep=""))
expQC
dev.off()

## process filtered full TSV file #######
tsv_files_to_read <- list.files(path=paste0(wd,"/Variants"), pattern = "-full.tsv$", recursive = TRUE)

#Read all tsv files in directory
tsv_files <- lapply(tsv_files_to_read,function(x) {
  read.csv(file = paste(paste0(wd,"/Variants"),x,sep="/"), quote="", sep = '\t', header = TRUE, skip=2)
})


#Combine file content list and file name list
content_lists <- mapply(c, tsv_files_to_read,tsv_files, SIMPLIFY = FALSE)

#Unlist all lists and add file name as new column
filtered_result <- rbindlist(content_lists, fill = T)

#Rename new column
names(filtered_result)[1] <- "Sample"

#Shorten Sample name by removing all characters after the first underscore
filtered_result$Sample<-gsub("_.*","",filtered_result$Sample)

#List of locus
snps <- read.csv(paste0(wd,"/data/ONCsolid.SNP_locus.csv"))
snps <- snps$Locus

#add the 2 CNV number replacing the VAF
combined_tsv <- distinct(filtered_result)
combined_tsv <- as.data.frame(sapply(filtered_result, function(x) gsub("\"", "", x)))
biallelicdata<-subset(combined_tsv, (grepl("?|p", protein, fixed = TRUE))
                      & (filter == "PASS") & (type %in% c("MNV","SNV","INDEL,MNV","INDEL","INDEL,SNV")))
write.csv(biallelicdata, file=paste("biallelicdata-final-",run_id,".csv",sep=""), row.names=FALSE)
#added snp check
set2 <- as.data.frame(pivot_wider(
  res.targeloci,
  names_from = "Sample",
  values_from = "X..Frequency",
  values_fn = list(X..Frequency= list)) %>%
  unchop(everything()))
# data.e<-spread(res.targeloci, Sample, X..Frequency)
set2[is.na(set2)] <- 0
#set2_processed <- as.data.frame((set2))
rownames(set2) <- set2$Locus
set2 <- set2[-1]
set2[] <- lapply(set2, function(x) as.numeric(as.character(x)))
reformat_data <- function(data1,data2=NULL) { #current run vs past run
  if(is.null(data2)) {corr_data <- cor(data1, use = "pairwise",method = "pearson")} else {
    # make two data frames compatible
    snps <- intersect(rownames(data1),rownames(data2))
    # if(dim(data2)[2] > dim(data1)[2]){corr_data <- cor(data2[snps,], data1[snps,], use = "pairwise",method = "pearson")} else {
    #   corr_data <- cor(data1[snps,], data2[snps,], use = "pairwise",method = "pearson")
    corr_data <- cor(data1[snps,], data2[snps,], use = "pairwise",method = "pearson")
    # }
  }
  # align the methods with ggcorr###
  m <- as.data.frame(corr_data*lower.tri(corr_data, diag = T))
  m$.ggally_ggcorr_row_names = rownames(m)
  m = reshape::melt(m, id.vars = ".ggally_ggcorr_row_names")
  names(m) = c("x", "y", "coefficient")
  return(m)
}
get_plot <- function(data, runID1, runID2) { # current run itself
  ggplot(data, aes(x, y, fill = coefficient)) + 
    theme(panel.background = element_blank()) + ggtitle(paste0(runID1," comparison with ", runID2)) +
    geom_tile(color = "white") + guides(fill = guide_legend(title = "Pearson\nCorrelation")) + 
    scale_y_discrete(position = "right") +
    xlab(runID1) +
    ylab(runID2) +
    scale_fill_gradient2(
      low = "darkblue",
      high = "tomato",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab"
    ) +
    theme(
      axis.title.x = element_text(size=5),
      axis.title.y = element_text(size=5),
      axis.text.x = element_text(angle = 45, vjust = 1, size = 5, hjust = 1),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(size=16, face="bold", hjust = 0.5),
      legend.position = "bottom", legend.key=element_rect(colour="black")
    ) +
    coord_fixed() 
}

pdf(paste0("snps_Comparison_",run_id,".pdf"))
get_plot(data=reformat_data(data1=set2),runID1 = run_id,runID2 = run_id)
dev.off()

# remove empty 
if (file.exists("./Rplots.pdf")) 
  file.remove("./Rplots.pdf")
