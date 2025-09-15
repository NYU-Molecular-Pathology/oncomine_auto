#
library(readxl)
library(tidyr)
library(ggplot2)
library(ini)
args = commandArgs(trailingOnly=TRUE)
runID <- args[1]
# read config

# read the config file
conf <- read.ini("/Users/yangy15/Documents/oncomine_automation/ion_config.conf")

# Access values like a list
output_path = conf$MYELOSEQ$DEST_PATH
work_dir = getwd()
dropout_dir = conf$MYELOSEQ$DROPOUT_DIR
#
PATH_Samplesheet <- paste0(output_path,"/worksheet.dropoffs/",
                           runID, ".xlsm")
# clean up samplesheet
samplesheet <- read_excel(PATH_Samplesheet, sheet = 1, skip = 5)
samplesheet <- as.data.frame(na.omit(samplesheet[,c("Bar code", "Accession #", "DNA #", "Chip #")]))
targets <- read.table(paste0(work_dir, "/QC_data/Oncomine_Myeloid.20170817.designed.DNA.bed"), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

chrX <- targets[targets$V1 == "chrX",]
samplesheet$`Bar code` <- ifelse(samplesheet$`Bar code` < 10,
                                 paste0("IonXpress_00", samplesheet$`Bar code`), paste0("IonXpress_0", samplesheet$`Bar code`))
samplesheet$sample_name <- paste0(samplesheet$`Accession #`,"-", samplesheet$`DNA #`)
chips <- unique(samplesheet$`Chip #`)
# get number of chips
get_one_plot <- function(samplesheet, runID, chip_num) {
  coverage_PATH <- paste0(dropout_dir,"/",runID)
  coverage_file <- list.files(coverage_PATH)[grepl(list.files(coverage_PATH),pattern = paste0(runID,"-", chip_num))]
  coverage_n <- read.table(paste0(coverage_PATH, "/",coverage_file), header = 1)
  chip <- samplesheet[samplesheet$`Chip #` == chip_num,]
  names(coverage_n)[-c(1,2)] <- chip$sample_name[match(names(coverage_n)[-c(1,2)], chip[,"Bar code"])]
  coverage_n$chrX <- ifelse(coverage_n$Target %in% chrX$V4, "chrX", "Not chrX")
  coverage_n <- coverage_n[!is.na(colnames(coverage_n))] # if we skip samples from samplesheet, this help prevent empty columns
  coverage_n_data <- coverage_n[,-c(1,2)]
  long_n <- coverage_n_data %>% gather(sample_name, coverage, -chrX)
  p <- ggplot(long_n, aes(x=chrX, y = coverage, fill=chrX)) + ylim(0,6000) +
    facet_wrap(~sample_name) +geom_dotplot(binaxis='y')+geom_boxplot(fill="white")+labs(title=paste0(runID, "-", chip_num))
  return(p)
}
output_PATH <- paste0("./")
pdf(paste0(output_PATH, "/", runID, "_chrX_coverage.pdf"),height = 10, width = 15)
lapply(chips, function(x) get_one_plot(samplesheet, runID, x))
dev.off()
