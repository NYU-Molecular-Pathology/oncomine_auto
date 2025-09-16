#!/usr/bin/env Rscript

## Script name: MGON-gene_coverage_plots.R
## Purpose of script: Read Oncosolid amplicon dropoff xls and save gene coverage PDF per sample chip
## Author: Jonathan Serrano
## Date Created: 2022-11-15

gb <- globalenv(); assign("gb", gb)
args <- commandArgs(TRUE)

specificRun <- args[1]
config_path <- args[2]
# read config

# read the config file
library("ini")
conf <- read.ini(config_path)

# Access values like a list
onco_InDir <- conf$SOLID$DEST_PATH
outFolder <- getwd()
xlsDropOff <- "amplicon.dropout.dropoff"

LoadLibrary <- function(pkgName) {
    suppressPackageStartupMessages(library(
        pkgName, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE, character.only = TRUE))
}

LoadPackageLibs <- function(){
    corePkgs <- c("readxl", "ggplot2", "tidyr", "stringr") # "qpdf"
    unlist(lapply(corePkgs, LoadLibrary))
}

#CreatePdfOutput <- function(newRunFolder){
 #   if(!dir.exists(newRunFolder)){
  #      message("Creating new Directory:\n", newRunFolder)
   #     dir.create(newRunFolder)
    #}
#}

XlsmChipParse <- function(shName, totalChips, colIdx = c(1:3, 10)) {
    samplesheet <- suppressMessages(readxl::read_excel(shName, sheet = 1, skip = 5))
    samplesheet <- as.data.frame(na.omit(samplesheet[, colIdx]))
    colnames(samplesheet) <- c("Barcode", "Accession #", "DNA #", "Chip #")
    samplesheet$Barcode <- ifelse(
        samplesheet$Barcode < 10,
        paste0("IonXpress_00", samplesheet$Barcode),
        paste0("IonXpress_0", samplesheet$Barcode)
    )
    samplesheet$sample_name <- paste0(samplesheet$`Accession #`, "-", samplesheet$`DNA #`)
    if(totalChips==2){
    chip1 <- samplesheet[samplesheet$`Chip #` == 1, ]
    chip2 <- samplesheet[samplesheet$`Chip #` == 2, ]
    return(list(chip1, chip2))
    }else{
        chip1 <- samplesheet[samplesheet$`Chip #` == 1, ]
        return(list(chip1))
    }
}

GetCoverage <- function(xlsFiles, chips){
    cov <- lapply(X=1:length(xlsFiles), function(idx){
        chipDat<- chips[[idx]]
	dim(chipDat)
        coverage <- read.table(xlsFiles[[idx]], header = 1)[-2]
        head(coverage)
        chipCols <- chipDat[, "Barcode"] %in% names(coverage)
        print(chipCols)
	print(head(chipDat))
	chipDat_f <- chipDat[chipCols,]
	print(dim(chipDat_f))
	#chipBars <- chipDat[, "Barcode"][chipCols]
        names(coverage)[-1] <-
		#chipDat$sample_name[match(names(coverage)[-1], chipBars)]
            chipDat_f$sample_name[match(names(coverage)[-1], chipDat_f$Barcode)]
        return(coverage)
    })
    return(cov)
}

GetFacetPlot <- function(coverage, lim=NULL) {
    long_1 <- coverage %>% tidyr::gather(sample_name, coverage, -Gene)
    p1 <-
        ggplot2::ggplot(long_1, aes(x = Gene, y = coverage, fill = Gene)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = 'none') +
        ggplot2::facet_wrap(~ sample_name, scales = 'free_x', ncol = 3) + ggplot2::geom_boxplot(outlier.fill = "red")
    if(!is.null(lim)){ p1 <- p1+coord_cartesian(ylim=c(0, 30000))}
    p1 <- p1 + ggplot2::theme(
        text = ggplot2::element_text(size = 20), plot.title = ggplot2::element_text(size = 24),
        strip.text.x = ggplot2::element_text(size = 30), plot.margin = ggplot2::margin(1, 1, 1, 1, "in"),
        panel.spacing = ggplot2::unit(.05, "lines"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = ggplot2::element_rect(color = "black", linewidth = 0.5,),
        panel.grid.minor.y = ggplot2::element_line(color = "#8ccde3", size = 0.25, linetype = 1),
        panel.grid.major.y = ggplot2::element_line(color = "#8ccde3", size = 0.5, linetype = 2)
        )
    return(p1)
}

GenerateBoxplots <- function(chp, coverage, runName){
    outFile <- paste(runName,"chip", chp, "coverage.csv", sep="_")
    #utils::write.csv(coverage[[chp]], file=outFile, quote = F, row.names=F)
    p1 <-  GetFacetPlot(coverage[[chp]])
    p2 <- GetFacetPlot(coverage[[chp]],lim = 30000)
    pdfOutput <- file.path(gb$outFolder)
    #CreatePdfOutput(pdfOutput)
    fiName <- paste0(runName, "_gene_coverage_chip", chp, ".pdf")
    message("Saving File... \n", fiName)
    ggplot2::ggsave(file.path(pdfOutput, fiName), plot = gridExtra::marrangeGrob(list(p1, p2), nrow = 1, ncol = 1), width = 11, height = 8.5, units = "in", device = "pdf", dpi = 350, scale = 3.5)
}

Merge_PDF_files <- function(runName){
    fiName <-  paste0(runName, "_gene_coverage.pdf")
    fiOutPath <- file.path(getwd(), fiName)
    pdfOutput <- file.path(gb$outFolder, runName)
    pdfFiles <- dir(pdfOutput, pattern=".pdf", all.files=T, full.names=T)
    message("Combining Files:\n", paste(capture.output(pdfFiles), collapse = "\n"))
    #qpdf::pdf_combine(input = pdfFiles, output = fiOutPath)
    message("Saving merged pdfs to directory:\n", fiOutPath)
    cmd1 <- paste("rm -rf", pdfOutput)
    #system(cmd1)
}

LoopChipOutput <- function(shFile, totalChips, xlsFiles, runName){
    chips <- XlsmChipParse(shFile, totalChips)
    coverage <- GetCoverage(xlsFiles, chips)
    for(chp in 1:totalChips){
        GenerateBoxplots(chp, coverage, runName)
    }
}

GetXlsFiles <- function(totalChips, xlsFiLi, xlsMatch) {
    message("Total Chips on run: ", totalChips)
    xlsFi1 <- xlsFiLi[xlsMatch][1]
    if (totalChips == 1) {
        xlsFiles <- list(xlsFi1)
    } else{
        xlsFi2 <- xlsFiLi[xlsMatch][2]
        xlsFiles <- list(xlsFi1, xlsFi2)
    }
    return(xlsFiles)
}


CheckRunFile <- function(runsToDo, specificRun) {
    if (!is.null(specificRun)) {
        runsToDo <- runsToDo[grepl(specificRun, runsToDo)]
        print(runsToDo)
        if (length(runsToDo) == 0) {
            msg1=paste0('Worksheet "', specificRun, '", was not found in the onco_InDir directory:\n')
            message(msg1, file.path(onco_InDir, "dropoff"))
            stopifnot(length(runsToDo) != 0)
        }
        return(runsToDo)
    }
}

GetRunsList <- function(dirToCheck, specificRun, exculsions= "test|xlsx|~\\$"){
    runsToDo <- dir(dirToCheck, specificRun, T, T, T)
    excludeFi <- grepl(exculsions, runsToDo, ignore.case = T)
    runsToDo <- CheckRunFile(runsToDo[!excludeFi], specificRun)
    return(runsToDo)
}

GetAllFilesMatch <- function(xlsMatch){
    realFile <- table(xlsMatch)[["FALSE"]]
    xLen <- length(xlsMatch)
    return(realFile >= xLen - 2 & realFile < xLen)
}

ParseChips <- function(xlsMatch, xlsFiLi, shFile){
    runName <- stringr::str_split_fixed(basename(shFile), ".xlsm", 2)[[1]]
    if (GetAllFilesMatch(xlsMatch)) {
        message("Running ", runName, "...")
        totalChips <- length(xlsFiLi[xlsMatch])
        xlsFiles <- GetXlsFiles(totalChips, xlsFiLi, xlsMatch)
        LoopChipOutput(shFile, totalChips, xlsFiles, runName)
    } else{message("Skipping ", runName, "...")}
}

ReadFileValues <- function(xlsFiLi, runsToDo, outFolder) {
    for (shFile in runsToDo) {
        runName <- stringr::str_split_fixed(basename(shFile), ".xlsm", 2)[[1]]
     #   newRunFolder <- file.path(outFolder, runName)
     #   CreatePdfOutput(newRunFolder)
        message("Reading '", runName, "' samplesheet:\n", shFile)
        xlsMatch <- grepl(pattern = runName, xlsFiLi, ignore.case = T)
        ParseChips(xlsMatch, xlsFiLi, shFile)
        #Merge_PDF_files(runName)
    }
}

GetFileList <- function(specificRun = NULL, onco_InDir, xlsDropOff, outFolder) {
    xlsFiLi <- dir(path = file.path(onco_InDir, xlsDropOff), pattern = ".xls", full.names = T)
    runsToDo <- GetRunsList(file.path(onco_InDir, "dropoff"), specificRun)
    ReadFileValues(xlsFiLi, runsToDo, outFolder)
}

LoadPackageLibs()

GetFileList(specificRun, onco_InDir, xlsDropOff, outFolder)
if (file.exists("./Rplots.pdf")) 
  file.remove("./Rplots.pdf")
