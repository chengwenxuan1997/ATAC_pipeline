work.path <- "i:/genomicdata/external/FigureYa/UCSCsnapshot";setwd(work.path)
# work.path <- "/share/genomics/cwx"
article.path <- file.path(work.path, "Literatures")
data.path <- file.path(work.path, "InputData")
pkg.path <- file.path(work.path, "Packages")
code.path <- file.path(work.path, "Codes")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

invisible(lapply(ls()[grep("path", ls())], function(x){
  if (!dir.exists(get(x))) dir.create(get(x), recursive = T)
}))

# BiocManager::install("rtracklayer")
# BiocManager::install("trackViewer")
# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("BRGenomics")


library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(trackViewer)
library(BRGenomics)

# put unzip file "GSE112813_RAW.tar" under data.path (./InputData)
# It looks like ./InputData/GSE112813_RAW/**.bw

# 1: Read bw and get peak Annotation -----------------------------------------
## time consuming, up to 12GB memory (for 350M BW file)
## in order to avoid unexpected error, including exceeding memory limit, connection time out
## do not start unnecessary software or downloading session
## run it *only* at the first round analysis (just once for each project)
## all outputs for the following usage are saved

## collect sample information
Sample_info <- list.files(file.path(data.path, "GSE112813_RAW"))
Sample_info <- gsub(pattern = ".bw", replacement = "", x = Sample_info)
### some sample loss annotation, fix it mannually
spec <- grep(pattern = paste(seq(3483994, 3484002), collapse = "|"), x = Sample_info)
Sample_info[spec] <- paste0(Sample_info[spec], "_", "ATAC")
Sample_info <- gsub(pattern = "CD56_Dim_CD57Neg", replacement = "CD56dimCD57neg", x = Sample_info)
Sample_info <- gsub(pattern = "CD56_Dim_CD57Pos", replacement = "CD56dimCD57pos", x = Sample_info)

Sample_info <- strsplit(x = Sample_info, split = "_")
Sample_info <- do.call(rbind, Sample_info)
Sample_info <- as.data.frame(Sample_info)
colnames(Sample_info) <- c("Sample", "Group", "DonarID", "SeqTech")
Sample_info$Species <- ifelse(test = Sample_info$DonarID %in% c("rep1", "rep2"),
                              yes = "mouse", no = "human")
Sample_info$Group[Sample_info$Group == "Bright"] <- "CD56bright"
Sample_info$Group[Sample_info$Group == "Dim"] <- "CD56dim"
### subset the table if some samples is not needed
### The code will process all sample in this table
write.table(Sample_info, file.path(res.path, "Sample_info.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)

## annotate the peak

for (filename in list.files(file.path(data.path, "GSE112813_RAW"))){
  sample <- substr(filename, 1, 10)
  
  ### if the sample has been preprocessed before, skip it automatically
  ### So if it broke down in the former analysis, run it directly with no need to modify the script
  if (!file.exists(file.path(res.path, "PeakAnnotation", paste0(sample, "_Annotated.Grange.rds")))){
    
    ### get annotation database and annotation range
    ### considering the stability and memory usage, set annotation range: rm all comtigs
    ### change the default settings manually if needed
    ### find annotation database for each sample
    species = Sample_info$Species[Sample_info$Sample == sample]
    if (species == "mouse") {
      txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
      anno = "org.Mm.eg.db"
      seq.range <- paste0("chr", c(1:19, "X", "Y"))
      seq.length <- seqlengths(txdb)[seq.range]
    }else if(species == "human"){
      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
      anno = "org.Hs.eg.db"
      seq.range <- paste0("chr", c(1:22, "X", "Y"))
      seq.length <- seqlengths(txdb)[seq.range]
    }
    
    cat(as.character(Sys.time()), sample, species, "\n")

    ### Read file
    #### Considering the code stability and memory limits, annotate each chromosome separately
    #### Annotate the chromosome, remove all contigs
    #### It is possible that some chromosome length does not match, especially for mouse
    #### Check the reference genome version when there are many (for example, more than 3) unmatched chromosomes
    Grange <- import(file.path(data.path, "GSE112813_RAW", filename))
    Grange <- Grange[seqnames(Grange) %in% seq.range]
    seqlevels(Grange) <- seq.range
    if (sum(seqlengths(Grange) != seq.length) >= 3) Sample_info$State[Sample_info$Sample == sample] = "warning"
    if (!all(seqlengths(Grange) == seq.length)){
      cat("The sequence length of", names(which(seqlengths(Grange) != seq.length)), "not match, try to fix it\n")
      seqlengths(Grange) = seq.length
    }

    ### Annotate peak using ChIPseeker
    cat(sample, ": annotate peak\n")
    
    #### Annotate all peaks in one time (run in server)
    # Annotated.Grange <- annotatePeak(peak = Grange,
    #                                  tssRegion=c(-3000, 3000),
    #                                  annoDb = anno,
    #                                  TxDb=txdb)

    #### Annotate peaks separately and merge them (run in PC or server)
    #### annotatePeak use data.frame to store the gene symbol annotation, which costs memory greatly
    #### split the annotation into two steps: annotate the peak, turn geneID to GeneSymbol
    Annotated.Grange <- split(x = Grange, f = seqnames(Grange))
    rm("Grange"); gc()

    Annotated.Grange <- lapply(seq.range, function(x){
      cat(x, "\t")
      return(as.GRanges(annotatePeak(peak = Annotated.Grange[[x]],
                                     tssRegion = c(-3000, 3000),
                                     # annoDb = anno,
                                     TxDb = txdb, 
                                     verbose = F)))
      gc()
    })
    Annotated.Grange <- unlist(as(Annotated.Grange, "GRangesList"))

    cat("Turn GeneID to GeneSymbol\n")
    geneAnno <- suppressMessages(ChIPseeker:::getGeneAnno(annoDb = anno, 
                                                          geneID = unique(Annotated.Grange$geneId), 
                                                          type = "Entrez Gene ID"))
    Annotated.Grange$ENSEMBL <- geneAnno$ENSEMBL[match(Annotated.Grange$geneId, geneAnno$ENTREZID)]
    Annotated.Grange$SYMBOL <- geneAnno$SYMBOL[match(Annotated.Grange$geneId, geneAnno$ENTREZID)]
    Annotated.Grange$GENENAME <- geneAnno$GENENAME[match(Annotated.Grange$geneId, geneAnno$ENTREZID)]

    ## The output will be saved in ./Results/PeakAnnotation/**_Annotated.Grange.rds
    cat("save results\n")
    saveRDS(Annotated.Grange, file.path(res.path, "PeakAnnotation", paste0(sample, "_Annotated.Grange.rds")))

    cat("---------------------------------------------------\n")
    ## clear intermediate variable to process next sample
    rm("Annotated.Grange");gc()
  }else{
    cat(as.character(Sys.time()), sample, "has been processed, skip it\n")
  }
}

### check whether some sample use false reference genome
Sample_info[Sample_info$State == "warning", ]

# 2: Get gene peak score -----------------------------------------------------

## run it once for each gene to plot
## if plot.genes change, start from here

plot.genes <- c("PRF1", "SH2D1B")
Sample_info <- read.table(file.path(res.path, "Sample_info.txt"),
                          sep = "\t", header = T)
Sample_info.bk <- Sample_info
# Sample_info <- Sample_info.bk
# Sample_info <- subset(Sample_info, Sample == "GSM3084203") ## whether to subset the sample
Filtered.Granges <- lapply(Sample_info$Sample, function(sample){
  ## sort out related gene from all samples
  ## need a little time
  cat(sample, "\t")
  Annotated.Grange <- readRDS(file.path(res.path, "PeakAnnotation", paste0(sample, "_Annotated.Grange.rds")))
  Annotated.Grange <- Annotated.Grange[Annotated.Grange$SYMBOL %in% plot.genes]
  gc()
  return(Annotated.Grange)
})
names(Filtered.Granges) <- Sample_info$Sample
## save the output for reproduce the plot
saveRDS(Filtered.Granges, file.path(res.path, "Filtered.Granges.rds"))

# 3: do track plot for PRF1  -----------------------------------------------------------
## take human sample as an example
## start from the output of part 2
## if the plot need to be polished, start from here
## citing: https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/changeTracksStyles.html

## Filter samples to plot
Sample_info <- read.table(file.path(res.path, "Sample_info.txt"),
                          sep = "\t", header = T)
Filtered.Granges <- readRDS(file.path(res.path, "Filtered.Granges.rds"))
sample.names <- Sample_info$Sample[Sample_info$Group %in% c("CD56bright", "CD56dim", "CD56dimCD57neg", "CD56dimCD57pos")]
plot.Grange <- Filtered.Granges[sample.names]
Sample_info <- subset(Sample_info, Sample %in% sample.names)
plot.Grange <- lapply(sample.names, function(x){
  plot.Grange[[x]] <- subset(plot.Grange[[x]], SYMBOL == "PRF1")
  plot.Grange[[x]]$Sample <- x
  plot.Grange[[x]]$Group <- Sample_info$Group[Sample_info$Sample == x]
  plot.Grange[[x]]$Tech <- Sample_info$SeqTech[Sample_info$Sample == x]
  return(plot.Grange[[x]])
})
names(plot.Grange) <- sample.names
plot.Grange <-  as(plot.Grange, "GRangesList")

## set plot range (make a template Grange for the plot)
### plot all related regions which was annotated with ChIPseeker (wide)
temp.Grange <- GRanges(seqnames = unique(unlist(seqnames(plot.Grange))), 
                       ranges = IRanges(max(min(start(plot.Grange))):
                                          min(max(start(plot.Grange)))))
### plot genomic coordinates (narrow)
# temp.Grange <- getLocation(symbol = "PRF1", 
#                            txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#                            org = "org.Hs.eg.db")
temp.Grange

### generate gene region function annotation
### subset the annotation if needed
trs <- geneModelFromTxdb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         orgDb = org.Hs.eg.db,
                         gr = temp.Grange)
names(trs)

viewTracks(trackList = trs, gr = temp.Grange)


## 3.1: plot by sample --------------------------------------------------------

### prepare data for plotting
plot.data <- lapply(plot.Grange, function(x){
  track <- new(Class = "track", 
               dat = trackViewer:::orderedGR(x),
               type = "data",
               format = "BED")
})

### optimize style
optSty <- optimizeStyle(c(plot.data, trs))
trackList <- optSty$tracks
viewerStyle <- optSty$style

### set parameter manually
setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .01)) # specify the margin

trackList[sample.names] <- lapply(trackList[sample.names], function(track){ # set ylab 
  setTrackStyleParam(track, "ylabgp", list(cex=.5, col="#007600"))
  setTrackStyleParam(track, "ylabpos", "topright") 
})

trackList[names(trs)] <- lapply(trackList[names(trs)], function(track){
  setTrackStyleParam(track, "ylabpos", "upstream") # set transcript name location
  setTrackStyleParam(track, "ylabgp", list(cex=.8)) # set transcript name size
})

names(trackList) <- c(paste0(Sample_info$Group, "_", Sample_info$DonarID, "_", Sample_info$SeqTech), rep("PRF1", 2))

setTrackXscaleParam(trackList[[1]], "draw", TRUE) # whether draw the scale legend
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8)) 
setTrackXscaleParam(trackList[[1]], attr="position",  # adjust the location, proportion
                    value = list(x = 72394000, y = 3, label = 10000))

### plot
pdf(file.path(fig.path, "PRF1_bySample.pdf"))
viewTracks(trackList = trackList, 
           gr = temp.Grange, 
           viewerStyle = viewerStyle)
dev.off()

## 3.2: plot by Group --------------------------------------------------------


### prepare data for plotting
### average score for all sample in the same group
Sample_info <- subset(Sample_info, Sample %in% sample.names)
Sample_info$plotgroup <- paste0(Sample_info$Group, "_", Sample_info$SeqTech)
plot.data <- lapply(unique(Sample_info$plotgroup), function(x){
  sample.to.merge <- Sample_info$Sample[Sample_info$plotgroup == x]
  merged.data <- mergeGRangesData(as(plot.Grange[sample.to.merge], "list"),
                                  field = "score",
                                  ncores = 1)
  merged.data$score <- merged.data$score/length(sample.to.merge)
  return(merged.data)
})
names(plot.data) <- unique(Sample_info$plotgroup)
plot.data <- lapply(plot.data, function(x){
  track <- new(Class = "track", 
               dat = trackViewer:::orderedGR(x),
               type = "data",
               format = "BED")
})

### optimize style
optSty <- optimizeStyle(c(plot.data, trs))
trackList <- optSty$tracks
viewerStyle <- optSty$style

### set parameter manually
setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "xgp", list(cex=.8, col = "black")) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "margin", c(.05, .05, .01, .01)) # specify the margin
# setTrackViewerStyleParam(viewerStyle, "flip", TRUE) # reverse the order, if the strand is "-"

names(trackList)
setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=.8, col = "#910000")) # set label color
setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=.8, col = "#910000"))
setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=.8, col = "#007600"))
setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8, col ="#007600"))
setTrackStyleParam(trackList[[5]], "ylabgp", list(cex=.8, col = "#007600"))
setTrackStyleParam(trackList[[6]], "ylabgp", list(cex=.8, col ="#007600"))

setTrackStyleParam(trackList[[1]], "ylabpos", "topright") # set label location
setTrackStyleParam(trackList[[2]], "ylabpos", "topright")
setTrackStyleParam(trackList[[3]], "ylabpos", "topright")
setTrackStyleParam(trackList[[4]], "ylabpos", "topright")
setTrackStyleParam(trackList[[5]], "ylabpos", "topright")
setTrackStyleParam(trackList[[6]], "ylabpos", "topright")

setTrackStyleParam(trackList[[7]], "ylabpos", "upstream") # set transcript name location
setTrackStyleParam(trackList[[8]], "ylabpos", "downstream")
setTrackStyleParam(trackList[[7]], "ylabgp", list(cex=.8)) # set transcript name size
setTrackStyleParam(trackList[[8]], "ylabgp", list(cex=.8))

setTrackStyleParam(trackList[[1]], "color", c("#910000")) # set track color
setTrackStyleParam(trackList[[2]], "color", c("#910000"))
setTrackStyleParam(trackList[[3]], "color", c("#007600"))
setTrackStyleParam(trackList[[4]], "color", c("#007600"))
setTrackStyleParam(trackList[[5]], "color", c("#007600"))
setTrackStyleParam(trackList[[6]], "color", c("#007600"))

setTrackStyleParam(trackList[[7]], "height", .1) # set track height
setTrackStyleParam(trackList[[8]], "height", .1)

names(trackList) <- c("CD56bright_K27ac", "CD56dim_K27ac", "CD56bright_ATAC", "CD56dim_ATAC",
                      "CD56dimCD57neg_ATAC", "CD56dimCD57pos_ATAC", rep("PRF1", 2))

setTrackXscaleParam(trackList[[1]], "draw", TRUE) # 
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8))
setTrackXscaleParam(trackList[[1]], attr="position", 
                    value = list(x = 72393000, y = 50, label = 10000))

### change plot range manually
temp.Grange <- GRanges(seqnames = unique(unlist(seqnames(plot.Grange))), 
                       ranges = IRanges(72350000:72380000))

### plot
pdf(file.path(fig.path, "PRF1_byGroup.pdf"))
viewTracks(trackList = trackList, 
           gr = temp.Grange, 
           viewerStyle = viewerStyle)
addGuideLine(guideLine = c(72318647,	72377415),
             col = "red",
             vp = viewTracks(trackList = trackList, 
                             gr = temp.Grange, 
                             viewerStyle = viewerStyle))
dev.off()



# 4: Another Example: SH2D1B ---------------------------------------------------------
## It is the same with part 3 completely, ignore it if it is no use

## Filter samples to plot
Sample_info <- read.table(file.path(res.path, "Sample_info.txt"),
                          sep = "\t", header = T)
Filtered.Granges <- readRDS(file.path(res.path, "Filtered.Granges.rds"))
sample.names <- Sample_info$Sample[Sample_info$Group %in% c("CD56bright", "CD56dim", "CD56dimCD57neg", "CD56dimCD57pos")]
plot.Grange <- Filtered.Granges[sample.names]
Sample_info <- subset(Sample_info, Sample %in% sample.names)
plot.Grange <- lapply(sample.names, function(x){
  plot.Grange[[x]] <- subset(plot.Grange[[x]], SYMBOL == "SH2D1B")
  plot.Grange[[x]]$Sample <- x
  plot.Grange[[x]]$Group <- Sample_info$Group[Sample_info$Sample == x]
  plot.Grange[[x]]$Tech <- Sample_info$SeqTech[Sample_info$Sample == x]
  return(plot.Grange[[x]])
})
names(plot.Grange) <- sample.names
plot.Grange <-  as(plot.Grange, "GRangesList")

## set plot range
### plot all related regions which was annotated with ChIPseeker (wide)
temp.Grange <- GRanges(seqnames = unique(unlist(seqnames(plot.Grange))), 
                       ranges = IRanges(max(min(start(plot.Grange))):
                                          min(max(start(plot.Grange)))))
### plot genomic coordinates (narrow)
# temp.Grange <- getLocation(symbol = "PRF1", 
#                            txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#                            org = "org.Hs.eg.db")
temp.Grange

### generate gene region function annotation
### subset the annotation if needed
trs <- geneModelFromTxdb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                         orgDb = org.Hs.eg.db,
                         gr = temp.Grange)
names(trs)

viewTracks(trackList = trs, gr = temp.Grange)


## 4.1: plot by sample --------------------------------------------------------

### prepare data for plotting

plot.data <- lapply(plot.Grange, function(x){
  track <- new(Class = "track", 
               dat = trackViewer:::orderedGR(x),
               type = "data",
               format = "BED")
})

### optimize style
optSty <- optimizeStyle(c(plot.data, trs))
trackList <- optSty$tracks
viewerStyle <- optSty$style

### set parameter manually
setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .01)) # specify the margin

trackList[sample.names] <- lapply(trackList[sample.names], function(track){
  setTrackStyleParam(track, "ylabgp", list(cex=.5, col="#007600"))
  setTrackStyleParam(track, "ylabpos", "topright") 
})

trackList[names(trs)] <- lapply(trackList[names(trs)], function(track){
  setTrackStyleParam(track, "ylabpos", "upstream") # set transcript name location
  setTrackStyleParam(track, "ylabgp", list(cex=.8)) # set transcript name size
})

names(trackList) <- c(paste0(Sample_info$Group, "_", Sample_info$DonarID, "_", Sample_info$SeqTech), rep("SH2D1B", 2))

setTrackXscaleParam(trackList[[1]], "draw", TRUE) # whether draw the scale legend
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8)) 
temp.Grange
setTrackXscaleParam(trackList[[1]], attr="position",  # adjust the location, proportion
                    value = list(x = 162410000, y = 3, label = 10000))

### plot
pdf(file.path(fig.path, "SH2D1B_bySample.pdf"))
viewTracks(trackList = trackList, 
           gr = temp.Grange, 
           viewerStyle = viewerStyle)
dev.off()

## 4.2: plot by Group --------------------------------------------------------


### prepare data for plotting
Sample_info <- subset(Sample_info, Sample %in% sample.names)
Sample_info$plotgroup <- paste0(Sample_info$Group, "_", Sample_info$SeqTech)
plot.data <- lapply(unique(Sample_info$plotgroup), function(x){
  sample.to.merge <- Sample_info$Sample[Sample_info$plotgroup == x]
  merged.data <- mergeGRangesData(as(plot.Grange[sample.to.merge], "list"),
                                  field = "score",
                                  ncores = 1)
  merged.data$score <- merged.data$score/length(sample.to.merge)
  return(merged.data)
})
names(plot.data) <- unique(Sample_info$plotgroup)
plot.data <- lapply(plot.data, function(x){
  track <- new(Class = "track", 
               dat = trackViewer:::orderedGR(x),
               type = "data",
               format = "BED")
})

### optimize style
optSty <- optimizeStyle(c(plot.data, trs))
trackList <- optSty$tracks
viewerStyle <- optSty$style

### set parameter manually
setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "xgp", list(cex=.8, col = "black")) # whether to plot x-axis
setTrackViewerStyleParam(viewerStyle, "margin", c(.05, .05, .01, .01)) # specify the margin
# setTrackViewerStyleParam(viewerStyle, "flip", TRUE) # reverse the order, if the strand is -

names(trackList)
setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=.8, col = "#910000")) # set label color
setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=.8, col = "#910000"))
setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=.8, col = "#007600"))
setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8, col ="#007600"))
setTrackStyleParam(trackList[[5]], "ylabgp", list(cex=.8, col = "#007600"))
setTrackStyleParam(trackList[[6]], "ylabgp", list(cex=.8, col ="#007600"))

setTrackStyleParam(trackList[[1]], "ylabpos", "topright") # set label location
setTrackStyleParam(trackList[[2]], "ylabpos", "topright")
setTrackStyleParam(trackList[[3]], "ylabpos", "topright")
setTrackStyleParam(trackList[[4]], "ylabpos", "topright")
setTrackStyleParam(trackList[[5]], "ylabpos", "topright")
setTrackStyleParam(trackList[[6]], "ylabpos", "topright")

setTrackStyleParam(trackList[[7]], "ylabpos", "upstream") # set transcript name location
setTrackStyleParam(trackList[[8]], "ylabpos", "downstream")
setTrackStyleParam(trackList[[7]], "ylabgp", list(cex=.8)) # set transcript name size
setTrackStyleParam(trackList[[8]], "ylabgp", list(cex=.8))

setTrackStyleParam(trackList[[1]], "color", c("#910000")) # set track color
setTrackStyleParam(trackList[[2]], "color", c("#910000"))
setTrackStyleParam(trackList[[3]], "color", c("#007600"))
setTrackStyleParam(trackList[[4]], "color", c("#007600"))
setTrackStyleParam(trackList[[5]], "color", c("#007600"))
setTrackStyleParam(trackList[[6]], "color", c("#007600"))

setTrackStyleParam(trackList[[7]], "height", .1) # set track height
setTrackStyleParam(trackList[[8]], "height", .1)

names(trackList) <- c("CD56bright_K27ac", "CD56dim_K27ac", "CD56bright_ATAC", "CD56dim_ATAC",
                      "CD56dimCD57neg_ATAC", "CD56dimCD57pos_ATAC", rep("SH2D1B", 2))

setTrackXscaleParam(trackList[[1]], "draw", TRUE) # 
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8))
setTrackXscaleParam(trackList[[1]], attr="position", 
                    value = list(x = 162410000, y = 50, label = 10000))

### change plot range manually

temp.Grange <- GRanges(seqnames = unique(unlist(seqnames(plot.Grange))), 
                       ranges = IRanges(162365056:162389000))

### plot
pdf(file.path(fig.path, "SH2D1B_byGroup.pdf"))
viewTracks(trackList = trackList, 
           gr = temp.Grange, 
           viewerStyle = viewerStyle)
addGuideLine(guideLine = c(162367000),
             col = "red",
             vp = viewTracks(trackList = trackList, 
                             gr = temp.Grange, 
                             viewerStyle = viewerStyle))
dev.off()

