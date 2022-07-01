# ATAC_pipeline
a preliminary pipeline to process the .bw file

**Data**  
Take [GSE112813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112813) for example.  

**Purpose**  
1: Annotate each peak as the gene symbols using R package [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).  
2: Draw gene-related peaks distribution based on their chromosome position using [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html).  
