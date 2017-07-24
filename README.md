# RDamIDSeq
All in R DamID sequencing pipeline.  
Includes adapter removal, Mapping and enrichment finding for DamID-Seq reads.

For more information see also: [Tools for DNA adenine methyltransferase identification analysis of nuclear organization during C. elegans development](http://onlinelibrary.wiley.com/doi/10.1002/dvg.22925/abstract)

install
-------

1) Install R (only 64 bit version!)  
[R Project](https://www.r-project.org/ "The R Project for Statistical Computing")


2) Download DamIDseq_0.1.3.tar.gz

3) install additional packages needed for the pipeline in R   

`install.packages(c('corrplot', 'aqfig'))`  
`source("https://bioconductor.org/biocLite.R")`  
`biocLite()`  
`biocLite(c("QuasR", "ShortRead", "GenomicRanges", "S4Vectors", "IRanges",`   
`"BiocInstaller", "Biobase", "Biostrings", "BSgenome", "GenomicFeatures", `    
`"GenomicAlignments", "BiocParallel", "GenomeInfoDb", "rtracklayer", `    
`"GenomicFiles", "Rbowtie"))`  

4) install DamIDSeq in R   
(exchange `PATH` with the path to the installation file on your disk)  

`install.packages("PATH/DamIDseq_0.1.3.tar.gz", repos = NULL, type = "source")`  




