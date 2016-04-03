\name{damid.seq}
\alias{damid.seq}
\title{DamID-Seq pipeline}
\usage{
damid.seq(raw.file.name, input.format="fastq", multi.core=F,
          exp.name="damid-gatc-sites",
          bin.len=100000, adapt.seq="CGCGGCCGAG", errors=1, qc=T,
          species="BSgenome.Celegans.UCSC.ce10", chr.names = NULL,
          restr.seq="GATC", normali=T, log10t=T, m.hits=1)
}
\description{
Full DamID-Seq pipeline; includes adapter removal, mapping and enrichment finding.
}

\arguments{
  \item{raw.file.name}{ A tab-separated file with 2 or 3 columns. See raw-file.txt as example!
    \itemize{
      \item first column: the path to the raw fastq (fastq.gz) or fasta file
            (if mapping = \code{FALSE}, the file must be a grange .RData object).
      \item Second column: name of the sample.
      \item Third column: single character ‘s’ for sample and ‘c’ for control (if not set the pipeline does not perform a log2 fold change output). The number of samples ("s") must be the same as
      the number of controls ("c"). }}
  \item{input.format}{ if \code{"fastq"} or \code{"fasta"}, the pipeline starts with the raw fastq/a(.gz) or fastq/a files after sequencing. If \code{"g"}, the pipeline starts with the GRange object made by the full pipeline to skip the adapter removal and rbowtie mapping. If \code{"b"} the pipeline starts with mapped .bam files}
  \item{multi.core}{ if \code{TRUE}, bowtie is performed with multiple cores}
  \item{exp.name}{ a generic name for your experiment }
  \item{bin.len}{ the bin length for the analysis }
  \item{errors}{ errors in adapter sequence allowed during adapter removal }
  \item{adapt.seq}{ the adapter sequence to be removed before mapping }
  \item{qc}{ if \code{TRUE}, a raw reads quality control is performed }
  \item{species}{ an object of class \code{"BSgenome"}, the name of the BSgenome object for rbowtie mapping and GATC sites extraction }
  \item{chr.names}{ a string vector of chromosome names convention of \code{\link[BSgenome]{BSgenome-utils}} or \code{NULL} to take all chromosomes from species }
  \item{restr.seq}{ the sequence of the restriction site }
  \item{normali}{ if \code{TRUE}, the reads are normalized by dividing the reads in a bin by the total reads numbers of the sample }
  \item{log10t}{ if \code{TRUE}, reads per bin are log transformed in qc plots }
  \item{m.hits}{Number of mapping positions allowed per read for bowtie mapping (1 = only unique mappable reads are mapped)}
  \item{\dots}{ not used. }
}

\details{
The tab-separated values in raw.file.name must be provided as shown in the example below to analyze a sample using controls.
Column names can be chosen freely.

  \tabular{lll}{
    FilePath/Name           \tab  sampleName  \tab  groupIndex \cr
    rawfile/sample.fastq.gz \tab  test.sample \tab  "s"        \cr
    rawfile/contro.fastq.gz \tab  test.control \tab  "c"
  }

The group index can be omited if no sample-control comparison is wished. (Log2 fold change between sample and control will not be plotted)


  \tabular{ll}{
    FilePath/Name           \tab  sampleName  \cr
    rawfile/sample.fastq.gz \tab  test.sample \cr
    rawfile/contro.fastq.gz \tab  test.contrl
  }

  several files are saved during the pipeline in several folders.
    \tabular{ll}{
    folder name \tab  saved files                                            \cr
    cutadapter  \tab  cut reads (.fastq/a)                                   \cr
    cutadapter  \tab  cut read length (.txt)                                 \cr
    cutadapter  \tab  cut read information (.txt)                            \cr
    bowtie      \tab  mapped reads (.bam)                                    \cr
    bowtie      \tab  bam information (.bam.bai)                             \cr
    bowtie      \tab  bam information (.txt)                                 \cr
    granges     \tab  reads per input (sample) file (.Rdata)                 \cr
    granges     \tab  reads per GATC site (.Rdata)                           \cr
    granges     \tab  reads per GATC site plus strand (.Rdata)               \cr
    granges     \tab  reads per GATC site minus strand (.Rdata)              \cr
    granges     \tab  reads per GATC fragment (.Rdata)                       \cr
    granges     \tab  GATC sites in genome                                   \cr
    qc-reports  \tab  fastqc report if \code{qc = T} (.txt)                  \cr
    qc-reports  \tab  correlation between samples: GATC and bin (.txt)       \cr
    qc-reports  \tab  number of reads lost in each step (.txt)               \cr
    results     \tab  correlation and read distribution: GATC and bin (.pdf) \cr
    results     \tab  read plots: GATC and bin (.pdf)                        \cr
    results     \tab  log2 fold change sample / control (.pdf)
  }
}

\value{
  A \code{\link[base]{list}} with 5 \code{\link[GenomicRanges]{GRanges-class}} objects.
  \tabular{lll}{
    1 \tab \code{$GATC} \tab sequencing reads per GATC site \cr
    2 \tab \code{$GATC.plus} \tab sequencing reads per GATC site on the plus strand \cr
    3 \tab \code{$GATC.minus} \tab sequencing reads per GATC site on the minus strand \cr
    4 \tab \code{$GATC.frag} \tab sequencing reads per GATC fragment \cr
    5 \tab \code{$bins} \tab sequencing reads per chosen bin
  }
}

\author{
  Dominic Ritler
}

\examples{
damid.seq(raw.file.name = "raw-file.txt", input.format = "fastq",
          multi.core = T, exp.name = "test", bin.len = 100000,
          adapt.seq="CGCGGCCGAG", errors = 1,  qc = F,
          species = "BSgenome.Celegans.UCSC.ce10",
          chr.names = NULL, restr.seq = "GATC", normali = T, log10t = F,
          m.hits = 1)



}
