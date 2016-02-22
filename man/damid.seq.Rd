\name{damid.seq}
\alias{damid.seq}
\title{DamID-Seq pipeline}
\usage{
damid.seq(raw.file.name, multi.core=F, exp.name="damid-gatc-sites",
          bin.len=100000, adapt.seq="CGCGGCCGAG", errors=1, mapping=T, qc=T,
          species="BSgenome.Celegans.UCSC.ce10", chr.names = NULL,
          restr.seq="GATC", normali=T, log10t=T, fasta.format="fastq", m.hits=1)
}
\description{
Full DamID-seq pipeline. Includes adapter removal, Mapping and enrichment finding.
}

\arguments{
  \item{raw.file.name}{ A tab-separated file with 2 or 3 columns. See raw-file.txt as example!
    \itemize{
      \item first column: the path to the raw fastq (fastq.gz) or fasta file
            (if mapping = \code{FALSE}, the file must be a grange .RData object).
      \item Second column: name of the sample.
      \item Third column: single character ‘s’ for sample and ‘c’ for control (if not set the pipeline does not perform a log2 fold change output). The number of samples ("s") must be the same as
      the number of cotrols ("c"). }}
  \item{multi.core}{ if \code{TRUE}, botwie is performed with multiple cores}
  \item{exp.name}{ a generic name for your experiment }
  \item{bin.len}{ the bin length for the analysis. }
  \item{errors}{ errors in adapter sequence allowed during adapter removal }
  \item{adapt.seq}{ the adapter sequence to be removed before mapping }
  \item{mapping}{ if \code{TRUE}, the pipeline starts with the raw fastq.gz or fastq files after sequencing. If \code{FALSE} the pipeline starts with the GRange object made by the full pipeline to skip the adapter removal and rbowtie mapping}
  \item{qc}{ if \code{TRUE}, a raw reads quality control is performed }
  \item{species}{ an object of class \code{"BSgenome"}. The name of the BSgenome object for rbowtie mapping and GATC sites extraction }
  \item{chr.names}{ a string vector of chromosome names convention from \code{\link[BSgenome]{BSgenome-utils}} or \code{NULL} to take all chromosomes from species }
  \item{restr.seq}{ the sequence of the restriction site. }
  \item{normali}{ if \code{TRUE}, the reads are normalized by dividing the reads in a bin by the total reads number of the sample. }
  \item{log10t}{ if \code{TRUE}, reads per bin are log transformed in qc plots }
  \item{fasta.format}{ input format either "fastq" (fastq, fastq.gz) or "fasta" }
  \item{m.hits}{Number of mapping positions allowed per read for bowtie mapping (1 = only unique mappable reads are mapped)}
  \item{\dots}{ not used. }
}

\details{
The tab-separated values in raw.file.name must be provided as shown in the example below to analyze a sample using controls.
Column names can be chosen freely.

  \tabular{lll}{
    FilePath/Name           \tab  sampleName  \tab  groupIndex \cr
    rawfile/sample.fastq.gz \tab  test.sample \tab  "s"        \cr
    rawfile/contro.fastq.gz \tab  test.contrl \tab  "c"
  }

The group index can be omit if no sample-control comparison is wished. (Log2 fold change between sample and control will not be plotted)


  \tabular{ll}{
    FilePath/Name           \tab  sampleName  \cr
    rawfile/sample.fastq.gz \tab  test.sample \cr
    rawfile/contro.fastq.gz \tab  test.contrl
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
damid.seq(raw.file.name = "raw-file.txt", multi.core = T,
          exp.name = "test", bin.len = 100000, errors = 1,
          mapping = T, qc = F, species = "BSgenome.Celegans.UCSC.ce10",
          chr.names = NULL, restr.seq = "GATC", normali = T, log10t = F,
          fasta.format = "fastq", m.hits = 1)



}
