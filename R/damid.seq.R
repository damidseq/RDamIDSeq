# Author:    Dominic Ritler (DR)
# Contact:   dominic.ritler@students.unibe.ch
# Version:   0.1.3
#
# Header:    DamID pipeline all in one. (Adapter removal, bowtie, analysis)


################################################################################
##                               functions                                    ##
################################################################################

fuse.chr <- function(query.gr) {
  # Fuse all chromosome listed in a grange object by adding the length of the
  # previous chromosomes to the GATC position of the current chromosome
  # resulting in a single string of GATC coordinates.
  #
  # Args:
  #   query.gr:  (grange) grange object
  #
  # Return:
  #   list with all GATC coordinates over all chromosomes added together

  # get length
  chr.leng <- seqlengths(query.gr)[seqnames(query.gr)@values]

  # make list with length of all listed chromosomes
  startplot <- rep(0, length(seqnames(query.gr)@values))

  # if more than 1 chromosome is chosen
  if (length(chr.leng) > 1) {

    # first position is 0 then sum up the sum of the lengths of the previous
    # chromosomes to the next chromosome.
    for (i in 2:length(seqnames(query.gr)@values)) {
      startplot[i] <- sum(chr.leng[1:i - 1])
    }

    # get the number of "hits" GATC sites for each chromosome
    coordinates <- Rle(startplot, runLength(seqnames(query.gr)))
    # add start value to the right position
    coordinates <- as.vector(coordinates) + as.vector(start(query.gr))

  } else {

    # make coordinates from 0 to length of the only chosen chromosome
    coordinates <- as.vector(start(query.gr))
  }

  return(coordinates)
}


read.sam.files.to.grange <- function(file.name) {
  # Read a bam/sam file and make a GRange with all reads from the sam/bam file.
  #
  # Args:
  #   file.name:  (String) the name of the bam file to import
  #
  # Return:
  #   new genomicRange with the reads from the sam/bam file

  # import sam/bam file
  mapps.galing <- readGAlignments(file.name)
  # make GRange
  mapps.grang  <- as(mapps.galing, "GRanges")

  return(mapps.grang)
}


find.restriction.sites.overlaps3 <- function(bam.grange,
                                             res.site.grages,
                                             strand="*", genom, chrom.nam,
                                             restsite="GATC") {
  #  Find all reads starting or ending with a GATC motif and return a grange
  #  listing all reads starting (plus strand) or/and ending (minus) with
  #  a GATC site. This step is to map the DNA seq bam/sam files to the list of
  #  all GATC sites.
  #
  # Args:
  #   bam.grange: (GRange) grange object with al the reads from sequencing
  #   res.site.grages (vector of grange) GATC grange (all,plus,minus,frag)
  #   strand: (chr) a string indicating the strand either "*", "+", "-"
  #                 or "all" to return all in a list.
  #   genom: (string) the name of the genome
  #   chrom.nam: (vector of strings) the chromosome names
  #   restsite: the restriction site string
  #
  # Return:
  #   vector with the number of matches of query in subject grange
  #   if strand == "all" a list with 3 vectors are returned (all, plus, minus)

  # get chromosome name
  gatc.plus  <- res.site.grages[[2]]
  gatc.minus <- res.site.grages[[3]]

  rest.sites.plus  <- gatc.plus[which(seqnames(gatc.plus) %in% chrom.nam)]
  rest.sites.minus <- gatc.minus[which(seqnames(gatc.minus) %in% chrom.nam)]

  # find overlaps for plus and minus strand
  olaps.start <- findOverlaps(bam.grange, rest.sites.plus, type = "start")
  olaps.end   <- findOverlaps(bam.grange, rest.sites.minus, type = "end")


  # make frequency (+ = plus strand only, - = minus strand only,
  # * = both strands, all = return list with plus, minus and all strand)
  if (strand == "+") {
    freq <- table(subjectHits(olaps.start))
  } else if (strand == "-") {
    freq <- table(subjectHits(olaps.end))
  } else if (strand == "*") {
    freq <- table(c(subjectHits(olaps.start), subjectHits(olaps.end)))
  } else if (strand == "all") {
    freq.all <- table(c(subjectHits(olaps.start), subjectHits(olaps.end)))
    freq.plu <- table(subjectHits(olaps.start))
    freq.min <- table(subjectHits(olaps.end))
  } else {# if not valid statement use * and give warning message.
    message(paste("The strand information", strand,
                  "is not valid: use '*', '+' or '-'. now '*' was used!"))
    freq <- table(c(subjectHits(olaps.start), subjectHits(olaps.end)))
  }

  # make vector with number of reads per GATC site
  if (strand == "all") {

    # initiate all GATC sites with 0 (if the site has no hits it will not be
    # overwritten in the next step therefore sett al to 0 in every iteration)
    scor.all <- rep(0,length(rest.sites.plus))#minus and plus GATC is palindrome
    scor.all[as.vector(as.integer(rownames(freq.all)))] <- freq.all

    scor.plu <- rep(0,length(rest.sites.plus))
    scor.plu[as.vector(as.integer(rownames(freq.plu)))] <- freq.plu

    scor.min <- rep(0,length(rest.sites.minus))
    scor.min[as.vector(as.integer(rownames(freq.min)))] <- freq.min

    scores <- list(scor.all, scor.plu, scor.min)

  } else {
    scores <- rep(0,length(rest.sites.plus))
    scores[as.vector(as.integer(rownames(freq)))] <- freq
  }

  return(scores)
}


add.metadata.to.grange <- function(metadata, grange, m.name) {
  # Add a metadata line to the end of the metadata in a grange file
  #
  # Args:
  #   metadata: (list) number of counts per GATC sites in the GRange
  #   grange: (GRange) the grange file to add the metadata
  #   m.name: (string) the name of the metadata column
  #
  # Return:
  #   GRange including the new metadata column at the end of the metadata.

  # get metadata from GRange
  meth.data <- elementMetadata(grange)

  # make data-frame with name
  new.meth <- data.frame(metadata)
  names(new.meth) <- m.name

  # bind new metadata column to the existing metadata
  final.meth  <- cbind(meth.data, new.meth)

  # add metadata back to the grange
  elementMetadata(grange) <- final.meth

  return(grange)
}

make.empty.genomicRange.bin3 <- function(animal, chrom.names,
                                         bin.size=100000, output.mesage=T) {
  # Make a empty genomicRange object with specific chosen bin size
  # (each bin one genomic range entry)
  #
  # Args:
  #   animal: (BSgenome) BSgenome object of a species
  #           species name as animal as species is a object "BioGenerics"
  #   chrom.names: (vector of strings) the names of the chromosomes
  #   bin.size: (int) length of the bins (length of each genomicRnage entry)
  #   output.mesage: (bool) if progress should be reported
  #
  # Return:
  #   new genomicRange object with specific genomic ranges
  #

  # set bin.grange zo NULL
  bin.grange <- NULL # old version is not working as exist search in global env. ****************************

  # go over all chromosomes chosen in chrom.names
  for (chr in chrom.names) {
    # get length of chromosome
    chr.length <- seqlengths(animal)[chr]

    # get last full range (bin.size=100, chrom.length=310, last full range=300)
    full.ranges1 <- floor(as.integer(chr.length)  / bin.size)
    full.ranges <- full.ranges1 # (minus 1 bin size)
    comp.length <- full.ranges * bin.size # the length to get only full ranges

    # test if the chromosome length is bigger or equal than the bin size.
    if (chr.length >= bin.size) {
      # indicator for clean up
      # 1 means all tmp.gr, tmp.gr.ful, tmp.gr.last have to be cleaned
      clean.all <- 1

      # make grange wit only full ranges
      tmp.gr.ful <- GRanges(seqnames = chr,
                            ranges = IRanges(start = seq(from = 1,
                                                         to = comp.length,
                                                         bin.size),
                                             width = bin.size,
                                             names = paste(chr,
                                                           seq(1, comp.length,
                                                               bin.size),
                                                           sep = ":")),
                            strand = "*", seqinfo = seqinfo(animal))


      # add the last range (not full)
      # but before test the extreme unlikely event that the chromosome has
      # exact the length of the last full range
      if (comp.length + bin.size != as.integer(chr.length)) {
        tmp.gr.last <- GRanges(seqnames = chr,
                               ranges = IRanges(start = comp.length + 1,
                                                width = chr.length - comp.length,
                                                names = paste(chr,
                                                              comp.length + 1,
                                                              sep = ":")),
                               strand = "*", seqinfo = seqinfo(animal))

        # add (concatenate) the last not full bin to the grange
        tmp.gr <- c(tmp.gr.ful, tmp.gr.last)
      } else {
        # indicator for chean up
        # 2 means all tmp.gr, tmp.gr.ful have to be cleaned
        clean.all <- 2

        # no last not full bin has to be concatenated
        tmp.gr <- tmp.gr.ful #if the range is completely and no must be appended
      }

    } else {# if the chromosome is smaller than the bin size.

      # indicator for clean up
      # 3 means tmp.gr have to be cleaned
      clean.all <- 3

      tmp.gr <- GRanges(seqnames = chr,
                        ranges = IRanges(start =  1,
                                         width = chr.length,
                                         names = paste(chr,
                                                       1,
                                                       sep = ":")),
                        strand = "*", seqinfo = seqinfo(animal))

    }

    # test if grange exist if not (first for loop iteration) make one
    # otherwise concatenate the reads from the chromosome
    if (is.null(bin.grange)) {
      bin.grange <- tmp.gr
    } else {
      bin.grange <- append(bin.grange, tmp.gr)
    }

    # clean-up
    if (clean.all == 1) {
      rm(tmp.gr, tmp.gr.last, tmp.gr.ful)
    } else if (clean.all == 2) {
      rm(tmp.gr, tmp.gr.ful)
    } else if (clean.all == 3) {
      rm(tmp.gr)
    }

    # writhe message if wanted
    if (output.mesage) {
      message(paste(format(Sys.time(), "%H:%M:%S:"),
                    "Chromosome ... ",chr," done!"))
    }

  }
  return(bin.grange)
}



################################################################################
##                            packages management                             ##
################################################################################
# This functions are obsolete when using the DamIDseq R package.
# Only BSgenomes when changing the species may need downloading packages
# from bioconductor

install.biolite.packages <- function(bc.install.names) {
  # install bioclite packages
  #
  # Args:
  #   install.names: (string) names of the packages to be installed from biocon
  #
  # Return:
  #   none (install packages if not installed)

  source("http://bioconductor.org/biocLite.R")
  biocLite()

  biocLite(bc.install.names)

}


install.r.packages <- function(r.install.names) {
  # install r packages
  #
  # Args:
  #   install.names: (string) names of the R packages to be installed from CRAN
  #
  # Return:
  #   none (install packages if not installed)

  install.packages(r.install.names)

}


test.installed.bc.packages <- function(bc.package.names) {
  # test if bioconductor package are installed and if not install it.
  #
  # Args:
  #   package.names: (vector of string) packages names for the pipeline
  #
  # Return:
  #   none (install packages if not installed)

  # find not installed packages
  installed <- bc.package.names %in% rownames(installed.packages())

  # install package if not installed but first test if list is empty
  if (length(bc.package.names[!installed]) != 0) {
    install.biolite.packages(bc.package.names[!installed])
  } else {
    print("all bioconductor packages installed")
  }
}


test.installed.r.packages <- function(r.package.names) {
  # test if r package are installed and if not install it.
  #
  # Args:
  #   package.names: (vector of string) packages to be needed for the analysis
  #
  # Return:
  #   none (install packages if not installed)

  # find not installed packages
  installed <- r.package.names %in% rownames(installed.packages())

  # install package if not installed but first test if list is empty
  if (length(r.package.names[!installed]) != 0) {
    install.r.packages(r.package.names[!installed])
  } else {
    print("all R packages installed")
  }
}


load.pakage.main <- function(chr.packages, bio.packages) {
  # all in one to load and install bioconductor and cran (r) packages
  #
  # Args:
  #   chr.packages: (vector of string) all package names for R cran
  #   bio.packages: (vector of string) all package names for bioconductor
  #
  # Return:
  #   none (install packages if not installed and load packages)

  if (bio.packages != "") {
    # test if all bioconductor packages are installed and if not install it
    test.installed.bc.packages(bio.packages)
    # load all bioconductor library’s
    lapply(bio.packages, require, character.only = TRUE)
  }

  if (chr.packages != "") {
    # test if all r packages are installed and install it if not
    test.installed.r.packages(chr.packages)
    # load all r library’s
    lapply(chr.packages, require, character.only = TRUE)
  }
}


################################################################################
##                         binning reads in GRanges                           ##
################################################################################

find.all.overlaps.and.sum.metadata3 <- function(query.gr, subject.gr,
                                                met.data=T, messag=T,
                                                sum.bin=T) {
  # Find genome overlaps in two GRange objects and return a list of hits per
  # grange of the query grange (use when binning granges).
  # either sum (sum.bin=T) metadata or mean metadata (sum.bin=F)
  #
  # Args:
  #   query.gr:   (GRange) query grange file
  #   subject.gr: (GRange) subject grange file
  #   meth.data:  (bool) F = matrix will be returned
  #                      T = query.gr will be returned with metadata
  #   messag: (bool) if a message should be made or not
  #   sum.bin: (bool) if T values in bin are summed if F the mean is taken
  #
  # Return:
  #   If met.data=F datata.frame with number of overlaps in metadata
  #   file and ranges from query (number of reads per bin)
  #   If met.data=T Bin Grange with metadata (GATC reads per bin)

  if (messag) {
    bin.siz <- width(subject.gr[1])
    message(paste(format(Sys.time(), "%H:%M:%S:"),
                  "find overlaps and sum up metadata entries for bin size:",
                  bin.siz, "bp"))
  }

  # calculate number of metadata rows
  met.rows <- dim(elementMetadata(query.gr))[2]

  # make empty matrix for new Grange (bin) metadata, each column is a sample
  bin.meta <- matrix(NA, length(subject.gr), met.rows) # sample matrix

  # Find overlaps (ov) between GATC Grange and bin Grange
  # This approach counts the GATC reads twice if the GATC site overlaps
  # bins (for example bin 1 ends with "GA", bin 2 starts with "TC").
  # type="within" cannot be used as it excludes the overlapping reads)
  overlaps <- findOverlaps(query.gr, subject.gr, ignore.strand = TRUE)

  # Get number of GATC sites for each separate bin.
  # Value is the position in the (bin grange),
  # length is the number of GATC site per bin
  overlaps.table <- Rle(subjectHits(overlaps))

  # get number of GATC for each bin separate and make list with start and end
  # of each bin (start, end = GATC grange index position)
  gatc.site.index <- PartitioningByWidth(width(overlaps.table))

  # decide if bin gets summed up or the mean of the bin is taken
  if (sum.bin) {
    sum.or.mean <- viewSums
  } else {
    sum.or.mean <- viewMeans
  }

  # go over all metadata and sum up reads per bin
  for (i in 1:met.rows) {
    # make list with 0 for each bin (0 reds per bin)
    reads.per.bin <- rep(0, length(subject.gr))

    # add the sum of GATC reads for each bin if bin has GATC reads
    reads.per.bin[runValue(overlaps.table)] <- sum.or.mean(Views(Rle(
      elementMetadata(query.gr)[queryHits(overlaps), i]), gatc.site.index))

    # add the sum of GATC reads per bin to the sample matrix
    bin.meta[,i] <- reads.per.bin
  }

  # assign the same column name from the GATC grange metadata to the bin Grange
  colnames(bin.meta) <- names(elementMetadata(query.gr))

  # return grange with metadata of matrix
  if (met.data) {
    elementMetadata(subject.gr) <- bin.meta
    return(subject.gr)
  } else {
    return(bin.meta)
  }
}


################################################################################
#     function to process already mapped reads starting with read Granges      #
################################################################################

make.full.gatc.grange.from.read.grange2 <- function(file.nam,
                                                    save.nam="gran.all",
                                                    messag=T, genom, chr.na,
                                                    restsite="GATC",
                                                    res.site.grages) {
  # Combine different function when starting the pipeline with read Granges
  # instead of a sam/bam file (sam reads as grange) to get reads per GATC grange
  #
  # Args:
  #   file.nam: (data.frame) data frame with 2 or 3 columns:
  #                          path: the path to the Grange file
  #                          name: the name of the sample
  #                          group: optional! 's' for sample or 'c' for control
  #   save.nam: (string) the name and path of the file to save the
  #                      resulting GATC grange
  #   messag: (bool) if a message should be returned
  #   genom: (string) the name of the BSgenome object
  #   chr.na (vector of string) the name of the chromosomes
  #   restsite: (string) the restriction site sequence
  #   res.site.grages (list of granges) all GATC granges (all,plus,minus,fragm)
  #
  # Return:
  #   grange with metadata of all reads

  # make empty GATC Granges
  gatc.all.m     <- res.site.grages[[1]] # for all GATC sites
  gatc.plus.m    <- res.site.grages[[2]] # for plus sites
  gatc.minus.m   <- res.site.grages[[3]] # for minus sites
  gatc.frag.m    <- res.site.grages[[4]]
  #rm(gatc.all)

  # print message if asked
  if (messag) {
    message(paste(format(Sys.time(), "%H:%M:%S:"),
                  "Load grange object and start mapping to restriction sites"))
  }

  # go over al samples and add read number to GATC GRange
  for (i in 1:length(file.nam[,1])) {

    # load grange object
    read.grange <- load.restriction.sites(toString(file.nam[i,1]))

    # map to GATC sites for
    mapped.to.gatc <- find.restriction.sites.overlaps3(read.grange,
                                                       res.site.grages,
                                                       strand = "all", genom,
                                                       chrom.nam = chr.na,
                                                       restsite = restsite)
    # add metadata to grange for all, plus and minus strand
    gatc.all.m <- add.metadata.to.grange(mapped.to.gatc[[1]], gatc.all.m,
                                         toString(file.nam[i,2]))
    gatc.plus.m <-  add.metadata.to.grange(mapped.to.gatc[[2]], gatc.plus.m,
                                           toString(file.nam[i,2]))
    gatc.minus.m <- add.metadata.to.grange(mapped.to.gatc[[3]], gatc.minus.m,
                                           toString(file.nam[i,2]))
  }

  # save grange object
  save.GRange.object(gatc.all.m, save.nam)
  save.nam.p <- paste(save.nam, "plus", sep = "-")
  save.GRange.object(gatc.plus.m, save.nam.p)
  save.nam.m <- paste(save.nam, "minus", sep = "-")
  save.GRange.object(gatc.minus.m, save.nam.m)

  # get fragment and save it (not reads per GATC site but per GATC fragment)
  frag.gran <- make.fragmet.counts(list(gatc.plus.m, gatc.minus.m),
                                   genom.data = genom, resr.site = restsite,
                                   gatc.fragm = gatc.frag.m)
  save.nam.f <- paste(save.nam, "fragment", sep = "-")
  save.GRange.object(frag.gran, save.nam.f)

  # return GATC sites grange with all metadata reads. 2016.01.20 also frag.gran
  return(list(gatc.all.m, gatc.plus.m, gatc.minus.m, frag.gran))
}

make.log.fold.change <- function(granges, ids, norm.tr=T) {
  # normalize reads and calculate fold change for every sample control pair
  # (sample1/sum(sample1) / control1/sum(control1))
  # then calculate the log2 of the mean fold change of all sample control pairs
  # if norm.tr=F don’t do the normalization per total number of reads
  #
  # Args:
  #   granges: (granges)  granges with all reads in metadata
  #   ids: (list of 2) list with two vectors sampl.ids and control.ids
  #   norm.tr: (bool) if norm by total number or reads should be made or not
  #
  # Return:
  #   matrix with mean log2 fold change

  # get sample and control ids
  sampl <- ids[[1]]
  contr <- ids[[2]]

  # get metadata
  meta.datas <- elementMetadata(granges)

  # build result matrix ncol n/2 as we have always sample and control
  res.matrix <- matrix(NA, nrow = length(meta.datas[,1]),
                       ncol = (length(meta.datas[1,])/2))

  if (norm.tr) {
    # loop over all sample control pairs (has been tested to have equal numbers)
    for (i in 1:length(sampl)) {


      # add one to each sample and control position (no 0 in the sample)
      #temp.samp.1 <- meta.datas[,sampl[i]] + 1
      #temp.cont.1 <- meta.datas[,contr[i]] + 1

      # sum the read number per sample and control for normalization
      norm.sam <- sum(meta.datas[,sampl[i]])
      norm.con <- sum(meta.datas[,contr[i]])

      # take the sum
      #norm.sam <- sum(temp.samp.1)
      #norm.con <- sum(temp.cont.1)

      # change all 0 values to 1 as we don’t want to divide by zero
      temp.cont.1 <- meta.datas[,contr[i]]
      temp.cont.1[which(temp.cont.1 == 0)] <- 1

      # make normalized reads reads/sum(reads)/control/sum(controls)
      temp.reads <- (meta.datas[,sampl[i]] / norm.sam) / (temp.cont.1 / norm.con)
      #temp.reads <- (temp.samp.1 / norm.sam) / (temp.cont.1 / norm.con)
      res.matrix[,i] <- temp.reads
    }
  } else { # this is for data who are not binned and therefore have lot of 0
    for (i in 1:length(sampl)) {
      # change all 0 values to 1 as we don’t want to divide by zero and also
      # don’t want to loos the value if 0 / x (0 as sample)
      t.cont.1 <- meta.datas[,contr[i]]
      t.samp.1 <- meta.datas[,sampl[i]]

      # search all possition who control and sample is 0
      sample.nul <- which(t.samp.1 == 0)
      contr.null <- which(t.cont.1 == 0)

      # find positions with sample and control as 0
      both.null <- intersect(sample.nul, contr.null)

      # the lowest value should be equal to 1 if we have 0.nnnn.
      # give all 0 in the controls the lowest value of the samples
      min.va <- min(t.samp.1)
      if (min.va == 0) {
       min.va <- min(t.samp.1[t.samp.1 != min(t.samp.1)])
      }
      t.cont.1[contr.null] <- min.va / 2

      # the lowest value should be equal to 1 if we have 0.nnnn.
      min.va <- min(t.cont.1)
      if (min.va == 0) {
        min.va <- min(t.cont.1[t.cont.1 != min(t.cont.1)])
      }
      # set only samples to a value where the control is not 0
      sample.nul.no.c.nul <- sample.nul[-both.null]
      t.samp.1[sample.nul.no.c.nul] <- min.va / 2
      # make normalized reads reads / control
      temp.reads <- t.samp.1  / t.cont.1
      res.matrix[,i] <- temp.reads
    }
  }

  # mean all rows, make log2 transformation and return values.
  lo.fo.ca <- log2(rowMeans(res.matrix))
  lo.fo.ca[which(lo.fo.ca == -Inf)] <- 0 # cange -inf to 0

  return(lo.fo.ca)
}

################################################################################
#                   make and load grange with restriction site                 #
################################################################################
find.gatc.helper <- function(chr, spec, res.seq) {
  # find restriction site overlap in one chromosome.
  #
  # Args:
  #   chr: (string) the name of the chromosome
  #   spec: (BSGenome) the genome
  #   res.seq: (string) restriction site string
  #
  # Return:
  #   list with chromosome string and restriction site matches.

  # make a empty list
  chr.res.match <- list()

  # Assign pointer string from chr to list and add restriction site positions
  chr.res.match <- matchPattern(res.seq, getSeq(spec, names = chr))

  return(chr.res.match)
}

make.restriction.site.grange.from.scratch <- function(species,
                                                      rest.seq="GATC",
                                                      save.names,
                                                      chro.nams) {
  # make empty GATC grange form scratch using BSgenome data
  #
  # Args:
  #   species: (BSgenomes) the genome of the species
  #   rest.seq: (string) the sequence of the restriction site
  #   save.names: the path to save the empty GATC Grange
  #   chro.nams: (string or NULL) the names of wished chromosomes or NUll
  #              for all chromosomes
  #
  # Return:
  #   The GATC granges for all, plus, minus, fragments

  # get restriction sites out of the genome
  if (!is.null(chro.nams)) { # for subset of chromosomes

    # call find.gatc.helper for all chosen chromosomes to extract GATC sites
    # use sapply approach as bsapply has error when using only one chrom.
    res.sit.all <- sapply(X = chro.nams, FUN = find.gatc.helper,
                          spec = species, res.seq = rest.seq)

    # make GRange out of the list.
    res.sit.all <- as(as(res.sit.all, "RangesList"), "GRanges")

    # set seqinfo for chosen chromosomes
    seqinfo(res.sit.all) <- seqinfo(species)[chro.nams]

    # make plus and minus
    res.sit.plus   <- res.sit.all
    strand(res.sit.plus) <- Rle("+", length(res.sit.plus))
    res.sit.minus  <- res.sit.all
    strand(res.sit.minus) <- Rle("-", length(res.sit.minus))

  } else {# when all chromosomes are chosen

    # get restriction sites out of genome
    res.sit <- vmatchPattern(rest.seq, species, fixed = FALSE)
    seqinfo(res.sit) <- seqinfo(species)

    # make plus and minus
    res.sit.plus   <- res.sit[which(strand(res.sit) == "+")]
    res.sit.minus  <- res.sit[which(strand(res.sit) == "-")]

    # Make Genomic Range without strand information
    res.sit.all         <- res.sit.plus
    strand(res.sit.all) <- Rle("*", length(res.sit.all))
  }

  # make fragment grange
  fragment.grang <- make.fragment.from.gatc.grange(res.sit.all)

  # save GRages (all, plus, minus, fragments)
  save.GRange.object(res.sit.all, file.name = save.names[1])
  save.GRange.object(res.sit.plus, file.name = save.names[2])
  save.GRange.object(res.sit.minus, file.name = save.names[3])
  save.GRange.object(fragment.grang, file.name = save.names[4])

  # load granges for global use. (maybe better way to do this!)
  return(list(res.sit.all, res.sit.plus, res.sit.minus, fragment.grang))
}


make.restriction.site.grange <- function(species, rest.seq="GATC",
                                         save.location="granges",
                                         chr.nam=NULL) {
  # Test if empty GATC granges file exist: if not make it, if yes load it
  #
  # Args:
  #   species: (BSgenomes) the genome of the species
  #   rest.seq: (string) the sequence of the restriction site
  #   save.location: (string) folder name to search for saved Granges
  #   chr.nam: (string or NULL) name of the chromosomes or for all NULL
  #
  # Return:
  #   list of 4 Granges, GATC.all, plus only, minus only, fragments

  # make name to load/save granges
  prov  <- providerVersion(species)
  relda <- releaseName(species)
  org   <- organism(species)

  # also make a string for chromnames to add at the end of the Grange filename
  if (is.null(chr.nam)) {
    cr.n <- "all-chr"
  } else {
    cr.n  <- toString(chr.nam)
  }

  # concatenate the name parts
  full.name <- paste(org, prov, relda, cr.n)

  # delete whitespace, points and commas
  full.name.nw <- gsub(" ", "", full.name, fixed = TRUE) # delete ws
  rele.name.n1 <- gsub(".", "", full.name.nw, fixed = TRUE) # delete .
  rele.name.nw <- gsub(",", "", rele.name.n1, fixed = TRUE) # delete ,

  # add the suffix to the name
  sufix <- c("all.RData", "plus.RData", "minus.RData", "fragment.RData")

  # concatenate name and suffix and add the path string
  names.all <- paste(rele.name.nw, sufix, sep = "-")
  path.all  <- file.path(save.location, names.all)

  # test if granges already exist and load it if yes
  exist.file <- file.exists(path.all)

  if (sum(exist.file) < 4) { # if not all empty GATC granges are found make it
    message(paste(format(Sys.time(), "%H:%M:%S:"), "Find all", rest.seq,
                  "restriction sites in genome for:", cr.n))
    gatc.grangs <- make.restriction.site.grange.from.scratch(species = species,
                                                             rest.seq,
                                                             path.all, chr.nam)
  } else {# otherwise load it
    message(paste(format(Sys.time(), "%H:%M:%S:"), "Load", rest.seq,
                  "restriction sites for", cr.n))
    gatc.all   <- load.restriction.sites(path.all[1])
    gatc.plus  <- load.restriction.sites(path.all[2])
    gatc.minus <- load.restriction.sites(path.all[3])
    gatc.fragm <- load.restriction.sites(path.all[4])

    gatc.grangs <- list(gatc.all, gatc.plus, gatc.minus, gatc.fragm)
  }
  return(gatc.grangs)
}


normalize.grange.md <- function(gran.data) {
  # Normalize reads per bin by dividing by total number of reads
  #
  # Args:
  #   gran.data: (grange) the grange containing all metadata
  #
  # Return:
  #   grange with normalized (by total number) metadata
  #

  # make matrix out of grange
  ful.met <- as.matrix(elementMetadata(gran.data))
  # make vector with all column sums
  normalizer <- colSums(ful.met)

  # make new matrix
  # (this step is not needed could use ful.met in for loop instead)
  new.met.norm <- ful.met

  # normalize
  for (i in 1:length(ful.met[1,])) {
    new.met.norm[,i] <- ful.met[,i] / normalizer[i]
  }

  # delete metadata from grange and write normalized read in it
  elementMetadata(gran.data) <- NA
  elementMetadata(gran.data) <- new.met.norm

  return(gran.data)
}


log10.transfom.grange.md <- function(gran.data) {
  # log10 transformation of all grange metadata values
  #
  # Args:
  #   gran.data: (grange) the grange containing all metadata
  #
  # Return:
  #   grange with log10 transformed data
  #

  # make matrix out of grange
  ful.met <- as.matrix(elementMetadata(gran.data))

  # make new matrix
  # (this step is not needed could use ful.met in for loop instead)
  new.met.log <- ful.met

  # log10 transformation (could also don with sapply)
  for (i in 1:length(ful.met[1,])) {
    new.met.log[,i] <- log10(ful.met[,i])
  }

  # delete metadata from grange and write normalized read in it
  elementMetadata(gran.data) <- NA
  elementMetadata(gran.data) <- new.met.log

  return(gran.data)
}


################################################################################
#                     make fragment not GATC site counts                       #
################################################################################

# GATC count = every sequencing read starting (plus) or ending (minus) with a
#              GATC counts for this GATC site.
#
# fragment count = (fragment = distance between two GATC sites)
#                  count per fragment means the reads staring with a GATC (plus)
#                  and the next ending with a GATC (minus) count
#                  for one fragment.


make.fragment.from.gatc.grange <- function(gran) {
  # Expand start and end position of a GATC grange to make fragment grange
  # HELPER FUNCTION FOR: make.fragmet.counts
  #
  # Args:
  #   gran: (grange) grange with grange entries (ranges with width of 4 bp)
  #
  # Return:
  #   new grange with fragments (distance between two GATC sites )

  # get name and length of chromosomes
  chr.nam <- seqnames(gran)
  chr.len <- seqlengths(gran)

  # loop over chromosomes
  for (chr in seqnames(gran)@values) {

    # get sub grange (one chromosome)
    chr.gran <- gran[which(seqnames(gran) == chr)]

    ## expand ranges
    # first form 1 bp (start of chromosome) to first GATC site
    # from 1 bp to first GATC must be added to the metadata
    # as the GATC site is expanded to the right.
    first.gatc.pos <- chr.gran[1]
    start(first.gatc.pos[1]) <- 1

    # expand the length from the first GATC (end of the next GATC)
    # to the second last GATC.
    # move end position one step to the right
    # end position 2 -> now end position 3
    new.end.pos <- end(chr.gran[-1])

    # add end position to grange
    end(chr.gran[-length(chr.gran)]) <- new.end.pos

    # last fragment from last GATC site to length(chromosome)
    end(chr.gran[length(chr.gran)]) <- chr.len[chr]

    # add first fragment to grange
    ful.chr <- append(first.gatc.pos, chr.gran)

    # append chromosome grange to full grange or in the first iteration make
    # the final grange
    if (exists("return.gran")) {
      return.gran <- append(return.gran, ful.chr)
    } else {
      return.gran <- ful.chr
    }

    # remove gran (cleanup)
    rm(first.gatc.pos, chr.gran, ful.chr)
  }
  return(return.gran)
}


make.fragmet.counts <- function(grangs, genom.data, resr.site="GATC",
                                gatc.fragm) {
  # Make read counts per fragment from GATC plus and minus grange instead of
  # read counts per GATC.
  #
  # Args:
  #   grangs:    (list of GRange) GATC sites plus and minus with metadata
  #   genome.da: (BSgenome)      genome of the species to use
  #   rest.site: (string)        restriction site string
  #   gatc.fragm: (GRange)       empty GATC fragment GRange
  #
  # Return:
  #   new grange with counts per fragment (plus and minus separated)

  # extract names
  p.nam <- names(elementMetadata(grangs[[1]]))
  m.nam <- names(elementMetadata(grangs[[2]]))

  # loop over all chromosomes
  for (chr in seqnames(grangs[[1]])@values) {

    # get GATC reads for one chromosome
    chr.gran.p <- grangs[[1]][which(seqnames(grangs[[1]]) == chr)]
    chr.gran.m <- grangs[[2]][which(seqnames(grangs[[2]]) == chr)]

    # extract metadata
    plus.gran <- elementMetadata(chr.gran.p)
    minu.gran <- elementMetadata(chr.gran.m)

    # make empty lists and add first (minus list) and last fragment (plus list)
    plus.frag <- rbind(rep(NA, length(plus.gran[1,])), data.frame(plus.gran))
    minu.frag <- rbind(data.frame(minu.gran), rep(NA, length(minu.gran[1,])))

    # make new matrix for plus and minus in one list and fill in first 2 column.
    # plus GATC site 1 minus GATC site 1, plus GATC site 2 minus GATC site 2, ..
    result.metadata <- data.frame(plus.frag[,1], minu.frag[,1])
    result.names <- c(paste(p.nam[1], "p", sep = "."),
                      paste(m.nam[1], "m", sep = "."))

    # fill in reads if more than one sample is listed in the grange metadata:
    # plus 2 minus 2, plus 3 minus 3 .... plus n minus n for each read
    if (length(plus.gran[1,]) > 1) {
      for (i in 2:length(plus.gran[1,])) {
        result.metadata <- data.frame(result.metadata, plus.frag[,i],
                                      minu.frag[,i])
        result.names <- c(result.names, paste(p.nam[i], "p", sep = "."),
                          paste(m.nam[i], "m", sep = "."))
        names(result.metadata) <- result.names
      }
    }

    # make result in the first iteration and append it for other iterations
    if (exists("result.d.f")) {
      result.d.f <- rbind(result.d.f, result.metadata)
    } else {
      result.d.f <- result.metadata
    }
  }

  # add fragment reads to metadata of empty fragment GRange
  return.md <- gatc.fragm # do not add stuff to gatc.fragm
  elementMetadata(return.md) <- NULL # just to be save
  elementMetadata(return.md) <- result.d.f # add metadata to result
  return(return.md)

}


################################################################################
##                           save and load objects                              ##
################################################################################

save.GRange.object <- function(grange.ob, file.name="test") {
  # save a GRange object to a .RData file
  #
  # Args:
  #   grange.ob: (GRange) the grange object you want to save.
  #   file.name: (string) the path and filename
  #
  # Return:
  #   none (it saves the object using the given path and name)


  # test if .RData is in the name and if not concatenate it
  test <- grep(".RData", file.name)

  # concatenate .Rdata if not already in name
  if (length(test) == 0) {
    file.name <- paste(file.name, ".RData", sep = '')
  }

  # save grange object
  save(grange.ob, file = file.name)
}


load.restriction.sites <- function(object.names) {
  # load sved R objects. (.RData)
  #
  # Args:
  #   object.names: list of R object names to be loaded.
  #
  # Return:
  #   imported R object


  # test if ending with .Rdata concatenate the .Rdata ending
  test <- grep(".RData", object.names)

  # concatenate .Rdata if not already in name
  if (length(test) != 0) {
    f.name <- object.names
  } else {
    f.name <- paste(object.names, ".RData", sep = '')
  }

  # import R object
  obj <- local(get(load(f.name)))

  return(obj)
}

################################################################################
# save metadata values of granges to a table                                   #
################################################################################
save.table.from.grange <- function(gran, file.pa) {
  # write table with metadata of granges GATC reads (strand less)
  #
  # Args:
  #   gran:    (GRange) grange with metadata
  #   file.pa: (string) path and name to save the table
  #
  # Return:
  #   none (save table to hard disk)

  # make data frame out of grange metadata and assign name
  met.df <- data.frame(elementMetadata(gran))
  colnames(met.df) <- names(elementMetadata(gran))

  # get start position and chromosome
  ind.df <- data.frame("chr" = seqnames(gran), "start.position" = start(gran))

  # concatenate both
  met.df.ful <- cbind(ind.df, met.df)

  # write table to hard disk
  write.table(met.df.ful, file.pa)
}


save.table.from.grange2 <- function(gran, gran.p, gran.m, file.pa) {
  # write table with metadata of granges
  # given in ob.name with plus and minus strand reads separately
  #
  # Args:
  #   gran:    (GRange) grange with metadata
  #   gran.p:  (GRange) grange with metadata plus strand
  #   gran.m:  (GRange) grange with metadata minus strand
  #   file.pa: (string) path and name to save the table
  #
  # Return:
  #    none (save table to hard disk)

  # all
  met.df <- data.frame(elementMetadata(gran))
  colnames(met.df) <- names(elementMetadata(gran))
  # plus
  met.p.df <- data.frame(elementMetadata(gran.p))
  colnames(met.p.df) <- names(elementMetadata(gran.p))
  # minus
  met.m.df <- data.frame(elementMetadata(gran.m))
  colnames(met.m.df) <- names(elementMetadata(gran.m))

  # test if more than one sample in dataframe and add it if true
  if (length(met.m.df[1,]) > 1) {

    # concatenate data frame for first sample (plus and minus)
    name.all  <- colnames(met.df[1,])
    name.plus <- paste(colnames(met.p.df[1,]), "plus", sep = ".")
    name.minu <- paste(colnames(met.m.df[1,]), "minus", sep = ".")
    met.all.in.one <- data.frame(met.df[,1], met.p.df[,1], met.m.df[,1])
    colnames(met.all.in.one) <- c(name.all[1], name.plus[1], name.minu[1])

    # loop from the 2. sample to the end
    for (i in 2:length(met.m.df[1,])) {
      temp.df <- data.frame(met.df[,i], met.p.df[,i], met.m.df[,i])
      colnames(temp.df) <- c(name.all[i], name.plus[i], name.minu[i])
      met.all.in.one <- cbind(met.all.in.one, temp.df)
    }

  } else {# if only one sample
    name.all  <- colnames(met.df)
    name.plus <- paste(colnames(met.p.df), "plus", sep = ".")
    name.minu <- paste(colnames(met.m.df), "minus", sep = ".")
    met.all.in.one <- data.frame(met.df[,1], met.p.df[,1], met.m.df[,1])
    colnames(met.all.in.one) <- c(name.all[1], name.plus[1], name.minu[1])
  }

  # get start position and chromosome
  ind.df <- data.frame("chr" = seqnames(gran), "start.position" = start(gran))

  # concatenate both
  met.df.ful <- cbind(ind.df, met.all.in.one)

  # write table
  write.table(met.df.ful, file.pa)
}

################################################################################
##                              function to plot                              ##
################################################################################
plot.all.chrom.in.one.2 <- function(grange, chr.nam, main.name="a plot",
                                    y.name="a y axis", y.lim=c(0,25000),
                                    file.name="resolution.png",
                                    plot.type="h",
                                    plot.size=c(12,8), met.dat.name="value",
                                    d.type="png", font.size=25, plot.col=1,
                                    log.t=F) {
  # Print .png / .pdf image of GRange restriction sites resolutions.
  # by plot all chromosome in one continuous plot.
  #
  # Args:
  #   grange:      (GRange) with metadata
  #   chr.nam:     (vector  of strings) name of the chromosomes to be plotted
  #   file.name:   (string) .png or .pdf file name
  #   plot.type:   (string) plot type e.g. h for histogram (r plot notation)
  #   plot.size:   (vector of two int) width and height of the plot
  #   y.lim        (vector of two int) y axis limitations
  #   y.name       (string) name of the y axis
  #   main.name    (string) name for the plot
  #   met.dat.name (string) name of the metadata column to be plotted
  #   d.type       (string or NA) name of the plot type (png or pdf)
  #                               if pdf w/h must be inches
  #                               (if NA then the plot will not be saved)
  #   font.size    (int): define the size of the font
  #   plot.col:    (int, string) the id or name for the plot color
  #   log.t:       (bool) if T values get log 10 transformed
  #
  # Return:
  #   none, safes a png file in the working directory

  # make a continuous list over all chromosomes
  coordinates <- fuse.chr(grange)

  # get metadata from grange either do a log10 transformation
  if (log.t) {
    metadata <- log10(elementMetadata(grange)[,met.dat.name])
    y.lim[2]  <- log10(y.lim[2]) # also log transform y axis limitations
    temp.min <- min(metadata)    # get minimum value
    y.lim[1]  <- floor(temp.min) # also log transform y axis limitations

  } else {
    metadata <- elementMetadata(grange)[,met.dat.name]
  }


  # choose plot data type
  if (d.type == "png") {
    png(file.name, width = plot.size[1], height = plot.size[2])
  }
  if (d.type == "pdf") {
    pdf(file.name, width = plot.size[1], height = plot.size[2])
  }

  # get maximal plot length (+1 to have a spacer, not really necessary)
  x.lim <- c(0, max(coordinates) + 1)

  # plot the graph
  opar = par(ps = font.size)
  plot(coordinates, metadata, type = plot.type, ylim = y.lim, xlim = x.lim,
       ylab = y.name, main = main.name, col = plot.col)

  # make vertical lines add also a red lie at the end of the last chromosome
  line.positions <- c(coordinates[start(seqnames(grange))],
                      max(coordinates[end(seqnames(grange))]))
  abline(v = line.positions, col = rgb(1,0,0,1/2)) # red half transparent

  # get chromosome names
  c.nams <- seqnames(grange)@values
  # write names to plot, (-last to exclude last line position)
  text(x = line.positions[-length(line.positions)],
       y = rep(y.lim[2], length(c.nams)), c.nams, pos = 4) # pos=4 = rigth

  opar

  # end ploting
  if (d.type == "png" | d.type == "pdf") {
    dev.off()
  }
}


plot.correlation <- function(grange, samp, cont, type="a", timesp=0,
                             bin.name="gatc", exp.name="test", qs.fold.name) {
  # plot the correlation of grange metadata (all columns vs. the others)
  # and save correlation in a file
  #
  # Args:
  #   grange: (GRange) with metadata
  #   samp:   (vector) index of sample  column positions in grange metadata
  #   cont:   (vector) index of control column positions in grange metadata
  #   type:   (string) if all, only sample and only control should be plot
  #                    3 values allowed "a" = all, "s" = sample, "c"=control
  #   bin.name (string) the name of the bin
  #   exp.name (string) the name of the experiment
  #   qs.fold.name (string) the name of the folder to save QC data
  #
  # Return:
  #   none, plot in current device

  # get metadata
  if (type == "s") {
    met.data <- elementMetadata(grange)
    met.data <- met.data[,samp]
  } else if (type == "c") {
    met.data <- elementMetadata(grange)
    met.data <- met.data[,cont]
  } else {
    met.data <- elementMetadata(grange)
  }

  # get correlation
  cor.gatc <- cor(as.matrix(met.data))

  # save correlation as file (only for type a)
  if (type == "a") {
    f.nam <- paste("correlation-",exp.name, "-", bin.name, "-",
                   timesp, ".txt", sep = "")
    rep.path <- file.path(qs.fold.name, f.nam)
    write.table(cor.gatc, file = rep.path, eol = "\r\n" )
  }

  # cl.lim (test if negative correlation and if true make cl.lim = -1,1)
  if (min(cor.gatc) < 0) {
    cor.lim <- c(-1,1)
    # color panel for the correlations -1 to 1
    colpal <- colorRampPalette (c("#7F0001", "#FFA501", "#808081", "#045A8E",
                                  "#00007E", "#045A8D", "lightgray",
                                  "orange", "#7F0000"))(300)
  } else {
    cor.lim <- c(0,1)
    # color panel for the correlations 0 to 1
    colpal <- colorRampPalette (c("#00007F","#045A8D","lightgray", "orange",
                                  "#7F0000","#00007F","#045A8D","lightgray",
                                  "orange", "#7F0000"))(300)
  }

  # first plot for the color plot
  corrplot(cor.gatc, method = "color", type = "lower", hclust = "complete",
           cl.lim = cor.lim, col = colpal, order = "hclust", outline = FALSE,
           addgrid.col = "black", plotCI = "n", tl.pos = "lt", diag = TRUE)
  # second plot for the numbers
  corrplot(cor.gatc,  method = "number", tl.pos = "n", type = "upper",
           diag = FALSE, hclust = "complete", order = "hclust", outline = FALSE,
           addgrid.col = "black", plotCI = "n", add = TRUE, cl.pos = "n",
           tl.cex = 8)

  # make list with both corplots in it
  #cor.plots <- c(first.plot, second.plot)
  #return(cor.plots)
}


plot.normalized.gatc.number.one.vs.another.log <- function(x.data, y.data,
                                                           nams=c("a", "b"),
                                                           bin.id) {
  # plot normalized numbers of GATC as smoothScatter and add correlation to it
  #
  # Args:
  #   x.data: (vector) value of sample 1
  #   y.data: (vector) value of sample 2
  #   nams:   (vector of 2) name of x and y sample
  #   bin.id: (int) length of the bin used for this data set
  #
  # Return:
  #   none, plot in current device

  # normalize
  nor.x <- x.data / sum(x.data)
  nor.y <- y.data / sum(y.data)

  # change 0 to NA for linear model
  nor.x.na <- nor.x
  nor.x.na[nor.x.na == 0] <- NA
  nor.y.na <- nor.y
  nor.y.na[nor.y.na == 0] <- NA

  # make linear model
  fit <- lm(log10(nor.y.na) ~ log10(nor.x.na), na.action = "na.omit")

  # make log fold scatterplot
  smoothScatter(log10(nor.x), log10(nor.y),
                main = paste("log10", nams[1],  "vs", nams[2], bin.id, "bins"),
                ylab = paste("log10", nams[2]), xlab = paste("log10", nams[1]),
                pch = ".", xlim = c(-6,0), ylim = c(-6,0))
  # plot linear model
  abline(fit, col = "red")
  # add legend with correlation value
  legend("topleft", paste("cor =", round(cor(nor.x, nor.y),3)))
}


################################################################################
##                      helper function to get genome                         ##
################################################################################
test.if.chrom.exist <- function(chr.names, bs.gen) {
  # test if chromosome names exist
  #
  # Args:
  #   chr.names: (vector of strings) chromosome names if NULL, all
  #   bs.gen:    (string) the name of the package (species)
  #
  # Return:
  #   the BSgenomeobject of the species if available


  # test if chromosome names match given values if not set NULL(NULL=take all)
  if (!is.null(chr.names)) {
    chr.n.test <- chr.names  %in%  seqnames(bs.gen)
    if (sum(chr.n.test) != length(chr.names)) {
      stop(paste("Chromosome name:",
                 toString(chr.names[which(!chr.n.test)]),
                 "not found in genome chromosome name set.",
                 "Choose only chromosome listed next:",
                 toString(seqnames(bs.gen)),
                 "from", toString(bs.gen@pkgname)))
    }
  }
}


################################################################################
##                              get genome form BSgenome                      ##
################################################################################
get.genome <- function(bs.genome.name="BSgenome.Celegans.UCSC.ce10",
                       chromosome.names=NULL) {
  # get genome form BSgenome, download, install, Xstring object
  #
  # Args:
  #   bs.genome.name:   (string) the name of the package (species)
  #   chromosome.names: (vector of strings) chromosome names if NULL, all
  #
  # Return:
  #   the BSgenomeobject of the species if available
  #
  # TODO: (DR) forge genome from fasta ... home made genome library.

  # test if package is installed
  if (bs.genome.name %in% rownames(installed.packages())) {
    # load the bsgenome package
    load.pakage.main("", bs.genome.name)

    # load the genome
    bs.genome <- getBSgenome(bs.genome.name, masked = FALSE)

    # test if chromosome names match given values if not set NULL(NULL=take all)
    test.if.chrom.exist(chromosome.names, bs.genome)

    return(bs.genome)
  }

  # test if package (genome) exist
  genome.list <- available.genomes()
  is.name.in.list <- bs.genome.name %in% genome.list
  if (!is.name.in.list) {
    stop(paste("Your genome", bs.genome.name, "is not supported"))
  } else {
    load.pakage.main("", bs.genome.name)
    bs.genome <- getBSgenome(bs.genome.name, masked = FALSE)

    # test if chromosome names match given values if not set NULL(NULL=take all)
    test.if.chrom.exist(chromosome.names, bs.genome)

    return(bs.genome)
  }
}

################################################################################
##                              read fastq.gz files                           ##
################################################################################
read.fastq <- function(file.path, fasta.for="fastq", seq.type="DNA") {
  # read a fastq (with quality string) or fasta file
  #
  # Args:
  #   file.path: (string) the path and filename of the fastq file
  #   fasta.for: (string) the format (fasta, fastq)
  #   seq.type:  (string) the type of the data, DNA, RNA, aa, b for anything
  #
  # Return:
  #   Xstrinng object with the reads and qc (fastq)

  # read fastq (with quality string) or fasta from file
  if (fasta.for == "fastq") {
    return(readFastq(file.path, withIds = T))
  } else {
    return(readDNAStringSet(filepath = file.path, format = fasta.for))
  }
}


################################################################################
##                              make QC report                                ##
################################################################################
make.qc <- function(file.list, time.stamp, exp.name, rep.dir) {
  # make QC (fastqc)
  #
  # Args:
  #   file.list: (data.frame) the war list data frame with 1=path, 2=name,
  #                           and optionally 3=group
  #                           the first column would be enough
  #   time.stamp: (Sys.time) the timestamp for the name of the saved file
  #   exp.name:   (string) the name of the experiment
  #   rep.dir:    (string) name of the directory to put in the reports
  #
  # Return:
  #    none (save report to file)

  # make filenames list
  file.names <- sapply(file.list[,1], as.character)

  # make pdf report path
  qc.file.name <- paste(exp.name, "-", time.stamp, ".pdf", sep = "")
  qc.full.path <- file.path(rep.dir, qc.file.name)

  # make qc (also be used with multicore)
  ret.report <- qQCReport(input = file.names, pdfFilename = qc.full.path)
}


################################################################################
##                              cutadapt homemade script                      ##
################################################################################
make.cutadapt.fastq <- function(read.list, adapter="CGCGGCCGAG", error.a=1,
                                out.name="test", compr=T, res.site="GATC",
                                dir.name.ca) {
  # cut adapter from the 5' end if the GATC sequence is following after
  # the adapter. This script cuts the adapter not only on the flanking regions.
  # this is for fastq files
  #
  # Args:
  #   read.list: (Xstring) all the reads form 1 sequencing sample
  #   adapter:   (string)  the adapter sequence
  #   errors.a:  (integer) the number of errors allowed
  #   out.name:  (string)  the name to save the cut sequence as fastq
  #   compr:     (bool)    save the file T=compressed (gzip), F=uncompressed
  #   res.site:  (string)  the motif of the restriction site
  #   dir.name.ca (stirng) the name of the folder for cut read files
  #
  # Return:
  #   (Xstrinng) only trimmed reads starting with GATC (saves cut fasta in file)

  # get timestamp to calculate the run time
  start.time <- Sys.time()

  # get number of raw reads
  raw.read.length <- length(read.list)

  # concatenate adapter + restriction site
  adapter.p <- paste(adapter, res.site, sep = "")

  # get quality strings
  quality.list <- quality(read.list)@quality
  # get sequence and also overwrite the read.list
  read.list    <- sread(read.list)

  # test if all sequences have the same length and give warning if not
  if (mean(width(read.list)) != width(read.list[1])) {
    warning("Your read list contains multiple read lengths!")
  }

  # get length of every read
  seq.length <- width(read.list)

  # find adapter + GATC position using vmatchpattern (no flanking regions)
  adapter.position <- vmatchPattern(adapter.p, read.list,
                                    max.mismatch = error.a)

  # get the end position of the matched adapter
  adapter.end.pos <- endIndex(adapter.position)

  # remove adapter.position (RAM friendly)
  rm(adapter.position)

  # change NULL to 0 as atomic vectors cannot have NULL
  # returns the full length of the sequence when no adapter is found.
  adapter.end.pos[unlist(lapply(adapter.end.pos, is.null))] <- as.integer(0)

  # vmatchpattern is able to return more than one element when
  # there are more matches in the DNA string.
  # get every first entry of each individual list
  # (could also exclude reads with multiple matches)
  adapter.end.pos.first <- lapply( adapter.end.pos, '[[', 1)

  # remove adapter.end.pos (RAM friendly)
  rm(adapter.end.pos)

  # make a atomic vector from list (faster)
  adapter.vect  <- as.vector(adapter.end.pos.first, mode = "integer")

  # remove adapter.end.pos.first (RAM friendly)
  rm(adapter.end.pos.first)

  # remove the restriction site form the adapter so that the read start with it
  adapter.vect.r <- adapter.vect - nchar(res.site)
  # set all adapter end positions to 0 if they are negative
  # (if no adapter was found 0 + -4 = -4 therefore set it back to 0)
  # 5' overlapping and the sequence starts with the restriction site(no adapter)
  adapter.vect.r[which(adapter.vect.r < 0)] <- 0

  # set all adapter end positions to 100 if they are 3' (right) overlap
  # this reads are not useful as we do not have any sequence starting with
  # a restriction site left. (only the adapter is overlapping the 3' end)
  adapter.vect.r[which(adapter.vect.r > seq.length)] <- seq.length

  # use the adapter end index - restriction site to truncate the read sequences
  # +1 as we have from 0 to 100 (or 50) and not from 1 to 100 in this object
  cut.read.list <- subseq(read.list, start = adapter.vect.r + 1) # for reads
  cut.qual.list <- subseq(quality.list, start = adapter.vect.r + 1) # quality

  # remove adapter.vec.r, read.list, and quality list   (RAM friendly)
  rm(adapter.vect.r, read.list, quality.list)

  # make a new list with only cut sequences and also cut the quality list
  cut.read.list.red <- cut.read.list[which(width(cut.read.list) != seq.length)]
  cut.qual.list.red <- cut.qual.list[which(width(cut.read.list) != seq.length)]

  # find adapter in 5' flanking region (overlapping)
  # get uncut sequences and also get quality list for uncut reads
  uncut.read.list <- cut.read.list[which(width(cut.read.list) == seq.length)]
  uncut.qual.list <- cut.qual.list[which(width(cut.read.list) == seq.length)]

  # make new seq.length list with only uncut reads
  seq.le.u.c <- seq.length[which(width(cut.read.list) == seq.length)]

  # use trimLTPattern to map adapter in flanking regions
  temp.trim <- trimLRPatterns(Lpattern = adapter.p, Rpattern = "",
                              uncut.read.list, ranges = T)

  # add the restriction site to the start of the remaining sequence
  start(temp.trim) <- start(temp.trim) - nchar(res.site)
  # remove all entries with negative starts and also remove adapter matching
  # only the first 3 or 4 bases of the read.
  start(temp.trim)[which(start(temp.trim) < 4)] <- 1
  # 4 as only take sequences with match 3 bases on the adapter
  # (increasing of this number truncates the max length of
  # reads 4 = max read length 97bp 5=96bp)

  # remove adapter sequence from sequences and also the quality data
  flanking.adapt <- subseq(uncut.read.list, start = start(temp.trim))
  flanking.quali <- subseq(uncut.qual.list, start = start(temp.trim))
  # get only sequences who are truncated
  flanking.adapt.red <- flanking.adapt[which(width(flanking.adapt)
                                             != seq.le.u.c)]
  flanking.quali.red <- flanking.quali[which(width(flanking.adapt)
                                             != seq.le.u.c)]

  # save fastq file compressed compr==T or uncompressed
  save.shre <- ShortReadQ(sread = c(cut.read.list.red, flanking.adapt.red),
                          quality = c(cut.qual.list.red, flanking.quali.red))
  if (compr) {
    export.name <- paste(out.name, "cut.fastq.gz", sep = "-")
    export.path <- file.path(getwd(), dir.name.ca, export.name, fsep = "/")

    # save compressed
    writeFastq(save.shre, export.path, mode = 'w', compress = compr)

  } else {
    export.name <- paste(out.name, "cut.fastq", sep = "-")
    export.path <- file.path(getwd(), dir.name.ca, export.name, fsep = "/")

    # save uncompressed
    writeFastq(save.shre, export.path, mode = 'w', compress = compr)
  }


  ## prepare report
  # get lengths of remaining sequences
  seq.len.c <- table(width(cut.read.list.red))  # from full adapter alignment
  seq.len.f <- table(width(flanking.adapt.red)) # from partial adapter alignment
  # concatenate both tables.
  seq.len <- c(seq.len.c, seq.len.f)

  # calculate number of reads with and without adapter
  no.adapter <- raw.read.length - (length(cut.read.list.red) +
                                     length(flanking.adapt.red))
  wi.adapter <- length(cut.read.list.red) + length(flanking.adapt.red)
  # calculate percentage
  adapt.presentage <- wi.adapter * 100 / raw.read.length


  ## save report file as .txt table.
  # make data.frame for easy report saving
  result.report <- as.data.frame(seq.len)
  colnames(result.report) <- "count"
  # make file name and set path
  export.name.txt <- paste(out.name, "cut-read-lengths.txt", sep = "-")
  export.path.txt <- file.path(getwd(),
                               dir.name.ca, export.name.txt, fsep = "/")
  # write header
  write.table(paste("length", "count", sep = "\t"),
              file = export.path.txt, eol = "\r\n", col.names = F,
              row.names = F, quote = F)
  # add length table data to the file
  write.table(seq.len, file = export.path.txt, sep = "\t",
              eol = "\r\n", quote = F, col.names = F, append = T)


  # get timestamp to calculate run time
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "auto")

  # write some stats (concatenate text and values for printing in file)
  total.r   <- paste("Total number of reads:", raw.read.length, sep = "\t")
  total.r.l <- paste("length of reads in bp:", toString(unique(seq.length))
                     , sep = "\t")
  total.r.a <- paste("Reads with adapter and GATC sequence:", wi.adapter,
                     "(%)", round(adapt.presentage, digits = 4),  sep = "\t")
  total.r.n <- paste("Remaining reads which didn't align to adapter and GATC",
                     no.adapter, sep = "\t")
  total.rut <- paste("Time used for file import and adapter removal:",
                     format(duration), sep = "\t")

  # concatenate report to single string
  some.stats <- paste(total.r, total.r.l, total.r.a, total.r.n, total.rut,
                      sep = "\r\n")

  # set name, path and write report (number and % of cut reads) to file
  export.name.txt2 <- paste(out.name, "cut-report.txt", sep = "-")
  export.path.txt2 <- file.path(getwd(), dir.name.ca, export.name.txt2,
                                fsep = "/")
  write(some.stats, file = export.path.txt2)

  # return the path of the fastq file and (wi.adapter) number of adapter found
  return(c(export.name, wi.adapter))
}

make.cutadapt.fasta <- function(read.list, adapter="CGCGGCCGAG", error.a=1,
                                out.name="test", compr=T, res.site="GATC",
                                dir.name.ca) {
  # cut adapter from the 5' end if the GATC sequence is following after
  # the adapter also when the adapter sequence is flanking the 5' end. for fasta
  #
  # Args:
  #   read.list: (Xstring) all the reads form 1 sequencing sample
  #   adapter:   (string)  the adapter sequence
  #   errors.a:  (integer) the number of errors allowed
  #   out.name:  (string)  the name to save the cut sequence as fastq
  #   compr:     (bool)    T for compressed (gzip), F for uncompressed
  #   res.site:  (string)  the string of the restriction site
  #   dir.name.ca(stirng)  the name of the folder for cut read files
  #
  # Return:
  #   (Xstrinng) only trimmed reads starting with GATC (saves cut fasta in file)

  # get timestamp to calculate the run time
  start.time <- Sys.time()

  # concatenate adapter + restriction site
  adapter.p <- paste(adapter, res.site, sep = "")

  # test if all sequences have the same length
  if (mean(width(read.list)) != width(read.list[1])) {
    warning("Your read list contains multiple read lengths!")
  }

  # get length of the sequences
  seq.length <- width(read.list)

  # find adapter position using vmatchpattern (no flanking regions)
  adapter.position <- vmatchPattern(adapter.p, read.list,
                                    max.mismatch = error.a)

  # get the end position of the matched adapter
  adapter.end.pos <- endIndex(adapter.position)

  # remove adapter.position (RAM friendly)
  rm(adapter.position)

  # change NULL to 0 as atomic vectors cannot have NULL
  # returns the full length of the sequence when no adapter is found.
  adapter.end.pos[unlist(lapply(adapter.end.pos, is.null))] <- as.integer(0)

  # vmatchpattern is able to return more than one element when
  # there are more matches in the DNA string.
  # get every first entry of each individual list
  adapter.end.pos.first <- lapply( adapter.end.pos, '[[', 1)

  # make a atomic vector from list (faster)
  adapter.vect  <- as.vector(adapter.end.pos.first, mode = "integer")

  # remove the restriction site form the adapter as this should not be removed
  adapter.vect.r <- adapter.vect - nchar(res.site)
  # set all adapter end positions to 0 if they are negative
  # (if no adapter was found 0 + -4 = -4 therefore set it back to 0)
  # 5' overlapping and the sequence starts with the restriction site(no adapter)
  adapter.vect.r[which(adapter.vect.r < 0)] <- 0

  # set all adapter end positions to 100 if they are 3' overlap
  # this reads are not useful as we do not have any sequence starting with
  # a restriction site left. (only the adapter is overlapping the 3' end)
  adapter.vect.r[which(adapter.vect.r > seq.length)] <- seq.length

  # use the adapter end index - restriction site to truncate the read sequences
  # +1 as we have from 0 to 100 (or 50) and not from 1 to 100 in this object
  cut.read.list <- subseq(read.list, start = adapter.vect.r + 1)

  # make a new list with only cut sequences
  cut.read.list.red <- cut.read.list[which(width(cut.read.list) != seq.length)]

  # find adapter in 5' flanking region (overlapping)
  # get uncut sequences
  uncut.read.list <- cut.read.list[which(width(cut.read.list) == seq.length)]

  # make new seq.length list with only uncut reads
  seq.le.u.c <- seq.length[which(width(cut.read.list) == seq.length)]

  # use trimLTPattern to map adapter to flanking regions
  temp.trim <- trimLRPatterns(Lpattern = adapter.p, Rpattern = "",
                              uncut.read.list, ranges = T)

  # add the restriction site to the start of the remaining sequence
  start(temp.trim) <- start(temp.trim) - nchar(res.site)
  # remove all entries with negative starts and also remove adapter matching
  # only the first 3 or 4 bases of the read.
  start(temp.trim)[which(start(temp.trim) < 4)] <- 1
  # 4 as only take sequences with match 3 bases on the adapter
  # (increasing of this number truncates the max length of
  # reads 4 = max read length 97bp 5=96bp)

  # remove adapter sequence from sequences
  flanking.adapt <- subseq(uncut.read.list, start = start(temp.trim))
  # get only sequences who are truncated
  flanking.adapt.red <- flanking.adapt[which(width(flanking.adapt)
                                             != seq.le.u.c)]

  # save fasta file compressed compr==T or uncompressed
  if (compr) {
    export.name <- paste(out.name, "cut.fasta.gz", sep = "-")
    export.path <- file.path(getwd(), dir.name.ca, export.name, fsep = "/")
    # save cut reads to compressed file
    writeXStringSet(c(cut.read.list.red, flanking.adapt.red), export.path,
                    format = "fasta", compress = T)
  } else {# uncompressed
    export.name <- paste(out.name, "cut.fasta", sep = "-")
    export.path <- file.path(getwd(), dir.name.ca, export.name, fsep = "/")
    # save cut reads to uncompressed file
    writeXStringSet(c(cut.read.list.red, flanking.adapt.red), export.path,
                    format = "fasta")
  }


  ## prepare report
  # get lengths of remaining sequence
  seq.len.c <- table(width(cut.read.list.red))    # from full adapter alignment
  seq.len.f <- table(width(flanking.adapt.red)) # from partial adapter alignment
  # concatenate both tables.
  seq.len <- c(seq.len.c, seq.len.f)

  # calculate number of reads with and without adapter
  no.adapter <- length(read.list) - (length(cut.read.list.red) +
                                       length(flanking.adapt.red))
  wi.adapter <- length(cut.read.list.red) + length(flanking.adapt.red)
  # calculate percentage
  adapt.presentage <- wi.adapter * 100 / length(read.list)


  ## save report file as .txt table.
  # make data.frame for easy report saving
  result.report <- as.data.frame(seq.len)
  colnames(result.report) <- "count"
  # make file name and set path
  export.name.txt <- paste(out.name, "cut-read-lengths.txt", sep = "-")
  export.path.txt <- file.path(getwd(), dir.name.ca, export.name.txt,
                               fsep = "/")
  # write header
  write.table(paste("length", "count", sep = "\t"),
              file = export.path.txt, eol = "\r\n", col.names = F,
              row.names = F, quote = F)
  # add cut read length table data to the file
  write.table(seq.len, file = export.path.txt, sep = "\t",
              eol = "\r\n", quote = F, col.names = F, append = T)

  # get timestamp for calculating runtime
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "auto")

  # write some stats (concatenate text and values for printing in file)
  total.r   <- paste("Total number of reads:", length(read.list), sep = "\t")
  total.r.l <- paste("length of reads in bp:", toString(unique(seq.length))
                     , sep = "\t")
  total.r.a <- paste("Reads with adapter and GATC sequence:", wi.adapter,
                     "(%)", round(adapt.presentage, digits = 4), sep = "\t")
  total.r.n <- paste("Remaining reads which didn't align to adapter and GATC",
                     no.adapter, sep = "\t")
  total.rut <- paste("Time used for file import and adapter removal:",
                     format(duration), sep = "\t")

  # concatenate stats to single string
  some.stats <- paste(total.r, total.r.l, total.r.a, total.r.n, total.rut,
                      sep = "\r\n")

  # set name, path and write report to file
  export.name.txt2 <- paste(out.name, "cut-report.txt", sep = "-")
  export.path.txt2 <- file.path(getwd(), dir.name.ca, export.name.txt2,
                                fsep = "/")
  write(some.stats, file = export.path.txt2)

  # return the path of the saved fastq file
  return(c(export.name, wi.adapter))
}



################################################################################
##                              bowtie                                        ##
################################################################################
make.bowtie <- function(fastq.files, multi.core=F, genome.name, pr.name,
                        dir.name.bt, m.hits=1) {
  # map all cut sequences to genome using bowtie (QuasR pipeline)
  #
  # Args:
  #   fastq.files: (string) path to the .txt folder containing all the cut
  #                         fastq files and a name (tap delimit)
  #   multi.core: (Boolean) T = enables multicore functionality
  #   genome.name: (string) the name of the BSgenome package
  #   pr.name:     (string) the name of the project
  #   dir.name.bt: (string) the name of the folder to store the .sam/.bam reads
  #   m.hits:         (int) the number of maximal hits allowed
  #
  # Return:
  #   save location of mapped .sam/.bam file

  # the actual mapping using rbowtie (either multicore or single core)
  if (multi.core) {
    # test how many cores are available
    num.of.cores <- detectCores()
    # make cluster for multicore
    clusts <- makeCluster(num.of.cores)
    # the bowtie mapping
    test <- qAlign(fastq.files, genome  =  genome.name,
                   aligner = "Rbowtie", maxHits = m.hits,
                   alignmentsDir = file.path(getwd(), dir.name.bt),
                   clObj =  clusts, projectName = pr.name)
    stopCluster(clusts) # stop cluster
  } else {
    # the bowtie mapping without multicore
    test <- qAlign(fastq.files, genome  =  genome.name,
                   aligner = "Rbowtie", maxHits = m.hits,
                   alignmentsDir = file.path(getwd(), dir.name.bt),
                   projectName = pr.name)
  }

  message(paste(format(Sys.time(), "%H:%M:%S:"), "All sequences mapped to",
                genome.name ))
  message(paste(format(Sys.time(), "%H:%M:%S:"),
                ".bam files are saved in 'bowtie' subfolder"))

  # return position of alignments (mapped .sam files)
  ret.d.f <- alignments(test)
  return(ret.d.f[1])
}

################################################################################
##                               test  input                                  ##
################################################################################
test.input.data <- function(input.file, inputs) {
  # test if the input data is valid
  #
  # Args:
  #   input.file: (data.frame) the path, names, and type of the files.
  #   inputs: (list with objects) all the objects to test
  #           multi.core, exp.name, bin.len, adapt.seq, errors, input.format,
  #           qc, species, chr.names, restr.seq, normali, log10t
  #
  # Return:
  #   atomic vector with indexes of sample.ids, controll.ids

  ## test if input are valid
  # test multi.core
  if (!is.logical(inputs[[1]])) {
    stop(paste("multi.core:", toString(inputs[[1]]),
               "must be Boolean (logical, TRUE or FALSE)"))
  }
  # test exp.name
  if (!is.character(inputs[[2]])) {
    stop(paste("exp.name:", toString(inputs[[2]]),
               "must be string: ('string', nostring)"))
  }
  # test bin.len
  if (!is.numeric(inputs[[3]])) {
    stop(paste("bin.len:", toString(inputs[[3]]),
               "must be int: (n,...,-3,-2,-1,0,1,2,3,..,n)"))
  }
  # test adapt.seq
  if (!is.character(inputs[[4]])) {
    stop(paste("adapt.seq:", toString(inputs[[4]]),
               "must be string: ('string', nostring)"))
  }
  # test errors
  if (!is.numeric(inputs[[5]])) {
    stop(paste("errors:", toString(inputs[[5]]),
               "must be int: (n,...,-3,-2,-1,0,1,2,3,..,n)"))
  }

  # test input.format
  if (inputs[[6]] != "fasta" & inputs[[6]] != "fastq") {
    if (substring(inputs[[6]], 1, 1) != "b") {
      if (substring(inputs[[6]], 1, 1) != "g") {
        stop(paste("input.format:", toString(inputs[[6]]),
                   "must be 'fastq', 'fasta', 'b' or 'g' of type",
                   "character string"))
      }
    }
  }

  # test qc
  if (!is.logical(inputs[[7]])) {
    stop(paste("qc:", toString(inputs[[7]]),
               "must be Boolean (logical, TRUE or FALSE)"))
  }
  # test species
  if (!is.character(inputs[[8]])) {
    stop(paste("species:", toString(inputs[[8]]),
               "must be string (character): ('string', nostring)"))
  }
  # test chr.names
  if (!is.null(inputs[[9]])) {
    if (!is.character(inputs[[9]])) {
      stop(paste("chr.names:", toString(inputs[[9]]),
                 "must be vector of character: c('string1', 'string2')",
                 "or NULL"))
    }
  }
  # test restr.seq
  if (!is.character(inputs[[10]])) {
    stop(paste("restr.seq:", toString(inputs[[10]]),
               "must be string (character): ('string', nostring)"))
  }
  # test normali
  if (!is.logical(inputs[[11]])) {
    stop(paste("normali:", toString(inputs[[11]]),
               "must be Boolean (logical, TRUE or FALSE)"))
  }
  # test log10t
  if (!is.logical(inputs[[12]])) {
    stop(paste("log10t:", toString(inputs[[12]]),
               "must be Boolean (logical, TRUE or FALSE)"))
  }

  #  test if raw.files have 3 or 2 columns
  if (length(colnames(input.file)) == 2) {
    with.contr <- F
  } else if (length(colnames(input.file)) >= 3) {
    with.contr <- T
  } else {
    stop("input list must have 3 or 2 columns: path to fastq file, name or",
         "id and group 'c' for control or 's' for sample if wished")
  }

  # test if name begins with number
  for (i in 1:length(input.file[,2])) {
    if (!is.na(suppressWarnings(as.numeric(substring(input.file[i,2],
                                                     first = 1, last = 1 ))))) {
      stop(paste("Name cannot start with a number! change",  input.file[i,2],
                 "first character to alphabetic: e.g. t.1234"))
    }
  }

  # count numbers of "s" sample and "c" control and test if they are equal
  if (with.contr) {
    sample.ids   <- NULL
    controll.ids <- NULL
    for (i in 1:length(input.file[,3])) {
      if (input.file[i,3] == "s") {
        sample.ids <- c(sample.ids, i)
      } else if (input.file[i,3] == "c") {
        controll.ids <- c(controll.ids, i)
      } else {# if not c or s is used
        stop(paste("group column '", input.file[i,3],
                   "' must be 's'  for sample and 'c' for control", sep = ""))
      }
    }
    if (length(sample.ids) != length(controll.ids)) {
      stop(paste("number of samples 's'", length(sample.ids),
                 "must be the same as number of controls 'c'",
                 length(controll.ids), "!"))
    }
  }

  # test if all cols have the same row lengths
  col1 <- length(input.file[,1])
  col2 <- length(input.file[,2])

  if (col1 != col2) {
    stop(paste("Your input .txt file has missing values: file paths, names:",
               col1, col2))
  }

  if (with.contr) {
    col3 <- length(input.file[,3])
    if (col1 != col3) {
      stop(paste("Your input .txt file has missing values: groups:",
                 "number of paths" ,col1, "versus number of group", col3))
    }
  }

  # return index position of sample.ids and controll.ids
  if (with.contr) {
    return(list(sample.ids, controll.ids))
  } else {
    return(NULL)
  }
}


################################################################################
##                             make empty bins                                ##
################################################################################
mapp.to.bins <- function(bin.len, genom.bs,  chromo.names) {
  # make empty bins (grange)
  #
  # Args:
  #   bin.len: (int) the length of the bins
  #   genom.bs: (BSGenome) genome
  #   chromo.names: (vector of string) chromosome names or NULL
  #
  # Return:
  #   grange with chosen bin size

  message(paste(format(Sys.time(), "%H:%M:%S:"), "Make", bin.len,
                "bp bins for chromosome",
                toString(chromo.names)))

  # get the empty bin grange
  return(make.empty.genomicRange.bin3(animal = genom.bs,
                                      chrom.names = chromo.names,
                                      bin.size = bin.len, output.mesage = F))
}


################################################################################
##                            export bed files                                ##
################################################################################
export.bed.file <- function(gran.data, exp.name, s.c.id, norm.width=10,
                            frag = F, bed.index="test", write.out = T) {
  # export wig/bed file for UCSC genomeBrowser import (https://genome.ucsc.edu/)
  # this script makes the rolling mean averaging after log2 fold change
  #
  # Args:
  #   gran.data: (grange) the Grange object containing the reads information
  #   exp.name: (string) the name of the wig file
  #   s.c.id: (vector of 2 int) index of sample and control
  #   norm.width: (integer) number of sliding region (standard = 10)
  #   frag: (bool) if T a GATC fragment GRange is expects if F a GATC site
  #   bed.index: (string) the name to put in the bed file
  #   write.out: (bool) if output should be saved
  #
  # Return:
  #   the grange with the log fold change in it if wirte.out = T (saves as bed)


  # make log 2 fold change
  # test if GATC fragment or site

  if (frag) {
    grange.em <- as.matrix(elementMetadata(gran.data))
    # change NA to 0
    grange.em[is.na(grange.em)] <- 0

    # sum forward and revers read hits for each fragment
    counter <- 1
    # make emty matrix
    new.met.dat <- matrix(nrow = length(grange.em[,1]),
                          ncol = (length(grange.em[1,]) / 2))
    colnames(new.met.dat) <- colnames(grange.em)[c(seq(1, length(grange.em[1,]),
                                                       2))]

    for (i in 1:(length(grange.em[1,]) / 2)) {
      new.met.dat[,i] <- rowSums(grange.em[,c(counter, counter + 1)])
      counter <- counter + 2
    }
    elementMetadata(gran.data) <- new.met.dat
  }


  # chrI is not allowed so must be changed to chr1 instead
  lfc.matrix <- make.log.fold.change(gran.data, s.c.id, norm.tr = F)

  # add log fold change to grange
  elementMetadata(gran.data) <- lfc.matrix

  # separate, split grange for chromosomes
#  chr.gran <- split(gran.data, seqnames(gran.data))

  #smoothing averaging with 10 GATC sites rolling over all.
  #library(zoo)

  # make empty result vector
#   res.smooth.lfc <- NULL
#
#   # go over all chromosomes
#   for (chr in seqnames(gran.data)@values) {
#     # get metadata for chromosome
#     temp.metdat <- elementMetadata(chr.gran[[chr]])[,1]
#
#     # make smoothing
#     smooth.lfc <- rollapply(temp.metdat, width = norm.width, FUN = mean,
#                             align = "center", partial = T) #, fill = NA)
#
#     # add chromosome to result vector
#     res.smooth.lfc <- append(res.smooth.lfc, smooth.lfc)
#   }

  # delete all metadata
#  elementMetadata(gran.data) <- NULL
  # make 1 line of log fold change metadata
#  elementMetadata(gran.data) <-  res.smooth.lfc


  # assinge score as name for wig file
  names(elementMetadata(gran.data)) <- "score"

  bed.name    <- paste(exp.name, ".bed", sep = "")

  if (write.out) {

    # write home made bedGRaph
    track.info <- paste('track type=bedGraph name="', bed.index,
                        '" description="BedGraph format"',
                        ' visibility=full color=200,100,0 altColor=0,100,200',
                        ' priority=20', sep = "")
    write.table(track.info, file = bed.name, row.names = F, col.names = F,
                quote = F, eol = "\n")
    bed.data.frame <- data.frame("chr" = seqnames(gran.data),
                                 "start" = start(gran.data),
                                 "stop" = end(gran.data),
                                 "score" = elementMetadata(gran.data))
    write.table(bed.data.frame, file = bed.name, row.names = F, col.names = F,
                quote = F, eol = "\n", append = T)
  }
  return(gran.data)
}

################################################################################
##                     combining multiple function to one                     ##
################################################################################
write.bed.file <- function(gran.data, exp.name, bed.index="test") {
    # write bed file.
    #
    # Args:
    #   gran.data: (grange) grange with 1 metadata column
    #   exp.name: (string) the name of the wig file
    #   bed.index: (string) the name to put in the bed file
    #
    # Return:
    #   saves a .bed file

  # give right name to metadata
  names(elementMetadata(gran.data)) <- "score"

  bed.name    <- paste(exp.name, ".bed", sep = "")

    # write home made bedGRaph
    track.info <- paste('track type=bedGraph name="', 'bin-',bed.index,
                        '" description="BedGraph format"',
                        ' visibility=full color=200,100,0 altColor=0,100,200',
                        ' priority=20', sep = "")
    write.table(track.info, file = bed.name, row.names = F, col.names = F,
                quote = F, eol = "\n")

    bed.data.frame <- data.frame("chr" = seqnames(gran.data),
                                 "start" = start(gran.data),
                                 "stop" = end(gran.data),
                                 "score" = elementMetadata(gran.data))
    write.table(bed.data.frame, file = bed.name, row.names = F, col.names = F,
                quote = F, eol = "\n", append = T)

}


################################################################################
##                     combining multiple function to one                     ##
################################################################################
from.fastq.to.grange <- function(raw.files, multi.core, exp.name, adapt.seq,
                                 errors, raw.file.name, timestampm, rep.dir,
                                 genome.name, chromosom.nams,
                                 restrict.site="GATC",
                                 dir.name.ca, dir.name.bt, dir.name.gr,
                                 res.site.grages, m.hits, bam.o.fastq) {
  # make the cutadapt, bowtie step and save GATC granges from sequencing reads
  #
  # Args:
  #   raw.files: (data.f)  data.frame with 3 or 2 cols:
  #                        col 1: path to the fasta/q sequencing read files
  #                        col 2: name of the sample
  #                        col 3: (optional), 's', 'c' for control or sample
  #   multi.core: (bool)   if T bowtie uses multiple cores
  #   exp.name:   (string) name of the experiment
  #   adapt.seq:  (string) adapter sequence
  #   errors:    (integer) number of errors allowed in the adapter removal step
  #   raw.file.name: (string) the name/path to the read index file (raw.files)
  #   timestampm: (Sys.time) the current start timestamp
  #   rep.dir:    (string)   name of the directory to put in the QS reports
  #   genome.name: (string)  the name of the genome (BSgenome)
  #   chromosom.nams: (vector or strings) the names of the chromosomes used
  #   restrict.site: (string) the restriction site sequence (mostly GATC :-)
  #   dir.name.ca:    (string) the name of the folder for cut read files
  #   dir.name.bt:    (string) the name of the folder for .bam/sam files
  #   dir.name.gr:    (string) the name of the folder for GRange Rdata files
  #   res.site.grages: (list of grange) all GATC granges (all,plus,minus,frag)
  #   m.hits: (int) the number of maximal hits in bowtie mapping
  #   bam.o.fastq: (string): 'fastq/a' = fastq/a, 'b' = bam imput file
  #
  # Return:
  #   grange with metadata of all sites reads
  #   save cut.fastq
  #   save mapped .bam files


  # test if input file format is fastq/a or bam.
  # if fastq/a = rawreads, if ,bam = already mapped
  f.o.b <- substring(bam.o.fastq, 1, 1)

  if (f.o.b == "f") {

    ##############
    ## read fastq and cutadapt step

    # make empty vector to store the names and location of the cut fasta/q files
    cut.file.name <- rep(NA, length(raw.files[,1]))
    # make empty vector to store how many reads are in raw files and cut files
    raw.file.reads <- rep(NA, length(raw.files[,1]))
    cut.file.reads <- rep(NA, length(raw.files[,1]))

    # read all fasta/q files (in raw.files) and cut adapter
    for (i in 1:length(raw.files[,1])) {
      message(paste(format(Sys.time(), "%H:%M:%S:"), "Reading file:", i, "of",
                    length(raw.files[,1]),"file:",  toString(raw.files[i,2])))

      # read raw file
      temp.raw <- read.fastq(toString(raw.files[i,1]), fasta.for = bam.o.fastq)

      # store and print number of reads in raw file
      raw.file.reads[i] <- length(temp.raw)
      message(paste(format(Sys.time(), "%H:%M:%S:"), "Cut adapter for",
                    toString(raw.files[i,2]), bam.o.fastq,
                    "; Total reads:", raw.file.reads[i]))

      # remove adapter either for fasta or fastq files
      if (bam.o.fastq == "fastq") {# for fastq
        cut.file.ret <- make.cutadapt.fastq(temp.raw,
                                            out.name = toString(raw.files[i,2]),
                                            adapter = adapt.seq,
                                            error.a = errors,
                                            res.site = restrict.site,
                                            dir.name.ca = dir.name.ca)
      } else {# for fasta
        cut.file.ret <- make.cutadapt.fasta(temp.raw,
                                            out.name = toString(raw.files[i,2]),
                                            adapter = adapt.seq,
                                            error.a = errors,
                                            res.site = restrict.site,
                                            dir.name.ca = dir.name.ca)
      }

      # store file name and number of reads in cut read file
      cut.file.name[i]  <- cut.file.ret[1]
      cut.file.reads[i] <- as.integer(cut.file.ret[2])
      # print number of cut reads
      message(paste(format(Sys.time(), "%H:%M:%S:"), cut.file.reads[i],
                    "out of",
                    raw.file.reads[i], "reads had adapter and start with",
                    restrict.site,
                    round(cut.file.reads[i] * 100 / raw.file.reads[i],
                          digits = 5), "%"))
    }

    # make tab separated list with cut read list for bowtie and save it
    bow.path <- file.path(getwd(), dir.name.ca, cut.file.name)
    bowtie.list <- data.frame("FileName" = bow.path, "SampleName" = raw.files[,2])
    # make file name and path
    bowtie.name <- paste("bowtie-import", timestampm, basename(raw.file.name),
                         sep = "-")
    bowtie.path <- file.path(getwd(), dir.name.bt, bowtie.name)
    # write cut read index file for QuasR import
    write.table(bowtie.list,
                file = bowtie.path, eol = "\r\n", col.names = T,
                row.names = F, quote = F, sep = "\t")


    ##############
    #### bowtie step

    # rune rbowtie
    alig.paths <- make.bowtie(bowtie.path, multi.core = multi.core,
                              genome.name = genome.name, pr.name = exp.name,
                              dir.name.bt = dir.name.bt, m.hits = m.hits)
    # get path to bam files
    alig.paths.d.f <- data.frame(alig.paths)
    # get only path to sam/bam file
    sam.directories <- alig.paths.d.f[,1]


  } else {# if the input files arealready mapped .bam files
    sam.directories <- as.vector(raw.files[,1])
    alig.paths.d.f <- data.frame(raw.files[,1], raw.files[,2])
    # set number of original reads to NA
    raw.file.reads <- 1
    # set number of cut reads to NA
    cut.file.reads <- 1
  }

  # copy GATC site grange to leave gatc.all untouched.
  gatc.all.m     <- res.site.grages[[1]] # for all GATC sites
  gatc.plus.m    <- res.site.grages[[2]] # for plus sites
  gatc.minus.m   <- res.site.grages[[3]] # for minus sites
  gatc.frag.m    <- res.site.grages[[4]] # for fragments

  # make empty vector to store number of reads in the .bam file
  mapped.reads.in.sam.file <- rep(NA, length(alig.paths.d.f[,1]))
  # make empty vector to store number of valide reads in .bam file
  gatc.star.reads.in.sam.file <- rep(NA, length(alig.paths.d.f[,1]))


  message(paste(format(Sys.time(), "%H:%M:%S:"),
                "Start mapping to restriction sites", sep = " "))

  # loop over all .bam files and make GRange.
  for (i in 1:length(alig.paths.d.f[,1])) {
    # import .bam file
    sam.grang <- read.sam.files.to.grange(sam.directories[i])

    # count number of reads in grange
    mapped.reads.in.sam.file[i] <- length(sam.grang)

    # save each reads grange separately as GRange
    save.name <- paste(alig.paths.d.f[i,2], exp.name, timestampm, sep = "-")
    save.GRange.object(sam.grang, file.path(getwd(), dir.name.gr, save.name))

    # map reads to GATC sites
    mapped.to.gatc <- find.restriction.sites.overlaps3(sam.grang,
                                                    res.site.grages,
                                                    strand = "all",
                                                    genom = genome.name,
                                                    chrom.nam = chromosom.nams,
                                                    restsite = restrict.site)

    # sum reads up to get number of valid reads (GATC matches)
    gatc.star.reads.in.sam.file[i] <- sum(mapped.to.gatc[[1]]) # for all reads

    # add metadata to all, plus, minus strand GRange
    gatc.all.m   <- add.metadata.to.grange(mapped.to.gatc[[1]],
                                           gatc.all.m, alig.paths.d.f[i,2])
    gatc.plus.m  <- add.metadata.to.grange(mapped.to.gatc[[2]],
                                           gatc.plus.m, alig.paths.d.f[i,2])
    gatc.minus.m <- add.metadata.to.grange(mapped.to.gatc[[3]],
                                           gatc.minus.m, alig.paths.d.f[i,2])
  }

  # save grange object (reads per GATC sites)
  save.name <- paste(exp.name, timestampm, sep = "-")
  message(paste(format(Sys.time(), "%H:%M:%S:"), "Save", restrict.site,
                "mapped GRange object as:", save.name))
  save.GRange.object(gatc.all.m, file.path(getwd(), dir.name.gr, save.name))
  # save also plus and minus strand as GRange
  save.name.p <- paste(exp.name,"plus", timestampm, sep = "-")
  save.GRange.object(gatc.plus.m, file.path(getwd(), dir.name.gr, save.name.p))

  save.name.m <- paste(exp.name,"minus", timestampm, sep = "-")
  save.GRange.object(gatc.minus.m, file.path(getwd(), dir.name.gr, save.name.m))
  # save fragment GRange
  frag.gran <- make.fragmet.counts(list(gatc.plus.m, gatc.minus.m),
                                   genom.data = genome.name,
                                   resr.site = restrict.site,
                                   gatc.fragm = gatc.frag.m )
  save.name.f <- paste(exp.name,"fragment", timestampm, sep = "-")
  save.GRange.object(frag.gran, file.path(getwd(), dir.name.gr, save.name.f))

  # calculate % of reads valide reads to export qc .txt file
  pres.usfull.reads <- gatc.star.reads.in.sam.file * 100 / raw.file.reads
  rep.d.f <- data.frame("Raw.reads" = raw.file.reads,
                        "Cut.reads" = cut.file.reads,
                        "Mapped.reads" = mapped.reads.in.sam.file,
                        "GATC.reads" = gatc.star.reads.in.sam.file,
                        "Percentage" = pres.usfull.reads)
  # give row names from sample
  row.names(rep.d.f) <- raw.files[,2]
  # save qc txt file
  qc.file.name <- paste(exp.name, "-", timestampm, ".txt", sep = "")
  qc.full.path <- file.path(rep.dir, qc.file.name)

  write.table(rep.d.f, file = qc.full.path, eol = "\r\n" )

  # return GATC sites grange with all metadata reads
  return(list(gatc.all.m,  gatc.plus.m, gatc.minus.m, frag.gran))
}


################################################################################
## plot log GATC reads of each sample against all other  ##
################################################################################

plot.normalized.gatc.number.one.vs.another.log.main <- function(ful.gra,
                                                                bin.name="gatc",
                                                                timest, main.n,
                                                                dir.name.res) {
  # plot the log10 read number per GATC and bin for each sample
  # one vs. all others in a single comparison plot (X vs. Y)
  #
  # Args:
  #   ful.gra: (grange)  grange with all metadata values
  #   bin.name: (string) the name of the bin (either GATC or bin length)
  #   timest: (timestamp) the timestamp as identifier for the plot file name
  #   main.n: (string) the name of the plot
  #   dir.name.res (string) the name of the result folder
  #
  # Return:
  #   plots are saved in the dir.name.res folder

  # test if ful.gra contains only one column
  # (required more than 1 sample for cor)
  if (length(elementMetadata(ful.gra)[1,]) <= 1) {
    return(NULL)
  }

  # get name and get metadata as matrix
  met.nam <- names(elementMetadata(ful.gra))
  ful.met <- as.matrix(elementMetadata(ful.gra))

  # make filename
  f.name <- paste(main.n, "-one-by-one-comp-", bin.name, "-", timest,
                  ".pdf", sep = "")
  f.path <- file.path(dir.name.res, f.name)

  # make pdf
  pdf(file = f.path, width = 8.267, height = 11.692)

  # overall file parameters for printing (sheet parameters)
  font.size.main <- 9
  border.layout  <- c(6,6,6,6)
  margin.layout  <- c(5.1,4.1,4.1,2.1)

  # count how many plots are needed: -1 as the first sample is compared against
  # all others and not with itself
  tot.metad    <- length(ful.met[1,]) - 1
  total.num    <- (tot.metad * (tot.metad + 1)) / 2
  num.of.plots <- 12 # number of maximal plots per page
  num.of.tog   <- 3  # number of plots in each row

  # calculate the number of pages needed to plot all comparisons
  num.of.pages <- ceiling(total.num / num.of.plots)

  # indexes to be able to compare all samples with each other.
  # get first sample and compare it with all others, then get second and
  # compare it with all except the first, and so on.
  up.nu <- tot.metad # the number when to change the x axis plot
  in.x <- 1          # the index for the x axis plot
  in.y <- 2          # the index for the y axis plot
  de.step <- 1       # the step to decrease the tot.metda number

  # go over all pages and make the plots
  for (i in 1:num.of.pages) {
    par(oma = border.layout)
    par(mar = margin.layout)

    par(mfrow = c(num.of.plots / num.of.tog, num.of.tog))

    # start plot for the next for loop (if first start is 1)
    if (i == 1) {
      j.start <- 1
    } else {
      j.start <- ((i - 1) * num.of.plots) + 1
    }

    # end position (plot) for next for loop (test also if last page)
    if ((j.start + num.of.plots - 1) < total.num) {
      j.end   <- j.start +  num.of.plots - 1
    } else {
      j.end   <- total.num
    }

    # loop to make num.of.plots plots
    for (j in j.start:j.end) {

      # plot
      plot.normalized.gatc.number.one.vs.another.log(ful.met[,in.x],
                                                     ful.met[,in.y],
                                                     nams = c(met.nam[in.x],
                                                              met.nam[in.y]),
                                                     bin.name)
      # test if x and y have to be changed
      if (j >= up.nu) {
        # increase up.nu
        up.nu <- up.nu + (tot.metad - de.step)
        de.step <- de.step + 1
        in.x <- in.x + 1 # increase in.x
        in.y <- in.x + 1 # make in.y is in.x + 1
      } else {
        # only increase in.y
        in.y <- in.y + 1 # always increase
      }
    }
    # plot title for this page
    title(main = "x vs y", outer = T)
  }
  dev.off() # close report pdf
}


################################################################################
##                                plot correlation                            ##
################################################################################

plot.results <- function(gran.data, exp.name, timest, bin.name = "gatc",
                         sap.cot.ids, norma = T, logt = T, lfc, dir.name.res,
                         qs.fold.name) {
  # plot reads distribution for each sample and correlation in one pdf file
  #
  # Args:
  #   gran.data: (grange) grange with metadata
  #   exp.name: (string)  the name of the experiment
  #   timest: (timestamp) the timestamp for the plot filename
  #   bin.name: (string)  the name of the bin ('gatc' or the bin size)
  #   sap.cot.ids: (list of 2 vectors) the index of control and sample
  #   norma: (bool) if data should be normalized (value / sum(all values))
  #   logt: (bool) if the plot should be log transformed
  #   lfc: (bool) if samples have controls or only samples
  #   dir.name.res (string) the name of the folder to store the results
  #
  # Return:
  #   a pdf file: with correlation for all files, correlation for samples and
  #               control separately if controls are available and the
  #               chromosome wide distribution of the reads for each sample

  # make normalization
  if (norma) {
    gran.data <- normalize.grange.md(gran.data)
  }

  # make  log10 transformation, corplot must be not log10 transformed! (-Int)
  if (logt) {
    gran.data.p <- log10.transfom.grange.md(gran.data)
  } else {
    gran.data.p <- gran.data
  }

  # make file name and path and start plotting
  f.name <- paste(exp.name, bin.name, timest, "docs.pdf", sep = "-")
  f.path <- file.path(dir.name.res, f.name)
  pdf(file = f.path, width = 8.267, height = 11.692)

  font.size.main <- 12
  border.layout  <- c(6,6,6,6)
  margin.layout  <- c(5.1,4.1,4.1,2.1) # default c(5,4,4,2)

  # plot correlation if more than one metadata col is available.
  if (length(elementMetadata(gran.data)[1,]) > 1 ) {

    # outer margin first page (correlation)
    par(oma = border.layout)
    par(mar = margin.layout)
    # page layout
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

    # plot correlation for all samples: type = 'a'
    # sap.cot.ids[[1]] = sample.ids, sap.cot.ids[[2]] = contorl.ids
    plot.correlation(gran.data, sap.cot.ids[[1]], sap.cot.ids[[2]],type = "a",
                     timest, bin.name, exp.name, qs.fold.name)
    # test if only one sample and control if only one sample and control
    # do not make correlation between samples only and control only
    if (length(sap.cot.ids[[1]]) > 1 && lfc == T) {
      plot.correlation(gran.data, sap.cot.ids[[1]], sap.cot.ids[[2]],
                       type = "s", qs.fold.name)
      plot.correlation(gran.data, sap.cot.ids[[1]], sap.cot.ids[[2]],
                       type = "c", qs.fold.name)

      # make title for more than one sample control pair
      title(main = paste("Correlations: all, sample only, control only for",
                         bin.name, "bins"), outer = T)
    } else {
      # make title if only one sample control pair
      title(main = paste("Correlations for",
                         bin.name, "bins"), outer = T)
    }
  }

  # plot second page (GATC methylation distribution for each sample separately)
  # plot only 4 plots each page
  num.of.plots <- 4 # number of plots each page
  total.num    <- length(elementMetadata(gran.data))
  # calculate number of pages
  num.of.pages <- ceiling(total.num/num.of.plots)

  # loop over all pages
  for (i in 1:num.of.pages) {
    par(oma = border.layout)
    par(mar = margin.layout)

    par(mfrow = c(num.of.plots,1))

    # start for the next for loop (if first start is 1)
    if (i == 1) {
      j.start <- 1
    } else {
      j.start <- ((i - 1) * num.of.plots) + 1
    }

    # end position for next for loop (also test if last page)
    if ((j.start + num.of.plots - 1) < total.num) {
      j.end   <- j.start + num.of.plots - 1
    } else {
      j.end   <- total.num
    }

    # loop to make num.of.plots plots
    for (j in j.start:j.end) {
      # get name, max and minimal values
      met.nam <- names(elementMetadata(gran.data.p)[j])
      max.value <- max(elementMetadata(gran.data.p)[,j])
      min.test <- elementMetadata(gran.data.p)[,j] # remove the -inf entries
      min.test <- min.test[which(min.test != -Inf)]
      min.value <- min(min.test)

      # define the space from the maximal value to the upper border
      if (norma) {
        # if normalized
        if (logt) {  # if log transformed
          overshoot <- 0.5
          p.type <- "l"
          y.na <- "norm log10 reads"
          min.value <- min.value - 0.5
        } else {# if normalized but not log transformed
          overshoot <- 0.0001
          p.type <- "l" #"hist"
          y.na <- "norm reads"
          min.value <- 0
        }

      } else {
        # if not normalized
        if (logt) { # if log transformed
          overshoot <- 2 # upper plot margin
          p.type <- "l"
          y.na <- "log10 reads"
        } else {# if not normalized and not log transformed
          overshoot <- 50
          p.type <- "l"
          y.na <- "reads"
          min.value <- 0
        }
      }
      # plot the chromosomal distribution
      plot.all.chrom.in.one.2(gran.data.p, main.name = met.nam,
                              y.name = y.na, #c(0, 0.01)
                              y.lim = c(min.value, max.value + overshoot),
                              file.name = "asdf.png", plot.type = p.type,
                              plot.size = c(12,8), met.dat.name = met.nam,
                              d.type = "NA", font.size = font.size.main)
    }
    # plot title for this page
    title(main = paste("distribution of methylated GATC sites for", bin.name,
                       "bins"), outer = T)
  }
  dev.off() # close report pdf
}

plot.log.fold.change <- function(gran.bin, sap.cot.ids, f.path, exp.name) {
  # plot log2 fold change between all samples and all controls
  #
  # Args:
  #   gran.bin: (grange) grange with metadata
  #   sap.cot.ids: (list of two) the ids for control and samples
  #   f.path: (string) the path and name of the plot after saving
  #   exp.name: (string) name of the experiment
  #
  # Return:
  #   pdf with log2 fold change

  # make the log2 fold transformation between samples and controls
  log.fold <- make.log.fold.change(gran.bin, sap.cot.ids)

  # add log fold values to grange, after removing all old metadata
  bin.lfc.grange <- gran.bin
  elementMetadata(bin.lfc.grange ) <- NULL # clean metadata
  bin.lfc.grange  <- add.metadata.to.grange(log.fold, bin.lfc.grange,
                                            m.name = "lfc")

  # define font size and start pdf plotting
  font.size.main <- 9
  pdf(file = f.path, width = 8.267, height = 4)

  # find maximal and minimal value
  max.reads <- max(elementMetadata(bin.lfc.grange)[,1])
  min.test <- elementMetadata(bin.lfc.grange)[,1] # remove the -Inf entries
  min.test <- min.test[which(min.test != -Inf)]
  min.value <- min(min.test)
  # plot log2 fold change  y.lim=c(0,(max.reads+0.5)
  plot.all.chrom.in.one.2(bin.lfc.grange,
                          main.name = paste("log2 fold change", exp.name),
                          y.name = "log2 fold change",
                          y.lim = c((min.value  - 0.5), (max.reads + 0.5)),
                          file.name = "resolution.png", plot.type = "h",
                          plot.size = c(12,8), met.dat.name = "lfc",
                          d.type = "NA", font.size = font.size.main)
  dev.off() # close pdf

  # return log2 fold change
  return(bin.lfc.grange)

}


################################################################################
##                                main  function                              ##
################################################################################
damid.seq <- function(raw.file.name, input.format="fastq", multi.core=F,
                      exp.name="damid-gatc-sites", bin.len=100000,
                      adapt.seq="CGCGGCCGAG", errors=1, qc=T,
                      species="BSgenome.Celegans.UCSC.ce10", chr.names = NULL,
                      restr.seq="GATC", normali=T, log10t=T, m.hits=1) {
  # DamID-seq pipeline: from sequencing reads to result: main, default function
  #
  # Args:
  #   raw.file.name: (string)  path to a .txt file with 3 or 2 columns.
  #                  columns must be tab separated
  #                  first row = header: you can choose your header freely
  #                  column 1: path and name to the sequencing files
  #                  column 2: name of the sample starting with a letter
  #                  column 3: (optional) group 's'=sample and 'c'=control
  #   input.format:  (string) if 'fa' / 'f' fastq sequencing reads are cut,
  #               mapped, if 'gr' / 'g' cut and mapping is skipped.
  #               Instead of fastq sequencing reads already mapped
  #               GRanges are listed in in raw.file.name
  #               if 'bam' / 'b' .bam file must be listed in the raw.file.name
  #               list as the mapping is scipped
  #   multi.core: (bool) if T all cores of the machine is used
  #               (some machine have trouble using multicore)
  #   exp.name:   (string) freely chosen name of the experiment
  #   bin.len:    (int) length of the bins
  #   adapt.seq:  (string) the sequence of the adapter
  #   errors:     (int) how many errors allowed in the cutadapt process.
  #   qc:   (bool) if quality of sequencing reads should be tested (fastqc)
  #   species: (string) the name of the BSgenome package
  #   chr.names: (vector of strings) name of the chromosomes to be considered
  #              in the analysis. if NULL all chromosome in the genome are used
  #   restr.seq: (string) the restriction site sequence (should be GATC)
  #   normali: (bool) if the data should be normalized for the chromosomal
  #            distribution plots.
  #   log10t:  (bool) if the number of reads should be log190 transformed in the
  #            chromosomal distribution plots.
  #   m.hits: (int) the number of maximal hits in bowtie mapping
  #           (1 = only unique hits)
  #
  # Return:
  #   list of GRanges (all, plus, minus, fragment, and also binned)

  # get timestamp of pipeline start time
  timestampm <-  as.integer(as.numeric(Sys.time()))

  # read raw.files .txt list
  raw.files <- read.table(raw.file.name, header = T, sep = "\t")

  # Test if the input file is valid (get the index of sample and control ids)
  # sap.cot.ids is NULL if no log2 fold change should be made
  # also test the rest of the input data
  sap.cot.ids <- test.input.data(raw.files, list(multi.core, exp.name, bin.len,
                                                 adapt.seq, errors,
                                                 input.format, qc,
                                                 species, chr.names, restr.seq,
                                                 normali, log10t))
  # set if log2 fold change should be made at the end of the pipeline
  if (length(sap.cot.ids) == 0) {
    lfc = F
  } else {
    lfc = T
  }

  # make integer if float was given for the number inputs (errors and bin size)
  errors  <- as.integer(errors)
  bin.len <- as.integer(bin.len)

  # test, load or install full genome data BSgenome
  # This function also test the correctness of the user input (chromosome names)
  # therefore this function has to be always first before make and work
  # further in the pipeline.
  bc.genome <- get.genome(species, chr.names)

  ## make all folders to store the resulting files
  # make cutadapter folder
  dir.name.ca <- "cutadapter"
  dir.create(file.path(getwd(), dir.name.ca), showWarnings = FALSE)
  # make directory to save alignments
  dir.name.bt <- "bowtie"
  dir.create(file.path(getwd(), dir.name.bt), showWarnings = FALSE)
  # make grange folder
  dir.name.gr <- "granges"
  dir.create(file.path(getwd(), dir.name.gr), showWarnings = FALSE)
  # qs folder name
  qs.fold.name <- "qc-reports"
  dir.create(file.path(getwd(), qs.fold.name), showWarnings = FALSE)
  # make directory for results
  dir.name.res <- "results"
  dir.create(file.path(getwd(), dir.name.res), showWarnings = FALSE)


  # make or load restriction sites (GATC) granges
  res.site.grages <- make.restriction.site.grange(species = bc.genome,
                                                  rest.seq = restr.seq,
                                                  save.location = dir.name.gr,
                                                  chr.nam = chr.names)

  # get all chromosome name if chroms = NULL
  if (is.null(chr.names)) {
    chr.names <- seqnames(bc.genome)
  }

  # make quality control of sequencing files
  if (qc) {
    make.qc(raw.files, timestampm, exp.name = exp.name, rep.dir = qs.fold.name)
  }

  # cut adapter, map reads using bowtie, save grange with reads per GATC
  fir.in.chr <- substring(input.format, 1, 1) # get fist character
  if (fir.in.chr == "f" | fir.in.chr == "b") {
    gatc.all.m.l <- from.fastq.to.grange(raw.files, multi.core, exp.name,
                                         adapt.seq, errors, raw.file.name,
                                         timestampm, rep.dir = qs.fold.name,
                                         genome.name = species,
                                         chromosom.nams = chr.names,
                                         restrict.site = restr.seq,
                                         dir.name.ca = dir.name.ca,
                                         dir.name.bt = dir.name.bt,
                                         dir.name.gr = dir.name.gr,
                                         res.site.grages,
                                         m.hits = m.hits,
                                         bam.o.fastq = input.format)
  } else {# or load granges with reads grange
    # prepare path and name to save results
    save.name.g <- paste(exp.name, timestampm, sep = "-")
    file.save.name <- file.path(getwd(), dir.name.gr, save.name.g)
    # import GRange with cut reads and make GATC Grange with metadata
    gatc.all.m.l <- make.full.gatc.grange.from.read.grange2(raw.files,
                                             save.nam = file.save.name,
                                             genom = species,
                                             restsite = restr.seq,
                                             chr.na = chr.names,
                                             res.site.grages = res.site.grages)
  }

  # make single grange out of grange list
  gatc.all.m   <- gatc.all.m.l[[1]]
  gatc.plus.m  <- gatc.all.m.l[[2]]
  gatc.minus.m <- gatc.all.m.l[[3]]
  gatc.fragm   <- gatc.all.m.l[[4]]

  ## sage GATC reads as table .txt
  # save GATC sites table for non R use of the data plus and minus separated
  file.naa <- paste("reads-", exp.name, "-gatc-reads-sep",  timestampm, ".txt",
                    sep = "")
  file.paa <- file.path(dir.name.res, file.naa)
  save.table.from.grange2(gatc.all.m, gatc.plus.m, gatc.minus.m, file.paa)

  # save table with GATC sites for non R use of the data
  file.na <- paste("reads-", exp.name, "-gatc-",  timestampm, ".txt", sep = "")
  file.pa <- file.path(dir.name.res, file.na)
  save.table.from.grange(gatc.all.m, file.pa)
  # save table with GATC fragments for non R use of the data
  file.na.f <- paste("fragment-reads", exp.name, "-gatc-",  timestampm, ".txt",
                     sep = "")
  file.pa.f <- file.path(dir.name.res, file.na.f)
  save.table.from.grange(gatc.fragm, file.pa.f)

  # make empty bin grange
  bin.grange <- mapp.to.bins(bin.len, genom.bs = bc.genome,
                             chromo.names =  chr.names)

  ############
  ## make reads per bin out of reads per GATC site and save table

  # bin gatc grange to bin grange
  bin.grange <- find.all.overlaps.and.sum.metadata3(gatc.all.m, bin.grange,
                                                    met.data = T, messag = T )

  # save table with GATC sites for non R use of the data
  file.na <- paste("reads-", exp.name, "-", toString(bin.len), "-",
                   timestampm, ".txt", sep = "")
  file.pa <- file.path(dir.name.res, file.na)
  save.table.from.grange(bin.grange, file.pa)


  ###########
  ## plot correlation, distribution of reads on the chromosome and lfc

  message(paste(format(Sys.time(), "%H:%M:%S:"),
                "plotting result and documentation"))

  # plot data for GATC sites
  plot.results(gatc.all.m, exp.name, timestampm, bin.name = "gatc", sap.cot.ids,
               norma = normali, logt = log10t, lfc, dir.name.res, qs.fold.name)

  # plot correlation if more than 1 samples is available.
  plot.normalized.gatc.number.one.vs.another.log.main(gatc.all.m,
                                                      bin.name = restr.seq,
                                                      timestampm, exp.name,
                                                      dir.name.res)


  # plot data for defined bin sites
  plot.results(bin.grange, exp.name, timestampm, bin.name = toString(bin.len),
               sap.cot.ids, norma = normali, logt = log10t, lfc, dir.name.res,
               qs.fold.name)
  # plot sample GATC reads comparison when multiple samples are available.
  plot.normalized.gatc.number.one.vs.another.log.main(bin.grange,
                                                  bin.name = toString(bin.len),
                                                  timestampm, exp.name,
                                                  dir.name.res)

  ##############
  ## plot log2 fold change and make bed file.

  if (lfc) {
    # make path and name for log2 fold change plot
    f.name.r <- paste(exp.name, timestampm, "res.pdf", sep = "-")
    f.path.r <- file.path(dir.name.res, f.name.r)

    # call function to do log fold change
    message(paste(format(Sys.time(), "%H:%M:%S:"), "Plot log2 fold change"))

    # plot result
    lof2.f <- plot.log.fold.change(bin.grange, sap.cot.ids, f.path.r, exp.name)


    ## print .bed file
    file.na.b <- paste(exp.name, "-", toString(bin.len), "-",
                     timestampm, sep = "")
    file.pa.b <- file.path(dir.name.res, file.na.b)


    export.bed.file(gatc.all.m, file.pa.b, sap.cot.ids, norm.width = 10,
                   frag = F, bed.index = exp.name, write.out = T)


    file.na.b <- paste(exp.name, "-", toString(bin.len), "-", bin.len,
                       "-bp-bin-", timestampm,  sep = "")
    file.pa.b <- file.path(dir.name.res, file.na.b)

    # export bed file using bins
    write.bed.file(lof2.f, file.pa.b, bed.index = paste(exp.name, "bin"))
  }

  # return bin and GATC grange
  if (lfc) {
    return(list("GATC" = gatc.all.m, "GATC.plus" = gatc.plus.m,
                "GATC.minus" = gatc.minus.m, "GATC.frag" = gatc.fragm,
                "bins" = bin.grange, "lfc" = lof2.f))
  } else {
    return(list("GATC" = gatc.all.m, "GATC.plus" = gatc.plus.m,
                "GATC.minus" = gatc.minus.m, "GATC.frag" = gatc.fragm,
                "bins" = bin.grange))
  }
}
