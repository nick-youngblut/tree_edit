#!/usr/bin/Rscript
mod <- "8/2/12 3:36 PM";
version <- "0.1";
author <- "Nick Youngblut";
#--------------------- version log ---------------------#
#
#
#-------------------------------------------------------#

### start-up
rm(list=ls())

### packages
suppressPackageStartupMessages(library(Hmisc))
pkgs <- Cs(
	optparse,
	ape
	)
for(i in 1:length(pkgs)){
	tmp <- pkgs[i]
	suppressPackageStartupMessages(library(pkgs[i], character.only=TRUE))
	}

### I/O
# initialize options
option_list <- list(
	make_option(c("-t", "--tree"), type="character", help="Tree file"),
	make_option(c("-f", "--format"), type="character", default="newick", help="Tree file format (newick or nexus)"),
	make_option(c("-o", "--outname"), type="character", help="Output file name. [default: modified input file name]"),
	make_option(c("-m", "--multi"), action="store_true", default=FALSE, help="Multiple trees in file? [FALSE]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output")
	)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


### functions 
ext.edit <- function(file, ext){
	file <- gsub("\\.[^\\.]+$|$", ext, file, perl=TRUE)
	return(file)
	}

### I/O error check
if(grepl("^new", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".nwk"} else
if(grepl("^nex", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".tre"} else{
	stop(" ERROR: tree format must be nexus or newick")
	}
if(is.null(opt$tree)){ stop(" ERROR: provide a tree file (-t)")}
if(is.null(opt$outname)){ opt$outname <- ext.edit(opt$tree, paste(c("_lad", ext), collapse="")) }


### data processing
# reading #
if(ext == ".nwk"){ tree <- read.tree(opt$tree) } else
if(ext == ".tre"){ tree <- read.nexus(opt$tree) }

# ladderizing
if(opt$multi==TRUE){ 
	tree <- lapply(tree, ladderize)
	tree <- do.call(c.phylo, tree) } else
	{ tree <- ladderize(tree) }

# writing tree #
if(ext == ".nwk"){ write.tree(tree, file=opt$outname) } else
if(ext == ".tre"){ write.nexus(tree, file=opt$outname) }


	
