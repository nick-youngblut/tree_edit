#!/usr/bin/Rscript
mod <- "10/27/12";
version <- "0.3";
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
	make_option(c("-f", "--format"), type="character", default="newick", help="Tree file format (newick or nexus). [newick]"),
	make_option(c("-n", "--name"), type="character", help="Name file; 2 column: old_names new_names   (no header!)"),
	make_option(c("-o", "--outname"), type="character", help="Output file name. [modified input file name]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output"),
	make_option(c("-z", "--Description"), action="store_false", default=TRUE, help="Script description: prune and/or rename the tips of a tree. Use 'delete' in 2nd column of rename file to drop a tip.")
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
if(is.null(opt$name)){ stop(" ERROR: provide a name file (-n)")}
if(is.null(opt$outname)){ opt$outname <- ext.edit(opt$tree, paste(c("_prn", ext), collapse="") ) }


### data processing
# loading files #
if(ext == ".nwk"){ nwk <- read.tree(opt$tree) } else
if(ext == ".tre"){ nwk <- read.nexus(opt$tree) }
name.tbl <- read.table(opt$name)
# sanity check #
#if(length(nwk$tip) != nrow(name.tbl)){
#	stop("Number of tips and number of names does not match!")
#}
name.tbl.ndel <- name.tbl[name.tbl[,2]!="delete", ]			
if(length(name.tbl.ndel[,2]) != length(unique(name.tbl.ndel[,2]))){
	print(table(as.vector(name.tbl.ndel[,2])))
	stop("The new name column of the rename file contains duplicate names!")
}

# order name table to phylogeny #
for(i in 1:nrow(name.tbl)){
	if(as.vector(name.tbl[i,2]) == "delete"){
		nwk <- drop.tip(nwk, as.vector(name.tbl[i,1]))
		
	} else {
		q <- paste(c("^", as.vector(name.tbl[i,1]), "$"), collapse="")
		nwk$tip.label[grep(q, nwk$tip)] <- as.vector(name.tbl[i,2])
	}
}


# writing tree #
if(ext == ".nwk"){ write.tree(nwk, file=opt$outname) } else
if(ext == ".tre"){ write.nexus(nwk, file=opt$outname) }
