#!/usr/bin/Rscript
mod <- "";
version <- "";
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
	make_option(c("-s", "--species"), action="store_true", default=FALSE, help="Species tree? [FALSE]"),
	make_option(c("-m", "--multi"), action="store_false", default=TRUE, help="Species tree? [TRUE]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output"),
	make_option(c("-x", "--Description"), action="store_false", default=TRUE, help="Script Description: 'cleans' trees for HiDe analysis (score.lua). No branch lengths, bifurcating nodes, no label punctuation")
	)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


### functions 
HiDe.clean.tree <- function(tree, in.file, species=FALSE){
	# in.file <- "RAxML_bootstrap.gtlenv_OTU2_LCB1"
# cleaning gene and species trees for HiDe #
	require(ape)
	
	## bifurcating, & leaf name punctuation #
	# all nodes to bifurcating #
	if(is.binary.tree(tree) == FALSE){
		tree <- multi2di(tree)
		}
	
	# rooting #
	tree <- root(tree, tree$tip.label[1], resolve.root=TRUE)

	# making output name #
	if(species==TRUE){
		out.file <- "species.newick"
		} else {
		out.file <- gsub("\\.nwk$|\\.newick$|$", ".newick", in.file, perl=TRUE)
			#out.file <- gsub("RAxML_bootstrap\\.", "", out.file, perl=TRUE)
		}
	
	message(out.file)
	write.tree(tree, out.file, digits=0)
	}

HiDe.clean.tree.multi <- function(tree, in.file, species=FALSE){
	# in.file <- "RAxML_bootstrap.gtlenv_OTU2_LCB1"
# cleaning gene and species trees for HiDe #
	require(ape)
	
	# bifurcating, & leaf name punctuation #
	tree <- lapply(tree, function(x){

		# all nodes to bifurcating #
		if(is.binary.tree(x) == FALSE){
			x <- multi2di(x)
			}
		
		return(x)
		})
	
	# rooting #
	tree <- lapply(tree, function(x){
		root(x, x$tip.label[1], resolve.root=TRUE)
		})
	
	# making multi-phylo object	
	tree.m <- do.call(c.phylo, tree)
	
	# making output name #
	if(species==TRUE){
		out.file <- "species.newick"
		} else {
		out.file <- gsub("\\.nwk$|\\.newick$|$", ".newick", in.file, perl=TRUE)
			#out.file <- gsub("RAxML_bootstrap\\.", "", out.file, perl=TRUE)
		}
	
	message(out.file)
	write.tree(tree.m, out.file, digits=0)
	}
	

### I/O error check
if(grepl("^new", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".nwk"} else
if(grepl("^nex", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".tre"} else{
	stop(" ERROR: tree format must be nexus or newick")
	}
if(is.null(opt$tree)){ stop(" ERROR: provide a tree file (-t)")}

### data processing
if(ext == ".nwk"){ tree <- read.tree(opt$tree) } else
if(ext == ".tre"){ tree <- read.nexus(opt$tree) }


if(opt$multi ==TRUE){ HiDe.clean.tree.multi(tree, opt$tree, opt$species) } else {
	HiDe.clean.tree(tree, opt$tree, opt$species) }


