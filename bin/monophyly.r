#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: monophyly.r [options] <tree> <taxonomy>

Options:
  <tree>      Newick tree file.
  <taxonomy>  A table of taxonomic classifications for each tree tip.
              The first column must match tree tip labels.
              Note: no table header.
  -r=<r>      Name of taxon to root tree on (if needed).
  -o=<o>      Output file name. [Default: monophyly_results]
  -h          Help

Description:
  Assess monophyly of clades in a tree.

  This is just a wrapper around the MonoPhy R package.
  Running MonoPhy in a Jupyter notebook can cause dead kernels,
  so this allows running outside of an R kernel.

  The taxonomy file can have extra taxa in it; these will be
  filtered out based on tree tip labels. Also, the root
  taxon will be added to the table if it is missing.

  OUTPUT:
    The R object produced by AssessMonophyly() will be saved
    for latter loading into an R/Jupyter session.

' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'MonoPhy', 'ape')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}


# functions
add_root_to_taxonomy = function(df.tax, root.id){
  n.col = ncol(df.tax)
  root = as.data.frame(c(root.id, 'Archaea', rep('unclassified', n.col - 2)))
  colnames(root) = 'root'
  root = t(root)
  colnames(root) = colnames(df.tax)
  df.tax = rbind(df.tax, root)
  return(df.tax)
}


# main
## loading tree
tree = read.tree(opts[['<tree>']])

## loading taxonomy
df.tax = read.delim(opts[['<taxonomy>']], sep='\t', header=FALSE) 

# filtering metadata to just OTUs in tree
df.tax = df.tax %>%
  filter(V1 %in% tree$tip)

# adding root to taxonomy if needed
root.id = opts[['-r']]
if (! is.null(root.id)){
  nroot = df.tax %>% filter(V1 == root.id) %>% nrow
  if (nroot == 0){
    df.tax = add_root_to_taxonomy(df.tax, root.id)
  }
}

## formatting tree
if(! is.null(root.id)){
  tree = unroot(tree)
  tree = root(tree, root.id)
}
if(! is.rooted(tree)){
  tree = unroot(tree)
  tree$root.edge = 0
}
tree = multi2di(tree)


# MonoPhy
monophy.res = AssessMonophyly(tree, df.tax)
outFile = file.path(workDir, opts[['-o']])
saveRDS(monophy.res, outFile)
cat('File saved to:', outFile, '\n')
