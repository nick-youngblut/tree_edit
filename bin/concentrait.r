#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'Usage: concentrait.r [Options] [-x | <tree> <trait>]

Options:
  <tree>       Newick tree file (multitree)
  <trait>      Trait table (no headers)
  -p=<p>       Percent shared trait cutoff.
               [Default: 90]
  -c=<c>       Cluster size file name prefix.
               [Default: cluster_size]
  -d=<d>       Cluster distance file name prefix.
               [Default: cluster_dist]
  -t=<m>       Tau_D table name.
               [Default: Tau_D.txt]
  -x=<x>       Create a test tree & trait file with `x` taxa.
  -b=<b>       Number of simulated trees to make.
               [Default: 10]
  -h           Help

Description:
  The script requires two input files as arguments:
  a Newick Tree and tab delimited text file with names
  of each taxon in the first column and then 0 or 1
  values for each trait in the following columns
  (and no headers).

  Tau_D is average consensus sequence distance (branch length)
  between trait values.

  Dependencies: data.table, adephylo, ape, docopt

  For more info: http://www.ess.uci.edu/group/amartiny/research/consentrait
' -> doc
opts = docopt(doc)


# packages
pkgs <- c('data.table', 'adephylo', 'ape')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# Test files
if(!is.null(opts[['-x']])){
  out.tree = 'concentrait_TEST.nwk'
  out.trait = 'concentrait_TEST.txt'
  # tree(s)
  n.taxa = as.numeric(opts['-x'])
  tree = rmtree(opts[['-b']], n.taxa)
  write.tree(tree, out.tree)
  
  ## traits
  taxa = sapply(1:n.taxa, function(x) paste0('t', x))
  df = data.frame('taxa' = taxa,
    'trait1' = sample(c(0,1), n.taxa, replace=TRUE),
    'trait2' = sample(c(0,1), n.taxa, replace=TRUE))
  write.table(df, out.trait, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
  msg = paste(c('Test files written: ', out.tree, ', ', out.trait, '\n'), collapse='')
  cat(msg)
  
  opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()
}


# Params
perc.share.cutoff = as.numeric(opts[['-p']])
stopifnot((perc.share.cutoff >= 0) &  (perc.share.cutoff <= 100))
perc.share.cutoff = perc.share.cutoff / 100


# Import
## Newick tree (multitree) - replace to read.nexus if using nexus tress
tree_all = read.tree(opts[['<tree>']],keep.multi = TRUE)
## Trait table w. no headers
table = read.table(opts[['<trait>']], sep = "\t", header=FALSE)

# Init
n.trees = length(tree_all)
Mean_all = matrix(nrow=ncol(table)-1,ncol=n.trees)

#loop through all trees
for (m in 1:length(tree_all)) {
  cat("Analyzing tree: ",m,"\n")
  
  # testing if table and tree contain the same entries - else drop tips
  tree = tree_all[[m]]  
  z = subset(tree$tip.label,!(tree$tip.label %in% table[,1]))
  if (length(z) > 0) {
    drop.tip(tree,z)
  }

  #rooting tree with first taxon - change if different root
  root_tree = root(tree,1,resolve.root=T)
  #replacing negative branch lengths - e.g., from PHYLIP
  root_tree$edge.length[root_tree$edge.length <= 0] =  0.00001
  subtree = subtrees(root_tree, wait=FALSE)

  cluster_mean = numeric(length=0)
  # loop through all traits
  for (j in 2:ncol(table)) {
     cat("  Analyzing trait: ",j-1,"\n")
     #Loading trait table
     table_tmp = table[,c(1,j)]
     colnames(table_tmp)[1] = "ID";
     colnames(table_tmp)[2] = "Trait";
     
     # removing all entries not in tree 
     table_tmp2 = data.table(table_tmp)
     setkey(table_tmp2,ID)
     table2 = table_tmp2[intersect(table_tmp2$ID,root_tree$tip.label)]
     setkey(table2,ID)

     #initializing result vectors and file names
     positives = vector(mode="list",length=0)
     cluster_size = numeric(length=0)
     cluster_size_file = paste(opts[['-c']],'_t',j-1,".txt",sep="")

     cluster_dist = numeric(length=0)
     cluster_dist_file = paste(opts[['-d']],'_t',j-1,".txt",sep="")

     # Init cluster size & distance files
     if (m == 1) {
         cat(c("trait","tree","subtree", "distance","cluster_size"), '\n', file = cluster_size_file, 
            sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
         cat(c("trait","tree","subtree", "cluster","distance"), '\n', file = cluster_dist_file, 
            sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
     }


     #loop through all subtrees and determining if any subtrees have >P% positives
     for (i in 1:length(subtree)){
       tip_names = subtree[[i]]$tip.label
       #change the value below if you want a new threshold
       if (mean(table2[tip_names][,Trait]) > perc.share.cutoff ) {
        match_test = match(tip_names,positives)
        if (all(is.na(match_test))) {
            positives = c(positives,tip_names)
            cluster_dist = distRoot(subtree[[i]],tip_names, method=c("p"))
            cluster_size = append(cluster_size,mean(cluster_dist))

            # printing to files###
            cat(j-1,m,i,mean(cluster_dist),length(cluster_dist), '\n', file = cluster_size_file,
                sep = "\t", fill = FALSE, labels = NULL,append = TRUE)

            for(cdl in 1:length(cluster_dist)){
              cat(j-1,m,i,cdl,cluster_dist[cdl], '\n', file = cluster_dist_file,
                  sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
            }
        }
        else if (any(is.na(match_test))) {
            print("Assertion error: NAs present")
        }
        else {
          
        }
      }
    }

    ##### find singletons ######
    a = table2[table2$Trait == 1,][,ID]
    g = as.character(a)

    singletons_names = setdiff(g,positives)
    if (length(singletons_names) > 0) {
       for (h in 1:length(singletons_names)){
           # weigh singletons with half
           we = which.edge(root_tree,singletons_names[h])
           singleton_edges = 0.5*root_tree$edge.length[we] 
           cluster_size = append(cluster_size,singleton_edges)

           cat(j-1,m,NA,singleton_edges,1, '\n', file = cluster_size_file, sep = "\t",
             fill = FALSE, labels = NULL,append = TRUE)
       }    
    }
  # means of cluster sizes
  Mean_all[j-1,m] = mean(cluster_size)
  }
}

#output file
## formatting
Mean_all = as.data.frame(t(Mean_all))
col.n = sapply(1:(ncol(table)-1), function(x) paste(c('t', x), collapse=''))
colnames(Mean_all) = col.n
Mean_all$tree = 1:nrow(Mean_all)
Mean_all = Mean_all[,c('tree', col.n)]
## writting
write.table(Mean_all,opts[['-t']], sep = "\t", quote=FALSE, row.names=FALSE)
