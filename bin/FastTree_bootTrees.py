#!/usr/bin/env python

"""
FastTree_bootTrees.py: make individual bootstrap trees with FastTree

Usage:
  FastTree_bootTrees.py [options] <alignment>
  FastTree_bootTrees.py -h | --help
  FastTree_bootTrees.py --version

Options:
  <alignment>   Alignment in fasta format.
  -b=<b>        Number of bootstrap replicates.
                [Default: 100]
  -p=<p>        Number of processors.
                [Default: 1]
  --version     Show version.
  --debug       Turn of parallel processing.
  -h --help     Show this screen.

Description:
  FastTree doesn't write out each bootstrap tree and instead
  just writes bootstrap values on the final tree.
  
  This script uses cogent to load/parse the alignment, 
  bootstrap the with alignment with replacement, and 
  call FastTree with no bootstrap replication. 

  The bootstrap trees will be written to stdout.
"""

from docopt import docopt
import sys,os
import tempfile
import shutil
from functools import partial
import multiprocessing as mp

import numpy as np
from cogent.app.fasttree import build_tree_from_alignment
from cogent import DNA, LoadSeqs    


def aln2array(aln):
    seqNames = [x.lstrip('>') for x in aln.toFasta().split('\n') 
                if x.startswith('>')]
    seqs = np.array([list(x.rstrip()) for x in aln.toFasta().split('\n') 
                     if not x.startswith('>')])
    return(seqNames, seqs)


def bootTree(i, aln_len, seqNames, seqs, outdir):
    # stats
    sys.stderr.write('Inferring bootstrap tree: {}\n'.format(i))

    # bootstrap index
    boot_idx = np.random.choice(range(aln_len), aln_len, replace=True)

    # conversion back to aln object
    aln_boot = seqs[:,boot_idx]
    aln_boot = {name:''.join(seq) for name,seq in zip(seqNames, aln_boot)}
    aln_boot = LoadSeqs(data=aln_boot, moltype=DNA)

    # inferring tree
    params = {'-boot' : 0}
    tree = build_tree_from_alignment(aln_boot, moltype=DNA, params=params)

    # writing out tree 
    outFile = os.path.join(outdir, 'boot{}.nwk'.format(i))
    tree.writeToFile(outFile)

    return(outFile)



if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')

    # loading sequences
    aln = LoadSeqs(uargs['<alignment>'], moltype=DNA)

    # converting sequences to numpy array
    seqNames,seqs = aln2array(aln)
    
    # making output directory
    tmpDir = tempfile.mkdtemp()

    # bootstrapping
    boots = [x + 1 for x in range(int(uargs['-b']))]
    
    if uargs['--debug']:
        treeFiles = [bootTree(i, len(aln), seqNames, seqs, tmpDir) 
                     for i in boots]
    else:
        p = mp.Pool(int(uargs['-p']))
        bootTree_p = partial(bootTree, 
                             aln_len=len(aln), seqNames=seqNames, 
                             seqs=seqs, outdir=tmpDir)
        treeFiles = p.map(bootTree_p, boots)
    

    # combining all bootstrap trees
    for f in treeFiles:
        with open(f, 'rb') as inFH:
            for x in inFH:
                print x
    shutil.rmtree(tmpDir)
