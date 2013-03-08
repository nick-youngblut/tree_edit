#!/usr/bin/env perl

use Bio::TreeIO;
use Bio::Tree::RandomFactory;
 
# initialize a TreeIO writer to output the trees as we create them
my $out = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">randomtrees.tre");
my @listoftaxa = qw(A B C D E F G H);
my $factory = new Bio::Tree::RandomFactory(-taxa => \@listoftaxa);
 
# generate 10 random trees
for( my $i = 0; $i < 10; $i++ ) {
    	print STDERR " tree $i written\n";
	 $out->write_tree($factory->next_tree);
}
# One can also just request a total number of taxa (8 here) and
# not provide labels for them
# In addition one can specify the total number of trees
# the object should return so we can call this in a while
# loop
#$factory = new Bio::Tree::RandomFactory(-num_taxa => 8
 #                                       -max_count=> 10);
#while( my $tree = $factory->next_tree) {
#  $out->write_tree($tree);
#}

