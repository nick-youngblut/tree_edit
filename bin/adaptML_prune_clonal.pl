#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $count_in, $format);
my $clonal_cut = 0;
GetOptions(
	   "tree=s" => \$tree_in,	       # tree file
	   "format=s" => \$format,	       # tree format
	   "cutoff=f" => \$clonal_cut,         # branch length cutoff for determining clonal
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage          # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);


### I/O error & defaults
(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_clps.nwk/;

### MAIN
my $treeo = tree_io($tree_in, $format);
$treeo = collapse_clonal($treeo, $clonal_cut);
tree_write($treeo, $outfile);



### Subroutines
sub tree_write{
  ### writting out a newick tree file ###
  my ($treeo, $outfile) = @_;
  my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
  $out->write_tree($treeo);
  print STDERR " Newick tree file written: $outfile\n";
}

sub write_root{
  my $treeo = shift;
  my $root = $treeo->get_root_node;
  print "Root_id = ", $root->id, "\n";
}

sub collapse_clonal{
  # expanding (replicating) taxa so that there is a rep (clonal) for each sample #
  my $treeo = shift or die $!;
  my $clonal_cut = shift or 0;
  
  my $leaf_cnt = 0;
  for my $leaf ($treeo->get_leaf_nodes){
    $leaf_cnt++;

    # filtering clonal
    if ($leaf->branch_length <= $clonal_cut){
      # getting all siblings and trimming redundant EcologyIDs
      my $anc = $leaf->ancestor;
      next unless defined $anc;

      my %ecoIDs; 
      for my $child ($anc->get_all_Descendents){
	my ($ecoID, $ID) = split /_/, $child->id, 2;
	if (exists $ecoIDs{$ecoID}){
	  # pruning taxa
	  print STDERR " Pruning: " . $child->id . "\n";
	  $treeo->remove_Node($child);
	}
	else{
	  $ecoIDs{$ecoID} = 1;
	}
      }
    }    
  }
  
  return $treeo;
}


sub tree_io{
  # loading tree object #
  my ($tree_in, $format) = @_;
  my $input = Bio::TreeIO -> new(-file => $tree_in,
				 -format => $format);
  my $treeio = $input->next_tree;	
  return $treeio;
}

sub check_format{
  # checking format input #
  my $format = shift;
  if($format =~ /^new/i){ $format = "newick"; }
  elsif($format =~ /^nex/i){ $format = "nexus"; }
  else{ die " ERROR: '$format' format not recognized (newick | nexus)\n"; }
  return $format;
}

sub print_node_ids{
  # error checking
  my ($tree, $leaf, $int) = @_;
  if($leaf =~ /true/i){
    foreach my $node ($tree->get_leaf_nodes){
      print "leaf: ", $node->id, "\n";
    }
  }
  if($int =~ /true/i){
    foreach my $node ($tree->get_nodes){
      if($node->is_Leaf == 0){
	print "internal: ", $node->id, "\n";
      }
    }
  }
}
  

__END__

=pod

=head1 NAME

adaptML_prune_clonal.pl -- Prune clonal taxa with the same EcologyID

=head1 SYNOPSIS

adaptML_prune_clonal.pl -t -c [-f]

=head2 options

=over

=item -t 	Tree file (newick | nexus)

=item -f 	Tree file format. [newick]

=item -c        Branch length cutoff for defining clonal taxa. [0]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc adaptML_prune_clonal.pl

=head1 DESCRIPTION

Remove clonal taxa from a tree with the same EcologyID.
This prevents spurious results when running AdaptML.

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

