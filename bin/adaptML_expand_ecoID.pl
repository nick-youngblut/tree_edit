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
my $split_str = "";
GetOptions(
	   "tree=s" => \$tree_in,	       # tree file
	   "format=s" => \$format,	       # tree format
	   "delim=s" => \$split_str,         #split string for ind. ecoIDs
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage          # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);


### I/O error & defaults
(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_exp.nwk/;

### MAIN
my $treeo = tree_io($tree_in, $format);
$treeo = expand_ecoID($treeo, $split_str);
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

sub expand_ecoID{
  # expanding (replicating) taxa so that there is a rep (clonal) for each sample #
  my $treeo = shift or die $!;
  my $split_str= shift;
  
  my $leaf_cnt = 0;
  for my $leaf ($treeo->get_leaf_nodes){
    $leaf_cnt++;

    my ($ecoID, $ID) = split /_/, $leaf->id, 2;
    if (! defined $ID){
      exit("ERROR: leaf name not formatted correctly: " . $leaf->ID . "\n")
    }

    my @ecoIDs = split /$split_str/, $ecoID;
    
    my @new_tips;
    if (scalar @ecoIDs > 1){
      my $anc = $leaf->ancestor;
      my $leaf_cnt = 0;
      foreach my $ind_ID (@ecoIDs){
	$leaf_cnt++;
	# adding taxon to tree
	my $new_tip = join("_", $ind_ID, join("x", $ID, $leaf_cnt));
	push @new_tips, $new_tip;
	$anc->add_Descendent($leaf->new(-id => $new_tip,
					-branch_length => 0));
      }
      # remove old node
      $treeo->remove_Node($leaf);
      # status
      print STDERR "Expanded " . $leaf->id . ":\n  " . join("\n  ", @new_tips) . "\n";
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

adaptML_expand_ecoID.pl -- Expand tips with multiple EcologyIDs into mutliple
clonal taxa

=head1 SYNOPSIS

adaptML_expand_ecoID.pl -t

=head2 options

=over

=item -t 	Tree file (newick | nexus)

=item -f 	Tree file format. [newick]

=item -d        Delimiter separating ecology IDs. [""}

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc adaptML_expand_ecoID.pl

=head1 DESCRIPTION



=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

