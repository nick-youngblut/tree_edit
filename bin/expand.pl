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
GetOptions(
	   "tree=s" => \$tree_in,		# tree file
	   "count=s" => \$count_in,		# habitat file (2column)
	   "format=s" => \$format,		# tree format
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);
die " ERROR: Provide a count file (Mothur format)." if ! $count_in;

### I/O error & defaults
(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_exp.nwk/;

### MAIN
my $treeo = tree_io($tree_in, $format);
my $count_ref = load_count($count_in);
$treeo = expand_taxa($treeo, $count_ref);
tree_write($treeo, $outfile);



### Subroutines
sub tree_write{
  ### writting out a newick tree file ###
  my ($treeo, $outfile) = @_;
  my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
  $out->write_tree($treeo);
  print STDERR " Newick tree file written:\n  $outfile\n";
}

sub write_root{
  my $treeo = shift;
  my $root = $treeo->get_root_node;
  print "Root_id = ", $root->id, "\n";
}

sub expand_taxa{
  # expanding (replicating) taxa so that there is a rep (clonal) for each sample #
  my ($treeo, $count_ref) = @_;
  
  my $leaf_cnt = 0;
  for my $leaf ($treeo->get_leaf_nodes){
    #print Dumper $leaf->id; 
    $leaf_cnt++;
    
    # finding ancestor #
    my $anc = $leaf->ancestor;
    
    # getting leaf info #
    my $leaf_id = $leaf->id;
    my $leaf_br = $leaf->branch_length;
    
    # adding clonal replicate taxa (for each sample) #
    die " ERROR: ", $leaf_id, " not found in count table\n" if ! exists $$count_ref{$leaf_id};
    foreach my $samp (@{$$count_ref{$leaf_id}}){
      # adding clonal taxa #
      $leaf->add_Descendent($leaf->new(-id => join("_", $samp, $leaf_cnt), 
				       -branch_length => 0));
      
      # renaming leaf #
      $leaf->id("");
    }
    
  }
  
  return $treeo;
}

sub load_count{
  # loading count table (mothur format #
  my $count_in = shift;
  
  open IN, $count_in or die $!;
  my %cnt_tbl;
  my @header;
  while(<IN>){
    chomp;
    if($.==1){	# header
      @header = split /\t/;
      @header = map{s/"//g; $_} @header;
    }
    else{
      my @line = split /\t/;
      @line = map{s/"//g; $_} @line;
      
      for my $i (2..$#line){
	push(@{$cnt_tbl{$line[0]}}, $header[$i]) if $line[$i] > 0;
      }
    }
  }
  close IN;
  #print Dumper %cnt_tbl; exit;
  return \%cnt_tbl;
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

tree_expand.pl -- Replicate (add clonal taxa) to a tree.

=head1 SYNOPSIS

tree_expand.pl -t -c [-f]

=head2 options

=over

=item -t 	Tree file (newick | nexus)

=item -c 	Count file (Mothur format)

=item -f 	Tree file format. [newick]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_expand.pl


=head1 DESCRIPTION

The script will duplicate taxa (clonal duplicates), in order to 
have a representative from each environment in the tree. Essentially,
tree tips are added, each without branch length to the original 
tip used for making the added tips. This is needed for AdaptML.

=head2 Example:

OTU1 is found in Sample1 & Sample2 (as defined in the cound file)

OTU1 would then be split in the tree to Sample1_1 Sample2_1

=head2 Mothur count file format:

http://www.mothur.org/wiki/Count_File


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

