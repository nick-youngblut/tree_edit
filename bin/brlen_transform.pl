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

my ($verbose, $tree_in, $find_max, $cutoff, $tips, $regex, $brlen_list);
my $format = "newick";
my $transform = "0.1";
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "transform=s" => \$transform,    # tranforming factor
	   "cutoff=s" => \$cutoff,	    # min cutoff for brlen modification
	   "brlens" => \$brlen_list, 
	   "max" => \$find_max,		    # just find max brlen?
	   "tips" => \$tips,		    # just tips? [TRUE]
	   "regex=s" => \$regex,	    # regex for ID of tip of interest
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die "Provide tree file (newick or nexus format)." if ! $tree_in;
$format = check_format($format);
(my $outfile = $tree_in) =~ s/\.[^.]+$|$/_tran.nwk/;
die " ERROR: cutoff must be numeric\n" if $cutoff && $cutoff !~ /^\d+/;
$regex = qr/$regex/ if $regex;

### MAIN
my $treeo = tree_io($tree_in, $format);
list_brlens($treeo, $tips) if $brlen_list;
find_max($treeo, $tips) if $find_max;

$treeo = transform_brlens($treeo, $tips, $transform, $cutoff, $regex);
tree_write($treeo, $outfile);


### Subroutines
sub tree_write{
  my ($tree, $outfile) = @_;
  my $out = new Bio::TreeIO(-file => "> $outfile",
			    -format => "newick");
  $out->write_tree($tree);
}


sub transform_brlens{
  # transforming branch lengths #
  my ($treeo, $tips, $transform, $cutoff, $regex) = @_;
  
  my $trans_cnt = 0;
  for my $node ($treeo->get_nodes){
    # screening nodes #
    next if ! $tips && ! $node->is_Leaf;
    next if $node->branch_length < $cutoff; # skipping if brlen is too small 
    next if $regex && $node->id !~ $regex;
    
    # transforming #
    my $new_brlen;
    if($transform =~ /^[\d.]+$/){
      $new_brlen = $node->branch_length * $transform;
      $trans_cnt++;
    }
    else{
      die " ERROR: transform must be numeric";
    }	
    
    $node->branch_length( $new_brlen )
  }	
  
  print STDERR " Number of branches transformed: $trans_cnt\n";
  
  return $treeo;
}


sub list_brlens{
  # list branch length values
  my ($treeo, $tips) = @_;
  
  for my $node ($treeo->get_nodes){
    next if ! $tips && ! $node->is_Leaf;
    
    my $brlen = $node->branch_length;
    if (defined $brlen and length $brlen > 0){
      print $brlen, "\n";
    }
  }  
  exit;
}


sub find_max{
  # find max brlen; then exiting #
  my ($treeo, $tips) = @_;
  
  my $max_brlen = 0;
  for my $node ($treeo->get_nodes){
    next if ! $tips && ! $node->is_Leaf;
    
    my $brlen = $node->branch_length;
    next if ! $brlen;
    $max_brlen = $brlen if $brlen > $max_brlen;
  }
  
  print "Maximum branch length: $max_brlen\n";
  exit;
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


__END__

=pod

=head1 NAME

brlen_transform.pl -- Alter branch lengths

=head1 SYNOPSIS

brlen_transform.pl [options]

=head2 Required Arguments

=over

=item -tree

Tree file (newick or nexus)

=back

=head2 Options Arguments

=over

=item -format

Tree format (newick or nexus). [newick]

=item -tips

Just select branch lengths of tips? [TRUE]

=item -transform

Multiplication factor for branch lengths of interest. [0.1]

=item -cutoff

Minimum branch length cutoff for transforming branch length. [0]

=item -regex

Regular expression for selecting tips (based on tip ID) for tranforming
branch lengths.

=item -brlens

Just list all branch lengths (one value per line)? [FALSE]

=item -max

Just write out the max branch length (just tips by default)? [FALSE]

=item -h	

This help message

=back

=head2 For more information:

perldoc tree_brlen_transform.pl

=head1 DESCRIPTION

Reduce or increase the branch lengths for certain taxa if the branch lengths
are too long or short. This can occur when using RAxML's EPA function.

=head1 EXAMPLES

=head2 Writing out max branch length

brlen_transform.pl -tree tree.nwk -m

=head2 Transforming branch lengths of just certain tips (regex)

brlen_transform.pl -tree tree.nwk -r "Methano"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

