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

my ($verbose, $species_in, $gene_in, $tformat);
my ($sregex, $gregex);
my $nclade = 4;
my $brlen;
GetOptions(
	   "species=s" => \$species_in,			# tree file
	   "gene=s" => \$gene_in,
	   "format=s" => \$tformat,				# tree format
	   "sregex=s" => \$sregex, 				# altering species tree names
	   "gregex=s" => \$gregex, 				# altering gene tree names
	   "clade=i" => \$nclade, 				# number of clades
	   "brlen=f" => \$brlen, 				# brlen jump for collapsing clades
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " Provide a species tree file (newick or nexus).\n" if ! $species_in;
#die " Provide a gene tree file (newick or nexus).\n" if ! $gene_in;
$tformat = check_tree_format($tformat);
$sregex = qr/$sregex/ if $sregex;
$gregex = qr/$gregex/ if $gregex;

### MAIN
# loading species tree #
my $streeo = tree_io($species_in, $tformat);

# collapsing tips #
## adding brlen if not present ##
#add_brlen($streeo);

## getting max branch length
#my $max_brlen = max_brlen($streeo);
my $max_brlen = max_height($streeo);
$brlen = $max_brlen / 50 unless $brlen;

## collapsing ##
$streeo = collapse_tips_brlen($species_in, $nclade, $brlen, $max_brlen);

#print Dumper $streeo; exit;

# writing modified tree #
tree_write($streeo);


### Subroutines
sub tree_write{
	### writting out a nexus tree file ###
	my ($treeo) = @_;
	
	my $out = new Bio::TreeIO(-fh => \*STDOUT, -format => "newick");
	$out->write_tree($treeo);
	}
	
sub add_brlen{
	my ($treeo) = @_;
	
	foreach my $node ($treeo->get_nodes){
		$node->branch_length($brlen) unless $node->branch_length();
		}
	}

sub collapse_tips_brlen{
# collapsing #
	my ($species_in, $nclade, $brlen, $max_brlen) = @_;

	my $nclade_brlen; 
	for (my $i=0; $i<=$max_brlen; $i+=$brlen){
		my $treeo_tmp = tree_io($species_in, $tformat);		# loading tree each time
		
		for my $node ($treeo_tmp->get_nodes){
			next if $node->is_Leaf;					# no collapsign of leaves directly

			# collapsing node #
			my @rm_list;
			if($node->height < $i){
				my $rm_num = 0;
				# removing leaf nodes 1st #
				for my $child ( $node->each_Descendent ){
					next unless $child->is_Leaf;
					
					push(@rm_list, $child->id);
					
						#print $child->id, "\n";
					$treeo_tmp->remove_Node($child);	
					$rm_num++;
					}
					
				# removing internal nodes #
				if ($rm_num){
					for my $child ( $node->each_Descendent ){	
						push(@rm_list, $child->id) if $child->id;
						$treeo_tmp->remove_Node($child);
						}
					}
				
				# new name for node #
				next unless @rm_list;		# only if something was removed
				
				# appending leaf names to node #
				if($node->id){ 
					$node->id( join("__", $node->id, @rm_list)); }
				else{ 
					$node->id( join("__", @rm_list)); 
					}
				}
			}
		
		# if N-nodes == $nclade, return nodes #
		my $nnodes = scalar $treeo_tmp->get_leaf_nodes;
		print STDERR "...Branch length cutoff: $i; Number of clades: $nnodes\n";
		return $treeo_tmp if $nnodes <= $nclade;
		}
		
	die " ERROR: the desired number of clades ($nclade) could not be found!\n";
	}
	
sub max_height{
	my ($treeo) = @_;
	
	# getting max branch length between nodes #
	my $max_brlen = 0;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		$max_brlen = $node->height if $node->height > $max_brlen;
		}
	print STDERR "...max branch length between leaves: $max_brlen\n";
	
	return $max_brlen;
	}	

sub max_brlen{
	my ($treeo) = @_;
	
	# getting max branch length between nodes #
	my $max_brlen = 0;
	for my $node1 ($treeo->get_leaf_nodes){
		for my $node2 ($treeo->get_leaf_nodes){
			next if $node1->id eq $node2->id;
			my $dist = $treeo->distance(-nodes => [$node1, $node2]);
			$max_brlen = $dist if $dist > $max_brlen;
			}
		}
	print STDERR "...max branch length between leaves: $max_brlen\n";
	
	return $max_brlen;
	}

sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO->new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		
	return $treeio;
	}

sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}


__END__

=pod

=head1 NAME

tree_collapse_tips.pl -- collapsing tips by a branch length cutoff

=head1 SYNOPSIS

tree_collapse_tips.pl [flags]

=head2 require flags

=over

=item -t

Tree file (newick or nexus).

=back

=head2 optional flags

=over

=item -format

Tree file format (newick or nexus). [newick]

=item -length

Branch length cutoff. [< 0]

=item -count

Count file in Mothur format

=item -meta

Metadata file in ITOL format

=item -color

Color file in ITOL format

=item -regex

Regular expression for excluding taxa from being collapsed.

=item -branch

Write out all branch lengths (internal nodes to most distant decendent).

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_collapse_tips.pl

=head1 DESCRIPTION

Collapse branches in a tree that have a branch length of < (-length)
from the ancestral node.

Collapsed nodes are labeled as: "collapsed-nodeID"_"number_taxa_collapsed"

If any metadata files are provided (count, ITOL-metadata, color), 
the taxon labels are updated and the abundances are summed. A random 
taxon in a collapsed clade will be used for the color info if a
color file is provided.

Output file names are based on input file names ("*br-col*")

=head2 WARNING

'-count' flag not fully tested!

=head1 EXAMPLES

=head2 Getting branch length distribution

tree_collapse_tips.pl -t tree.nwk -branch

=head2 Collapse tree with a metadata file

tree_collapse_tips.pl -t tree.nwk -meta meta.txt 

=head2 Collapse tree with a count file

tree_collapse_tips.pl -t tree.nwk -count count.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

