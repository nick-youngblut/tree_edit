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

my ($verbose, $tree_in, $full_in, $format, $abund_cut);
GetOptions(
	   "tree=s" => \$tree_in,		# tree file
	   "full=s" => \$full_in,		# habitat file (2column)
	   "x=s" => \$format,			# tree format
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
$format = "newick" if ! $format;
$format = check_format($format);
die " ERROR: Provide a 'full' file (from AdaptML)." if ! $full_in;


### MAIN
my $treeo = tree_io($tree_in, $format);
my ($full_r,$hab_index_r) = load_full($full_in);
#$node_list_r = get_node_tips($treeo, $full_r)
$treeo = label_lca($treeo, $full_r);
transverse_tree($treeo, $hab_index_r);

### Subroutines
sub transverse_tree{
# transversing a tree and finding lineage histories #
	my ($treeo, $hab_index_r) = @_;
	
	# header #
	print join("\t", qw/tip_label bifur_num habitat color/), "\n";
	
	# getting lineages #
	my %trans;
	for my $child ($treeo->get_nodes){
		next if ! $child->is_Leaf;		# just tips
		my @nodes = $treeo->get_lineage_nodes($child);	# getting nodes from root to tip
		next if ! @nodes;
		
		$trans{$child->id} = \@nodes;
		}
		
	# writing out table; sort by number of bifurcations #
	foreach my $id (sort{scalar @{$trans{$a}} <=> scalar @{$trans{$b}}} keys %trans){
		for my $i (1..$#{$trans{$id}}){
			print join("\t", $id, $i, ${$trans{$id}}[$i]->id, $$hab_index_r{${$trans{$id}}[$i]->id}), "\n"
			}
		}
	}

sub label_lca{
# labeling the lca of all nodes by habitat #
	my ($treeo, $full_r) = @_;

	foreach my $taxa (@$full_r){	# foreach taxon set
		# getting nodes #
		my @noi;		# nodes of interest
		foreach my $taxon (@{$$taxa{"taxa"}}){
			my @nodes = $treeo->find_node(-id => $taxon);
			if(! @nodes){
				print STDERR " WARNING: $taxon not found! Skipping\n" if ! @nodes;
				next;
				}
			die " ERROR: >1 hit for $taxon!\n" if scalar @nodes > 1;
			push(@noi, $nodes[0]);
			}
		next if scalar @noi != 2;

		# lca ID #
		my $lca = $treeo->get_lca(-nodes => \@noi);
		die " ERROR: LCA: $lca already has an ID!\n" if $lca->id;
		$lca->id($$taxa{"habitat"});
		
		# status #
		print STDERR " Processing taxa: ", join(", ", @{$$taxa{"taxa"}}), "\n" if $verbose;
		}
	
	return $treeo;
	}


sub load_full{
# loading full file from adaptML #
	my $full_in = shift;
	open IN, $full_in or die $!;
	
	my (@full, @header, @colors);
	while(<IN>){
		chomp;
		s/^ +//;
		my @line = split /,/;
		if(/^LABELS/){
			@header = @line;
			}
		elsif(/^COLORS/){	
			@colors = @line;
			}
		elsif(! /^[^,]+\|/){ 	# if tip
			next;
			}
		else{
			my %tax_data;
			@{$tax_data{"taxa"}} = split /\|/, $line[0];
			for my $i (1..$#line){
				next if $line[$i] !~ /^[\d\.]+$/;
				if($line[$i] > 0){
					$tax_data{"habitat"} = $header[$i-1];
					#$tax_data{"color"} = $colors[$i-1];
					}
				}
			push(@full, \%tax_data);
			}
		}
	close IN;
	
	# making habitat-color index #
	my %hab_index;
	for my $i (1..$#header){
		$hab_index{$header[$i]} = $colors[$i];
		}
		#print Dumper @full; exit;
	return (\@full, \%hab_index);
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

adaptML_niche_transition.pl -- Plot niche transitions inferred by AdaptML

=head1 SYNOPSIS

adaptML_niche_transition.pl -t -f [-x] > out.txt

=head2 options

=over

=item -t 	Tree file (newick or nexus)

=item -x 	Tree file format [newick]

=item -f 	Full.file (AdaptML output)

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc adaptML_niche_transition.pl

=head1 DESCRIPTION

The script uses the full file to label the nodes by habitat.
The phylogeny is then transversed from root to tip for all tips and
the habitat at each node is recorded.

The script helps show evolutionary history of transitions among niches.

-head3 The output is in 4 columns:

=over

=item 1) Tip label

=item 2) The number of bifurcations from the root

=item 3) The inferred habitat

=item 4) The habitat's color (from full.file)

=back

=head1 EXAMPLES

=head2 Normal usage

adaptML_niche_transition.pl -t itol.tree -f full.file > niche_trans.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

