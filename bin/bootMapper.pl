#!/usr/bin/env perl

=pod

=head1 NAME

tree_bootMapper.pl -- script template

=head1 SYNOPSIS

tree_bootMapper.pl [options] -ref has_boot.nwk -query needs_boot.nwk > needs_boot_with-boot.nwk

=head2 Required flags

=over

=item -ref

Reference tree with bootstrap values.

=item -query

Query tree that bootstrap values will be mapped onto.

=back

=head2 Optional flags

=over

=item -format

Tree formats (2 values; 'newick' or 'nexus'). ['newick' 'newick']

=item -missing

Internal nodes with no bootstrap values after mapping will get this value. []

=item -internal

Internal node IDs are bootstrap values. [TRUE]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_bootMapper.pl

=head1 DESCRIPTION

Map the bootstrap values on internal nodes of a reference
tree to the internal nodes of a query tree containing
the same taxa.

This is needed if bootstrap values are stripped by a program
and need to be added back.

=head1 EXAMPLES

=head2 Basic usage:

tree_bootMapper.pl -ref ref.nwk -query query.nwk -missing 0 > query_wboot.nwk

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $rtree_in, $qtree_in, $internal_b, $missing_boot);
my @format = ("newick", "newick");
GetOptions(
	"reference=s" => \$rtree_in,
	"query=s" => \$qtree_in,
	"format=s{2,2}" => \@format,
	"missing=i" => \$missing_boot,
	"internal" => \$internal_b,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide a reference tree file!\n" unless $rtree_in;
die "ERROR: provide a query tree file!\n" unless $qtree_in;
map{ check_format($_) } @format;

#--- MAIN ---#
my $rtreeo = tree_io($rtree_in, $format[0], $internal_b);
my $qtreeo = tree_io($qtree_in, $format[1], $internal_b);

# getting all lca bootstrap values in reference #
my $lca_r = get_all_lca($rtreeo);

# mapping on bootstrap values #
map_bootstrap($qtreeo, $lca_r);

# adding bootstrap to nodes missing a value #
add_missing_bootstrap($qtreeo, $missing_boot) if defined $missing_boot;

# moving bootstrap to id #
move_bootstrap_to_id($qtreeo);

# writing out editting query tree #
tree_write($qtreeo);

#--- Subroutines ---#
sub tree_write{
	my ($treeo) = @_;
	my $out = new Bio::TreeIO(-fh => \*STDOUT,
							  -format => "newick");
	$out->write_tree($treeo);
	}
	
sub move_bootstrap_to_id{
	my $treeo = shift;

	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		
		$node->id($node->bootstrap);
		}
	}

sub add_missing_bootstrap{
	my ($qtreeo, $missing_boot) = @_;
	
	for my $node ($qtreeo->get_nodes){
		next if $node->is_Leaf;
		unless(defined $node->bootstrap){
			$node->bootstrap($missing_boot);
			}
		}
	}

sub map_bootstrap{
# mapping bootstrap values onto query tree #
	my ($qtreeo, $lca_r) = @_;
	
	# status #
	print STDERR "Adding bootstrap support values to query tree\n";
	
	my @lnodes = $qtreeo->get_leaf_nodes;
	for my $i (0..$#lnodes){
		for my $ii (0..$#lnodes){
			next if $i <= $ii;		# lower triangle
			my $lca_boot;
			
			# if query LCA pair exists in reference tree #
			if(exists $lca_r->{$lnodes[$i]->id}{$lnodes[$ii]->id}){
				$lca_boot = $lca_r->{$lnodes[$i]->id}{$lnodes[$ii]->id};				
				}
			elsif(exists $lca_r->{$lnodes[$ii]->id}{$lnodes[$i]->id}){
				$lca_boot = $lca_r->{$lnodes[$ii]->id}{$lnodes[$i]->id};
				}
			
			# if LCA reference bootstrap, add to query #
			if(defined $lca_boot){
				# adding bootstrap support to LCA on query tree #
				my $lca = $qtreeo->get_lca(-nodes => [$lnodes[$i], $lnodes[$ii]]);
				$lca->bootstrap($lca_boot);
				
				print STDERR "LCA of: ", $lnodes[$i]->id, "<->", $lnodes[$ii]->id,
					" bootstrap assigned: ", $lca->bootstrap, "\n" if $verbose_b;
				#print Dumper $lnodes[$i]->id, $lnodes[$ii]->id, $lca->bootstrap;
				}
			}
		}	

	}

sub get_all_lca{
# getting all lca bootstrap for each LCA of each strain pair #
	my ($treeo) = @_;
	
	my @lnodes = $treeo->get_leaf_nodes;
	# status #
	print STDERR "Number of leaves in reference tree: ", 
		scalar @lnodes, "\n";
	print STDERR "Number of LCAs to check: ", 
		(scalar @lnodes) * (scalar @lnodes -1) / 2, "\n";
	
	my %lca;
	for my $i (0..$#lnodes){
		for my $ii (0..$#lnodes){
			next if $i <= $ii;		# lower triangle
			my $lca = $treeo->get_lca(-nodes => [$lnodes[$i], $lnodes[$ii]]);
			if(defined $lca && defined $lca->bootstrap){
				#push @lca, [$lnodes[$i]->id, $lnodes[$ii]->id, $lca->bootstrap];
				$lca{$lnodes[$i]->id}{$lnodes[$ii]->id} = $lca->bootstrap;
				}
			else{
				print STDERR "No bootstrap values for: ", $lnodes[$i]->id, 
							" <-> ", $lnodes[$ii]->id, "\n" if $verbose_b;
				}
			}
		}
		
	return \%lca;
	}



sub check_format{
# checking format input #
	my $format = shift;
	if($format =~ /^new/i){ $format = "newick"; }
	elsif($format =~ /^nex/i){ $format = "nexus"; }
	else{ die " ERROR: '$format' format not recognized (newick | nexus)\n"; }
	return $format;
	}

sub tree_io{
# loading tree object #
	my ($tree_in, $format, $internal_b) = @_;
	my $input;
	if(! $internal_b){			# on by default
		$input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format,
								-internal_node_id => 'bootstrap');	
		}
	else{
		$input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
		}
	my $treeo = $input->next_tree;	
	return $treeo;
	}
	

