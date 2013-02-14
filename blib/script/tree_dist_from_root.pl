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

my ($verbose, $tree_in, $habitat_in, $full_in, $format);
GetOptions(
	   "tree=s" => \$tree_in,		# tree file
	   #"habitat=s" => \$habitat_in,	# habitat file (2column)
	   "x=s" => \$full_in,		# full file (from adaptML)
	   "format=s" => \$format,		# tree format
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);


### MAIN
my $treeo = tree_io($tree_in, $format);
#my $habs_ref = load_habitat($habitat_in) if $habitat_in;
my $full_ref = load_full_file($full_in);
my $depths_ref = get_node_distances($treeo, $full_ref);
write_depths($depths_ref);

### Subroutines
sub write_depths{
# writing out branch length from not to root (parsed by habitat #
	my ($depths_ref) = @_;
	
	print join("\t", qw/Habitat Root_distance/), "\n";
	foreach my $hab (keys %$depths_ref){
		foreach my $depth (@{$$depths_ref{$hab}}){
			print join("\t", $hab, $depth), "\n";
			}
		}
	}

sub get_node_distances{
# getting node distances from root for each not specified in full file #
	my ($treeo, $full_ref) = @_;
	
	my $root = $treeo->get_root_node;
	
	#getting nodes #
	my %depths;
	my $cnt = 0;
	foreach (@$full_ref){
		$cnt++;
		my @nodes1 = $treeo->find_node(-id => $$_[0]);
		my @nodes2 = $treeo->find_node(-id => $$_[1]);
		next if ! @nodes1 || ! @nodes2;
		die " ERROR: $_ is found 2x in the tree file\n"
			if scalar @nodes1 > 1 || scalar @nodes2 > 1;
		my $lca = $treeo->get_lca(-nodes => [$nodes1[0], $nodes2[0]] );
		
		# loading distances from root #
		push(@{$depths{$$_[2]}}, $lca->depth);
		
		# status #
			#last if $cnt == 100;
		print STDERR "Nodes processed: $cnt of ", scalar @$full_ref, "\n" 
			if $cnt % 100 == 0 and $verbose;
		}

		#print Dumper %depths;
	return \%depths;
	}

sub load_full_file{
	my $full_in = shift;
	
	open IN, $full_in or die $!;
	my @full;
	my @labels;
	while(<IN>){
		chomp;
		if(/^LABELS/){
			@labels = split /,/;
			}
		elsif(/^COLORS/){ next; }
		else{
			next unless /\|/;		# must be a node
			my @line = split /,/;
			my @taxa = split /\|/, $line[0];
			for my $i (2..$#line){
				push(@full, [$taxa[0], $taxa[1], $labels[$i-1]]) if $line[$i] > 0;
				}
			}
		}
	close IN;

		#print Dumper scalar @full; exit;
	return \@full;
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
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;	
	return $treeio;
	}
	
sub load_habitat{
## loading habitat file ##
	my $habitat_in = shift;
	open IN, $habitat_in or die $!;
	my %habs;
	while(<IN>){
		chomp;
		s/#.+//;
		next if /^\s$/;
		my @line = split /\t/;
		die " ERROR: habitat file not formated correctly\n" if scalar @line != 2;
		$habs{$line[0]} = $line[1];
		}
	close IN;	
		
	return \%habs;
	}

__END__

=pod

=head1 NAME

tree_dist_from_root.pl -- Get the branch length from root to each node (for AdaptML habitats)

=head1 SYNOPSIS

tree_dist_from_root.pl [-f] [-v] -t -x > output

=head2 options

=over

=item -t 	Tree file (newick or nexus)

=item -f 	Tree file format [newick]

=item -x 	Full file (from AdaptML)

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_dist_from_root.pl

=head1 DESCRIPTION

The script will find the branch lengths from a particular node to the root.
The purpose is to determine whether particular habitats are ancestral.

=head1 EXAMPLES

tree_dist_from_root.pl -t itol.tree -full full.file -v > node-root_dist.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

