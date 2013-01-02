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

my ($verbose, $tree_in, $habitat_in, $full_in, $format, $abund_cut);
GetOptions(
	   "tree=s" => \$tree_in,		# tree file
	   "habitat=s" => \$habitat_in,	# habitat file (2column)
	   "format=s" => \$format,		# tree format
	   "cutoff=s" => \$abund_cut,	# abundance cutoff for filtering habitat
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);
die " ERROR: Provide a habitat file (from AdaptML)." if ! $habitat_in;
$abund_cut = 20 if ! $abund_cut;

### MAIN
my $treeo = tree_io($tree_in, $format);
my $habs_ref = load_habitat($habitat_in);
$treeo = find_pops($treeo, $habs_ref, $abund_cut);
#print_node_ids($treeo, "true", "true");
#$habs_ref = filter_habs($habs_ref, $abund_cut);				# filtering by abundance
exit;
#get_leaf_branch_lengths($treeo);

### Subroutines
sub find_pops{
	my ($treeo, $habs_ref, $abund_cut) = @_;
	
	my $root = $treeo->get_root_node;
	#my %groups;
	#my $groups_ref = groups_of_ten($treeo, $root, \%groups);
	
	my %pops;
	my $pops_ref = check_monophyl($root, $treeo, \%pops, $habs_ref, $abund_cut);
	
		#print Dumper $pops_ref; exit;
	
	}

sub groups_of_ten{
	my ($treeo, $node, $groups_ref) = @_;
	
	if(scalar $node->get_all_Descendents <= 20){
		$node->id(scalar keys %$groups_ref);
		$$groups_ref{$node->id} = 1;
		}
	else{
		for my $child ($node->each_Descendent){
			$groups_ref = groups_of_ten($treeo, $child, $groups_ref);
			}
		}
	return $groups_ref;
	}


sub check_monophyl{
	my ($node, $treeo, $pops_ref, $habs_ref, $abund_cut) = @_;
	
	# checking for monophyl #
	my %node_habs;
	for my $child ($node->get_all_Descendents){			# getting all leaves
		next if ! $child->is_Leaf;
		if(! exists $$habs_ref{$child->id}){
			print STDERR " WARNING: ", $child->id, " was not found in habitat file; skipping\n";
			next;
			}
		else{	
			$node_habs{$$habs_ref{$child->id}}{"count"}++;
			push(@{$node_habs{$$habs_ref{$child->id}}{"taxa"}}, $child->id);
			}
		}
	
	# filtering #
	my $node_habs_sum = 0;
	my @taxa;
	foreach my $key (keys %node_habs){
		if ($node_habs{$key}{"count"} < $abund_cut){
			print STDERR "Filtering: ", join(",", @{$node_habs{$key}{"taxa"}}),"\n";
			delete $node_habs{$key};
			}
		else{
			$node_habs_sum += $node_habs{$key}{"count"};
			push(@taxa, @{$node_habs{$key}{"taxa"}});
			}
		}

	# if monophyl, note node; else continue down tree #
	if(scalar keys %node_habs == 1){
		# labeling node #
		$node->id(scalar keys %$pops_ref);

		# getting node characteristics #
		$$pops_ref{$node->id}{"taxa_count"} = $node_habs_sum;			# number of taxa in population
		$$pops_ref{$node->id}{"taxa"} = \@taxa;							# taxa in node
		$$pops_ref{$node->id}{"root_dist"} = $node->depth;				# number of taxa in population
		
		}
	else{					# continue down tree if not monophyletic
		for my $child ($node->each_Descendent){
			$pops_ref = check_monophyl($child, $treeo, $pops_ref, $habs_ref, $abund_cut);
			}
		}
	return $pops_ref;
	}

sub check_monophyl_OLD{
	my ($node, $habs_ref, $abund_cut) = @_;
	
	my %hab_cnt;
	for my $child ($node->each_Descendent){
		if(! $child->is_Leaf){		# if node
			if($child->id){
				$hab_cnt{$child->id}++; 		# getting habitat of node
				}
			else{							# getting all descendents
				for my $node_child ($node->get_all_Descendents){
				
					}
				}
			} 
		else{						# if leaf
			if(! exists $$habs_ref{$child->id}){
				print STDERR " WARNING: ", $child->id, " not found in habitat file\n";
				}
			else{
				$hab_cnt{ $$habs_ref{$child->id} }++;
				}
			}
		}
	
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

sub filter_habs{
# filtering by abundance cutoff #
	my ($habs_ref, $abund_cut) = @_;
	
	# making a hab count table #
	my %hab_cnt;
	map{$hab_cnt{$_}++ } values %$habs_ref; 

	print Dumper %hab_cnt; exit;
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
	return \%habs;
	}


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

