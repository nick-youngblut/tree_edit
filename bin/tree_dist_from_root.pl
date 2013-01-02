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
	   "habitat=s" => \$habitat_in,	# habitat file (2column)
	   "full=s" => \$full_in,		# full file (from adaptML)
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
my $habs_ref = load_habitat($habitat_in) if $habitat_in;
$treeo = internal_node_label($treeo, $habs_ref);
#print_node_ids($treeo, "true", "true");
exit;
#get_leaf_branch_lengths($treeo);

### Subroutines
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

sub internal_node_label{
	# using tip labels to label internal nodes #
	my ($treeo, $habs_ref) = @_;
	my @leaf_names = keys %$habs_ref;
	foreach my @nodes(
		for my $ii (0..$#leaf_names){
			next if $i >= $ii;
			my $lca = $treeo->get_lca(-nodes => $leaf_names[$i], $leaf_names[$ii]);
			print Dumper $lca; exit;
			}
		}
		#exit;
	#return $tree;
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

