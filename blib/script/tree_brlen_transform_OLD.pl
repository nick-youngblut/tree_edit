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

my ($verbose, $tree_in, $find_max, $cutoff, $tips);
my $format = "newick";
my $transform = "log";
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "transform=s" => \$transform,		# tranforming factor
	   "cutoff=i" => \$cutoff,				# min cutoff for brlen modification
	   "max" => \$find_max,					# just find max brlen?
	   "tips" => \$tips,					# just tips? [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die "Provide tree file (newick or nexus format)." if $tree_in;
$format = check_format($format);

### MAIN
my $treeo = tree_io($tree_in, $format);
find_max($treeo, $tips) if $find_max;


### Subroutines
sub find_max{
# finding max brlen; then exiting #
	my ($treeo, $tips) = @_;
	
	for my $node ($treeo->get_nodes){
		print Dumper $node; exit;
		}
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

