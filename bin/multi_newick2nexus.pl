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

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );


### I/O error & defaults
### MAIN
# checking names #
#(my $infile = $ARGV[0])
my $in = new Bio::TreeIO(-file => $ARGV[0],
						-format => "newick");
						
my $names_r = check_names($in);

# converting newick to merged nexus #
convert_newick($ARGV[0], $names_r);


### Subroutines
sub convert_newick{
	my ($infile, $names_r) = @_;
	open IN, $infile or die $!;

	# writing header #
	print "#NEXUS\n\n";
	print "BEGIN taxa;\n";
	print "DIMENSIONS ntax=", scalar keys %$names_r, ";\n";
	print "TAXLABELS\n";

	my @names = keys %$names_r;
	for my $i (0..$#names){
		print "[",$i+1,"] '$names[$i]'\n";
		}
	print ";\nEND [taxa];\n";

	# writing trees #
	print "BEGIN trees;\n";
	while(<IN>){
		print "tree t$.=$_";
		}
		
	# writing footer #
	print "END [trees];\n";
	}

sub check_names{
# checking to make sure all trees have the same taxa names #
	my $in = shift;

	my %names;
	my $tree_cnt = 0;

	# getting leaf names #
	while(my $tree = $in->next_tree){
		$tree_cnt++;
		for my $node ($tree->get_nodes){
			next unless $node->is_Leaf;
			$names{$node->id}++;
			}
		}
		
	# checking that all leaf name counts = N-trees 
	foreach my $leaf (keys %names){
		die " ERROR: $tree_cnt trees, but $leaf found $names{$leaf} times\n"
			unless $names{$leaf} == $tree_cnt;
		}
	#	print Dumper %names; exit;
	return \%names;
	}


__END__

=pod

=head1 NAME

tree_multi_newick2nexus.pl -- convert a multi-tree newick to (new) nexus format

=head1 SYNOPSIS

tree_multi_newick2nexus.pl input > output

=head2 options

=over

=item -h	This help message

=back

=head2 For more information:

perldoc tree_multi_newick2nexus.pl

=head1 DESCRIPTION

Convert a newick tree file containing multiple trees to a nexus file
readable by Splitstree.

=head1 EXAMPLES

=head2 Basic Usage:

tree_multi_newick2nexus.pl file.nwk > file.nex

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

