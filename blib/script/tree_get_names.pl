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

my ($tree_in, $format, $verbose, $multi);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "multi" => \$multi,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$format = check_tree_format($format);

### MAIN
my $treeio = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
my %names;
for my $treeo ($treeio->next_tree){
	for my $node ($treeo->get_nodes){
		next unless $node->is_Leaf;
		$names{ $node->id }++;
		}
	last unless $multi;
	}
write_names(\%names, $multi);

### Subroutine
sub write_names{
# writing out names from tree file #
	my ($names_r, $multi) = @_;
	
	foreach my $name (sort{$names_r->{$b} <=> $names_r->{$a}} keys %$names_r){
		if($multi){ print join("\t", $name, $names_r->{$name}), "\n";}
		else{ print $name, "\n"; }
		}

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

tree_get_names.pl -- get leaf names from a newick/nexus tree file

=head1 SYNOPSIS

tree_get_names.pl -t [-f] [-m]  > names.txt

=head2 options

=over

=item -t 	Tree file.

=item -f 	Format (newick | nexus). [newick]

=item -m 	Multiple trees. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc tree_get_names.pl

=head1 DESCRIPTION

Get a list of leaf names in a newick or nexus tree.

If the file contains multiple trees, the first tree will be
used unless the '-m' flag is used, which will count the 
number of trees that each name appears in and that count 
will be written as a second column (good for checking that concatentated
trees all have the same names.

=head1 EXAMPLES

=head2 Usage: standard

tree_get_names.pl -t file.nwk > names.txt

=head2 Usage: count of name prevalence among trees in file

tree_get_names.pl -t file.nwk -multi > names.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

