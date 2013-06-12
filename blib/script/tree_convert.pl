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

my ($informat, $outformat, $verbose);
GetOptions(
		"in=s" => \$informat,
		"out=s" => \$outformat,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$informat = check_format($informat);
$outformat = check_format($outformat);

### MAIN
my $in = new Bio::TreeIO(-fh => \*STDIN,
						-format => $informat);
my $out = new Bio::TreeIO(-fh => \*STDOUT,
						-format => $outformat);
						
while(my $tree = $in->next_tree){
	$out->write_tree($tree);
	}


### subroutines
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

tree_convert.pl -- convert tree format (newick <-> nexus)

=head1 SYNOPSIS

tree_convert.pl -i -o < tree_in > tree_out

=head2 options

=over

=item -i 	Input tree format (newick|nexus)

=item -o 	Output tree format (newick|nexus)

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_convert.pl

=head1 DESCRIPTION

Simple Bio::TreeIO script for converting the format of a tree file.
Accepted formats: newick or nexus.
Files can contain multiple trees.

=head1 EXAMPLES

=head2 Usage method 1

tree_convert.pl -i newick -o nexus < in.nwk > out.tre

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

