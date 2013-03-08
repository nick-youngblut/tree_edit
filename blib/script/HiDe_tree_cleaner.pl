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

my ($verbose, @gene_in, $species_in, $tformat, $species, $multi);
GetOptions(
	   "gene=s{,}" => \@gene_in,
	   "format=s" => \$tformat,
	   "species=s" => \$species_in,		# species tree?
	   #"multi" => \$multi, 			# multiple trees?
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " Provide a species & >=1 gene tree file (newick or nexus).\n" 
	if ! @gene_in || ! $species_in;
$tformat = check_tree_format($tformat);

#my $outfile;
#if($species){ $outfile = "species.newick"; }
#else{ ($outfile = $gene_in) =~ s/\.[^\.]+$|$/_cln.nwk/; } 
#unlink($outfile) if -e $outfile;

### MAIN
# loading species tree names #
#my $treeio = tree_io($gene_in, $tformat);
check_names($species_in, \@gene_in, $tformat);
	exit;

# editting and writing #
foreach my $gene_in (@gene_in){
my $treeio = tree_io($gene_in, $tformat);
my $out = new Bio::TreeIO(-fh => \*STDOUT, -format => "newick");
for my $treeo ($treeio->next_tree){

	$treeo = edit_leaf_labels($treeo);
	$treeo = remove_node_labels($treeo);
	$treeo = remove_bootstrap($treeo);
	$treeo = remove_branch_length($treeo);
	$treeo = check_tree_cleaning($treeo) if $verbose;

		#print Dumper $treeo; exit;
	$out->write_tree($treeo);
	}
}

#print STDERR " Newick tree file written:\n  $outfile\n";

### Subroutines
sub check_names{
	my ($species_in, $gene_in_r, $tformat) = @_;
	open IN, "tree_get_names.pl -t $species_in -f $tformat | " or die $!;
	while(<IN>){
		print Dumper $_;
		}
	close IN;
	

	}

sub remove_branch_length{
	my ($treeo) = @_;
	for my $node ($treeo->get_nodes){
		$node->branch_length("");
		}
	return $treeo;
	}

sub edit_leaf_labels{
	my ($treeo) = @_;
	for my $node ($treeo->get_nodes){
		next unless $node->is_Leaf;
		(my $id = $node->id) =~ s/[_.+-]/X/g;
		$node->id($id);
		}
	return $treeo;
	}

sub check_tree_cleaning{
# checking to see if cleaning worked #
	my ($treeo) = @_;
	
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		print "boot: ", $node->bootstrap, "\n"; # if $node->bootstrap;
		print "id: ", $node->id, "\n";	# if $node->id;
		print "brlen: ", $node->branch_length, "\n";
		}
	return $treeo;
	}

sub remove_bootstrap{
# removing bootstrap values #
	my ($treeo) = @_;
	for my $node ($treeo->get_nodes){
		$node->bootstrap("");
		}
	return $treeo;
	}

sub remove_node_labels{
# removing any node labels; replacing w/ "" #
	my ($treeo) = @_;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;		# skipping leaf nodes
		$node->id("");
		}
	return $treeo;
	}

sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}

sub tree_io{
	# loading tree object #
	my ($gene_in, $format) = @_;
	my $treeio = Bio::TreeIO -> new(-file => $gene_in,
								-format => $format);
	return $treeio;
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

