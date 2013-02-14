#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use List::Util qw/sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $tformat, $count_in, $count_header, $mothur);
my $abund_cut = 5;
GetOptions(
	   "tree=s" => \$tree_in,			# tree file
	   "format=s" => \$tformat,			# tree format
	   "count=s" => \$count_in,			# count table file
	   "xheader" => \$count_header,		# header in count table? [T]
	   "mothur" => \$mothur,			# mothur formatted count file? [F]
	   "abundance=i" => \$abund_cut, 	# abundance cutoff 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " Provide a tree file (newick or nexus).\n" if ! $tree_in;
$tformat = check_tree_format($tformat);

### MAIN
# load tree & count files #
my $treeo = tree_io($tree_in, $tformat);
my ($count_r, $header) = load_count($count_in, $count_header, $mothur);

# sum abundances #
my $abund_r = sum_count($count_r, $abund_cut, $mothur);

# checking taxon existence in tree #
check_names($treeo, $abund_r, $count_r);

# writing out prune list for R script & calling script #
my $prune_file = write_prune_list($count_in, $abund_r);			
call_prune_tree($tree_in, $prune_file);							# calling R script

# writting out pruned count file #
write_pruned_count_file($count_r, $count_in, $header);

#$treeo = prune_by_abundance($treeo, $count_r);
#tree_write($treeo, $tree_in);

### Subroutines
sub write_pruned_count_file{
	my ($count_r, $count_in, $header) = @_;
	
	(my $outfile = $count_in) =~ s/\.[^.]+$|$/_prn.txt/;
	open OUT, ">$outfile" or die $!;
	
	print OUT $header, "\n" if $header;
	foreach my $row (keys %$count_r){
		print OUT join("\t", $row, @{$$count_r{$row}}), "\n";
		}
	close OUT;	
		
	print STDERR " Pruned count file written: '$outfile'\n";
	print STDERR "  Number of taxa in pruned count file: ", scalar keys %$count_r, "\n";
	}

sub call_prune_tree{
# calling prune_tree.r #
	my ($tree_in, $prune_file) = @_;
	
	my $cmd = "Rscript ../../bin/tree_PruneByAbundance.r -t $tree_in -n $prune_file";
#	my $cmd = "tree_PruneByAbundance.r -t $tree_in -n $prune_file";
	print STDERR "\n$cmd\n" if $verbose;
	system($cmd);
	
	}

sub write_prune_list{
# writing out prune list #
	my ($count_in, $count_r) = @_;
	
	(my $outfile = $count_in) =~ s/\.[^.]+$|$/_prn-list.txt/;
	open OUT, ">$outfile" or die $!;
	
	foreach my $taxon (keys %$count_r){
		print OUT join("\t", $taxon, $$count_r{$taxon}), "\n";
		}
	close OUT;
	
	return $outfile;
	}

sub tree_write{
	### writting out a nexus tree file ###
	my ($treeo, $tree_in) = @_;

	(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_prn.nwk/;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	print STDERR " Newick tree file written: '$outfile'\n";

	}

sub prune_by_abundance{
# pruning the phylogeny by abundance #
	my ($treeo, $count_r) = @_;
	
	
	for my $node ($treeo->get_leaf_nodes){
		next if ! exists $$count_r{$node->id};
		
		$treeo->remove_Node($node);
		}
	
	return $treeo;
	}

sub check_names{
# checking names #
	my ($treeo, $abund_r, $count_r) = @_;

	my $count_rows = scalar keys %$count_r;

	my %nodes = map{$_->id, 1} $treeo->get_leaf_nodes;	

	my $prune_cnt = 0;
	foreach my $taxon (keys %$abund_r){
		if (! exists $nodes{$taxon}){
			print STDERR " WARNING! not found in tree file: '$taxon'\n";
			delete $$abund_r{$taxon};			# removing from count file because taxon not in tree
			delete $$count_r{$taxon};			# removing from count file because taxon not in tree
			}
		else{
			$prune_cnt++ if $$abund_r{$taxon} eq "delete";
			}
		}
	
	print STDERR "\n Number of taxa in count file: $count_rows\n";
	print STDERR " Number of leaves: ", scalar keys %nodes, "\n";
	print STDERR " Number of taxa to be pruned: $prune_cnt\n";
	}

sub sum_count{
# summing abundances in count file #
	my ($count_r, $abund_cut, $mothur) = @_;
	
	my %abund;
	foreach my $row (keys %$count_r){
		my $rowsum;
		if($mothur){ $rowsum = sum(@{$$count_r{$row}}[1..$#{$$count_r{$row}}]); }
		else{ $rowsum = sum( @{$$count_r{$row}} ); }
		
		
		if($rowsum >= $abund_cut){
			$abund{$row} = $row;
			}
		else{
			$abund{$row} = "delete";
			}
		}
		
		#print Dumper %abund; exit;
	return \%abund;
	}

sub load_count{
# loading count file #
# 
	my ($count_in, $count_header, $mothur) = @_;
	open IN, $count_in or die $!;
	
	my %count;
	my $header;
	while(<IN>){
		chomp;
		
		# header #
		if ($.==1 && ! $count_header && ! $mothur){
			$header = $_;
			next;
			}
		
		my @line = split /\t/;
		die " ERROR: the count file must be at least 2 columns (rownames, count)\n"
			if scalar @line < 2;
			
		$count{$line[0]} = [@line[1..$#line]];
		}
	close IN;
			
		#print Dumper %count; exit;
	return \%count, $header;			# returning taxa for pruning
	}

sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		
		#for my $node ($treeio->get_nodes){ print "nodeID: ", $node->id, "\n"; }
		#exit;	
	return $treeio;
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

tree_PruneByAbundance.pl -- Prune tree by taxon abundances

=head1 SYNOPSIS

tree_PruneByAbundance.pl -t -c [-f] [-x] [-m] [-a]

=head2 options

=over

=item -t

Tree file (newick or nexus).

=item -f

Tree file format (newick or nexus). [newick]

=item -c

Count file (tab-delimited, 1st row = rownames)

=item -x

Header in 1st line of count file? [FALSE]

=item -m

Mothur-formatted count file? [FALSE]

=item -a

Abundance cutoff for pruning (>=). [5]

=item -v

Verbose output

=item -h

This help message

=back

=head2 Requirements:

prune_TreeByAbundance.r

=head2 For more information:

perldoc tree_PruneByAbundance.pl

=head1 DESCRIPTION

Prune taxa from a tree based on the abundances of the taxa found.

=head1 EXAMPLES

=head2 Basic usage

tree_PruneByAbundance.pl -t test.nwk -c count.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

