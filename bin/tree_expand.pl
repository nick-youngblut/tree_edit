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

my ($verbose, $tree_in, $count_in, $format);
GetOptions(
	   "tree=s" => \$tree_in,		# tree file
	   "count=s" => \$count_in,		# habitat file (2column)
	   "format=s" => \$format,		# tree format
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide tree file (newick or nexus format)." if ! $tree_in;
if(! $format){ $format = "newick"; }
$format = check_format($format);
die " ERROR: Provide a count file (Mothur format)." if ! $count_in;

### I/O error & defaults
(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_exp.nwk/;

### MAIN
my $treeo = tree_io($tree_in, $format);
my $count_ref = load_count($count_in);
$treeo = expand_taxa($treeo, $count_ref);
	#write_root($treeo);
tree_write($treeo, $outfile);

### Subroutines
sub tree_write{
	### writting out a newick tree file ###
	my ($treeo, $outfile) = @_;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	print STDERR " Newick tree file written:\n  $outfile\n";
	}

sub write_root{
	my $treeo = shift;
	my $root = $treeo->get_root_node;
	print "Root_id = ", $root->id, "\n";
	}

sub expand_taxa{
# expanding (replicating) taxa so that there is a rep (clonal) for each sample #
	my ($treeo, $count_ref) = @_;
	
	my $leaf_cnt = 0;
	for my $leaf ($treeo->get_leaf_nodes){
		#print Dumper $leaf->id; 
		$leaf_cnt++;
		
		# finding ancestor #
		my $anc = $leaf->ancestor;
		
		# getting leaf info #
		my $leaf_id = $leaf->id;
		my $leaf_br = $leaf->branch_length;
		
		# adding clonal replicate taxa (for each sample) #
		die " ERROR: ", $leaf_id, " not found in count table\n" if ! exists $$count_ref{$leaf_id};
		my $samp_cnt = 64;			 # naming
		my $tojoin = "";			# naming
		foreach my $samp (@{$$count_ref{$leaf_id}}){
			$samp_cnt++;
			
			# adding clonal node #
			#$anc->add_Descendent($anc->new(-branch_length => $leaf->branch_length));
			
			# adding clonal taxa #
			$leaf->add_Descendent($leaf->new(-id => join("_", $samp, join("", $leaf_cnt, $tojoin, chr $samp_cnt)), 
					-branch_length => 0));
			
			# naming counter (AA, AAA, etc) #
			if($samp_cnt == 90){
				$tojoin .= chr $samp_cnt;
				$samp_cnt = 64;
				}
			
			# renaming leaf #
			$leaf->id("");
			}

		}
		
	return $treeo;
	}

sub load_count{
# loading count table (mothur format #
	my $count_in = shift;
	
	open IN, $count_in or die $!;
	my %cnt_tbl;
	my @header;
	while(<IN>){
		chomp;
		if($.==1){	# header
			@header = split /\t/;
			}
		else{
			my @line = split /\t/;
			for my $i (2..$#line){
				push(@{$cnt_tbl{$line[0]}}, $header[$i]) if $line[$i] > 0;
				}
			}
		}
	close IN;
		#print Dumper %cnt_tbl; exit;
	return \%cnt_tbl;
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

