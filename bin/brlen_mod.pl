#!/usr/bin/perl
my $mod = "8/10/11 8:07 PM";
my $version = "0.4";
my $author = "Nick Youngblut";
#--------------------- version log ---------------------#
#
#
#-------------------------------------------------------#

### packages/perl_flags
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::TreeIO;

### global variables
my ($error);

### I/O
if ($#ARGV < 0){
	&usage;
	}
my ($tree_in, $format, $brlen_in, $verbose, $outfile);
GetOptions(
	   "tree=s" => \$tree_in,
	   "brlen=s" => \$brlen_in,
	   "verbose" => \$verbose,
	   "outfile=s" => \$outfile,
	   "help|?" => \&usage # Help
	   );

### Input error check
if(! $tree_in){
	$error = "Provide tree file (newick or nexus format).";
	error_routine($error, 1);
	}
if(! $format){ $format = "newick"; }
if(! $brlen_in){
	$error = "Provide a file of all branch lengths to include in tree.";
	error_routine($error, 1);
	}
if(! $outfile){
	($outfile = $tree_in) =~ s/\.[^\.]+$|$/_br-cor.nwk/;
	}
	
### Routing main subroutines
my $treeo = tree_io($tree_in, $format);
my $brlen_tbl = load_brlen_tbl($brlen_in);
$treeo = internal_node_label($treeo);
	#print_node_ids($treeo,"true", "true"); exit;
$treeo = mod_brlen($treeo, $brlen_tbl);
	#check_brlen($treeo); exit;
tree_write($treeo, $outfile);

#----------------------Subroutines----------------------#
sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;	
	return $treeio;
	}

sub load_brlen_tbl{
	# loading table of branch lengths #
	my $brlen_in = shift;
	my %brlen_tbl; # %
	open(FILE, $brlen_in) or die $!;
	my $line1 = <FILE>;
	while(<FILE>){
		chomp;
		my @tmp = split(/\t/);
		$tmp[1] =~ s/.+;//;
		$brlen_tbl{$tmp[1]} = 1 - abs($tmp[3]);
		}
		#print Dumper(%brlen_tbl); exit;
	return \%brlen_tbl;
	}
	
sub internal_node_label{
	# using tip labels to label internal nodes #
	my $tree = shift;
	foreach my $node ($tree->get_leaf_nodes){
			#print $node->id, "\n";
		my $rent = $node->ancestor;
		my @tmp = split(/-/, $node->id);
		foreach(@tmp){ $_ =~ s/(.+)_.+$/$1/; }
			#print $node->ancestor->id, "\n";
			#print "tmp: $tmp[4]\n";
		if(! $node->ancestor->id){
			$node->ancestor->id($tmp[4]);
			#print "rent: ", $node->ancestor->id, "\n";
			}
			#print $node->ancestor->id, "\n";
		if($tmp[0] eq $tmp[4]){
			#print $node->ancestor->id, "\n";
			$node->ancestor->branch_length(0);
			#print $node->ancestor->branch_length, "\n";
			}
		$node->id($tmp[4]);
		}
		#exit;
	return $tree;
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

sub mod_brlen{
	# adding branch lengths based on correlation values #
	my ($tree, $brlen_tbl) = @_;	#brlne_tbl(%)
		#print Dumper($brlen_tbl); exit;
	my @int_nodes;
	foreach my $node ($tree->get_nodes){		#brlen to internal nodes
		#(my $rent_name = $node->ancestor->id) =~ s/_i$//;
		#my $q = join("__", $rent_name, $node_name);
		#print "q: $q\n";
		#print "len: $$brlen_tbl{$q}\n";
		if($node->is_Leaf == 1){ next;}
		if(! exists($$brlen_tbl{$node->id})){ 
				#die "can't find node id -> ", $node->id, "\n"; 
			$node->branch_length(0);
				#print "node br ", $node->id, "\n";
			}
		else{ $node->branch_length($$brlen_tbl{$node->id});}	
		push(@int_nodes, $node->id);
		}
	#print Dumper(@int_nodes); exit;
	foreach my $node ($tree->get_leaf_nodes){	#leaf nodes
		my $q = $node->id;
		if(scalar(grep(/^$q$/, @int_nodes)) > 0){ 
			$node->branch_length(0);
				#print "added\n";
			}
		elsif(! $node->branch_length){ 
			if(exists($$brlen_tbl{$node->id})){
				$node->branch_length($$brlen_tbl{$node->id});
				}
			else{
				$node->branch_length(0); 
				}
			}
		elsif($node->branch_length == 0){ next; }	#skipping phyla
		else{
			$node->branch_length($$brlen_tbl{$node->id});
			}
		}
		#exit;
	return $tree;
	}

sub check_brlen{
	# checking branch lengths of parent-daughter #
	my $tree = shift;
	print "\nBranch lengths\n";
	foreach my $node ($tree->get_leaf_nodes){
		print "leaf: ", $node->id, ": ", $node->branch_length, "\n";
		}
	foreach my $node ($tree->get_nodes){
		if($node->is_Leaf == 1){next;}
		else{ print "int: ", $node->id, ": ", $node->branch_length, "\n"; }
		}

	}

sub tree_write{
	my ($tree, $outfile) = @_;
	my $out = new Bio::TreeIO(-file => "> $outfile",
							  -format => "newick");
	$out->write_tree($tree);
	}

sub error_routine{
	my $error = $_[0];
	my $exitcode = $_[1];
	print STDERR "ERROR: $error\nSee help: [-h]\n";
	exit($exitcode);
	}

sub usage {
 my $usage = <<HERE;
Usage:
  brlen_mod.pl -tree -brlen [-format] [-o]
  * used with output from cor_network_nest_mthr.R
Options:
  -tree		Tree file name. (newick or nexus)
  -brlen	Correlation file used for branch length values
  -format	Tree format ('newick' or 'nexus')
Description:
  Program takes a tree file and correlation file created by
  cor_network_nest_mthr.R and applies the correlation values
  to the branch lengths.
Requires:
  1) Bioperl
Notes:
	Version: $version
	Last Modified: $mod
	Author: $author
Categories:
	Phylogeny
		
HERE
	print $usage;
    exit(1);
}
