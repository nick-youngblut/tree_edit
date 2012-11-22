#!/usr/bin/perl
my $mod = "8/1/12 3:47 PM";
my $version = "0.2.1";
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
use File::Spec;
use Bio::TreeIO;
use List::Util qw/maxstr/;

### global variables
my ($error);

### I/O
if ($#ARGV < 0){
	&usage;
	}
my ($tree_in, $tformat, $attrib_in, $att_header, $verbose, $outfile);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$tformat,
	   "attribute=s" => \$attrib_in,
	   "header" => \$att_header,
	   "verbose" => \$verbose,
	   "outfile=s" => \$outfile,
	   "help|?" => \&usage # Help
	   );

### Input error check
die " Provide a tree file (newick or nexus).\n" if ! $tree_in;
$tree_in = File::Spec->rel2abs($tree_in);
$tformat = check_tree_format($tformat);
$attrib_in = File::Spec->rel2abs($attrib_in) if $attrib_in;
($outfile = $tree_in) =~ s/\.[^\.]+$|$/_abnd.nwk/;

### Routing main subroutines
$attrib_in ? my $attrib_ref = load_attrib($attrib_in, $att_header) : die " Provide an attribute file (make_abund.pl & rep_pop_overlap.pl)\n";
my $treeo = tree_io($tree_in, $tformat);
$treeo = node_attribute($treeo, $attrib_ref);
	#check_additivity($treeo);
tree_write($treeo, $outfile);
tree_attrib_edit($treeo, $outfile);

#----------------------Subroutines----------------------#
sub tree_attrib_edit{
	### re-formatting newick output tree into nexus tree; 
	### taking attribute in names and moving to attribute portion
	my ($treeo, $outfile) = @_;
	
	my $nodes_ref = get_tip_names($treeo);	# getting all nodes
		
	# tree I/O #
	open (FILE, $outfile) or die $!;
	my @nex_tree = <FILE>;
	close FILE;
	
	# modding input tree external node labels #
	$nex_tree[0] =~ s/"(.*?)"/$1/g;
		
 	# modding tree file #
	$outfile =~ s/\.nwk/.tree/;
	
	open(FILE, "> $outfile") or die $!;
	print FILE "#NEXUS\nbegin taxa;\n\tdimensions ntax=", scalar(@$nodes_ref), ";\n\ttaxlabels\n";
	foreach(@$nodes_ref){ print FILE "\t$_\n"; }
	print FILE ";\nend;\n\nbegin trees;\n\ttree TREE1 = [&R] $nex_tree[0]\nend;\n\n";
	close FILE;
 	
 	print STDERR " Nexus tree file written:\n $outfile\n";
 	
	sub get_tip_names{
		my $tree = shift;
		my @nodes;
		for my $node ($tree -> get_nodes){
			if($node->is_Leaf){
				(my $tmp = $node->id) =~ s/\s*\[.+//; 
				push(@nodes, $tmp);
				}
			}
		return \@nodes;
		}
	}

sub tree_write{
	### writting out a newick tree file ###
	my ($treeo, $outfile) = @_;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	print STDERR " Newick tree file written:\n  $outfile\n";
	}

sub node_attribute{
	### adding abundances to tree nodes ###
	my ($treeo, $attrib_ref) = @_;
	
	# finding totals for columns for calc of proportions #
	my @totals;
	foreach(keys %$attrib_ref){ 
		next if $_ eq "header" || $_ eq "config";
		for (my $i=0; $i<=$#{$$attrib_ref{$_}}; $i++){
			if(${$$attrib_ref{"config"}}[$i+1] eq "int"){ $totals[$i] += ${$$attrib_ref{$_}}[$i];  }
			else{ $totals[$i] = 0; }
			}
		}
	
	# changing node label ID #
	for my $node ($treeo -> get_nodes){
		if($node->is_Leaf){ $node = label_tip_node($node, $attrib_ref, \@totals); }	# getting abundance for tip
		else{ $node = label_internal_node($node, $attrib_ref, \@totals); }	# adding up abundance for node
		}
	return $treeo;
	
	#---- node_abund subroutines ----#
	sub label_tip_node{
		### transforming attribance & adding attribance to node id ###
		my ($node, $attrib_ref, $total_ref) = @_;
		
		if(exists $$attrib_ref{$node->id}){ 	# if found in attribute table
			my $vals_ref = $$attrib_ref{$node->id};
			for(my $i=0; $i<=$#$vals_ref; $i++){ 
				if(${$$attrib_ref{"config"}}[$i+1] eq "int"){ $$vals_ref[$i] = $$vals_ref[$i]/$$total_ref[$i]; }	# richnes/abundance over total 
				elsif(${$$attrib_ref{"config"}}[$i+1] eq "group"){ 
					my $count = 1;
					$count++ while $$vals_ref[$i] =~ /,/g;
					$$vals_ref[$i] = $count;
					} # number of groups
				elsif(${$$attrib_ref{"config"}}[$i+1] eq "char"){ $$vals_ref[$i] = $$vals_ref[$i]; }				# no change
				else{ die $!; }
				$$vals_ref[$i] = ${$$attrib_ref{"header"}}[$i+1] . "=$$vals_ref[$i]";
				}
			$node->id(join("", $node->id, " [&", join(", ", @$vals_ref), "]") ); 
			} 
		else{ print " WARNING: ", $node->id, " not found in attribute table\n"; }
		return $node;
		}
	sub label_internal_node{
		### transforming attribance & adding attribance to node id ###
		my ($node, $attrib_ref, $total_ref) = @_;		
		my @attribs;
		for my $child ($node -> get_all_Descendents){	#summing attributes of all tips of internal node
			if($child->is_Leaf){ sum_children($$attrib_ref{$child->id}, \@attribs, $$attrib_ref{"config"});  }
			}
		
		total_attribs(\@attribs, $$attrib_ref{"config"}, $total_ref);	# dividing by total (if int); counting number of groups (group); using most prevalent (char);

		for (my $i=0; $i<=$#attribs; $i++){ $attribs[$i] = ${$$attrib_ref{"header"}}[$i+1] . "=$attribs[$i]"; }	# adding category info to attributes in tree
		
		if( ! $node->id){
			$node->id(join("", " [&", join(", ", @attribs), "]") );
			}
		else{ $node->id(join("", $node->id, " [&", join(", ", @attribs), "]") ); }
		return $node;
		
		sub sum_children{
			### summing children based on config
			my ($child_attrib_ref, $attrib_ref, $config_ref) = @_;
			for(my $i=0; $i<=$#$child_attrib_ref; $i++){ 
				if($$config_ref[$i+1] eq "int"){ $$attrib_ref[$i] += $$child_attrib_ref[$i]; }
				elsif($$config_ref[$i+1] =~ /group|char/){ push(@{$$attrib_ref[$i]}, split(/,/, $$child_attrib_ref[$i])); }
				else{ die $!; }
				}
			return $attrib_ref;
			}
		sub total_attribs{
			### dividing by total (if int); counting number of groups (group); using most prevalent (char) -> maxstr();
			my ($attribs_ref, $config_ref, $totals_ref) = @_;
			for(my $i=0; $i<=$#$attribs_ref; $i++){ 
				if($$config_ref[$i+1] eq "int"){
					$$attribs_ref[$i] = $$attribs_ref[$i]/$$totals_ref[$i]; 
					$$attribs_ref[$i] = 0.9999999999999 if $$attribs_ref[$i] eq "1"; ### VERY IMPORTANT for FigTREE to work!! # making root attribance < 1
					}
				elsif($$config_ref[$i+1] eq "group"){ 
					my @tmp = unique( @{$$attribs_ref[$i]}); 
					$$attribs_ref[$i] = scalar(@tmp); }	# counting number of groups
				elsif($$config_ref[$i+1] eq "char"){ 
					my %count; my $most_prev;
					foreach( @{$$attribs_ref[$i]}){ $count{$_}++; }
					$$attribs_ref[$i] = (sort{$count{$b} <=> $count{$a}} keys %count)[0];	# should give key w/ highest value
					}
				else{ die $!; }
				}
				#print Dumper($attribs_ref); #exit;
			return $attribs_ref;
			}
		}
	
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

sub load_attrib{
	###
	my ($attrib_in, $att_header) = @_;
	open IN, $attrib_in or die $!;
	my %attrib; my %col_format;	# attibute table and format of column placed in attribute table
	while(<IN>){
		chomp;
		$_ =~ s/#.+//;
		next if $_ =~ /^\s*$/;
		my @tmp = split(/\t/);
		die " Group file not formatted correctly\n" if scalar(@tmp) < 2;	# if not at least 2 column
		die " Sequence names in group file are not unique!\n" if exists($attrib{$tmp[0]});
		$attrib{"header"} = \@tmp, next if $.==1 && ! $att_header;
	
		for(my $i=0; $i<=$#tmp; $i++){ 	#figuring out column format
			if(exists $col_format{$i} && $col_format{$i} eq "group"){ next; }		# group
			elsif($tmp[$i] =~ /,/){ $col_format{$i} = "group"; }
			elsif(exists $col_format{$i} && $col_format{$i} eq "char"){ next; }
			elsif($tmp[$i]  =~ /[A-Za-z]/){  $col_format{$i} = "char"; }
			elsif($tmp[$i] =~ /^\d+\.*\d*$/){ $col_format{$i} = "int"; }
			else{ die " ERROR: can't figure out colume format!\n", $!; }
			}
		$attrib{shift @tmp} = \@tmp;
		}
	close IN;	
	$attrib{"config"} = [@col_format{ sort keys %col_format}]; 
		#print Dumper($attrib{"config"}); exit;
	return \%attrib;		# %@
	}

sub log10{
	my $n = shift;
	return log($n)/log(10);
	}

sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}

sub unique {
	# version: 1.0
	# usage: finding uniques in an array, returns an array of uniques
	my %utmp = map{$_, 1} @_;  #these lines pull out unique lake names and sort them
	@_ = keys %utmp;
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
 tree_add_attribue.pl -t -a [-f] [-h]
Options:
 -t 	Tree file.
 -f 	Tree file format (newick or nexus).
 -a 	Attribute file (names must match tip labels!)
 -h 	Header in attribute file? [default: TRUE]
Description:
 Program adds attributes to tree file based
 on a tab-delimited table of tip labels and attributes.
 Attribute formats:
 	integer, group (comma-delimited), character
 Attribute formatting:
 	integer: value/total
 	group: unique_groups_in_all_children
 	char: unique_groups_in_all_children (need to update)
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