#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $tformat, $count_in, $meta_in, $regex, $color_in, $brlen_write);
my ($use_rep);
my $brlen_cut = 0;
GetOptions(
	   "tree=s" => \$tree_in,			# tree file
	   "format=s" => \$tformat,			# tree format
	   "count=s" => \$count_in,			# count table file
	   "color=s" => \$color_in, 		# color range file
	   "regex=s" => \$regex, 
	   "meta=s" => \$meta_in,
	   "length=f" => \$brlen_cut, 
	   "branch" => \$brlen_write, 		# write out all brlens
	   "x" => \$use_rep, 				# use a representative taxon for name
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " Provide a tree file (newick or nexus).\n" if ! $tree_in;
$tformat = check_tree_format($tformat);
$regex = qr/$regex/ if $regex;

### MAIN
# loading metadata #
my $count_r = load_count($count_in) if $count_in;
my $meta_r = load_metadata($meta_in) if $meta_in;
my $color_r = load_color_range($color_in) if $color_in;

# loading tree #
my $treeo = tree_io($tree_in, $tformat);

# writing br-lengths #
write_brlen($treeo) if $brlen_write;

# collapsing tips #
my $brlen_r = make_node_brlen_index($treeo);
collapse_tips_brlen($treeo, $brlen_r, $brlen_cut, $count_r, $meta_r, $regex, $color_r);

# writing modified tree #
tree_write($treeo, $tree_in);

# writing metadata #
write_metadata($meta_r,$meta_in) if $meta_in;
write_color_range($color_r, $color_in) if $color_in;
write_count($count_r, $count_in) if $count_in;


### Subroutines
sub write_count{
# writing metadata #
	my ($count_r, $count_in) = @_;
	
	(my $count_out = $count_in) =~ s/\.[^\.]+$|$/_br-col.txt/;
	open OUT, ">$count_out" or die $!;
	print OUT join("\t", @{$count_r->{"HEADER"}{"HEADER"}}), "\n";
	foreach my $taxon (sort keys %{$count_r->{"BODY"}} ){
		print OUT join("\t", $taxon, @{$count_r->{"BODY"}{$taxon}}), "\n";
		}
	close OUT;

	print STDERR " Count file written: '$count_out'\n";
	}

sub write_color_range{
# writting color range #
	my ($color_r, $color_in) = @_;
	
	(my $color_out = $color_in) =~ s/\.[^\.]+$|$/_br-col.txt/;
	open OUT, ">$color_out" or die $!;
	
	foreach my $taxon (sort keys %$color_r){ 
		print OUT join("\t", $taxon, @{$color_r->{$taxon}}), "\n";
		}
	close OUT;
	
	print STDERR " Color range file written: '$color_out'\n";
	}

sub write_metadata{
# writing metadata #
	my ($meta_r, $meta_in) = @_;
	
	(my $meta_out = $meta_in) =~ s/\.[^\.]+$|$/_br-col.txt/;
	open OUT, ">$meta_out" or die $!;
	print OUT join("\t", @{$meta_r->{"LABELS"}}), "\n";
	print OUT join("\t", @{$meta_r->{"COLORS"}}), "\n";
	foreach my $taxon (sort keys %{$meta_r->{"BODY"}} ){
		print OUT join("\t", $taxon, @{$meta_r->{"BODY"}{$taxon}}), "\n";
		}
	close OUT;

	print STDERR " Metadata file written: '$meta_out'\n";
	}

sub tree_write{
	### writting out a nexus tree file ###
	my ($treeo, $tree_in) = @_;

	(my $outfile = $tree_in) =~ s/\.[^\.]+$|$/_br-col.nwk/;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	print STDERR " Newick tree file written: '$outfile'\n";

	}

sub write_brlen{
# writing brlens of all nodes #
	my ($treeo) = @_;
	
	#for my $node ($treeo->get_nodes){
	#	next if $node->is_Leaf;
	#	$node->id($node->height);
	#	my @brlens;
	#	for my $child ($node->get_all_Descendents){
	#		push @brlens, $treeo->distance(-nodes => [$node, $child]);
	#		}
	#	$node->id(max(@brlens));
	#	}
	#tree_write($treeo, $tree_in); exit;
	
	
	my $node_cnt = 0;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		
		my @brlens;
		for my $child ($node->get_all_Descendents){
			push @brlens, $treeo->distance(-nodes => [$node, $child]);
			}
		print max(@brlens), "\n";
		}

	exit;
	}

sub make_node_brlen_index{
	my ($treeo) = @_;
	
	my %brlen;
	my $node_cnt = 0;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		# labeling internal nodes w/ unique ID #
		$node_cnt++;
		$node->description($node_cnt);
		next unless $node->ancestor;
		
		# getting distance from #
		my @brlens;
		for my $child ($node->get_all_Descendents){
			push @brlens, $treeo->distance(-nodes => [$node, $child]);
			}
		$brlen{$node->description} = max @brlens; #$treeo->distance(-nodes => [$node->ancestor, $node]);
		}

		#print Dumper %brlen; exit;
	return \%brlen;
	}

sub collapse_tips_brlen{
# collapsing tips based on branch length #
	my ($treeo, $brlen_r, $brlen_cut, $count_r, $meta_r, $regex, $color_r) = @_;
		
	# removing nodes from < height to > height #
	my $col_cnt = 0;
	foreach my $node (sort{$brlen_r->{$a}<=>$brlen_r->{$b}} keys %$brlen_r){
		
		# removing taxa w/ height < brlen cutoff #
		next if $brlen_r->{$node} >= $brlen_cut;
			
		# getting internal node #
		my @nodes = $treeo->find_node(-description => $node);
		die " LOGIC ERROR: no node found! $!\n" unless @nodes;
		die " LOGIC ERROR: >1 node found! $!\n" if scalar @nodes > 1;
		
		# removing all leaves 1st #
		#my @rm_list;
		my %rm_list;
		for my $child ( $nodes[0]->get_all_Descendents ){
			next unless $child->is_Leaf;
			next if $child->is_Leaf && $regex && $child->id =~ /$regex/;
					
			#push(@rm_list, $child->id) if $child->is_Leaf;
			$rm_list{$child->id} = $treeo->distance(-nodes => [$nodes[0],$child])
				if $child->is_Leaf;
			$treeo->remove_Node($child);
			}
		
		# next, removing all internal nodes < cutoff #
		for my $child ( $nodes[0]->get_all_Descendents ){
			next if $child->is_Leaf && $regex && $child->id =~ /$regex/;
			
			#push(@rm_list, $child->id) if $child->is_Leaf;
			#$rm_list{$child->id} = $treeo->distance(-nodes => [$nodes[0],$child])
			#	if $child->is_Leaf;		
			$treeo->remove_Node($child);
			}	
			
		# new name for node #
		next unless keys %rm_list;		# only if something was removed
		
		
		# adding collapsed node or changing current node #
		$col_cnt++;
		my $col_node_name;
		if($use_rep){
			foreach my $child_id (sort{$rm_list{$b}<=>$rm_list{$a}} keys %rm_list){
				#print Dumper $child_id;
				$col_node_name = $child_id;
				last;			
				}
			}
		else{
			$col_node_name = join("_", $col_cnt, scalar keys %rm_list);
			}
		
		my $new_node = new Bio::Tree::Node(-id => $col_node_name, 
								-branch_length => $brlen_r->{$node});			
		$nodes[0]->add_Descendent($new_node);
	
		# editing metadata #
		sum_count($count_r, [keys %rm_list], $col_node_name) if $count_r;
		sum_meta($meta_r, [keys %rm_list], $col_node_name) if $meta_r;
		sum_color_range($color_r, [keys %rm_list], $col_node_name) if $color_r;
		}
	
	print STDERR " Number of nodes collapsed: $col_cnt\n";
	}

sub collapse_tips_brlen_OLD{
# collapsing tips based on branch length #
	my ($treeo, $brlen_cut, $count_r, $meta_r, $regex, $color_r) = @_;
	
	my $col_cnt = 1;		# number of collapsed clades
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;					# no collapsign of leaves directly
		
		# collapsing node #
		my @rm_list;
		if($node->height < $brlen_cut){
			my $rm_num = 0;
			# removing leaf nodes 1st #
			for my $child ( $node->each_Descendent ){
				next unless $child->is_Leaf;
				next if $regex && $child->id =~ /$regex/;		# skipping if match
				
				push(@rm_list, $child->id);
				
					#print $child->id, "\n";
				$treeo->remove_Node($child);	
				$rm_num++;
				}
				
			# removing internal nodes #
			if ($rm_num){
				for my $child ( $node->each_Descendent){
					next if $child->is_Leaf && $child->id =~ /$regex/;
					
					push(@rm_list, $child->id) if $child->is_Leaf;
					$treeo->remove_Node($child);
					}
				}
			
			# new name for node #
			next unless @rm_list;		# only if something was removed
			
			# adding collapsed node or changing current node #
			$col_cnt++;
			my $col_node_name = "$col_cnt\_$rm_num";
			if($node->is_Leaf){
				my $new_node = new Bio::Tree::Node(-id => $col_node_name, 
											-branch_length => $node->height);
				$node->add_Descendent($new_node);
				}
			else{
				$node->id($col_node_name);
				}
		
			# editing metadata #
			sum_count($count_r, \@rm_list, $col_node_name) if $count_r;
			sum_meta($meta_r, \@rm_list, $col_node_name) if $meta_r;
			sum_color_range($color_r, \@rm_list, $col_node_name) if $color_r;
			}
		}
	}

sub sum_count{
# summing values in count file for collapsed clades #
	my ($count_r, $rm_list_r, $col_node_name) = @_;
	
	# summing #
	my @sums;
	foreach my $taxon (@$rm_list_r){
		die " ERROR: $taxon not found in the count file!\n"
			unless exists $count_r->{"BODY"}{$taxon};
		for my $i (0..$#{$count_r->{"BODY"}{$taxon}}){
			$sums[$i] += ${$count_r->{"BODY"}{$taxon}}[$i];
			}
		delete $count_r->{"BODY"}{$taxon};
		}
		
	$count_r->{"BODY"}{$col_node_name} = \@sums;
	}

sub sum_color_range{
# pruning color range; replacing with random line from deleted #
	my ($color_r, $rm_list_r, $col_node_name) = @_;
	
	my $replace;
	foreach my $taxon (@$rm_list_r){
		die " ERROR: $taxon not in color range file!\n"
			unless exists $color_r->{$taxon};
		
		$replace = $color_r->{$taxon};
		delete $color_r->{$taxon};
		}

	$color_r->{$col_node_name} = $replace;
	}

sub sum_meta{
# summning for collapsed clades #
	my ($meta_r, $rm_list_r, $col_node_name) = @_;
	
	# summing #
	my @sums;
	foreach my $taxon (@$rm_list_r){
		die " ERROR: $taxon not found in the metadata file!\n"
			unless exists $meta_r->{"BODY"}{$taxon};
		for my $i (0..$#{$meta_r->{"BODY"}{$taxon}}){
			$sums[$i] += ${$meta_r->{"BODY"}{$taxon}}[$i];
			}
		delete $meta_r->{"BODY"}{$taxon};
		}
		
	$meta_r->{"BODY"}{$col_node_name} = \@sums;
	
		#print Dumper $meta_r->{"BODY"}; 
	}

sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO->new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		
	return $treeio;
	}

sub load_color_range{
	my $color_in = shift;
	
	open IN, $color_in or die $!;

	my %color;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		$color{$line[0]} = [@line[1..$#line]];
		}
		#print Dumper %color; exit;
	return \%color;
	}

sub load_metadata{
# loading metadata table (itol format) #
	my ($meta_in ) = @_;
	
	open IN, $meta_in or die $!;
	my %meta;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		if($. == 1){
			die " ERROR: 1st line of metadata table should start with 'LABELS'\n"
				unless /^LABELS/;
			$meta{"LABELS"} = \@line;
			}
		elsif($. == 2){
			die " ERROR: 2nd line of metadata table should start with 'COLORS'\n"
				unless /^COLORS/;
			$meta{"COLORS"} = \@line;
			}
		else{
			$meta{"BODY"}{$line[0]} = [@line[1..$#line]];
			}
		}
	close IN;

		#print Dumper %meta; exit;
	return \%meta;
	}

sub load_count{
# loading count file #
	my ($count_in) = @_;
	open IN, $count_in or die $!;
	
	my %count;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		die " ERROR: the count file must be at least 2 columns (rownames, count)\n"
			if scalar @line < 2;
			
		# header #
		if($. == 1){ $count{"HEADER"}{"HEADER"} = \@line; }
		else{ $count{"BODY"}{$line[0]} = [@line[1..$#line]]; }
		}
	close IN;
			
		#print Dumper %count; exit;
	return \%count;			# returning taxa for pruning
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

tree_collapse_tips.pl -- collapsing tips by a branch length cutoff

=head1 SYNOPSIS

tree_collapse_tips.pl [flags]

=head2 require flags

=over

=item -tree  <char>

Tree file (newick or nexus).

=back

=head2 optional flags

=over

=item -format  <char>

Tree file format (newick or nexus). [newick]

=item -length  <float>

Branch length cutoff. [< 0]

=item -x  <bool>

Use a repesentative sequnece for collapsed clades? [FALSE]

=item -count  <char>

Count file in Mothur format

=item -meta  <char>

Metadata file in ITOL format

=item -color  <char>

Color file in ITOL format

=item -regex  <char>

Regular expression for excluding taxa from being collapsed.

=item -branch  <bool>

Write out all branch lengths (internal nodes to most distant decendent).

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_collapse_tips.pl

=head1 DESCRIPTION

Collapse branches in a tree that have a branch length of < (-length)
from the ancestral node.

=head2 Collapsed node labeling 

=head3 Default

Collapsed nodes are labeled as: "collapsed-nodeID"_"number_taxa_collapsed"

=head3 '-x'

Leaf with the greatest branch length from the collapsed node will be used
as the representative taxon.

=head2 Metadata

If any metadata files are provided (count, ITOL-metadata, color), 
the taxon labels are updated and the abundances are summed. A random 
taxon in a collapsed clade will be used for the color info if a
color file is provided.

Output file names are based on input file names ("*br-col*")

=head2 WARNING

'-count' flag not fully tested!

=head1 EXAMPLES

=head2 Getting branch length distribution

tree_collapse_tips.pl -t tree.nwk -branch

=head2 Collapse tree and use repesentative taxon names

tree_collapse_tips.pl -t tree.nwk -l 0.001 -x

=head2 Collapse tree with a metadata file

tree_collapse_tips.pl -t tree.nwk -meta meta.txt -l 0.001

=head2 Collapse tree with a count file

tree_collapse_tips.pl -t tree.nwk -count count.txt -l 0.001

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

