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

my ($verbose, $tree_in, $tformat, $count_in, $meta_in, $regex, $color_in);
my $brlen_cut = 0;
GetOptions(
	   "tree=s" => \$tree_in,			# tree file
	   "format=s" => \$tformat,			# tree format
	   "count=s" => \$count_in,			# count table file
	   "color=s" => \$color_in, 		# color range file
	   "regex=s" => \$regex, 
	   "meta=s" => \$meta_in,
	   "length=f" => \$brlen_cut, 
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

# collapsing tips #
#root2tips($treeo, $treeo->get_root_node(), $brlen_cut, $count_r, $meta_r);
collapse_tips_brlen($treeo, $brlen_cut, $count_r, $meta_r, $regex, $color_r);

# writing modified tree #
tree_write($treeo, $tree_in);

# writing metadata #
write_metadata($meta_r,$meta_in) if $meta_in;
write_color_range($color_r, $color_in) if $color_in;

### Subroutines
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

sub root2tips{
# moving from the root to the tips #
	my $treeo = shift;
	my $anc = shift;


	# starting w/ root #
	#$anc = $tree->get_root_node() unless $anc;
	
	# I/O #
	die " ERROR: no ancestor!\n" unless $anc;
	

	# if max length to tip < brlen collapse #
	if($anc->height < $brlen_cut){
		print Dumper "here"; exit;
		collapse_tips_brlen($treeo, $anc, @_);
		}
	else{ 	# else: next node #
		for	my $child ($anc->each_Descendent){
			next if $child->is_Leaf;
			root2tips($treeo, $anc, @_);
			}
		}
	}

sub collapse_tips_brlen{
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
				next if $child->id =~ /$regex/;		# skipping if match
				
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
			#sum_count($count_r, \@rm_list, $col_node_name) if $count_r;
			sum_meta($meta_r, \@rm_list, $col_node_name) if $meta_r;
			sum_color_range($color_r, \@rm_list, $col_node_name) if $color_r;
			}
		}
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

sub collapse_tips_brlen_OLD{
# collapsing tips based on branch length #
	my ($treeo, $node, $brlen_cut, $count_r, $meta_r) = @_;
			
		#print Dumper " here"; exit;	
			
	for my $child ( $node->each_Descendent ){
		my $brlen = $treeo->distance(-nodes => [$node, $child]);
		#print Dumper $brlen; 
		if( ! $brlen || $brlen < $brlen_cut){
			$treeo->remove_Node($child);			
			}
		}
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

