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

my ($verbose, $tree_in, $tformat, $count_in, $outfile);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$tformat,
	   "count=s" => \$count_in,
	   "outfile=s" => \$outfile,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: Provide a tree file (newick or nexus).\n" if ! $tree_in;
die " ERROR: Provide a count file (mothur format).\n" if ! $count_in;
$tformat = check_tree_format($tformat);

### MAIN
# loading files #
my $treeo = tree_io($tree_in, $tformat);
my $count_r = load_count($count_in);

# moving node_ids to bootstrap #
id_2_bootstrap($treeo);

# finding lca and annotating #
$treeo = find_lca_annotate($treeo, $count_r);
my $tree_out = tree_write($treeo, $outfile, $tree_in);

# editing tree file (convert to nexus) #
tree_attrib_edit($treeo, $tree_out);

# finding lca taxa and making a color definition file #
my $lca_taxa_r = find_lca_taxa($treeo, $count_r);
write_color_definitions_file($lca_taxa_r, $outfile, $tree_in);



### Subroutines
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
	$outfile =~ s/\.nwk$|$/.tree/;
	
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

sub id_2_bootstrap{
# moving internal node ids to bootstrap #
	my $treeo = shift;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		next if ! $node->id;
		if($node->id =~ /\d{1,2}/){
			$node->bootstrap($node->id);
			$node->id("");
			}
		}

	return $treeo;
	}

sub write_color_definitions_file{
# writing out color definitions file for itol #
	my ($lca_taxa_r, $outfile, $tree_in) = @_;
	
	($outfile = $tree_in) =~ s/\.[^\.]+$|$/_col-dev.txt/;
	open OUT, ">$outfile" or die $!;
	
	print OUT join("\t", qw/NODE_ID TYPE COLOR LABEL/), "\n";
	
	my @colors = ("#0000ff","#00ff00","#ff0000","#aaffaa","#ff00ff","#CCCCCC","#66ffff");
	my $col_cnt = 0;
	foreach my $samp (keys %$lca_taxa_r){
		print OUT join("\t", join("|", @{$$lca_taxa_r{$samp}}), "clade", $colors[$col_cnt], $samp), "\n";
		
		$col_cnt++;
		$col_cnt = 0 if $col_cnt > 6;		# cycling through  colors
		}
	
	close OUT;
	print STDERR " Color definition file written:\n  $outfile\n";
	}

sub tree_write{
	### writting out a nexus tree file ###
	my ($treeo, $outfile, $tree_in) = @_;

	($outfile = $tree_in) =~ s/\.[^\.]+$|$/_lca.nwk/;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	print STDERR " Newick tree file written:\n  $outfile\n";

	return $outfile;
	}

sub find_lca_annotate{
# finding the lca for each group and annotationing the lca nodes #
	my ($treeo, $count_r) = @_;
	
	print STDERR " Annotating tree\n";
	
	foreach my $samp (keys %$count_r){
		my @samp_nodes;
		foreach my $node (@{$$count_r{$samp}}){
			my @found_nodes = $treeo->find_node(-id => $node);
			if (! @found_nodes){
				print STDERR " WARNING: $node not found in tree file!\n";
				next;
				}
			else{
				push @samp_nodes, $found_nodes[0];
				}
			}
		
		die " ERROR: no nodes found for $samp\n" if ! @samp_nodes;

		my $lca = $treeo->get_lca(-nodes => \@samp_nodes);
		
		if( ! $lca->id){
			$lca->id( join("", " [&sample={'", $samp, "'}]") );
			}
		else{ 
			(my $lca_id = $lca->id) =~ s/\}]/,/;
			$lca->id(join("", $lca_id, "'", $samp, "'}]") ); 
			}
		}
	return $treeo;
	}


sub find_lca_taxa{
# finding lca for each sample #
	my ($treeo, $count_r) = @_;
	
	
	my %lca_taxa;
	foreach my $samp ( keys %$count_r){
		# status #
		print STDERR "Finding lca for: $samp\n";
		
		# finding all of the tip nodes for a sample #
		my @samp_nodes;
		foreach my $taxon ( @{$$count_r{$samp}} ){
			my @found_nodes = $treeo->find_node(-id => $taxon);
			print STDERR " WARNING: $taxon not found in tree file!\n" if ! @found_nodes;
			die " ERROR: $taxon hit multiple taxa in tree file!\n" if scalar @found_nodes > 1;
			push(@samp_nodes, $found_nodes[0]);
			}
		
		
		# finding the lca for all of the found nodes #
		my $lca = $treeo->get_lca(-nodes => \@samp_nodes);
		
		# finding most divergent taxa for lca #
		my @taxa_rm;
		while(1){
			#last if scalar @samp_nodes == 2; 
			#my $rm_len = int (scalar @samp_nodes / 10); 	# removing 10% at a time
			my %rands;
			while(1){
				my $rand = int rand( scalar @samp_nodes); 
				$rands{$rand} = 1;
				last if scalar keys %rands == scalar @samp_nodes - 1
				}
				
			@taxa_rm = @samp_nodes[keys %rands];
			my $rm_lca = $treeo->get_lca(-nodes => \@taxa_rm);
			
			if($lca == $rm_lca){	# if 
				@samp_nodes = @taxa_rm; 				# reducing taxa
				last if scalar @samp_nodes == 2;		# if down to 2 taxa
				}
			else{
				next;									# try again
				}
			}
		
		# loading sample node ids #
		my @samp_node_ids;
		map{ push(@samp_node_ids, $_->id) } @samp_nodes;
		
		$lca_taxa{$samp} = \@samp_node_ids;
		
		}
		#print Dumper %lca_taxa; exit;
	return \%lca_taxa;
	}

sub load_count{
# loading count file (mothur format) #
	my $count_in = shift;
	open IN, $count_in or die $!;
	
	my (%count, @header);
	while(<IN>){
		chomp;
		my @line = split /\t/;
		if($.==1){	# header
			@header = @line;
			next;
			}
		else{
			for my $i (2..$#line){		# load sample-specific taxa
				push (@{$count{$header[$i]}}, $line[0]) if $line[$i] > 0;
				}
			
			}
		}
	close IN;
		#print Dumper %count; exit;
	return \%count;
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

lca_finder.pl -- Find LCA for a group of taxa and make a color definition file for iTOL

=head1 SYNOPSIS

lca_finder.pl -t -c [-f] [-o]

=head2 options

=over

=item -t 	Tree file (newick or nexus format)

=item -f 	Tree file format. [newick]

=item -c 	Count file (Mothur format).

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc lca_finder.pl

=head1 DESCRIPTION

Finding last common ancestor for the taxa in a sample.

A nexus tree file contains the annotated notes ('sample').

The 'col-dev.txt' file contains a color definition file for iTOL

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

