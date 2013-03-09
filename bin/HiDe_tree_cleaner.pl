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

my ($verbose, @gene_in, $species_in, $tformat, $species, $multi, $index_in, $HiDe_in);
GetOptions(
	   "gene=s{,}" => \@gene_in,		# gene trees
	   "format=s" => \$tformat,
	   "species=s" => \$species_in,		# species tree
	   "index=s" => \$index_in, 		# index file
	   "out=s" => \$HiDe_in, 			# hide input
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# if HiDe output provided #
HiDe_out_change_names($HiDe_in, $index_in) if $HiDe_in;

die " ERROR: Provide a species & >=1 gene tree file (newick or nexus).\n" 
	if ! @gene_in || ! $species_in;
$tformat = check_tree_format($tformat);


### MAIN
# loading species tree names #
my $names_r = check_names($species_in, \@gene_in, $tformat);

# making name index #
my $nameI_r = make_name_index($names_r, $species_in);


# editting and writing #
clean_tree($species_in, $tformat, $nameI_r, "species.newick");
foreach my $gene_in (@gene_in){
	clean_tree($gene_in, $tformat, $nameI_r);
	}


# main subroutines #
sub clean_tree{
	my ($gene_in, $tformat, $nameI_r, $outfile) = @_;
	
	# I/O #
	($outfile = $gene_in) =~ s/\.[^\.]+$|$/_cln.newick/ unless $outfile;
	unlink $outfile if -e $outfile;
	
	my $treeio = Bio::TreeIO -> new(-file => $gene_in,
								-format => $tformat);
	my $out = new Bio::TreeIO(-file => ">>$outfile", -format => "newick");
	
	# cleaning trees #
	my $tree_cnt = 0;
	while (my $treeo = $treeio->next_tree){
		$treeo = edit_leaf_labels($treeo, $nameI_r);
		$treeo = remove_node_labels($treeo);
		$treeo = remove_bootstrap($treeo);
		$treeo = remove_branch_length($treeo);
		$treeo = check_tree_cleaning($treeo) if $verbose;
		$out->write_tree($treeo);
		$tree_cnt++;
		}
	
	# removing all colons (from branch lengths #
	`perl -pi -e 's/://g; s/;/;\n/g' $outfile`;	
	
	# calling HiDe_clean_trees.r #
	my $cmd = "HiDe_clean_trees.r -t $outfile";
	$cmd .= " -s" if $outfile eq "species.newick";		# species flag if species
	$cmd .= " -m" if $tree_cnt == 1;						# turning on multi-tree flag in R script
	print STDERR "CMD: $cmd\n";
	`$cmd`;
	
	print STDERR "...Cleaned tree written: $outfile\n";
	}

sub HiDe_out_change_names{
	my ($HiDe_in, $index_in) = @_;
	die " ERROR: provide the HiDe output and the index file (*NI.txt).\n"
		if ! $HiDe_in || ! $index_in;
	
	my $nameI_r = load_name_index($index_in);
	
	(my $outfile = $HiDe_in) =~ s/\.[^\.]+$|$/_rn.txt/;
	open IN, $HiDe_in or die $!;
	open OUT, ">$outfile" or die $!;
	while(<IN>){
		chomp;
		my @line = split /\t| --> /;
		foreach my $cat (@line[1..2]){
			if($cat =~ /LCA/){					# if LCA of 2 taxa
				my @parts = split /[(,)]/, $cat;
				foreach my $index (@parts[1..2]){
					die " ERROR: '$index' not found in name index!\n" 
						unless exists $nameI_r->{$index};
					$index = $nameI_r->{$index}
					}
				$cat = "$parts[0]($parts[1],$parts[2])";
				}
			else{					# if taxon
				my @parts = split /-/, $cat;
				die " ERROR: '$parts[1]' not found in name index!\n" 
						unless exists $nameI_r->{$parts[1]};
				$parts[1] = $nameI_r->{$parts[1]};
				$cat = "$parts[0]-$parts[1]";
				}
			}
		if($line[3]){
			print OUT join("\t", $line[0], join(" --> ", @line[1..2]), $line[3]), "\n";
			}
		else{
			print OUT join("\t", $line[0], join(" --> ", @line[1..2])), "\n";
			}
		}
	close IN;
	close OUT;
	
	
	print STDERR "...New HiDe file written: $outfile\n";
	exit;
	}

### Subroutines
# names #

sub load_name_index{
	my ($index_in) = @_;
	open IN, $index_in or die $!;
	my %nameI;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		$nameI{$line[0]} = $line[1];
		}
	close IN;
	return \%nameI;		# %{index}=>{name}
	}

sub make_name_index{
# making a name index for conversion of names post HiDe #
	my ($names_r, $species_in) = @_;

	# making index #
	my %nameI;
	for my $i (0..$#$names_r){
		$nameI{$$names_r[$i]} = $i;
		}

	# writing out index #
	(my $outfile = $species_in) =~ s/\.[^\.]+$|$/_NI.txt/;
	open OUT, ">$outfile" or die $!;
	foreach my $name (sort{$nameI{$a}<=>$nameI{$b}} keys %nameI){
		print OUT join("\t", $nameI{$name} + 1, $name), "\n";
		}
	close OUT;
	
	print STDERR "...name index file written: $outfile\n";
	
	return \%nameI;
	}

sub check_names{
	my ($species_in, $gene_in_r, $tformat) = @_;

	print STDERR "...checking for matching names in species and gene trees\n";
	my %names;

	# species tree #
	open IN, "tree_get_names.pl -t $species_in -f $tformat | " or die $!;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		$names{$line[0]}++;
		die " ERROR in species tree: $line[0] found twice\n" if $names{$line[0]} > 1;
		}
	close IN;
	
	# gene trees #
	foreach my $gene_in (@$gene_in_r){
		open IN, "tree_get_names.pl -t $gene_in -f $tformat | " or die $!;
		while(<IN>){
			chomp;
			my @line = split /\t/;
			die " ERROR in $gene_in: $line[0] not found in species tree\n" unless
				exists $names{$line[0]};
			}
		close IN;	
		}

	print STDERR "...names match!\n";
	
	return [keys %names];
	}

# tree editting #

sub root_tree{
# rooting tree if not rooted #
	my ($treeo) = @_;
	if(! $treeo->get_root_node->id){
		for my $node ($treeo->get_leaf_nodes){
			$treeo->reroot($node);			# re-rooting on first leaf node
			last;
			}
		}
	return $treeo; 
	}

sub remove_branch_length{
	my ($treeo) = @_;
	for my $node ($treeo->get_nodes){
		$node->branch_length("");
		}
	return $treeo;
	}

sub edit_leaf_labels{
	my ($treeo, $nameI_r) = @_;
	for my $node ($treeo->get_nodes){
		next unless $node->is_Leaf;
		my $node_id = $node->id;
		die " ERROR: $node_id not found in name index!\n" unless
			exists $nameI_r->{$node_id};
		$node->id($nameI_r->{$node_id} + 1);			# converting IDs to numbers
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


__END__

=pod

=head1 NAME

HiDe_tree_cleaner.pl -- Cleaning trees for HiDe analysis (score.lua)

=head1 SYNOPSIS

=head2 Pre-HiDe

HiDe_tree_cleaner.pl -s -g [-f] [-v]

=head2 Post-HiDe

HiDe_tree_cleaner.pl -i -o

=head2 options

=over

=item -s 	Species tree (newick or nexus).

=item -g 	>=1 gene/LCB tree (newick or nexus).

=item -f 	Format of trees (newick or nexus). [newick]

=item -i 	Index file.

=item -o 	score.lua output.

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc HiDe_tree_cleaner.pl

=head1 DESCRIPTION

HiDe (eg score.lua) requires very simple trees for input.
This script removes:

=over

=item * 	Branch lengths

=item * 	Internal node labels

=item * 	Bootstrap values

=back

It also coververts all leaf names to numbers. The index is saved to "*_NI.txt"

Leaf names must match among species and gene trees!

=head2 Requires:

=over

=item * 	tree_get_names.pl

=back

=head1 EXAMPLES

=head2 Usage: pre-HiDe

HiDe_tree_cleaner.pl -s species.nwk -g gene*.newick

=head2 Usage: HiDe

score.lua . 80 > score.lua_out.txt

=head2 Usage: post-HiDe

HiDe_tree_cleaner.pl -o score.lua_out.txt -i species_NI.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

