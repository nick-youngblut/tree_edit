#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $org_in, $runID);
my $format = "newick";
my $options = "-a";	
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "organism=s" => \$org_in,
	   "runID=s" => \$runID,
	   "x=s" => \$options,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if($tree_in){
	die " ERROR: $tree_in file not found!\n" unless -e $tree_in;
	$format = format_check($format);
	}
die " ERROR: $org_in not found!\n" if $org_in && ! -e $org_in;
die " ERROR: provide a runID!\n" unless $runID;


### MAIN
# loading leaf labels #
my $labels_r;
if($tree_in){ $labels_r = get_leaves_from_tree($tree_in, $format); }
else{ $labels_r = get_leaves_from_stdin(); }

# converting leaf names #
my $rev_org_r = convert_names($labels_r, $org_in) if $org_in;

# getting COG categories #
my %COG; my $cnt = 0;
for my $name (@$labels_r){
	print STDERR "...getting COGs for $name\n" unless $verbose;
	get_COG_cat($name, $runID, $options, \%COG);
	
	$cnt++; last if $cnt > 2;
	}

# writing out metadata table #
write_metadata(\%COG, $labels_r, $rev_org_r);


### Subroutines
sub write_metadata{
	my ($COG_r, $labels_r, $rev_org_r) = @_;
	
	# getting all COG categories #
	my %COG_cat;
	foreach my $org (keys %$COG_r){
		foreach my $cat (keys %{$COG_r->{$org}}){
			$COG_cat{$cat} = 1;
			}
		}
	my @COG_cat = keys %COG_cat;
	
	# getting hexideximal colors for the categories #
	my @hex = get_hex_colors($#COG_cat);
	
	# writing out header #
	print join("\t", "LABELS", @COG_cat), "\n";
	print join("\t", "COLORS", @hex), "\n";
	
	# writing out body #
	foreach my $org (keys %$COG_r){
		my @line;
		foreach my $cat (@COG_cat){
			if(exists $COG_r->{$org}{$cat}){
				push(@line, $COG_r->{$org}{$cat});
				}
			else{ push(@line, "NA"); }
			}
		# changing organisms names back #
		if($rev_org_r){
			die " LOGIC ERROR: $!\n" unless exists $rev_org_r->{$org};
			print join("\t", $rev_org_r->{$org}, @line), "\n"; 
			}
		else{
			print join("\t", $org, @line), "\n";
			}
		}
	}



sub get_hex_colors{
# getting hexidecimal colors for metadata #
	my $n_col = shift;
	my @hex = qw/FF0000 FF6600 33FF00 0000FF FF00FF FF0099 33CCFF 990000 CC0099 000066 006600 CC6600/;
	map{$_ =~ s/^/#/} @hex;
	
	if($n_col > scalar @hex){
		for my $i (0..int( $n_col/ scalar @hex)){
			push @hex, @hex;
			}
		}
	
	return @hex[0..$n_col];
	}

sub get_COG_cat{
# getting the COG categories of all the of genes specified #
	my ($name, $runID, $options, $COG_r) = @_;
	
	my $cmd = "echo \"$name\" | db_findClustersByOrganismList.py $options $runID | db_getGenesInClusters.py | db_getExternalClusterGroups.py -d cog -g 3 | db_getExternalClustersById.py -c 13 | ";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, $cmd or die $!;
	
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		next unless /COG\d+/;
		
		s/^.+\[|\].*$|\/.+//g;		# removing all but COG category
		s/[-, ;:]+$//g;				# removing junk from end 
		s/[-, ;:]+/_/g;				# sanitizing
		
		$COG_r->{$name}{$_}++;
		}
	close PIPE;
	
		#print Dumper %$COG_r; exit;
	}

sub convert_names{
# converting names to FIG ID if organism file provided #
	my ($labels_r, $org_in) = @_;
	
	# making organism hash #
	my %org;
	my %rev_org;
	open IN, $org_in or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		(my $sanitized = $line[0]) =~ tr/-. /_/;
		$org{$sanitized} = $line[0];
		$rev_org{$line[0]} = $sanitized;
		}
	close IN;
		#print Dumper %org; exit;
		
	# changing names #
	foreach my $name (@$labels_r){
		die " ERROR: $name not found in organism file!\n"
			unless exists $org{$name};
		$name = $org{$name};
		}
		
		#print Dumper @$labels_r; exit;
	return \%rev_org;
	}

sub get_leaves_from_stdin{
# loading leaves from stdin #
	print STDERR " No tree provided; getting leaf names from STDIN\n";
	my @labels;
	while(<>){
		chomp;
		next if /^\s*$/;
		push @labels, $_;
		}
	return \@labels;
	}

sub get_leaves_from_tree{
# loading leaves from tree file #
	use Bio::TreeIO;
	my ($tree_in, $format) = @_;
	my $treeio = Bio::TreeIO -> new(-file => $tree_in,
							-format => $format);	
	my @labels;
	while (my $tree = $treeio->next_tree){
		for my $node ($tree->get_leaf_nodes){
			push @labels, $node->id;
			}
		}
		#print Dumper @labels; exit;
	return \@labels;
	}

sub format_check{
	my $format = shift;
	if ($format =~ /nwk|newick/i){ $format = "newick" if $format =~ /nwk|newick/i; }
	elsif ($format =~ /nex|nexus|tre/i){ $format = "nexus"; }
	else{ die " ERROR: format not recognized [newick or nexus accepted]\n"; }
	return $format;
	}

__END__

=pod

=head1 NAME

tree_makeCOGMetadata.pl -- COG category info to ITOL metadata table

=head1 SYNOPSIS

=head2 Input: tree file

tree_makeCOGMetadata.pl [options] -t tree.nwk -r

=head2 Input: leaf names from STDIN

nw_labels -I tree.nwk | tree_makeCOGMetadata.pl -r 

=head2 Required flags (depending on input)

=over

=item -t 	Tree file (newick or nexus).

=item -r 	ITEP cluster runID

=back

=head2 Options

=over

=item -f 	Tree format (only needed if '-t' provided; newick or nexus). [newick] 

=item -o 	ITEP organism file (for associating tree names to ITEP names).

=item -x 	Inclusiveness of clusters; flag(s) for db_findClustersByOrganismList.py. [-a]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_makeCOGMetadata.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

tree_makeCOGMetadata.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

tree_makeCOGMetadata.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

