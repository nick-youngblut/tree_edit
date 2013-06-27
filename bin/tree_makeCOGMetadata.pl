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
my $COG_cutoff = 0.01;			# COG category must be >= 1% of PEGs w/ COG
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "organism=s" => \$org_in,
	   "runID=s" => \$runID,
	   "x=s" => \$options,
	   "cutoff=f" => \$COG_cutoff,
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
die " ERROR: provide an organism file!\n" unless $org_in;


### MAIN
# loading leaf labels #
my $labels_r;
if($tree_in){ $labels_r = get_leaves_from_tree($tree_in, $format); }
else{ $labels_r = get_leaves_from_stdin(); }

# converting leaf names #
my $name_index_r = make_name_index($labels_r, $org_in);

# getting COG categories #
my %COG; 		#my $cnt = 0;
for my $name (keys %$name_index_r){
	print STDERR "...getting COGs for $name\n" unless $verbose;
	get_COG_cat($name, $runID, $name_index_r, $options, \%COG);
				#$cnt++; last if $cnt > 2;
	}

# normalizing each org by total COG #
normalize_COG(\%COG);

# writing out metadata table #
write_metadata(\%COG, $labels_r);


### Subroutines
sub write_metadata{
	my ($COG_r, $labels_r) = @_;
	
	# summing by category #
	my %COG_sum;
	foreach my $org (keys %$COG_r){
		foreach my $cat (keys %{$COG_r->{$org}}){
			$COG_sum{$cat} += $COG_r->{$org}{$cat};
			}
		}
	
	# ordering by COG category sum #
	my @COG_cat = sort{$COG_sum{$b}<=>$COG_sum{$a}} keys %COG_sum;
	
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
			else{ push(@line, "0"); }	
			}
		# changing organisms names back #
		print join("\t", $org, @line), "\n";
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
	
sub normalize_COG{
	my ($COG_r) = @_;
	foreach my $org (keys %$COG_r){
		my $sum = 0;
		map{ $sum += $COG_r->{$org}{$_} } keys %{$COG_r->{$org}};		# getting organism sum #
		map{ $COG_r->{$org}{$_} /= $sum } keys %{$COG_r->{$org}};		# normalizing
		map{ delete $COG_r->{$org}{$_} if $COG_r->{$org}{$_} < $COG_cutoff } keys %{$COG_r->{$org}};	# removing any categories below cutoff 
		}
		#print Dumper %$COG_r; exit;
	}

sub get_COG_cat{
# getting the COG categories of all the of genes specified #
	my ($name, $runID, $name_index_r, $options, $COG_r) = @_;

	my $grep_q =  join("", "fig|", $name_index_r->{$name}{'fig'}, ".peg");
	my $cmd = "echo \"$name_index_r->{$name}{'org'}\" | db_findClustersByOrganismList.py $options $runID | db_getGenesInClusters.py | grep \"$grep_q\"  | db_getExternalClusterGroups.py -d cog -g 3 | db_getExternalClustersById.py -c 13 | ";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, $cmd or die $!;
	
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		next unless /COG\d+/;
		my @line = split /\t/;
		
		s/^.+\[|\].*$|\/.+//g;		# removing all but COG category
		s/[-, ;:]+$//g;				# removing junk from end 
		s/[-, ;:]+/_/g;				# sanitizing
		
		$COG_r->{$name}{$_}++;
		}
	close PIPE;

		#print Dumper %$COG_r; exit;
	}

sub make_name_index{
# making an index that connects the provided labels to the fig & orgID the org file #
	my ($labels_r, $org_in) = @_;
	
	# making organism hash #
	my %sani_org;
	open IN, $org_in or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		#$org{$line[0]} = $line[1];
		
		(my $sanitized = $line[0]) =~ tr/-. /_/;
		
		$sani_org{$sanitized}{"org"} = $line[0]; 
		$sani_org{$sanitized}{"fig"} = $line[1];
		}
	close IN;
		
	# changing names #
	my %index;
	foreach my $name (@$labels_r){
		(my $sani_name = $name) =~ tr/-. /_/;
		if(exists $sani_org{$sani_name}){
			$index{$name} = $sani_org{$sani_name};
			}
		else{
			die " ERROR: $name not found in organism file (even with sanitized organism names)!\n	"
			}
		}
		
		#print Dumper %index; exit;
	return \%index;
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

cat leaf_labels.txt | tree_makeCOGMetadata.pl -r 

=head2 Required flags (depending on input)

=over

=item -t 	Tree file (newick or nexus).

=item -r 	ITEP cluster runID

=item -o 	ITEP organism file (for associating tree names to ITEP names).

=back

=head2 Options

=over

=item -f 	Tree format (only needed if '-t' provided; newick or nexus). [newick] 

=item -c 	COG category must be '-c' fraction of total COGs for the organism. [0.01]

=item -x 	Inclusiveness of clusters; flag(s) for db_findClustersByOrganismList.py. [-a]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc tree_makeCOGMetadata.pl

=head1 DESCRIPTION

Make an ITOL metadata table with the number of COG categories
associated with PEGs from each genome.

The COG categories are normalized by total COGs hits for the 
genome. The COGs are also sorted by category totals.

=head1 EXAMPLES

=head2 Tree file provided:

tree_makeCOGMetadata.pl -t tree.nwk -org organisms -r  mazei_I_2.0_c_0.4_m_maxbit > COG_metadata.txt

=head2 Just leaf labels provided:

cat leaf_labels.txt | tree_makeCOGMetadata.pl -org organisms -r  mazei_I_2.0_c_0.4_m_maxbit > COG_metadata.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

