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

my ($verbose, @range, @names, $append);
my ($delimiter, $punct, $taxa_in);
my $column = 2;
GetOptions(
	   "taxa=s" => \$taxa_in,
	   "range=i{,}" => \@range,
	   "names=s{,}" => \@names,
	   "delimit=s" => \$delimiter,
	   "append" => \$append,
	   "column=i" => \$column,
	   "w" => \$punct, 					# non-characters can be anything [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a list of taxon names (-t)!\n"
	unless $taxa_in;
die " ERROR: cannot find $taxa_in!\n"
	unless -e $taxa_in;
die " ERROR: provide a range designating clades!\n"
	unless @range;
$column--;
if(! $delimiter){
	if($append){  $delimiter = "\t"; }
	else{ $delimiter = "__"; }
	}

### MAIN
my $taxon_names_r = load_taxon_names($taxa_in);
my $index_r = make_name_clade_index($taxon_names_r, \@range, \@names);


if($append){
	append_clade($index_r, $delimiter, $column);
	}
else{
	rename_taxa_add_clade($index_r, $delimiter);
	}


### Subroutines
sub rename_taxa_add_clade{
	my ($index_r, $delimiter) = @_;
	
	while(<>){
		chomp;
		next if /^\s*$/;
		
		foreach my $taxon_regex (keys %$index_r){
			my $r = join("", $delimiter, $index_r->{$taxon_regex});
			$_ =~ s/($taxon_regex)/$1$r/g;
			#last if $_ =~ /$taxon_regex/;
			}
		
		print "$_\n";
		}
	}

sub append_clade{
# using regex to determine clade; appending clade name as final column #
	my ($index_r, $delimiter, $column) = @_;
	
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /$delimiter/;
		
		foreach my $taxon_regex (keys %$index_r){
			my $regex = join("", "^", $taxon_regex, "\$");
			if($line[$column] =~ /$regex/){
				push @line, $index_r->{$taxon_regex};
				last;
				}
			}
		
		print join($delimiter, @line), "\n";
		}

	}

sub make_name_clade_index{
# making hash: name => clade #
	my ($taxon_names_r, $range_r, $names_r) = @_;
	
	# filling out range if missing last value #
	if(scalar(@$range_r) % 2){	# odd number
		push @$range_r, scalar @$taxon_names_r;
		}
	
	# 0-indexing range values #
	map{$_--} @$range_r;
	
	# making index #
	my %index;
	for (my $i=0;$i<=$#$range_r;$i+=2){
		die " ERROR: ranges should be ascending!\n" 
			unless $$range_r[$i] <= $$range_r[$i+1];
		
		foreach my $taxon (@$taxon_names_r[@$range_r[$i]..$$range_r[$i+1]]){
			# making regex #
			$taxon =~ s/[^A-Za-z0-9]/./g unless $punct;
			
			# loading hash #
			if($$names_r[$i/2]){
				$index{$taxon} = $$names_r[$i/2];
				}
			else{
				$index{$taxon} = $i/2 + 1;
				}
			}
		}
	
	# warning if range does not extend to all taxa #
	print STDERR " WARNING: ranges provided do not include all taxa!\n"
		unless scalar keys %index == scalar @$taxon_names_r;
	
		#print Dumper %index; exit;
	return \%index;
	}

sub load_taxon_names{
# loading names from tree; using nw_display #
	my $taxa_in = shift;
	
	open PIPE, "perl -p -e 's/[)]\\d+/)/g' $taxa_in | nw_display -S - |" or die $!;
	
	my %names;
	my @names;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		next unless $. % 2; 	# odd line numbers
		s/.+[|+]//;
		s/^ +| +$//g;
		s/ /_/g;
		
		print STDERR " WARNING: $_ found multiple times in names list!\n"
		if exists $names{$_};
		$names{$_} = $.;
		push @names, $_;
		}
	close PIPE;
	
	die " ERROR: no taxon names provided!\n"
		unless @names;
		
		#print Dumper @names; exit;
	return \@names;
	}

sub load_taxon_names_OLD{
# loading names 
	my $taxa_in = shift;
	open IN, $taxa_in or die $!;
	
	my %names;
	my @names;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		print STDERR " WARNING: $_ found multiple times in names list!\n"
			if exists $names{$_};
		$names{$_} = $.;
		push @names, $_;
		}
	close IN;
	
	die " ERROR: no taxon names provided!\n"
		unless @names;

	return \@names;
	}


__END__

=pod

=head1 NAME

rename_addCladeName.pl -- appending clade names to organism names or as a last column in a table

=head1 SYNOPSIS

cat (table or other file) | rename_addCladeName.pl [flags] -t tree.nwk > clades_appended_file

=head2 Required flags

=over

=item -tree  <char>

Tree file (newick format).

=item -range  <int>

Ranges delimiting clades as orderd by tree (index by 1; see DESCRIPTION)

=back

=head2 Optional flags

=over

=item -names  <char>

Clade names (see DESCRIPTION)

=item -delimit  <char>

Delimiter for appending clade name to taxon name 
OR parsing table and appeanding clade as last column. ["__"] or ["\t"] if '-append'

=item -append  <bool>

Append clade names as column on table instead of adding directly to names? [FALSE]

=item -column  <int>

Column in table for taxon name matching (index by 1; only if '-append'). [2]

=item -w  <bool>

Non-alpha-numeric characters count as anything ('.') in regular expression? [TRUE]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc rename_addCladeName.pl

=head1 DESCRIPTION

Add clade names to taxa names in a tree file, table, or other file.

=head2 -range

The range should designate a set of taxa names as ordered in the tree (provided 
with '-taxa'). Tree ordering is from top of the tree to the bottom as displayed
by nw_display. For example: '-range 1 5' designates the 1st 5 taxa in the tree as
the same clade. Example2: '-range 1 3 4 5' will split 5 taxa into 2 clades.

Providing an odd number of arguments will cause the last range
to extend to the end of the taxa in the tree. For example: '-range 1 3 4' is the 
same as '-range 1 3 4 5' for a tree of 5 taxa.  

=head2 -name

If no names are provided, clades are named as the order found in the tree.

A name should be provided for each range. For example: '-range 1 3 4 6' should have
2 names '-name clade1 clade2'. 

=head1 EXAMPLES

=head2 Renaming taxa names in gene tree by clades determined from species tree.

cat gene.nwk | rename_addCladeName.pl -t species.nwk -r 1 3 4 6

=head2 Renaming taxa names in ITEP gene info table 

printf "all_I_2.0_c_0.4_m_maxbit\t1" | db_getClusterGeneInformation.py |
rename_addCladeName.pl -t tree.nwk -r 1 3 4 6

=head2 Appending clade names on end of ITEP gene info table 

printf "all_I_2.0_c_0.4_m_maxbit\t1" | db_getClusterGeneInformation.py |
rename_addCladeName.pl -t tree.nwk -a -r 1 3 4 6

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

