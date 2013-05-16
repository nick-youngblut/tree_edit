#!/usr/bin/env perl

### packages/perl_flags
use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use File::Spec;

### I/O
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $name_in, $regex);
GetOptions(
	   "tree=s" => \$tree_in,
	   "name=s" => \$name_in,
	   "regex=s" => \$regex,
	   "help|?" => \&pod2usage 		# Help
	   );

### Input error check
die " ERROR: you must provide a regex in format of 'find/replace'  (must be in quotes).\n"
	if ! $regex; 
my @regs = split(/\/|\|/, $regex);
$regs[0] = qr/$regs[0]/;

die " ERROR: provide both and newick tree file and a name file (mothur format)!\n"
	if ! $tree_in || ! $name_in;
$tree_in = File::Spec->rel2abs($tree_in);
$name_in = File::Spec->rel2abs($name_in);

### Routing main subroutines
my $name_ref = load_name($name_in);
my $sub_list = make_sub_list($name_ref, \@regs);
expand_names($sub_list, $tree_in);

#----------------------Subroutines----------------------#
sub expand_names{
	### expanding names in tree based on sub_list (from names file) ###
	my ($sub_list, $tree_in) = @_;
	open IN, $tree_in or die $!;
	while (<IN>){
		s/#.+//;
		next if /^\s*$/;
		foreach my $sub (keys %$sub_list){
			if(scalar @{$sub_list->{$sub}} == 1){		# if just changing name
				s/$sub/${$sub_list->{$sub}}[0]/;
				}
			else{							# if adding clonal taxa with different habitat label
				my $let = "a";
				my @rep;
				foreach(@{$sub_list->{$sub}}){ 
					push (@rep, join("", $_, "$let:0"));
					$let++; 
					}
				my $rep = join(",", @rep);
				s/$sub/($rep)/;
				}
			}
		print $_;
		}
	close IN;
	}

sub make_sub_list{
	### making list for substituting in names ###
	my ($name_ref, $regs_ref) = @_;
	my %sub_list;
	my $count = 0;
	foreach my $rep (keys %$name_ref){
		$count++;
		my %cnt;
		foreach my $tax (@{$$name_ref{$rep}}){
			(my $tmp = $tax) =~ s/$$regs_ref[0]/"$$regs_ref[1]"/ee;
			$cnt{"$tmp\_$count"}++;
			}
		#print Dumper %cnt; exit;
		$sub_list{$rep} = [keys %cnt];
		}
		#print Dumper %sub_list; exit;	
	return \%sub_list;
	}

sub load_name{
	### loading name file ###
	my $name_in = shift;
	open IN, $name_in or die $!;
	my %name;
	while(<IN>){
		chomp;
		if($_ =~ /^\s+#/){next;}
		my @tmp = split(/\t|,/);
		$name{shift(@tmp)} = \@tmp;
		}
	close IN;
	die " ERROR: name file has nothing in it!\n" if scalar keys %name == 0;
	return \%name;
	}

sub error_routine{
	my $error = $_[0];
	my $exitcode = $_[1];
	print STDERR "ERROR: $error\nSee help: [-h]\n";
	exit($exitcode);
	}

=pod

=head1 NAME

adaptML_name_tree.pl -- 'expand' the rep. taxa names in a tree to include clonal reps.

=head1 SYNOPSIS

=head2 Getting leaf names

adaptML_name_tree.pl -t -n -r > expanded.nwk

=head2 options

=over

=item -tree

Tree file (newick format)

=item -name

Name file (Mothur format)

=item -regex

Regex for getting traits from names (see description)

=item -help

This help message

=back

=head2 For more information:

perldoc adaptML_name_tree.pl

=head1 DESCRIPTION
Renames and/or removes leaves from a
nexus or newick tree. Renaming/removal is based
on a file of old names and new names (or 'delete'
to remove the leaf).

=head2 General example

OTU1 in the tree is found in Sample1 & Sample2.
It would then become (S1_1a:0, S2_1b:0) in the
newick file.

=head2 REGEX ('-r')

Regex for getting traits from names.
The 'trait' protion of the name should be all that remains.

=head3 Example:

-r '\\d{6}_([^_]+)_.+/\$1'

changes '090531_SSBS_IN_59_E11' to 'SSBS' 

=head1 EXAMPLES

=head2 Usage

adaptML_name_tree.pl -t tree.nwk -n mothur_file.name -r '\\d{6}_([^_]+)_.+/\$1' > out.nwk

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

