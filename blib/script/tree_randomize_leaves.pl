#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
my $mod = "8/23/12 1:12 PM";
my $version = "0.1";
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

### global variables
my ($error);

### I/O
my ($delim, $col);
GetOptions(
	   "delimiter=s" => \$delim,
	   "col=i" => \$col,
	   "help|?" => \&usage # Help
	   );

### Input error check
$col = 1 if ! $col;

### Routing main subroutines
my $tree_ref = get_tree();
foreach(@$tree_ref){
	my $names_ref = get_names($_);
	if ($delim){
		$names_ref = shuffle_column($names_ref, $delim, $col);
		}
	else{
		$names_ref = shuffle_names($names_ref);
		}
	$_ = change_names($_, $names_ref);
	}

#----------------------Subroutines----------------------#
sub change_names{
	### changing names based on key ###
	my ($tree, $names_ref) = @_;
	my @tree = split /:/, $tree;
	
	my @new_tree = @tree;
	foreach my $chg (keys %$names_ref){
			#print Dumper $chg;
		my $hit_cnt = 0;
		for my $i (0..$#tree){
				#print Dumper $tree[$i];
			if ($tree[$i] =~ /$chg$/){
				$hit_cnt++;
				(my $tmp = $tree[$i]) =~ s/$chg$/$$names_ref{$chg}/;
				$new_tree[$i] = $tmp;
				}
			}
		die " ERROR: 0 hits for $chg." if $hit_cnt < 1;
		print STDERR " WARNING: >1 hit for $chg.\n" if $hit_cnt > 1;
		}
		#print Dumper @new_tree; exit;
	print join(":", @new_tree), "\n"; 
	}

sub shuffle_names{
	### simple shuffling of names ###
	my $names_ref = shift;
	my @toshuf = @$names_ref;
	my $shuf_ref = fisher_yates_shuffle(\@toshuf);

	my %new_names;
	for my $i (0..$#$names_ref){
		$new_names{$$names_ref[$i]} = $$shuf_ref[$i]; 
		}
		#print Dumper %new_names;  exit;
	return \%new_names;
	}

sub shuffle_column{
	### shuffling column after delimiting names ###
	my ($names_ref, $delim, $col) = @_;
	my $delimg = qr/$delim/;
	
	# parsing names by delim #
	my (@toshuf, @other);
	foreach(@$names_ref){
		my @tmp =  split /$delimg/;
		push(@toshuf, splice(@tmp, $col-1, 1));
		push(@other, join("", @tmp));
		}
		
	# shuffling by fisher-yates method #
	my $shuf_ref = fisher_yates_shuffle(\@toshuf);
	
	# recombining #
	my %new_names;
	for my $i (0..$#$names_ref){
		if($col == 1){
			$new_names{$$names_ref[$i]} = join($delim, $$shuf_ref[$i], $other[$i]);
			}
		else{
			$new_names{$$names_ref[$i]} = join($delim, $other[$i], $$shuf_ref[$i]);
			}
		}
		#print Dumper %new_names; exit;
	return \%new_names;
	}

sub get_names{
	### getting names from tree ###
	my $tree = shift;
	$tree =~ s/[(),;]//g;
	my @tmp = split /:[\d\.]+/, $tree;
	
	my @tmp2;
	foreach(@tmp){ push(@tmp2, $_) if $_ !~ /^\s*$/; }
	
	return \@tmp2;
	}

sub get_tree{
	### getting tree file from stdin ###
	my @tree;
	while(<>){
		next if /^\s*$/;
		next if ! /^\s*\(.+;\s*$/;
		push(@tree, $_);
		}
	die " ERROR: newick tree does not seem to be formatted correcly!\n" if ! $tree[0];
		#print Dumper @tree; exit;
	return \@tree;
	}

sub error_routine{
	my $error = $_[0];
	my $exitcode = $_[1];
	print STDERR "ERROR: $error\nSee help: [-h]\n";
	exit($exitcode);
	}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
        }
   	return $array;
	}

sub usage {
 my $usage = <<HERE;
Usage:
  tree_randomize_leaves.pl [-d] [-c] < newick_file > newick_file
Options:
  -d 	Delimiter if only randomizing part of the name.
 		Example: "_" for only randomizing part of the 
 			names: "HH_01" & "SS_02"
  -c 	Column in name for randomization. Only needed if using
 		the -d flag. Index by 1.
 		[Default: 1]
Description:
  The program shuffles leaf names or parts of leaf names
  by the Fisher-Yates method.
  * multiple newick trees can be randomized
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