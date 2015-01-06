#!/usr/bin/perl
my $mod = "2/16/12 8:57 AM";
my $version = "0.9";
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
if ($#ARGV < 0){
	&usage;
	}
my ($verbose, $format_in, $first_flag);
GetOptions(
	   "in_format=s" => \$format_in,	# input format
	   #"out_format=s" => \$format_out,	# output format
	   "first" => \$first_flag,
	   "verbose" => \$verbose,
	   "help|?" => \&usage # Help
	   );

### Input error check
@ARGV = sort{		# sorting by cluster
	my @a = $a =~ /cluster(\d+)/;
	my @b = $b =~ /cluster(\d+)/;
	$a[0] <=> $b[0]} @ARGV if $ARGV[0] =~ /cluster\d+/; 
foreach(@ARGV){
	$_ = File::Spec->rel2abs($_);
	if(! -e $_ && ! -d $_){ die " ERROR: cannot find $_\n"; }
	}
$format_in = "nexus" if ! $format_in;
$format_in = check_tree_format($format_in);

### main subroutines ###
cat_newick_trees($first_flag) if $format_in eq "newick";
cat_nexus_trees($first_flag) if $format_in eq "nexus";

#----------------------Subroutines----------------------#
sub cat_newick_trees{
	### simply concatenating trees ###
	my ($first_flag) = @_;
	foreach(@ARGV){
		open IN, $_ or die $!;
		while(<IN>){
			print $_;
			last if $first_flag;	# only using 1st tree if $first_flag
			}
		close IN;
		}
	print STDERR " Concatenated tree file written\n";
	exit;
	}

sub cat_nexus_trees{
	### concatenating nexus trees ###
	my ($first_flag) = @_;
	
	# getting taxa names #
	my %node_names;
	foreach(@ARGV){ $_ = get_names($_, \%node_names);  }
	
	# setting up names file #
	my $namesfile = unique_file("", "cat_tree.names");
	open my $namesfh, ">$namesfile", or die $!;
	
	# concatenating trees #
	for(my $i=0; $i<=$#ARGV; $i++){ name_and_cat_trees($ARGV[$i], \%node_names, $i, $first_flag, $namesfh);  }	# writting out all trees
	
	close $namesfh;
	}

sub get_names{
	### removing duplicates; checking for multiple genes in 1 tree ###
	my ($tree_file, $names_ref) = @_;
		#print STDERR "...processing $parts[2]\n" if $verbose;
	open my $treefh, $tree_file or die $!;
	my @names_last = keys %$names_ref;
	while(<$treefh>){
		chomp;
		if($_ =~ /Translate/){
			while(<$treefh>){
				$_ =~ s/^\s+//;
				last if $_ =~ /^\s*;/;	# last delimiter;
				my @tmp = split(/\s/);
				(my $base = $tmp[1]) =~ s/(_peg)*_[^_]+$//;
				die " ERROR name cannot begin with '_'\n" if ! $base;
				$$names_ref{$base} = $tmp[0];
				last if $_ =~ /;/;		# last delimiter;
				}
			}
		}
	die " ERROR: extra/less tip names in tree (paralogs?):\n $tree_file\n"  if scalar(keys %$names_ref) != scalar(@names_last) && $names_last[0]; 
	close $treefh;
	return $tree_file;
	}

sub name_and_cat_trees{
	### parsing out trees from each file
	### naming trees & writing out nexus file
	my ($tree_file, $node_names_ref, $tree_file_cnt, $first_flag, $namefh) = @_;
		#print Dumper();
	open IN, $tree_file or die $!;
	my $tree_cnt = 1; my %tree_names;	# tree_names = name order for this tree
	
	while (my $line = <IN>){
		if($line =~ /Translate/){ 	# getting name order for this tree #
			while (<IN>){
				$_ =~ s/^\s+//;
				last if $_ =~ /^\s*;/;	# last delimiter;
				my @tmp = split(/\s+/);
				(my $base = $tmp[1]) =~ s/(_peg)*_[^_]+$//;
				die " ERROR name cannot begin with '_'\n" if ! $base;
				$tree_names{$base} = $tmp[0];
				last if $_ =~ /;/;		# last delimiter;
				}
			}
		next if $line !~ /^\s*tree/;	# passing this means line is tree

		chomp $line;
		$line =~ s/^.+=.*?\(/(/;		# removing tree name
		
		$line = number2name($line, \%tree_names, "$tree_file\_tree$tree_file_cnt");	# checking label order

		if($tree_file =~ /cluster(\d+)/){ print $namefh "file$1\_tree$tree_cnt\n";  }	# name file corresponding to tree names (for renaming distance matrix)
		else{ print $namefh "file$tree_file\_tree$tree_file_cnt\n"; }	# name file corresponding to tree names (for renaming distance matrix)
		
		print $line,"\n";			#writting file
		last if $first_flag;	#one 1st tree if $first_flag
		$tree_cnt++;
		}
	close IN;
	}

sub number2name{
	### changing numbers in trees to name and back to number of 'main taxon key' #
	my ($line, $tree_names_ref, $file_tree_number) = @_;		# tree_names = this tree; node_names = final node order
		
	# changing order; numbers2genes (specific tree);
	print STDERR " Correcting order for $file_tree_number\n" if $verbose;
	foreach (keys %$tree_names_ref){ 
		die " ERROR: $$tree_names_ref{$_} not found in tree:\n $line\n" if $line !~ /([\(,])$$tree_names_ref{$_}/;
		$line =~ s/([\(,])$$tree_names_ref{$_}:/$1$_:/; 	
		}
	return $line;
	}

sub check_tree_format{
	if($_[0] =~ /new|newick/i){ $_[0] = "newick"; }
	elsif($_[0] =~ /nex|nexus/i){ $_[0] = "nexus"; }
	else{ die " ERROR: tree format not recognized\n  (must be newick or nexus)\n";}
	return $_[0];
	}
	
sub unique_file{
	# making unique log file
		# requires File::Spec
		# IN: directory, file_name
	my ($dir, $end) = @_;
	if(! $dir){ $dir = File::Spec->rel2abs(File::Spec->curdir()); }
	else{ $dir = File::Spec->rel2abs($dir); }
	my $ufile;
	my $cnt = 1;
	my @ends = split(/\./, $end);
	while(1){
		my $loopend;
		if(scalar(@ends) > 1){ $loopend = join("", @ends[0..$#ends-1], $cnt, ".", $ends[$#ends]); }
		else{ $loopend = $end . "$cnt"; }
		$ufile = join("", $dir, "/", file_date(), "_", $loopend);
		if(-e $ufile || -d $ufile){ $cnt++; next;}
		else{ last;}
		}
		#print Dumper($ufile); exit;
	return $ufile;
	}

sub file_date{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	$year = substr(($year + 1900),2);
	$mon++;
	##Correcting for single digits
	if ($mon !~ /\d\d/){
		$mon = "0$mon";
		}
	if ($mday !~ /\d\d/){
		$mday = "0$mday";
		}
	my $date = join("",$year,$mon,$mday);
	return $date;
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
 cat.pl [-i] [-o] [-f] [-w] ...
Options:
 -i 	Input tree file format (newick or nexus)
 			[default: nexus]
 -o 	Output tree file format (newick or nexus)
 			[default: newick]
 -f 	Use only first tree in file.
 			[default: FALSE]
 -w 	Written file name.
 ... 	Tree files (1 or more)
Description:
 The script formats and concatenates trees in
 multiple files into 1 tree file.
Notes:
	Version: $version
	Last Modified: $mod
	Author: $author
Required: 
	Bioperl
Categories:
	Phylogeny

HERE
	print $usage;
    exit(1);
}
