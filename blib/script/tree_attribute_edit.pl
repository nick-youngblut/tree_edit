#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
my $mod = "4/27/12 9:28 AM";
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
use POSIX;

### global variables
my ($error);

### I/O
if ($#ARGV < 0){
	&usage;
	}
my ($treein, $attribute, $change, $verbose, $outfile);
GetOptions(
	   "tree=s" => \$treein,
	   "attribute=s" => \$attribute,
	   "change=s" => \$change,
	   "verbose" => \$verbose,
	   "outfile=s" => \$outfile,
	   "help|?" => \&usage # Help
	   );

### Input error check
die " ERROR provide a tree file!" if ! $treein;
$treein = rel_2_abs($treein);
($outfile = $treein) =~ s/\.[^\.]+$/_e.tre/;
$change = "round" if ! $change;
$attribute = "Abund" if ! $attribute;

### Routing main subroutines
open my $infh, $treein or die $!;
open my $outfh, ">$outfile" or die $!;
parse_nexus_tree($infh, $outfh, $attribute, $change);

close $infh; close $outfh;

#----------------------Subroutines----------------------#
sub parse_nexus_tree{
	my ($infh, $outfh, $attribute, $change) = @_;
	while (<$infh>){
		if($_ =~ /^\s*begin trees/i){
			my @line = split(/\[|\]/, <$infh>);
			my $line_ref = change_attribute(\@line, $attribute, $change);
			print $outfh "$_$line_ref";
			}
		else{ print $outfh $_; }
		}
	}

sub change_attribute{
	my ($line_ref, $attrib, $change) = @_;
	my @newline;
	foreach(@$line_ref){
		if($_ =~ /\s*&/ && $_ !~ /&[RU]$/){	#if attirbute
			my %tmp = split(/=|\s*,\s*/, $_);
			foreach(keys %tmp){
				if ($_ =~ /^[&]*$attrib$/ && $change =~ /round/i){
					$tmp{$_} = sprintf("%.2f", $tmp{$_}) if $tmp{$_} =~/^[\d\.e-]+$/;
					}
				}
			# remaking attribute #
			my @newattrib;
			foreach(sort keys %tmp){  push(@newattrib, join("=", $_, $tmp{$_}));  }
			push(@newline, join("", "[", join(", ", @newattrib), "]"));
			}
		elsif($_ =~ /\s*&/ && $_ =~ /&[RU]$/){	#if rooting
			push(@newline, join("", "[", $_, "]"));
			}
		else{ push(@newline, $_); }
		}
		#print Dumper(@newline); exit;
	return join("", @newline);
	}

sub rel_2_abs{
	### usage = filename, "f,d, or fd" ###
	### f = file, d = directory, fd = either [default = both] ###
	$_[0] = File::Spec->rel2abs($_[0]);
	$_[1] = "df" if ! $_[1];
	die " ERROR: directory not found:\n $_[0]\n" if ! -d $_[0] && $_[1] eq "d";
	die " ERROR: file not found:\n $_[0]\n" if ! -f $_[0] && $_[1] eq "f";
	die " ERROR: directory/file not found:\n $_[0]\n" if ! -d $_[0] && ! -e $_[0] && $_[1] =~ /[df]+/;
	return $_[0];
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
 tree_attribute_edit.r -t [-a] [-c]
Options:
 -t 	Tree file name [nexus format only!].
 -a 	Attribute name (same as in tree file).
 			default = 'Abund'
 -c 	How to change the attribute:
 			round = round integer to 2 decimal places.
 			default = 'round'
Description:
 The program edits the attributes in a nexus tree file.
 Currently, it can only edit the 1st tree in the file.
 Also, only a rounding feature is offered.
Notes:
	# only nexus file with 1 tree #
	Version: $version
	Last Modified: $mod
	Author: $author
Categories:
	Phylogeny
		
HERE
	print $usage;
    exit(1);
}