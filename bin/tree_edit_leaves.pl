#!/usr/bin/env perl

### packages/perl_flags
use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;


### global variables
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

### I/O
my ($tree_in, $format, $name_in, $get_names, $verbose, $outfile, $bootstrap);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "name=s" => \$name_in,
	   "get" => \$get_names,
	   "bootstrap" => \$bootstrap, 		# bootstrap tree? [TRUE]
	   "verbose" => \$verbose,
	   "outfile=s" => \$outfile,
	   "help|?" => \&pod2usage # Help
	   );

### Input error check
$tree_in = rel_2_abs($tree_in);
$name_in = rel_2_abs($name_in);
die " ERROR: provide both a tree file and a name file!\n" if (! $tree_in || ! $name_in && ! $get_names) || (! $tree_in && $get_names);
$format = guess_tree_format($tree_in) if ! $format;
format_check($format);

### Routing main subroutines
my $treeio = tree_io($tree_in, $format);
get_leaf_names($treeio, $tree_in) if $get_names;
my $name_ref = load_name($name_in);

my $out = new Bio::TreeIO(-fh => \*STDOUT,
						-format => $format);


while(my $tree = $treeio->next_tree){
	$tree = remove_rename_leaves($tree, $name_ref);
	$tree = reroot_tree($tree);
	$out->write_tree($tree);
	print "\n";
	}


#----------------------Subroutines----------------------#
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

sub get_leaf_names{
	### printing out leaf names ###
	my ($treeio, $infile) = @_;
	my %names;
	while(my $tree = $treeio->next_tree){		
		my @leaves = $tree->get_leaf_nodes;
		foreach(@leaves){
			if(exists($names{$_->id})){
				print STDERR " ERROR: different or repeditive leaf names in tree file!\n";
				}
			else{ 
				$names{$_->id} = 1; 
				}
			}
		last unless $bootstrap; 		# just using 1st tree
		}

	foreach (keys %names){  print "$_\n"; }
	exit;
	}

sub format_check{
	my $format = shift;
	if ($format =~ /nwk|newick/i){ $format = "newick" if $format =~ /nwk|newick/i; }
	elsif ($format =~ /nex|nexus|tre/i){ $format = "nexus"; }
	else{ die " ERROR: format not recognized [newick or nexus accepted\n"; }
	return $format;
	}

sub guess_tree_format{
	### guessing tree format by file extension (nwk or tre or nex) ###
	my $errormsg = " ERROR: cannot guess tree format based on file extension!\n";
	(my $ext = $_[0]) =~ s/.+\.//;
	die $errormsg if ! $ext;
	my $format;
	if($ext =~ /nwk/i){ $format = "newick"; }
	elsif($ext =~ /tre|nex/i){ $format = "nexus"; }
	else{ die $errormsg; }
	return $format;
	}

sub tree_io{
	### loading tree object ###
	my ($tree_in, $format) = @_;
	my $treeio = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);	
	return $treeio;
	}

sub load_name{
	### loading name file ###
	my $name_in = shift;
	open(NAMEFILE, $name_in) or die "\nERROR: name file not found!\n";
	my %name;
	while(<NAMEFILE>){
		chomp;
		if($_ =~ /^\s+#/){next;}
		my @tmp = split(/\t|,/);
		$name{shift(@tmp)} = \@tmp;
		}
		#print Dumper(%name);exit;
	return \%name;
	}

sub remove_rename_leaves{
	my ($tree, $name_ref) = @_;
	# checking for leaf existence in tree #
	my %nodeids;
	for my $leaf ($tree->get_leaf_nodes){
		$nodeids{$leaf->id} = 1;
		}
	foreach (keys %$name_ref){  
 		if(! exists($nodeids{ $_ })){
 			print STDERR " ERROR: $_ not found in tree!\n";
 			}
		}
	# remove/rename leaf #
	for my $node ($tree->get_nodes){
		if($node->id && exists($$name_ref{$node->id})){
			if(${$$name_ref{$node->id}}[0] =~ /delete|remove/i){
				$tree->remove_Node($node);
				if($verbose){ print $node->id, " deleted!\n"; }
				}
			else{
				$node->id(${$$name_ref{$node->id}}[0]);
				}
			}
		}
		my $root = $tree->get_root_node;
		#print Dumper($root->id); exit;
	return $tree;
	}

sub reroot_tree{
	### moving root to more inclusive root if all divergent taxa have benn removed ###
	my $tree = shift;
	my @nodes_with_ids;
	for my $child ($tree->get_leaf_nodes){		# finding all leaves not deleted
		push(@nodes_with_ids, $child) if $child->id;
		}		
	if(scalar(@nodes_with_ids) < scalar($tree->get_leaf_nodes)){		# if nodes have been deleted
		my @leaves = $tree->get_leaf_nodes;
		for my $leaf (@leaves){		#random rerooting
			if ($leaf->id){
				$tree->reroot($leaf->ancestor);
				print STDERR " tree rerooted on: ", $leaf->id, " ancestor because other leaves were deleted.\n";
				last;
				}
			}
 		for my $leaf (@leaves){		#removing blank leaves
 			$tree->remove_Node($leaf->ancestor) if ! $leaf->id;
 			}
		}
	return $tree;
	}

sub error_routine{
	my $error = $_[0];
	my $exitcode = $_[1];
	print STDERR "ERROR: $error\nSee help: [-h]\n";
	exit($exitcode);
	}

__END__

=pod

=head1 NAME

tree_edit_leaves.pl -- remove and/or rename leaves

=head1 SYNOPSIS

=head2 Getting leaf names

tree_edit_leaves.pl -t -g [-f] > names.txt

=head2 Renaming/removing leaves

tree_edit_leaves.pl -t -n [-f] [-o] [-v] > editted_tree

=head2 options

=over

=item -tree

Tree file (nexus or newick). [newick]

=item -format

Tree file format 

=item -get

Get leaf names and quit? [FALSE]

=item -bootstrap

Bootstrap trees (if multi-tree file)? [TRUE]

=item -name

Name file; 2 columns: old_name, new_name (or 'delete')

=item -verbose

Verbose output

=item -help

This help message

=back

=head2 For more information:

perldoc tree_edit_leaves.pl

=head1 DESCRIPTION
Renames and/or removes leaves from a
nexus or newick tree. Renaming/removal is based
on a file of old names and new names (or 'delete'
to remove the leaf).

=head1 WARNING

Deletion of leaves can cause rerooting!

=head1 EXAMPLES

=head2 Getting leaf names

tree_edit_leaves.pl -t file.nwk -g > names.txt

=head2 edit names.txt with favorite text editor (adding 2nd column)

vim names.txt

=head2 Renaming/removing leaves

tree_edit_leaves.pl -t file.nwk -n names.txt > file_rn.nwk

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

