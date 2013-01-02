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

my ($verbose, $sum_bool);
GetOptions(
	   "sum" => \$sum_bool,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $tbl_ref = hab_matrix_2_table();
if($sum_bool){
	$tbl_ref = sum_by_label($tbl_ref) if $sum_bool;			# summing by label position
	write_hab_table_summed($tbl_ref);
	}
else{
	write_hab_table($tbl_ref); 
	}

### Subroutines
sub write_hab_table_summed{
# writing out the habitat table in long format (for ggplot2) #
	my $tbl_ref = shift;

	print join("\t", qw/Habitat Label_Position Label Proportion/), "\n";
	foreach my $hab (sort {$a<=>$b} keys %$tbl_ref){
		foreach my $pos (sort {$a<=>$b} keys %{$$tbl_ref{$hab}}){
			foreach my $lab (sort keys %{$$tbl_ref{$hab}{$pos}}){
				print join("\t", $hab, $pos, $lab, $$tbl_ref{$hab}{$pos}{$lab}), "\n";
				}
			}
		}

	}

sub write_hab_table{
# writing out the habitat table in long format (for ggplot2) #
	my $tbl_ref = shift;

	print join("\t", qw/Habitat Label Proportion/), "\n";
	foreach my $hab (sort {$a<=>$b} keys %$tbl_ref){
		foreach my $lab (sort keys %{$$tbl_ref{$hab}}){
			print join("\t", $hab, $lab, $$tbl_ref{$hab}{$lab}), "\n";
			}
		}

	}
	
sub sum_by_label{
# summing proportions by each label position (and habitat) #
	my $tbl_ref = shift;

	# making an index of what to sum #
	my %hab_idx;
	foreach my $hab (keys %$tbl_ref){
		foreach my $label (keys %{$$tbl_ref{$hab}}){
			my @lab_pos = split //, $label;
			for my $i (0..$#lab_pos){
				push(@{$hab_idx{$hab}{$i}{$lab_pos[$i]}}, $label);
				}
			}
		}
	
	# summing by index #
	my %hab_sum;
	foreach my $hab (keys %hab_idx){
		foreach my $pos (keys %{$hab_idx{$hab}}){
			foreach my $arrkey (keys %{$hab_idx{$hab}{$pos}}){
				foreach my $lab (@{$hab_idx{$hab}{$pos}{$arrkey}}){
					$hab_sum{$hab}{$pos}{$arrkey} += $$tbl_ref{$hab}{$lab};
					}
				}
			}
		}
		#print Dumper %hab_sum; exit;
	return \%hab_sum;
	}

sub hab_matrix_2_table{
# converting a habitat matrix to a table that R can read (long ggplot format) #
# columns: habitat#, label, proportion
	my %hab;
	my $hab_name;
	while(<>){
		chomp;
		next if /^\s*$/;
		die " ERROR: the habitat matrix should be only 1 line long\n" if ! eof && $.>1; 
		my @habs = split /\s*[{}]\s*/;
			#print Dumper @habs; exit;
		foreach my $i (1..$#habs){
			$habs[$i] =~ s/'//g ;
			my @tmp = split /,/, $habs[$i];
			if(scalar (grep (! /^\s*$/, @tmp)) == 1){
				shift @tmp if ! $tmp[0];
				$hab_name = $tmp[0];
				$hab_name =~ s/^ +|:|habitat //gi;
				}
			else{
				for my $ii (0..$#tmp){
					my @lab_port = split / *: */, $tmp[$ii];
					$lab_port[0] =~ s/^ +//;
					$hab{$hab_name}{$lab_port[0]} = $lab_port[1];
					}
				}
			}
		}
	return \%hab;
	}

__END__

=pod

=head1 NAME

adaptml_hab-mtx_convert.pl -- Convert an AdaptML habitat matrix to a table for R plotting

=head1 SYNOPSIS

adaptml_hab-mtx_convert.pl < habitat.matrix > habitat.mtx.txt

=head2 options

=over

=item -h	This help message

=back

=head2 For more information:

perldoc adaptml_hab-mtx_convert.pl

=head1 DESCRIPTION

Convert an AdaptML habitat matrix to a table for R plotting. The table is in 'long' format
for easy plotting in R. The table is tab-delimited.

=head2 3 Column output:

=over

=item 1) AdaptML inferred habitat

=item 2) Taxon label

=item 3) Proportion of label in habitat (inferred by AdaptML)

=back

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/tree_edit

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

