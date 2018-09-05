#!/usr/bin/perl

use warnings;
use strict;

# include module global
use lib '/home/clara/Projekte/CMPERL/'; # search for moduls into cwd
use CM::gncttr_global qw($genom $seq $minsize $tmp $model $header $limit %starts %stops);
use CM::codons qw(finds codon translate);
use CM::hmm qw(nofltest compare);
use CM::cds qw(cds);


### find & extract genes from genomes


## variables
# parameter
=stuff
my $genom = "";
my $seq = 0;
my $minsize = 100; # minsize of genes
my $tmp = "/tmp";
my $model = "";
my $header = "";
my $limit = 100000;

# data
#my @starts = ('ATG', 'TTG', 'CTG', 'GTG');
my %starts = (
	'ATG' => 1,
	'TTG' => 1,
	'CTG' => 1,
	'GTG' => 1
);
#my @stops = ('TAA', 'TAG', 'TGA');
my %stops = (
	'TAA' => 1,
	'TAG' => 1,
	'TGA' => 1
);

# statistic
my $plus = 0;
my $minus = 0;
my $hit = 0; # HMM hits
my $loose = 0; # no hits of HMM

=cut


## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?g(?:enom)?=(.+)/)					{ $genom = $1; }					# genome for processing					-g | --genome			for $genom
	elsif ($para =~m/--?m(?:insize)?=(\d+)/)			{ $minsize = $1; }					# minsize for extracted sequences		-m | --minsize			for $minsize
	elsif ($para =~m/--?h(?:mm)?=(.*)/)					{ $model = $1; }					# HMM Model vor testing					-h | --hmm				for $model
}

# read file containing genome # read_genom ($genom)
sub main {
	my $file = $genom;
	my $seq = "";
	my $first = 0;
	open (FILE, "<$file") || die ($!); # open genom file 
	print STDERR "[INFO]\tExtraction of potential gene sequences from $file..\n";
	while (<FILE>) {
		chomp;
		if ($_ =~/^[ATCGN]/) {
			# save sequence window of length $limit into string
			if (length $seq >= $limit) {
				# cut first 81 (3*N) nt of seq of sequence ($seq) 
				$seq = substr($seq , 81, (length $seq) - 81);
				$first += 81;
			}
			$seq .= $_;
			#print "$seq\n";
			my $codons = codon $seq; # split $seq into codons
			my ($starts, $stops) = finds $codons; # search stopcodon -> stopcodon? -> search for start codons

			# save all cds -> all starts - stop
			foreach (sort {$a <=> $b} keys %$stops) {
				my $cds = cds ($starts, $_, $first, $codons);

				# split string $seq into hash -> ">xxx" => [ATCG]+
				my %result;
				foreach (sort keys %$cds) {
					#print "$_\n";
					
					# translate cds into AA
					my $AA = translate($$cds{$_});
					#print "$AA\n";
			
					# test HMM
					my $result = nofltest ("$_\n$AA", $model);
					unless ($result == 0) {
						$result{"$_\n$AA"} = @$result[2];
					}
				}
				compare \%result;
			}
		}
	}
}


## programm

# read genome
if (@ARGV) {
	main;
}
