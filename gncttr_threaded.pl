#!/usr/bin/perl

use warnings;
use strict;

# include module global
use lib '/home/clara/Projekte/CCApredictor/gncttr/'; # search for moduls into cwd
use CM::GNCTTR::global qw($genom $seq $minsize $tmp $model $header $limit %starts %stops);
use CM::GNCTTR::codons qw(finds codon translate);
use CM::GNCTTR::hmm qw(nofltest compare);
use CM::GNCTTR::cds qw(cds);

use Thread::Pool; # Pool for multithreading


### find & extract genes from genomes


## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?g(?:enom)?=(.+)/)					{ $genom = $1; }					# genome for processing					-g | --genome			for $genom
	elsif ($para =~m/--?m(?:insize)?=(\d+)/)			{ $minsize = $1; }					# minsize for extracted sequences		-m | --minsize			for $minsize
	elsif ($para =~m/--?h(?:mm)?=(.*)/)					{ $model = $1; }					# HMM Model vor testing					-h | --hmm				for $model
}


# initialize thread pool
my $pool = Thread::Pool->new(
	{
		optimize => 'cpu',
		workers => 10, # threads
		do => sub {
			my ($starts, $top, $first, $codons) = @_;
			my $cds = cds ($starts, $top, $first, $codons);

			# split string $seq into hash -> ">xxx" => [ATCG]+
			my %result;
			foreach (sort keys %$cds) {
				
				# translate cds into AA
				my $AA = translate($$cds{$_});

				# test HMM
				my $result = nofltest ("$_\n$AA", $model);
				unless ($result == 0) {
					$result{"$_\n$AA"} = @$result[2];
					print "#$_\n$AA => @$result[2]\n";
				}
			}
			compare \%result;	
		}
	}
);

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
				# define & start jobs
				$pool->job(($starts, $_, $first, $codons));
			}
			$pool->join; #wait for all worker threads to finish
		}
	}
	$pool->shutdown;
}


## programm

# read genome
if (@ARGV) {
	main;
}
