#!/usr/bin/perl

use warnings;
use strict;

# include module global
use lib '/home/clara/Projekte/CMPERL/'; # search for moduls into cwd
use CM::gncttr_global qw($minsize $model $genom $score $name);
use CM::hmm qw(nofltest);


### programm for splitting translated AA sequences

## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?g(?:enom)?=(.+)/)					{ $genom = $1; }					# genome for processing					-g | --genome			for $genom
	elsif ($para =~m/--?m(?:insize)?=(\d+)/)			{ $minsize = $1; }					# minsize for extracted sequences		-m | --minsize			for $minsize
	elsif ($para =~m/--?h(?:mm)?=(.+)/)					{ $model = $1; }					# HMM Model vor testing					-h | --hmm				for $model
	elsif ($para =~m/--?s(?:core)?=(.+)/)				{ $score = $1; }					# score limit for output				-s | --score			for $score
	elsif ($para =~m/--n?(?:ame)?=(.+)/)				{ $name = $1; }						# name for output						-n | --name				for $name		
}


# read seq 
if ($genom) {
	my $i = 0;
	open (my $FAA, $genom) || die($!);
	while (<$FAA>) {
		chomp;
		# find '*' -> split
		my @cds = split (/\*/, $_);
		foreach (@cds) {
			$i++;
			if (length($_) >= $minsize) {
				my $result = nofltest (">$i\n$_", $model);
				if ($result ne 0 && @$result[2] >= $score) {
					print ">$name\_CDS$i@$result\t\n$_\n"; # print out "CDS"
				}
			}
		}
	}
	close $FAA;
}
