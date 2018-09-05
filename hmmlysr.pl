#!/usr/bin/perl

use warnings;
use strict;


# include module
use lib '/home/clara/Projekte/CMPERL/'; # search for moduls into cwd
use CM::file qw(readvz cleanvz wrt wrta);
use CM::hmmlysr_global qw($in $out $score $evalue $filter @scores @evalues);


### Programm for analysis of predicted CCAI/CCAII/A-add/CC-add from predicted CDS


## parameters
foreach my $para (@ARGV){
	if ($para =~m/--?i(?:n)?=(.+)/)					{ $in = $1; }				# input directory for processing						-i | --in			for $in
	elsif ($para =~m/--?o(?:ut)?=(.+)/) 				{ $out = $1; } 				# output path of results								-o | --out			for $out
	elsif ($para =~m/--?s(?:core)?=(.+)/)				{ $score = $1; }			# minimum score for output								-s | --score		for $score
	elsif ($para =~m/--?e(?:value)?=(.+)/)				{ $evalue = $1; }			# minimum evalue for output								-e | --evalue		for $evalue
	elsif ($para =~m/--?f(?:ilter)?=(.+)/)				{ $filter = $1; }			# filter for inpput files								-f | --filter		for $filter
}


## Programm

# read directory & clean up empty files
if ($in) {
	my $dir = cleanvz $in;

	## read files -> get scores of hits
	foreach (@$dir) {
		unless ($_ =~/$filter/) { next; } # filter files

		open (my $LOG, "<$in/$_") || die ($!);
		while (<$LOG>) {
			chomp;
			if ($_ =~/^\[HIT\]/) { #  grep "\[HIT\]"
				my @result = split (/\s+/, $_); # my @result = split (/\t/, $line) -> @result[1] = evalue, @result[2] = score
				#print "$result[1]\t$result[2]\n";
				if ($result[2] >= $score && $result[1] <= $evalue) {
					push @scores, $result[2];
					push @evalues, $result[1];
				}
			}			
		}
		close $LOG;
	}
}


# write scores into $ARGV[1].scores
wrta ("$out.scores", \@scores);

# write evalues into $ARGV[1].evalues
wrta ("$out.evalues", \@evalues);
