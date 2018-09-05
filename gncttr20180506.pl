#!/usr/bin/perl

use warnings;
use strict;



### find & extract genes from genomes


## variables
# parameter
my $genom = "";
my $seq = 0;
my $minsize = 5; # minsize of genes
my $header = "";

# data
my %starts = (
	'CAC' => '+',
	'CAA' => '+',	
	'CAT' => '+',	
	'ATG' => '-',	
	'GTG' => '-',	
	'TTG' => '-',	
);
my %stops = (
	'CTA' => '+',
	'TCA' => '+',
	'TTA' => '+',
	'TAG' => '-',
	'TGA' => '-',
	'TAA' => '-',
);
my %scoord; # POS => '+' || '-'
my %ocoord; # POS => '+' || '-'

# statistic
my $plus = 0;
my $minus = 0;


## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?g(?:enom)?=(.+)/)					{ $genom = $1; }					# genome for processing					-g | --genome			for $genom
	elsif ($para =~m/--?m(?:insize)?=(\d+)/)			{ $minsize = $1; }					# minsize for extracted sequences		-m | --minsize			for $minsize
	elsif ($para =~m/--?n(?:ame)?=(.+)/)				{ $header = $1; }					# name for extraction					-n | --name				for $header
}

## subroutine

# get reverse complement of sequence # $rev_comp_sequence = rev_comp ($sequence)
sub rev_comp {
	my $rev = reverse($_[0]); #String rueckwaerts einlesen
	$rev =~tr/ACTGactg/TGACtgac/; #Nucleotide uebersetzen
	return $rev;
}


# read file containing genome # read_genom ($genom)
sub read_genom {
	my $file = shift;
	my $seq = "";
	open (FILE, "<$file") || die ($!); # open genom file 
	print STDERR "[INFO]\tExtraction of potential gene sequences from $file..\n";
	while (<FILE>) {
		chomp;
		if ($_ =~/^>(.*)/ && $header eq "") {
			$header = $1; # set header
			$header =~s/ /_/g;
		}
		if ($_ =~/[ATCGN]/) {
		 $seq .= $_;
		}
	}
	return $seq;
}

# find codons inside genom sequence # code ($seq, $codon)
sub code {
	my $seq = shift;
	my $codon = shift;
	my $pos = 0;
	my @pos; # array containing codon position
	while ($pos != -1) { # as long as codon found
		my $offset = $pos + 1;
		$pos = index($seq, $codon, $offset);
		push @pos, $pos; # remember position of codon
	}
	return \@pos;
}

# find position of certain set of codons # codpos ($seq, \%codons, \%)
sub codpos {
	my $seq = shift; # input sequence
	my $codons = shift; # hash containing codons
	my $hit = shift; # output hash of coordinates
	foreach my $codon (sort {$$codons{$a} cmp $$codons{$b}} keys %{ $codons }) {
		my $pos = 0;
		while ($pos != -1) {
			my $offset = $pos + 1;
			$pos = index ($seq, $codon, $offset);
			if ($pos != -1) { $$hit{$pos} = $$codons{$codon}; }
		}
	}
}

# extract part of sequence
sub xtrc {
	my $seq = shift;
	my $start = shift;
	my $end = shift;
	# shorten sequence
	my $new = substr($seq, $start, $end - $start);
	return $new;
}

# extract nt sequence into string
$seq = read_genom ($genom);

# find codon position inside genom
print STDERR "[INFO]\tSearching for Codons..\n";
codpos ($seq, \%starts, \%scoord); # find all startcodons
codpos ($seq, \%stops, \%ocoord); # find all stopcodons

# found of start -> if + -> stop codon downstream? -> if - -> stop codon upstream?
print STDERR "[INFO]\tProcessing Data..\n";
foreach my $tart (sort {$a <=> $b} keys %scoord) {
	if ($scoord{$tart} eq '+' && (my @tmp = grep {$_ > $tart} sort {$a <=> $b} keys %ocoord)) { # '+' -> find next downstream coordinate inside %ocoords (STOP codons)
		$plus++; # count + strand
		if (length (my $gene = xtrc ($seq, $tart, $tmp[0])) > $minsize) {
			print ">$header--$tart-$tmp[0] ($scoord{$tart})\n";
			print "$gene\n";
		}
	} elsif ($scoord{$tart} eq '-' && (@tmp = grep {$_ < $tart} sort {$a <=> $b} keys %ocoord)) { # '-' -> find next upstream coordinate inside %ocoords (STOP codons)
		$minus++; # count - strand
		if (length (my $gene = xtrc ($seq, $tmp[0], $tart)) > $minsize){ # test if length of reverse complement sequences is long enough
			$gene = rev_comp ($gene);
			print ">$header--$tmp[-1]-$tart ($scoord{$tart})\n";
			print "$gene\n";
		}
	}
}

# print statistic
print STDERR "[STAT]\tnumber of genes: ".($plus + $minus)."\n[STAT]\tnumber of (+) genes: $plus (".($plus/($plus + $minus)*100)."%)\n[STAT]\tnumber of (-): $minus (".($minus/($plus + $minus)*100)."%)\n";


