#!/usr/bin/perl

use warnings;
use strict;

# include module global
use lib '/home/clara/Projekte/CMPERL/'; # search for moduls into cwd
use CM::rev_comp qw(rev_comp); # load reverse complement function
use CM::codons qw(codon translate);

### create all 6 reading frames from one genome sequence

## variables
# parametesr
my $name = "";

## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?n(?:ame)?=(.+)/)					{ $name = $1; }					# set name for output					-n	   | --name		for $name
}



# read genome from STDIN/PARAMETER # multiple chromosomes? -> hash?? # $chro {$header} = $genom 
sub rdgenom {
	my $genom = "";
	my $header = "";
	while (<STDIN>) {
		chomp;
		if ($_ =~/^(>.)/) { $header = $_; } # save header
		elsif ($_ =~/^[ACTGN]/) {
			$genom .= $_;
		}
	}
	return ($header, $genom);
}

# print out fasta containing lines of 80 nt
sub print80 {
	my $seq = shift;
	my $out = "";
	for (my $i = 0; $i <= (length $seq); $i+=80) {
		$out .= substr($seq, $i, 80);
		$out .= "\n";
	} 
	return $out;
}

sub savfl {
	my $out = shift;
	my $path = shift;
	open (OUT, ">$path") || die ($!);
	print OUT $out;
}

# sense
my ($header, $genom) = rdgenom; # frame #1 -> no changes
my $out = print80($genom);
#my $out = $genom;
my $out2 = print80((substr $genom, 1)); # frame #2 -> remove first nt
#my $out2 = substr $genom, 1;
my $out3 = print80((substr $genom, 2)); # frame #3 -> remove first + second nt
#my $out3 = substr $genom, 2;

# antisense
my $antigenom = rev_comp ($genom); # frame #1 -> revcomp seq
my $antiout = print80($antigenom);
#my $antiout = $antigenom; 
my $antiout2 = print80((substr $antigenom, 1)); # frame #2 -> revcomp seq - first nt
#my $antiout2 = substr $antigenom, 1;
my $antiout3 = print80((substr $antigenom, 2)); # frame #3 -> revcomp seq - first + second nt
#my $antiout3 = substr $antigenom, 2;


# print out different frames # each line -> 80 nt
# sense frame 1
savfl("$header--sense_frame1\n$out\n", "$name\_sense.fna");
#print "$header--sense_frame1\n$out";

# translate frames into AA
my $protS1 = translate $out;

savfl("$header--sense_frame1\n$protS1\n", "$name\_sense.faa");
#print "$header--sense_frame1\n$protS1";

# sense frame 2
savfl("$header--sense_frame2\n$out2\n", "$name\_sense2.fna");
#print "$header--sense_frame2\n$out2";

# translate frames into AA
my $protS2 = translate $out2;
savfl("$header--sense_frame2\n$protS2\n", "$name\_sense2.faa");
#print "$header--sense_frame1\n$protS2";


# sense frame 3
savfl("$header--sense_frame3\n$out3\n", "$name\_sense3.fna");
#print "$header--sense_frame3\n$out3";

# translate frames into AA
my $protS3 = translate $out3;
savfl("$header--sense_frame3\n$protS3\n", "$name\_sense3.faa");
#print "$header--sense_frame1\n$protS1";


# antisense frame 1
savfl("$header--antisense_frame1\n$antiout\n", "$name\_antisense.fna");
#print "$header--antisense_frame1\n$antiout";

# translate frames into AA
my $protA1 = translate $antiout;
savfl("$header--antisense_frame1\n$protA1\n", "$name\_antisense.faa");
#print "$header--sense_frame1\n$protA1";


# antisense frame 2
savfl("$header--antisense_frame2\n$antiout2\n", "$name\_antisense2.fna");
#print "$header--antisense_frame2\n$antiout2";

# translate frames into AA
my $protA2 = translate $antiout2;
savfl("$header--antisense_frame2\n$protA2\n", "$name\_antisense2.faa");
#print "$header--sense_frame1\n$protA2";


# antisense frame 3
savfl("$header--antisense_frame3\n$antiout3\n", "$name\_antisense3.fna");
#print "$header--antisense_frame3\n$antiout3";

# translate frames into AA
my $protA3 = translate $antiout3;

savfl("$header--antisense_frame3\n$protA3\n", "$name\_antisense3.faa");
#print "$header--sense_frame1\n$protA3";
