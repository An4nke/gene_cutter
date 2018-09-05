#!/usr/bin/perl

use warnings;
use strict;



### find & extract genes from genomes


## variables
# parameter
my $genom = "";
my $seq = 0;
my $minsize = 100; # minsize of genes
my $tmp = "/tmp";
my $model = "";
my $header = "";

# data
#my $start = 'ATG'; # ATG, TTG, CTG
my @starts = ('ATG', 'TTG', 'CTG', 'GTG');
#my $stop = 'TAG'; #TAA, TAG, TGA
my @stops = ('TAA', 'TAG', 'TGA');

# statistic
my $plus = 0;
my $minus = 0;
my $hit = 0; # HMM hits
my $loose = 0; # no hits of HMM


# else
my $CCAII = "MNKILTKFRPHALKIKDAGGTLYIVGGAVRRFLMCQTPHDVDFCVCGLTVDTFKDLFPDARQQGNQFPVFVVDDCEFAMARTEKKVCAGYNGFEINSDPSVTIEEDLCRRDLTINSIAIDVITGMIVDPFNGAEDLINGIITPTSDAFKEDPVRVIRAARFMCEYPTFRASPTLIMYMLDLVHEIKLISDDHKFSELKKVFSSPKPSRYFNILNLGHALYITFPEIYSLIGIPQAHHSDGDAFEHTMRVLDDCRELTDDPVCLFAALTHDLGKATTPTEILPAHHDHETRSVEIIDKIDWVPNEWKYFAKVFAADHMRGHRFREMRRGKRVSLLERIHKSNRGLEGFCKVLYADKPTPGTMRDIALMHATYAKIYSISGDDLPATTPKGEAFGKVLHQKRVEMI";


## NCBI standard AA code
# TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys  
# TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys  
# TTA L Leu      TCA S Ser      TAA * Ter      TGA * Ter  
# TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp  

# CTT L Leu      CCT P Pro      CAT H His      CGT R Arg  
# CTC L Leu      CCC P Pro      CAC H His      CGC R Arg  
# CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg  
# CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg  

# ATT I Ile      ACT T Thr      AAT N Asn      AGT S Ser  
# ATC I Ile      ACC T Thr      AAC N Asn      AGC S Ser  
# ATA I Ile      ACA T Thr      AAA K Lys      AGA R Arg  
# ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg  

# GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly  
# GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly  
# GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly  
# GTG V Val      GCG A Ala      GAG E Glu      GGG G Gly 


## NCBI Bacterial, Archaeal and Plant Plastid Code
# TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys  
# TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys  
# TTA L Leu      TCA S Ser      TAA * Ter      TGA * Ter  
# TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp  

# CTT L Leu      CCT P Pro      CAT H His      CGT R Arg  
# CTC L Leu      CCC P Pro      CAC H His      CGC R Arg  
# CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg  
# CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg  

# ATT I Ile i    ACT T Thr      AAT N Asn      AGT S Ser  
# ATC I Ile i    ACC T Thr      AAC N Asn      AGC S Ser  
# ATA I Ile i    ACA T Thr      AAA K Lys      AGA R Arg  
# ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg  

# GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly  
# GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly  
# GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly  
# GTG V Val i    GCG A Ala      GAG E Glu      GGG G Gly  

my %code = (
	'TTT' => 'F', # Phe
	'TTC' => 'F', # Phe
	'TTA' => 'L', # Leu
	'TTG' => 'L', # Leu
	'TCT' => 'S', # Ser
	'TCC' => 'S', # Ser
	'TCA' => 'S', # Ser
	'TCG' => 'S', # Ser
	'TAT' => 'Y', # Tyr 
	'TAC' => 'Y', # Tyr 
	'TAA' => '*', # Ter
	'TAG' => '*', # Ter 
	'TGT' => 'C', # Cys
	'TGC' => 'C', # Cys 
	'TGA' => '*', # Ter
	'TGG' => 'W', # Trp 
	'CTT' => 'L', # Leu 
	'CTC' => 'L', # Leu
	'CTA' => 'L', # Leu 
	'CTG' => 'L', # Leu
	'CCT' => 'P', # Pro 
	'CCC' => 'P', # Pro
	'CCA' => 'P', # Pro
	'CCG' => 'P', # Pro
	'CAT' => 'H', # His
	'CAC' => 'H', # His
	'CAA' => 'Q', # Gln
	'CAG' => 'Q', # Gln
	'CGT' => 'R', # Arg
 	'CGC' => 'R', # Arg
	'CGA' => 'R', # Arg
	'CGG' => 'R', # Arg
	'ATT' => 'I', # Ile
	'ATC' => 'I', # Ile
	'ATA' => 'I', # Ile
	'ATG' => 'M', # Met
	'ACT' => 'T', # Thr
	'ACC' => 'T', # Thr
	'ACA' => 'T', # Thr
	'ACG' => 'T', # Thr
	'AAT' => 'N', # Asn
	'AAC' => 'N', # Asn
	'AAA' => 'K', # Lys
	'AAG' => 'K', # Lys
	'AGT' => 'S', # Ser
	'AGC' => 'S', # Ser
	'AGA' => 'R', # Arg
	'AGG' => 'R', # Arg
	'GTT' => 'V', # Val
	'GTC' => 'V', # Val
	'GTA' => 'V', # Val
	'GTG' => 'V', # Val
	'GCT' => 'A', # Ala
	'GCC' => 'A', # Ala
	'GCA' => 'A', # Ala
	'GCG' => 'A', # Ala
	'GAT' => 'D', # Asp
	'GAC' => 'D', # Asp
	'GAA' => 'E', # Glu
	'GAG' => 'E', # Glu
	'GGT' => 'G', # Gly
	'GGC' => 'G', # Gly
	'GGA' => 'G', # Gly
	'GGG' => 'G', # Gly
);


## parameters
foreach my $para (@ARGV) {
	if ($para =~m/--?g(?:enom)?=(.+)/)					{ $genom = $1; }					# genome for processing					-g | --genome			for $genom
	elsif ($para =~m/--?m(?:insize)?=(\d+)/)			{ $minsize = $1; }					# minsize for extracted sequences		-m | --minsize			for $minsize
	elsif ($para =~m/--?h(?:mm)?=(.*)/)					{ $model = $1; }					# HMM Model vor testing					-h | --hmm				for $model
	elsif ($para =~m/--?n(?:ame)?=(.+)/)				{ $header = $1; }					# name for extraction					-n | --name				for $header
	elsif ($para =~m/--?t(?:mp)?=(.+)/)					{ $tmp = $1; }						# directory for temp file saving		-t | --tmp				for $tmp
}

## subroutine

# write sequences into file # seqout ($file, $string)
sub seqout {
	my $file = shift;
	my $string = shift;
	print STDERR "[INFO]\topen $file for writing extrated random protein sequence..\n";
	open (OUT, ">$file") || die($!); # open file for reading
	print OUT $string; # write string into file
	close OUT;
}

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
			#$header = $1; # set header
			#$header =~s/ /_/g;
			next;
		}
		if ($_ =~/[ATCGN]/) {
			$seq .= $_;
		}
	}
	return $seq;
}

# testing hmm of cds # test (\@cds, $tmp_file)
sub test {
	my $tmp = shift;
	my $input = shift;
	seqout ($tmp, $input); # save cds into tmp file
	my $cmd = "hmmsearch --cpu 8 --noali $model $tmp";
	print STDERR "$cmd\n";
	my @tmp = split (/\n/, qx{$cmd}); # do hmmsearch, save and split output
	my @good = (grep (/^\s+\d+/, @tmp)); # extract lines containing results
	if ($good[0]) {
		$hit++;
		print "$good[0]\n";
		# save seq
	} else { $loose++; } # no hit detected
	unlink $tmp; # remove tmp file
}

# testing hmm without writing cds into file?
sub nofltest {
	my $string = shift;
	my $cmd = "echo \"$string\" | hmmsearch --noali $model -"; # use hmmsearch from stdin
	print STDERR "$cmd\n";
	my @tmp = split (/\n/, qx{$cmd}); # do hmmsearch, save and split output
	my @good = (grep (/^\s+\d+/, @tmp)); # extract lines containing results
	if ($good[0]) {
		$hit++;
		print "[HIT]\t$good[0]\n";
		# seqout ("$outfile", $string); # save seq
	} else { $loose++; } # no hit detected
}

# split string into codon
sub codon {
	my $seq = shift;
	my $name = shift;
	my $b = 0; # button for remembering START
	my $aa; # translated cds AA sequence
	my $s = 0; # startcodon nr of cds
	for (my $i = 0; $i < (((length $seq)/3) -1); $i++) {
		my $codon = substr($seq, $i*3, 3); # split into codons
		#if ($code{$codon} eq "M" && $b == 0) { # test start codon
		if ((grep $codon, @starts) && $b == 0) {
			$b = 1;
			$s = $i; # remember startcodon nr
		}
		#if ($codon eq $stop) { # test stop codon
		if ($code{$codon} eq "*" || $i >= (((length $seq)/3)) -2) {
		#if ((grep $codon, @starts) && $b == 1) {
			$b = 0;
			if ($aa && length $aa >= $minsize) { # save sequences with right size only
				#test ("$tmp/$s\_$i.faa", ">cds $s - $i\n$aa\n"); # write seq into tmp file -> test HMM -> rmv file after processing
				nofltest(">$name $s - $i\n$aa\n"); # testing HMM without writing into file
				print ">$name $s - $i\n$aa\n";
			}
		}
		if ($b == 1) { # write into cds
			$aa .= $code{$codon}; # save translated protein sequence
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


## programm

# test
#test ("$tmp/test.faa", ">test\n$CCAII");

# extract nt sequence into string
$seq = read_genom ($genom); # 1st reading frame
my $seq2 = substr $seq, 1; # 2nd reading frame
my $seq3 = substr $seq, 2; # 3rd reading frame

# split sequence into codons -> translate into AA -> test HMM
my $frame1 = codon ($seq, "sense 1");
my $frame2 = codon ($seq2, "sense 2");
my $frame3 = codon ($seq3, "sense 3");


# create reverse complement of seq
my $revseq = rev_comp ($seq); # 1st reading frame
my $revseq2 = substr $revseq, 1; # second reading frame
my $revseq3 = substr $revseq, 2; # third reading frame

# split reverse complement sequence into codons -> translate into AA -> test HMM
my $revcfrm1 = codon ($revseq, "antisense 1");
my $revcfrm2 = codon ($revseq2, "antisense 2");
my $revcfrm3 = codon ($revseq3, "antisense 3");


# print statistics
print STDERR "[STAT]\tHMM hits\t$hit\n[STAT]\tHMM blanks\t$loose\n";
