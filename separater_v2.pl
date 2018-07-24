#!/bin/env perl -w 
use strict; 

##################################################     separater.pl      ####################################################################
######	Separater the reading into three parts by using the adapter sequences.
######	Run Sample: 
#############################################################################################################################################

#command line inputs
my ($fasta,$len,$p1c,$p2c) = @ARGV;

#Open Fasta file and output file
open (FA,"./$fasta") or die "can't open gene file $fasta\n";
open (P1,'>',"part1_$fasta") or die "can't write output file $!\n";
open (P2,'>',"part2_$fasta") or die "can't write output file $!\n";

my $num = 1;
while (<FA>) {
	chomp;
	if ($_ =~ /(ACGTGGA).*(TCCACTGT)/) {
		my ($p_1,$p_2) = ($`,$');
		my $len1 = length($p_1);
		my $len2 = length($p_2);
		if ( $len1 > $len and $len2 > $len ) {
			$p_1 = substr($p_1,0,$len1-$p1c);
			$p_2 = substr($p_2,$p2c-1,$len2-$p2c);
			print P1 "\>$num\n$p_1\n";
			print P2 "\>$num\n$p_2\n";
			$num++;
		}
	}
}
close P1;
close P2;
system ("~/bin/fastx_reverse_complement -i part1_$fasta -o part1_RC_$fasta");
print "part1_$fasta is reversed!\n";
system ("~/bin/fastx_reverse_complement -i part2_$fasta -o part2_RC_$fasta");
print "part2_$fasta is reversed!\n";
$num--;
print "Program \"Separater\" on $fasta Completed!\nThere are $num results can be matched\n";
