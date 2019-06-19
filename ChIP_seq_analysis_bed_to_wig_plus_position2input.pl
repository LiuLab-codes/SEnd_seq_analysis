#!/bin/env perl -w 
use strict; 





my ($bed_file_1, $bed_file_2) = @ARGV; #Open bed file and sam file

use File::Basename;
my $bed_name = basename ("$bed_file_2",  ".bedpe");
#my %start_count_positive;
my %total_signal_1;
my %total_signal_2;

my $genome;
my @genome = 1..4649560;
foreach $genome (@genome ) {
	$total_signal_1{$genome} =0;
	$total_signal_2{$genome} =0;
}



open (bedpe_file_1,"$bed_file_1") or die "can't open input file $!\n";  #input file



	while (<bedpe_file_1>) {
		chomp;
		my @original_mapping_bed_line = "$_";
		my @line_mapping = split(/\t/);
		my ($chrom_mapping,$read_start,$read_end,$read_name, $read_direction) = @line_mapping[0,1,2,3,5];
		my $nucleotide_mapping;
		for $nucleotide_mapping ($read_start..$read_end) {
				$total_signal_1{$nucleotide_mapping} +=1;
			    }
			}
close bedpe_file_1;

open (bedpe_file_2,"$bed_file_2") or die "can't open input file $!\n";  #Chip-seq file 



	while (<bedpe_file_2>) {
		chomp;
		my @original_mapping_bed_line = "$_";
		my @line_mapping = split(/\t/);
		my ($chrom_mapping,$read_start,$read_end,$read_name, $read_direction) = @line_mapping[0,1,2,3,5];
		my $nucleotide_mapping;
		for $nucleotide_mapping ($read_start..$read_end) {
				$total_signal_2{$nucleotide_mapping} +=1;
			    }
			}
close bedpe_file_2;



open (bed_wig,'>',"./$bed_name\_bed_2_wig_log2_input.wig") or die "can't write output file $!\n";
print bed_wig "track name=\"$bed_name\" color=98,0,234 altColor=255,0,0 graphType=bar viewLimits=0:1
fixedStep chrom=gi|556503834|ref|NC_000913.3| start=1 step=1\n";
open (bed_wig_txt,'>',"./$bed_name\_bed_2_wig_plus_position_log2input.txt") or die "can't write output file $!\n";
print bed_wig_txt "position\tinput_100nt\tIP_signal\t$bed_name\_ratio\n";

foreach $genome (@genome ) {
	
	
	my $input_value = 1;
	
	for (-50..50) {
		my $number = $_;
		$input_value += $total_signal_1{$genome+$number};
	}
	
	
	 
	my $ratio_value = ($total_signal_2{$genome}+1)/$input_value*101;
	
	print bed_wig_txt "$genome\t$total_signal_1{$genome}\t$total_signal_2{$genome}\t$ratio_value\n";
	print bed_wig "$ratio_value\n";
}

close bed_wig_txt;

close bed_wig;

