#!/bin/env perl -w 
use strict; 
##################################################     transcript_End_analysis.pl      ####################################################################
######	This program is designed to detect how many TSS before each annoated gene and how about the distance before translation site.
##################################################################################################################################################

#command line inputs
my ($bed_file_1) = @ARGV; #Open bed file and sam file

use File::Basename;
my $bed_name_1 = basename ("$bed_file_1",  ".bedpe");

unless (-d "./result_TSS_number/") {
	mkdir "./result_TSS_number/", 0755 or die "can't write output file $!\n";
}

my $n;




my %strand_information;
my %gene_positive_region;
my %gene_negative_region;
my %gene_5end_500_plus_inteval;
my %reads_1_positive_count;
my %reads_1_negative_count;

my %start_1_count_positive;
my %start_1_count_negative;

my %TSS_known_positive;
my %TSS_known_negative;
my @genome = 1..4649560;
my $genome;

foreach $genome (@genome ) {
	#$start_1_count_positive{$genome} =0;
	#$start_1_count_negative{$genome} =0;

	
	
	#$reads_1_positive_count{$genome} =0;
	#$reads_1_negative_count{$genome} =0;
	
	$gene_positive_region{$genome} =["0", "", "","",""];
	$gene_negative_region{$genome} =["0", "", "","",""];
	
	$strand_information{$genome}=0;
	$gene_5end_500_plus_inteval{$genome}=0;
	
	
	$TSS_known_positive{$genome}=["0", "", "","","","0"];
	$TSS_known_negative{$genome}=["0", "", "","","","0"];
	#
}

open (genome_bed,"./NC_000913.bed") or die "can't open bed file NC_000913.bed bed file\n";

while (<genome_bed>) {
	chomp;
	my @original_gene_bed_line = "$_";
	my @line = split(/\t/);
	my ($chrom,$gene_start,$gene_end,$gene_direction, $gene_name) = @line[0,3,4,6,8];
	my $nucleotide;
	my $nucleotide_inteval;
	my $gene_end_around;
	$gene_name =~ /\=([\s\S]+?)\;/;
	my $gene_name_s= $1;
	
	
	if ($gene_direction =~ /\+/ ) {
		 
		for $nucleotide ($gene_start..$gene_end) {
			#$gene_positive_region{$nucleotide} =["1", "\+", "$gene_start","$gene_end","$gene_name_s)"];
			
		}
		#print "1\t\+\t$gene_start\t$gene_end\t$gene_name_s\n";	
	} elsif ($gene_direction =~ /-/)  {
		
		for $nucleotide ($gene_start..$gene_end) {
			#$gene_negative_region{$nucleotide} =["-1", "\-", "$gene_start","$gene_end","$gene_name_s)"];
			
		}
		#print "-1\t\-\t$gene_start\t$gene_end\t$gene_name_s\n";	
	}
}
close genome_bed;
open (TSS_all,"$bed_file_1") or die "can't open bed file all tss file \n";
	while (<TSS_all>) {
		chomp;

	my @original_mapping_bed_line = "$_";

	my @line_mapping = split(/\t/);

	my ($TSS_site, $strand, $TSS_count, $TSS_left_signal, $TSS_right_signal) = @line_mapping[0,1,2,3,4];


	if ($strand =~ /\+/ && $TSS_site>0) {
		$TSS_known_positive{$TSS_site} =["$TSS_site","$strand", "$TSS_count","$TSS_left_signal", "$TSS_right_signal", "known"];
		
	} elsif ($strand =~ /-/ && $TSS_site>0) {
			$TSS_known_negative{$TSS_site} =["$TSS_site","$strand", "$TSS_count","$TSS_left_signal", "$TSS_right_signal", "known"];
			
	}

}


close TSS_all;

	my @TSS_500_head;
	open (genome_bed,"./NC_000913.bed") or die "can't open bed file NC_000913.bed bed file\n";
while (<genome_bed>) {
	chomp;
	my @original_gene_bed_line = "$_";
	my @line = split(/\t/);
	my ($chrom,$gene_start,$gene_end,$gene_direction, $gene_name) = @line[0,3,4,6,8];
	my $nucleotide;
	my $nucleotide_inteval;
	my $gene_end_around;
	$gene_name =~ /\=([\s\S]+?)\;/;
	my $gene_name_s= $1;
	
	if ($gene_direction =~ /\+/ ) {
		my $gene_start_500_head=0;
		
		my $TSS_500_head;
		 $gene_start >=300 ? ($gene_start_500_head=$gene_start-300) : ($gene_start_500_head=0);
		 my $count =0;
		 my $distance =0;
		  my $position =0;
		 foreach  ($gene_start_500_head..($gene_start+50)){
			  $position = $_;
			 
			if ($TSS_known_positive{$position}[0] !=0) {
				$count+=1;
				$distance =$gene_start-$position+1;
				#push (@TSS_500_head, "$gene_start:+:$count:$distance");
				#print "$gene_start\t\+\t$position\t$count\t$distance\n";
				
			} 
			
		 }  
      print "$gene_start\t$gene_name_s\t\+\t$position\t$count\t$distance\n";
	} elsif ($gene_direction =~ /\-/ ) {
		my $gene_start_500_head=0;
		
		my $TSS_500_head;
		 $gene_end < 4649060 ? ($gene_start_500_head=$gene_end+300) : ($gene_start_500_head=4649560);
		 my $count =0;
		 my $distance =0;
		 my $position =0;
		 foreach  (($gene_end-50)..$gene_start_500_head){
			  $position = $_;
			 
			if ($TSS_known_negative{$position}[0] !=0) {
				$count+=1;
				$distance =$position+1-$gene_end;
				#push (@TSS_500_head, "$gene_end:+:$count:$distance");
				#print "$gene_end\t\-\t$position\t$count\t$distance\n";
				
			} 
			
		 }
		 print "$gene_end\t$gene_name_s\t\-\t$position\t$count\t$distance\n";  
	}
}
close genome_bed;