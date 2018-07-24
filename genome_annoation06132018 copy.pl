#!/bin/env perl -w 
use strict; 

##################################################     transcript_End_analysis.pl      ####################################################################
######	compare the 3End sites in two samples
######	Run Sample: Bed, S5.sorted.bam
##################################################################################################################################################

#command line inputs
my ($bed_file_1, $bed_file_2, $bed_file_3) = @ARGV; #Open bed file and sam file


my $n;

my %strand_information;
my %gene_end_50_plus_inteval;
my @genome = 1..4649560;
my $genome;

foreach $genome (@genome ) {
	$strand_information{$genome}=["0", "", "","",""];
}


open (genome_bed,"./NC_000913.bed") or die "can't open bed file NC_000913.bed bed file\n";
my $previous_gene = "null:+";
my $previous_gene_end = 0;
my $next_gene_start = "0";
my $previous_gene_start= 0;
my $next_gene = "";
my $gap_start = 1;
my $gap_end = 4641652;
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
			my $postition = sprintf ("%.f", 100*($nucleotide-$gene_start+1)/($gene_end-$gene_start+1));
			my $distance_to_previous_gene_end = $nucleotide - $previous_gene_end  ;
			$strand_information{$nucleotide} =["\+", "$gene_name_s:+:$gene_start:$gene_end", "$postition","$previous_gene:$distance_to_previous_gene_end", ""];			
		}
		
	} elsif ($gene_direction =~ /-/)  {		
		for $nucleotide ($gene_start..$gene_end) {
			my $postition = sprintf ("%.f", 100*($gene_end-$nucleotide+1)/($gene_end-$gene_start+1));
			my $distance_to_previous_gene_end = $nucleotide - $previous_gene_end ;
			$strand_information{$nucleotide} =["-", "$gene_name_s:-:$gene_end:$gene_start", "$postition","$previous_gene:$distance_to_previous_gene_end", ""];		
		}		
	}
	
	
	my $nucletide_gap;
	foreach  $nucletide_gap ($gap_start..($gene_start-1) ) {
		if ($strand_information{$nucletide_gap}[0] =~ 0 ) {
			my $distance_to_previous_gene_end =  ($nucletide_gap - $previous_gene_end );
			my $postition = sprintf ("%.f", 100*($nucletide_gap -$previous_gene_end)/($gene_start-1-$previous_gene_end));
			my $gap_left_site = $previous_gene_end+1;
			my $gap_wright_site = $gene_start-1;
			my $distance_to_next_gene = $gene_start - $nucletide_gap;
			$strand_information{$nucletide_gap} = ["0", "gap:0:$gap_left_site:$gap_wright_site", "$postition","$previous_gene:$distance_to_previous_gene_end","$distance_to_next_gene:$gene_direction:$gene_name_s"]
		}
	}
	my $nucletide_previous;
	$previous_gene_end =0 if $gene_start == 190;
	foreach  $nucletide_previous ($previous_gene_start..$previous_gene_end) {
		my $distance_to_next_gene = $gene_start - $nucletide_previous;
		$strand_information{$nucletide_previous}[4]= "$distance_to_next_gene:$gene_direction:$gene_name_s";
	}
	
	$previous_gene ="$gene_name_s:$gene_direction";
	$previous_gene_end = $gene_end;
	$gap_start =  $gene_end+1;
	$previous_gene_start= $gene_start;
	}
	
	
	# to generate the next gene information
	
	
	
	
	foreach $genome (@genome ) {
		
	 print "$genome\t$strand_information{$genome}[0]\t$strand_information{$genome}[1]\t$strand_information{$genome}[2]\t$strand_information{$genome}[3]\t$strand_information{$genome}[4]\n";
		#print "$genome\t$strand_information{$genome}\t$gene_end_50_plus_inteval{$genome}\n";
	}
	