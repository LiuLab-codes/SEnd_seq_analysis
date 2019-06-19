#!/bin/env perl -w 
use strict; 

# This script is used to extract TSSs from the bedfile. 
my ($bed_file, $bam_file) = @ARGV; #Open bed file and sam file

use File::Basename;
my $bed_name = basename ("$bed_file",  ".bedpe");
my $bam_file_name = basename ("$bam_file",  ".bam");
unless (-d "./result/") {
	mkdir "./result/", 0755 or die "can't write output file $!\n";
}
my $n;

my %strand_information;
my %gene_positive_region;
my %gene_negative_region;
my %gene_end_50_plus_inteval;
my %reads_positive_count;
my %reads_negative_count;
my %reads_positive_negative_overlap;
my %peak_start_positive_sharp;
my %peak_end_positive_sharp;
my %peak_start_negative_sharp;
my %peak_end_negative_sharp;
my %start_count_positive;
my %end_count_positive;
my %start_count_negative;
my %end_count_negative;

my @genome = 1..4649560;
my $genome;

foreach $genome (@genome ) {
	$reads_positive_count{$genome} =0;
	$reads_negative_count{$genome} =0;
	#$strand_information{$genome}=0;
	#$gene_end_50_plus_inteval{$genome}=3;
	#$reads_positive_negative_overlap{$genome} =0;
	$peak_start_positive_sharp{$genome} =0;
	#$peak_end_positive_sharp{$genome} =0;
	$peak_start_negative_sharp{$genome} =0;
	#$peak_end_negative_sharp{$genome} =0;
	$start_count_positive{$genome} = 0;
	#$end_count_positive{$genome} =0;
	$start_count_negative{$genome} = 0;
	#$end_count_negative{$genome} =0;
	
}




open (mapping_bed,"./$bed_file") or die "can't open bed file $bed_file\n";



	while (<mapping_bed>) {
		chomp;
	
	my @original_mapping_bed_line = "$_";
	my @line_mapping = split(/\t/);
	
	my ($chrom_mapping,$read_start,$read_end,$read_name, $read_direction) = @line_mapping[0,1,2,3,5];
	my $nucleotide_mapping;
	
	if ($read_direction =~ /\+/ ) {
		
		#$reads_positive_end{$read_end} +=1;
		
		#$reads_positive_start{$read_start} +=1;
		$start_count_positive{$read_start} +=1;
		$end_count_positive{$read_end} +=1;
	
		    for $nucleotide_mapping ($read_start..$read_end) {
			$reads_positive_count{$nucleotide_mapping} +=1;
		    }
		} elsif ($read_direction =~ /-/ ) {
			
			$start_count_negative{$read_end} +=1;
			$end_count_negative{$read_end}  +=1;
			#$reads_negative_end{$read_start} +=1;
			
			#$reads_negative_start{$read_end} +=1;
			
	       for $nucleotide_mapping ($read_start..$read_end) {
			$reads_negative_count{$nucleotide_mapping} +=1;
	}
	
}
}

close mapping_bed;
my @peak_start_positive;
my @peak_end_positive;
my @peak_start_negative;
my @peak_end_negative;
my $peak_0 =0;
my $sum;
my $now;
my $max;

#for the start site of positive strand filter
foreach $genome (@genome ) {

	$sum = $peak_0 +10;
	if (($start_count_positive{$genome} +$start_count_positive{$genome+1} + $start_count_positive{$genome+2}) >=10 && $genome >$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome +$next;
			$max = $now  if $start_count_positive{$now} > $start_count_positive{$max};
		}
		$peak_0 = $max;
		
		push(@peak_start_positive,$max);
		
	}
	
}
#for the end site of positive strand filter
$peak_0 =0;
$sum =0 ;


#for the start site of negative strand filter
$peak_0 =4649570;

my @genome_r = reverse @genome;
foreach $genome (@genome_r) {
	
	$sum = $peak_0 -10;
	if (($start_count_negative{$genome} +$start_count_negative{$genome+1} + $start_count_negative{$genome+2}) >=10 && $genome <$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome -$next;
			$max = $now  if $start_count_negative{$now} > $start_count_negative{$max};
		}
		$peak_0 = $max;
		push(@peak_start_negative,$max)
		
	}
	
}

#for the end site of negative strand filter
$peak_0 =0;




my $n_peak_start_positive = @peak_start_positive;
my $n_peak_end_positive = @peak_end_positive;
my $n_peak_start_negative = @peak_start_negative;
my $n_peak_end_negative =@peak_end_negative;
print "$n_peak_start_positive, $n_peak_end_positive, $n_peak_start_negative, $n_peak_end_negative\n";

my @peak_start_positive_soaring;
my @peak_end_positive_soaring;
my @peak_start_negative_soaring;
my @peak_end_negative_soaring;
my $peak_start_positive;
my $peak_end_positive;
my $peak_start_negative;
my $peak_end_negative;

open (TSS,'>',"./result/$bed_name\_TSS.txt") or die "can't write output file overlap.txt\n";






foreach my $end_site (@peak_start_positive) {
	
	
	my $direction = "\+";
	my $end_site_analyis_result = &end_site_analysis_positive($end_site, $direction );

	my @var_positive = split (/\t/, $end_site_analyis_result);
	pop @var_positive;
	my %positive_gene_coverage;
	my ($positve_coverage_left,$positve_coverage_right ) = (($end_site-1),($end_site+100));
	for my $read_infor_s (@var_positive) {
		#print "1 is $read_infor_s is 1\n";
			my @read_infor = split (/:/, $read_infor_s) ;
		if ($read_infor[2] =~ /\+/  ){
		
			#$positive_gene_coverage{$read_infor[0]}[1]+=1;
			for my $read_position (($read_infor[0])..$read_infor[1]){
				$positive_gene_coverage{$read_position}[0] +=1 ;
				#$positive_gene_coverage{$read_position}[0] =1 unless $positive_gene_coverage{$read_position}[0]=~ /^[0-9,.E]+$/;
			}
			$positve_coverage_left = $read_infor[0] if ($read_infor[0] < $positve_coverage_left);
			$positve_coverage_right = $read_infor[1] if ($read_infor[1] > $positve_coverage_right);
			}
	
		}
		my $peak_end_positive_sum_left = 0;
		my $peak_end_positive_sum_right =0;
		for (0..10) {
			my $number = $_;
			$peak_end_positive_sum_left += $positive_gene_coverage{$end_site-$number-1}[0] if $positive_gene_coverage{$end_site-$number-1}[0] ;
			$peak_end_positive_sum_right += $positive_gene_coverage{$end_site+$number}[0] if $positive_gene_coverage{$end_site+$number}[0];
		}
	
		my $total_5Ends_count_positive = 0;
	
		for (-5..5) {
			my $number= $_;
			$total_5Ends_count_positive += $start_count_positive{$end_site+$number} if $start_count_positive{$end_site+$number} >0 ;
		}
	
	
		if ($peak_end_positive_sum_right > (1.3*$peak_end_positive_sum_left)) {
			print TSS "$end_site\t\+\t$total_5Ends_count_positive\t$peak_end_positive_sum_left\t$peak_end_positive_sum_right\t$bam_file_name\n";
		}
}	
	



foreach my $end_site (@peak_start_negative) {
	
	my $direction = "\-";
	my $end_site_analyis_result = &end_site_analysis_negative($end_site, $direction );
	my @var_negative = split (/\t/, $end_site_analyis_result);
	pop @var_negative;
	my %negative_gene_coverage;
	my ($negative_coverage_left,$negative_coverage_right ) = (($end_site-100),($end_site+1));
	for my $read_infor_s (@var_negative) {
			my @read_infor = split (/:/, $read_infor_s) ;
		if ($read_infor[2] =~ /-/  ){
		
			#$negative_gene_coverage{$read_infor[0]}[1]+=1;
			for my $read_position ($read_infor[0]..($read_infor[1])){
				$negative_gene_coverage{$read_position}[0] +=1 ;
				#$negative_gene_coverage{$read_position}[0] =1 unless $negative_gene_coverage{$read_position}[0]=~ /^[0-9,.E]+$/;
			}
			$negative_coverage_left = $read_infor[0] if ($read_infor[0] < $negative_coverage_left);
			$negative_coverage_right = $read_infor[1] if ($read_infor[1] > $negative_coverage_right);
			}
	
		}
		my $peak_end_negative_sum_upstream = 0;
		my $peak_end_negative_sum_downstream =0;
		for (0..10) {
			my $number = $_;
			$peak_end_negative_sum_upstream += $negative_gene_coverage{$end_site+$number+1}[0] if $negative_gene_coverage{$end_site+$number+1}[0] ;
			$peak_end_negative_sum_downstream += $negative_gene_coverage{$end_site-$number}[0] if $negative_gene_coverage{$end_site-$number}[0];
		}
	
		my $total_5Ends_count_negative;
	
		for (-5..5) {
			my $number= $_;
			$total_5Ends_count_negative += $start_count_negative{$end_site+$number} if $start_count_negative{$end_site+$number} > 0  ;
		}
	
	
		if ($peak_end_negative_sum_downstream > (1.3*$peak_end_negative_sum_upstream)  ) {
			print TSS "$end_site\t-\t$total_5Ends_count_negative\t$peak_end_negative_sum_upstream\t$peak_end_negative_sum_downstream\t$bam_file_name\n";
		}
	
	
}


close TSS;





sub get_genome_sequence {
	
	my $extract_start = shift (@_);
	
	my $extract_end = shift (@_);
	
	my $extract_length = $extract_end -$extract_start +1;
	$extract_start -=1;
	open(genome_sequence_fa, "./NC_000913_spikeRNA.fna")   or die "Could not open file genome file $!";

	my $genome_sequence_input;
	my $count_sub = 1;
	while (<genome_sequence_fa>) {
	  if ($count_sub == 2) {
	    #print line, then get out of here
	    $genome_sequence_input =$_;
		chomp $genome_sequence_input;
	   
	    
	  }
	  $count_sub +=1;
	}
	my $target_sequence = substr($genome_sequence_input, $extract_start , $extract_length);
	return $target_sequence;
	 close genome_sequence_fa;
}



sub end_site_analysis_positive {
	my $site = shift (@_);
	my $direction = shift (@_);
	if ($direction =~ /\+/ ){
		my $site_check_start = $site+60;
		my $site_check_end = $site+100;
		
		my $bed_each_gene = "\"gi|556503834|ref|NC_000913.3|:$site_check_start-$site_check_end\"";
		my $var_left_original = `samtools view  -b $bam_file $bed_each_gene | bedtools bamtobed -i stdin |cut -f 2,3,6  | grep \+ |tr "\t" ":"|tr "\n" "\t" `;
		chomp $var_left_original;
		return "$var_left_original\n";
	}
}

sub end_site_analysis_negative {
	my $site = shift (@_);
	my $direction = shift (@_);
	if ($direction =~ /\-/ ){
		my $site_check_start = $site-100;
		my $site_check_end = $site-60;
		
		my $bed_each_gene = "\"gi|556503834|ref|NC_000913.3|:$site_check_start-$site_check_end\"";
		my $var_left_original = `samtools view  -b $bam_file $bed_each_gene | bedtools bamtobed -i stdin |cut -f 2,3,6  | grep \- |tr "\t" ":"|tr "\n" "\t" `;
		chomp $var_left_original;
		return "$var_left_original\n";
	}
}

	















