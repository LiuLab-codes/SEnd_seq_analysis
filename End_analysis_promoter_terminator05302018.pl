#!/bin/env perl -w 
use strict; 

##################################################     transcript_End_analysis.pl      ####################################################################
######	Locate experiment data into reference genom.
######	Run Sample: Bed, S5.sorted.bam
##################################################################################################################################################

#command line inputs
my ($bed_file) = @ARGV; #Open bed file and sam file

use File::Basename;
my $bed_name = basename ("$bed_file",  ".bedpe");

unless (-d "./result/") {
	mkdir "./result/", 0755 or die "can't write output file $!\n";
}

my $n;

open (genome_bed,"./NC_000913.bed") or die "can't open bed file NC_000913.bed bed file\n";


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

my @genome = 1..4649560;
my $genome;
foreach $genome (@genome ) {
	#$gene_positive_region{$genome} = 1;
	#$gene_negative_region{$genome} =-1;
	#$gene_interval{$genome} = 0;
	
	$reads_positive_count{$genome} =0;
	$reads_negative_count{$genome} =0;
	
	
	
	$strand_information{$genome}=0;
	$gene_end_50_plus_inteval{$genome}=3;
	$reads_positive_negative_overlap{$genome} =0;
	
	$peak_start_positive_sharp{$genome} =0;
	$peak_end_positive_sharp{$genome} =0;
	$peak_start_negative_sharp{$genome} =0;
	$peak_end_negative_sharp{$genome} =0;
	
	#$reads_positive_end{$genome} =0;
	#$reads_negative_end{$genome} =0;
	#$reads_positive_start{$genome} =0;
	#$reads_negative_start{$genome} =0;
	#$n +=1;
	#print "line $n is $start_count{$genome} and the end is $end_count{$genome}\n";
}
while (<genome_bed>) {
	chomp;
	my @original_gene_bed_line = "$_";
	my @line = split(/\t/);
	my ($chrom,$gene_start,$gene_end,$gene_direction) = @line[0,3,4,6];
	my $nucleotide;
	my $nucleotide_inteval;
	my $gene_end_around;
	if ($gene_direction =~ /\+/ ) {
		 
		for $nucleotide ($gene_start..$gene_end) {
			$strand_information{$nucleotide} =1;
			
		}
		
		$gene_end_around = $gene_end -50;
		for $nucleotide_inteval ($gene_end_around ..$gene_end+50) {
			$gene_end_50_plus_inteval{$nucleotide_inteval} =0;
		}
		
		
		
	
	} elsif ($gene_direction =~ /-/)  {
		
		for $nucleotide ($gene_start..$gene_end) {
			$strand_information{$nucleotide} =-1;
			
		}
		$gene_end_around = $gene_start +50;
		for $nucleotide_inteval ($gene_start-50 ..$gene_end_around) {
			$gene_end_50_plus_inteval{$nucleotide_inteval} =0;
		}
	}
	
	}
	
	foreach $genome (@genome ) {
		
		if ($strand_information{$genome} ==0) {
			$gene_end_50_plus_inteval{$genome} =0
		}
		#print "$genome\t$strand_information{$genome}\t$gene_end_50_plus_inteval{$genome}\n";
	}

close genome_bed;



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
		
	
		    for $nucleotide_mapping ($read_start..$read_end) {
			$reads_positive_count{$nucleotide_mapping} +=1;
		    }
		} elsif ($read_direction =~ /-/ ) {
			
			
			#$reads_negative_end{$read_start} +=1;
			
			#$reads_negative_start{$read_end} +=1;
			
	       for $nucleotide_mapping ($read_start..$read_end) {
			$reads_negative_count{$nucleotide_mapping} +=1;
	}
	
}
}



my %start_count_positive;
my %end_count_positive;
my %start_count_negative;
my %end_count_negative;

foreach $genome (@genome ) {
	$start_count_positive{$genome} = 0;
	$end_count_positive{$genome} =0;
	$start_count_negative{$genome} = 0;
	$end_count_negative{$genome} =0;
	#$n +=1;
	#print "line $n is $start_count{$genome} and the end is $end_count{$genome}\n";
}
close mapping_bed;

open (mapping_bed,"./$bed_file") or die "can't open bed file $bed_file\n";
while (<mapping_bed>) {
	chomp;

	my @original_line_mapping_bed_peak = "$_";
	my @line_mapping_bed_peak = split(/\t/);
	my ($chrom_mapping_peak,$gene_start_peak,$gene_end_peak,$reads_name_peak,$gene_direction_peak) = @line_mapping_bed_peak[0,1,2,3,5];
	#print "start site is $gene_start and the end site is $gene_end and the direction is $gene_direction \n";
    if ($gene_direction_peak =~ /\+/ ) {
    	
		$start_count_positive{$gene_start_peak} +=1;
		$end_count_positive{$gene_end_peak} +=1;
    } elsif ($gene_direction_peak =~ /-/) {
		$start_count_negative{$gene_end_peak} +=1;
		$end_count_negative{$gene_start_peak}  +=1;
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
	if ($start_count_positive{$genome} > 10 && $genome >$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome +$next;
			$max = $now  if $start_count_positive{$now} > $start_count_positive{$max};
		}
		$peak_0 = $genome;
		push(@peak_start_positive,$max)
		
	}
	
}
#for the end site of positive strand filter
$peak_0 =0;
$sum =0 ;
foreach $genome (@genome ) {
	
	
	$sum = $peak_0 +10;
	if ($end_count_positive{$genome} >=5 && $genome >$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome +$next;
			$max = $now  if $end_count_positive{$now} > $end_count_positive{$max};
		}
		$peak_0 = $genome;
		push(@peak_end_positive,$max)
		
	}
	
}

#for the start site of negative strand filter
$peak_0 =0;

foreach $genome (@genome ) {
	
	
	$sum = $peak_0 +10;
	if ($start_count_negative{$genome} > 10 && $genome >$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome +$next;
			$max = $now  if $start_count_negative{$now} > $start_count_negative{$max};
		}
		$peak_0 = $genome;
		push(@peak_start_negative,$max)
		
	}
	
}

#for the end site of negative strand filter
$peak_0 =0;

foreach $genome (@genome ) {
	
	
	$sum = $peak_0 +10;
	if ($end_count_negative{$genome} >3 && $genome >$sum) {
		
		
		 $max =$genome;
		
		for (0..10) {
			 my $next = $_;
			 $now = $genome +$next;
			$max = $now  if $end_count_negative{$now} > $end_count_negative{$max};
		}
		$peak_0 = $genome;
		push(@peak_end_negative,$max)
		
	}
	
}


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
open (Termi,'>',"./result/$bed_name\_terminator.txt") or die "can't write output file overlap.txt\n";





foreach $peak_start_positive (@peak_start_positive) {
	my $peak_start_positive_sum_left = 0;
	my $peak_start_positive_sum_right =0;
	for (1..10) {
		my $number = $_;
		$peak_start_positive_sum_left += $reads_positive_count{$peak_start_positive-$number};
		$peak_start_positive_sum_right += $reads_positive_count{$peak_start_positive+$number-1};
	}
	if ($peak_start_positive_sum_right > (1.2*$peak_start_positive_sum_left) or $start_count_positive{$peak_start_positive} >30 ) {
		push (@peak_start_positive_soaring, $peak_start_positive);
		$peak_start_positive_sharp{$peak_start_positive} =$start_count_positive{$peak_start_positive};
		print TSS "$peak_start_positive\t\+\t$start_count_positive{$peak_start_positive}\t$peak_start_positive_sum_left\t$peak_start_positive_sum_right\n";
	}
		
}	
	
foreach $peak_end_positive (@peak_end_positive) {
	my $peak_end_positive_sum_left = 0;
	my $peak_end_positive_sum_right =0;
	for (1..10) {
		my $number = $_;
		$peak_end_positive_sum_left += $reads_positive_count{$peak_end_positive-$number+1};
		$peak_end_positive_sum_right += $reads_positive_count{$peak_end_positive+$number};
	}
	if ($peak_end_positive_sum_left > (1.2*$peak_end_positive_sum_right)  or $end_count_positive{$peak_end_positive} >30 ) {
		push (@peak_end_positive_soaring, $peak_end_positive);
		$peak_end_positive_sharp{$peak_end_positive} = $end_count_positive{$peak_end_positive};
		print Termi "$peak_end_positive\t\+\t$end_count_positive{$peak_end_positive}\t$peak_end_positive_sum_left\t$peak_end_positive_sum_right\n";
}
}


foreach $peak_start_negative (@peak_start_negative) {
	my $peak_start_negative_sum_left = 0;
	my $peak_start_negative_sum_right =0;
	for (1..10) {
		my $number = $_;
		$peak_start_negative_sum_left += $reads_negative_count{$peak_start_negative-$number+1};
		$peak_start_negative_sum_right += $reads_negative_count{$peak_start_negative+$number};
	}
	
	if ($peak_start_negative_sum_left  > (1.2*$peak_start_negative_sum_right) or $start_count_negative{$peak_start_negative} >30 ) {
		push (@peak_start_negative_soaring, $peak_start_negative );
		$peak_start_negative_sharp{$peak_start_negative} = $start_count_negative{$peak_start_negative};
		print TSS "$peak_start_negative\t\-\t$start_count_negative{$peak_start_negative}\t$peak_start_negative_sum_right \t$peak_start_negative_sum_left\n";
	
	}
}

foreach $peak_end_negative (@peak_end_negative) {
	my $peak_end_negative_sum_left = 0;
	my $peak_end_negative_sum_right =0;
	for (1..10) {
		my $number = $_;
		$peak_end_negative_sum_left += $reads_negative_count{$peak_end_negative-$number};
		$peak_end_negative_sum_right += $reads_negative_count{$peak_end_negative+$number+1};
	}
	
	if ((1.2*$peak_end_negative_sum_left)  < $peak_end_negative_sum_right or $end_count_negative{$peak_end_negative} >30) {
		push (@peak_end_negative_soaring, $peak_end_negative );
		$peak_end_negative_sharp{$peak_end_negative} = $end_count_negative{$peak_end_negative-1};
		print Termi "$peak_end_negative\t\-\t$end_count_negative{$peak_end_negative-1}\t$peak_end_negative_sum_left\t$peak_end_negative_sum_right\t$peak_end_negative_sharp{$peak_end_negative}\t$reads_negative_count{$peak_end_negative}\n";
	
	}
}

close TSS;
close Termi;

open (Overlap,'>',"./result/$bed_name\_overlap.txt") or die "can't write output file overlap.txt\n";
open (Overlap_sequence,'>',"./result/$bed_name\_overlap_sequence.fa") or die "can't write output file overlap.txt\n";

my $peak_end_positive_soaring;
my $peak_end_negative_soaring;
my $overlap_count =1;
my $previous_overlap_start=0;
foreach $peak_end_positive_soaring (@peak_end_positive_soaring) {
	
	if ($reads_negative_count{$peak_end_positive_soaring}>3  ) {
		my $overlap_start_site;
		my $overlamp_end_site = $peak_end_positive_soaring;
		for (10..90) {
			my $number_back = $_;
			if ($peak_end_negative_sharp{($peak_end_positive_soaring-$number_back)} >0 && ($peak_end_positive_soaring-$number_back) > $previous_overlap_start) {
				$overlap_start_site = ($peak_end_positive_soaring - $number_back) ;
				my $length_overlap = $overlamp_end_site-$overlap_start_site+1;
				
				my $overlap_positive_nucletide_count=0;
				my $overlap_negative_nucletide_count=0;
				my $overlap_unannoted_nucletide_count=0;
				for ($overlap_start_site..$overlamp_end_site){
					if ($strand_information{$_} ==0) {
						$overlap_unannoted_nucletide_count +=1;
						
					} elsif ($strand_information{$_} ==1) {
						$overlap_positive_nucletide_count +=1;
					} elsif ($strand_information{$_} ==-1) {
						$overlap_negative_nucletide_count +=1;
					}
				}
				
				my $positive_signal_sum=0;
				my $negative_signal_sum=0;
				my $ration =0;
				for ($overlap_start_site..$overlamp_end_site){
					my $signal_number = $_;
					 $positive_signal_sum += $reads_positive_count{$signal_number};
					$negative_signal_sum += $reads_negative_count{$signal_number};
					$ration= $positive_signal_sum/$negative_signal_sum;
				}
				
				print  Overlap "$overlap_count\t$overlap_start_site\t$overlamp_end_site\t$length_overlap\t$overlap_positive_nucletide_count\t$overlap_negative_nucletide_count\t$overlap_unannoted_nucletide_count\n";
				my $overlap_sequence = &get_genome_sequence($overlap_start_site, $overlamp_end_site);
				print Overlap_sequence "\>$overlap_count\_$length_overlap\_$overlap_start_site\_$overlamp_end_site\n$overlap_sequence\n";
				print "$overlap_count\t$overlap_start_site\t$overlamp_end_site\t$length_overlap\t$positive_signal_sum\t$negative_signal_sum\t$ration\n";
				#print "$length_overlap\n";
				$overlap_count +=1;
				$previous_overlap_start = $overlap_start_site;
			}
			
		}
		
	}
	
	
	
	
	
	
}

close Overlap;
close Overlap_sequence;




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





	















