#!/bin/env perl -w 
use strict; 

##################################################     3A_5A analysis.pl      ####################################################################
######	Locate experiment data into reference genom.
######	Run Sample: E.coli_K12_MG1655.gtf--gene_file, S5.sam--exp_file
##################################################################################################################################################

#command line inputs
my ($gene_file,$exp_file) = @ARGV;

#collect experiment and genom data
my ($forward_gene_start_array_ref,$forward_gene_end_array_ref,$forward_gene_name_array_ref,$forward_gene_len,$backward_gene_start_array_ref,$backward_gene_end_array_ref,$backward_gene_name_array_ref,$backward_gene_len) = &get_gene($gene_file);
my ($read_name_array_ref,$read_start_array_ref,$read_end_array_ref,$read_distance_array_ref,$read_strand_array_ref,$exp_file_len) = &get_exp($exp_file);

#open output file
open (OUTPUT,'>',"./output") or die "can't write output file ...\n";

#map experiment data to genom data
for (my $i=0; $i<$exp_file_len/2; $i++){	
	my ($start_match,$end_match,$end_match_value,$gene_strand) = (0,0,0,0);
	my $start_temp = 0;
	my ($start,$end) = (0,0);
	my ($gene_start,$gene_end) = (0,0);
	my $cov = 0;
	my @cov_gene;
	if (${$read_start_array_ref}[2*$i] < ${$read_start_array_ref}[2*$i+1]) {
		$start = ${$read_start_array_ref}[2*$i];
		$end = ${$read_start_array_ref}[2*$i+1];
	}
	else {
		$start = ${$read_start_array_ref}[2*$i+1];
		$end = ${$read_start_array_ref}[2*$i];
	}
	next unless (abs($start - $end) < 10000);
	for (my $j=0; $j<$forward_gene_len; $j++) {
		my $com_sta;
		if ($j eq 0) {
			$com_sta = 0;
		}
		else {
			$com_sta = ${$forward_gene_end_array_ref}[$j-1];
		}					
		my $com_end = ${$forward_gene_end_array_ref}[$j];
		next unless ($com_sta < $start and $start < $com_end);
		$start_match = $start - ${$forward_gene_start_array_ref}[$j];
		$gene_start = ${$forward_gene_start_array_ref}[$j];
		$gene_strand = '+';
		$cov = 1;
		@cov_gene = ();
		for (my $k = $j; $k <= $forward_gene_len; $k++) {
			if ($end <= ${$forward_gene_end_array_ref}[$k]) {
				$end_match_value = $end - ${$forward_gene_end_array_ref}[$k];
				$end_match = "End_IN_gene";
				push @cov_gene,${$forward_gene_name_array_ref}[$k];
				$gene_end = ${$forward_gene_end_array_ref}[$k];
				last;
			}
			elsif (${$forward_gene_end_array_ref}[$k] < $end and $end < ${$forward_gene_start_array_ref}[$k+1]) {
				$end_match_value = $end - ${$forward_gene_end_array_ref}[$k];
				$end_match = "End_BETWEEN_gene";
				push @cov_gene,(${$forward_gene_name_array_ref}[$k],"GAP");
				$gene_end = ${$forward_gene_end_array_ref}[$k];
				last;
			}
			else {
				push @cov_gene,(${$forward_gene_name_array_ref}[$k],"GAP ;");
				$cov++;
			}
		}
	}
	($start,$end) = ($end,$start);
	for (my $j=0; $j<$backward_gene_len;$j++) {
		my $com_sta = ${$backward_gene_end_array_ref}[$j];
		my $com_end;
		if ($j eq $backward_gene_len-1) {
			$com_end = 0;
		}
		else {
			$com_end = ${$backward_gene_end_array_ref}[$j+1];
		}
		next unless ($com_sta < $start and $start < $com_end);
		$start_temp = $start - ${$backward_gene_start_array_ref}[$j];
		if (abs($start_temp) < abs($start_match)){
			$start_match = $start_temp;
			$gene_start = ${$backward_gene_start_array_ref}[$j];
			$gene_strand = '-';
			$cov = 1;
			@cov_gene = ();
			for (my $k = $j; $k >= 0; $k--) {
				if ($end >= ${$backward_gene_end_array_ref}[$k]) {
					$end_match_value = $end - ${$backward_gene_end_array_ref}[$k];
					$end_match = "End_IN_gene";
					push @cov_gene,${$backward_gene_name_array_ref}[$k];
					$gene_end = ${$backward_gene_end_array_ref}[$k];
					last;
				}
				elsif (${$backward_gene_start_array_ref}[$k-1] < $end and $end < ${$backward_gene_end_array_ref}[$k]) {
					$end_match_value = $end - ${$backward_gene_end_array_ref}[$k];					
					$end_match = "End_BETWEEN_gene";
					push @cov_gene,(${$backward_gene_name_array_ref}[$k],"GAP");
					$gene_end = ${$backward_gene_end_array_ref}[$k];
					last;
				}
				else {
					push @cov_gene,(${$backward_gene_name_array_ref}[$k],"GAP ;");
					$cov++;
				}
			}
		}
	}
	if ($gene_strand eq '+'){
		($start,$end) = ($end,$start);
	}
	print OUTPUT "${$read_name_array_ref}[2*$i]\t$start\t$end\t$gene_strand\t$start_match\t$end_match_value\t$end_match\t$cov\t$gene_start\t@cov_gene\t$gene_end\n";
}
close OUTPUT;
print "Program Complete!\n";	
	
	
###########################################################################################
############################# Additional Functions ########################################
###########################################################################################

#extract gene data
sub get_gene {
	my ($gene_file) = @_;
	my @forward_gene_start_array;
	my @forward_gene_end_array;
	my @forward_gene_name_array;
	my @backward_gene_start_array;
	my @backward_gene_end_array;
	my @backward_gene_name_array;
	my $c = 1;
	my $f = 0;
	my $b = 0;
	open(GENE,"./$gene_file") or die "can't open gene file $gene_file\n";
	open (FORWARD,'>',"./output_forwardgene") or die "can't write output file output_forwardgene\n";
	open (BACKWARD,'>',"./output_backwardgene") or die "can't write output file output_backwardgene\n";
	while (<GENE>) {
		chomp;
		my @line = split(/\t/);
		my ($style,$gene_start,$gene_end,$gene_strand,$gene_name) = @line[2,3,4,6,8];
		next unless ($style eq 'exon');
		if ($gene_strand eq '+') {
			push @forward_gene_start_array,$gene_start;
			push @forward_gene_end_array,$gene_end;
			push @forward_gene_name_array,$gene_name;
			print FORWARD "$c -> $gene_start\t$gene_end\t$gene_strand\t$gene_name\n";
			$f++;
		}
		else {
			push @backward_gene_start_array,$gene_end;
			push @backward_gene_end_array,$gene_start;
			push @backward_gene_name_array,$gene_name;
			print BACKWARD "$c -> $gene_end\t$gene_start\t$gene_strand\t$gene_name\n";
			$b++;
		}
		$c++;
	}
	close FORWARD;
	close BACKWARD;
	return (\@forward_gene_start_array,\@forward_gene_end_array,\@forward_gene_name_array,$f,\@backward_gene_start_array,\@backward_gene_end_array,\@backward_gene_name_array,$b);
}

#extract experiment data and earse invaild data
sub get_exp {
	my ($exp_file) = @_;
	my @read_name_array;
	my @read_start_array;
	my @read_end_array;
	my @read_distance_array;
	my @read_strand_array;
	my $read_strand;
	my $c = 0;
	open(EXP,"./$exp_file") or die "can't open experiment file $exp_file\n";
	open (OUTPUT,'>',"./output_exp") or die "can't write output file output_exp\n";
	while (<EXP>) {
		chomp;
		my @line = split(/\t/);
		next unless ($#line > 8);	
		my ($read_name,$read_start,$read_end,$read_distance) = @line[0,3,7,8];	
		next unless ($read_distance != 0);			
		$read_strand = ($read_distance > 0) ? '+':'-';		
		print OUTPUT "$c -> $read_name\t$read_start\t$read_end\t$read_distance\t$read_strand\n";
		push @read_name_array,$read_name;
		push @read_start_array,$read_start;
		push @read_end_array,$read_end;
		push @read_distance_array,$read_distance;
		push @read_strand_array,$read_strand;
		$c++;
	}
	close OUTPUT;
	return (\@read_name_array,\@read_start_array,\@read_end_array,\@read_distance_array,\@read_strand_array,$c);
}


			
		
	
	
	
	
	
	
	
	
	