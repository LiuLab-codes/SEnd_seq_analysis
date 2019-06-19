#!/bin/env perl -w 
use strict; 


my ($TTS_file) = @ARGV; #Open bed file and sam file

use File::Basename;
my $TTS_file_name = basename ("$TTS_file",  ".txt");


my %TTS_known;

	
open (result_analysis,'>',"./$TTS_file_name\_annoation.txt") or die "can't write output file $!\n";


open (TTS_all,"$TTS_file") or die "can't open bed file all TTs file \n";
			while (<TTS_all>) {
				chomp;
			my @original_mapping_bed_line = "$_";
			my @line_mapping = split(/\t/);
			my ($TTS_n, $TTS_site, $strand, $log_infor, $stationary_infor) =[" ", " ", " ", " ", " "];
			($TTS_n, $TTS_site, $strand, $log_infor, $stationary_infor) = @line_mapping[0,1,2,3,4];
			#print "$TTS_site \t 5\t$strand, $log_infor, $stationary_infor\n";
			next unless $line_mapping[0]=~ /^[0-9,.E]+$/;
	  if ($strand =~ /\+/ ) {
		   
		my $sequence_3End = &get_genome_sequence(($TTS_site-50), ($TTS_site+15));
		my $end_sequence_analysis = &TTS_sequence_analysis($sequence_3End);
		print result_analysis "$TTS_n\t$TTS_site\t$strand\t$end_sequence_analysis\t$log_infor\t$stationary_infor\n";
	
	} elsif ($strand =~ /-/ && $TTS_site>0) {
				
				my $sequence_3End_n = &get_genome_sequence(($TTS_site-15), ($TTS_site+50));
				my $sequence_3End_revcom = &reverse_complement_IUPAC($sequence_3End_n);
		
				my $end_sequence_analysis_n = &TTS_sequence_analysis($sequence_3End_revcom);
				
				print result_analysis "$TTS_n\t$TTS_site\t$strand\t$end_sequence_analysis_n\t$log_infor\t$stationary_infor\n";
				
		}
		   
		   
		   
}		   
close result_analysis;



sub TTS_sequence_analysis {
	my ($seq_DNA_o)= @_;
	my $seq_DNA_45 = substr ($seq_DNA_o, 9, 45);
	#print "original_trim\t$seq_DNA_45\n";
	my $DNA_length_seq_DNA = length ($seq_DNA_45);
	my %nucletide_fold;
	
	
	my @rna_fold_output_seq_DNA = `echo $seq_DNA_45| RNAfold -p`;
	my $rna__seq_DNA = $rna_fold_output_seq_DNA[0]; chomp $rna__seq_DNA;

	
	$rna_fold_output_seq_DNA[1] =~ /^(.+)\s+\((.+)\)/;
	my ($fold_seq_DNA,$energy_seq_DNA) = ($1,$2); 
	
	my $seq_DNA=$seq_DNA_45;
	if ($fold_seq_DNA =~ m/\)\.{0,8}+\(/) {
		print "left ($`) ($&) right($')\n";
		my ($left_fold, $middle, $right_fold) = ($`,$&,$');
		my $left_folding_length = () = $left_fold =~ /\(/g;
		my $right_folding_length = () = $right_fold =~ /\(/g;
		my $right_folding_position = index ($fold_seq_DNA, $middle);
		if ($left_folding_length <=$left_folding_length) {
		
		#my $seq_DNA = substr($seq_DNA_o, (9+$right_folding_position-$left_folding_length), 45);
		 $seq_DNA = substr($seq_DNA_45, ($right_folding_position-$left_folding_length));
		 $DNA_length_seq_DNA = length ($seq_DNA);
	 } else {
	 	
		 $seq_DNA = substr($seq_DNA_45,0, ($right_folding_position+$right_folding_length));
		$DNA_length_seq_DNA = length ($seq_DNA);
	 }
	 
	 
		#print "$left_folding_length\t$right_folding_length\t$right_folding_position\t$seq_DNA\t$DNA_length_seq_DNA\n";;
		
	}
	
	my @seq_DNA_base= split(//,$seq_DNA); 
	
	
	my @rna_fold_output_seq_DNA_after_2fold = `echo $seq_DNA| RNAfold -p`;
	my $rna__seq_DNA_after_2fold = $rna_fold_output_seq_DNA_after_2fold[0]; chomp $rna__seq_DNA_after_2fold;

	
	$rna_fold_output_seq_DNA_after_2fold[1] =~ /^(.+)\s+\((.+)\)/;
	my ($fold_seq_DNA_45,$energy_seq_DNA_45) = ($1,$2); 
	
	my @RNA_fold_base_seq_DNA= split(//,$fold_seq_DNA_45);
	
	
	
	for (0..($DNA_length_seq_DNA-1)){
		my $sequence_position = $_;
		$nucletide_fold{$sequence_position}[0] = $seq_DNA_base[$sequence_position];
		$nucletide_fold{$sequence_position}[1] = $RNA_fold_base_seq_DNA[$sequence_position];
		#print "$sequence_position\t$nucletide_fold{$sequence_position}[0]\t$nucletide_fold{$sequence_position}[1]\n";
		}
	
	
	
		
	
	
	my $new_sequence;
	my $add_tail;
	
	my $add_tail_dot;
	
	my $seq_length = length($seq_DNA);
	my $later = substr($seq_DNA, -13, 13);
	my $first_T =index ($later, "TTT");
	my $position_TTT = $seq_length-13+$first_T;
	my $U_tract_sequence = substr($later, $first_T-2, 13);
	
	#print "U_tract\t$U_tract_sequence\n";
	my $U_tract_count_to_RNAfold =0;
	if ($first_T != -1) {
		for (($position_TTT-1)..($DNA_length_seq_DNA-1)){
			my $U_position =$_;
		$U_tract_count_to_RNAfold ++ if (($nucletide_fold{$U_position}[0] eq "T") && ($nucletide_fold{$U_position}[1] eq "\)"));
		
	}
		
	#print 	"Utract\t$first_T\t\t$position_TTT\t$U_tract_count_to_RNAfold\n";
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	my $concetive_T_site= ($first_T==-1)? $seq_length : ($DNA_length_seq_DNA-13+$first_T);
	my $last_fold_r =rindex ($fold_seq_DNA, "))");
	
	
	
	my $warning = " ";
	
	if (($U_tract_count_to_RNAfold>3) && ($first_T != -1 )){
		
		 $new_sequence = substr($seq_DNA, 0, $concetive_T_site);
		 
		 $add_tail = $seq_length - $concetive_T_site;
		 $add_tail_dot =  "\." x $add_tail;
		 $warning = "U-tract is need to be removed";
		 
		# print "U-tract is need to be removed\n";
		# print "$later\t$last_fold_r\t$first_T\t$seq_length\t$concetive_T_site\t$position_TTT\t$rna_fold_output_seq_DNA[0]\n$rna_fold_output_seq_DNA[1]\t$new_sequence\n";
		 #system " echo $seq_DNA| RNAfold -p ";
		 #print "$U_tract_output";
	} else{
		
		 $new_sequence = $seq_DNA;
		 $add_tail =0;
		 $add_tail_dot = "";
		
	}
	
	
	
	
	
	
	
	my @rna_fold_output = `echo $new_sequence| RNAfold -p`;
	
	
	
	my $rna = $rna_fold_output[0]; chomp $rna;
	
	
	$rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
	my ($fold,$energy) = ($1,$2); 
	
	print "$rna\n$fold\n";
	
	
	
	my $new_folding_length = () = $fold =~ /\(/g;
	
	print "ew_folding_length\t$new_folding_length\n";
	
	
	if (($U_tract_count_to_RNAfold>3) && ($first_T != -1 )&&($new_folding_length <7)){
		
		 $new_sequence = substr($seq_DNA, 0, $concetive_T_site+(7-$new_folding_length));
		
	}
	
	
	 @rna_fold_output = `echo $new_sequence| RNAfold -p`;
	
	
	
	 $rna = $rna_fold_output[0]; chomp $rna;
	
	
	$rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
	 ($fold,$energy) = ($1,$2); 
	
	
	
	
	
	
	
	
	#$fold =~ m/(\((.+)\)),(\((.+)\))/;
	#$fold =~ m/(\((.+)\))(\((.+)\))/;
	my $warning_2 = " " ;
	$warning_2 ="two_fold" if $fold =~ m/(\(\.{0,8}+\)+\.{0,15}\()/;
	#$fold =~ m/(\(*+\)*+)(\.{0,10}\(*+\)*+.{0,20}$/; 
#$fold =~ m/(\(+\.{0,3}\(+\.{0,2}\(+)(\.+)(\)+\.{0,3}\)+\.{0,2}\)+).{0,20}$/; 
	print "two_fold\n" if $fold =~ m/(\(\.{0,8}+\)+\.{0,15}\()/;
	
	my $first_fold =index ($fold, "(");

	my $last_fold = rindex($fold, ")");
	
	
	print "last sequence \n$new_sequence\n$fold\n$first_fold\t$last_fold\n";
	
	my $last_RNA_fold_only_sequence= substr($new_sequence, $first_fold, ($last_fold-$first_fold+1));
	
	print "$last_RNA_fold_only_sequence\n";
	my ( $last_fold_5_end, $last_middle, $last_fold_3_end) ;
	
	if ($seq_DNA_o =~ m/$last_RNA_fold_only_sequence/) {
		
		( $last_fold_5_end, $last_middle, $last_fold_3_end)  = ($`, $&, $');
		print "left ($`) ($&) right($')\n";
		
	}
	
	my $flank_5end = substr($last_fold_5_end, -8);
	
	my $flank_3end = substr($last_fold_3_end, 0, 8);
	
	
	my $last_only_fold = substr($fold,  $first_fold, ($last_fold-$first_fold+1));
	
	
	my $last_TTS_sequence_fold =  convert_DNA_2_RNA($last_RNA_fold_only_sequence);
	
	my $flank_5end_RNA =  convert_DNA_2_RNA($flank_5end);
	my $flank_3end_RNA =  convert_DNA_2_RNA($flank_3end);
	my $seq_RNA_45 =  convert_DNA_2_RNA($seq_DNA_45);
	
	
	
	my $last_TTS_sequence = $flank_5end_RNA.$last_TTS_sequence_fold.$flank_3end_RNA;
	
	my $last_rna_fold = "\." x8 . "$last_only_fold"."\." x8 ;
	print "last \t $last_TTS_sequence\t$last_rna_fold\n";
	my @nucs_head = split(//,$flank_5end_RNA); 
	my $number_of_A_residues = 0;
	foreach my $nuc_head (@nucs_head) { 
		$number_of_A_residues++ if ($nuc_head eq 'A');
	}
	my @nucs = split(//,$flank_3end_RNA); 
	my $number_of_U_residues = 0;
	
	foreach my $nuc (@nucs) { 
		$number_of_U_residues++ if ($nuc eq 'U');
	}
	my $stem_loop =0;
	 $stem_loop = &is_stem_loop($last_rna_fold,$last_TTS_sequence) if ($last_rna_fold);
	
	 my $TTS_length_fold= length ($last_TTS_sequence_fold);
	 return "$seq_RNA_45\t$last_TTS_sequence\t$last_TTS_sequence_fold\t$TTS_length_fold\t$energy\t$flank_5end_RNA\t$number_of_A_residues\t$flank_3end_RNA\t$number_of_U_residues\t$stem_loop\t$warning\t$warning_2";
		
	 
		
}



sub count_U_stretch { 
	my ($seq) = @_;
	my $last_T = rindex($seq, "TT");
	my $tail = ($last_T < (length($seq)-10))? substr($seq,length($seq)-8,8) : substr($seq,$last_T-6,8); 
	my $first_A =index ($seq, "A");
	my $head = ($first_A >10)? substr($seq, 0, 8) : substr($seq, $first_A, 8); 
	my @nucs_head = split(//,$head); 
	my $number_of_A_residues = 0;
	foreach my $nuc_head (@nucs_head) { 
		$number_of_A_residues++ if ($nuc_head eq 'A');
	}
	my @nucs = split(//,$tail); 
	my $number_of_T_residues = 0;
	
	foreach my $nuc (@nucs) { 
		$number_of_T_residues++ if ($nuc eq 'T');
	}
	
	return "$head\t$number_of_A_residues\t$tail\t$number_of_T_residues";
}


sub is_stem_loop { 
	my ($fold, $DNA_stem_loop) = @_;
	
	$fold =~ m/(\(+\.{0,3}\(+\.{0,2}\(+)(\.+)(\)+\.{0,3}\)+\.{0,2}\)+).{0,20}$/; 

	if ($1 and $2 and $3 ) { 
		my $length_1 =length($1);
		my $length_2 =length($2);
		my $length_3 =length($3);
		my $stem_start_site = index($fold, $1);
		my $stem_sequence_left = substr ($DNA_stem_loop, $stem_start_site, $length_1);
		my $stem_sequence_left_RNA = convert_DNA_2_RNA($stem_sequence_left);
		my $stem_GC_content = 0;
		$stem_GC_content = &GC_content_calculation($stem_sequence_left) if ($stem_sequence_left) ;
		my $stem_diff = abs($length_1 - $length_3);
		my $loop_stem_diff_1 = $length_2 / $length_1;
		my $loop_stem_diff_2 = $length_2/ $length_3;
		if ($stem_diff <= 6 and $loop_stem_diff_1 < 2 and $loop_stem_diff_2 < 2) { 
			my $stem_loop = "$fold\t$length_1\t$length_2\t$length_3\t$stem_sequence_left_RNA\t$stem_GC_content";
			return $stem_loop;
		} else {
			my $stem_loop = "$fold\t$length_1\t$length_2\t$length_3\tno_stem_loop\t0";
			return $stem_loop;
		}
	}
	return "$fold\t0\t0\t0\t0\t0";
}

sub GC_content_calculation {
	my ($sequence_DNA) = @_;
	my $a=($sequence_DNA=~tr/A//);
	my $b=($sequence_DNA=~tr/C//);
	my $c=($sequence_DNA=~tr/G//);
	my $d=($sequence_DNA=~tr/T//);
	my $Total=$a+$b+$c+$d;
	my $GCper=sprintf("%.2f", (($b+$c)/($Total)*100));
	return $GCper;
	
}

sub reverse_complement_IUPAC {
        my $dna = shift;
        # reverse the DNA sequence
        my $revcomp = reverse($dna);
        # complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
	}
	

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
	
sub convert_DNA_2_RNA {
		$_ = shift (@_);
		if  (/^[ATCG]/){
			my $RNA_total =$_;
			$RNA_total =~ s/T/U/g;
			return $RNA_total;
		}
		
	}
