#!/bin/env perl -w 


# This script is used to split the merged fasta file to two paired-end files based on the 3' adatpors. If different adapotrs were used, the sequences in this script are need to be updated.

# this script is based on the 3' adapotr of "/5Phos/rNrNrN rNrArA rCrCrU rGrCrU rArUrC rArArC rUrG/3ddC/ ".

use strict; 

use File::HomeDir;
my $fileSpec = File::HomeDir->my_home;


my $direcotry = `pwd`;
chomp $direcotry;
#print "$direcotry";
opendir(DIR, $direcotry) or die $!;

my $file_name_1;
while (my $file = readdir(DIR)) {
  next if  $file =~ /\.$/;
 if ($file =~ /\.n\.fa/){
 	$file=~ /\.n\.fa/;
	 $file_name_1 = $`;
	print "$file_name_1\n";
	#&pair_end_fq_check("$file_name_1\_R1_001.fastq","$file_name_1\_R2_001.fastq");
	&fasta_split;
 };
  
 

}
sub fasta_split{

my $fa_name = " ";
my $fa_sequence = " ";


system("cat $file_name_1.n.fa | $fileSpec/miniconda2/bin/seqkit locate -p CAGTTGATAGCAGG -G > $file_name_1.txt ");

my %position_infor;

open( position_infor, "./$file_name_1.txt") or die "can't read postion information file $!\n";
	while (<position_infor>){
		chomp;
		my $line = $_;
		my @line_infor = split(/\t/,$line); 
		$position_infor{($line_infor[0])}[1] = "2" if $position_infor{($line_infor[0])};
		$position_infor{($line_infor[0])} = ["1", "$line_infor[3]", "$line_infor[4]", "$line_infor[5]" ]unless $position_infor{($line_infor[0])};
		
		
	}

close position_infor;

open(fasta_p1, '>', "./$file_name_1\_split_R1.fa") or die "can't write output file $!\n";
open(fasta_p2, '>', "./$file_name_1\_split_R2.fa") or die "can't write output file $!\n";

open (merged_fa," ./$file_name_1.n.fa") or die "can't open gene_annotation file $!\n";
my $adptor_exist_check=0;
my $str=0;
my ($strand_3a, $start_3a, $end_3a)=("0", "0", "0");
while (<merged_fa>) {
	chomp;
	my $line = $_;
if ($line=~ /^\>/){
	$adptor_exist_check =0;
	$fa_name = $line;   
	$str =$fa_name;
	$str=~  tr/\>//d;
	print "$str\n" if (($str%10000)==0);
	next unless $position_infor{$str};
	$adptor_exist_check =1 if ($position_infor{$str}[0]==1);
	($strand_3a, $start_3a, $end_3a) = ($position_infor{$str}[1], $position_infor{$str}[2], $position_infor{$str}[3]);
	#print "($strand_3a, $start_3a, $end_3a)\n";
	
}else {
	next if $adptor_exist_check ==0;
	$fa_sequence = $line;	
	#my $seqence = "$fa_name\n$fa_sequence";
	#print "$seqence\n";
	#my $detect_string = `echo -e '$seqence' | seqkit locate -p CAGTTGATAGCAGG -G | wc -l `;
	#my $detect_string = `echo '$seqence' | $fileSpec/miniconda2/bin/seqkit locate -p CAGTTGATAGCAGG -G | wc -l `;
	#print "number is $detect_string";
	
	#my $detect_infor = `echo '$seqence' | $fileSpec/miniconda2/bin/seqkit locate -p CAGTTGATAGCAGG -G | grep -vE "seqID"  `;
	
	
	#chomp;
	#print "$detect_infor";
	#my $nr_of_lines = $detect_infor =~ tr/\n//;
	#print "$nr_of_lines\n" if ( $nr_of_lines >1);
	#next if ($nr_of_lines !=1);
	#my @adapotr_position = split(/\t/,$detect_infor); 
	next if (($start_3a<15) or ($end_3a >(length($fa_sequence)-15)));
	#print "\t$fa_name\t$adapotr_position[3]\t$adapotr_position[4]\t$adapotr_position[5]\n";
	my ($P1, $P2) = (" ", " ");
	
	if ($strand_3a =~ /\+/){
		
		$P1 = substr($fa_sequence, 0,($start_3a-4));
		#print "$P1\n";
		next unless $P1;
		#next if (length($P1) < 15);
		$P1 = reverse $P1;
		$P1 =~ tr/ACGTacgt/TGCAtgca/;
		
		if ($P1 =~ /CTCTCCACGT/){
			#print "$P1\n";
			my $posiiton_5A = index ($P1, "CTCTCCACGT");
			
			$P1 = (length($P1)>($posiiton_5A+16))? (substr($P1, ($posiiton_5A+14))) : " "; 
			#print "$posiiton_5A\n$P1\n";
		}
		next unless $P1;
		next if (length($P1) < 15);
		
		
		$P1 = substr($P1, 0, 75) if (length($P1) >75); 
		
		
		my $adaptor_check = substr($fa_sequence, $end_3a,2);
		#print "fist two is $adaptor_check\n";
		next unless $adaptor_check eq "TT";
		
		$P2 = substr($fa_sequence, ($end_3a+6),75);
		
		next unless $P2;
		next if (length($P2) < 15);
		$P2 = reverse $P2;
		$P2 =~ tr/ACGTacgt/TGCAtgca/;
		print fasta_p1 "$fa_name\n$P1\n";
		print fasta_p2 "$fa_name\n$P2\n";
		#print fasta_p2 "$P2\n";
	} elsif ($strand_3a =~ /-/){
		
		$P1 = substr($fa_sequence, ($start_3a+3));
		next unless $P1;
		
		#print fasta_p1 "$P1\n";
		if ($P1 =~ /CTCTCCACGT/){
			#print "$P1\n";
			my $posiiton_5A = index ($P1, "CTCTCCACGT");
			
			$P1 = (length($P1)>($posiiton_5A+16))? (substr($P1, ($posiiton_5A+14))) : " "; 
			#print "$posiiton_5A\n$P1\n";
		}
		next unless $P1;
		next if (length($P1) < 15);
		$P1 = substr($P1, 0, 75) if (length($P1) >75); 
		my $adaptor_check = substr($fa_sequence, ($start_3a-3),2);
		
		next unless $adaptor_check eq "AA";
		#print "$adaptor_check\n";
		
		$P2 = substr($fa_sequence, 0, ($start_3a-7));
		
		#print "$fa_sequence\n$end_3a\t$adaptor_check\n$P2\n";
		$P2 =substr ($P2, (length($P2)-74), 75) if (length ($P2) >75);
		next unless $P2;
		next if (length($P2) < 15);
		#print "$fa_sequence\n$end_3a\t$adaptor_check\n$P2\n";
		#print fasta_p2 "$P2\n";
		print fasta_p1 "$fa_name\n$P1\n";
		print fasta_p2 "$fa_name\n$P2\n";
		
	}
	
	
}


	#print "$fa_name\n" if $_ =~ /\>/;
	
	
   # print "$fa_name\n$fa_sequence\n";



   
}

close merged_fa;

close fasta_p1;
close fasta_p2;


system ("rm $file_name_1.txt ");
#system ("cutadapt ");



}


sub pair_end_fq_check {
my ($paired_fq_file_1, $paired_fq_file_2) = @_; #Open bed file and sam file
print "$paired_fq_file_1\t$paired_fq_file_2\n";

my %fasta_name;





open (fastq_file_1, "./$paired_fq_file_1") or die "can't open the input fastq_file\n";
my $previous_line_number =-3;
my $line_number =0;
my @lines_file_1;
my $target_line =4;
while (<fastq_file_1>){
	$line_number +=1;
	$target_line +=1;
chomp;

my $line_1 = "$_";
if (($line_number == ($previous_line_number+4)) && ($line_1 =~ /^\Q@\E/)){
	my @line_name= split(/ /);
	#print "$line_name[0]\n";
	$fasta_name{$line_number}=$line_name[0];
		$target_line =1;
		
	$previous_line_number +=4;
}




}

close fastq_file_1;


open (fastq_file_2, "./$paired_fq_file_2") or die "can't open the input fastq_file\n";
 $previous_line_number =-3;
 $line_number =0;
my @lines_file_2;
$target_line =4;
while (<fastq_file_2>){
	$line_number +=1;
	$target_line +=1;
chomp;

my $line_2 = "$_";
if (($line_number == ($previous_line_number+4)) && ($line_2 =~ /^\Q@\E/)){
	my @line_name= split(/ /);
	#print "$line_name[0]\n";
	#print "$fasta_name{$line_number}\t$line_name[0]\n";
	print "$paired_fq_file_2 is wrong\n" unless ($fasta_name{$line_number} eq $line_name[0]);
	last unless $fasta_name{$line_number} eq $line_name[0];	#print "$line_name[0]\n";
	}
	$previous_line_number +=4;
}






close fastq_file_2;
}
