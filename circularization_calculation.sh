set -ue
# To calculate the efficiency of circularization

# get name inforamtion

fa_file=$1


# get the file name before.sam


filename_1=$(basename "$fa_file")

filename_1="${filename_1%.*}"



grep -vE "(AGTTGATAGCAGG).*(AGTTGATAGCAGG)|(CCTGCTATCAACT).*(CCTGCTATCAACT)" $fa_file | grep -E "AGTTGATAGCAGGTT" -B 1 | grep -v "^--$" > $filename_1.f.fa
grep -vE "(AGTTGATAGCAGG).*(AGTTGATAGCAGG)|(CCTGCTATCAACT).*(CCTGCTATCAACT)" $fa_file | grep -E "AACCTGCTATCAACT" -B 1 | grep -v "^--$" > $filename_1.r.fa

fastx_reverse_complement -i $filename_1.r.fa -o $filename_1.r_RC.fa
cat $filename_1.f.fa $filename_1.r_RC.fa > $filename_1.c2.fa 
rm $filename_1.f.fa $filename_1.r_RC.fa $filename_1.r.fa


#seperate the forward adaptor and reverse adaptor, and then reverse complementary change the reverse adaptor congtaing reads. and then combine the two files togehter.

awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $filename_1.c2.fa  | grep -vE "^>" > $filename_1.length.txt



 ~/bin/fastaRegexFinder.py -f $filename_1.c2.fa  --noreverse  -r AGTTGATAGCAGGTT  -q |awk '{print $1, $2, $3}' > $filename_1.3L_position.txt 
 paste $filename_1.3L_position.txt $filename_1.length.txt | awk '{print $1,$2,$3,$4-$3+1,$4}' > $filename_1.3L_position_2.txt 
 
line_number="$(wc -l $filename_1.3L_position.txt   | awk '{print $1}')" 
 
 circulized_number="$(cat $filename_1.3L_position.txt |awk '{print $2}'| awk '{if($1==$1+0 && $1>4)print $1}' |wc -l )" 
 
 circulization_efficiency=`echo "scale=4; $circulized_number/$line_number*100" | bc`
 echo  "$fa_file \t $line_number \t $circulized_number \t $circulization_efficiency \t file_name:total_3L_nmeber:circularized_number:efficiency\n" >> ./$filename_1.circ_efficiency.txt
 #$line_number $circulized_number $circulization_efficiency 
 rm $filename_1.c2.fa $filename_1.3L_position.txt
#printf "%s\t%s\t%s\t%s\n" "$fa_file" "$line_number" " $circulized_number" "$circulization_efficiency" >> ./$filename_1.circ_efficiency.txt



#for x in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16  
#do


#fastx_collapser -v -i S$x.n.fa -o S$x\_collapse.n.fa
##cat legnth_S$x\_3L.bowtie2mapping.rRNA_filtered_sorted.bam.txt | cut -f 17,15 |grep -vE "\t0$" >legnth_S$x\_3L.bowtie2mapping.rRNA_filtered_sorted.bam_cut.txt
#awk '{ total += $2; count++ } END { print total/count }' legnth_S$x\_3L.bowtie2mapping.rRNA_filtered_sorted.bam_cut.txt
#done




#bedtools bamtobed -i S$x\_3L.bowtie2mapping.rRNA_filtered_sorted.bam > S$x\_3L_rRNA_filtered_SE.bedpe
#bedtools bamtobed -i S$x._5L_3L.bowtie2mapping.rRNA_filtered_sorted.bam > S$x\_5L_3L_rRNA_filtered_SE.bedpe
	#perl gene_length_5endplus100nt11142017.pl NC_000913.bed S$x._5L_3L.bowtie2mapping.rRNA_filtered_sorted.bam
	#perl gene_length_3endplus100nt11142017.pl NC_000913.bed S$x._5L_3L.bowtie2mapping.rRNA_filtered_sorted.bam
#perl gene_length_3endplus100nt11142017.pl NC_000913.bed S$x._5L_3L.bowtie2mapping_PE_2_S.bowtie2mapping.sorted.bam



#~/bin/bin/samtools sort S$x._5L_3L.bowtie2mapping.bam > S$x._5L_3L.bowtie2mapping_sorted.bam

#~/bin/bin/samtools index S$x._5L_3L.bowtie2mapping_sorted.bam
 #sort -n -k2 S$x._5L_3L.bowtie2mapping_PE_2.bedpe > S$x._5L_3L.bowtie2mapping_PE_2_sorted.bedpe
 #samtools view -F 0x904 -c   S$x._5L_3L.bowtie2mapping.rRNA_filtered_sorted.bam


