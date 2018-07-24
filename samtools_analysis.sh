set -ue

# To analyzed the length distribution of ribosome RNA and non-ribosome RNA.

# get input and name inforamtion

sam_file=$1
#name=$2


# get the file name before.sam


filename_1=$(basename "$sam_file")

filename_1="${filename_1%.*}"


#calculate the insert length distribution of total reads 

cat $sam_file  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_total_mapping_1.txt


# convert the sam file to the index bam file

~/bin/bin/samtools view -Sb  $sam_file >$filename_1.bam
#~/bin/bin/samtools sort $filename_1.bam> $filename_1.sorted.bam
#~/bin/bin/samtools index $filename_1.sorted.bam

#~/bin/bin/samtools view $filename_1.sorted.bam  "gi|556503834|ref|NC_000913.3|:4643500-4645606" > $filename_1.C_D.mapping.sam
#~/bin/bin/samtools view $filename_1.sorted.bam  "gi|556503834|ref|NC_000913.3|:4642788-4643494" > $filename_1.H.mapping.sam
#~/bin/bin/samtools view $filename_1.sorted.bam "gi|556503834|ref|NC_000913.3|:4645646-4647500" > $filename_1.Fluc.mapping.sam
#cat $filename_1.H.mapping.sam  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_H_spikinRNA.txt
#cat $filename_1.Fluc.mapping.sam  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_Fluc_spikinRNA.txt
#cat $filename_1.C_D.mapping.sam |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt

#line_C_D_mapping_sam=`wc -l < $filename_1.C_D.mapping.sam`
#echo $line_C_D_mapping_sam
#if [ "$line_C_D_mapping_sam" -ge "10" ]
#then 
	
#  cat $filename_1.C_D.mapping.sam |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt
#fi
#rm $filename_1.C_D.mapping.sam $filename_1.H.mapping.sam $filename_1.Fluc.mapping.sam

# to make the pair end to single end fasta file then do further analysis

~/bin/bin/samtools sort -n $filename_1.bam  -o $filename_1.sorted_n.bam

~/bin/bin/samtools view -bf 0x2 $filename_1.sorted_n.bam | bedtools bamtobed -i stdin -bedpe > $filename_1\_PE.bedpe
cut -f 1,2,6,7,8,9 $filename_1\_PE.bedpe > $filename_1\_PE_2.bedpe
rm $filename_1\_PE.bedpe
bedtools getfasta -fi  ~/ref_E.coli/NC_000913_spikeRNA.fna -bed $filename_1\_PE_2.bedpe  -s  -name  > $filename_1\_PE_2_S.fa

#
~/bin/bowtie2 -p 5   -f --very-sensitive-local    --ff -x ~/ref_E.coli/NC_000913_spikeRNA $filename_1\_PE_2_S.fa -S $filename_1\_PE_2_S.bowtie2mapping.sam 
cat $filename_1\_PE_2_S.bowtie2mapping.sam  |cut -f 6 |sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_total_mapping.txt
~/bin/bin/samtools view -Sb $filename_1\_PE_2_S.bowtie2mapping.sam >$filename_1\_PE_2_S.bowtie2mapping.bam
~/bin/bin/samtools sort $filename_1\_PE_2_S.bowtie2mapping.bam > $filename_1\_PE_2_S.bowtie2mapping.sorted.bam
~/bin/bin/samtools index $filename_1\_PE_2_S.bowtie2mapping.sorted.bam 


bedtools intersect -abam  $filename_1\_PE_2_S.bowtie2mapping.sorted.bam -b  ~/ref_E.coli/RNA/NC_000913.rRNA_tRNA.bed -v > $filename_1.rRNA_filtered.bam
~/bin/bin/samtools index $filename_1.rRNA_filtered.bam
~/bin/bin/samtools view  $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered.sam
cat $filename_1.rRNA_filtered.sam   |cut -f 6 | sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1.RNA_filtered_length.txt


~/bin/bin/samtools view $filename_1\_PE_2_S.bowtie2mapping.sorted.bam "gi|556503834|ref|NC_000913.3|:4643500-4645606" > $filename_1.C_D.mapping.sam
~/bin/bin/samtools view $filename_1\_PE_2_S.bowtie2mapping.sorted.bam "gi|556503834|ref|NC_000913.3|:4642788-4643494" > $filename_1.H.mapping.sam
~/bin/bin/samtools view $filename_1\_PE_2_S.bowtie2mapping.sorted.bam "gi|556503834|ref|NC_000913.3|:4645646-4647500" > $filename_1.Fluc.mapping.sam
cat $filename_1.H.mapping.sam  |cut -f 6 | sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_H_spikinRNA.txt
cat $filename_1.Fluc.mapping.sam  |cut -f 6 |sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_Fluc_spikinRNA.txt
cat $filename_1.C_D.mapping.sam |cut -f 6 |sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt

line_C_D_mapping_sam=`wc -l < $filename_1.C_D.mapping.sam`
echo $line_C_D_mapping_sam
if [ "$line_C_D_mapping_sam" -ge "10" ]
then

  cat $filename_1.C_D.mapping.sam |cut -f 6 |sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt
fi
rm $filename_1.C_D.mapping.sam $filename_1.H.mapping.sam $filename_1.Fluc.mapping.sam













rm $filename_1.rRNA_filtered.sam 

~/bin/bin/samtools sort $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered_sorted.bam
~/bin/bin/samtools index $filename_1.rRNA_filtered_sorted.bam
# calculate the rRNA length distribution

~/bin/bin/samtools view -L ~/ref_E.coli/RNA/NC_000913.rRNA_tRNA.bed  $filename_1\_PE_2_S.bowtie2mapping.sorted.bam |cut -f 6 |sed "s/M//g" | grep -vE "^-" | grep -vE "^0" > $filename_1.rRNA_length.txt

file="txt"
if [ -d "$file" ]
then
	echo "$file exists and the new files will be added into txt file."
else
	mkdir txt
fi

mv $filename_1*.txt txt/



#remove the adaptor duplicated reads.
