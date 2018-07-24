set -ue

# To analyzed the length distribution of ribosome RNA and non-ribosome RNA.

# get input and name inforamtion

sam_file=$1
#name=$2


# get the file name before.sam



filename_1=$(basename "$sam_file")

filename_1="${filename_1%.*}"

#calculate the insert length distribution of total reads 

cat $sam_file  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_total_mapping.txt


# convert the sam file to the index bam file

~/bin/bin/samtools view -Sb  $sam_file >$filename_1.bam
~/bin/bin/samtools sort $filename_1.bam> $filename_1.sorted.bam
~/bin/bin/samtools index $filename_1.sorted.bam

~/bin/bin/samtools view $filename_1.sorted.bam  "gi|556503834|ref|NC_000913.3|:4643500-4645606" > $filename_1.C_D.mapping.sam
~/bin/bin/samtools view $filename_1.sorted.bam  "gi|556503834|ref|NC_000913.3|:4642788-4643494" > $filename_1.H.mapping.sam
~/bin/bin/samtools view $filename_1.sorted.bam "gi|556503834|ref|NC_000913.3|:4645646-4647500" > $filename_1.Fluc.mapping.sam
cat $filename_1.H.mapping.sam  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_H_spikinRNA.txt
cat $filename_1.Fluc.mapping.sam  |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_Fluc_spikinRNA.txt
#cat $filename_1.C_D.mapping.sam |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt

line_C_D_mapping_sam=`wc -l < $filename_1.C_D.mapping.sam`
echo $line_C_D_mapping_sam
if [ "$line_C_D_mapping_sam" -ge "10" ]
then 
        
  cat $filename_1.C_D.mapping.sam |cut -f 9 | grep -vE "^-" | grep -vE "^0" >$filename_1\_C_D_spikinRNA.txt
fi


bedtools intersect -abam $filename_1.sorted.bam -b  ~/ref_E.coli/RNA/NC_000913.rRNA_tRNA.bed -v > $filename_1.rRNA_filtered.bam

~/bin/bin/samtools view  $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered.sam
cat $filename_1.rRNA_filtered.sam   |cut -f 9 |  grep -vE "^-" | grep -vE "^0" >$filename_1.RNA_filtered_length.txt

rm $filename_1.rRNA_filtered.sam 
~/bin/bin/samtools sort $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered_sorted.bam
~/bin/bin/samtools index $filename_1.rRNA_filtered_sorted.bam

~/bin/bin/samtools view -L ~/ref_E.coli/RNA/NC_000913.rRNA_tRNA.bed  $filename_1.sorted.bam |cut -f 9 | grep -vE "^-" | grep -vE "^0" > $filename_1.rRNA_length.txt

file="txt"
if [ -d "$file" ]
then
        echo "$file exists and the new files will be added into txt file."
else
        mkdir txt
fi

mv $filename_1*.txt txt/

rm $filename_1.C_D.mapping.sam $filename_1.H.mapping.sam $filename_1.Fluc.mapping.sam
