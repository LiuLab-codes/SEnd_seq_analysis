set -ue

# To convert the merged fasta sequence data to 3' adatpor labled single end data, this program will remove the 5'end adator at last

# get name inforamtion

fa_file=$1
name=$2


grep -E "AGTTGATAGCAGGTT" -B 1 $1 | grep -v "^--$" > $name.f.fa
grep -vE "AGTTGATAGCAGGTT" $1| grep -E  "AACCTGCTATCAACT" -B 1 | grep -v "^--$" > $name.r.fa

fastx_reverse_complement -i $name.r.fa -o $name.r_RC.fa

cat $name.f.fa $name.r_RC.fa > $name.c2.fa 
rm $name.f.fa $name.r_RC.fa $name.r.fa

# seperate the adaptor containing reads into two parts, one is for the 5' end and one is for the 3'end

sed "s/AGTTGATAGCAGGTT/XXXXX\\`echo -e '\n\r'`YYYYY/g"  $name.c2.fa  > $name.sep.c2.fa

grep -vE "XXXXX"  $name.sep.c2.fa | grep -E "YYYYY" -B1 | sed "s/YYYYY//g" |grep -v '^--$' > $name.sep.c2.P2.fa

dos2unix $name.sep.c2.P2.fa

# use the cutadapt software to remove possible 5' adaptor and remove too short reads


~/.local/bin/cutadapt  -e 0.1  -u 4 -m 18  --trim-n -o $name.sep.c2.P2_cut.fa   $name.sep.c2.P2.fa

~/.local/bin/cutadapt  -l 150 -o $name.sep.c2.P2_cut_75bp.fa $name.sep.c2.P2_cut.fa
fastx_reverse_complement -i $name.sep.c2.P2_cut_75bp.fa -o $name.sep.c2.P2_cut_75bp_RC.fa

# save the final files to the fold of final.

file="3L_single_end_final"
if [ -d "$file" ]
then
	echo "$file exists and the new files will be added into final."
else
	mkdir 3L_single_end_final
fi

 cp  *75bp_RC.fa 3L_single_end_final
 rm $name.c2.fa  $name.sep.c2.fa $name.sep.c2.P2.fa $name.sep.c2.P2_cut.fa  $name.sep.c2.P2_cut_75bp.fa   $name.sep.c2.P2_cut_75bp_RC.fa

~/bin/bowtie2 -p8 -f --very-sensitive-local  --ff -x ~/ref_E.coli/NC_000913_spikeRNA  ./3L_single_end_final/$name.sep.c2.P2_cut_75bp_RC.fa -S ./3L_single_end_final/$name.3L_single_end.sam

~/bin/bin/samtools view -Sb ./3L_single_end_final/$name.3L_single_end.sam > ./3L_single_end_final/$name.3L_single_end.bam

~/bin/bin/samtools sort  ./3L_single_end_final/$name.3L_single_end.bam > ./3L_single_end_final/$name.3L_single_end.sorted.bam
~/bin/bin/samtools index ./3L_single_end_final/$name.3L_single_end.sorted.bam