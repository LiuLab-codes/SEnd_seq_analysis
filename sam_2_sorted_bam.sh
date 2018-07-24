set -ue

# To convert the merged fasta sequence data to 3' adatpor labled single end data, this program will remove the 5'end adator at last

# get name inforamtion

sam_file=$1
#name=$2

filename_1=$(basename "$sam_file")

filename_1="${filename_1%.*}"



~/bin/bin/samtools view -Sb $sam_file > $filename_1.bam

~/bin/bin/samtools sort  $filename_1.bam > $filename_1.sorted.bam
~/bin/bin/samtools index $filename_1.sorted.bam