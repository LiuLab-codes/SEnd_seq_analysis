set -ue

####################################### bam2bed.sh #################################
####	convert a sorted bam format file to bed format file
####    apply the alignment files for further quantitative analysis
####   require to install the bed tools http://bedtools.readthedocs.io/en/latest/
####	run example: sh bam2bed.sh  *.bam
############################################################################################################









bam_file=$1



filename_1=$(basename "$bam_file")

filename_1="${filename_1%.*}"

bedtools bamtobed -i $bam_file> $filename_1.bedpe