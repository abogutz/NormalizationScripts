#! /bin/bash
#SBATCH --account=def-mlorincz            # required (format def-name)
#SBATCH --cpus-per-task=20                        # number of cpus
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=02-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL
# Check loss of mappability upon Spike-In

# 1 = Host .fa file
# 2 = Host .bed file of all chr
# 3 = Spike-In .fa file

module load bedtools/2.30.0
module load seqtk/1.3
module load samtools/1.12
module load bwa/0.7.17


RUN_THREAD=20
SCRATCH="$HOME/scratch/"

REF=$1
BED=$2
SPIKE=$3
REF_BED="$SCRATCH/Ref.bed"
REF_FASTA="$SCRATCH/Ref.fa"
REF_FASTQ="$SCRATCH/Ref.fq"
SAM="$SCRATCH/temp.sam"
REF_BAM=${REF//.fa/}".bam"
CAT="$SCRATCH/Cat.fa"
CAT_BAM=${REF//.fa/}"+"${SPIKE//.fa/}".bam"

#bedops --chop 100 $2 > $REF_BED
bedtools makewindows -w 100 -s 5 -b $BED > $REF_BED
bedtools getfasta -fi $REF -fo $REF_FASTA -bed $REF_BED
rm $REF_BED
seqtk seq -F '#' $REF_FASTA > $REF_FASTQ
rm $REF_FASTA
gzip $REF_FASTQ
REF_FASTQ=$REF_FASTQ".gz"

bwa index $REF
bwa mem -t $RUN_THREAD $REF $REF_FASTQ > $SAM
samtools view -bhS -@ $RUN_THREAD $SAM > $REF_BAM

cat $REF $SPIKE > $CAT
bwa index $CAT
bwa mem -t $RUN_THREAD $CAT $REF_FASTQ > $SAM
samtools view -bhS -@ $RUN_THREAD $SAM > $CAT_BAM
rm $SAM

samtools view $REF_BAM | cut -f 5 | awk 'BEGIN{for(x=0; x<256; x++){a[x]=0}{
	a[$1]++
} END {for(x=0; x<256; x++){print x, a[x]}' > "Ref.txt"

samtools view $CAT_BAM | cut -f 5 | awk 'BEGIN{for(x=0; x<256; x++){a[x]=0}{
	a[$1]++
} END {for(x=0; x<256; x++){print x, a[x]}' > "Cat.txt"


