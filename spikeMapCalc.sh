#! /bin/bash
#SBATCH --account=def-mlorincz					# required (format def-name)
#SBATCH --cpus-per-task=20							# number of cpus
#SBATCH --mem-per-cpu=4G								# memory; default unit is megabytes
#SBATCH --time=00-12:00									# time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL

# Check loss of mappability upon inclusion of Spike-In

# 1 = Ref .fa file
# 2 = Ref .bed file of all chr
# 3 = Spike-In .fa file

module load bedtools/2.30.0
module load seqtk/1.3
module load samtools/1.12
module load bwa/0.7.17


RUN_THREAD=20
SCRATCH="$HOME/scratch/$SLURM_JOB_ID/"
mkdir $SCRATCH

READ_LEN=66
SLIDING_WIN=5

REF=$1
SIZES=${REF//.fa/.sizes}
BED=$2
SPIKE=$3

REF_NAME=$(basename $REF)
REF_NAME=${REF_NAME//.fa/}
SPIKE_NAME=$(basename $SPIKE)
SPIKE_NAME=${SPIKE_NAME//.fa/}
echo "Testing "$REF_NAME" and "$SPIKE_NAME"..."

REF_BED="$SCRATCH/Ref.bed"
REF_FASTA="$SCRATCH/Ref.fa"
REF_FASTQ="$SCRATCH/Ref_R"$READ_LEN"-w"$SLIDING_WIN".fq"
REF_FASTQ_GZ="$REF_FASTQ.gz"
SAM="$SCRATCH/temp.sam"
REF_BAM=$REF_NAME"_R"$READ_LEN"-w"$SLIDING_WIN".bam"
REF_TXT=$REF_NAME"_R"$READ_LEN"-w"$SLIDING_WIN".txt"
CAT="$SCRATCH/${REF//.fa/}"+"${SPIKE//.fa/}.fa"
CAT_BAM=$REF_NAME"+"$SPIKE_NAME"_R"$READ_LEN"-w"$SLIDING_WIN".bam"
CAT_TXT=$REF_NAME"+"$SPIKE_NAME"_R"$READ_LEN"-w"$SLIDING_WIN".txt"

if [[ ! -f $REF_FASTQ_GZ ]]; then
	echo "Creating windowed fastq..."
	bedtools makewindows -w $READ_LEN -s $SLIDING_WIN -b $BED > $REF_BED
	bedtools getfasta -fi $REF -fo $REF_FASTA -bed $REF_BED
	rm $REF_BED
	seqtk seq -F '#' $REF_FASTA > $REF_FASTQ
	rm $REF_FASTA
	gzip $REF_FASTQ
fi

if [[ ! -f $REF_TXT ]]; then
#	bwa index $REF
	bwa mem -t $RUN_THREAD $REF $REF_FASTQ_GZ > $SAM
	samtools view -bhS -@ $RUN_THREAD $SAM > $REF_BAM

	samtools view $REF_BAM | cut -f 5 | awk 'BEGIN{for(x=0; x<256; x++){a[x]=0}}{
		a[$1]++
	 } END {for(x=0; x<256; x++){print x, a[x]}}' > $REF_TXT
fi

if [[ ! -f $CAT ]]; then # Create concatenated genome
	CHR_PREFIX=">$REF_NAME_"
	cat $REF | sed "s/>/$CHR_PREFIX/g" > temp
	CHR_PREFIX=">$SPIKE_NAME_"
	cat $SPIKE | sed "s/>/$CHR_PREFIX/g" > temp2
	cat temp temp2 > $CAT
	rm temp temp2
	bwa index $CAT
fi

bwa mem -t $RUN_THREAD $CAT $REF_FASTQ_GZ > $SAM
samtools view -bhS -@ $RUN_THREAD $SAM > $CAT_BAM
rm $SAM

samtools view $CAT_BAM | cut -f 5 | awk 'BEGIN {
	for(x=0; x<256; x++) {
		a[x]=0;
	}
}{
	a[$1]++;
} END {for(x=0; x<256; x++){print x, a[x]}}' > $CAT_TXT

for FILE in $REF_BAM $CAT_BAM
do
	samtools view $FILE | awk -v win=$SLIDING_WIN 'OFS="\t"{
		split($1, a, /[:,]/);
		print a[1], a[2], a[2]+win-1, $5/60;
	}' > ${FILE//.bam/.bedgraph}
	bedGraphToBigWig ${FILE//.bam/.bedgraph} $SIZES ${FILE//.bam/_mappability.bw}
done

rm -r $SCRATCH

