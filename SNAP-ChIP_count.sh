#! /bin/bash
#SBATCH --account=def-mlorincz       		   # required
#SBATCH --cpus-per-task=8                        # number of MPI processes
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=00-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL


# SNAP-ChIP Normalization
# Run in folder with bam files to be normalized
# Normalization = RPM * (#Barcode Reads WT / #Mm Reads WT) / (#Barcode Reads Sample / #Mm Reads Sample)
# Will iterate through all marks covered in Barcode fasta file (excluding wt)
# Setup for Graham ( but can be adapted for other servers)
# Will produce un-normalized bigwigs, normalization of experimental samples should be done using output from this script

module load samtools/1.12
module load bedtools/2.30.0
module load bwa/0.7.17

THREADS=8
MEM=4G
SORT_MEM=3G
BIN_SIZE=1
SMOOTH_WIN=0
NORMALIZE="CPM"
MAPQ=5

BARCODE="/project/def-mlorincz/reference_genomes/snapChIP/KmetStat_barcoded.fasta"
BAM_COVERAGE="/project/def-mlorincz/scripts/utilities/miniconda3/bin/bamCoverage"
BWA="bwa mem -v 0 -t $THREADS $BARCODE"

OUTPUT="Barcode_Reads.txt"

SECONDS=$(date +%s)
TEMP_DIR="$HOME/scratch/$SECONDS"
TEMP_BAM="$TEMP_DIR/temp.bam"
TEMP_SAM="$TEMP_DIR/temp.sam"
TEMP_FQ="$TEMP_DIR/temp.fastq"
TEMP_FQ1="$TEMP_DIR/temp_1.fastq"
TEMP_FQ2="$TEMP_DIR/temp_2.fastq"
mkdir $TEMP_DIR

for FILE in *.bam
do
	BAM_BARCODE=${FILE//.bam/_barcode.bam}

	# Figure out if paired end or single end
	PAIRED=$(($(samtools view $FILE | head -n 1 | cut -f 2) % 2)) # 1 if paired, 0 if not

	if [[ $PAIRED == 1 ]]; then

		echo "Sorting $FILE by read name..."
		samtools sort -n -o $TEMP_BAM -@ $THREADS -m $SORT_MEM $FILE
		
		echo "Extracting Fastqs..."
		bamToFastq -i $TEMP_BAM -fq $TEMP_FQ1 -fq2 $TEMP_FQ2
		rm $TEMP_BAM
		
		echo "Aligning to Barcode..."
		$BWA $TEMP_FQ1 $TEMP_FQ2 > $TEMP_SAM
		rm $TEMP_FQ1 $TEMP_FQ2
		
	else

		echo "Extracting Fastqs from $FILE..."
		bamToFastq -i $FILE -fq $TEMP_FQ
		
		echo "Aligning to Barcode..."
		$BWA $TEMP_FQ > $TEMP_SAM
		rm $TEMP_FQ
			
	fi
		
	echo "Filtering Barcode Reads..."
	samtools view -b -@ $THREADS -q 5 -o $TEMP_BAM $TEMP_SAM
	rm $TEMP_SAM

	echo "Sorting bam..."
	samtools sort -o $BAM_BARCODE -@ $THREADS -m $MEM $TEMP_BAM
	rm $TEMP_BAM

	echo "Indexing Barcode Reads..."
	samtools index $BAM_BARCODE
	samtools idxstats $BAM_BARCODE > ${BAM_BARCODE//.bam/_stats.txt}
	MM_READS=$(samtools view -c -F 4 $FILE)
	echo -e "Mm_Reads\ta\t$MM_READS" >> ${BAM_BARCODE//.bam/_stats.txt}
	MM_READS=$(samtools view -c -q $MAPQ $FILE)
	echo -e "Mm_Reads_MapQ\ta\t$MM_READS" >> ${BAM_BARCODE//.bam/_stats.txt}
	MM_READS=$(samtools view -c -q $MAPQ -F 1024 $FILE)
	echo -e "Mm_Reads_MapQ_NoDup\ta\t$MM_READS" >> ${BAM_BARCODE//.bam/_stats.txt}
	MM_READS=$(samtools view -c -f 4 $FILE)
	echo -e "Unmapped\ta\t$MM_READS" >> ${BAM_BARCODE//.bam/_stats.txt}

	samtools index $FILE
	$BAM_COVERAGE -b $FILE --outFileName ${FILE//.bam/.bw} --binSize $BIN_SIZE -p $THREADS --normalizeUsing $NORMALIZE --smoothLength $SMOOTH_WIN --outFileFormat bigwig --minMappingQuality $MAPQ --ignoreDuplicates
done

HEADER="Barcode\t"
for FILE in *.txt
do
	HEADER=$HEADER${FILE//_stats.txt/}"\t"
done

let N=1
	
for FILE in *.txt
do
	if [[ $N == 1 ]] ; then
		echo "Header +" $FILE
		echo -e $HEADER > $OUTPUT
		awk 'BEGIN{OFS="\t"} {print $1, $3}' $FILE >> $OUTPUT
	else
		echo "File" $FILE
		awk 'BEGIN{getline; OFS="\t"; print $0}{
			if ( NR == FNR ) {
				a[NR]=$0;
			} else {
				print a[FNR+1], $3
			}
		}' $OUTPUT $FILE > temp.txt
		mv temp.txt $OUTPUT
	fi
	((N++))
done

rm -r $TEMP_DIR
