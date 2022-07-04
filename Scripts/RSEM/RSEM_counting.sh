

: " SCRIPT TO COUNT FEATURE FROM A BAM ALIGNMENT FILE.
 It requires an indexed genome (required by RSEM) using 'rsem-prepare-reference'
 This script does the following feature couting using rsem-calculate-expression. "

# Saving of the arguments into variables:
# Strandedness of the RNA-Seq experiment
STRND=$1
echo "$STRND"

# Name of the indexed genome.
INDX_GENOME_NAME=$2
INDX_GENOME_NAME=./${INDX_GENOME_NAME}
echo "Index genome: $INDX_GENOME_NAME"


# Input BAM file.
BAM=$3
echo "Input BAM: $BAM"
# Sample name
SAMPLE=${BAM%%.*}
SAMPLE=${SAMPLE/BAMs/COUNTS}  # Changing the directory to save the output counts. THIS TO CHANGE
echo "Sample: $SAMPLE"

# Number of threads to be used. Parallelism <3
N_THREADS=$4
echo "Nº of threads: $N_THREADS"

########################################################################################################################
: " TESTING OF THE INPUT BAM USING 'rsem-sam validator'"

# In memory of @Miki:
echo "Testing if $BAM is comptaible with RSEM"

rsem-sam-validator \
$BAM
########################################################################################################################
: " CALCULATION OF THE FEATURE COUNTS using 'rsem-calculate-expression' "

echo "Calculation of the expression counts of $SAMPLE..."

rsem-calculate-expression \
--strandedness $STRND \
--append-names \
--no-bam-output \
--num-threads $N_THREADS \
--alignments \
--paired-end \
$BAM \
$INDX_GENOME_NAME \
$SAMPLE
