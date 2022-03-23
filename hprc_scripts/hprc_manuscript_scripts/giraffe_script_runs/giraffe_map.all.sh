#!/usr/bin/env bash
# giraffe_map.sh: map HG001, NA12891, NA12892, HG002, HG003, HG004, HG005, HG006, HG007 150bp paired reads with VG GIRAFFE to the HS38d1-based graph references
# giraffe_map.sh /data/Udpbinfo/usr/markellocj/hprc_feb3_giraffe_runs HG002 GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7 .novaseq.pcr-free.30x
module load singularity

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

WORKDIR=$1
SAMPLE_NAME=$2
INDEX_BASENAME=$3
REF_FASTA_BASENAME=$4
READ_BASENAME=$5
REF_FASTA="${REF_FASTA_BASENAME}.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp_${SAMPLE_NAME}"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR

SEQ_DICT="${REF_FASTA_BASENAME}.dict"
XG_IDX="${INDEX_BASENAME}.xg"
GBWT_IDX="${INDEX_BASENAME}.gbwt"
MIN_IDX="${INDEX_BASENAME}.min"
GG_IDX="${INDEX_BASENAME}.gg"
DIST_IDX="${INDEX_BASENAME}.dist"
GBZ_IDX="${INDEX_BASENAME}.giraffe.gbz"
INPUT_BAM=giraffe_${SAMPLE_NAME}.bam

# run giraffe mapper
cd $WORKDIR
echo ${SAMPLE_NAME}${READ_BASENAME}
singularity exec \
--env SEQ_DICT=${SEQ_DICT},GBZ_IDX=${GBZ_IDX},XG_IDX=${XG_IDX},MIN_IDX=${MIN_IDX},DIST_IDX=${DIST_IDX},SAMPLE_NAME="${SAMPLE_NAME}",READ_BASENAME="${READ_BASENAME}" \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/vgteam/vg:v1.38.0 \
vg giraffe \
-t 16 \
-x ${XG_IDX} \
-m ${MIN_IDX} \
-Z ${GBZ_IDX} \
-d ${DIST_IDX} \
-f "${SAMPLE_NAME}${READ_BASENAME}.R1.fastq.gz" \
-f "${SAMPLE_NAME}${READ_BASENAME}.R2.fastq.gz" \
--output-format BAM \
--ref-paths ${SEQ_DICT} > ${INPUT_BAM}

cd $WORKDIR
INPUT_BAM=giraffe_${SAMPLE_NAME}.bam
singularity exec \
--env INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
    samtools sort -@ 32 ${INPUT_BAM} -O BAM > positionsorted.${INPUT_BAM} && rm ${INPUT_BAM}

singularity exec \
--env INPUT_BAM=${INPUT_BAM},SEQ_DICT=${SEQ_DICT} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://broadinstitute/picard:2.21.9 \
  java -Xmx20g -XX:ParallelGCThreads=16 -jar /usr/picard/picard.jar \
  ReorderSam \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=positionsorted.${INPUT_BAM} \
  OUTPUT=reordered.positionsorted.${INPUT_BAM} \
  SEQUENCE_DICTIONARY=${SEQ_DICT} && rm positionsorted.${INPUT_BAM}

# Indel realign input BAM
cd $WORKDIR
singularity exec \
--env INPUT_BAM=${INPUT_BAM},SAMPLE_NAME=${SAMPLE_NAME} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
samtools addreplacerg \
-@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 \
reordered.positionsorted.${INPUT_BAM} > gatk_ready.reordered.positionsorted.${INPUT_BAM}

singularity exec \
--env INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
samtools index -@ 32 gatk_ready.reordered.positionsorted.${INPUT_BAM}

singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME},REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://broadinstitute/gatk3:3.8-1 \
  java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
  -R ${REF_FASTA} \
  -I gatk_ready.reordered.positionsorted.${INPUT_BAM} -o ${SAMPLE_NAME}.intervals



singularity exec \
--env REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/freebayes-samtools:1.2.0_1.10 \
bamleftalign < gatk_ready.reordered.positionsorted.${INPUT_BAM} \
> left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM} \
--fasta-reference ${REF_FASTA} \
--compressed

singularity exec \
--env REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/freebayes-samtools:1.2.0_1.10 \
samtools index -b left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM} left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM}.bai

singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME},REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0 \
java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
--remove_program_records \
-drf DuplicateRead \
--disable_bam_indexing \
-nt 16 \
-R reference.fa \
-I left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM} \
--out forIndelRealigner.intervals.${SAMPLE_NAME}

awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals.${SAMPLE_NAME} > forIndelRealigner.intervals.${SAMPLE_NAME}.bed

singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME},REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0 \
bedtools slop -i forIndelRealigner.intervals.${SAMPLE_NAME}.bed -g "${REF_FASTA}.fai" -b "160" > ${SAMPLE_NAME}.intervals.widened.bed

mv ${SAMPLE_NAME}.intervals.widened.bed forIndelRealigner.intervals.${SAMPLE_NAME}.bed

singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME},INPUT_BAM=${INPUT_BAM},REF_FASTA=${REF_FASTA} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba \
java -Xmx20G -jar /opt/abra2/abra2.jar \
  --targets forIndelRealigner.intervals.${SAMPLE_NAME}.bed \
  --in left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM} \
  --out indel_realigned.left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM} \
  --ref ${REF_FASTA} \
  --index \
  --threads 16

singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
    samtools index -@ 32 indel_realigned.left_shifted.gatk_ready.reordered.positionsorted.${INPUT_BAM}

