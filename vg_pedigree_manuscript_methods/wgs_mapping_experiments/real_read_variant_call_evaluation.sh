#!/usr/bin/env bash
# real_read_deepvariant_call.sh: indel realign input BAMs and run the DeepVariant Genotyper

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

USE_DEFAULT_MODEL=true
WORKDIR=${HOME}/run_deeptrio_genotyping
CHILD_NAME="${1}"
VCF_FILE_CHILD="${2}"
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna"


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
copy ${VCF_FILE_CHILD} "${WORKDIR}/${CHILD_NAME}.vcf.gz"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna "${WORKDIR}/${REF_FASTA}"
if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
fi

docker run \
-e REF_FASTA=${REF_FASTA} \
-v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
    rtg format \
    ${REF_FASTA} \
    -o ${REF_FASTA}.sdf

# Evaluate Called File
for REGION in "HG002_GRCh38_CMRG_smallvar_v1.00.bed" "HG002_GRCh38_v4.2.1.MHC.bed" "HG002_GRCh38_v4.2.1.alllowmapandsegdupregions.bed" "HG002_GRCh38_v4.2.1.all_difficult_regions.bed" "HG002_GRCh38_CHROM1-22_v4.2.1.highconf.NO_SNP1KG.specific.bed" "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" ; do


cd $WORKDIR
rm -r $WORKDIR/happy_vcfeval_output_${CHILD_NAME}_${REGION}
TRUTH_VCF="${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
TRUTH_BED="${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
CALLED_VCF="${CHILD_NAME}.vcf.gz" \
docker run \
-e TRUTH_VCF=${TRUTH_VCF} \
-e TRUTH_BED=${TRUTH_BED} \
-e CALLED_VCF=${CALLED_VCF} \
-e REF_FASTA=${REF_FASTA} \
-v ${PWD}:${HOME} -w ${HOME} jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    ${TRUTH_VCF} \
    ${CALLED_VCF} \
    -f ${TRUTH_BED} \
    --reference ${REF_FASTA} \
    --threads 16 \
    --engine=vcfeval \
    -o happy_vcfeval_output_${CHILD_NAME}_${REGION}

rm -r $WORKDIR/rtg_vcfeval_output_${CHILD_NAME}_${REGION}
docker run \
-e TRUTH_VCF=${TRUTH_VCF} \
-e TRUTH_BED=${TRUTH_BED} \
-e CALLED_VCF=${CALLED_VCF} \
-e REF_FASTA=${REF_FASTA} \
-e OUT_NAME=rtg_vcfeval_output_${CHILD_NAME}_${REGION} \
-v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
    rtg vcfeval \
    --baseline ${TRUTH_VCF} \
    --calls ${CALLED_VCF} \
    --evaluation-regions ${TRUTH_BED} \
    --output ${OUT_NAME} \
    --template ${REF_FASTA}.sdf \
    --threads 16

if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    REGION="CMRG"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.bed"
    cd $WORKDIR
    rm -r $WORKDIR/happy_vcfeval_output_${CHILD_NAME}_${REGION}
    TRUTH_VCF="${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    TRUTH_BED="${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    CALLED_VCF="${CHILD_NAME}.vcf.gz" \
    docker run \
    -e TRUTH_VCF=${TRUTH_VCF} \
    -e TRUTH_BED=${TRUTH_BED} \
    -e CALLED_VCF=${CALLED_VCF} \
    -e REF_FASTA=${REF_FASTA} \
    -v ${PWD}:${HOME} -w ${HOME} jmcdani20/hap.py:v0.3.12 \
        /opt/hap.py/bin/hap.py \
        ${TRUTH_VCF} \
        ${CALLED_VCF} \
        -f ${TRUTH_BED} \
        --reference ${REF_FASTA} \
        --threads 16 \
        --engine=vcfeval \
        -o happy_vcfeval_output_${CHILD_NAME}_${REGION}

    rm -r $WORKDIR/rtg_vcfeval_output_${CHILD_NAME}_${REGION}
    docker run \
    -e TRUTH_VCF=${TRUTH_VCF} \
    -e TRUTH_BED=${TRUTH_BED} \
    -e CALLED_VCF=${CALLED_VCF} \
    -e REF_FASTA=${REF_FASTA} \
    -e OUT_NAME=rtg_vcfeval_output_${CHILD_NAME}_${REGION} \
    -v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
        rtg vcfeval \
        --baseline ${TRUTH_VCF} \
        --calls ${CALLED_VCF} \
        --evaluation-regions ${TRUTH_BED} \
        --output ${OUT_NAME} \
        --template ${REF_FASTA}.sdf \
        --threads 16
    
    # TODO
    cd /data/markellocj/benchmark_data/HG002_cohort
    bedtools intersect \
    -a /data/markellocj/benchmark_data/HG002_cohort/HG002_GRCh38_CMRG_smallvar_v1.00.bed \
    -b /data/markellocj/benchmark_data/HG002_cohort/HG002_GRCh38_CHROM1-22_v4.2.1.highconf.NO_SNP1KG.specific.bed \
    > HG002_GRCh38_v4.2.1.NO_SNP1KG.specific.CMRG.bed
fi
