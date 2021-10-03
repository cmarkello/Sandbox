#!/usr/bin/env bash
# real_read_variant_call_evaluation.sh: Run happy and rtg vcfevals on input experiment VCF variant called file against various regions.

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

function make_bedfile() {
    if [ ! -e "${3}" ] ; then
        docker run \
        -e HIGH_CONF_BED=${1} \
        -e REGION_BED=${2} \
        -e HIGH_CONF_REGION_BED=${3} \
        -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1
        bedtools intersect \
        -a ${HIGH_CONF_BED} \
        -b ${REGION_BED} \
        > ${HIGH_CONF_REGION_BED}
    fi
}

function make_nosnp1kg_bedfile() {
    if [ ! -e "${3}" ] ; then
        docker run \
        -e HIGH_CONF_BED=${1} \
        -e SNP1KG_VCF_SITES_FILE=${2} \
        -e HIGH_CONF_NOSNP1KG_BED=${3} \
        -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1
        bedtools subtract \
        -a ${HIGH_CONF_BED} \
        -b ${SNP1KG_VCF_SITES_FILE} \
        > ${HIGH_CONF_NOSNP1KG_BED}
    fi
}

function run_happy() {
    if [ ! -d "${1}" ] ; then
        OUT_DIR=${1}
        TRUTH_VCF=${2}
        TRUTH_BED=${3}
        CALLED_VCF=${4}
        REF_FASTA=${5}
        docker run \
        -e OUT_DIR=${OUT_DIR} \
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
            -o ${OUT_DIR}/happy_output
    fi
}

function run_rtgvcfeval() {
    if [ ! -d "${1}" ] ; then
        OUT_DIR=${1}
        TRUTH_VCF=${2}
        TRUTH_BED=${3}
        CALLED_VCF=${4}
        REF_FASTA=${5}
        docker run \
        -e TRUTH_VCF=${TRUTH_VCF} \
        -e TRUTH_BED=${TRUTH_BED} \
        -e CALLED_VCF=${CALLED_VCF} \
        -e REF_FASTA=${REF_FASTA} \
        -v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
            rtg vcfeval \
            --baseline ${TRUTH_VCF} \
            --calls ${CALLED_VCF} \
            --evaluation-regions ${TRUTH_BED} \
            --output ${OUT_DIR} \
            --template ${REF_FASTA}.sdf \
            --threads 16
    fi
}

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
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/ALL.GRCh38.genotypes.20170504.no_segdups_gt10kb.sites.sorted.vcf.gz "${WORKDIR}/ALL.GRCh38.genotypes.20170504.no_segdups_gt10kb.sites.sorted.vcf.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alldifficultregions.bed.gz "${WORKDIR}/GRCh38_alldifficultregions.bed.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alllowmapandsegdupregions.bed.gz "${WORKDIR}/GRCh38_alllowmapandsegdupregions.bed.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz "${WORKDIR}/GRCh38_MHC.bed.gz"
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

# Evaluate Called File in All High Confident Regions
make_nosnp1kg_bedfile ${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ALL.GRCh38.genotypes.20170504.no_segdups_gt10kb.sites.sorted.vcf.gz ${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.specific.bed
run_happy happy_vcfeval_output_${CHILD_NAME}_allhighconfregions "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
run_happy happy_vcfeval_output_${CHILD_NAME}_allhighconfregions_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.specific.bed" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.specific.bed" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}

# Evaluate Called File in Difficult Regions
for REGION in "MHC.bed" "alllowmapandsegdupregions.bed" "alldifficultregions.bed" ; do
    cd $WORKDIR
    make_bedfile ${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed GRCh38_${REGION}.gz ${CHILD_NAME}_GRCh38_v4.2.1.${REGION}
    make_bedfile ${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.specific.bed GRCh38_${REGION}.gz ${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}

    run_happy happy_vcfeval_output_${CHILD_NAME}_${REGION} "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_happy happy_vcfeval_output_${CHILD_NAME}_${REGION}_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_${REGION} "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_${REGION}_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
done

if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_v4.2.1.CMRG.bed"
    REGION="CMRG.bed"
    cd $WORKDIR
    make_nosnp1kg_bedfile ${CHILD_NAME}_GRCh38_v4.2.1.${REGION} ALL.GRCh38.genotypes.20170504.no_segdups_gt10kb.sites.sorted.vcf.gz ${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}
    run_happy happy_vcfeval_output_${CHILD_NAME}_allhighconfregions "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_happy happy_vcfeval_output_${CHILD_NAME}_allhighconfregions_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
    run_rtgvcfeval rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions_NOSNP1KG "${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.CMRG.vcf.gz" "${CHILD_NAME}_GRCh38_v4.2.1.NO_SNP1KG.specific.${REGION}" "${CHILD_NAME}.vcf.gz" ${REF_FASTA}
fi

