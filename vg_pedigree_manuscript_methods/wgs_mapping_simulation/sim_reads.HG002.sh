#!/bin/bash
# sim_reads.HG002.sh: Simulate reads from the HG002 sample graph and stratify by various regions.

set -ex
set -o pipefail

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


WORKDIR=${HOME}/run_sim_reads


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"

# Format the reads
zcat HG002.R1.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG002.R1-1m.fq.gz
zcat HG002.R2.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG002.R2-1m.fq.gz
paste <(zcat HG002.R1-1m.fq.gz) <(zcat HG002.R2-1m.fq.gz) | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 > "HG002.R1-shuffled-1m.fq"; print $2,$4,$6,$8 > "HG002.R2-shuffled-1m.fq"}'
gzip HG002.R1-shuffled-1m.fq
gzip HG002.R2-shuffled-1m.fq
docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bbmap:38.93--he522d1c_0 \
reformat.sh t=2 in=HG002.R1-shuffled-1m.fq.gz in2=HG002.R2-shuffled-1m.fq.gz out=HG002_merged_interleaved-shuffled-1m.fastq.gz

## SIMULATE BASELINE READS

# Simulate 1M reads for all high confident region stratification
docker run \
-e NREADS=1000000 \
-e FASTQ=HG002_merged_interleaved-shuffled-1m.fastq.gz \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
/bin/bash -c 'vg sim -r -I -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.1m.raw.gam'

# Simulate 100M reads for difficult region stratification
docker run \
-e NREADS=100000000 \
-e FASTQ=HG002_merged_interleaved-shuffled-1m.fastq.gz \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
/bin/bash -c 'vg sim -r -I -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.100m.raw.gam'

declare -a REGION_LIST=( "high_conf_hg002_v4.2.1_regions" "all_difficult_regions_hg002_v4.2.1_regions" "alllowmapandsegdupregions_hg002_v4.2.1_regions" "mhc_hg002_v4.2.1_regions" "cmrg_hg002_v4.2.1_regions" "high_conf_NO1000GP_hg002_v4.2.1_regions" )
declare -a BED_FILE_LIST=( "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "HG002_GRCh38_v4.2.1.all_difficult_regions.bed" "HG002_GRCh38_v4.2.1.alllowmapandsegdupregions.bed" "HG002_GRCh38_v4.2.1.MHC.bed" "HG002_GRCh38_CMRG_smallvar_v1.00.bed" "HG002_GRCh38_CHROM1-22_v4.2.1.highconf.NO_SNP1KG.bed" )

for index in "${!REGION_LIST[@]}"; do
    REGION="${REGION_LIST[index]}"
    BED_FILE="${BED_FILE_LIST[index]}"
    mkdir -p ${REGION}
    cp sim.raw.gam hg002_sample_grch38.vg ${REGION}/
    cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG002_sim_reads/${REGION}
    singularity exec --env REGION=${REGION},BED_FILE=${BED_FILE} -H ${PWD}:${HOME} --pwd ${HOME} \
    -B /data/markellocj/benchmark_data/HG002_cohort:${HOME}/HG002_cohort \
    docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
    /bin/bash -xc 'set -eo pipefail && \
    echo annotate with paths && \
    time vg annotate \
    -p -x hg002_sample_grch38.vg -a sim.raw.gam \
    --bed-name HG002_cohort/${BED_FILE} \
    --threads 32 > sim.${REGION}.gam && \
    time vg filter -i -U -F "" sim.${REGION}.gam > sim.gam && \
    rm -f sim.${REGION}.gam sim.fq.gz sim.fq && \
    echo convert to fastq && \
    time vg view -X -a sim.gam | gzip > sim.fq.gz && \
    gunzip sim.fq.gz && \
    sed "s/_1\$//g" sim.fq | sed "s/_2\$//g" > sim.paired.fq && \
    echo format true position information && \
    time vg view -a sim.gam | jq -c -r "[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(\",\")] else [\".\"] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv" > true.pos'
done




