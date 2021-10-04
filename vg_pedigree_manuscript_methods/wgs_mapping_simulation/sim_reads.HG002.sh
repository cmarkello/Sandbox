#!/bin/bash
module load singularity
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG002_sim_reads

## SIMULATE BASELINE READS
singularity shell -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG002:${HOME}/HG002 \
docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
NREADS=1000000 \
FASTQ=HG002/HG002_merged_interleaved-shuffled-1m.fastq.gz
echo simulate reads
vg sim -r -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.raw.gam

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




