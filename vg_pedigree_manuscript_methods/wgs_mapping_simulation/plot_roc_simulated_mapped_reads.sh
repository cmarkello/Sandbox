#!/bin/bash
module load R
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG002_sim_reads

SCRIPT_DIR="/data/markellocj/git/giraffe-sv-paper/scripts/plotting"

# Replace all names of mappers with human-readable ones
function humanize_names() {
    sed -e 's/[a-zA-Z0-9_.]*bwamem[a-zA-Z0-9_.]*/BWA/' -e 's/[a-zA-Z0-9_.]*giraffe_parental_default[a-zA-Z0-9_.]*/GiraffeParental/' -e 's/[a-zA-Z0-9_.]*giraffe_snp1kg_default[a-zA-Z0-9_.]*/Giraffe1000GP/' -e 's/[a-zA-Z0-9_.]*giraffe_primary_default[a-zA-Z0-9_.]*/GiraffePrimary/'
}

#REGION="high_conf_hg002_v4.2.1_regions"
REGION="all_difficult_regions_hg002_v4.2.1_regions"
GRAPH="hg002_sample_grch38"
READS="HG002"
PAIRING="paired"
SPECIES="human"
PE_OPTS="-- -pe"
GBWT='sampled\.64'
LINEAR_GRAPH="NOTAPPLICABLE"
echo "Extracting toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv for region ${REGION}"
cat bwamem_roc_stats.${REGION}.tsv | head -n1 > toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv

# Grab linear BWA for this graph and these reads
tail -q -n +2 bwamem_roc_stats.${REGION}.tsv | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | grep ${PE_OPTS} | humanize_names >> toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv


# Grab giraffe and map and graphaligner non-linear, and other linear mappers
tail -q -n +2 giraffe_parental_roc_stats.${REGION}.tsv giraffe_snp1kg_roc_stats.${REGION}.tsv giraffe_primary_roc_stats.${REGION}.tsv | grep ${PE_OPTS} | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | sed 's/null/0/g' | humanize_names >> toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv

time Rscript ${SCRIPT_DIR}/plot-roc.R toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv roc-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.${REGION}.png

