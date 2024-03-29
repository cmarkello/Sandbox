version 1.0

workflow vus_cooccurrence {
    input {
        File SAMPLE_BCF
        File SAMPLE_BCF_INDEX
        File? SITES_VCF
        File? SITES_VCF_INDEX
        String OUTPUT_NAME
        String GENE
        Boolean QC_FILTER = false
    }
    
    if (QC_FILTER) {
        call hail_preprocess_qc {
            input:
                in_sample_bcf=SAMPLE_BCF,
                in_sample_bcf_index=SAMPLE_BCF_INDEX
        }
    }
    call setup_brcaexchange_data {
        input:
            in_gene=GENE,
            outname=OUTPUT_NAME
    }
    call normalize_brcaexchange_data as normalize_brca_path {
        input:
            in_vcf=setup_brcaexchange_data.pathogenic_vcf,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    call normalize_brcaexchange_data as normalize_brca_vus {
        input:
            in_vcf=setup_brcaexchange_data.vus_vcf,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    call normalize_brcaexchange_data as normalize_brca_all {
        input:
            in_vcf=setup_brcaexchange_data.all_vcf,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    
    File preprocessed_sample_bcf = select_first([hail_preprocess_qc.hail_qc_bcf, SAMPLE_BCF])
    File preprocessed_sample_bcf_index = select_first([hail_preprocess_qc.hail_qc_bcf_index, SAMPLE_BCF_INDEX])
    
    call setup_sample_data {
        input:
            in_sample_bcf=preprocessed_sample_bcf,
            in_sample_bcf_index=preprocessed_sample_bcf_index,
            in_gene=GENE,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    
    call intersect_variants as intersect_path_variants {
        input:
            in_base_vcf=normalize_brca_path.normalized_vcf,
            in_base_vcf_index=normalize_brca_path.normalized_vcf_index,
            in_query_vcf=setup_sample_data.normalized_vcf,
            in_query_vcf_index=setup_sample_data.normalized_vcf_index,
            complement=false
    }
    call intersect_variants as intersect_vus_variants {
        input:
            in_base_vcf=normalize_brca_vus.normalized_vcf,
            in_base_vcf_index=normalize_brca_vus.normalized_vcf_index,
            in_query_vcf=setup_sample_data.normalized_vcf,
            in_query_vcf_index=setup_sample_data.normalized_vcf_index,
            complement=false
    }
    call intersect_variants as intersect_all_variants {
        input:
            in_base_vcf=normalize_brca_all.normalized_vcf,
            in_base_vcf_index=normalize_brca_all.normalized_vcf_index,
            in_query_vcf=setup_sample_data.normalized_vcf,
            in_query_vcf_index=setup_sample_data.normalized_vcf_index,
            complement=true
    }
    
    call concat_vcfs {
        input:
            in_a_vcf=intersect_vus_variants.intersected_vcf,
            in_a_vcf_index=intersect_vus_variants.intersected_vcf_index,
            in_b_vcf=intersect_all_variants.intersected_vcf,
            in_b_vcf_index=intersect_all_variants.intersected_vcf_index,
    }
    
    call detect_vus_benign {
        input:
            in_intersect_vus_vcf=concat_vcfs.concatenated_vcf,
            in_intersect_vus_vcf_index=concat_vcfs.concatenated_vcf_index,
            in_intersect_path_vcf=intersect_path_variants.intersected_vcf,
            in_intersect_path_vcf_index=intersect_path_variants.intersected_vcf_index,
            in_sites_vcf=SITES_VCF,
            in_sites_vcf_index=SITES_VCF_INDEX,
            outname=OUTPUT_NAME
    }
    
    output {
        File cooccurrence_report = detect_vus_benign.cooccurrence_report
        File complete_cooccurrence_report = detect_vus_benign.complete_cooccurrence_report
        File hwe_hom_vus_report = detect_vus_benign.hwe_report
        File hwe_stat_hist = detect_vus_benign.hwe_stat_hist
        File hwe_pvalue_hist = detect_vus_benign.hwe_pvalue_hist
        File hwe_freq_hist = detect_vus_benign.hwe_freq_hist
        File apparent_benign_vus_vcf = detect_vus_benign.apparent_benign_vus_vcf
        File apparent_benign_vus_vcf_index = detect_vus_benign.apparent_benign_vus_vcf_index
    }
}

task setup_brcaexchange_data {
    input {
        String in_gene
        String outname
    }
    command <<<
        set -exu -o pipefail
        
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
        gzip -d hg38.p12.fa.gz
        gzip -d ncbiRefSeq.txt.gz
        
        python3 /usr/src/app/extract_brcaexchange_data.py -g ~{in_gene} -o ~{outname}
    >>>
    output {
        File reference_file = "hg38.p12.fa"
        File ref_seq_file = "ncbiRefSeq.txt"
        File pathogenic_vcf = glob("*.pathogenic.vcf")[0]
        File vus_vcf = glob("*.vus.vcf")[0]
        File all_vcf = glob("*.all.vcf")[0]
    }
    runtime {
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task normalize_brcaexchange_data {
    input {
        File in_vcf
        File in_ref_file
        File in_seq_file
    }
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_vcf} input.vcf
        bgzip input.vcf
        tabix -p vcf input.vcf.gz
       
        mkdir -p ${PWD}/seqrepo
        seqrepo -r ${PWD}/seqrepo pull -i 2019-06-20
        export HGVS_SEQREPO_DIR=${PWD}/seqrepo/2019-06-20
        python3 /usr/src/app/hgvs_normalize.py \
            -i input.vcf.gz \
            -o brcaexchange.norm.vcf \
            -r ~{in_ref_file} \
            -g ~{in_seq_file}
        
        vcf-sort -p 8 brcaexchange.norm.vcf > brcaexchange.norm.sorted.vcf
        bgzip brcaexchange.norm.sorted.vcf
        tabix -p vcf brcaexchange.norm.sorted.vcf.gz
        rm -f brcaexchange.norm.vcf
    >>>
    output {
        File normalized_vcf = "brcaexchange.norm.sorted.vcf.gz"
        File normalized_vcf_index = "brcaexchange.norm.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task hail_preprocess_qc {
    input {
        File in_sample_bcf
        File in_sample_bcf_index
    }
    command <<<
    set -exu -o pipefail
    bcftools view -O v ~{in_sample_bcf} | bgzip > input.vcf.gz
    HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
    export PYSPARK_SUBMIT_ARGS="
    --jars $HAIL_HOME/backend/hail-all-spark.jar
    --conf spark.driver.extraClassPath="$HAIL_HOME/backend/hail-all-spark.jar"
    --conf spark.executor.extraClassPath=./hail-all-spark.jar
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator
    --conf spark.driver.memory=5G
    --conf spark.executor.memory=30G
    pyspark-shell"
    python3 <<CODE
    import hail as hl
    hl.init(default_reference = "GRCh38", log = 'Hail_QC.log')
    recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    mt = hl.import_vcf("input.vcf.gz", force_bgz = True, contig_recoding=recode, min_partitions=150)
    mt = mt.filter_rows(hl.len(mt.filters) == 0, keep = True)
    
    ## sample and variant QC ##
    # Run sample QC
    mt = hl.sample_qc(mt)
    # Run variant QC
    mt = hl.variant_qc(mt)
    # Summarize each field
    mt.summarize()
    
    ## genotype QC ##
    #ab = mt.AD[1] / hl.sum(mt.AD)
    #filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
    #                        (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
    #                        (mt.GT.is_hom_var() & (ab >= 0.9)))
    #
    #fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition_ab))
    #print(f'Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.')
    #mt = mt.filter_entries(filter_condition_ab)
    
    ## variant QC and filtering ##
    mt = mt.filter_rows(hl.min(mt.variant_qc.AC) > 0)
    #mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
    mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-2)
    # Bin variants by frequency
    mt = mt.annotate_rows(frq_bin =(hl.case()
        .when(hl.min(mt.variant_qc.AC) == 1, "singleton")
        .when((hl.min(mt.variant_qc.AC) >= 2) & (hl.min(mt.variant_qc.AC) <= 4), "mac2_4")
        .when((hl.min(mt.variant_qc.AC) >= 4) & (hl.min(mt.variant_qc.AC) <= 19), "mac5_19")
        .when((hl.min(mt.variant_qc.AC) >= 20) & (hl.min(mt.variant_qc.AF) < 0.001), "mac20_maf01")
        .when((hl.min(mt.variant_qc.AF) >= 0.001) & (hl.min(mt.variant_qc.AF) < 0.01), "maf01_maf1")
        .when((hl.min(mt.variant_qc.AF) >= 0.01) & (hl.min(mt.variant_qc.AF) < 0.05), "maf1_5")
        .when(hl.min(mt.variant_qc.AF) >= 0.05, "maf5")
        .default("None")))
    # Count the number of variants in each bin
    mt.aggregate_rows(hl.agg.counter(mt.frq_bin))
    hl.export_vcf(mt, "raw_sample.hail_qc.vcf")
    CODE
    bcftools view -Ob raw_sample.hail_qc.vcf > raw_sample.hail_qc.bcf
    bcftools index raw_sample.hail_qc.bcf
    >>>
    output {
        File hail_qc_bcf = "raw_sample.hail_qc.bcf"
        File hail_qc_bcf_index = "raw_sample.hail_qc.bcf.csi"
    }
    runtime {
        cpu: 1
        docker: 'quay.io/cmarkello/hail:latest'
    }
}

task setup_sample_data {
    input {
        File in_sample_bcf
        File in_sample_bcf_index
        String in_gene
        File in_ref_file
        File in_seq_file
    }
    command <<<
        set -exu -o pipefail
        
        # Extract specific region from vcf
        ln -s ~{in_sample_bcf} input_raw_sample.bcf
        ln -s ~{in_sample_bcf_index} input_raw_sample.bcf.csi
        
        mkdir -p ${PWD}/cache
        export PYENSEMBL_CACHE_DIR=${PWD}/cache
        pyensembl install --release 99 --species homo_sapiens
        gene_region=$(
        python <<CODE
        import pyensembl
        ensembl = pyensembl.EnsemblRelease(release=99)
        gene_data = ensembl.genes_by_name("~{in_gene}")
        print('chr{}:{}-{}'.format(gene_data[0].contig, gene_data[0].start, gene_data[0].end))
        CODE
        )
        
        bcftools filter --threads 8 --regions ${gene_region} input_raw_sample.bcf > raw_sample.vcf
        bgzip raw_sample.vcf
        tabix -p vcf raw_sample.vcf.gz
        
        # Separate vcf variant data from genotype data
        bcftools view -G raw_sample.vcf.gz -O v -o raw_sample.no_genotypes.vcf
        bcftools query -f 'GT\t[%GT\t]\n' raw_sample.vcf.gz > raw_sample.genotypes.vcf
        bgzip raw_sample.no_genotypes.vcf
        tabix -p vcf raw_sample.no_genotypes.vcf.gz
        
        # Normalize only the variant data
        mkdir -p ${PWD}/seqrepo
        seqrepo -r ${PWD}/seqrepo pull -i 2019-06-20
        export HGVS_SEQREPO_DIR=${PWD}/seqrepo/2019-06-20
        python3 /usr/src/app/hgvs_normalize.py \
            -i raw_sample.no_genotypes.vcf.gz \
            -o raw_sample.no_genotypes.norm.vcf \
            -r ~{in_ref_file} \
            -g ~{in_seq_file}
        
        # Merge normalized variant data with genotype data
        bcftools view -H raw_sample.no_genotypes.norm.vcf > raw_sample.no_genotypes.norm.no_header.vcf
        paste raw_sample.no_genotypes.norm.no_header.vcf raw_sample.genotypes.vcf > raw_sample.norm.paste.vcf
        bcftools view -h raw_sample.vcf.gz >> raw_sample.norm.vcf
        cat raw_sample.norm.paste.vcf >> raw_sample.norm.vcf
        vcf-sort -p 8 raw_sample.norm.vcf > raw_sample.norm.sorted.vcf
        bgzip raw_sample.norm.sorted.vcf
        tabix -p vcf raw_sample.norm.sorted.vcf.gz
        rm -f raw_sample.norm.vcf raw_sample.norm.paste.vcf raw_sample.no_genotypes.norm.no_header.vcf raw_sample.no_genotypes.norm.vcf raw_sample.no_genotypes.vcf.gz raw_sample.no_genotypes.vcf.gz.tbi raw_sample.genotypes.vcf
    >>>
    output {
        File normalized_vcf = "raw_sample.norm.sorted.vcf.gz"
        File normalized_vcf_index = "raw_sample.norm.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task intersect_variants {
    input {
        File in_base_vcf
        File in_base_vcf_index
        File in_query_vcf
        File in_query_vcf_index
        Boolean complement
    }
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_base_vcf} base.vcf.gz
        ln -s ~{in_base_vcf_index} base.vcf.gz.tbi
        ln -s ~{in_query_vcf} query.vcf.gz
        ln -s ~{in_query_vcf_index} query.vcf.gz.tbi
        
        if [ ~{complement} == false ]; then
            bcftools isec \
                -O v \
                -n =2 -w 1 \
                -o intersected_variants.vcf \
                query.vcf.gz \
                base.vcf.gz
        else
            bcftools isec \
                -O v \
                -w 1 \
                -C \
                -o intersected_variants.vcf \
                query.vcf.gz \
                base.vcf.gz
        fi
        vcf-sort -p 8 intersected_variants.vcf > intersected_variants.sorted.vcf
        bgzip intersected_variants.sorted.vcf
        tabix -p vcf intersected_variants.sorted.vcf.gz
        rm -f intersected_variants.vcf
    >>>
    output {
        File intersected_vcf = "intersected_variants.sorted.vcf.gz"
        File intersected_vcf_index = "intersected_variants.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task concat_vcfs {
    input {
        File in_a_vcf
        File in_a_vcf_index
        File in_b_vcf
        File in_b_vcf_index
    }
    command <<<
        ln -s ~{in_a_vcf} a.vcf.gz
        ln -s ~{in_a_vcf_index} a.vcf.gz.tbi 
        ln -s ~{in_b_vcf} b.vcf.gz
        ln -s ~{in_b_vcf_index} b.vcf.gz.tbi
        
        bcftools concat \
            -O v \
            --threads 8 \
            -o "concatenated.vcf" \
            a.vcf.gz \
            b.vcf.gz
        
        vcf-sort -p 8 concatenated.vcf > concatenated.sorted.vcf
        bgzip concatenated.sorted.vcf
        tabix -p vcf concatenated.sorted.vcf.gz
        rm -f concatenated.vcf
    >>>
    output {
        File concatenated_vcf = "concatenated.sorted.vcf.gz"
        File concatenated_vcf_index =  "concatenated.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task detect_vus_benign {
    input {
        File in_intersect_vus_vcf
        File in_intersect_vus_vcf_index
        File in_intersect_path_vcf
        File in_intersect_path_vcf_index
        File? in_sites_vcf
        File? in_sites_vcf_index
        String outname
    }
    
    Boolean in_sites_vcf_available = defined(in_sites_vcf)
    Boolean in_sites_vcf_index_available = defined(in_sites_vcf_index)
    
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_intersect_vus_vcf} vus.vcf.gz
        ln -s ~{in_intersect_vus_vcf_index} vus.vcf.gz.tbi
        ln -s ~{in_intersect_path_vcf} path.vcf.gz
        ln -s ~{in_intersect_path_vcf_index} path.vcf.gz.tbi
        if [[ ~{in_sites_vcf_available} == true && ~{in_sites_vcf_index_available} == true ]]; then
            ln -s ~{in_sites_vcf} sites.vcf.gz
            ln -s ~{in_sites_vcf_index} sites.vcf.gz.tbi
            python3 /usr/src/app/detect_vus_benign.py \
                -i vus.vcf.gz \
                -j path.vcf.gz \
                -s sites.vcf.gz \
                -o ~{outname}.cooccurrence_report \
                -v ~{outname}.apparent_benign_vus_list.vcf
        else
            python3 /usr/src/app/detect_vus_benign.py \
                -i vus.vcf.gz \
                -j path.vcf.gz \
                -o ~{outname}.cooccurrence_report \
                -v ~{outname}.apperent_benign_vus_list.vcf
        fi
    
        vcf-sort -p 8 ~{outname}.apparent_benign_vus_list.vcf > ~{outname}.apparent_benign_vus_list.sorted.vcf
        bgzip ~{outname}.apparent_benign_vus_list.sorted.vcf
        tabix -p vcf ~{outname}.apparent_benign_vus_list.sorted.vcf.gz
        rm -f ~{outname}.apparent_benign_vus_list.vcf
    >>>
    output {
        File cooccurrence_report = "~{outname}.cooccurrence_report.txt"
        File complete_cooccurrence_report = "complete_~{outname}.cooccurrence_report.txt"
        File hwe_report = "hom_vus_hwe_~{outname}.cooccurrence_report.txt"
        File hwe_stat_hist = "hom_vus_chi_square_stat.~{outname}.cooccurrence_report.png"
        File hwe_pvalue_hist = "hom_vus_chi_square_pvalue.~{outname}.cooccurrence_report.png"
        File hwe_freq_hist = "hom_vus_allele_frequencies.~{outname}.cooccurrence_report.png"
        File apparent_benign_vus_vcf = "~{outname}.apparent_benign_vus_list.sorted.vcf.gz"
        File apparent_benign_vus_vcf_index = "~{outname}.apparent_benign_vus_list.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}





