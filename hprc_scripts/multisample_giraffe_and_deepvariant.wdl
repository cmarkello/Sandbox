version 1.0

### multisample_giraffe_and_deepvariant.wdl ###
## Author: Charles Markello
## Description: Multi-sample workflow running the core VG Giraffe mapping and DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

import "./giraffe_and_deepvariant.wdl" as vgGiraffeDeepVariantWorkflow

workflow vgTrioPipeline {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        description: "Multi-sample workflow running the core VG Giraffe mapping and DeepVariant calling workflow for single sample datasets."
    }
    input {
        Array[File]+ SAMPLE_INPUT_READ_FILE_1_LIST      # List of input 1st read pair fastq.gz
        Array[File]+ SAMPLE_INPUT_READ_FILE_2_LIST      # List of input 2nd read pair fastq.gz where the order follows the same sample order of SAMPLE_INPUT_READ_FILE_1_LIST
        Array[String]+ SAMPLE_NAME_LIST                 # List of sample names where the order follows the same sample order of SAMPLE_INPUT_READ_FILE_1_LIST
        # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        String VG_CONTAINER = "quay.io/vgteam/vg:ci-3272-f6ec6f200ef9467ec1b62104a852acb187d4d3d8"
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        String? GIRAFFE_OPTIONS                         # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        File XG_FILE                                    # Path to .xg index file
        File GBWT_FILE                                  # Path to .gbwt index file
        File GGBWT_FILE                                 # Path to .gg index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
        File? TRUTH_VCF                                 # Path to .vcf.gz to compare against
        File? TRUTH_VCF_INDEX                           # Path to Tabix index for TRUTH_VCF
        File? EVALUATION_REGIONS_BED                    # BED to restrict comparison against TRUTH_VCF to
        File? DV_MODEL_META                             # .meta file for a custom DeepVariant calling model
        File? DV_MODEL_INDEX                            # .index file for a custom DeepVariant calling model
        File? DV_MODEL_DATA                             # .data-00000-of-00001 file for a custom DeepVariant calling model
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 16
        Int MAP_DISK = 10
        Int MAP_MEM = 50
        Int CALL_CORES = 8
        Int CALL_DISK = 40
        Int CALL_MEM = 50
    }
    
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            call extractPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_extract_disk=MAP_DISK,
                    in_extract_mem=MAP_MEM
            }
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call subsetPathNames {
                input:
                    in_path_list_file=extractPathNames.output_path_list_file
            }
        }
    }
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, subsetPathNames.output_path_list_file, written_path_names_file])

    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    call extractReference {
        input:
            in_xg_file=XG_FILE,
            in_vg_container=VG_CONTAINER,
            in_extract_disk=MAP_DISK,
            in_extract_mem=MAP_MEM
    }
    File reference_file = extractReference.reference_file

    call indexReference {
        input:
            in_reference_file=reference_file,
            in_index_disk=MAP_DISK,
            in_index_mem=MAP_MEM
    }
    File reference_index_file = indexReference.reference_index_file
    File reference_dict_file = indexReference.reference_dict_file
    
    
    Array[Pair[File,File]] read_pair_files_list = zip(SAMPLE_INPUT_READ_FILE_1_LIST, SAMPLE_INPUT_READ_FILE_2_LIST)
    scatter (read_pair_set in zip(read_pair_files_list, SAMPLE_NAME_LIST)) {
        Pair[File,File] read_pair_files = read_pair_set.left
        call vgGiraffeDeepVariantWorkflow.vgMultiMap as sample_giraffe_deepvariant_run {
            input:
                INPUT_READ_FILE_1=read_pair_files.left,
                INPUT_READ_FILE_2=read_pair_files.right,
                SAMPLE_NAME=read_pair_set.right,
                VG_CONTAINER=VG_CONTAINER,
                READS_PER_CHUNK=READS_PER_CHUNK,
                GIRAFFE_OPTIONS=GIRAFFE_OPTIONS,
                CONTIGS=CONTIGS,
                PATH_LIST_FILE=pipeline_path_list_file,
                XG_FILE=XG_FILE,
                GBWT_FILE=GBWT_FILE,
                GGBWT_FILE=GGBWT_FILE,
                DIST_FILE=DIST_FILE,
                MIN_FILE=MIN_FILE,
                TRUTH_VCF=TRUTH_VCF,
                TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
                EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
                DV_MODEL_META=DV_MODEL_META,
                DV_MODEL_INDEX=DV_MODEL_INDEX,
                DV_MODEL_DATA=DV_MODEL_INDEX,
                MIN_MAPQ=MIN_MAPQ,
                REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
                SPLIT_READ_CORES=SPLIT_READ_CORES,
                SPLIT_READ_DISK=SPLIT_READ_DISK,
                MAP_CORES=MAP_CORES,
                MAP_DISK=MAP_DISK,
                MAP_MEM=MAP_MEM,
                CALL_CORES=CALL_CORES,
                CALL_DISK=CALL_DISK,
                CALL_MEM=CALL_MEM,
                REFERENCE_FILE=reference_file,
                REFERENCE_INDEX_FILE=reference_index_file,
                REFERENCE_DICT_FILE=reference_dict_file
        }
    }
    output {
        Array[File]? output_vcfeval_evaluation_archive_samples = sample_giraffe_deepvariant_run.output_vcfeval_evaluation_archive
        Array[File]? output_happy_evaluation_archive_samples = sample_giraffe_deepvariant_run.output_happy_evaluation_archive
        Array[File] output_vcf_file_samples = sample_giraffe_deepvariant_run.output_vcf
        Array[File] output_vcf_index_file_samples = sample_giraffe_deepvariant_run.output_vcf_index
        Array[Array[File]] output_indelrealigned_bams_samples = sample_giraffe_deepvariant_run.output_indelrealigned_bams
        Array[Array[File]] output_indelrealigned_bam_indexes_samples = sample_giraffe_deepvariant_run.output_indelrealigned_bam_indexes
    }
}

########################
### TASK DEFINITIONS ###
########################

task extractPathNames {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt
    }
    output {
        File output_path_list_file = "path_list.txt"
    }
    runtime {
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task subsetPathNames {
    input {
        File in_path_list_file
    }

    command <<<
        set -eux -o pipefail

        grep -v _decoy ~{in_path_list_file} | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM > path_list.txt
    >>>
    output {
        File output_path_list_file = "path_list.txt"
    }
    runtime {
        memory: "1 GB"
        disks: "local-disk 10 SSD"
        docker: "ubuntu:20.04"
    }
}

task extractReference {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --extract-fasta \
            --xg ${in_xg_file} > ref.fa
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem
        Int in_index_disk
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_reference_file} ref.fa

        samtools faidx ref.fa

        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict
    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}


