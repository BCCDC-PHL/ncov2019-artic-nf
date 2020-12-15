#!/bin/bash
set -eo pipefail

export PATH=/opt/conda/bin:$PATH

ref_storage_dir="${GITHUB_WORKSPACE}/composite_ref"
mkdir ${ref_storage_dir}

# get the GRCh38 human genome
# as per https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" > ${ref_storage_dir}/GRC38_no_alt_analysis_set.fna.gz
gunzip ${ref_storage_dir}/GRC38_no_alt_analysis_set.fna.gz

curl -s "https://raw.githubusercontent.com/BCCDC-PHL/artic-ncov2019/master/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta" > ${ref_storage_dir}/nCoV-2019.reference.fasta

# create composite reference of human and virus for competitive bwt mapping 
# based host removal
cat ${ref_storage_dir}/GRC38_no_alt_analysis_set.fna ${ref_storage_dir}/nCoV-2019.reference.fasta > ${ref_storage_dir}/composite_human_viral_reference.fna
bwa index ${ref_storage_dir}/composite_human_viral_reference.fna

