#!/bin/bash

# Set the variables
feature_file_name="/data/illumina_pipeline/scripts/feature_files/yeast/replication_origins/eaton_origins/"
feature_file_name="${feature_file_name}2010_genome_research_geo_submission/nrACS/"
feature_file_name="${feature_file_name}eaton_nrACS_GSM424494_sacCer2_sgd_r61_feature_file.csv"
output_file_name="${feature_file_name:0:${#feature_file_name}-4}_2.csv"

pwm_file_name="/data/illumina_pipeline/scripts/feature_files/yeast/replication_origins/pwm/ACS.pwm"
bg_model_file_name="/data/illumina_pipeline/scripts/feature_files/yeast/genome/sgd_R61_ucsc_sacCer2/yeast_1mer_to_6mer_background_nucleotide_freq.tsv"
genome_fasta_header="/data/illumina_pipeline/scripts/feature_files/yeast/genome/sgd_R61_ucsc_sacCer2/sacCer2_sgdR61_chr"

# Set the parameters
program="FeatureMotifScore"

# Set the jar file path
jar_file_path="/data/illumina_pipeline/scripts/java_scripts/jar_files/motif-processing.jar"

# Execute the script
java -jar $jar_file_path $program $feature_file_name $output_file_name $pwm_file_name $bg_model_file_name $genome_fasta_header
