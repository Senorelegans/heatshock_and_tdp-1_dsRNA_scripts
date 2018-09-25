#!/usr/bin/env bash



GTF=input_files/gtf/c_elegans.PRJNA13758.WS258.canonical_geneset.gtf
IN=OK_LFC_FLIPPED_DOG_2ndtime     #This is the output folder from part II
COMPARISON_PATH=${IN}/N2_vs_OK_2ndtime

BAMPATH=BAMS
BEDPATH=BEDS/BEDS_OK
CPUS=4

mkdir -p ${IN}
mkdir -p ${COMPARISON_PATH}


STAROUT=output_from_star


#
##python 1.0_Dogcatcher_flatten_gtf.py --annotation_file_with_path ${GTF}
##
###############Loop through POS and MIN Files ###################
declare -a MIN=($(ls -d ${BEDPATH}/*min.BedGraph | sort))
declare -a POS=($(ls -d ${BEDPATH}/*plu.BedGraph | sort))
L=${#MIN[@]}
X=0
while [  $X -le $L ]; do
    echo "**********"
    echo ${MIN[X]}
    echo ${POS[X]}
    python 2.0_Dogcatcher.py \
    --cpus $CPUS \
    --BedGraph_input_min_strand ${MIN[X]} \
    --BedGraph_input_plu_strand ${POS[X]} \
    --output_prefix  ${IN}/ \
    --annotation_file_with_path ${GTF} \
    --window_size 100
    let X=X+1
done
#
######################################################
python 2.5_Dogcatcher_filter.py \
--filter longest \
--input_prefix ${IN}/ \
--Dogcatcher_plu_strand_list \
N2WT_INP_split_plu.BedGraph \
N2-J2-1__split_plu.BedGraph \
N2-J2-2__split_plu.BedGraph \
N2-J2-3__split_plu.BedGraph \
OK80-INP_split_plu.BedGraph \
OK8-J2-1_split_plu.BedGraph \
OK8-J2-2_split_plu.BedGraph \
OK8-J2-3_split_plu.BedGraph \
--Dogcatcher_min_strand_list \
N2WT_INP_split_min.BedGraph \
N2-J2-1__split_min.BedGraph \
N2-J2-2__split_min.BedGraph \
N2-J2-3__split_min.BedGraph \
OK80-INP_split_min.BedGraph \
OK8-J2-1_split_min.BedGraph \
OK8-J2-2_split_min.BedGraph \
OK8-J2-3_split_min.BedGraph \
--output_prefix ${COMPARISON_PATH}/
#
#python 3.2_Create_R_subread_ENRICHMENT.py \
#--annotation_file_with_path ${GTF} \
#--control_TOTAL_BAM_list \
#${BAMPATH}/N2WT_INP_.bam \
#--treatment_TOTAL_BAM_list \
#${BAMPATH}/OK803-INP__.bam \
#--control_J2_BAM_list \
#${BAMPATH}/N2-J2-1_.bam \
#${BAMPATH}/N2-J2-2_.bam \
#${BAMPATH}/N2-J2-3_.bam \
#--treatment_J2_BAM_list \
#${BAMPATH}/OK803-J2-1_.bam \
#${BAMPATH}/OK803-J2-2_.bam \
#${BAMPATH}/OK803-J2-3_.bam \
#--input_R_template_file initial_R_subread_TEMPLATE.R \
#--input_prefix ${IN} \
#--output_prefix ${IN}/initial_Rsubread \
#--cpus ${CPUS} \
#--padj 0.05
#
##Run the created R script
#R CMD BATCH ${IN}/initial_Rsubread/Rsubread_initial.R ${IN}/initial_Rsubread/Rsubread_initial.R.out
#
#python 3.3_INP_J2_normalization_library_rRNA_mapped.py \
#--control_TOTAL_BAM_list \
#${BAMPATH}/N2WT_INP_.bam \
#--treatment_TOTAL_BAM_list \
#${BAMPATH}/OK803-INP__.bam \
#--control_J2_BAM_list \
#${BAMPATH}/N2-J2-1_.bam \
#${BAMPATH}/N2-J2-2_.bam \
#${BAMPATH}/N2-J2-3_.bam \
#--treatment_J2_BAM_list \
#${BAMPATH}/OK803-J2-1_.bam \
#${BAMPATH}/OK803-J2-2_.bam \
#${BAMPATH}/OK803-J2-3_.bam \
#--control_TOTAL_STAROUT_list \
#${STAROUT}/N2WT_INP_TLog.final_rRNA.out \
#--treatment_TOTAL_STAROUT_list \
#${STAROUT}/OK803-INP__TLog.final_rRNA.out \
#--control_J2_STAROUT_list \
#${STAROUT}/N2-J2-1__TLog.final_rRNA.out \
#${STAROUT}/N2-J2-2__TLog.final_rRNA.out \
#${STAROUT}/N2-J2-3__TLog.final_rRNA.out \
#--treatment_J2_STAROUT_list \
#${STAROUT}/OK803-J2-1_TLog.final_rRNA.out \
#${STAROUT}/OK803-J2-2_TLog.final_rRNA.out \
#${STAROUT}/OK803-J2-3_TLog.final_rRNA.out \
#--input_prefix ${IN}/initial_Rsubread
#
#
#python 3.4_Create_R_DESeq2_ENRICHMENT.py \
#--input_R_template_file initial_DEseq2_TEMPLATE.R \
#--output_prefix ${IN}/initial_Rsubread \
#--cpus ${CPUS} \
#
##Run the created R script
#R CMD BATCH ${IN}/initial_Rsubread/DESeq2_initial.R ${IN}/initial_Rsubread/DESeq2_initial.R.out
#echo runningR
#
#
#python 4.0_Dogcatcher_Rsubread_DESeq2_ENRICHMENT.py \
#--annotation_file_with_path ${GTF} \
#--input_prefix ${COMPARISON_PATH} \
#--input_prefix_initial ${IN}/initial_Rsubread \
#--output_prefix ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes \
#--padj 0.05
#echo finished nonsig
#
#
####Run the Rsubread scripts from DOG gtf's with non-significant genes
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsigRsubread.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsigRsubread.R.out
#echo finished R scripts plu DOG
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsigRsubread.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsigRsubread.R.out
#echo finished R scripts min DOG
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsigRsubread.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsigRsubread.R.out
#echo finished R scripts plu ADOG
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsigRsubread.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsigRsubread.R.out
#echo finished R scripts min ADOG
#
###Normalize DOGS
#python 4.1_INP_J2_normalization_library_rRNA_mapped_2ndtime.py \
#--control_TOTAL_BAM_list \
#${BAMPATH}/N2WT_INP_.bam \
#--treatment_TOTAL_BAM_list \
#${BAMPATH}/OK803-INP__.bam \
#--control_J2_BAM_list \
#${BAMPATH}/N2-J2-1_.bam \
#${BAMPATH}/N2-J2-2_.bam \
#${BAMPATH}/N2-J2-3_.bam \
#--treatment_J2_BAM_list \
#${BAMPATH}/OK803-J2-1_.bam \
#${BAMPATH}/OK803-J2-2_.bam \
#${BAMPATH}/OK803-J2-3_.bam \
#--control_TOTAL_STAROUT_list \
#${STAROUT}/N2WT_INP_TLog.final_rRNA.out \
#--treatment_TOTAL_STAROUT_list \
#${STAROUT}/OK803-INP__TLog.final_rRNA.out \
#--control_J2_STAROUT_list \
#${STAROUT}/N2-J2-1__TLog.final_rRNA.out \
#${STAROUT}/N2-J2-2__TLog.final_rRNA.out \
#${STAROUT}/N2-J2-3__TLog.final_rRNA.out \
#--treatment_J2_STAROUT_list \
#${STAROUT}/OK803-J2-1_TLog.final_rRNA.out \
#${STAROUT}/OK803-J2-2_TLog.final_rRNA.out \
#${STAROUT}/OK803-J2-3_TLog.final_rRNA.out \
#--input_prefix ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes
#
#
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsigDESeq2.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsigDESeq2.R.out
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsigDESeq2.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsigDESeq2.R.out
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsigDESeq2.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsigDESeq2.R.out
#R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsigDESeq2.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsigDESeq2.R.out
#echo finished R scripts

#Filter out significant runon/runin from DESeq2 and match with the csv so you can get biotypes etc.
python 5.0_filter_sig_DESeq2.py \
--annotation_file_with_path ${GTF} \
--input_prefix ${COMPARISON_PATH} \
--output_prefix ${COMPARISON_PATH}/FINAL_OUT \
--input_DOG_DESeq2_plu_sense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv \
--input_DOG_DESeq2_min_sense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv \
--input_ADOG_DESeq2_plu_antisense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv \
--input_ADOG_DESeq2_min_antisense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv \
--padj 0.05
#wait
