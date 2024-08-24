#!/usr/bin/bash

############以下为对结果进行初步的格式处理得到linear_result_all.txt文件的最终格式为:"CpG  Estimate  SE  Pval"
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD//02.Diffrential_CpG_analysis/step3.annotation.results.files/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD//02.Diffrential_CpG_analysis/step2.logistic.model.result.files/"

rm ${wd}linear_result_all.txt;
cat ${wd2}residual_result_block* >>  ${wd}linear_result_all.txt;
sed -i '/Estimate/d' ${wd}linear_result_all.txt; 
sed -i "s/\"/ /g"  ${wd}linear_result_all.txt;
sed -i "1iCpG   Estimate   SE   Pval" ${wd}linear_result_all.txt;
sed -i "s/^ //g" ${wd}linear_result_all.txt;

#############以下为可选操作
##########排序查看P值(可选操作）
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD//02.Diffrential_CpG_analysis/step3.annotation.results.files/"

sort -k 4  -g  ${wd}linear_result_all.txt > ${wd}GDPH.Sort.results &&##排序

############筛选过滤得到与兴趣变量相关的位点（可选操作）

#awk -F " " '{if ($4 <= 9e-08) print $0}' ${wd}linear_result_all.txt > ${wd}site_with_ESRD.txt;#筛选与Sentrix_ID相关的位点

#############################以下为曼哈顿图所需格式处理
########提取曼哈顿图所需的前三行注释文件
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD//02.Diffrential_CpG_analysis/step3.annotation.results.files/"

#awk -F " " '{print $1,$2,$3}' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/EWAS_20210812/annotation_EPIC.txt > 3_col_annotation_EPIC.txt;sed -i '1d' 3_col_annotation_EPIC.txt;sed -i "1iSNP CHR BP" 3_col_annotation_EPIC.txt;sed -i 's/\"//g' 3_col_annotation_EPIC.txt


####合并整理为可以做曼哈顿图的文件"file_for_manhattan_plot.txt"
qsub -cwd -l vf=50g,num_proc=4 -P P20Z10200N0041_super -binding linear:8 -q st_supermem.q /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD//02.Diffrential_CpG_analysis/Important_files/merge.sh &&

#merge.sh的目的：perl /zfssz2/ST_MCHRI/BIGDATA/USER/zhouxiaohong/EWAS_20210812/SG_20211110/SG_719samples_ESRD_DM/02.correct_methylation_value/tmp_correct_mval/04.annotation_all_results/merge.pl  linear_result_all.txt 3_col_annotation_EPIC.txt  file_for_manhattan_plot.txt;生成文件file_for_manhattan_plot.txt

sed -i "s/chr//g" ${wd}file_for_manhattan_plot.txt;

###记得把file_for_manhattan_plot.txt重命名

