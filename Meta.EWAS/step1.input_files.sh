
#####################This manuscript was aimed to generate input files for step2(step2.Bacon_adjusted.R), "AHPL.GDPH.logistic.result.all.txt", of which the first column is EWAS results in Singapore, the last column is EWAS results in GDPH.

wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis//00.SGandSY_EWAS_results/";

cp /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartB.AHPL_SG_EWAS/1.Raw_QC_Normalized/02.Diffrential_CpG_analysis/step3.annotation.results.files/linear_result_all.txt ${wd}SG_linear_result_all.txt;
sed -i "s/\s+/\t/g" ${wd}SG_linear_result_all.txt;

cp  /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step3.annotation.results.files/linear_result_all.txt ${wd}SY_linear_result_all.txt;
sed -i "s/\s+/\t/g" ${wd}SY_linear_result_all.txt;

########################################去除结果中缺失值大于5%的探针
sh /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/Remove_0.05Missing_Probe.sh;
sed -i "s/\s+/\t/g" ${wd}SG_linear_result_all.txt;
sed -i "s/\s+/\t/g" ${wd}SY_linear_result_all.txt;

####################合并两个队列的数据
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis//00.SGandSY_EWAS_results/";
perl ${wd}merge_SY_SG.pl  ${wd}SG_linear_result_all.txt ${wd}SY_linear_result_all.txt ${wd}AHPL.GDPH.logistic.result.all.txt;

sed -i "1d" ${wd}AHPL.GDPH.logistic.result.all.txt;
sed -i "1icgname Estimate.SY Std.SY P.SY Estimate.SG Std.SG P.SG" ${wd}AHPL.GDPH.logistic.result.all.txt;
