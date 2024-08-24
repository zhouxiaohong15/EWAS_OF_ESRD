#!usr/bin/bash
for i in $(seq 1 100)
do 
sed -n "`echo $[7776*($i-1)+2]`,`echo  $[7776*($i-1)+1+7766]p`" /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.A1.residual_mvals.txt  > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step1.mval.files/residual_mvals_block$i.txt
done;

tail -n 32 /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.A1.residual_mvals.txt > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step1.mval.files/residual_mvals_block101.txt;
