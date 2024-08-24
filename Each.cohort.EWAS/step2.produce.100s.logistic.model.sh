###############Step1 生成批量的101个R文件：logistic模型
#!/usr/bin/bash
for i in $(seq 1 101)
do 
echo " #!/usr/bin/env Rscript

###############读入residual mval
row$i=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step1.mval.files/residual_mvals_block$i.txt\",header=F,sep=\" \")
row.names=row$i[,1]
data$i=t(row$i[,-1])
colnames(data$i)=row.names
file=row$i


##############读入协变量
covariance=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.A1.pheno.txt\",sep=\" \",header=T)

###########################细胞类型
cell=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.A1.Cell.Type.txt\",sep=\" \",header=T)
cell=cell[,1:6]



#################################读入ASA芯片基因型信息主成分
ASA=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.A1.genotype_PCA.txt\",header=T)
genotypePC=ASA[,2:11]


############################性别和年龄
Age=covariance\$Age
Gender=covariance\$Gender


############读入分组
Sample_group=covariance\$Sample_Group


################读入吸烟
Smoking=covariance\$smokingScore


######logistic回归分析
result=c()
for ( i in 1:nrow(file)) { fit = glm(Sample_group~data$i[,i] +Age+Smoking+cell[,1]+cell[,2]+cell[,3]+cell[,4]+cell[,5]+cell[,6],family=binomial); result=as.data.frame(rbind(result,c(colnames(data$i)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c(\"CpG\",\"Estimate\" ,\"SE\" ,\"Pvalue\")}
write.table(result,file=\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/residual_result_block$i\",sep = \" \", row.names = F)" > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/block_residual_logistic_$i.R
done;


##################################Step2.生成批量处理R的shell文件(R to shell)
for i in $(seq 1 101)
do
 echo "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/Rscript /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/block_residual_logistic_$i.R" > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/block_residual_logistic_$i.sh 
done;

#####################################Step3.将所有shell文件生成批量的qusb脚本投递(sh to qsub)
#!/usr/bin/bash

for i in $(seq 1 101)
do 
echo  "qsub -cwd -l vf=50g,num_proc=4 -P P20Z10200N0041_super -binding linear:8 -q st_supermem.q /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/block_residual_logistic_$i.sh" >> /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/batch_shell_to_qusb.sh
done;

sh /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/batch_shell_to_qusb.sh;


rm /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step2.logistic.model.result.files/batch_shell_to_qusb.sh;
