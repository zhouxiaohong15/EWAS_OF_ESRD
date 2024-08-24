
#################################变量“data”为/jdfssz1/ST_HEALTH/P20H10200N0033/Zhouxiaohong/Meta_analysis/00.SGandSY_EWAS_results/AHPL.GDPH.logistic.result.all.txt。它是来自脚本“step1.input_files.sh”的结果。格式为
#cgname Estimate Std P Estimate Std P
#cg07364396 2.06680686930211 0.392231496288757 1.36904208043052e-07 2.65628849064727 0.443357129814798 2.08161792941148e-09
#前三列是省医的结果，后三列是新加坡的结果。

data=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/AHPL.GDPH.logistic.result.all.txt",header=T,sep=" ")
len=nrow(data)
data$ID=c(1:len)

######将data赋值为另一变量，用于结尾作QQ图
dataqq=data


#############es为SG和SY队列的estimate,se为为SG和SY队列的standard error
es=cbind(data[,2],data[,5])
se=cbind(data[,3],data[,6])

#############进行bacon校正统计量
library(bacon)
bc <- bacon(NULL, es, se)
##################提取校正后的统计量（可选部分，不一定非要执行）
Bacon.se=se(bc)
Bacon.estimate=es(bc)
Bacon.pval=pval(bc)
###############重构"data"变量，列名为cgname    Bacon.Estimate_sg  Bacon.se_sg   Bacon.pval_sg  Bacon.Estimate_sy  Bacon.se_sy  Bacon.pval_sy ID;可以将data保存起来（可选部分，不一定非要执行）

data=data.frame(cgname=data$cgname,Bacon.Estimate_sy=Bacon.estimate[,1],Bacon.se_sy=Bacon.se[,1],Bacon.pval_sy=Bacon.pval[,1],Bacon.Estimate_sg=Bacon.estimate[,2],Bacon.se_sg=Bacon.se[,2],Bacon.pval_sg=Bacon.pval[,2],ID=data$ID)

row=nrow(data)
row
write.table(data,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/Bacon.AHPL.GDPH.logistic.result.txt")



########################制作A1和A3队列QQ Plot 
data=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/AHPL.GDPH.logistic.result.all.txt",header=T,sep=" ")
len=nrow(data)
data$ID=c(1:len)
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results")
#########A1 QQ Plot
es.A1=data[,2]
se.A1=data[,3]

bc.A1=bacon(NULL,es.A1,se.A1)

P.alone=plot(bc.A1,type="qq") +theme_classic()
ggsave("GDPH.A1.QQplot.jpg", plot = P.alone, device = "jpeg", dpi = 600)


#########A3 QQ Plot
es.A3=data[,5]
se.A3=data[,6]

bc.A3=bacon(NULL,es.A3,se.A3)

P.alone=plot(bc.A3,type="qq") +theme_classic()
ggsave("AHPL.A3.QQplot.jpg", plot = P.alone, device = "jpeg", dpi = 600)



#####保存SY和SG经过bacon校正后的inflation factor;因为输入bc，自动会输出相关的iflation bias等信息（参考https://www.bioconductor.org/packages/devel/bioc/vignettes/bacon/inst/doc/bacon.html#4_Fixed-effect_meta-analysis）
bc
save(bc,file="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon_adjusted_inflation.Rdata")
#############################Perform fixed-effect meta-analysis and the inspection of results. 该结果包含
bcm <- meta(bc)
save(bcm,file="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Adjusted.FE.Meta.results.Rdata")

######meta结果 QQ plot
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results")
P.meta=plot(bcm,type="qq") + theme_classic()
ggsave("Meta.QQplot.jpg", plot = P.meta, device = "jpeg", dpi = 600)

load("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Adjusted.FE.Meta.results.Rdata")
str(bcm)
row=dim(bcm@ teststatistics)[1] ######获得行数
inflation(bcm)
####打印出所需结果变量“n"，格式为数据框："eff.size.meta","std.err.meta","pval.adj.meta","pval.org.meta","tstat.meta","eff.size.bacon.SY","eff.size.bacon.SG","std.err.bacon.SY","std.err.bacon.SG","pval.bacon.SY","pval.bacon.SG","tstat.SY","tstat.SG"
n=topTable(bcm,n=row) 
n=as.data.frame(n)
names(n)=c("eff.size.meta","std.err.meta","pval.adj.meta","pval.org.meta","tstat.meta","eff.size.bacon.SY","eff.size.bacon.SG","std.err.bacon.SY","std.err.bacon.SG","pval.bacon.SY","pval.bacon.SG","tstat.SY","tstat.SG")

##################################由于n数据框没有CpG位点的名称，因此需要借助文件“/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/Bacon.AHPL.GDPH.logistic.result.txt”，将CpG名称加入到n数据框中。
ba=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/00.SGandSY_EWAS_results/Bacon.AHPL.GDPH.logistic.result.txt",header=T)
ba=ba[order(ba$Bacon.Estimate_sy),]
n=n[order(n$eff.size.bacon.SY),]
identical(round(ba$Bacon.Estimate_sy,3),round(n$eff.size.bacon.SY,3))

###########如果为identical(round(ba$Bacon.Estimate_sy,3),round(n$eff.size.bacon.SY,3))结果TRUE，则如下
new=cbind(ba$cgname,n)
new=new[order(new$pval.adj.meta),]
write.table(new,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Adjusted.FE.Meta.results.txt",row.names=F)

#########################将bacon校正后的SY SG Meta EWAS结果整理为曼哈顿图的文件格式
anno=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartB.AHPL_SG_EWAS/1.Raw_QC_Normalized/02.Diffrential_CpG_analysis/Important_files/3_col_annotation_EPIC.txt",header=T)
EWAS=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Adjusted.FE.Meta.results.txt",header=T)
names(EWAS)[1]="SNP"
EWAS_new=merge(EWAS,anno,by="SNP")

man_Meta=EWAS_new[,c(1,15,16,5)]
man_SY=EWAS_new[,c(1,15,16,11)]
man_SG=EWAS_new[,c(1,15,16,12)]
names(man_Meta)[4]="P"
names(man_SY)[4]="P"
names(man_SG)[4]="P"
man_Meta$CHR=gsub('[chr]','',man_Meta$CHR)
man_SY$CHR=gsub('[chr]','',man_SY$CHR)
man_SG$CHR=gsub('[chr]','',man_SG$CHR)

man_SG=man_SG[order(man_SG$CHR,man_SG$BP),]
man_SY=man_SY[order(man_SY$CHR,man_SY$BP),]
man_Meta=man_Meta[order(man_Meta$CHR,man_Meta$BP),]

write.table(man_SG,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.SG.manhattan.txt",row.names=F)
write.table(man_SY,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.SY.manhattan.txt",row.names=F)
write.table(man_Meta,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Meta.manhattan.txt",row.names=F)

##########################制作SY SG Meta的QQ图
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/")
load("Bacon.Adjusted.FE.Meta.results.Rdata")
library(ggplot2)

########################################################################GDPH QQ
Pval=as.data.frame(pval(bcm))
head(Pval)
names(Pval)=c("GDPH","AHPL","Meta")
pvals=Pval$GDPH
rate <- 2.30258509299405 #####该数值是通过命令“a=plot(bcm,type="qq")；a$layer;”所提取

#########提取期望值，也就是QQ图的横坐标
theoretical_quantiles <- qexp(ppoints(length(pvals)), rate = 2.30258509) ##期望值
pvals=Pval$GDPH
empirical_quantiles <- -(log10(pvals))

#######形成数据框，横坐标期望的P值，纵坐标是实际P值，并对他们进行排序
data = data.frame(var1=theoretical_quantiles[order(theoretical_quantiles)], var2=empirical_quantiles[order(empirical_quantiles)]) ###实际值
head(data)

#################保存GDPH结QQ图
GDPH= ggplot(data, aes(x = var1, y = var2))+ geom_point() + geom_abline(intercept = 0, slope = 1, color = "red",linewidth=2)+ theme_classic() +labs(x = "Expected -log10(Pvalue)",   y = "Observed -log10(Pvalue)",     title = "QQ Plot of GDPH") + theme(    text = element_text(size=15,face="bold"),  axis.title = element_text(size=20, face="bold"), axis.text = element_text(size=15,face="bold"),  axis.line = element_line(linewidth=1.5)  ) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=25))+ coord_cartesian(ylim = c(0, 30))
ggsave("QQ.GDPH.jpg", plot = AHPL, device = "jpeg", dpi = 600)



########################################################################GDPH QQ
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/")
load("Bacon.Adjusted.FE.Meta.results.Rdata")
library(ggplot2)

Pval=as.data.frame(pval(bcm))
head(Pval)
names(Pval)=c("GDPH","AHPL","Meta")
pvals=Pval$AHPL
rate <- 2.30258509299405 #####该数值是通过命令“a=plot(bcm,type="qq")；a$layer;”所提取

#########提取期望值，也就是QQ图的横坐标
theoretical_quantiles <- qexp(ppoints(length(pvals)), rate = 2.30258509) ##期望值
pvals=Pval$AHPL
empirical_quantiles <- -(log10(pvals))

#######形成数据框，横坐标期望的P值，纵坐标是实际P值，并对他们进行排序
data = data.frame(var1=theoretical_quantiles[order(theoretical_quantiles)], var2=empirical_quantiles[order(empirical_quantiles)]) ###实际值
head(data)

#################保存AHPL结QQ图
AHPL= ggplot(data, aes(x = var1, y = var2))+ geom_point() + geom_abline(intercept = 0, slope = 1, color = "red",linewidth=2)+ theme_classic() +labs(x = "Expected -log10(Pvalue)",   y = "Observed -log10(Pvalue)",     title = "QQ Plot of AHPL") + theme(    text = element_text(size=15,face="bold"),  axis.title = element_text(size=20, face="bold"), axis.text = element_text(size=15,face="bold"),  axis.line = element_line(linewidth=1.5)  ) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=25))+ coord_cartesian(ylim = c(0, 30))
ggsave("QQ.AHPL.jpg", plot = AHPL, device = "jpeg", dpi = 600)


#############################################################Meta QQ图
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/")
load("Bacon.Adjusted.FE.Meta.results.Rdata")
library(ggplot2)

pvals=Pval$Meta
rate <- 2.30258509299405 #####该数值是通过命令“a=plot(bcm,type="qq")；a$layer;”所提取

#########提取期望值，也就是QQ图的横坐标
theoretical_quantiles <- qexp(ppoints(length(pvals)), rate = 2.30258509) ##期望值
pvals=Pval$Meta
empirical_quantiles <- -(log10(pvals))

#######形成数据框，横坐标期望的P值，纵坐标是实际P值，并对他们进行排序
data = data.frame(var1=theoretical_quantiles[order(theoretical_quantiles)], var2=empirical_quantiles[order(empirical_quantiles)]) ###实际值
head(data)

#################保存Meta结QQ图
Meta= ggplot(data, aes(x = var1, y = var2))+ geom_point() + geom_abline(intercept = 0, slope = 1, color = "red",linewidth=2)+ theme_classic() +labs(x = "Expected -log10(Pvalue)",   y = "Observed -log10(Pvalue)",     title = "QQ Plot of Meta") + theme(    text = element_text(size=15,face="bold"),  axis.title = element_text(size=20, face="bold"), axis.text = element_text(size=15,face="bold"),  axis.line = element_line(linewidth=1.5)  ) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=25))+ coord_cartesian(ylim = c(0, 30))
ggsave("QQ.Meta.jpg", plot = Meta, device = "jpeg", dpi = 600)

#################################提取显著性位点的结果
man=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Adjusted.FE.Meta.results.txt",header=T)
head(man)
dim(man)
0.01/771388
DMP=subset(man,pval.bacon.SY < 0.01 & pval.bacon.SG < 0.01 & pval.org.meta < 1.296364e-08)
DMP$direction=DMP$eff.size.bacon.SY*DMP$eff.size.bacon.SG
DMP0=subset(DMP, direction > 0)
dim(DMP0)
write.table(DMP0,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.meta.DMP.txt",row.names=F)




#########################对SY SG Meta制作QQ plot，注意：不可以用正态分布的方式作图，一定要用bacon的先验分布。
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/")
plot(bcm, type="qq")+ theme_classic()
dev.off()


############################对DMP 的meta结果进行异质性分析（固定效应 模型）
library(metafor)
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.meta.DMP.txt",header=T,sep=" ")
data=DMP[,c(1,7,9,11,8,10,12)]
row=nrow(data)
row
data$ID=c(1:row)

colnames(data)=c("cgname", "Estimate.sy", "Std.sy", "P.sy", "Estimate.sg","Std.sg",  "P.sg", "ID")
library(dplyr)
data=mutate(data,Estimate.sy=as.numeric(Estimate.sy),Std.sy=as.numeric(Std.sy),P.sy=as.numeric(P.sy),Estimate.sg=as.numeric(Estimate.sg),Std.sg=as.numeric(Std.sg),P.sg=as.numeric(P.sg))

row=nrow(data)
row


###############################meta: Fix effect model
result=c()
es=t(data[,c(2,5)])
se=t(data[,c(3,6)])
########file的格式为 
###                  ID    Estimate       Std
#################### SG -0.85690053 0.4220631
#################### SY  0.07684989 0.6519630

for ( i in 1:row) { file=data.frame(ID=c("SG","SY"),Estimate=es[,i],Std=se[,i]);meta=rma(yi=Estimate,data=file,sei=Std,method="FE",weighted=TRUE);result=rbind(result,c(summary(meta)$beta,summary(meta)$se,summary(meta)$pval,summary(meta)$I2,summary(meta)$ QE,summary(meta)$ QEp))}

result=as.data.frame(result)
str(result)

colnames(result)=c("meta.beta","meta.se","meta.pval","meta.I2","meta.Q","meta.Qpval")
all=cbind(data,result)
all[1:2,]
str(all)

bacon_pval=DMP[,c(1,4,5)]
names(all)[1]="cgname"
names(bacon_pval)[1]="cgname"
all=merge(all,bacon_pval,by="cgname")
all=all[,c(1:11,15,16,12:14)]
all=all[order(all$pval.adj.meta),]
Fix=all

write.table(Fix,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/MetaBacon.SYSG.FixModel.Result.txt",row.names=F)



############################对DMP 的meta结果进行异质性分析（随机效应 模型）
library(metafor)
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.meta.DMP.txt",header=T,sep=" ")
data=DMP[,c(1,7,9,11,8,10,12)]
row=nrow(data)
row
data$ID=c(1:row)

colnames(data)=c("cgname", "Estimate.sy", "Std.sy", "P.sy", "Estimate.sg","Std.sg",  "P.sg", "ID")
library(dplyr)
data=mutate(data,Estimate.sy=as.numeric(Estimate.sy),Std.sy=as.numeric(Std.sy),P.sy=as.numeric(P.sy),Estimate.sg=as.numeric(Estimate.sg),Std.sg=as.numeric(Std.sg),P.sg=as.numeric(P.sg))

row=nrow(data)
row


###############################meta
result=c()
es=t(data[,c(2,5)])
se=t(data[,c(3,6)])
########file的格式为 
###                  ID    Estimate       Std
#################### SG -0.85690053 0.4220631
#################### SY  0.07684989 0.6519630

for ( i in 1:row) { file=data.frame(ID=c("SG","SY"),Estimate=es[,i],Std=se[,i]);meta=rma(yi=Estimate,data=file,sei=Std,method="REML");result=rbind(result,c(summary(meta)$beta,summary(meta)$se,summary(meta)$pval,summary(meta)$I2,summary(meta)$ QE,summary(meta)$ QEp))}

result=as.data.frame(result)
str(result)

colnames(result)=c("meta.beta","meta.se","meta.pval","meta.I2","meta.Q","meta.Qpval")
all=cbind(data,result)
all[1:2,]
str(all)

bacon_pval=DMP[,c(1,4,5)]
names(all)[1]="cgname"
names(bacon_pval)[1]="cgname"
all=merge(all,bacon_pval,by="cgname")
all=all[,c(1:11,15,16,12:14)]
all=all[order(all$pval.adj.meta),]

Random=all

#########################################################整合固定效应和随机效应模型的结果
Fix=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/MetaBacon.SYSG.FixModel.Result.txt",header=T,sep=" ")

Fix$meta.re.pval=Random$meta.pval

Fix=Fix[,-11]

for(i in 1:nrow(Fix)) {
  # 如果meta.I2 大于等于50
  if(Fix$meta.Qpval[i] >= 0.1) {
    Fix$P.final[i] <- Fix$pval.org.meta[i]
  }
  # 如果meta.I2 小于50
  else if(Fix$meta.Qpval[i] < 0.1) {
    Fix$P.final[i] <- Fix$meta.re.pval[i]
  }
}


write.table(Fix,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Meta.DMP.FEandRE.Result.txt",row.names=F)
