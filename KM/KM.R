setwd("../mRNA")
###整理数据
rt=read.table("fpkm.txt",sep="\t",header=T,check.names=F)                #读入并整理表达矩阵
rt=t(rt)                                                                 #转置
write.csv(rt,"fpkm1.csv")
rownames(rt)=rt[,1]
rt=read.table("fpkm1.txt",sep="\t",header=T,check.names=F)
cli=read.table("time.txt",sep="\t",check.names=F,header=T,row.names=1)   #读入并整理生存数据
cli[,3]=rownames(cli)
colnames(cli)=c("os","stata","id")
merge_os=merge(rt,cli,by="id",all=F)                                     #合并表达数据与生存数据
write.csv(merge_os,"OS.csv")

###开始分析
setwd("./")
os=read.table("OS.txt",sep="\t",header=T,check.names=F,row.name=1)       #带有表达量,生存时间和状态
os$futime=os$futime/30    #将时间坐标改为月
library(survival)
library(survminer)
col=c("red","blue")    
outTab=data.frame()
{
  group <- numeric(length(os[,3]))   
  group
  group[os[,3]>58.3062075]<-"high" 
  group[os[,3]<44.0878696]<-"low"    #以平均值上下15%浮动为分组标准
  diff=survdiff(Surv(futime, fustat) ~group,data = os)
  pValue=1-pchisq(diff$chisq,df=1)
  outVector=cbind(gene,pValue)
  outTab=rbind(outTab,outVector)
  if(pValue<0.0001){
    pValue="p<0.0001"
  }else{
    pValue=paste0("p=",sprintf("%.04f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = os)
  surPlot=ggsurvplot(fit,            #画图               
                     data=os,
                     pval=pValue,
                     pval.size=6,
                     legend.labs=c("high","low"),
                     legend.title=paste0(gene," levels"),
                     font.legend=12,
                     xlab="Time(mouths)",
                     palette=col,
                     break.time.by = 20,
                     conf.int=T,
                     fontsize=4.5,
                     risk.table=TRUE,
                     ylab="Overall survival",
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0(gene,".pdf"),onefile = FALSE,
      width = 6,         #
      height =5)         #
  print(surPlot)
  dev.off()
}
write.table(outTab,file="OS.result.xls",sep="\t",row.names=F,quote=F)
