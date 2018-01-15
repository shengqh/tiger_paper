library(reshape2)
library(ggplot2)

splitmirna<-function(dat, pattern){
  multi<-grepl(pattern, rownames(dat))
  multiv<-dat[multi,]
  dat<-dat[!multi,]
  idx<-1
  for(idx in c(1:nrow(multiv))){
    names<-unlist(strsplit(rownames(multiv)[idx], pattern))
    name<-names[1]
    for(name in names){
      df<-multiv[idx,]
      df<-df / length(names)
      rownames(df)<-name
      if(name %in% rownames(dat)){
        dat[name,] = dat[name,] + df
      }else{
        dat<-rbind(dat, df)
      }
    }
  }
  return(dat)
}

chimira<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.chimira.miRNA.csv", header=T, row.names=1)
rownames(chimira)<-gsub("miR","mir",rownames(chimira))

excerptrMirnaFile="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.excerptr.miRNA.fixed.csv"
if(!file.exists(excerptrMirnaFile)){
  excerptr<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.excerptr.miRNA.csv", header=T, row.names=1)
  colnames(excerptr)<-gsub("sample_","",colnames(excerptr))
  colnames(excerptr)<-gsub("_fastq","",colnames(excerptr))
  excerptr<-splitmirna(excerptr, "\\|")
  excerptr["mmu-miR-28c",] = excerptr["mmu-miR-28c",] + excerptr["mmu-mir-28c",]
  excerptr["mmu-miR-378d",] = excerptr["mmu-miR-378d",] + excerptr["mmu-mir-378d",]
  enames<-c("mmu-mir-28c", "mmu-mir-378d")
  excerptr<-excerptr[!(rownames(excerptr) %in% enames),]
  rownames(excerptr)<-gsub("miR","mir",rownames(excerptr))
  write.csv(excerptr, file=excerptrMirnaFile)
}else{
  excerptr<-read.csv(excerptrMirnaFile, row.names=1)
}

mirgeMirnaFile<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRge.miRNA.fixed.csv"
if(!file.exists(mirgeMirnaFile)){
  mirge<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRge.miRNA.csv", header=T, row.names=1)
  mirge<-splitmirna(mirge, "/")
  rownames(mirge)<-gsub("miR","mir",rownames(mirge))
  write.csv(mirge, file=mirgeMirnaFile)
}else{
  mirge<-read.csv(file=mirgeMirnaFile, header=T, row.names=1)
}

oasisMirnaFile<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.oasis.miRNA.fixed.csv"
if(!file.exists(oasisMirnaFile)){
  oasis<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.oasis.miRNA.csv", header=T, row.names=1)
  oasis<-splitmirna(oasis, "__")
  rownames(oasis)<-gsub("miR","mir",rownames(oasis))
  write.csv(oasis, file=oasisMirnaFile)
}else{
  oasis<-read.csv(file=oasisMirnaFile, header=T, row.names=1)
}

tigerMirnaFile<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.tiger.miRNA.fixed.csv"
if(!file.exists(tigerMirnaFile)){
  tiger<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.miRNA.count", sep="\t", header=T, row.names=1)
  tiger<-tiger[,colnames(mirge)]
  tiger<-splitmirna(tiger, ";")
  rownames(tiger)<-gsub("miR","mir",rownames(tiger))
  tiger<-tiger[apply(tiger, 1, sum) > 0,]
  write.csv(tiger, file=tigerMirnaFile)
}else{
  tiger<-read.csv(file=tigerMirnaFile, header=T, row.names=1)
}

x<-chimira
y<-excerptr
calcCorrelation <- function(x, y){
  commonNames<-rownames(x)[rownames(x) %in% rownames(y)]
  samples<-sort(colnames(x))
  cx<-x[commonNames, samples]
  cy<-y[commonNames, samples]
  sample<-samples[1]
  do.call(rbind, lapply(samples,function(sample){
    ccx<-cx[,sample,drop=T]
    ccy<-cy[,sample,drop=T]
    sumxy<-(ccx + ccy) > 0
    ccx<-ccx[sumxy]
    ccy<-ccy[sumxy]
    ct<-cor.test(ccx, ccy, method="spearman")
    data.frame("Sample"=sample, "SpearmanR"=ct$estimate, "SpearmanP"=ct$p.value)
  }))
}

xname<-"chimira"
yname<-"excerptr"
getPoint <- function(x, y, xname, yname){
  commonNames<-rownames(x)[rownames(x) %in% rownames(y)]
  samples<-sort(colnames(x))
  cx<-x[commonNames, samples]
  cy<-y[commonNames, samples]
  sample<-samples[1]
  do.call(rbind, lapply(samples,function(sample){
    ccx<-cx[,sample,drop=T]
    ccy<-cy[,sample,drop=T]
    sumxy<-(ccx + ccy) > 0
    ccx<-ccx[sumxy]
    ccy<-ccy[sumxy]
    data.frame("Sample"=sample, "X"=ccx, "Y"=ccy, "Xname" = xname, "Yname" = yname)
  }))
}

datalistnames<-list("Chimira", "excerptr", "miRge", "OASIS", "TIGER")
datalists<-list(chimira, excerptr, mirge, oasis, tiger)

finalpoint<-NULL
finaltable<-NULL
ix<-1
iy<-2
for (ix in c(1:(length(datalistnames) - 1))){
  ixname<-datalistnames[[ix]]
  ixdata<-datalists[[ix]]
  for (iy in c((ix+1):length(datalistnames))){
    iyname<-datalistnames[[iy]]
    iydata<-datalists[[iy]]
    
    corrs<-calcCorrelation(ixdata, iydata)
    colnames(corrs)<-c("Sample", paste0(ixname, "_", iyname, c("_R", "_pvalue")))
    if(is.null(finaltable)){
      finaltable<-corrs
    }else{
      finaltable<-cbind(finaltable, corrs[,c(2,3),drop=F])
    }
    
    points<-getPoint(ixdata, iydata, ixname, iyname)
    if(is.null(finalpoint)){
      finalpoint<-points
    }else{
      finalpoint<-rbind(finalpoint, points)
    }
  }
}

write.csv(finaltable, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_r_pvalue.csv", row.names=F)

#common<-rownames(ixdata)[rownames(ixdata) %in% rownames(iydata)]
#commonx<-ixdata[common, 1]
#commony<-iydata[common, 1]

#cd<-data.frame(name=common, x=commonx, y=commony)
#sumxy<-(commonx + commony) > 0
#cd<-cd[sumxy, ]
#cd<-cd[order(cd$name),]
#cor.test(cd$x, cd$y, method="spearman")

finaltableR = finaltable[,grepl("_R", colnames(finaltable))]
rownames(finaltableR)<-finaltable$Sample
colnames(finaltableR) = gsub("_R", "", colnames(finaltableR))
colnames(finaltableR) = gsub("_", " vs ", colnames(finaltableR))
finaltableRGroup<-gsub("_WT_.*", " WT", rownames(finaltableR))
finaltableR.mean <- aggregate(x = finaltableR, 
                              by = list(finaltableRGroup), 
                              FUN = mean)
colnames(finaltableR.mean)[1]<-"Group"
write.csv(finaltableR.mean, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_group_R_mean_new.csv", row.names=F)

finaltableR.stdev <- aggregate(x = finaltableR, 
                               by = list(finaltableRGroup), 
                               FUN = sd)
colnames(finaltableR.stdev)[1]<-"Group"
write.csv(finaltableR.stdev, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_group_R_stdev.csv", row.names=F)

meltdata<-melt(finaltable, id="Sample")
meltdata<-meltdata[grepl("_R", meltdata$variable),]
meltdata$variable<-gsub("_R", "", meltdata$variable)
meltdata$variable<-gsub("_", " vs ", meltdata$variable)
meltdata$Group<-gsub("_WT.*", " WT", meltdata$Sample)

pdf("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_group_boxplot.pdf", height=7, width=7)
g<-ggplot(meltdata, aes(x=variable, y=value)) + 
  geom_boxplot(width=0.5) + 
  ylab("Spearman correlation") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5),
        strip.background =element_rect(fill=NA)) + 
  facet_grid(Group~.)
print(g)
dev.off()

pdf("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_group_violin.pdf", height=7, width=7)
g<-ggplot(meltdata, aes(x=variable, y=value)) + 
  geom_violin(width=0.7, color="red", fill="red") + 
  geom_jitter(width=0.2, size=1) +
  ylab("Spearman correlation") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5),
        strip.background =element_rect(fill=NA)) + 
  facet_grid(Group~.)
print(g)
dev.off()

# finalpoint$X <- log2(finalpoint$X+1)
# finalpoint$Y <- log2(finalpoint$Y+1)
# finalpoint<-finalpoint[,c("Sample", "X", "Y", "Xname", "Yname")]
# pcopy<-finalpoint[,c("Sample", "Y", "X", "Yname", "Xname")]
# colnames(pcopy)<-c("Sample", "X", "Y", "Xname", "Yname")
# pdata<-rbind(finalpoint, pcopy)
# pdata$Xname<-as.factor(as.character(pdata$Xname))
# pdata$Yname<-as.factor(as.character(pdata$Yname))
# pdata$Group<-gsub("_WT.*", " WT", pdata$Sample)
# 
# pdf("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_scatter.pdf", height=7, width=7)
# g<-ggplot(pdata, aes(x=X, y=Y)) + 
#   geom_point() + 
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) + 
#   facet_grid(Xname~Yname)
# print(g)
# dev.off()
# 
# pdata$Sample<-as.character(pdata$Sample)
# sample<-unique(pdata$Sample)[1]
# for(sample in unique(pdata$Sample)){
#   outfile = paste0("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/pipeline_correlation/KCV_3018_77_78_79.miRNA.correlation_scatter_", sample, ".pdf")
#   cdata<-subset(pdata, Sample==sample)
#   pdf(outfile, height=7, width=7)
#   g<-ggplot(cdata, aes(x=X, y=Y)) + 
#     geom_point() + 
#     theme_bw() +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) + 
#     facet_grid(Xname~Yname)
#   print(g)
#   dev.off()
# }
# 

getMeltData<-function(x, name){
  x$Feature<-rownames(x)
  result<-melt(x, id="Feature")
  result$SampleFeature<-paste0(result$variable, ":", result$Feature)
  result$Pipeline<-name
  return(result[c("SampleFeature", "Pipeline", "value")])
}

mChimira<-getMeltData(chimira, "Chimira")
mExcerpt<-getMeltData(excerptr, "exceRpt")
mMirge<-getMeltData(mirge, "miRge")
mOasis<-getMeltData(oasis, "OASIS")
mTiger<-getMeltData(tiger, "TIGER")
mFinal<-rbind(mChimira, mExcerpt)
mFinal<-rbind(mFinal, mMirge)
mFinal<-rbind(mFinal, mOasis)
mFinal<-rbind(mFinal, mTiger)
ftable<-dcast(mFinal, SampleFeature ~ Pipeline)
rownames(ftable)<-ftable$SampleFeature
ftable<-ftable[,c(2:ncol(ftable))]
ftable<-log2(ftable+1)

corTestWithoutZero <- function(x, y, method="spearman") {
  sumxy<-!is.na(x) & !is.na(y) & (x != 0) & (y != 0)
  ccx<-x[sumxy]
  ccy<-y[sumxy]
  cor.test(ccx, ccy, method=method)
}

panel.cor <- function(x, y, digits = 2, prefix = "", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  test <- corTestWithoutZero(x, y)
  txt <- format(c(test$estimate, 0.123456789), digits = digits)[1]
  txt <- paste("R = ", txt, sep = "")
  Signif <- symnum(test$p.value, 
                   corr = FALSE, 
                   na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("p<0.0001", "p<0.001", "p<0.01", " ", ""))
  text(0.5, 0.4, txt, cex = 2)
  text(0.5, 0.7, Signif, cex = 2)
}

#t1<-ftable[grepl("APOB_WT_01", rownames(ftable)),]
#t1<-t1[!is.na(t1$Chimira) & !is.na(t1$exceRpt) & (t1$Chimira + t1$exceRpt > 0),]
#cor.test(t1$Chimira, t1$exceRpt, method="spearman")

for (group in c("APOB_WT", "HDL_WT", "Liver_WT")){
  #pdf(paste0("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_scatter_", group, ".pdf"), height=7, width=7)
  png(paste0("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_scatter_", group, ".png"), height=4000, width=4000, res=600)
  pairs(ftable[grepl(group, rownames(ftable)),], gap = 0, upper.panel = panel.cor, main="log2(miRNA count + 1)", pch=19, cex=0.5)
  dev.off()
}
#pdf(paste0("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_scatter_all.pdf"), height=7, width=7)
png(paste0("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRNA.correlation_scatter_all.png"), height=4000, width=4000, res=600)
pairs(ftable, gap = 0, upper.panel = panel.cor, main="log2(miRNA count + 1)", pch=19, cex=0.5)
dev.off()
