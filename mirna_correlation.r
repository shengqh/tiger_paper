library(reshape2)
library(ggplot2)

chimira<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.chimira.miRNA.csv", header=T, row.names=1)
rownames(chimira)<-gsub("miR","mir",rownames(chimira))

excerptr<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.excerptr.miRNA.csv", header=T, row.names=1)
colnames(excerptr)<-gsub("sample_","",colnames(excerptr))
colnames(excerptr)<-gsub("_fastq","",colnames(excerptr))
multi<-grepl("\\|", rownames(excerptr))
excerptr<-excerptr[!multi,]
mnames<-rownames(excerptr)[multi]
enames<-unlist(lapply(mnames, function(x){
  strsplit(x, "\\|")
}))
enames<-c("mmu-mir-28c", "mmu-mir-378d", enames)
excerptr<-excerptr[!(rownames(excerptr) %in% enames),]
rownames(excerptr)<-gsub("miR","mir",rownames(excerptr))

mirge<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRge.miRNA.csv", header=T, row.names=1)
rownames(mirge)<-gsub("/.+","",rownames(mirge))
rownames(mirge)<-gsub("miR","mir",rownames(mirge))

tiger<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.miRNA.count", sep="\t", header=T, row.names=1)
tiger<-tiger[,colnames(mirge)]
rownames(tiger)<-gsub(";.+","",rownames(tiger))
rownames(chimira)<-gsub("miR","mir",rownames(chimira))

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
    data.frame("Sample"=sample, "Spearman"=cor(ccx, ccy, method="spearman"))
  }))
}

datalistnames<-list("Chimira", "excerptr", "miRge", "TIGER")
datalists<-list(chimira, excerptr, mirge, tiger)

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
    colnames(corrs)<-c("Sample", paste0(ixname, "_", iyname))
    if(is.null(finaltable)){
      finaltable<-corrs
    }else{
      finaltable<-cbind(finaltable, corrs[,2,drop=F])
    }
  }
}

meltdata<-melt(finaltable, id="Sample")
meltdata$X<-gsub("_.*", "", meltdata$variable)
meltdata$Y<-gsub(".*_", "", meltdata$variable)
meltdata<-meltdata[,c("X", "Y", "value")]
mcopy<-meltdata[,c("Y", "X", "value")]
colnames(mcopy)<-c("X", "Y", "value")
mdata<-rbind(meltdata, mcopy)

ggplot(mdata, aes(y=value)) + 
  geom_violin(aes(x=1)) + 
  geom_point(aes(x=1)) +
  xlab("") +
  ylab("Spearman correlation") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + facet_grid(X~Y)
