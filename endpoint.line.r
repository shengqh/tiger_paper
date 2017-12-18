rm(list=ls()) 
setwd("/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/host_endpoint_vis")
outFile='KCV_3018_77_78_79.EndpointVis.Host'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='KCV_3018_77_78_79.endpoint.txt'
parFile2='/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.mapped.count'
parFile3=''

options(bitmapType='cairo')

# TODO: Add comment
# 
# Author: zhaos
###############################################################################

maxFeature=100

groupFileList=parSampleFile1
visLayoutFileList=parSampleFile2
positionFile = parFile1
totalCountFile<-parFile2

countName = "TotalReads"

#load Rcpp package first because of the error with reshape2 package
library(Rcpp)
library(grid)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)

position<-read.delim(positionFile, header=T,as.is=T)
colnames(position)[colnames(position)=="Sample"]<-"File"
colnames(position)[colnames(position)=="RelativeSampleEndpoint"]<-"Position"
colnames(position)[colnames(position)=="EndpointCount"]<-"PositionCount"
colnames(position)[colnames(position)=="Group"]<-"Category"

samples<-unique(position$File)

groupInfo<-read.delim(groupFileList, header=F, as.is=T)
groups<-unique(groupInfo$V2)
groupSize<-table(groupInfo[,2])

totalCount<-read.delim(totalCountFile,as.is=T,header=T,row.names=1,check.names=FALSE)
totalCount<-unlist(totalCount[countName,])

position$FeatureLabel<-paste0(position$Feature,"(",round(position$TotalCount,0),")")

allCategoryPosition<-NULL
for(group in groups){
  groupSamples<-groupInfo$V1[groupInfo$V2==group]
  groupPositions<-position[position$File %in% groupSamples,]
  groupPositions$Group<-group
  allCategoryPosition<-rbind(allCategoryPosition, groupPositions)
}

allPositionByGroup<-aggregate(x = allCategoryPosition, by = list(allCategoryPosition$Feature, allCategoryPosition$Group, allCategoryPosition$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
allPositionByGroup$GroupPositionCountFraction<-as.vector(allPositionByGroup$PositionCountFraction/groupSize[allPositionByGroup$Group])
allPositionByGroup$Position<-allPositionByGroup$Group.3

if (visLayoutFileList!="") {
  visLayout<-read.delim(visLayoutFileList,as.is=T,header=F)
  visLayout<-sapply(split(visLayout[,1],visLayout[,2]),function(x) x)
  row.names(visLayout)<-visLayout[,"Groups"]
  visLayout<-data.frame(visLayout[,-which(colnames(visLayout)=="Groups")])
  visLayout$Col_Group<-factor(visLayout$Col_Group,levels=unique(visLayout$Col_Group))
  visLayout$Row_Group<-factor(visLayout$Row_Group,levels=unique(visLayout$Row_Group))
  
  allPositionByGroup<-data.frame(allPositionByGroup,visLayout[allPositionByGroup[,"Group"],])
} else {
  allPositionByGroup$Group<-factor(allPositionByGroup$Group,levels=groups)
}

if (visLayoutFileList!="") {
  rowGroupCount=length(unique(allPositionByGroup$Row_Group))
  colGroupCount=length(unique(allPositionByGroup$Col_Group))
  height=max(rowGroupCount*1500,2000)
  width=max(colGroupCount*1500,2000)
} else {
  height=3000
  groupCount=length(unique(allPositionByGroup$Group))
  width=max(groupCount*2000, height)
}

allPositionSlim<-allPositionByGroup[,c("Category", "Group", "Position", "GroupPositionCountFraction", "Col_Group", "Row_Group")]
allPositionMean<- aggregate(x = allPositionSlim, by = list(allPositionSlim$Category, allPositionSlim$Group, allPositionSlim$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
allPositionMean$Position<-allPositionMean$Group.3
allPositionMean$GroupPositionCountFractionMean<-as.numeric(allPositionMean$GroupPositionCountFraction / groupSize[allPositionMean$Group])
allPositionMean<-subset(allPositionMean, Position >= -10)
allPositionMean<-subset(allPositionMean, Position <= 10)

colors<-primary.colors(length(unique(allPositionMean$Category)))

png(paste0(outFile,".allPositionLine.png"),width=width,height=height,res=300)
m<-ggplot(allPositionMean, aes(x = Position,y=GroupPositionCountFractionMean, color=Category)) +
  geom_line(size=1.5) +
  theme_bw()+
  ylab("average read coverage (per million total reads)")+
  xlim(-10, 10) +
  theme(text = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1.5, "lines"), 
        legend.position="right") +
  scale_color_manual(values=colors) +
  guides(fill= guide_legend(ncol=2,keywidth=1, keyheight=1.5))

if (visLayoutFileList!="") {
  m<-m+facet_grid(Row_Group~Col_Group,space = "free",scale="free")
} else {
  m<-m+facet_grid(.~Group,space = "free",scale="free")
}
print(m)
dev.off()
