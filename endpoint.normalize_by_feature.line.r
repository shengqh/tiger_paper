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

groupFileList=parSampleFile1
visLayoutFileList=parSampleFile2
positionFile = parFile1
totalCountFile<-parFile2

#load Rcpp package first because of the error with reshape2 package
library(Rcpp)
library(grid)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

position<-read.delim(positionFile, header=T,as.is=T)
position<-subset(position, SampleRank < 11)

colnames(position)[colnames(position)=="Sample"]<-"File"
colnames(position)[colnames(position)=="RelativeSampleEndpoint"]<-"Position"
colnames(position)[colnames(position)=="EndpointCount"]<-"PositionCount"
colnames(position)[colnames(position)=="Group"]<-"Category"

samples<-unique(position$File)

groupInfo<-read.delim(groupFileList, header=F, as.is=T)
groups<-unique(groupInfo$V2)
groupSize<-table(groupInfo[,2])

countTable<-read.delim(totalCountFile,as.is=T,header=T,row.names=1,check.names=FALSE)

position$FeatureLabel<-paste0(position$Feature,"(",round(position$TotalCount,0),")")

allGroupPosition<-NULL
for(group in groups){
  groupSamples<-groupInfo$V1[groupInfo$V2==group]
  groupPositions<-position[position$File %in% groupSamples,]
  groupPositions$Group<-group
  allGroupPosition<-rbind(allGroupPosition, groupPositions)
}

categories<-unique(allGroupPosition$Category)

allCategoryPosition<-allGroupPosition
groupData<-allCategoryPosition[,c("Group", "Category", "File", "Feature")]
groupData<-groupData[!duplicated(groupData),]
groupTable<-table(groupData[,c(1,2)])

allPositionByGroup<-aggregate(x = allCategoryPosition, by = list(allCategoryPosition$Group, allCategoryPosition$Category, allCategoryPosition$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
allPositionByGroup$Position<-allPositionByGroup$Group.3
gsize<-apply(allPositionByGroup, 1, function(x){
  return (groupTable[x["Group"], x["Category"]])
})

allPositionByGroup$GroupPositionCountFraction<-as.vector(allPositionByGroup$Percentage/gsize)

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

allPositionMeanMiRNA<-subset(allPositionMean, Category == "miRNA")

colors<-c("grey", brewer.pal(length(categories)-1, "Set1"))
names(colors)<-c("miRNA", categories[categories != "miRNA"])

png(paste0(outFile,".allPositionLine.normalizedByAnchor.png"),width=width,height=height,res=300)
m<-ggplot() +
  geom_area(data=allPositionMeanMiRNA, aes(x = Position, y=GroupPositionCountFractionMean), fill="grey") +
  geom_line(data=allPositionMean, aes(x = Position,y=GroupPositionCountFractionMean, color=Category), size=1.5) +
  facet_grid(Row_Group~Col_Group,space = "free",scale="free") +
  theme_bw()+
  ylab("Average read coverage percentage (per feature)")+
  xlim(-10, 10) +
  theme(text = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1.3, "lines"), 
        legend.position="right") +
  scale_color_manual(values=colors) +
  guides(fill= guide_legend(ncol=2,keywidth=1, keyheight=1.5))
print(m)
dev.off()
