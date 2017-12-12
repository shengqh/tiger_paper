summary_file<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/miRge/report.txt"
summary<-read.table(summary_file, sep="\t", header=T)
summary<-summary[summary$File.name.s. != "",]
summary<-summary[sort(rownames(summary)),]
summary$File.name.s. = gsub("/gpfs23/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/preprocessing/cutadapt/result/", "", summary$File.name.s.)
summary$File.name.s. = gsub("_clipped.fastq.gz", "", summary$File.name.s.)
rownames(summary)<-summary$File.name.s.
summary$All.miRNA.Reads...Filtered.miRNA.Reads = gsub("^.+\\s","",summary$All.miRNA.Reads...Filtered.miRNA.Reads)
summary<-t(summary)
summary<-gsub(" ","", summary)
summary<-summary[c(2,4,7,8),sort(colnames(summary))]
class(summary)<-"numeric"
rownames(summary)<-c("Total Reads", 
                   "miRNA total reads", 
                   "other ncRNA total reads", 
                   "mRNA total reads")

miRNA<-read.csv("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/miRge/miR.Counts.csv", header=T, row.names=1)
colnames(miRNA)<-gsub("_clipped","",colnames(miRNA))
miRNA<-miRNA[c(2:nrow(miRNA)),colnames(summary)]

write.csv(miRNA, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRge.miRNA.csv")


miRNArpm<-t(t(miRNA * 1000000) / unlist(summary[1,]))
miRNArpmm<-t(t(miRNA * 1000000) / unlist(summary[2,]))

miRNArpm10<-apply(miRNArpm, 2, function(x){
  sum(x > 10)
})
miRNArpmm100<-apply(miRNArpmm, 2, function(x){
  sum(x > 100)
})
miR22_3p_rpm<-miRNArpm["mmu-miR-22-3p",]
miR22_3p_rpmm<-miRNArpmm["mmu-miR-22-3p",]
miR92a_3p_rpm<-miRNArpm["mmu-miR-92a-3p",]
miR92a_3p_rpmm<-miRNArpmm["mmu-miR-92a-3p",]

finalTable<-rbind(summary[1:2,],
            "miRNA > 10RPM" = miRNArpm10,
            "miRNA > 100RPMM" = miRNArpmm100,
            "miR-22-3p RPM" = miR22_3p_rpm,
            "miR-22-3p RPMM" = miR22_3p_rpmm,
            "miR-92a-3p RPM" = miR92a_3p_rpm,
            "miR-92a-3p RPMM" = miR92a_3p_rpmm,
            summary[3:nrow(summary),,drop=F])

write.csv(finalTable, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.miRge.summary.csv")
