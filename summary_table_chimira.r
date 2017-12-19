tiger<-read.csv(file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.tiger.summery.csv", header=T, row.names=1)

miRNA<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/chimira/4d94b105-ff45-634f-b562-a56c01b00213.plain_counts_per_file.counts", sep="\t", header=T, row.names=1)
colnames(miRNA)<-gsub("X","",colnames(miRNA))
colnames(miRNA)<-gsub("\\.","",colnames(miRNA))
samples<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/chimira/samples.txt", header=T, row.names=1)
colnames(miRNA)<-samples[colnames(miRNA), "Name"]
colnames(miRNA)<-gsub(".gz","", colnames(miRNA))
miRNA<-miRNA[,colnames(tiger)]

write.csv(miRNA, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.chimira.miRNA.csv")

totalMiRNA<-apply(miRNA, 2, sum)

summary<-rbind("Total Reads" = tiger["Total reads",],
               "miRNA total reads" = totalMiRNA)


miRNArpm<-t(t(miRNA * 1000000) / unlist(summary[1,]))
miRNArpmm<-t(t(miRNA * 1000000) / unlist(summary[2,]))

miRNArpm10<-apply(miRNArpm, 2, function(x){
  sum(x > 10)
})
miRNArpmm100<-apply(miRNArpmm, 2, function(x){
  sum(x > 100)
})
miR22_3p_rpm<-miRNArpm["mmu-mir-22-3p",]
miR22_3p_rpmm<-miRNArpmm["mmu-mir-22-3p",]
miR92a_3p_rpm<-miRNArpm["mmu-mir-92a-3p",]
miR92a_3p_rpmm<-miRNArpmm["mmu-mir-92a-3p",]

annotated<-apply(summary, 2, function(x){
  round(x["miRNA total reads"] / x["Total reads"],2)
})

finalTable<-round(rbind(summary[1:2,],
            "miRNA > 10RPM" = miRNArpm10,
            "miRNA > 100RPMM" = miRNArpmm100,
            "miR-22-3p RPM" = miR22_3p_rpm,
            "miR-22-3p RPMM" = miR22_3p_rpmm,
            "miR-92a-3p RPM" = miR92a_3p_rpm,
            "miR-92a-3p RPMM" = miR92a_3p_rpmm,
            "% assigned" = annotated),2)

write.csv(finalTable, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.chimira.summary.csv")


