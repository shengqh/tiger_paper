summary_file<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/excerptr/postProcessedResults_v4.6.3/PPR-2017-11-27-16%3A00%3A55_exceRpt_readMappingSummary.txt"
summary<-read.table(summary_file, sep="\t", header=T, row.names=1)
summary<-summary[sort(rownames(summary)),]
summary<-t(summary)

reads_file<-"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/excerptr/postProcessedResults_v4.6.3/PPR-2017-11-27-16%3A00%3A55_exceRpt_biotypeCounts.txt"
reads<-read.table(reads_file, sep="\t", header=T, row.names=1)
reads<-reads[,colnames(summary)]

final<-rbind(summary[c("input","genome"),],
             reads[c("miRNA", "tRNA","snRNA"),],
             summary["rRNA",,drop=F],
             reads[c("lincRNA","misc_RNA", "snoRNA", "piRNA"),])

rownames(final)<-c("Total Reads", 
                   "Host Reads", 
                   "miRNA total reads", 
                   "tDR total reads", 
                   "snDR total reads", 
                   "rDR total reads", 
                   "lncDR total reads", 
                   "misc_RNA total reads", 
                   "snoDR total reads", 
                   "piRNA total reads")

miRNA<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/excerptr/postProcessedResults_v4.6.3/PPR-2017-11-27-16%3A00%3A55_exceRpt_miRNA_ReadCounts.txt", sep="\t", header=T, row.names=1)
miRNA<-miRNA[,colnames(final)]
write.csv(miRNA, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.excerptr.miRNA.csv")

miRNArpm<-t(t(miRNA * 1000000) / unlist(final[1,]))
miRNArpmm<-t(t(miRNA * 1000000) / unlist(final[3,]))

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

tRNA<-read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/excerptr/postProcessedResults_v4.6.3/PPR-2017-11-27-16%3A00%3A55_exceRpt_tRNA_ReadCounts.txt", sep="\t", header=T, row.names=1)
tRNA<-tRNA[,colnames(final)]
tRNArpmt<-t(t(tRNA * 1000000) / unlist(final[4,]))
tRNArpmt100<-apply(tRNArpmt, 2, function(x){
  sum(x > 100)
})

finalTable<-rbind(final[1:3,],
            "miRNA > 10RPM" = miRNArpm10,
            "miRNA > 100RPMM" = miRNArpmm100,
            "miR-22-3p RPM" = miR22_3p_rpm,
            "miR-22-3p RPMM" = miR22_3p_rpmm,
            "miR-92a-3p RPM" = miR92a_3p_rpm,
            "miR-92a-3p RPMM" = miR92a_3p_rpmm,
            final[4,,drop=F],
            "tDR > 100 RPMtDR" = tRNArpmt100,
            final[5:nrow(final),,drop=F])
colnames(finalTable)<-gsub("sample_","",colnames(finalTable))
colnames(finalTable)<-gsub("_fastq","",colnames(finalTable))

write.csv(finalTable, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.excerptr.summery.csv")
