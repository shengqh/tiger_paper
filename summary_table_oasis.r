host_file<-"/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.mapped.count"
final_unmapped_file<-"/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/final_unmapped/final_unmapped_reads_summary/result/KCV_3018_77_78_79.count"
nonhost_file<-"/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/reads_in_tasks/result/KCV_3018_77_78_79.TaskReads.csv"

host<-read.table(host_file, sep="\t", header=T, row.names=1)
unannotated<-read.table(final_unmapped_file, sep="\t", header=T, row.names=1)
nonhost<-read.csv(nonhost_file, row.names=1)

combined<-rbind(host[c("TotalReads", "MappedReads", "FeatureReads", "miRNA", "tRNA", "snRNA", "rRNA", "lincRNA", "misc_RNA", "snoRNA"),],
                unannotated["FeatureReads",],
                nonhost[c("bacteria_group1", "bacteria_group2", "fungus_group4", "miRBase", "tRNA", "rRNA"),])

annotated<-apply(combined, 2, function(x){
  round((x["MappedReads"] + x["FeatureReads1"]) / x["TotalReads"],2)
})

final<-rbind(combined, annotated)
rownames(final)<-c("Total reads", 
                   "Host reads", 
                   "Host smallRNA reads", 
                   "miRNA total reads", 
                   "tDR total reads", 
                   "snDR total reads", 
                   "rDR total reads",
                   "lncDR total reads",
                   "misc_RNA total reads",
                   "snoDR total reads",
                   "Nonhost reads",
                   "HMB (Bac) reads",
                   "Environment (Bac) reads",
                   "Fungi reads",
                   "Non-host miRNA reads",
                   "Non-host tDR reads",
                   "Non-host rDR reads",
                   "% assigned")

miRNA<-read.table("/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.miRNA.count", sep="\t", header=T, row.names=1)
miRNA<-miRNA[,colnames(final)]
miRNArpm<-t(t(miRNA * 1000000) / unlist(final[1,]))
miRNArpmm<-t(t(miRNA * 1000000) / unlist(final[4,]))

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

tRNA<-read.table("/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.tRNA.count", sep="\t", header=T, row.names=1)
tRNA<-tRNA[,colnames(final)]
tRNArpmt<-t(t(tRNA * 1000000) / unlist(final[5,]))
tRNArpmt100<-apply(tRNArpmt, 2, function(x){
  sum(x > 100)
})

snRNA<-read.table("/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.snRNA.count", sep="\t", header=T, row.names=1)
snRNA<-snRNA[,colnames(final)]
snRNArpms<-t(t(snRNA * 1000000) / unlist(final[6,]))
snRNArpms100<-apply(snRNArpms, 2, function(x){
  sum(x > 100)
})

rRNA<-read.table("/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV_3018_77_78_79.rRNA.count", sep="\t", header=T, row.names=1)
rRNA<-rRNA[,colnames(final)]
rRNArpmr<-t(t(rRNA * 1000000) / unlist(final[7,]))
rRNArpmr100<-apply(rRNArpmr, 2, function(x){
  sum(x > 100)
})

finalTable<-rbind(final[1:4,],
            "miRNA > 10RPM" = miRNArpm10,
            "miRNA > 100RPMM" = miRNArpmm100,
            "miR-22-3p RPM" = miR22_3p_rpm,
            "miR-22-3p RPMM" = miR22_3p_rpmm,
            "miR-92a-3p RPM" = miR92a_3p_rpm,
            "miR-92a-3p RPMM" = miR92a_3p_rpmm,
            final[5,],
            "tDR > 100 RPMtDR" = tRNArpmt100,
            final[6,],
            "snDR > 100 RPMsnDR" = snRNArpms100,
            final[7,],
            "rDR > 100 RPMrDR" = rRNArpmr100,
            final[8:nrow(final),])

write.csv(finalTable, file="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/KCV_3018_77_78_79.summery.csv")
