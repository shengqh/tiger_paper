#Oasis
data <- read.table("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/Oasis/59e909b321b68/data/summary/images/data.tsv", sep = '\t', header = TRUE)
colnames(data) <- colnames(data[2:17])
data <- t(data[1:(ncol(data)-1)])

#summary files per sample
sample <- dir("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/Oasis/59e909b321b68/data/counts/species",
    pattern = "smallrna.table",full.names=TRUE)

summary <- data.frame(matrix(NA,nrow=5, ncol=length(sample)))
for (i in 1:length(sample)){
  summary[,i]<-fread(sample[i],select = 2, data.table=FALSE)
  }

rownames(summary) <- c("miRNA", "piRNA", "snoRNA", "snRNA", "rRNA")
colnames(summary) <- gsub(".smallrna.table", "", basename(sample))

#miRNA files
sample_miRNA <- dir("T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/other_pipelines/Oasis/59e909b321b68/data/counts/species",
            pattern = "mirnaCounts.txt",full.names=TRUE)

summary_miRNA <- list()
for (i in 1:length(sample_miRNA)){
  summary_miRNA[[i]]<-read.table(sample_miRNA[i])
}

summary_miRNA_final <- Reduce(function(...) merge(..., all=TRUE, by= "V1"), summary_miRNA)
summary_miRNA_final <- data.frame(summary_miRNA_final[,2:ncol(summary_miRNA_final)],row.names = summary_miRNA_final[,1])
colnames(summary_miRNA_final) <- gsub(".smallrna.table", "", basename(sample))
write.csv(summary_miRNA_final,"T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.oasis.miRNA.csv")

miRNArpm <- t(t(summary_miRNA_final* 1000000) / unlist(data[1,]))
miRNArpmm <- t(t(summary_miRNA_final * 1000000) / unlist(data[9,]))

miRNArpm10 <- apply(miRNArpm, 2, function(x){
  sum(x > 10)
})
miRNArpmm100 <- apply(miRNArpmm, 2, function(x){
  sum(x > 100)
})

miR22_3p_rpm<-miRNArpm["mmu-miR-22-3p",]
miR22_3p_rpmm<-miRNArpmm["mmu-miR-22-3p",]
miR92a_3p_rpm<-miRNArpm["mmu-miR-92a-3p",]
miR92a_3p_rpmm<-miRNArpmm["mmu-miR-92a-3p",]

host_reads<- (data["Genome",]+ colSums(summary))

annotated <- (data["Genome",]+ colSums(summary))/ data["Initial",]

finalTable <- rbind(data[1,],
                    host_reads,
                    summary[1,],
                    miRNArpm10,
                    miRNArpmm100,
                    miR22_3p_rpm,
                    miR22_3p_rpmm,
                    miR92a_3p_rpm,
                    miR92a_3p_rpmm,
                    summary[c(4,5,3,2),],
                    annotated)
rownames(finalTable)<-c("Total Reads","Host Reads", "miRNA total reads",
                        "miRNA > 10RPM", "miRNA > 100RPMM", "miR-22-3p RPM", "miR-22-3p RPMM",
                        "miR-92a-3p RPM", "miR-92a-3p RPMM", "snDR total reads", "rDR total reads",
                        "snoDR total reads", "piRNA total reads", "% assigned")

write.csv(finalTable, file="T:/Shared/Labs/Vickers Lab/Tiger/projects/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/forpaper/KCV_3018_77_78_79.oasis.summary.csv")

