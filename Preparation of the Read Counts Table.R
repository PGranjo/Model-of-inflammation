###########################################
#######Create Read Counts table ###########
###########################################

#Download RNA-seq Files from each sample in different conditions
#A working directory solely with this files was selected
ReadCount<- c()


for (i in list.files()){
  read_count_single_file <- as.data.frame(read.csv(file = i, header = FALSE, sep = "\t"))
  colnames(read_count_single_file )<-c("GeneID","Read Counts")
  rownames(read_count_single_file )<-read_count_single_file$GeneID
  read_count_single_file$GeneID<-NULL
  ReadCount<- as.matrix(cbind(ReadCount,read_count_single_file$"Read Counts"))
}

colnames(ReadCount) <- c("1WTStim","2WTStim","3WTStim",
                               "4WTNStim","5WTNStim")
rownames(ReadCount) <- rownames(read_count_single_file)


########################################
#### Change Gene Names Annotation ######
########################################

#install.packages("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://uswest.ensembl.org")

attributes = listAttributes(ensembl)
#attributes[1:50,]

geneIDs <- rownames(ReadCount)

allInfo <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters = "ensembl_gene_id", 
                 values = geneIDs,
                 mart = ensembl)

geneNames <- allInfo[match(geneIDs, allInfo[,2]),1]


#substitute the gene codes by gene names 
ReadCounts <- ReadCount[match(allInfo$ensembl_gene_id,rownames(ReadCount)),]

genenames1<-allInfo$hgnc_symbol
rownames(ReadCounts) <- genenames1
ReadCounts <- as.data.frame(ReadCounts)


#######################################################
############## Prepare Metadata File ##################
#######################################################
sample <- colnames(ReadCounts)
sampletype <- gsub("^[0-9]+","",sample)
Metadata <- data.frame(sample,sampletype)


##########################################################################
########Save Key Variables as RData for further analysis##################
##########################################################################

save(ReadCounts,Metadata, file="SamplesExpression.RData")



