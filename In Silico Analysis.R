#install.packages(readxl)
#install.packages("tidyr")
#install.packages("corrplot")
#install.packages("RColorBrewer")
#install.packages("edgeR")
#install.packages("PCATools")
#install.packages("openxlsx")
#install.packages("gplots")
#install.packages("Hmisc")

#Select the Directory where you can find the RData from the script: ReadCounts_TablePrep and Supplementary Material 1

##################################################################################
############# Exploratory analysis of transcriptome profiles #####################
##################################################################################

load("SamplesExpression.RData")

group <- factor(Metadata[,"sampletype"], levels=unique(Metadata[,"sampletype"]))

# Creating a DGEList object for use in edgeR
library(edgeR)

y <- DGEList(counts=ReadCounts,group=group)

# Filtering low expressed genes

keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

# Calculate normalizing factors

y <- calcNormFactors(y)

#Creating a CPMS matrix (normalized values for PCA and clustering)

logCPMs<- edgeR::cpm(y, offset = y$offset, log = TRUE)


#########
#Figure 2 C - Principal component analysis biplots of human skin fibroblasts before and upon TNF-α stimulus
########

library(PCAtools)

#PCA

metadata <- data.frame(group)
rownames(metadata) <- colnames(logCPMs)

#Run PCA

pca.res <- pca(logCPMs, metadata=metadata)

#Plot variance explained by each component

# Plot PC1 versus PC2

biplot(pca.res, colby="group",title = "Biplot PC1 vs PC2") # Biplot with colors by sample type
biplot(pca.res, x="PC1", y="PC2",lab="",colby="group",legendPosition = 'right', labSize = 200,legendLabSize = 20,pointSize =6 ) # Biplot with PC1 and PC2



##########################################################################################
####################### Differential expression analysis #################################
##########################################################################################

library(edgeR)


load("SamplesExpression.RData")


group <- factor(Metadata[,"sampletype"], levels=unique(Metadata[,"sampletype"]))


y <- DGEList(counts=ReadCounts,group=group)


# Filtering low expressed genes
keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

# Calculate normalizing factors
y <- calcNormFactors(y)

# Estimating the dispersion and Fit model

###Creation of the Design and Contrast Matrix

design_matrix <- model.matrix(~0+group); colnames(design_matrix) <- gsub("group", "", colnames(design_matrix))

contrasts <- makeContrasts("WTStim_vs_WTNStim"=WTStim-WTNStim,levels=design_matrix)

y <- estimateDisp(y, design_matrix)

fit <- glmQLFit(y,design_matrix)

# Identify DEGs - WTStim_vs_WTNSTim

#Only cut-off criteria was FDR ≤ 0.05 

qlf <- glmQLFTest(fit, contrast = contrasts)

DEG <-topTags(qlf, n=nrow(qlf$table),p.value=0.05)

#-1   1 
#29 130

# Calculate the Log2(Fold-Change) of all genes regardless of the significant change

LogFC <- topTags(qlf, n=nrow(qlf$table),p.value=1)

#############################################################################################
###################Writing Output Table with the Differently Expression Genes################
#############################################################################################



################################ Primary Biological Questions################################
###########Separated DEGs based on upregulation and downregulation upon stimulus##############
#############################################################################################

#Follow the prior Pipeline

library(openxlsx)

separate_positive_negativeDEG <- function(df){
  # Function that creates variables with DEG upregulated and downregulated
  DEG_positive <- vector()
  DEG_negative <- vector()
  for (i in 1:length(rownames(df))){
    ifelse(sign(df[i,1]) == +1,DEG_positive <- c(DEG_positive, rownames(df)[i]), DEG_negative <- c(DEG_negative,rownames(df)[i]) )
  }
  return (list(pos = as.data.frame(DEG_positive),neg = as.data.frame(DEG_negative)))
}



# write output table for WT VS PMM2-CDG
xlsx_logCPM <- cbind("GeneName"=rownames(logCPMs), logCPMs)
LogFoldChange <- cbind("GeneName"=rownames(LogFC$table), LogFC$table)
Discriminatory_list <- list("Norm_Expression" = xlsx_logCPM ,"WTStim_vs_WTNStim"=LogFoldChange,separetadedDEGs["pos"],separetadedDEGs["neg"])
names(Discriminatory_list) <- c("Norm_Expression","Log(Fold-Change)","DEGs_upregulated","DEGs_downregulated")
write.xlsx(Discriminatory_list, file = "Expression.xlsx")



###########################################################################################
################################# Figure 2 D ##############################################
#################################Volcano Plot##############################################
###########################################################################################

# Produce plots
library(gplots)


plotData <- as.data.frame(topTags(qlf, n=nrow(qlf$table)))

# Assign colors to upregulated and downregulated genes
plotData$colour <- ifelse(plotData$logFC > 0 & plotData$FDR < 0.05,"springgreen4",
                         ifelse(plotData$logFC < 0 & plotData$FDR < 0.05, "gray20","tomato3"))




#Select DEGs that should be Highlighted in the Plot
DEGs_Higlights <- rownames(DEG$table[DEG$table[,1]> 5 | DEG$table[,1] < -2.2,])


# Create the volcano plot

ggplot(plotData, aes(x = logFC, y = -log10(FDR), color = colour)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("tomato3",  "gray20", "springgreen4"),
                     labels = c("Downregulated", "Non Significant", "Upregulated")) +
  xlim(c(-6, 15)) +
  ylim(c(0, 3.25)) +
  geom_hline(yintercept = -log10(max(DEG$table$FDR)), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  labs(title = "Volcano Plot", x = "Log2(Fold Change)", y = "-log10(FDR)") +
  guides(fill = 'transparent',colour = guide_legend(override.aes = list(size=6),label.theme = element_text(size = 15))) +
  geom_text_repel(data = plotData[DEGs_Higlights , ], 
                  aes(x = logFC, y = -log10(FDR), label = DEGs_Higlights), hjust = 0, fontface= "bold", vjust = 0, color = "grey27")+
  theme( panel.background = element_rect(fill = "white"),
         panel.spacing = unit(1.5, "lines"),
         legend.position = "top",
         legend.key = element_blank(),
         legend.box.background = element_blank(),
         panel.grid.major = element_line(colour = "gray87"), panel.grid.minor = element_line(colour = "gray87"),
  axis.line = element_line(colour = "black"), legend.background = element_rect(fill='transparent')
)


dev.off()



###########################################################################################
################################# Figure 2 E ##############################################
##############################Correlation Plot#############################################
###########################################################################################

library(corrplot)
library(openxlsx)
library(Hmisc)

LogFC <- topTags(qlf, n=nrow(qlf$table),p.value=1)


#Download From Cytosig the LogFC of the genes from every experiment that was conducted in Immune cells and Fibroblasts upon TNF-α stimulus
#Cytosig database: https://cytosig.ccr.cancer.gov/
TNAStim_cell <- as.data.frame(read_excel("TNFA_stim.xlsx", col_names = TRUE, sheet =2))

rownames(TNAStim_cell) <- TNAStim_cell[,1]

TNAStim_cell <- TNAStim_cell[,-1]

#Match Rows Between Our Fibroblasts, Immune Cells and other Fibroblasts to have a complete dataframe of LogFC for each experiment

#Remove NA values
TNAStim_cell <- TNAStim_cell[!rowSums(TNAStim_cell== "NA"),]

id <- match(rownames(TNAStim_cell),rownames(LogFC$table))

SK_Stim <- LogFC$table[id[!is.na(id)],]

TNAStim_cell <- TNAStim_cell[rownames(SK_Stim),]

skin_Fibroblast <- SK_Stim[,1]



#Creation of a final Matrix with all cells that were stimulated with fibroblasts and our skin fibroblats

AllConditions <- data.frame(TNAStim_cell,skin_Fibroblast)

AllConditions <- as.data.frame(sapply(AllConditions,as.numeric))

rownames(AllConditions) <-rownames(SK_Stim)



#Calculation of the pearson coefficient correlation score and their significance

cor_matrix <- rcorr(as.matrix(AllConditions), type = "pearson")$r

p_values <- rcorr(as.matrix(AllConditions), type = "pearson")$P

adj_p <- p.adjust(p_values, method = "fdr")


#Matrix of the pearson coefficient correlation scores

skin_Fibroblast <- as.matrix(cor_matrix[,"skin_Fibroblast"][-10])
colnames(skin_Fibroblast) <- "Skin Fibroblast (This Study)"
rownames(skin_Fibroblast) <- colnames(AllConditions)[-10]
rownames(skin_Fibroblast) <- gsub("[_.]", " ", rownames(skin_Fibroblast))
rownames(skin_Fibroblast)[1] <- "CD4+ Memory T Cells"


#Matrix of the adjusted p-value

adj_p<- as.matrix(adj_p[80:88])
colnames(adj_p) <- "Skin Fibroblast (This Study)"
rownames(adj_p) <- colnames(Final)[-10]
rownames(adj_p) <- gsub("[_.]", " ", rownames(skin_Fibroblast))
rownames(adj_p)[1] <- "CD4+ Memory T Cells"


#Color Gradient to measure the correlation plot

scalebluered <- colorRampPalette(c("#2166AC","white","#D6604D"))(200)

corrplot(skin_Fibroblast, 
         method="color", 
         type="full", 
         cl.pos = "n",
         tl.srt=65,
         order="original", 
         p.mat = adj_p, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.75,
         diag=TRUE, 
         tl.col="black", 
         tl.cex =1,
         col=colorRampPalette(c("#2166AC","white","#D6604D"))(200),
         addgrid.col="black",
)

#Legend for the Gradient that It is used as a form of measurement of the Correlation Plot
colorlegend(
  scalebluered,
  c(seq(-1,1,0.25)),
  xlim=c(1.65,2.25),
  ylim=c(0.5,9.5),
  align="l",
  vertical=TRUE,
  addlabels=TRUE,
  cex=1
)


#############################################################################################
################################## Figure 3 E ###############################################
########################Heatmap Molecular Interactions#######################################
#############################################################################################


library(readxl)
library(tidyr)
library(corrplot)
library(RColorBrewer)

#Upload Supplementary Material 1 where there is the different interactions betweens fibroblasts DEGs and immune cells receptors and secreted modulators
Interactions <- as.data.frame(read_excel("Supplementary Material 1.xlsx", col_names = TRUE, sheet = "Fibroblast_immune interactions"))


#Create a vector with all the DEGs identified in an interaction with an immune cell
selDEG <- as.vector(Interactions[,"DEG_Fibroblasts"])
split_elements <-as.vector(unlist(strsplit(selDEG, ",")))
updated_string <- gsub(" ", "", split_elements)
uniqueDEG <- unique(updated_string)

          
#Selection of the positive Differently Expressed Genes within the interaction list
Sel_positive_Deg<- DEG$table[uniqueDEG,"logFC"] >0


#Log2(Fold-Change) of Positive DEGs
expgenes <- as.matrix(as.numeric(DEG$table[as.vector(unlist(uniqueDEG)),"logFC"]))
rownames(expgenes) <- as.vector(unlist(uniqueDEG))
expgenes <- as.matrix(expgenes[Sel_positive_Deg,])
colnames(expgenes)<- " "

#FDR of the Positive DEGs
adj_p <- as.matrix(as.numeric(DEG$table[as.vector(unlist(uniqueDEG)),"FDR"]))
rownames(adj_p) <- as.vector(unlist(uniqueDEG))
adj_p <- as.matrix(adj_p[Sel_positive_Deg,]) 
colnames(adj_p)<- " "


# Subset the gene expression data for the genes of interest
gene_data <- rownames(expgenes)

# Create a vector of presence/absence data for each gene

Immune_Int <- Interactions[,c("DEG_Fibroblasts","Monocytes", "B cells", "T cells", "DCs", "Eosinophils", "Macrophages", "NK cells", "Neutrophils")]

#Removal of string with values separated with a comma
Immune_Int_sep<-as.data.frame(Immune_Int %>%               
  separate_rows(DEG_Fibroblasts, sep=","))

#Replace NA with Absent
Immune_Int_sep[is.na(as.data.frame(Immune_Int_sep))] <- "Absent"

#Creation of an empty Matrix where we will classify the interactions between Fibroblasts and Immune Cells

Int <- matrix(nrow=length(rownames(expgenes)), ncol=8)
rownames(Int) <- gene_data

for(i in 1:length(gene_data)){
  for(n in 1:length(rownames(Immune_Int_sep))){
    if(gene_data[i] == Immune_Int_sep[n,"DEG_Fibroblasts"]){
      for(l in 2:length(colnames(Immune_Int_sep))){
        if(Immune_Int_sep[n,l] == 1){
          col <- l - 1
          Int[i,col] <- "Present"
        }
      }
    }
    
  }
}

#Replace NA with Absent
Int[is.na(as.data.frame(Int))] <- "Absent"

colnames(Int) <-c("Monocytes", "B cells", "T cells", "DCs", "Eosinophils", "Macrophages", "NK cells", "Neutrophils")


# creation of a data frame with the gene names, log fold change values, and interaction information
df <- data.frame(gene = rownames(expgenes),
                 logfc = expgenes[,1],
                 interaction = Int)

#Heatmap with Predicted Interactions
corrplot(expgenes, 
         method="color", 
         type="full", 
         cl.pos = "n",
         is.corr = FALSE,
         #cl.length = length(rownames(skin_Fibroblast)),
         #cl.cex = 1,
         #cl.ratio =1,
         #cl.align.text ="r",
         p.mat = adj_p, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         tl.srt=65,
         order="original",
         pch.cex = 0.75,
         diag=FALSE, 
         tl.col="black", 
         tl.cex =0.75,
         col=brewer.pal(8,"Greens"),
         addgrid.col="black",
)


# Addition of a color gradient legend based on the Log2FC of each gene upon stimulus

colorlegend(
  brewer.pal(8,"Greens"),
  c(seq(0,round(max(expgenes),1), by =2.5)),
  xlim=c(-3.2,-1.3),
  ylim=c(18.1,22.6),
  align="l",
  vertical=TRUE,
  addlabels=TRUE,
  cex=0.75
)


#Addition of Predicted Interactions based on the downloaded excel file

#create a vector of colors to use for the points

col <- c("tomato","goldenrod1",brewer.pal(7,"YlGnBu")[2:7])

# add coloured squares based on the interaction with each immune cell

for (i in 1:nrow(Int)) {
  for (j in 1:ncol(Int)) {
    x <- j + ncol(expgenes)+0.5
    y <- nrow(df) +1 - i
    if(Int[i,j] == "Present"){
    symbols(x, y, squares = 0.75, inches = FALSE, add = TRUE, fg = col[j], bg = col[j], lwd = 1, ylim = c(0, ncol(expgenes)+ncol(Int)), xaxt = "n", yaxt = "n", ann = FALSE, pch = 15)
    }
  }
}





# add text for different scales

text(-2.6,23.2, "Log2(FC)", cex = 0.85, adj =0.5, family = "Arial")

x  <- (ncol(expgenes) + 0.5 + ncol(Int))/2+1
y <- nrow(df)+1

text(x,y,"Predicated Interactions", cex = 0.95, adj = 0.5)

for(j in 1:ncol(Int)){
  x<-  j + ncol(expgenes) + 0.5 
  y <- -3.5
  text(x,y,colnames(Int)[j], cex = 0.85, adj = c(0, 0.5),srt =90)
  
}


dev.off()



###############################################################################
############################## Figure 4A ######################################
#############################Heatmap Plot######################################
#######H#######   Heatmap of immune-related gene expression ###################
###############################################################################

#You will need to obtain the LogCPMs variable which is at the top of the script

library(RColorBrewer); library(gplots)

plotCol_exp <- brewer.pal(9, "Greens")
plotColSamples <- c(rep("grey50",3),rep("grey70", 2))
par(mar = rep(2, 4))

#Creation of a Vector with the Expression of the genes of interest
Genes_interest <- c("IL1B", "IL6", "IL15", "CXCL1", "CXCL5", "CXCL8", "CCL2", "CCL5", "NFKB1","NFKB2", "NFKBIA")

#Used HeatMap
heatmap.2(logCPMs[Genes_interest,], Rowv=FALSE,  scale="row", dendrogram = "none", cexCol = 1, cexRow = 0.57, trace="none", density.info="none", ColSideColors = plotColSamples, col=plotCol_exp,main="Immunological Response",margins = c(0.01,25))
