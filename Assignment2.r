#install.packages("PoiClaClu")
#install.packages("glmpca")
#install.packages("ggbeeswarm")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("pheatmap")
#BiocManager::install("apeglm")
#BiocManager::install("genefilter")
#BiocManager::install("ReportingTools")
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("msigdbr")
#BiocManager::install("ggnewscale")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("OmnipathR")

# Load libraries
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")  # for annotations
#install.packages("ggridges")

library(GenomicFeatures)
library(tximport)
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db") # for annotations
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("glmpca")
library("ggbeeswarm")
library("apeglm")
library("genefilter")
library("ReportingTools")
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)
library("pathview")
library(EnsDb.Hsapiens.v86)


organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

## ----------------------------------
##      
## To run the R script create the following folder tree:
##  
## - Intitial folder
##    - Assignment.R
##    - Data     # contains data files.
##        - Homo_sapiens.GRCh38.101.gtf.gz
##        - NB_SY5Y_MYCN_data_salmon_counts   # contains salon files
##            - ERR571495     
##            - ERR571497     
##            ...
##            - sampleslist.txt
##    - Pathways    # it will contain the produced pathways
##
## ----------------------------------  

## ----------------------------------  
##
##    SET WORKING DIRECTORY to initial folder
##
## ----------------------------------  

# 1. Import Data via tximport

# Load the samples table
setwd("Data\\NB_SY5Y_MYCN_data_salmon_counts")
samples <- read.table("samplelist.txt", header = FALSE)
colnames(samples) <- c("condition","run")
# set as the reference the untreated samples
samples$condition <- factor(samples$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))
samples$condition <- relevel(samples$condition, "TOT_0h")
samples

# Prepare a list with all salmon files
files <- file.path("NB_SY5Y_MYCN_data_salmon_counts", samples$run, "quant.sf")
names(files) <- samples$run
all(file.exists(files))
files

# Create map from tx to gene
# This is time consuming process thus
# If it is 1st run, skip step A and run the next B to create tx2gene.csv
# If it is NOT 1st run, just run step A to read tx2gene.csv
# !!!! Skip the next step after the 1st run !!!!
# Step A
setwd("..")
tx2gene <- read.csv("tx2gene.csv", header = TRUE)

# If it is 1st run, RUN it
# otherwise SKIP IT!!
# Step B
setwd("..")
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.101.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
write.csv(tx2gene,"tx2gene.csv", row.names = FALSE)

#Note, you can query the txdb by keys of a keytype and get columns back. 
#To see the available keytypes: keytypes(txdb) and coluns: columns(txdb)
#keytypes(txdb)
#columns(txdb)

# 2. Import Data via tximport
# tx import salmon files
# BE CAREFUL!!!! Current Directory should be: Data
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
head(txi.salmon$counts)

# 3. Construct DESeqDataSet object via DESeqDataSetFromTximport
# DESeq2 object
dds <- DESeqDataSetFromTximport(txi.salmon, colData =  samples, ~ condition)
dds$condition <- factor(dds$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))
dds
# TOT_0Η should be the 1st level at condition design formula
# DEseq will use the first level as reference. 
# It makes more sense if this is the untreated control
dds$condition # check the levels of the design formula

# 4.0 Filtering & Visualization
# 4.1. Pre-filtering the dataset
# 
nrow(dds)
# find out how many rows have 0 counts
zero_rows <- rowSums(counts(dds)) == 0
sum(zero_rows)
# find out how many rows have 1 counts
one_rows <- rowSums(counts(dds)) == 1
sum(one_rows)
# find out how many rows have <5 counts
sum(rowSums(counts(dds)) < 5)
# delete the rows that have count=0 and count=1 
keep <- rowSums(counts(dds)) > 1
sum(keep)
dds <- dds[keep,]
nrow(dds)


# 4.2 The variance stabilizing transformation and the rlog
# VST transformation
# For visualization purposes is it usually beneficial to somehow normalize or standardize 
# the data. 
# variance shrinkage algorithm 
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
# rlog transformation
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# plot the first sample against the second, first simply using 
# the log2 function (after adding 1, to avoid taking the log of zero), 
# and then using the VST and rlog-transformed values. 
ddsE <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(ddsE, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


# 4.3 Sample distances
# Euclidean distance based on VST values
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$run, vsd$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Calculate sample distances by using Poisson Distance
# that takes the original count matrix 
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vsd$run, vsd$condition, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# 4.4 PCA Plot

plotPCA(vsd, intgroup = c("run"))
plotPCA(vsd, intgroup = c("condition")) 
plotPCA(vsd, intgroup = c("run","condition"))
plotPCA(vsd, intgroup = c("condition","run")) # used for reporting!

# Another way of doing it that offer more plotting flexibility:
# is extracting the data and using theggplot2 package
pcaData <- plotPCA(vsd, intgroup = c( "run", "condition"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = run, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# 4.5 PCA plot using Generalized PCA

gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$run <- dds$run
gpca.dat$condition <- dds$condition

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = run, shape = condition)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# 4.6 MDS Plot
# based on vst
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = run, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# based on Poisson Distance 
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = run, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

# 5. Differential expression analysis
# 5.1 Running the differential expression pipeline
dds <- DESeq(dds)
res <- results(dds)
res
mcols(res, use.names = TRUE) # lists the meaning of the columns
#head(res)
summary(res)
table(res$padj < 0.1)
# 5.2 Filter the results for significance
# run results for adjusted p-value < 0.05
res.05 <- results(dds, alpha = 0.05)
summary(res.05)
table(res.05$padj < 0.05)
#head(res.05)


# Filter the results by log2 fold change threshold
resLFC1SE <- results(dds, lfcThreshold=1)
summary(resLFC1SE)
table(resLFC1SE$padj < 0.1)
# Filter the results false discovery and by log2 fold change (LFC)
resLFC1SE_a <- results(dds, alpha = 0.05, lfcThreshold=1)
summary(resLFC1SE_a)
table(resLFC1SE_a$padj < 0.05)

# 5.3 Other comparisons

other_res_0_1<- results(dds, contrast = c("condition", "TOT_1h", "TOT_0h"))
summary(other_res_0_1)

other_res_1_4<- results(dds, contrast = c("condition", "TOT_4h", "TOT_1h"))
summary(other_res_1_4)

other_res_4_24<- results(dds, contrast = c("condition", "TOT_24h", "TOT_4h"))
summary(other_res_4_24)

other_res_0_4<- results(dds, contrast = c("condition", "TOT_4h", "TOT_0h"))
summary(other_res_0_4)

# 5.4 Multiple testing !!!!!!!!
# compute the number of genes with a p value < 0.05 
# for which the test succeeded, Ho is rejected
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
# there are 1785 genes with a p value < 0.05 among 24548 genes 
# for which the test succeeded, Ho is rejected

# If Ho is true, then  by the definition of the p value, we expect 
# up to 5% of the genes to have a pvalue<0.05
sum(!is.na(res$pvalue)) * 0.05
# this amounts to 1227 genes
# if we consider that the list of genes with a pvalue<0.05 as 
# differentially expresd, this list should contain up to 
# 1227 / 1785 = 68.7% False Positive
sum(!is.na(res$pvalue)) * 0.05 / sum(res$pvalue < 0.05, na.rm=TRUE)
# Consider that 10% false positives are acceptable, then all genes 
# with an adjusted pvalue <10% = 0.1 are considered significant.
sum(res$padj < 0.1, na.rm=TRUE)

# the significant genes with the strongest down-regulation
resSig <- subset(res, padj < 0.1)
dim(resSig)
down<-head(resSig[ order(resSig$log2FoldChange), ],10)
down

# the significant genes with the strongest up-regulation:
up<-head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ],10)
up

# 6. Plotting results

# 6.1 Counts plot
# Normalized counts for a single gene over treatment group.
# visualize the counts for a particular gene
# Have a look at the tip gene:
# the gene with the smallest padj value
topGene <- rownames(res)[which.min(res$padj)]
topGene   # it is the MYCN gene - ENSG00000134323

plotCounts(dds, gene = topGene, intgroup=c("condition"))
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("run","condition"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = run, y = count, color = condition)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)


# # Normalized counts for a single gene over treatment group.
# DESeq test actually takes into account the cell line effect, 
# so this figure more closely depicts the difference being tested
ggplot(geneCounts, aes(x = run, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# 6.2 MA-plot
# Provides overview for the distribution of the estimated coefficients in the model, 
# e.g. the comparisons of interest, across all genes.
# An MA-plot of changes induced by treatment.
# Comparison of interest: TOT_24h vs TOT_0h
# Shrink the noisy log2 fold changes
resultsNames(dds)
shrink_res <- lfcShrink(dds, coef="condition_TOT_24h_vs_TOT_0h", type="apeglm")
plotMA(shrink_res, ylim = c(-2, 2))

# If not shrink the noisy log2 fold changes
res.noshr <- results(dds, name="condition_TOT_24h_vs_TOT_0h")
plotMA(res.noshr, ylim = c(-5, 5))

# label individual points on the MA-plot as well
# topGene = MYCN
# "ENSG00000135679" = MDM2
# "ENSG00000141510" = PT53
plotMA(shrink_res, ylim = c(-4,4))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
mdm2gene <- "ENSG00000135679"   # MDM2
with(res[mdm2gene, ], {
  points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
  text(baseMean, log2FoldChange, mdm2gene, pos=4, col="red")
})
pt53gene <- "ENSG00000141510"   # PT53
with(res[pt53gene, ], {
  points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
  text(baseMean, log2FoldChange, pt53gene, pos=4, col="green")
})

# Histogram of the p values
# genes with very small counts excluded, otherwise generate spikes appear
hist(shrink_res$pvalue[shrink_res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# Fit curve to gene-wise dispersion estimates
plotDispEsts(dds)


# 6.3 Clustering !!!
# see this web page for a nice illustration of functionality
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

# clustering is only relevant for genes that actually carry a signal
# so cluster a subset of the most highly variable genes
# select the 20 genes with the highest variance across samples
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

# center each genes’ values across samples
# clustering by run and condition
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("run", "condition")])
pheatmap(mat, annotation_col = anno)

# clustering by condition
anno1 <- as.data.frame(colData(vsd)[, c("condition")])
row.names(anno1) <- colnames(mat)
colnames(anno1)[1] = "condition"
pheatmap(mat, annotation_col = anno1,cluster_cols=F)

# heatmap of gene under investigation
myGenes <- c(which(rownames(res)=="ENSG00000134323"),
             which(rownames(res)=="ENSG00000141510"),
             which(rownames(res)=="ENSG00000135679")
             )
mat1  <- assay(vsd)[ myGenes, ]
mat1  <- mat1 - rowMeans(mat1)
anno2 <- as.data.frame(colData(vsd)[, c("run", "condition")])
pheatmap(mat1, annotation_col = anno2)

# Save the results to be ready for Pathway Analysis
# order by log2FoldChange
resOrderedByFC <- res[order(res$log2FoldChange),]
resDF <- as.data.frame(resOrderedByFC)
write.csv(resDF, file = "deseq_results.csv")

# Although the DESeq results do not report any outlier
# Just for verification, check if outliers exist
# Boxplot of Cook's distances shows that the samples have similar forms
# for scaling purposes Cook's distances transformed by log10
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, ylab="log10 of Cook's distances")
par(mar=c(1,1,1,1))

# 7.0 Annotating and exporting results
# Annotation is important to make the data and result interpretable. 
# People usually know genes by their gene symbol (abbreviated gene name).
# the result table so far only contains the Ensembl gene IDs
# We can get the gene symbols (and many other types of identifiers and information)
# using the AnnotationDbi package and corresponding database for our organism, 
# here human: ogr.Hs.eg.db. 

#columns(org.Hs.eg.db)
#columns(EnsDb.Hsapiens.v86)
ens.str <- rownames(res)
ens.str
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name <- mapIds(org.Hs.eg.db, 
                   keys=ens.str,
                   column="GENENAME",
                   keytype="ENSEMBL", 
                   multiVals="first")

# Check the entries having NA at ENTREZID
# Most of NAs are pseudogenes !!!! 
sum(is.na(res$entrez))
mapIds(EnsDb.Hsapiens.v86, keys="ENSG00000288393", keytype="GENEID", column="SYMBOL")
mapIds(EnsDb.Hsapiens.v86, keys="ENSG00000288393", keytype="ENTREZID", column="SYMBOL")
mapIds(EnsDb.Hsapiens.v86, keys="ENSG00000288393", keytype="GENENAME", column="SYMBOL")
#head(res,20)


## 7.1 Exporting the results
# It is nice to have the most significant results on the top. 
# Sort the results by to pvalue
resOrdered <- res[order(res$pvalue),]
head(resOrdered,10)
# save the results into ordered_results.csv
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
resOrderedDF <- apply(resOrderedDF,2,as.character)
write.csv(resOrderedDF, file = "ordered_results.csv")

# get those with padj < 0.1 and sort by log2FoldChange
resSig <- subset(res, padj < 0.1)
down<-head(resSig[ order(resSig$log2FoldChange), ],20)
# save the results into down_regulated.csv
downDF <- as.data.frame(down)
write.csv(downDF, file = "down_regulated.csv")

# sort by log2FoldChange in descending order
up<-head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ],20)
up
# save the results into up_regulated.csv
upDF <- as.data.frame(up)
write.csv(upDF, file = "up_regulated.csv")

# 8.0 PATHWAY ANALYSIS

# reading in data from deseq2 analysis
df = read.csv("deseq_results.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X
any(is.na(df$log2FoldChange))

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# ----------------------
# Gene set enrichment - Gene Ontology
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", #BP: Biological Processes, #ALL: BP, CC, MF
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# Dot plot visualizes enriched terms. 
# It depicts the gene ratio as bar height and color
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
# Enrichment map organizes enriched terms into a network with edges 
# connecting overlapping gene sets. 
x1<-pairwise_termsim(gse)
emapplot(x1, showCategory = 20)
# cnetplot depicts the linkages of genes and biological concepts as a network
cnetplot(gse, categorySize="pvalue", foldChange=gene_list)
# The ridgeplot helps in interpreting up/down-regulated pathways.
ridgeplot(gse) + labs(x = "enrichment distribution")


# ----------------------
# Gene set enrichment - KEGG
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df2 = df[df$X %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
# create the gene list
kegg_gene_list <- df2$log2FoldChange
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
# Dot plot visualizes enriched terms. 
# It depicts the gene ratio as bar height and color
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
# Enrichment map organizes enriched terms into a network with edges 
# connecting overlapping gene sets. 
x2 <- pairwise_termsim(kk2)
emapplot(x2)
# cnetplot depicts the linkages of genes and biological concepts as a network
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list,showCategory = 3)
# The ridgeplot helps in interpreting up/down-regulated pathways.
ridgeplot(kk2) + labs(x = "enrichment distribution")

# Produce KEGG pathways
# Test: 
# Neuroblastoma - Related Pathway - hsa04722
# !!!! MYCN related pathway - 	hsa05202 
# Transcriptional misregulation in cancer: hsa05202
# Pathways in cancer: hsa05200
# Produce the native and regulated KEGG plot (PNG)
setwd("..\\Pathways")
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05202", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05202", species = kegg_organism, kegg.native = F)

knitr::include_graphics("hsa04722.pathview.png")

# --------------------------
## GSEA using msigdbr
# create gen list
gene <- names(kegg_gene_list)[abs(kegg_gene_list) > 2]
geneList <- kegg_gene_list
head(gene)

msigdbr_species()

# Get all human gene sets:
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

# Alternatively, get a specific collection, e.g. Hallmarks (H)
# Note: dplyr selects the columns we need for enrichemnt analysis: 
# gene set names (gs_name) and gene identifiers (entrez_gene). 
m_df_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
head(m_df_H_t2g)

# View what collections there are
collections <- msigdbr_collections()
collections

# Select (subset) the Hallmark collection
m_df_H <- m_df[m_df$gs_cat=="H",]
head(m_df_H)

# GSEA using msigdbr
#For enrichemnt analysis we need a dataframe with two columns
# 1st column gene set name (gs_name), 2nd column gene identifier (entrez_gene)
#
#msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Enrichment by hypergeometric test implemented in the enricher function
res_enricher <- enricher(gene = gene, TERM2GENE = msigdbr_t2g)

# Gene set enrichment analysis is implemented in the GSEA funcion
res_GSEA <- GSEA(geneList, TERM2GENE = msigdbr_t2g)

# The GSEA result is a special data type: gseaResult. 
# We can convert it into a dataFrame:
res_GSEA_df <- as.data.frame(res_GSEA)

# Visualization
dotplot(res_GSEA, showCategory=22)

# Network and heatmap visualization
# This is sometimes usefull, but not here for the hallmark genes set
# because there are too may genes in the genesets with too much
# overlap
#
# convert gene ID to Symbol
edox <- setReadable(res_GSEA, 'org.Hs.eg.db', 'ENTREZID')
# visualize 1. network, then heatmap:
cnetplot(edox, foldChange=geneList)
heatplot(edox, foldChange=geneList)

# -----------------------
## Omnipath allows us to query information from many database resources.
## Here we use it to get the target genes of the transcription factors 
## MYCN

library(OmnipathR)
library(igraph)

## We check some of the different interaction databases
get_interaction_resources()

## We query and store the interactions into a dataframe
interactions <- import_all_interactions(
  #resources = c('HPRD', 'BioGRID'),
  organism = 9606
)

## We select the most confident interactions for a given TF and we print
## the interactions to check the way it regulates its different targets
interactions_A_MYCN  <- dplyr::filter(
  interactions, #dorothea_level==c("A","B"),
  source_genesymbol == "MYCN" & type == "transcriptional")
print_interactions(interactions_A_MYCN)

## than for MYCN
dim(interactions_A_MYCN)

## It is often useful to interpret data in terms of the underlying network

##  Convert interaction data into a graph (network)
I_g <- interaction_graph(interactions = interactions)


## Find and print shortest paths on the directed network between proteins
## of interest (requires igraph package):
print_path_es(igraph::shortest_paths(I_g,from = "MYCN",to = "MDM2",
                                     output = 'epath')$epath[[1]],I_g)

## Find and print all shortest paths between proteins of interest:
print_path_vs(all_shortest_paths(I_g,from = "MYCN",
                                 to = "MDM2")$res,I_g)
