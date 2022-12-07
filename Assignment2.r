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

# SET WORKING DIRECTORY to intial folder

# 1. Import Data via tximport

setwd("Data\\NB_SY5Y_MYCN_data_salmon_counts")
samples <- read.table("samplelist.txt", header = FALSE)
colnames(samples) <- c("condition","run")
# set as the reference the untreated samples
samples$condition <- factor(samples$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))
samples$condition <- relevel(samples$condition, "TOT_0h")
samples

files <- file.path("NB_SY5Y_MYCN_data_salmon_counts", samples$run, "quant.sf")
names(files) <- samples$run
all(file.exists(files))
files

# If 1st run, skip it and run the next step
# !!!! Skip the next step after the 1st run !!!!
# Next run is time consuming
setwd("..")
tx2gene <- read.csv("tx2gene.csv", header = TRUE)

# If it is 1st run, RUN it
# otherwise SKIP IT!!
# prepare txdb, create map from tx to gene
setwd("..")
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.101.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
#write.csv(tx2gene,"tx2gene.csv", row.names = FALSE)

#Note, you can query the txdb by keys of a keytype and get columns back. 
#To see the available keytypes: keytypes(txdb) and coluns: columns(txdb)
#keytypes(txdb)
#columns(txdb)

# 2. Import Data via tximport
# tx import salmon files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
head(txi.salmon$counts)

# 3. Construct DESeqDataSet object via DESeqDataSetFromTximport

# run the data through DESeq2
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
zero_rows <- rowSums(counts(dds)) == 0
sum(zero_rows)
one_rows <- rowSums(counts(dds)) == 1
sum(one_rows)
sum(rowSums(counts(dds)) < 5)
keep <- rowSums(counts(dds)) > 1
sum(keep)
dds <- dds[keep,]
nrow(dds)

# boxplots
# Boxplot of Cooks distances shows that the samples have similar forms
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
par(mar=c(1,1,1,1))

# 4.2 The variance stabilizing transformation and the rlog
# VST transformation
# For visualization purposes is it usually beneficial to somehow normalize or standardize 
# the data. Here we use the variance shrinkage algorithm 
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
# based on VST values
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
rownames(samplePoisDistMatrix) <- paste( dds$run, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# 4.4 PCA Plot

par(mfrow = c(2,2))
plotPCA(vsd, intgroup = c("run"))
plotPCA(vsd, intgroup = c("condition"))
plotPCA(vsd, intgroup = c("run","condition"))
plotPCA(vsd, intgroup = c("condition","run"))
par(mfrow = c(1,1))


# Copy visualizing steps from the provided DESeq2.r 
## Principal components plot
plotPCA(vsd, intgroup = c("condition"))

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
table(res$padj < 0.05)
# 5.2 Filter the results for significance

# Filter the results by false discovery rate
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
table(res.05$padj < 0.1)
#head(res.05)
summary(res.05)

# Filter the results by log2 fold change threshold
resLFC1SE <- results(dds, lfcThreshold=1)
summary(resLFC1SE)
table(resLFC1SE$padj < 0.1)
table(res.05$padj < 0.05)

# Filter the results false  discovery and by log2 fold change (LFC)
LFC1SE <- results(dds, lfcThreshold=1, alpha=0.05)
summary(LFC1SE)

# 5.3 Other comparisons - TO BE IMPLEMENTED
# Complicated contrasts
# Refer to https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html#pathway-analysis
other_res_0_1<- results(dds, contrast = c("condition", "TOT_1h", "TOT_0h"))
other_res_0_1
summary(other_res_0_1)
other_res_1_4<- results(dds, contrast = c("condition", "TOT_4h", "TOT_1h"))
other_res_1_4
summary(other_res_1_4)
other_res_4_24<- results(dds, contrast = c("condition", "TOT_24h", "TOT_4h"))
other_res_4_24
summary(other_res_4_24)
other_res_0_4<- results(dds, contrast = c("condition", "TOT_4h", "TOT_0h"))
other_res_0_4
summary(other_res_0_4)
other_res_<- results(dds, contrast = c("condition", "TOT_24h", "TOT_0h"))
other_res_
summary(other_res_)

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
head(resSig[ order(resSig$log2FoldChange), ],10)

# the significant genes with the strongest up-regulation:
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ],10)

# 6. Plotting results
# 6.1 Counts plot
# Normalized counts for a single gene over treatment group.
# visualize the counts for a particular gene
#Have a look at the tip gene:
# the gene with the smallest padj value
topGene <- rownames(res)[which.min(res$padj)]
topGene
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

resultsNames(dds)
shrink_res <- lfcShrink(dds, coef="condition_TOT_24h_vs_TOT_0h", type="apeglm")
plotMA(shrink_res, ylim = c(-2, 2))

# If not shrink the noisy log2 fold changes
res.noshr <- results(dds, name="condition_TOT_24h_vs_TOT_0h")
plotMA(res.noshr, ylim = c(-5, 5))

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
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("run" , "condition")])
pheatmap(mat, annotation_col = anno)

# Save the results to be ready for Pathway Analysis
# order by log2FoldChange
resOrderedByFC <- res[order(res$log2FoldChange),]
resDF <- as.data.frame(resOrderedByFC)
write.csv(resDF, file = "deseq_results.csv")

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
# We sort the results by to pvalue
resOrdered <- res[order(res$pvalue),]
head(resOrdered,10)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
resOrderedDF <- apply(resOrderedDF,2,as.character)
write.csv(resOrderedDF, file = "ordered_results.csv")

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
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
x1<-pairwise_termsim(gse)
emapplot(x1, showCategory = 10)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list)
ridgeplot(gse) + labs(x = "enrichment distribution")



# ----------------------
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]


df2 = df[df$X %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
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
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
x2 <- pairwise_termsim(kk2)
emapplot(x2)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")

# Produce KEGG pathways
# Test: hsa05202, hsa00043, hsa04722, hsa04974
# Produce the native KEGG plot (PNG)
# Neuroblastoma - Related Pathway - hsa04722
# !!!! MYCN related pathway - 	hsa05202 
# Transcriptional misregulation in cancer: hsa05202
# 	Pathways in cancer: hsa05200
setwd("..\\Pathways")
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05200", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05200", species = kegg_organism, kegg.native = F)

knitr::include_graphics("hsa04722.pathview.png")

# --------------------------
## GSEA using msigdbr
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
# barplot(res_GSEA, showCategory=2)
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

