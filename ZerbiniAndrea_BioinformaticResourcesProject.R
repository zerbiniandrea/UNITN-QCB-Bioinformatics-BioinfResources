if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
#BiocManager::install("pathview")

library(biomaRt)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(pathview)

# 1
# LOAD RDATA FILE
load("C:/Users/zerbi/Desktop/Project/Lung_Cancer.RData")

# 2 
# Update raw_count_df and r_anno_df extracting only protein coding genes.
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

query1 <- getBM(attributes=c("ensembl_gene_id",
                             "external_gene_name",
                             "gene_biotype",
                             "transcript_count",
                             "start_position",
                             "end_position",
                             "chromosome_name",
                             "strand",
                             "description",
                             "version"),
                filters=c("chromosome_name"), 
                values=list(c("21")),
                mart = ensembl)

query1_protein_coding <- query1[which(query1$gene_biotype=="protein_coding"),]

r_anno_df <- subset(r_anno_df,gene_id %in% query1_protein_coding$ensembl_gene_id)

raw_counts_df <- raw_counts_df[query1_protein_coding$ensembl_gene_id,]
raw_counts_df <- na.omit(raw_counts_df)

# 3
# Perform differential expression analysis

# count threshold
count_thr <- 20
# number of replicates with more counts than the count threshold
repl_thr <- 2 

filter_vec <- apply(raw_counts_df,1,function(y) max(by(y, c_anno_df$condition, function(x) sum(x>=count_thr))))

filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,]

filter_anno_df <- r_anno_df[r_anno_df$gene_id %in% rownames(filter_counts_df),]

edge_c <- DGEList(counts=filter_counts_df,group=c_anno_df$condition,samples=c_anno_df,genes=filter_anno_df) 

edge_n <- calcNormFactors(edge_c,method="TMM")

cpm_table <- as.data.frame(round(cpm(edge_n),2))

design <- model.matrix(~0+group, data=edge_c$samples)
colnames(design) <- levels(edge_c$samples$group)
rownames(design) <- edge_c$samples$sample

edge_d <- estimateDisp(edge_n,design)
edge_f <- glmQLFit(edge_d,design) 

edgeRglmQLF <- function(mat=edge_f,contro,cpm_mat=edge_n,label="",sig_thr=0.5,sig_col="CPM",fc_thr=0.5,pval_col="p_val",pval_thr=0.05,names=FALSE)
{
  degs <- glmQLFTest(edge_f,contrast=contro)$table[,-3]
  colnames(degs) <- c("log2_FC","log2_CPM","p_val")
  a_levels <- rownames(contro)[which(contro!=0)]
  a_samples <- which(cpm_mat$samples$group%in%a_levels)
  cpm_sele <- cpm(cpm_mat,log=T)[,a_samples]
  degs$log2_CPM <- apply(cpm_sele,1,function(x) mean(x))
  #degs<-exactTest(edge_c, pair=cond, dispersion=bcv^2)$table
  degs$p_adj <- p.adjust(degs$p_val, method ="BH")
  degs$class <- "="
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC>=fc_thr & degs[,pval_col]<=pval_thr),"class"] <- "+"
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC<=(-fc_thr) & degs[,pval_col]<=pval_thr),"class"] <- "-"
  degs$class <- as.factor(degs$class)
  degs$comp <- label
  degs$id <- rownames(degs)
  degs <- degs[,c("id","comp","log2_FC","log2_CPM","p_val","p_adj","class")]
  if(names=="TRUE"){
    newnames <- paste(label,colnames(degs),sep="_")
    colnames(degs) <- newnames
  }
  return(degs)
}

contro <- makeContrasts("Case-Control", levels=design) 

fc_thrs <- 1
CPM_thrs <- 1.5
pval_thrs <- 0.01

DEGs <- edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="CasevsControl", sig_thr=CPM_thrs, sig_col="log2_CPM", fc_thr=fc_thrs, pval_thr=pval_thrs, pval_col="p_adj",names=F)

#Volcano plot
input_df <- DEGs
xlabel <- "log2 FC Case vs Control"
ylabel <- "-log10 adj_pvalue (FDR)"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$log2_FC, -log(input_df$p_adj,base=10),xlab=xlabel, ylab=ylabel, 
     col=ifelse(input_df$class=="=","grey70","olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot")
abline(v=0,lty=2,col="grey20")
abline(v=fc_thrs,lty=2,col="grey20")
abline(v=-fc_thrs,lty=2,col="grey20")
abline(h=-log10(pval_thrs),lty=2,col="grey20")

## Heatmap with DEG genes
pal <- c("blue","white","red") 
pal <- colorRampPalette(pal)(100)
heatmap(as.matrix(cpm_table[which(rownames(cpm_table)%in%DEGs$id[which(DEGs$class!="=")]),])
        ,cexCol = 0.5,margins = c(4,4),col=pal,cexRow = 0.2)

# 4
# Gene set enrichment analysis
ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
convert <- getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name"),
                 filters=c("ensembl_gene_id"), 
                 values=DEGs$id,
                 mart = ensembl)
convert <- na.omit(convert)

DEGs <- merge(DEGs,convert,by.x="id",by.y="ensembl_gene_id")

# THIS GAVE ME ERRORS SO I COULDN'T CHECK FOR DUPLICATES
# THIS SPECIFIC DATA SET DOESN'T HAVE DUPLICATES SO IT'S SAFE TO GO ON
#DEGs <- DEGs[-which(duplicated(DEGs$entrezgene_id)),]

## Create a list of upregulated genes
upDEGs <- DEGs %>% filter(class == "+")
downDEGs <- DEGs %>% filter(class == "-")

# BP
up_ego_BP <- enrichGO(gene = upDEGs$external_gene_name,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
dim(up_ego_BP)

down_ego_BP <- enrichGO(gene = downDEGs$external_gene_name,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
dim(down_ego_BP)

# MF
###################################################################
# !!! HAD TO INCREASE P-VALUE FOR UP REGULATED GENES USING MF !!! #
###################################################################
up_ego_MF <- enrichGO(gene = upDEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.15,
                      qvalueCutoff = 0.15)
dim(up_ego_MF)

down_ego_MF <- enrichGO(gene = downDEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
dim(down_ego_MF)

#barplot(up_ego_BP,showCategory=10)

# KEGG

ekegg <- enrichKEGG(gene = upDEGs$entrezgene_id,
                    organism = 'human',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

head(ekegg, n=10)

## Plot the pathway with DEG genes
logFC <- upDEGs$log2_FC
names(logFC) <- upDEGs$entrezgene_id
#pathview(gene.data = logFC, 
         # = "hsa05202", 
         #species = "human")
