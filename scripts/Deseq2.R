#!/usr/bin/Rscript --vanilla
######
# Code written by Anais Bardet (2023)
######

## Retrieve arguments
args=commandArgs(TRUE)

## Help
help <- function(){
    cat("\nDESeq2 pipeline \n")
    cat("Usage: DESeq2_pipeline.R -n STUDY_DESIGN_NAME -d STUDY_DESIGN_FILE [-c COUNTS_DIR] [-o OUTPUT_DIR]\n")
    cat("-n : Study design name used in output files\n")
    cat("-d : Study design file: sample / condition / replicate (column after sample is used for design) (at least 3 columns required)\n")
    cat("-c : Directory containing count files (default: ./expression/deseq2/counts)\n")
    cat("-o : Output directory (default: ./results/expression/deseq2)\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    for(ii in 1:length(args)){
	if(grepl("^-",args[ii]) && args[ii] != "-"){
	    if(ii+1<=length(args) && (!grepl("^-",args[ii+1]) || args[ii+1]=="-")){
	        assign(gsub("-","",args[ii]),args[ii+1])
	    } else { assign(gsub("-","",args[ii]),1) }
	}
}}

## Set defaults for optional arguments
if(exists("c")) { countsdir <- c; rm("c") } else countsdir <- "./expression/deseq2/counts"
if(!exists("o")) o <- "./results/expression/deseq2"
outdir <- o
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

## Load data from input into a table
if(exists("d")){
    design=n
    if (file.exists(d)){
        colData=as.data.frame(read.table(d,header=T,row.names=1))
	samples=rownames(colData)
	cat(paste0("Design using ",colnames(colData)[1],"\n"))
    } else { cat("Study design file does not exist\n"); q() }
} else { cat("No input specified\n"); q() }

suppressMessages(library(DESeq2))
suppressMessages(library(NMF))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

# CountData
d=vector()
for(sample in samples){
    x=read.table(paste0(countsdir,"/",sample,"_genes_count.txt"))
    d=cbind(d,x[,3])
    colnames(d)[ncol(d)]=sample
    rownames(d)=x[,1]
    genenames=x[,2]
    genelength=x[,4]
}

## DESeq matrix
dds=suppressMessages(DESeqDataSetFromMatrix(countData=d,colData=colData,design=as.formula(paste0("~",colnames(colData)[[1]]))))
# Filter out genes with no counts
dds=dds[rowSums(counts(dds))>1,]

cat("Data loaded\n")

## Transform data for variance analysis
# changed to use varianceStabilizingTransformation for small datasets (less than 1000 genes) to avoid errors
rld=suppressMessages(
  if(nrow(dds) >= 1000) vst(dds, blind=FALSE)
  else varianceStabilizingTransformation(dds, blind=FALSE)
)

## Principal components analysis (PCA)
pcaData=plotPCA(rld,intgroup=colnames(colData),returnData=T)
percentVar=round(100*attr(pcaData,"percentVar"))
pdf(paste0(outdir,"/",design,"_deseq2_pca.pdf"))
par(bg="white")
ggplot(pcaData,aes(x=PC1,y=PC2,color=as.character(get(colnames(colData)[1])),shape=as.character(get(colnames(colData)[2])))) +
  ggtitle("Gene expression") +
  geom_point(size=3) +
  scale_shape_manual(values=seq(0,length(unique(pcaData$condition)))) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color=colnames(colData)[1])  +
  labs(shape=colnames(colData)[2])  +
  coord_fixed()
null=dev.off()

cat("PCA plot generated\n")
cat(paste0("=> ",outdir,"/",design,"_deseq2_pca.pdf\n"))

## Differential expression (on non transformed data)

# Contrasts
conds=unique(colData[[1]])
res=list()
contrasts=c()
for(i in 1:(length(conds)-1)){
    for(j in (i+1):(length(conds))){
    	  contrast=paste0(conds[i],"_vs_",conds[j])
    	  sel=which(colData[,1]==conds[i]|colData[,1]==conds[j])
	  selColData=colData[sel,]
	  seld=d[,sel]
	  dds=suppressMessages(DESeqDataSetFromMatrix(countData=seld,colData=selColData,design=~condition))
	  dds=dds[rowSums(counts(dds))>1,]
	  design(dds) <- ~ condition
	  dds=suppressMessages(DESeq(dds))

        res[[contrast]]=suppressMessages(results(dds,contrast=c(colnames(colData)[1],as.character(conds[i]),as.character(conds[j])),cooksCutoff=FALSE))
        res[[contrast]]$log2FoldChange[which(is.na(res[[paste0(conds[i],"_vs_",conds[j])]]$padj))]=NA
        res[[contrast]]$padj[which(res[[contrast]]$padj==0)]=10^-30
	contrasts=c(contrasts,contrast)
	cat(contrast," comparison performed\n")

	mean=list()
	for(cond in c(as.character(conds[i]),as.character(conds[j]))){
	    mean[[cond]]=apply(counts(dds,normalized=T)[,which(as.character(dds[[colnames(selColData)[1]]])==cond)],1,mean)
	}

	# Result table
	write.table(cbind("id"=rownames(res[[contrast]]),"gene"=genenames[match(rownames(res[[contrast]]),rownames(seld))],data.frame(sapply(mean,c)),res[[contrast]][,c(2,6)]),paste0(outdir,"/",design,"_",contrast,"_deseq2_diff_table.txt"),quote=F,sep="\t",row.names=F,col.names=T)

	cat("Summary table generated\n")
	cat(paste0("=> ",outdir,"/",design,"_",contrast,"_deseq2_diff_table.txt\n"))

}}

max_pval=ceiling(max(-log10(unlist(sapply(res,"[[",6))),na.rm=T))
max_fc=ceiling(max(abs(unlist(sapply(res,"[[",2))),na.rm=T))

# Volcano plot for each condition
for(contrast in contrasts){
    pdf(paste0(outdir,"/",design,"_",contrast,"_deseq2_diff_volcano.pdf"),useDingbats=FALSE)
    par(bg="white")
    sel=head(order(-log10(res[[contrast]]$padj),decreasing=T),20)
    labels=genenames[match(rownames(res[[contrast]]),rownames(d))]
    labels[setdiff(1:length(rownames(res[[contrast]])),sel)]=NA
p=ggplot(as.data.frame(res[[contrast]]),aes(x=log2FoldChange,y=-log10(padj),label=labels)) +
    ggtitle(contrast) +
    geom_point(size=1) +
    xlab("Fold Change (log2)") +
    ylab("P-value (-log10)") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_vline(xintercept=c(-1, 1), color="red", linetype="dashed") +
    geom_text_repel(size=3)
    n_up   <- sum(res[[contrast]]$log2FoldChange >  1 & res[[contrast]]$padj < 0.05, na.rm=T)
    n_down <- sum(res[[contrast]]$log2FoldChange < -1 & res[[contrast]]$padj < 0.05, na.rm=T)

    p = p + annotate("text", x =  max_fc * 0.8, y = max_pval * 0.95, 
                    label = paste0("Up: ", n_up), color = "red", size = 4, hjust = 1) +
            annotate("text", x = -max_fc * 0.8, y = max_pval * 0.95, 
                    label = paste0("Down: ", n_down), color = "blue", size = 4, hjust = 0)
    suppressWarnings(print(p))
    null=dev.off()
    cat("Volcano plot generated\n")
    cat(paste0("=> ",outdir,"/",design,"_",contrast,"_deseq2_diff_volcano.pdf\n"))
   }


