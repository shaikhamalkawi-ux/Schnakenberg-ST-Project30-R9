#!/usr/bin/env Rscript
# GSE97358 donor-blocked Negative Binomial sensitivity; NO substitution or raw-data redistribution.
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE)
dir.create("results",showWarnings=FALSE)
files <- c("GSE97358_fibroblasts_samples_info.txt.gz","GSE97358_168_fibroblast_samples_gene_counts.txt.gz")
expected <- c("0c1bb7564e7c943ab68a4ec9a486da473a44d92a89081b9a5cd001510afcf747",
              "f07c0230c8824b91cb73c28c0d544df989073a38172a2506167692cb6fb460d1")
for(f in files) if(!file.exists(file.path("raw",f))) stop("Missing GEO raw file: ",f)
observed <- vapply(files,function(f) strsplit(system2("sha256sum",file.path("raw",f),stdout=TRUE),"\\s+")[[1]][1],character(1))
if(!identical(unname(observed),unname(expected))) stop("Official GEO source hash changed")
read_tab <- function(name) read.delim(gzfile(file.path("raw",name)),check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE)
meta <- read_tab(files[1])
if(!"Sample_ID" %in% names(meta) || anyDuplicated(meta$Sample_ID)) stop("Sample info integrity failure")
raw <- read_tab(files[2])
if(!"Ensembl_geneID" %in% names(raw) || anyDuplicated(raw$Ensembl_geneID)) stop("Count row ID integrity failure")
genes <- raw$Ensembl_geneID
raw$Ensembl_geneID <- NULL
counts <- as.matrix(raw)
mode(counts) <- "numeric"
rownames(counts) <- genes
if(anyNA(counts) || any(counts<0) || any(abs(counts-round(counts))>1e-9)) stop("Raw counts are not nonnegative integers")
samples <- colnames(counts)
if(ncol(counts)!=168 || length(meta$Sample_ID)!=168 || !setequal(samples,meta$Sample_ID) || anyDuplicated(samples)) stop("Exactly 168 source-aligned samples required")
if(!all(grepl("^Pa[0-9]+(BS|ST)$",samples))) stop("Unrecognized sample-name encoding")
donor <- sub("(BS|ST)$","",samples)
group <- factor(sub("^Pa[0-9]+","",samples),levels=c("BS","ST"))
if(length(unique(donor))!=84 || any(table(donor,group)!=1)) stop("Expected 84 donor-linked basal/stimulated pairs")
ids <- c(CTHRC1="ENSG00000164932",POSTN="ENSG00000133110",FN1="ENSG00000115414",
         LOX="ENSG00000113083",COL1A1="ENSG00000108821",COL3A1="ENSG00000168542",
         ASPN="ENSG00000106819")
if(!all(ids %in% rownames(counts))) stop("Missing predefined Ensembl target")
coldata <- data.frame(donor=factor(donor),condition=group,row.names=samples)
design <- model.matrix(~donor+condition,data=coldata)
if(qr(design)$rank!=ncol(design) || tail(colnames(design),1)!="conditionST") stop("Non-identifiable paired design")
# Independent primary-endpoint replay: original median-ratio normalized log2(count+0.5).
nonzero <- rowSums(counts>0)==168
if(sum(nonzero)<100) stop("Insufficient genes to estimate original median-ratio normalization")
geom <- exp(rowMeans(log(counts[nonzero,,drop=FALSE])))
sf <- apply(counts[nonzero,,drop=FALSE],2,function(v) median(v/geom))
if(any(!is.finite(sf)) || any(sf<=0)) stop("Invalid original count scaling")
b <- paste0(unique(donor),"BS")
t <- paste0(unique(donor),"ST")
stopifnot(all(b %in% samples),all(t %in% samples))
paired <- log2(counts[ids["CTHRC1"],t]/sf[t]+0.5)-log2(counts[ids["CTHRC1"],b]/sf[b]+0.5)
replay <- list(mean=unname(as.numeric(mean(paired))),median=unname(as.numeric(median(paired))),
               positive=as.integer(sum(paired>0)),n=as.integer(length(paired)))
# edgeR QL Negative Binomial GLM with TMM offsets; donor as a fixed blocking factor.
dge <- DGEList(counts=counts,group=group)
keep <- filterByExpr(dge,group=group)
if(sum(keep)<1000 || !keep[match(ids["CTHRC1"],rownames(dge))]) stop("Insufficient expressed genes")
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge,method="TMM")
dge <- estimateDisp(dge,design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
test <- glmQLFTest(fit,coef="conditionST")
all <- test$table
if(any(!is.finite(all$logFC))) stop("Nonfinite model coefficients")
all$q_BH_all <- p.adjust(all$PValue,"BH")
ix <- match(ids,rownames(all))
if(anyNA(ix)) stop("Target filtered out by expression rule")
sub <- all[ix,,drop=FALSE]
result <- data.frame(gene=names(ids),EnsemblID=unname(ids),n_pairs=84,
  log2FC_NB_TMM=unname(sub$logFC),F_QL=unname(sub$F),p_unadjusted=unname(sub$PValue),
  q_BH_all_filtered=unname(sub$q_BH_all),q_BH_7_prespecified=p.adjust(sub$PValue,"BH"),
  row.names=NULL)
write.csv(result,"results/GSE97358_NB_7genes.csv",row.names=FALSE)
status <- if(replay$positive==82 && abs(replay$mean-0.882)<0.015) "PRIMARY_REPLAY_APPROX_PASS" else "PRIMARY_REPLAY_MISMATCH_HOLD"
if(!is.finite(result$log2FC_NB_TMM[1]) || !is.finite(result$p_unadjusted[1])) status <- "NB_NONFINITE_HOLD"
report <- list(status=status,source="Official NCBI GEO GSE97358",source_gzip_sha256=as.list(unname(observed)),
 n_samples=168,n_control=84,n_treated=84,n_complete_donor_pairs=84,
 original_replay=replay,original_expected=c(mean=0.882,positive=82),
 model="edgeR robust quasi-likelihood negative-binomial GLM, donor fixed effect + groupST, TMM normalization",
 filter="edgeR filterByExpr defaults, grouped BS/ST",filtered_genes=sum(keep),
 original_genes=nrow(counts),rank=qr(design)$rank,design_columns=ncol(design),
 R_version=R.version.string,edgeR_version=as.character(packageVersion("edgeR")),
 CTHRC1=list(log2FC_NB_TMM=as.numeric(result$log2FC_NB_TMM[1]),
 p_unadjusted=as.numeric(result$p_unadjusted[1]),q_BH_all_filtered=as.numeric(result$q_BH_all_filtered[1])),
 estimand_boundary="NB model coefficient is not directly equal to the original mean donor-paired log2 normalized expression",
 raw_counts_redistributed=FALSE)
jsonlite::write_json(report,"results/GSE97358_NB_summary.json",auto_unbox=TRUE,pretty=TRUE)
writeLines(c(paste("STATUS",status),paste("ORIGINAL_MEAN",signif(replay$mean,8)),
 paste("ORIGINAL_POSITIVE",replay$positive,"of",replay$n),
 paste("NB_LOG2FC",signif(result$log2FC_NB_TMM[1],8)),
 paste("NB_PVALUE",signif(result$p_unadjusted[1],8)),
 paste("NB_Q_ALL",signif(result$q_BH_all_filtered[1],8)),
 paste("FILTERED",sum(keep)),"RAW_DATA_NOT_REDISTRIBUTED"),
 "results/console_summary.txt")
pdf("results/CTHRC1_84_pair_changes.pdf",width=7.5,height=4.2,useDingbats=FALSE)
par(mar=c(4,4.7,1,0.4))
plot(seq_along(paired),sort(as.numeric(paired)),xlab="84 donor pairs, ordered by contrast",
 ylab="Treated - basal (log2 normalized counts + 0.5)",
 pch=19,col="gray35",cex=0.6)
abline(h=0,lty=2,lwd=1.2)
abline(h=mean(paired),lwd=1.3)
legend("topleft",legend=c("Zero","Mean contrast"),lty=c(2,1),bty="n",cex=.8)
dev.off()
cat(paste(readLines("results/console_summary.txt"),collapse="\n"),"\n")
if(status!="PRIMARY_REPLAY_APPROX_PASS") stop("Scientific HOLD; original primary analysis did not match")
