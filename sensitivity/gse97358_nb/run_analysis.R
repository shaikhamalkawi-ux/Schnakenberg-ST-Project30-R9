#!/usr/bin/env Rscript
# Project 30: sensitivity of the existing seven-gene panel, not new discovery.
# Raw GEO gene counts are not FASTQ read alignment. No count imputation/replacement.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
rawdir <- args[1]; outdir <- args[2]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
options(stringsAsFactors=FALSE, digits=17)
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
writeLines(capture.output(sessionInfo()),file.path(outdir,'R_sessionInfo.txt'))
versions <- sapply(c('DESeq2','S4Vectors','SummarizedExperiment','BiocParallel','jsonlite'),function(x)as.character(packageVersion(x)))
write.csv(data.frame(package=names(versions),version=versions),file.path(outdir,'software_versions.csv'),row.names=FALSE)
protocol <- list(dataset='GSE97358',role='post-review sensitivity of previously specified targets',primary='CTHRC1',secondary=c('POSTN','FN1','LOX','COL1A1','COL3A1','ASPN'),design='~ donor + treatment',reference='Control',contrast='TGFB1 minus Control at 24 h',normalization='DESeq2 median-ratio size factors estimated on full input before filtering',filter='at least 10 counts in at least 84 samples; fixed before model fit',fitType='parametric',test='two-sided Wald',betaPrior=FALSE,lfcShrink=FALSE,minReplicatesForReplace='Inf (no count replacement)',independentFiltering=FALSE,cooksCutoff='default; report explicit influence diagnostic',multiplicity='primary raw P; BH across six predefined secondary genes; BH across all tested rows separately',decision='report all seven genes regardless of result; no primary-endpoint replacement or clinical inference')
write_json(protocol,file.path(outdir,'ANALYSIS_SPECIFICATION_BEFORE_FIT.json'),pretty=TRUE,auto_unbox=TRUE)
files <- c('GSE97358_168_fibroblast_samples_gene_counts.txt.gz','GSE97358_fibroblasts_samples_info.txt.gz')
expected <- c('f07c0230c8824b91cb73c28c0d544df989073a38172a2506167692cb6fb460d1','0c1bb7564e7c943ab68a4ec9a486da473a44d92a89081b9a5cd001510afcf747')
observed <- vapply(file.path(rawdir,files),function(x)strsplit(system2('sha256sum',shQuote(x),stdout=TRUE),' ')[[1]][1],character(1))
stopifnot(identical(unname(observed),expected))
cm <- read.delim(gzfile(file.path(rawdir,files[1])),check.names=FALSE,row.names=1)
meta <- read.delim(gzfile(file.path(rawdir,files[2])),check.names=FALSE)
stopifnot(ncol(cm)==168L,nrow(cm)==64253L,nrow(meta)==168L,!anyDuplicated(rownames(cm)),!anyDuplicated(meta$Sample_ID),setequal(colnames(cm),meta$Sample_ID))
cm <- as.matrix(cm)
stopifnot(all(is.finite(cm)),all(cm>=0),all(cm==floor(cm)),max(cm)<.Machine$integer.max)
storage.mode(cm)<-'integer'
meta<-meta[match(colnames(cm),meta$Sample_ID),,drop=FALSE]
stopifnot(identical(meta$Sample_ID,colnames(cm)))
stopifnot(all(grepl('^Pa[0-9]{3}(BS|ST)$',meta$Sample_ID)))
meta$donor<-factor(sub('(BS|ST)$','',meta$Sample_ID))
stopifnot(all(meta$Stimulation %in% c('Control-24h','TGFB1-24h')))
meta$treatment<-factor(ifelse(meta$Stimulation=='Control-24h','Control','TGFB1'),levels=c('Control','TGFB1'))
stopifnot(all((grepl('BS$',meta$Sample_ID))==(meta$treatment=='Control')))
stopifnot(nlevels(meta$donor)==84L,all(table(meta$donor,meta$treatment)==1))
rownames(meta)<-meta$Sample_ID
batch_per_donor<-tapply(meta$LibraryPrep_Batch,meta$donor,function(x)length(unique(x)))
stopifnot(all(batch_per_donor==1L))
design<-model.matrix(~donor+treatment,meta)
stopifnot(qr(design)$rank==85L,nrow(design)-qr(design)$rank==83L)
write.csv(meta,file.path(outdir,'validated_sample_manifest.csv'),row.names=FALSE)
write.csv(design,file.path(outdir,'design_matrix.csv'))
targets<-data.frame(gene=c('CTHRC1','POSTN','FN1','LOX','COL1A1','COL3A1','ASPN'),ensembl=c('ENSG00000164932','ENSG00000133110','ENSG00000115414','ENSG00000113083','ENSG00000108821','ENSG00000168542','ENSG00000106819'))
stopifnot(all(targets$ensembl %in% rownames(cm)))
write.csv(targets,file.path(outdir,'target_gene_map.csv'),row.names=FALSE)
dds_all<-DESeqDataSetFromMatrix(countData=cm,colData=meta,design=~donor+treatment)
dds_all<-estimateSizeFactors(dds_all,type='ratio')
sf<-sizeFactors(dds_all)
stopifnot(all(is.finite(sf)),all(sf>0))
write.csv(data.frame(sample=names(sf),size_factor=unname(sf),library_total=colSums(cm)),file.path(outdir,'size_factors.csv'),row.names=FALSE)
keep<-rowSums(cm>=10L)>=84L
write.csv(data.frame(ensembl=rownames(cm),kept=keep,nonzero_samples=rowSums(cm>0),samples_ge10=rowSums(cm>=10),total_count=rowSums(cm)),file.path(outdir,'gene_filter_manifest.csv'),row.names=FALSE)
stopifnot(all(keep[targets$ensembl]))
# Before any model inference, export donor-level transformed observations.
norm<-sweep(cm[targets$ensembl,,drop=FALSE],2,sf,'/')
logexpr<-log2(norm+0.5)
ctrl<-match(paste0(levels(meta$donor),'BS'),colnames(cm));stim<-match(paste0(levels(meta$donor),'ST'),colnames(cm))
diff<-logexpr[,stim,drop=FALSE]-logexpr[,ctrl,drop=FALSE]
paired<-do.call(rbind,lapply(seq_len(nrow(targets)),function(i)data.frame(gene=targets$gene[i],ensembl=targets$ensembl[i],donor=levels(meta$donor),control_sample=colnames(cm)[ctrl],treated_sample=colnames(cm)[stim],control_count=cm[targets$ensembl[i],ctrl],treated_count=cm[targets$ensembl[i],stim],control_size_factor=sf[ctrl],treated_size_factor=sf[stim],control_normalized=norm[i,ctrl],treated_normalized=norm[i,stim],paired_log2_difference=diff[i,])))
write.csv(paired,file.path(outdir,'seven_gene_donor_pairs.csv'),row.names=FALSE)
summary<-data.frame(gene=targets$gene,mean=rowMeans(diff),median=apply(diff,1,median),positive=rowSums(diff>0),negative=rowSums(diff<0),zero=rowSums(diff==0),min=apply(diff,1,min),max=apply(diff,1,max),loo_min=apply(diff,1,function(x)min((sum(x)-x)/83)),loo_max=apply(diff,1,function(x)max((sum(x)-x)/83)))
write.csv(summary,file.path(outdir,'transformed_endpoint_reconstruction.csv'),row.names=FALSE)
write_json(list(input_genes=nrow(cm),samples=ncol(cm),donors=nlevels(meta$donor),pairs_complete=TRUE,design_rank=qr(design)$rank,residual_df=83L,all_pairs_same_libraryprep_batch=all(batch_per_donor==1L),genes_positive_all_samples=sum(rowSums(cm>0)==168L),tested_genes=sum(keep),excluded_low_counts=sum(!keep),all_target_genes_kept=TRUE,zero_imputed_counts=TRUE),file.path(outdir,'INPUT_QA.json'),pretty=TRUE,auto_unbox=TRUE)
cat('INPUT GATES PASS. Retained genes:',sum(keep),'rank:',qr(design)$rank,'\n')
dds<-dds_all[keep,]
rm(dds_all)
dds<-DESeq(dds,test='Wald',fitType='parametric',betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE,quiet=FALSE)
stopifnot(identical(counts(dds),cm[keep,,drop=FALSE]))
res<-results(dds,contrast=c('treatment','TGFB1','Control'),independentFiltering=FALSE,alpha=0.05)
allres<-data.frame(ensembl=rownames(res),as.data.frame(res),beta_converged=mcols(dds)$betaConv,dispersion=dispersions(dds))
write.csv(allres,file.path(outdir,'DESeq2_all_tested_genes.csv'),row.names=FALSE)
rr<-allres[match(targets$ensembl,allres$ensembl),]
rr<-data.frame(gene=targets$gene,rr)
rr$CI95_lower<-rr$log2FoldChange-qnorm(0.975)*rr$lfcSE
rr$CI95_upper<-rr$log2FoldChange+qnorm(0.975)*rr$lfcSE
rr$fold_change<-2^rr$log2FoldChange
rr$q_BH_secondary6<-NA_real_;rr$q_BH_secondary6[-1]<-p.adjust(rr$pvalue[-1],method='BH')
rr$disp_gene_est<-mcols(dds)$dispGeneEst[match(targets$ensembl,rownames(dds))]
rr$disp_fitted_trend<-mcols(dds)$dispFit[match(targets$ensembl,rownames(dds))]
cooks<-assays(dds)[['cooks']][targets$ensembl,,drop=FALSE]
rr$max_cooks<-apply(cooks,1,max,na.rm=TRUE)
rr$cooks_F99_cutoff<-qf(.99,ncol(design),nrow(design)-ncol(design))
rr$cooks_above_cutoff<-rr$max_cooks>rr$cooks_F99_cutoff
rr$dispersion_outlier<-mcols(dds)$dispOutlier[match(targets$ensembl,rownames(dds))]
write.csv(rr,file.path(outdir,'DESeq2_seven_gene_results.csv'),row.names=FALSE)
write.csv(coef(dds)[targets$ensembl,,drop=FALSE],file.path(outdir,'target_coefficients_log2.csv'))
write.csv(assays(dds)[['mu']][targets$ensembl,,drop=FALSE],file.path(outdir,'target_fitted_counts.csv'))
write.csv(cooks,file.path(outdir,'target_cooks_distance.csv'))
write.csv(data.frame(ensembl=rownames(dds),gene_dispersion=mcols(dds)$dispGeneEst,trend=mcols(dds)$dispFit,MAP=dispersions(dds)),file.path(outdir,'dispersion_diagnostics.csv'),row.names=FALSE)
write_json(list(dispersion_fit_type=attr(dispersionFunction(dds),'fitType'),dispersion_fit_coefficients=attr(dispersionFunction(dds),'coefficients'),target_all_beta_converged=all(rr$beta_converged),target_all_p_finite=all(is.finite(rr$pvalue)),target_all_intervals_finite=all(is.finite(rr$CI95_lower)&is.finite(rr$CI95_upper)),target_influence_flags=sum(rr$cooks_above_cutoff),nonconverged_genomewide=sum(!allres$beta_converged,na.rm=TRUE),NA_p_genomewide=sum(is.na(allres$pvalue)),count_replacement=FALSE,primary_endpoint_not_overwritten=TRUE),file.path(outdir,'MODEL_QA.json'),pretty=TRUE,auto_unbox=TRUE)
print(rr)
stopifnot(all(rr$beta_converged),all(is.finite(rr$pvalue)),all(is.finite(rr$CI95_lower)),all(is.finite(rr$CI95_upper)))
# Compact full audit object: source counts stay in raw/, no full DESeq object needed.
writeLines(capture.output(sessionInfo()),file.path(outdir,'R_sessionInfo.txt'))
cat('SEVEN-GENE COUNT-MODEL SENSITIVITY COMPLETED\n')
