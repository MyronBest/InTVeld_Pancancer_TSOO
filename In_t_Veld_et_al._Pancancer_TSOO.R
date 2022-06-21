# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.g.best@amsterdamumc.nl; g.intveld1@amsterdamumc.nl
# Summary       : Script to reproduce the main results (pan-cancer classifier and tumor-site-of-origin classifier)
#                 as presented in the manuscript by In 't Veld et al.
# Date          : 15th of June 2022
# Note          : This script builds on the previously published pipeline by Best et al. Nature Protocols 2019,
#                 and requires the same input scripts and R-packages as described.

# load the required functions by running the following commands: 
source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO.R')
nCores = 50

########################
#      PANCANCER       #
########################

# load and prepare the input data provided in GEO (GSE183635). Contains the raw read counts
# from all included samples (n=2351), and all (pre-filtered) counts (n=5440)
load("TEP_Count_Matrix.RData")
colnames(TEP_Count_Matrix) <- sub(".*?-", "", colnames(TEP_Count_Matrix))

# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('sampleInfo_libSize.csv', sep=",", row.names=1)
sampleInfo <- sampleInfo[order(match(rownames(sampleInfo), colnames(TEP_Count_Matrix))), ]

# load gene info and create edgeR object
load('bin/dgeGenesEnsembl75.RData')
library(edgeR)
dge <- DGEList(counts = TEP_Count_Matrix,
               group = sampleInfo$Group,
               genes = genes[which(rownames(genes) %in% rownames(TEP_Count_Matrix)),]
)
dge$samples <- cbind(dge$samples, sampleInfo)
dge$samples$lib.size <- sampleInfo$lib.size

# dichotomize the groups (i.e. Cancer and asymptomatic controls)
dge$samples$group <- factor(dge$samples$group, levels = c(levels(dge$samples$group),"Cancer","asymptomaticControls"))
dge$samples$group[which(dge$samples$group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                 "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                 "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                 "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer"))] <- "Cancer"
# here all non-cancer groups are names 'asymptomaticControls' to enable validation purposes
dge$samples$group[which(dge$samples$group %in% c("Asymptomatic controls","Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                 "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                 "Pancreatic diseases","Pulmonary Hypertension"))] <- "asymptomaticControls"
dge$samples <- droplevels(dge$samples)
summary(dge$samples$group)

# Perform thromboSeq ANOVA differential expression analysis of splice junctions
thromboSeq.anova <- thromboSeqANOVA(dge = dge[,c(colnames(dge)[dge$samples$Training.series==1],
                                                 colnames(dge)[dge$samples$Evaluation.series==1])],
                                    k.variables = 2,
                                    variable.to.assess = c("lib.size"),
                                    variable.threshold = c(0.8),
                                    ruvg.pvalue.threshold.group = 1e-2,
                                    ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                    training.series.only = FALSE,
                                    select.biomarker.FDR = TRUE,
                                    plot = FALSE,
                                    iteration = NULL,
                                    figureDir = "figureOutputFolder",
                                    number.cores = nCores,
                                    verbose = TRUE)

# Perform PSO-enhanced thromboSeq classifier development.
thromboPSO <- thromboSeqPSO(dge = dge[,colnames(dge)[dge$samples$group %in% c("asymptomaticControls","Cancer")]], 
                            training.samples.provided = colnames(dge)[dge$samples$Training.series==1],
                            evaluation.samples.provided = colnames(dge)[dge$samples$Evaluation.series==1],
                            swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"),
                            swarm.boundaries = c(-0.1, 1.0, 50, as.numeric(summary(thromboSeq.anova$FDR < 0.005)[3]), 0.5, 1.0, 50, as.numeric(summary(thromboSeq.anova$FDR < 0.005)[3])),
                            k.variables = 2,
                            variable.to.assess = c("lib.size"),
                            variable.threshold = c(0.8),
                            ruvg.pvalue.threshold.group = 1e-8,
                            ruvg.pvalue.threshold.strongest.variable = 1e-2,
                            select.biomarker.FDR = TRUE,
                            minimum.n.transcripts.biomarkerpanel = 2,
                            svm.gamma.range = 2^(-20:0),
                            svm.cost.range = 2^(0:20),
                            number.cross.splits = 2,
                            n.particles.gamma.cost.optimization = 50,
                            n.iterations.gamma.cost.optimization = 4,
                            n.particles = 60,
                            n.iterations = 8,
                            figureDir = "figureOutputFolder",
                            number.cores = nCores, 
                            rule.in.optimization = TRUE,
                            verbose = TRUE
)

# Validate the developed thromboSeq algorithm.
thromboPSOreadout <- thromboSeqPSO.readout(dge = dge,
                                           replace.counts.validation = 0, 
                                           filter.clinical.characteristics.validation = NA,
                                           filter.clinical.characteristics.group = NA,
                                           filter.clinical.characteristics.specified = NA,
                                           readout.training = T, 
                                           readout.evaluation = T, 
                                           readout.validation = T, 
                                           apply.rule.in.readout = T, 
                                           rule.in.setting = 99, 
                                           apply.rule.out.readout = F, 
                                           rule.out.setting = NA,
                                           additional.dge.to.predict = NA,
                                           number.cores = nCores, 
                                           clinical.data.in.output = c('Sample.ID','Group','Sample.supplying.institution','Stage')
)


# Perform control experiments for thromboSeq classifier development.
control <- thromboSeqPSO.controls(dge = dge,
                                  filter.clinical.characteristics.validation = NA,
                                  filter.clinical.characteristics.group = NA,
                                  filter.clinical.characteristics.specified = NA,
                                  thromboSeqPSO.shuffled = T,
                                  thromboSeqPSO.iterations = F,
                                  number.cores = nCores, 
                                  n.shuffled = 1000
)

control <- thromboSeqPSO.controls(dge = dge,
                                  filter.clinical.characteristics.validation = NA,
                                  filter.clinical.characteristics.group = NA,
                                  filter.clinical.characteristics.specified = NA,
                                  thromboSeqPSO.shuffled = F,
                                  thromboSeqPSO.iterations = T,
                                  number.cores = nCores, 
                                  n.iterations = 1000
)

# create figures from the output data
load('results.classification.training.RData')
load('results.classification.evaluation.rule.in.RData')
load('results.classification.validation.rule.in.RData')

# ROC curve
library(ROCR)
library(pROC)

## training series
rocra.training <- prediction(as.numeric(as.character(results.classification.training$svm.summary[, 5])), 
                             results.classification.training$svm.summary[, 3], 
                             label.ordering = rev(levels(results.classification.training$svm.summary$real.group)))
perfa.training <- performance(rocra.training, "tpr", "fpr")
print(paste("AUC Training Series:", attributes(performance(rocra.training, 'auc'))$y.values[[1]]))

roc.summary <- data.frame(
  cutOffs = unlist(attributes(rocra.training)$cutoffs),
  tp = unlist(attributes(rocra.training)$tp),
  tn = unlist(attributes(rocra.training)$tn),
  fp = unlist(attributes(rocra.training)$fp),
  fn = unlist(attributes(rocra.training)$fn),
  accuracy = (unlist(attributes(rocra.training)$tp) + unlist(attributes(rocra.training)$tn)) /
    (unlist(attributes(rocra.training)$fp) + unlist(attributes(rocra.training)$tn) + unlist(attributes(rocra.training)$tp) + unlist(attributes(rocra.training)$fn)),
  xValues = unlist(attributes(perfa.training)$x.values),
  yValues = unlist(attributes(perfa.training)$y.values)
)
roc.optimal.accuracy <- max(roc.summary$accuracy)

# calculate confidence interval
roc.95ci <- roc(results.classification.training$svm.summary$real.group,
                results.classification.training$svm.summary$asymptomaticControls, 
                ci = TRUE
)


## evaluation series
rocra.evaluation <- prediction(as.numeric(as.character(results.classification.evaluation.rule.in$svm.summary[, 4])), 
                               results.classification.evaluation.rule.in$svm.summary[, 3], 
                               label.ordering = rev(levels(results.classification.evaluation.rule.in$svm.summary$real.group)))
perfa.evaluation <- performance(rocra.evaluation, "tpr", "fpr")
print(paste("AUC Evaluation Series: ", attributes(performance(rocra.evaluation, 'auc'))$y.values[[1]], sep = ""))

roc.summary <- data.frame(
  cutOffs = unlist(attributes(rocra.evaluation)$cutoffs),
  tp = unlist(attributes(rocra.evaluation)$tp),
  tn = unlist(attributes(rocra.evaluation)$tn),
  fp = unlist(attributes(rocra.evaluation)$fp),
  fn = unlist(attributes(rocra.evaluation)$fn),
  accuracy = (unlist(attributes(rocra.evaluation)$tp) + unlist(attributes(rocra.evaluation)$tn)) /
    (unlist(attributes(rocra.evaluation)$fp) + unlist(attributes(rocra.evaluation)$tn) + unlist(attributes(rocra.evaluation)$tp) + unlist(attributes(rocra.evaluation)$fn)),
  xValues = unlist(attributes(perfa.evaluation)$x.values),
  yValues = unlist(attributes(perfa.evaluation)$y.values)
)
roc.optimal.accuracy <- max(roc.summary$accuracy)

# calculate confidence interval
roc.95ci <- roc(results.classification.evaluation.rule.in$svm.summary$real.group,
                results.classification.evaluation.rule.in$svm.summary$asymptomaticControls, 
                ci = TRUE
)

## validation series
# pancancer ROC
samples <- rownames(results.classification.validation.rule.in$svm.summary)[!results.classification.validation.rule.in$svm.summary$Group %in% 
                                                                             c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                               "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                               "Pancreatic diseases","Pulmonary Hypertension")]
results.classification.validation.rule.in$svm.summary <- results.classification.validation.rule.in$svm.summary[samples,]

rocra.validation <- prediction(as.numeric(as.character(results.classification.validation.rule.in$svm.summary[, 4])), 
                               results.classification.validation.rule.in$svm.summary[, 3], 
                               label.ordering = rev(levels(droplevels(results.classification.validation.rule.in$svm.summary$real.group))))
perfa.validation <- performance(rocra.validation, "tpr", "fpr")
print(paste("AUC Validation Series: ", attributes(performance(rocra.validation, 'auc'))$y.values[[1]], sep = ""))

roc.summary <- data.frame(
  cutOffs = unlist(attributes(rocra.validation)$cutoffs),
  tp = unlist(attributes(rocra.validation)$tp),
  tn = unlist(attributes(rocra.validation)$tn),
  fp = unlist(attributes(rocra.validation)$fp),
  fn = unlist(attributes(rocra.validation)$fn),
  accuracy = (unlist(attributes(rocra.validation)$tp) + unlist(attributes(rocra.validation)$tn)) /
    (unlist(attributes(rocra.validation)$fp) + unlist(attributes(rocra.validation)$tn) + unlist(attributes(rocra.validation)$tp) + unlist(attributes(rocra.validation)$fn)),
  xValues = unlist(attributes(perfa.validation)$x.values),
  yValues = unlist(attributes(perfa.validation)$y.values)
)
roc.optimal.accuracy <- max(roc.summary$accuracy)

# calculate confidence interval
roc.95ci <- roc(droplevels(results.classification.validation.rule.in$svm.summary$real.group),
                results.classification.validation.rule.in$svm.summary$asymptomaticControls, 
                ci = TRUE
)

pdf("ROC_pancancer.pdf")
plot(perfa.training, 
     lwd = 4, 
     col = "#C0C0C0", 
     lty=2, 
     ylim = c(0,1))
par(new=T)
plot(perfa.evaluation, 
     lwd = 6, 
     col = "#969696",
     ylim = c(0,1))
par(new=T)
plot(perfa.validation, 
     lwd = 6, 
     col = "#B03B3D",
     ylim = c(0,1)
)
dev.off()

## validation series
# pancancer ROC per tumor type
tumorTypes <- levels(results.classification.validation.rule.in$svm.summary$Group)
tumorTypes <- tumorTypes[which(!tumorTypes %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                  "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                  "Pancreatic diseases","Pulmonary Hypertension","Asymptomatic controls"))]
matrix <- matrix(NA, nrow = length(tumorTypes), ncol = 22)
rownames(matrix) <- tumorTypes
colnames(matrix) <- c('nTraining','nEvaluation','nValidation','AUCtraining','AUCevaluation','AUCvalidation','95AUCtrainingLower','95AUCtrainingUpper','95AUCevaluationLower','95AUCevaluationUpper',
                      '95AUCvalidationLower','95AUCvalidationUpper','correctStageI','correctStageII','correctStageIII','correctStageIV','correctStageNA','totalStageI','totalStageII','totalStageIII','totalStageIV','totalStageNA')
for(tumor in tumorTypes){
  if (!tumor %in% c("Esophageal carcinoma","Lymphoma")) {
    rocra.training <- prediction(as.numeric(as.character(results.classification.training$svm.summary[which(results.classification.training$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 5])), 
                                 results.classification.training$svm.summary[which(results.classification.training$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 3], 
                                 label.ordering = rev(levels((results.classification.training$svm.summary$real.group))))
    perfa.training <- performance(rocra.training, "tpr", "fpr")
    roc.95ci.tr <- roc(results.classification.training$svm.summary[which(results.classification.training$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$real.group,
                       results.classification.training$svm.summary[which(results.classification.training$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$asymptomaticControls, 
                       ci = TRUE
    )
    
    rocra.evaluation <- prediction(as.numeric(as.character(results.classification.evaluation.rule.in$svm.summary[which(results.classification.evaluation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4])), 
                                   results.classification.evaluation.rule.in$svm.summary[which(results.classification.evaluation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 3], 
                                   label.ordering = rev(levels(results.classification.evaluation.rule.in$svm.summary$real.group)))
    perfa.evaluation <- performance(rocra.evaluation, "tpr", "fpr")
    roc.95ci.ev <- roc(results.classification.evaluation.rule.in$svm.summary[which(results.classification.evaluation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$real.group,
                       results.classification.evaluation.rule.in$svm.summary[which(results.classification.evaluation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$asymptomaticControls, 
                       ci = TRUE
    )
    
    rocra.validation <- prediction(as.numeric(as.character(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4])), 
                                   results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 3], 
                                   label.ordering = rev(levels(droplevels(results.classification.validation.rule.in$svm.summary$real.group))))
    perfa.validation <- performance(rocra.validation, "tpr", "fpr")
    roc.95ci.val <- roc(droplevels(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$real.group),
                        results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$asymptomaticControls, 
                        ci = TRUE
    )
    
    save(perfa.validation, file = paste('perfa.validation-',tumor,".RData",sep=""), version = 2)
    
    pdf(paste("ROC-",tumor,".pdf",sep=""))
    plot(perfa.training, 
         lwd = 4, 
         col = "#C0C0C0", 
         lty=2)
    par(new=T)
    plot(perfa.evaluation, 
         lwd = 6, 
         col = "#969696")
    par(new=T)
    plot(perfa.validation, 
         lwd = 6, 
         col = "#B03B3D"
     )
    dev.off()
    
    matrix[which(rownames(matrix) == tumor),] <- c(
      length(as.numeric(as.character(results.classification.training$svm.summary[which(results.classification.training$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 5]))),
      length(as.numeric(as.character(results.classification.evaluation.rule.in$svm.summary[which(results.classification.evaluation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4]))),
      length(as.numeric(as.character(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4]))),
      attributes(performance(rocra.training, 'auc'))$y.values[[1]],
      attributes(performance(rocra.evaluation, 'auc'))$y.values[[1]],
      attributes(performance(rocra.validation, 'auc'))$y.values[[1]],
      roc.95ci.tr$ci[1],
      roc.95ci.tr$ci[3],
      roc.95ci.ev$ci[1],
      roc.95ci.ev$ci[3],
      roc.95ci.val$ci[1],
      roc.95ci.val$ci[3],
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","")),])
    )
    
  } else {
    rocra.validation <- prediction(as.numeric(as.character(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4])), 
                                   results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 3], 
                                   label.ordering = rev(levels(droplevels(results.classification.validation.rule.in$svm.summary$real.group))))
    perfa.validation <- performance(rocra.validation, "tpr", "fpr")
    roc.95ci.val <- roc(droplevels(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$real.group),
                        results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)),]$asymptomaticControls, 
                        ci = TRUE
    )
    
    save(perfa.validation, file = paste('perfa.validation-',tumor,".RData",sep=""), version = 2)
    
    pdf(paste("ROC-",tumor,".pdf",sep=""))
    plot(perfa.validation, 
         lwd = 6, 
         col = "#B03B3D"
    )
    dev.off()
    
    matrix[which(rownames(matrix) == tumor),] <- c(
      NA,
      NA,
      length(as.numeric(as.character(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls",tumor)), 4]))),
      NA,
      NA,
      attributes(performance(rocra.validation, 'auc'))$y.values[[1]],
      NA,
      NA,
      NA,
      NA,
      roc.95ci.val$ci[1],
      roc.95ci.val$ci[3],
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$predicted.group=="Cancer"),])/nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]),
      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group==tumor & results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","")),])
    )
  }
}
write.csv(matrix, file = "ROCtumorTypes.csv")

## validation series
# pancancer barplots per stage
Healthy <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),])
Healthy.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]$predicted.group)))[1], 
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]))

StageI <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),])
StageI.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]$predicted.group)))[2], 
                        nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I")),]))

StageII <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),])
StageII.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]$predicted.group)))[2], 
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("II")),]))

StageIII <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),])
StageIII.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]$predicted.group)))[2], 
                          nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III")),]))

StageI_II <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II")),])
StageI_II.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II")),]$predicted.group)))[2], 
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II")),]))

StageIII_IV <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III","IV","n.a.","n.i.") &
                                                                                                   results.classification.validation.rule.in$svm.summary$real.group == "Cancer"),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III","IV","n.a.","n.i.") &
                                                                     results.classification.validation.rule.in$svm.summary$real.group == "Cancer"),])
StageIII_IV.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III","IV","n.a.","n.i.") &
                                                                                                                 results.classification.validation.rule.in$svm.summary$real.group == "Cancer"),]$predicted.group)))[2], 
                             nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("III","IV","n.a.","n.i.") &
                                                                                                results.classification.validation.rule.in$svm.summary$real.group == "Cancer"),]))

StageI_III <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II","III")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II","III")),])
StageI_III.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II","III")),]$predicted.group)))[2], 
                            nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("I","II","III")),]))


StageIV <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),])
StageIV.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]$predicted.group)))[2],
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("IV")),]))

StageNA <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                                                                                                       "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                                                                                                       "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                                                                                                       "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                                                                                        "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                                                                                        "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                                                                                        "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),])
StageNA.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                                                                                                                                "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                                                                                                                                "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                                                                                                                                "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]$predicted.group)))[2],
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Stage %in% c("n.a.","n.i.","") & results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                                                                                                    "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                                                                                                    "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                                                                                                    "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]))

StageAll <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                            "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                            "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                            "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                 "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                 "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                 "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),])
StageAll.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                               "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                               "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                               "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]$predicted.group)))[2],
                          nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer","Endometrial cancer",
                                                                                                                                                                         "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                                                                                                                                         "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                                                                                                                                         "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer")),]))

pdf("Barplot-detectionPerStage.pdf")
par(mai   = c(1.5, 1, 0.8, 1))
barCenters <- barplot(height = c(StageI*100, StageII*100, StageIII*100, StageIV*100, StageNA*100, StageAll*100, Healthy*100),
                      beside = true, las = 2,
                      ylim = c(0, 100),
                      cex.names = 0.75, xaxt = "n",
                      border = "black", col= c(rep("#B03B3D", times = 6), "#1F497D") ,
                      axes = TRUE)
segments(barCenters, c(StageI.ci$conf.int[1]*100,
                       StageII.ci$conf.int[1]*100,
                       StageIII.ci$conf.int[1]*100,
                       StageIV.ci$conf.int[1]*100,
                       StageNA.ci$conf.int[1]*100,
                       StageAll.ci$conf.int[1]*100,
                       Healthy.ci$conf.int[1]*100)
         , barCenters,
         c(StageI.ci$conf.int[2]*100,
           StageII.ci$conf.int[2]*100,
           StageIII.ci$conf.int[2]*100,
           StageIV.ci$conf.int[2]*100,
           StageNA.ci$conf.int[2]*100,
           StageAll.ci$conf.int[2]*100,
           Healthy.ci$conf.int[2]*100), lwd = 3)
dev.off()

## validation series
# pancancer barplots per tumor type
# in manuscript presented as coxcombplot (created in MATLAB)
BrCa <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer")),])
BrCa.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer")),]$predicted.group)))[2], 
                      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Breast cancer")),]))
esophagus <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Esophageal carcinoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Esophageal carcinoma")),])
esophagus.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Esophageal carcinoma")),]$predicted.group)))[2],
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Esophageal carcinoma")),]))
PDAC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic cancer")),])
PDAC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic cancer")),]$predicted.group)))[2],
                      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic cancer")),]))
endometrial <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Endometrial cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Endometrial cancer")),])
endometrial.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Endometrial cancer")),]$predicted.group)))[2],
                             nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Endometrial cancer")),]))
CRC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Colorectal cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Colorectal cancer")),])
CRC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Colorectal cancer")),]$predicted.group)))[2],
                     nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Colorectal cancer")),]))
GBM <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Glioma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Glioma")),])
GBM.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Glioma")),]$predicted.group)))[2],
                     nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Glioma")),]))
HNC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Head and neck cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Head and neck cancer")),])
HNC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Head and neck cancer")),]$predicted.group)))[2],
                     nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Head and neck cancer")),]))
Melanoma <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Melanoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Melanoma")),])
Melanoma.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Melanoma")),]$predicted.group)))[2],
                          nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Melanoma")),]))
OvCa <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Ovarian cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Ovarian cancer")),])
OvCa.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Ovarian cancer")),]$predicted.group)))[2],
                      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Ovarian cancer")),]))
Cholangio <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Cholangiocarcinoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Cholangiocarcinoma")),])
Cholangio.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Cholangiocarcinoma")),]$predicted.group)))[2],
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Cholangiocarcinoma")),]))
RCC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Renal cell cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Renal cell cancer")),])
RCC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Renal cell cancer")),]$predicted.group)))[2],
                     nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Renal cell cancer")),]))
HL <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Lymphoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Lymphoma")),])
HL.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Lymphoma")),]$predicted.group)))[2],
                    nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Lymphoma")),]))
NSCLC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Non-small-cell lung cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Non-small-cell lung cancer")),])
NSCLC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Non-small-cell lung cancer")),]$predicted.group)))[2],
                       nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Non-small-cell lung cancer")),]))
Sarcoma <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Sarcoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Sarcoma")),])
Sarcoma.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Sarcoma")),]$predicted.group)))[2],
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Sarcoma")),]))
HCC <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hepatocellular carcinoma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hepatocellular carcinoma")),])
HCC.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hepatocellular carcinoma")),]$predicted.group)))[2],
                     nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hepatocellular carcinoma")),]))
UrothelCa <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Urothelial cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Urothelial cancer")),])
UrothelCa.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Urothelial cancer")),]$predicted.group)))[2],
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Urothelial cancer")),]))
MM <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple Myeloma")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple Myeloma")),])
MM.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple Myeloma")),]$predicted.group)))[2],
                    nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple Myeloma")),]))
PrCa <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Prostate cancer")),]$predicted.group)))[2] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Prostate cancer")),])
PrCa.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Prostate cancer")),]$predicted.group)))[2],
                      nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Prostate cancer")),]))

Healthy <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),])
Healthy.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]$predicted.group)))[1],
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Asymptomatic controls")),]))

Sympt <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                                              "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                                              "Pancreatic diseases","Pulmonary Hypertension")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                      "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                      "Pancreatic diseases","Pulmonary Hypertension")),])

pdf("Barplot-detectionPerGroup.pdf")
par(mai   = c(4, 1, 0.8, 0.1))
barCenters <- barplot(height = c(
  PDAC*100,BrCa*100,esophagus*100, endometrial*100,  CRC*100,
  GBM*100,Melanoma*100, HNC*100, OvCa*100, Cholangio*100, RCC*100,
  HL*100, NSCLC*100, Sarcoma*100, HCC*100, UrothelCa*100, MM*100, 
  PrCa*100, 0,Healthy*100
),
beside = true, las = 2,
ylim = c(0, 100),
cex.names = 0.75, # xaxt = "n",
border = "black", col="#387EB9", axes = TRUE,
names.arg = c("PDAC","BrCa","esophagus","endometrial","CRC",
              "GBM","Melanoma","HNC","OvCa","Cholangio","RCC",
              "HL","NSCLC","Sarcoma","HCC","UrothelCa",
              "MM","PrCa","","Healthy"))
segments(barCenters, c(PDAC.ci$conf.int[1]*100,
                       BrCa.ci$conf.int[1]*100,
                       esophagus.ci$conf.int[1]*100,
                       endometrial.ci$conf.int[1]*100,
                       CRC.ci$conf.int[1]*100,
                       GBM.ci$conf.int[1]*100,
                       Melanoma.ci$conf.int[1]*100,
                       HNC.ci$conf.int[1]*100,
                       OvCa.ci$conf.int[1]*100,
                       Cholangio.ci$conf.int[1]*100,
                       RCC.ci$conf.int[1]*100,
                       HL.ci$conf.int[1]*100,
                       NSCLC.ci$conf.int[1]*100,
                       Sarcoma.ci$conf.int[1]*100,
                       HCC.ci$conf.int[1]*100,
                       UrothelCa.ci$conf.int[1]*100,
                       MM.ci$conf.int[1]*100,
                       PrCa.ci$conf.int[1]*100,
                       0,
                       Healthy.ci$conf.int[1]*100)
         , barCenters,
         c(PDAC.ci$conf.int[2]*100,
           BrCa.ci$conf.int[2]*100,
           esophagus.ci$conf.int[2]*100,
           endometrial.ci$conf.int[2]*100,
           CRC.ci$conf.int[2]*100,
           GBM.ci$conf.int[2]*100,
           Melanoma.ci$conf.int[2]*100,
           HNC.ci$conf.int[2]*100,
           OvCa.ci$conf.int[2]*100,
           Cholangio.ci$conf.int[2]*100,
           RCC.ci$conf.int[2]*100,
           HL.ci$conf.int[2]*100,
           NSCLC.ci$conf.int[2]*100,
           Sarcoma.ci$conf.int[2]*100,
           HCC.ci$conf.int[2]*100,
           UrothelCa.ci$conf.int[2]*100,
           MM.ci$conf.int[2]*100,
           PrCa.ci$conf.int[2]*100,
           0,
           Healthy.ci$conf.int[2]*100), lwd = 3)
dev.off()

## validation series
# pancancer barplots of symptomatic non-cancer groups
load('results.classification.validation.rule.in.RData')
Hematuria <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hematuria")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hematuria")),])
Hematuria.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hematuria")),]$predicted.group)))[1],
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Hematuria")),]))

PulmHyper <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pulmonary Hypertension")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pulmonary Hypertension")),])
PulmHyper.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pulmonary Hypertension")),]$predicted.group)))[1],
                           nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pulmonary Hypertension")),]))


formerSarcoma <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Former sarcoma")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Former sarcoma")),])
formerSarcoma.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Former sarcoma")),]$predicted.group)))[1],
                               nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Former sarcoma")),]))

anginaP <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris")),])
anginaP.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris")),]$predicted.group)))[1],
                         nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris")),]))

epilepsy <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Medically-intractable epilepsy")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Medically-intractable epilepsy")),])
epilepsy.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Medically-intractable epilepsy")),]$predicted.group)))[1],
                          nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Medically-intractable epilepsy")),]))

pancreaticDis <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic diseases")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic diseases")),])
pancreaticDis.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic diseases")),]$predicted.group)))[1],
                               nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Pancreatic diseases")),]))

bowel <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Bowel disease")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Bowel disease")),])
bowel.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Bowel disease")),]$predicted.group)))[1],
                       nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Bowel disease")),]))

nSTEMI <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("nSTEMI")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("nSTEMI")),])
nSTEMI.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("nSTEMI")),]$predicted.group)))[1],
                        nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("nSTEMI")),]))

MS <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple sclerosis")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple sclerosis")),])
MS.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple sclerosis")),]$predicted.group)))[1],
                    nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Multiple sclerosis")),]))

Sympt <- summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                                              "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                                              "Pancreatic diseases","Pulmonary Hypertension")),]$predicted.group)))[1] / 
  nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                      "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                      "Pancreatic diseases","Pulmonary Hypertension")),])
Sympt.ci <- binom.test(summary(factor(unlist(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                                                            "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                                                            "Pancreatic diseases","Pulmonary Hypertension")),]$predicted.group)))[1],
                       nrow(results.classification.validation.rule.in$svm.summary[which(results.classification.validation.rule.in$svm.summary$Group %in% c("Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                                                                                                                           "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                                                                                                                           "Pancreatic diseases","Pulmonary Hypertension")),]))

pdf("Barplot-detectionPerGroupSymptomatics.pdf")
par(mai   = c(4, 1, 0.8, 0.1))
barCenters <- barplot(height = c(
  Hematuria*100,PulmHyper*100, formerSarcoma*100, anginaP*100, 
  pancreaticDis*100, bowel*100, epilepsy*100, nSTEMI*100, MS*100, 0, Sympt*100
),
beside = true, las = 2,
ylim = c(0, 100),
cex.names = 0.75,
border = "black", col= "#1F497D", #"#387EB9", 
axes = TRUE,
names.arg = c("Hematuria","PukmHyper","formerSarcoma","anginaP","epilepsy",
              "pancreaticDis","bowel","nSTEMI","MS","","Sympt"))
segments(barCenters, c(Hematuria.ci$conf.int[1]*100,
                       PulmHyper.ci$conf.int[1]*100,
                       formerSarcoma.ci$conf.int[1]*100,
                       anginaP.ci$conf.int[1]*100,
                       pancreaticDis.ci$conf.int[1]*100,
                       bowel.ci$conf.int[1]*100,
                       epilepsy.ci$conf.int[1]*100,
                       nSTEMI.ci$conf.int[1]*100,
                       MS.ci$conf.int[1]*100,
                       0,
                       Sympt.ci$conf.int[1]*100
)
, barCenters,
c(Hematuria.ci$conf.int[2]*100,
  PulmHyper.ci$conf.int[2]*100,
  formerSarcoma.ci$conf.int[2]*100,
  anginaP.ci$conf.int[2]*100,
  pancreaticDis.ci$conf.int[2]*100,
  bowel.ci$conf.int[2]*100,
  epilepsy.ci$conf.int[2]*100,
  nSTEMI.ci$conf.int[2]*100,
  MS.ci$conf.int[2]*100,
  0,
  Sympt.ci$conf.int[2]*100), lwd = 3)
dev.off()

########################
# TUMOR SITE OF ORIGIN #
########################

source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO_multiclass.R')
# fixed analyses parameters
nIter = 1
nCores = 50
percentage.for.training = 40
percentage.for.evaluation = 40
percentage.for.validation = 20
library(foreach)

# load and prepare the input data provided in GEO (GSE183635). Contains the raw read counts
# from all included samples in the TSOO analysis (n=1025), and all annotated counts (n=57736)
load("TEP_Count_Matrix_TSOO.RData")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('sampleInfo_libSize.csv', sep=",", row.names=1)
sampleInfo <- sampleInfo[which(rownames(sampleInfo) %in% colnames(TEP_Count_Matrix_TSOO)),]
sampleInfo <- sampleInfo[order(match(rownames(sampleInfo), colnames(TEP_Count_Matrix_TSOO))), ]

# load gene info and create edgeR object
load('bin/dgeGenesEnsembl75.RData')
library(edgeR)
dge <- DGEList(counts = TEP_Count_Matrix_TSOO,
               group = sampleInfo$Group,
               genes = genes
)
dge$samples <- cbind(dge$samples, sampleInfo) # add the sample info to the object
dge$samples$lib.size <- sampleInfo$lib.size

# rename the included five tumor types
dge$samples$group <- factor(dge$samples$group, levels = c(levels(dge$samples$group),"GLIO","HNSCC","NSCLC","PDAC","OVCAR"))
dge$samples$group[which(dge$samples$group %in% c("Glioma"))] <- "GLIO"
dge$samples$group[which(dge$samples$group %in% c("Head and neck cancer"))] <- "HNSCC"
dge$samples$group[which(dge$samples$group %in% c("Non-small-cell lung cancer"))] <- "NSCLC"
dge$samples$group[which(dge$samples$group %in% c("Ovarian cancer"))] <- "PDAC"
dge$samples$group[which(dge$samples$group %in% c("Pancreatic cancer"))] <- "OVCAR"
dge$samples <- droplevels(dge$samples)
summary(dge$samples$group)

# create new directory to separate both analyses
dir.create('TSOO/')
setwd('TSOO')

# randomly select the training-evaluation-validation series in a five-fold cross validation setting
# it requires to five times split the dataset into these series while prohibiting that samples
# are included twice in the validation series. The set.seed is locked to 1 to reproduce the 
# settings as presented by In 't Veld et al.
set.seed(nIter)
series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.training <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
  ) 
  
  training.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i]],
    size = n.samples.training,
    replace = F
  )
  
  # container
  series <- list()
  series[["training.samples.subset"]] <- training.samples.subset
  series
}
training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  

series.evaluation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.evaluation <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.evaluation / 100)
  ) 
  
  evaluation.samples.subset <- sample(
    colnames(dge[, dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% training.samples]),
    size = n.samples.evaluation,
    replace = F
  )
  
  # container
  series <- list()
  series[["evaluation.samples.subset"]] <- evaluation.samples.subset
  series
}
evaluation.samples <- unlist(lapply(series.evaluation,  function(x){x[["evaluation.samples.subset"]]}))  

training.samples.1 <- training.samples
evaluation.samples.1 <- evaluation.samples
validation.samples.1 <- colnames(dge[,which(!colnames(dge) %in% c(training.samples, evaluation.samples))])

set.seed(nIter)
series.validation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.validation <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.validation / 100)
  ) 
  
  validation.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% validation.samples.1],
    size = n.samples.validation,
    replace = F
  )
  
  # container
  series <- list()
  series[["validation.samples.subset"]] <- validation.samples.subset
  series
}
validation.samples <- unlist(lapply(series.validation, function(x){x[["validation.samples.subset"]]}))  

series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.training <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
  ) 
  
  training.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% validation.samples],
    size = n.samples.training,
    replace = F
  )
  
  # container
  series <- list()
  series[["training.samples.subset"]] <- training.samples.subset
  series
}
training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  

training.samples.2 <- training.samples
validation.samples.2 <- validation.samples
evaluation.samples.2 <- colnames(dge)[which(!colnames(dge) %in% c(training.samples.2, validation.samples.2))]

set.seed(nIter)
series.validation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.validation <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.validation / 100)
  ) 
  
  validation.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% c(validation.samples.1, validation.samples.2)],
    size = n.samples.validation,
    replace = F
  )
  
  # container
  series <- list()
  series[["validation.samples.subset"]] <- validation.samples.subset
  series
}
validation.samples <- unlist(lapply(series.validation, function(x){x[["validation.samples.subset"]]}))  

series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.training <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
  ) 
  
  training.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% validation.samples],
    size = n.samples.training,
    replace = F
  )
  
  # container
  series <- list()
  series[["training.samples.subset"]] <- training.samples.subset
  series
}
training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  

training.samples.3 <- training.samples
validation.samples.3 <- validation.samples
evaluation.samples.3 <- colnames(dge)[which(!colnames(dge) %in% c(training.samples.3, validation.samples.3))]

set.seed(nIter)
series.validation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.validation <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.validation / 100)
  )
  
  validation.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% c(validation.samples.1, validation.samples.2, validation.samples.3)],
    size = n.samples.validation,
    replace = F
  )
  
  # container
  series <- list()
  series[["validation.samples.subset"]] <- validation.samples.subset
  series
}
validation.samples <- unlist(lapply(series.validation, function(x){x[["validation.samples.subset"]]}))

series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.training <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
  )
  
  training.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% validation.samples],
    size = n.samples.training,
    replace = F
  )
  
  # container
  series <- list()
  series[["training.samples.subset"]] <- training.samples.subset
  series
}
training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))

training.samples.4 <- training.samples
validation.samples.4 <- validation.samples
evaluation.samples.4 <- colnames(dge)[which(!colnames(dge) %in% c(training.samples.4, validation.samples.4))]

set.seed(nIter)
validation.samples.5 <- colnames(dge)[which(!colnames(dge) %in% c(validation.samples.1, validation.samples.2, validation.samples.3, validation.samples.4))]

series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
  n.samples.training <- round(length(which(
    dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
  )
  
  training.samples.subset <- sample(
    colnames(dge)[dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% validation.samples.5],
    size = n.samples.training,
    replace = F
  )
  
  # container
  series <- list()
  series[["training.samples.subset"]] <- training.samples.subset
  series
}
training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))

training.samples.5 <- training.samples
evaluation.samples.5 <- colnames(dge)[which(!colnames(dge) %in% c(training.samples.5, validation.samples.5))]
evaluation.samples <- colnames(dge)[which(!colnames(dge) %in% c(training.samples.5, validation.samples.5))]

# Filter the dataset for low abundant RNAs.
dgeIncludedSamples <- filter.for.platelet.transcriptome.TrainEval.group(dge=dge, 
                                                                        minimum.read.counts = 30,
                                                                        minimum.prct.cohort = 90,
                                                                        training.series.only = TRUE,
                                                                        training.series.only.samples = c(training.samples, evaluation.samples),
                                                                        groupPlateletomeSelection = TRUE,
                                                                        verbose = TRUE)

# Perform thromboSeq ANOVA differential expression analysis of splice junctions
thromboSeq.anova <- thromboSeqANOVA(dge = dgeIncludedSamples[,c(training.samples, evaluation.samples)],
                                    k.variables = nlevels(dgeIncludedSamples$samples$group)+2, ## seems groups is highly represented
                                    variable.to.assess = c("lib.size"),
                                    variable.threshold = c(0.8),
                                    ruvg.pvalue.threshold.group = 1e-2,
                                    ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                    training.series.only = FALSE,
                                    select.biomarker.FDR = TRUE,
                                    plot = F,
                                    clinical.info.heatmap = c("group","isolationlocation","patientgroup"),
                                    swarm.optimization = F,
                                    n.particles = 60,
                                    n.iterations = 8,
                                    iteration = NULL,
                                    figureDir = "figureOutputFolder",
                                    number.cores = nCores,
                                    verbose = TRUE)

# Perform PSO-enhanced thromboSeq classifier development.
# For this, run the following command:
thromboPSO <- thromboSeqPSO(dge = dgeIncludedSamples, 
                            training.samples.provided = training.samples,
                            evaluation.samples.provided = evaluation.samples,
                            swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"),
                            swarm.boundaries = c(-0.1, 1.0, 50, as.numeric(summary(thromboSeq.anova$FDR < 1e-20)[3]), 0.5, 1.0, 50, as.numeric(summary(thromboSeq.anova$FDR < 1e-20)[3])),
                            k.variables = nlevels(dgeIncludedSamples$samples$group)+2,
                            variable.to.assess = c("lib.size"),
                            variable.threshold = c(0.8),
                            ruvg.pvalue.threshold.group = 1e-8,
                            ruvg.pvalue.threshold.strongest.variable = 1e-2,
                            select.biomarker.FDR = TRUE,
                            minimum.n.transcripts.biomarkerpanel = 2,
                            svm.gamma.range = 2^(-10:0),
                            svm.cost.range = 2^(0:10),
                            number.cross.splits = 2,
                            n.particles.gamma.cost.optimization = 50,
                            n.iterations.gamma.cost.optimization = 4,
                            n.particles = 60,
                            n.iterations = 8,
                            figureDir = "figureOutputFolder",
                            number.cores = nCores, 
                            verbose = TRUE
)

# Validate the developed thromboSeq algorithm.
thromboPSOreadout <- thromboSeqPSO.readout(dge = dgeIncludedSamples,
                                           replace.counts.validation = 0, 
                                           filter.clinical.characteristics.validation = NA,
                                           filter.clinical.characteristics.group = NA,
                                           filter.clinical.characteristics.specified = NA,
                                           readout.training = F, 
                                           readout.evaluation = T, 
                                           readout.validation = T, 
                                           apply.rule.in.readout = F, 
                                           rule.in.setting = 99, 
                                           apply.rule.out.readout = F, 
                                           rule.out.setting = NA,
                                           additional.dge.to.predict = NA,
                                           number.cores = nCores, 
                                           clinical.data.in.output = c('group','Sample.ID','Sample.supplying.institution','Stage','Sex')
)

# create figures from the output data
load('results.classification.evaluation.RData')
load('results.classification.validation.RData')

# count first and second algorithm choice.
results.classification.validation$svm.summary <- cbind(results.classification.validation$svm.summary, matrix(NA, nrow = nrow(results.classification.validation$svm.summary), ncol=6))
results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Sex=="F"),c('prostateCancer')] <- 0
results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Sex %in% c("M")),c('breastCancer','ovarianCancer','endometrialCancer')] <- 0
samples <- rownames(results.classification.validation$svm.summary)
nCol <- which(colnames(results.classification.validation$svm.summary)=='group')-1
for(sample in samples){
  results.classification.validation$svm.summary[sample,'1'] <- colnames(results.classification.validation$svm.summary[sample,c(4:nCol),])[which(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),]) == max(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),])))]
  results.classification.validation$svm.summary[sample,colnames(results.classification.validation$svm.summary[sample,c(4:nCol),])[which(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),]) == max(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),])))]] <- 0
  
  results.classification.validation$svm.summary[sample,'2'] <- colnames(results.classification.validation$svm.summary[sample,c(4:nCol),])[which(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),]) == max(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),])))]
  results.classification.validation$svm.summary[sample,colnames(results.classification.validation$svm.summary[sample,c(4:nCol),])[which(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),]) == max(as.numeric(results.classification.validation$svm.summary[sample,c(4:nCol),])))]] <- 0
  
  if (results.classification.validation$svm.summary[sample,'real.group'] == results.classification.validation$svm.summary[sample,'1']){
    results.classification.validation$svm.summary[sample,'4'] <- "1st correct"
    results.classification.validation$svm.summary[sample,'6'] <- "1st/2nd"
  }
  if (results.classification.validation$svm.summary[sample,'real.group'] == results.classification.validation$svm.summary[sample,'2']){
    results.classification.validation$svm.summary[sample,'5'] <- "2nd correct"
    results.classification.validation$svm.summary[sample,'6'] <- "1st/2nd"
  }
}

# total first correct
summary(factor(results.classification.validation$svm.summary$`4`))
# total second correct
summary(factor(results.classification.validation$svm.summary$`5`))
# total 1st/2nd correct
summary(factor(results.classification.validation$svm.summary$`6`))
# total classifications (ie size of validation series)
nrow(results.classification.validation$svm.summary)

# summarize accuracy per stage for only first classification
# accuracy stage I-III
stageI_III <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]$`4`=="1st correct")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),])
stageI_III.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]$`4`=="1st correct")[2]),
                            nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]))

# accuracy stage IV
stageIV <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]$`4`=="1st correct")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),])
stageIV.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]$`4`=="1st correct")[2]),
                         nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]))

# accuracy stage NA
stageNA <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]$`4`=="1st correct")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),])
stageNA.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]$`4`=="1st correct")[2]),
                         nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]))

pdf("Barplot-detectionPerStage_1stprediction.pdf")
par(mai   = c(4, 1, 0.8, 5))
barCenters <- barplot(height = c(stageI_III*100, stageIV*100, stageNA*100),
                      beside = true, las = 2,
                      ylim = c(0, 100),
                      cex.names = 0.75, xaxt = "n",
                      border = "black", col=c("#fb6a4a","#de2d26","#a50f15"), axes = TRUE)
segments(barCenters, c(stageI_III.ci$conf.int[1]*100,
                       stageIV.ci$conf.int[1]*100,
                       stageNA.ci$conf.int[1]*100
)
, barCenters,
c(stageI_III.ci$conf.int[2]*100,
  stageIV.ci$conf.int[2]*100,
  stageNA.ci$conf.int[2]*100
), lwd = 3)
dev.off()

group = "GLIO"
Glioma_1 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                     & results.classification.validation$svm.summary$`1` == group),])
Glioma_2 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                     & results.classification.validation$svm.summary$`2` == group),])
Glioma <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group),])

group = "HNSCC"
headAndNeck_1 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                          & results.classification.validation$svm.summary$`1` == group),])
headAndNeck_2 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                          & results.classification.validation$svm.summary$`2` == group),])
headAndNeck <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group),])

group = "NSCLC"
NSCLC_1 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                    & results.classification.validation$svm.summary$`1` == group),])
NSCLC_2 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                    & results.classification.validation$svm.summary$`2` == group),])
NSCLC <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group),])

group = "OVCAR"
ovarianCancer_1 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                            & results.classification.validation$svm.summary$`1` == group),])
ovarianCancer_2 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                            & results.classification.validation$svm.summary$`2` == group),])
ovarianCancer <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group),])

group = "PDAC"
PDAC_1 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                   & results.classification.validation$svm.summary$`1` == group),])
PDAC_2 <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group
                                                                   & results.classification.validation$svm.summary$`2` == group),])
PDAC <- nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$real.group == group),])

data <- matrix(NA, nrow = 1, ncol = 5)
data[,1] <- c(round((NSCLC_1/(NSCLC))*100, digits = 0))
data[,2] <- c(round((Glioma_1/(Glioma))*100, digits = 0))
data[,3] <- c(round((headAndNeck_1/(headAndNeck))*100, digits = 0))
data[,4] <- c(round((ovarianCancer_1/(ovarianCancer))*100, digits = 0))
data[,5] <- c(round((PDAC_1/(PDAC))*100, digits = 0))

NSCLC.bin <- binom.test(NSCLC_1, NSCLC)
Glioma.bin <- binom.test(Glioma_1, Glioma)
PDAC.bin <- binom.test(PDAC_1, PDAC)
HNC.bin <- binom.test(headAndNeck_1, headAndNeck)
OvCa.bin <- binom.test(ovarianCancer_1, ovarianCancer)

pdf('StackedBarplotTOO.pdf')
par(mai   = c(3, 1, 0.8, 2))
barPlot <- barplot(data,
                   col = c("#d7191c","#2b83ba"),
                   border="white", 
                   space=0.04, 
                   font.axis=2, 
                   ylim  = c(0,100))
segments(barPlot, 
         c(NSCLC.bin$conf.int[1]*100,
           Glioma.bin$conf.int[1]*100,
           HNC.bin$conf.int[1]*100,
           OvCa.bin$conf.int[1]*100,
           PDAC.bin$conf.int[1]*100
         )
         , barPlot,
         c(NSCLC.bin$conf.int[2]*100,
           Glioma.bin$conf.int[2]*100,
           HNC.bin$conf.int[2]*100,
           OvCa.bin$conf.int[2]*100,
           PDAC.bin$conf.int[2]*100
         ), 
         lwd = 3)
dev.off()

# create confusion matrix of only first classifications counted
library(reshape)
confusion.matrix <- as.data.frame(
  cast(results.classification.validation$svm.summary,
       1 ~ real.group,
       length,
       value = "sampleName"))
rownames(confusion.matrix) <- confusion.matrix$`1`
confusion.matrix <- confusion.matrix[,-1]
lev <- sort(unique(c(colnames(confusion.matrix), rownames(confusion.matrix))))
confusion.matrix <- confusion.matrix[lev, lev]
colnames(confusion.matrix) <- lev
rownames(confusion.matrix) <- lev
confusion.matrix <- as.matrix(confusion.matrix)
confusion.matrix[is.na(confusion.matrix)] <- 0
melted.confusion.matrix <- as.data.frame(confusion.matrix)
melted.confusion.matrix$predicted <- rownames(melted.confusion.matrix)
melted.confusion.matrix <- melt(melted.confusion.matrix, id.vars = "predicted")
colnames(melted.confusion.matrix) <- c("predicted", "real", "frequency")
library(ggplot2)
tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
  geom_tile(aes(fill = frequency)) +
  scale_fill_continuous(low = "white", high = "red") +
  geom_text(aes(label = frequency)) +
  ggtitle("SVM classification results") +
  labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
  theme_bw() +
  theme(legend.position = "top")

pdf("ConfusionMatrix.pdf", paper = "a4")
print(tiles)
dev.off()

# summarize accuracy per stage for first and second classification
data <- matrix(NA, nrow = 2, ncol = 5)
data[,1] <- c(round((NSCLC_1/(NSCLC))*100, digits = 0), round((NSCLC_2/(NSCLC))*100, digits = 0))
data[,2] <- c(round((Glioma_1/(Glioma))*100, digits = 0), round((Glioma_2/(Glioma))*100, digits = 0))
data[,3] <- c(round((headAndNeck_1/(headAndNeck))*100, digits = 0), round((headAndNeck_2/(headAndNeck))*100, digits = 0))
data[,4] <- c(round((PDAC_1/(PDAC))*100, digits = 0), round((PDAC_2/(PDAC))*100, digits = 0)) 
data[,5] <- c(round((ovarianCancer_1/(ovarianCancer))*100, digits = 0), round((ovarianCancer_2/(ovarianCancer))*100, digits = 0))

NSCLC.bin <- binom.test(NSCLC_1+NSCLC_2, NSCLC)
Glioma.bin <- binom.test(Glioma_1+Glioma_2, Glioma)
PDAC.bin <- binom.test(PDAC_1+PDAC_2, PDAC)
HNC.bin <- binom.test(headAndNeck_1+headAndNeck_2, headAndNeck)
OvCa.bin <- binom.test(ovarianCancer_1+ovarianCancer_2, ovarianCancer)

pdf('StackedBarplotTOO_1st2nd.pdf')
par(mai   = c(3, 1, 0.8, 2))
barPlot <- barplot(data,
                   col = c("#99000d","#ef3b2c"),
                   border="white", 
                   space=0.04, 
                   font.axis=2, 
                   ylim  = c(0,100))
segments(barPlot, 
         c(NSCLC.bin$conf.int[1]*100,
           Glioma.bin$conf.int[1]*100,
           HNC.bin$conf.int[1]*100,
           PDAC.bin$conf.int[1]*100,
           OvCa.bin$conf.int[1]*100)
         , barPlot,
         c(NSCLC.bin$conf.int[2]*100,
           Glioma.bin$conf.int[2]*100,
           HNC.bin$conf.int[2]*100,
           PDAC.bin$conf.int[2]*100,
           OvCa.bin$conf.int[2]*100
         ), 
         lwd = 3)
dev.off()

# accuracy stage I-III, NA
stageI_III <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]$`6`=="1st/2nd")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),])
stageI_III.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]$`6`=="1st/2nd")[2]),
                            nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("I","II","III")),]))

# accuracy stage IV
stageIV <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]$`6`=="1st/2nd")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),])
stageIV.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]$`6`=="1st/2nd")[2]),
                         nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("IV")),]))

# accuracy stage NA
stageNA <- as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]$`6`=="1st/2nd")[2])/
  nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),])
stageNA.ci <- binom.test(as.numeric(summary(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]$`6`=="1st/2nd")[2]),
                         nrow(results.classification.validation$svm.summary[which(results.classification.validation$svm.summary$Stage %in% c("n.a.","n.i.")),]))

pdf("Barplot-detectionPerStage_1st2ndprediction.pdf")
par(mai = c(2, 3, 1, 2))
barCenters <- barplot(height = c(stageI_III*100, stageIV*100, stageNA*100),
                      beside = true, las = 2,
                      ylim = c(0, 100),
                      cex.names = 0.75, xaxt = "n",
                      border = "black", col=c("#fb6a4a","#de2d26","#a50f15"), axes = TRUE)
segments(barCenters, c(stageI_III.ci$conf.int[1]*100,
                       stageIV.ci$conf.int[1]*100,
                       stageNA.ci$conf.int[1]*100
)
, barCenters,
c(stageI_III.ci$conf.int[2]*100,
  stageIV.ci$conf.int[2]*100,
  stageNA.ci$conf.int[2]*100
), lwd = 3)
dev.off()

# end of script


