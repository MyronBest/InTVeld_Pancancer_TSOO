# R-tools for thromboSeq
# Functions included for PSO-enhanced classification algorithm development
# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Date          : 1st of September 2018
# Revision      : 19th of November 2018
# Revision v1.1 : 28th of November 2018
# Revision v1.2 : 29th of December 2018
# Revision v1.3 : 20th of March 2019
# Revision v1.4 : 29th of August 2019

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thromboSeqPSO ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thromboSeqPSO <- function(dge = dgeIncludedSamples,
                          percentage.for.training = 40,
                          percentage.for.evaluation = 30,
                          training.samples.provided = NA,
                          evaluation.samples.provided = NA,
                          swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"),
                          swarm.boundaries = c(-0.1, 1.0, 0.1, 1.0, 0.5, 1.0, 200, "all"),
                          k.variables = 3,
                          variable.to.assess = c("Age","lib.size"),
                          variable.threshold = c(0.2, 0.8),
                          ruvg.pvalue.threshold.group = 1e-2,
                          ruvg.pvalue.threshold.strongest.variable = 1e-2,
                          select.biomarker.FDR = FALSE,
                          minimum.n.transcripts.biomarkerpanel = 2,
                          svm.gamma.range = 2^(-20:0),
                          svm.cost.range = 2^(0:20),
                          number.cross.splits = 2,
                          n.particles.gamma.cost.optimization = 50,
                          n.iterations.gamma.cost.optimization = 4,
                          class.weights = TRUE,
                          rule.in.optimization = FALSE,
                          n.particles = 100,
                          n.iterations = 10,
                          figureDir = "figureOutputFolder",
                          number.cores = 2, 
                          verbose = TRUE){
  # Perform PSO-enhanced thromboSeq classification algorithm development. For this, first
  # select the training and evaluation series. Next, perform PSO-optimization. Finally, summarize data
  # and output the PSO progression plot.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   percentage.for.training: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the training series.
  #   percentage.for.evaluation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the evaluation series.
  #   training.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           training series.
  #   evaluation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           evaluation series.
  #   swarm.parameters: Vector with parameters to be optimized by swarm intelligence. By default the parameters FDR,
  #                     stable genes based on lib size correlation, correlated genes, and ranked genes are included.
  #                     Additional clinical parameters may be added by adding the sample info column names to this vector
  #                     plus the 'variable.to.assess'-vector and a default value to the 'variable.threshold'-vector.
  #   swarm.boundaries: Vector with lower and upper boundaries per variable in swarm parameters, separated
  #                     by a comma. Number of input values should match with the number of swarm parameters provided.
  #   k.variables: Number of (expected) variables/axes in the dataset.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
  #                                be reached by the correlation between the counts and  any variable
  #                                in order to bypass the wanted variable 'group' to be selected.
  #   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
  #                                be reached by the correlation between the counts and the specific 
  #                                variable in order to this variable to be assigned to the RUVg axis.
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
  #                                         biomarker panel.
  #   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
  #   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
  #   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   class.weights: TRUE/FALSE whether class weights should be included into the tune.svm process.
  #   rule.in.optimization: TRUE/FALSE whether or not to optimize the training process towards 99% specificity.
  #   n.particles: Numeric-value with number of PSO particles per iteration for classifier development.
  #   n.iterations: Numeric-value with number of PSO iterations in total for classifier development.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Returns the best settings selected by PSO.
  #  Also produces output files of the training and optimization process, including the ppso.log, ppso.log, and individual RData
  #  files in the folder outputPSO that contain the trained support vectors.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!percentage.for.training %in% seq(1,100,by=1e-3)){
    stop("percentage.for.training should be within 1 and 100")
  }
  
  if (!percentage.for.evaluation %in% seq(1,100,by=1e-3)){
    stop("percentage.for.evaluation should be within 1 and 100")
  }
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) > 1){
    print("Both percentage for training and a separate training list provided. The provided training list will be used.")
  }
  
  if (length(training.samples.provided) > 1 & length(evaluation.samples.provided) == 1){
    stop("list of training samples provided but no evaluation samples specified. Please specify evaluation samples.")
  }
  
  if (length(swarm.parameters)*2 != length(swarm.boundaries)){
    stop("number of swarm.parameters and provided swarm.boundaries do not match.")
  }
  
  if (!is.numeric(k.variables)){
    stop("Provide numeric value for k.variables")
  }
  
  if (k.variables < 1){
    stop("k.variables should be 1 or more")
  }
  
  if (all(variable.to.assess %in% colnames(dge$samples)) != TRUE){
    stop("Inputted variables do not match column names of the sample info table.")
  }
  
  if (length(variable.to.assess) != length(variable.threshold)){
    stop("Number of variables in variable.to.assess should match number of thresholds")
  }
  
  if (!is.numeric(ruvg.pvalue.threshold.group)){
    stop("Provide numeric value for ruvg.pvalue.threshold.group")
  }
  
  if (!is.numeric(ruvg.pvalue.threshold.strongest.variable)){
    stop("Provide numeric value for ruvg.pvalue.threshold.strongest.variable")
  }
  
  if (!is.numeric(minimum.n.transcripts.biomarkerpanel)){
    stop("Provide numeric value for minimum.n.transcripts.biomarkerpanel")
  }
  
  if (!is.numeric(number.cross.splits)){
    stop("Provide numeric value for number.cross.splits")
  }
  
  if (!is.numeric(n.particles.gamma.cost.optimization)){
    stop("Provide numeric value for n.particles.gamma.cost.optimization")
  }
  
  if (!is.numeric(n.iterations.gamma.cost.optimization)){
    stop("Provide numeric value for n.iterations.gamma.cost.optimization")
  }
  
  if (!is.numeric(n.particles)){
    stop("Provide numeric value for n.particles")
  }
  
  if (!is.numeric(n.iterations)){
    stop("Provide numeric value for n.iterations")
  }
  
  library(parallel)
  if (number.cores > detectCores()){
    print(paste("Please note that more computer threads have been selected than available with this machine. This may cause problems with start of the PSO process. It is advised to reduce the number of cores to ", detectCores(), sep = ""))
  }
    
  # load required packages
  if (verbose == TRUE){
    print("Load required packages ppso, edgeR, e1071, and RUVSeq")
  } 
  suppressMessages(library(ppso))
  suppressMessages(library(edgeR))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(e1071))
  suppressMessages(library(foreach))
  
  # create subdirectory for PSO-process
  if (!file.exists("pso-enhanced-thromboSeq")){
    dir.create("pso-enhanced-thromboSeq", recursive = T)
  }
  # store current directory and subdirectory
  workDir_main <- getwd()
  setwd("pso-enhanced-thromboSeq/") # transfer to subdirectory
  workDir <- getwd()
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) == 1){
  # randomly select samples for the training and evaluation series
  # here it is assumed the group size and potential confounding factors 
  # (e.g. age of the individuals) are similar among both groups
  set.seed(1000) # lock random number generator
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
  write.csv(training.samples, "trainingSamples.csv")
  
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
  write.csv(evaluation.samples,"evaluationSamples.csv")
  } else {
    training.samples <- training.samples.provided
    evaluation.samples <- evaluation.samples.provided
  }
  
  if (verbose == TRUE){
    print("Samples selected for training series:")
    print(summary((dge[, training.samples])$samples$group))
    print("Samples selected for evaluation series:")
    print(summary((dge[, evaluation.samples])$samples$group))
  } 
  
  # store these particular variables
  k.variables = k.variables
  variable.to.assess = variable.to.assess
  ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group
  ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable
  class.weights = class.weights
  rule.in.optimization = rule.in.optimization
  
  # input parameter boundaries for ppso:
  if ('all' %in% swarm.boundaries){
    # in case rankedTranscripts supplied with 'all', replace by numeric value of all RNAs detected
    swarm.boundaries[swarm.boundaries == 'all'] <- nrow(dge$counts)
  }
  
  # prepare matrix with per variable to optimize the lower and upper boundary value
  parameter.bounds <- matrix(ncol = 2, nrow = length(swarm.parameters))
  rownames(parameter.bounds) <- swarm.parameters
  parameter.bounds[,1] <- swarm.boundaries[seq(1, length(swarm.boundaries), by = 2)]
  parameter.bounds[,2] <- swarm.boundaries[seq(2, length(swarm.boundaries), by = 2)]
  parameter.bounds <- apply(parameter.bounds, 2, as.numeric) # ensure the matrix contains numeric values
  
  # create an RData file with all data in this R session, 
  # employed by ppso for analyzing all individual particles
  save(list = ls(envir = environment(), all.names = TRUE), 
       file = "Pre-PSO-snapshot.RData", envir = environment()) 
  if (!file.exists("outputPSO")){
    dir.create("outputPSO", recursive = T)
  } # separate output folder for individual particles
  
  # run ppso
  # Of note, this ppso-function may require days to complete. The algorithm will read each 
  # second whether the 'slaves' have produced new results and updates automatically the swarming process.
  # Continuous output can be seen in the 'outputPSO'-folder, the individual slave-log-files, 
  # and following each iteration of n particles the ppso.log-file.
  set.seed(1000) # lock random number generator
  resultPSO <- optim_ppso_robust(objective_function        = thrombo.algo,
                                 nslaves                   = number.cores - 1,
                                 number_of_parameters      = nrow(parameter.bounds),
                                 plot_progress             = FALSE,
                                 number_of_particles       = n.particles,
                                 max_number_of_iterations  = n.iterations,
                                 max_number_function_calls = n.particles * n.iterations,
                                 parameter_bounds          = parameter.bounds,
                                 tryCall                   = TRUE, 
                                 verbose                   = if (verbose == TRUE){TRUE} else {FALSE},
                                 lhc_init                  = TRUE, 
                                 wait_complete_iteration   = TRUE,
                                 logfile                   = paste(workDir, "/ppso.log", sep = ""),
                                 projectfile               = paste(workDir, "/ppso.pro", sep = ""),
                                 break_file                = "stopPPSO.txt"
  )
  # end of ppso
  
  # Select the best parameter setting
  logged.PSO.distribution <- read.csv("ppso.log", sep = "\t")
  logged.PSO.distribution <- logged.PSO.distribution[
    order(logged.PSO.distribution$objective_function), 
    ]
  logged.PSO.distribution.Index <- logged.PSO.distribution[
    which(logged.PSO.distribution$objective_function == min(logged.PSO.distribution$objective_function)), 
    ]
  
  set.seed(1000)
  # in case more particles have the same AUC output value, select randomly one as the particle for readout
  if (nrow(logged.PSO.distribution.Index) > 1){
    logged.PSO.distribution.Index <- logged.PSO.distribution.Index[
      sample(1 : nrow(logged.PSO.distribution.Index), size = 1),
      ]
  }
  
  # plot the pso-progress
  # add swarm.parameters to columns of ppso-overview log-file
  colnames(logged.PSO.distribution) <- c("time",swarm.parameters,"objective_function","worker")
  colnames(logged.PSO.distribution.Index) <- c("time",swarm.parameters,"objective_function","worker")
  # for each PSO-optimization parameter plot the input value to 1-AUC-result
  # highlight the ultimate best particle in red
  # generate and store plot in output folder
  if (!file.exists(paste(workDir_main,"/",figureDir,sep=""))){
    dir.create(paste(workDir_main,"/",figureDir,sep=""), recursive = T)
  }
  
  pdf(paste(workDir_main,"/",figureDir,"/PSO_optimizationPlots.pdf",sep=""), paper="a4")
  par(mfrow=c(2,2))
  for (parameter in swarm.parameters){
    plot(logged.PSO.distribution[,parameter],logged.PSO.distribution$objective_function,
         xlim = c(min(logged.PSO.distribution[,parameter]), max(logged.PSO.distribution[,parameter])), 
         ylim = if (max(logged.PSO.distribution$objective_function) != Inf) {
           c(0, max(logged.PSO.distribution$objective_function)) } else {
             c(0, 1)
           },
         pch = 20, 
         ylab = "1-AUC", 
         xlab = parameter)
    points(logged.PSO.distribution.Index[,parameter], logged.PSO.distribution.Index$objective_function, 
           col = "red", 
           pch = 19)
  }
  dev.off()
  
  # return the string with the best settings selected by PSO.
  return(logged.PSO.distribution.Index)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo <- function(x){
  # thromboSeq classification algorithm to be optimized by PSO. In brief, the classification 
  # algorithm performs data correction and normalization employing only the training series. 
  # Next, it performs ANOVA differential expression analysis of splice junctions. Then, it 
  # identifies highly correlated RNAs and removes those according to the PSO-proposed threshold. 
  # Following, it trains the first SVM algorithm of which the most contributing RNAs are 
  # maintained, as suggested by PSO. Then, an updated SVM-model is trained of which the gamma 
  # and cost setting are also PSO-optimized. Finally, the classification model is stored and the 
  # evaluation series are classified.
  # 
  # Args
  #   x: Vector with values proposed by PSO (ppso-function)
  #
  # Returns:
  #   Inverse AUC-value of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings
  
  # load R environment data
  load(paste(getwd(), "/Pre-PSO-snapshot.RData", sep = ""))
  # load packages
  suppressMessages(library(edgeR, warn.conflicts = F, quietly = T))
  suppressMessages(library(e1071, warn.conflicts = F, quietly = T))
  suppressMessages(library(foreach, warn.conflicts = F, quietly = T))
  suppressMessages(library(ROCR, warn.conflicts = F, quietly = T))
  suppressMessages(library(RUVSeq, warn.conflicts = F, quietly = T))
  suppressMessages(library(caret, warn.conflicts = F, quietly = T))
  suppressMessages(library(ppso, warn.conflicts = F, quietly = T))
  # load functions
  source(paste(workDir_main, "/bin/thromboSeqTools_PreProcessing_2.R", sep = ""))
  source(paste(workDir_main, "/bin/thromboSeqTools_ANOVA.R", sep = ""))
  source(paste(workDir_main, "/bin/thromboSeqTools_PSO.R", sep = ""))
  
  # prepare data according to training and evaluation series
  dgeTraining <- dge[, training.samples]
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  
  # ruv.correction
  # if variable.to.assess in PSO optimization, replace default value by PSO-proposed value
  for (variable in variable.to.assess){
    if (variable %in% swarm.parameters){
      variable.threshold[which(variable.to.assess==variable)] <- x[which(variable==swarm.parameters)]
    }
  }
  
  dgeTraining <- perform.RUVg.correction(dge = dge[, c(training.samples, evaluation.samples)], 
                                         k.variables = k.variables, 
                                         variable.to.assess = variable.to.assess,
                                         variable.threshold = variable.threshold, 
                                         ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                         ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                         training.series.only = T,
                                         training.series.only.samples = training.samples,
                                         verbose = verbose)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$raw.lib.size <- dgeTraining$samples$lib.size
  dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
  
  # perform TMM normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
  
  # Likelihood-ratio test modified for thromboSeq (ANOVA)
  dgeTraining$ruv.counts <- dgeTraining$counts
  dgeTraining$counts <- dgeTraining$raw.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$raw.lib.size
  thromboSeq.anova <- anovaLikeEdgeRthromboSeq(dgeTraining[, training.samples], 
                                               method = "TMM",
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.column = dgeTraining$refSample)
  
  # select FDR threshold. When FDR/fdr was included as a PSO optimization variable, select PSO-proposed value
  if (any(c("fdr","FDR") %in% swarm.parameters)){
    fdr <- x[which(swarm.parameters %in% c("fdr","FDR"))]
  } else {
    fdr <- nrow(thromboSeq.anova[which(thromboSeq.anova$FDR < 0.05),]) # in case fdr/FDR was not provided by user
  }
  
  # select biomarker panel using either FDR or p-value statistics
  if (select.biomarker.FDR == TRUE){
    selected.transcripts <- rownames(thromboSeq.anova)[
      thromboSeq.anova$logCPM > 3 & 
        thromboSeq.anova$chromosome_name %in% c(1:22, "X")
      ][c(0:fdr)][!is.na(rownames(thromboSeq.anova)[
        thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22, "X")
        ][c(0:fdr)])]
  } else {
    selected.transcripts <- rownames(thromboSeq.anova)[
      thromboSeq.anova$logCPM > 3 & 
        thromboSeq.anova$chromosome_name %in% c(1:22, "X")
      ][c(0:fdr)][!is.na(rownames(thromboSeq.anova)[
        thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22, "X")
        ][c(0:fdr)])]
  }
  
  # only continue when biomarker panel size is more than provided threshold
  if (length(selected.transcripts) > minimum.n.transcripts.biomarkerpanel){
    # remove transcripts with a Pearson's correlation to any other transcripts over the PSO-proposed threshold
    correlation.matrix <- cor(t(normalized.counts[selected.transcripts, training.samples]))
    if ("correlatedTranscripts" %in% swarm.parameters){
      correlated.transcripts <- x[which(swarm.parameters == "correlatedTranscripts")]
    } else {
      correlated.transcripts <- 1.0 # in case this variable was not provided by user
    }
    highly.correlated <- colnames(correlation.matrix)[findCorrelation(correlation.matrix, cutoff = correlated.transcripts)]
    
    # select biomarker panel using either FDR or p-value statistics
    selected.transcripts <- selected.transcripts[which(!selected.transcripts %in% highly.correlated)]
    
    # only continue when biomarker panel size is more than provided threshold
    if (length(selected.transcripts) > minimum.n.transcripts.biomarkerpanel){
      # first SVM-model, with grid search for gamma and cost
      if(class.weights == TRUE) {
      tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                            y           = real.groups.training,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            class.weights = 1/summary(dgeTraining$samples$group),
                            probability = TRUE
      ) } else {
        tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                              y           = real.groups.training,
                              gamma       = svm.gamma.range,
                              cost        = svm.cost.range,
                              tunecontrol = tune.control(cross = number.cross.splits),
                              probability = TRUE
        ) 
      }
      
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      # rank transcripts (features) of the biomarker panel according to their relative contribution to the SVM model
      if (nlevels(dgeTraining$samples$group) == 2){
        # binary feature ranking algorithm
        svm.ranking <- svmrfeFeatureRanking(
          t(normalized.counts[selected.transcripts, training.samples]), 
          real.groups.training, tuned.svm.model$gamma, tuned.svm.model$cost)
      } else {
        # multiclass feature ranking algorithm
        svm.ranking <- svmrfeFeatureRankingForMulticlass(
          t(normalized.counts[selected.transcripts, training.samples]), 
          real.groups.training, tuned.svm.model$gamma, tuned.svm.model$cost)
      }
      
      # employ PSO-proposed threshold and set final spliced RNA biomarker panel of this PSO particle
      if ("rankedTranscripts" %in% swarm.parameters){
        ranked.transcripts <- round(x[which(swarm.parameters == "rankedTranscripts")], digits = 0)
      } else {
        ranked.transcripts <- nrow(dgeTraining$counts) # in case rankedTranscripts was not provided by user
      }
      selected.transcripts <- (cbind(colnames(t(normalized.counts[selected.transcripts, training.samples])), svm.ranking)[
        which(as.numeric(cbind(colnames(t(normalized.counts[selected.transcripts, training.samples])), svm.ranking)[,2]) <= as.numeric(ranked.transcripts)),
        ])[,1]
      
      # second SVM-model, re-train with adjusted biomarker panel and employ a grid search for gamma and cost
      if(class.weights == TRUE) {
        tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                              y           = real.groups.training,
                              gamma       = svm.gamma.range,
                              cost        = svm.cost.range,
                              tunecontrol = tune.control(cross = number.cross.splits),
                              class.weights = 1/summary(dgeTraining$samples$group),
                              probability = TRUE
        ) } else {
          tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                                y           = real.groups.training,
                                gamma       = svm.gamma.range,
                                cost        = svm.cost.range,
                                tunecontrol = tune.control(cross = number.cross.splits),
                                probability = TRUE
          ) 
        }
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      ## optimize gamma and cost by a second PSO-optimization algorithm
      # employ the training series for SVM algorithm training, and the evaluation series for optimization
      # make sure the values for cost and gamma are not infinite
      if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
        # if necessary - broaden the range of cost and gamma
        if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
          svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
        }
        if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
          svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
        }
        svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                             svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
        svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                            svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])

        # input parameters:
        parameter.bounds.gamma.cost <- matrix(ncol = 2, nrow = 2)
        rownames(parameter.bounds.gamma.cost) <- c("gamma", "cost")
        parameter.bounds.gamma.cost[, 1] <- c(svm.gamma.range[1], svm.cost.range[1])
        parameter.bounds.gamma.cost[, 2] <- c(svm.gamma.range[2], svm.cost.range[2])

        # PSO
        set.seed(2000)

        result.internal.gamma.cost <- optim_pso(objective_function        = if (nlevels(dgeTraining$samples$group) == 2){
                                                                                  thrombo.svm.gamma.cost
                                                                            } else {
                                                                                  thrombo.svm.gamma.cost.multiclass
                                                                            },
                                                number_of_parameters      = nrow(parameter.bounds.gamma.cost),
                                                plot_progress             = FALSE,
                                                number_of_particles       = n.particles.gamma.cost.optimization,
                                                max_number_of_iterations  = n.iterations.gamma.cost.optimization,
                                                max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
                                                parameter_bounds          = parameter.bounds.gamma.cost,
                                                tryCall                   = TRUE,
                                                verbose                   = FALSE,
                                                lhc_init                  = FALSE,
                                                wait_complete_iteration   = TRUE,
                                                logfile                   = NULL,
                                                projectfile               = NULL,
                                                break_file                = "stopPPSO.txt"
        )

        # employ PSO proposed gamma and cost parameters
        if(class.weights == TRUE) {
          tuned.svm <- svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                           y           = real.groups.training,
                           gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                           cost        = as.numeric(result.internal.gamma.cost$par[2]),
                           tunecontrol = tune.control(cross = number.cross.splits),
                           probability = TRUE,
                           class.weights = 1/summary(dgeTraining$samples$group)
          )
        } else {
          tuned.svm <- svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                           y           = real.groups.training,
                           gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                           cost        = as.numeric(result.internal.gamma.cost$par[2]),
                           tunecontrol = tune.control(cross = number.cross.splits),
                           probability = TRUE
          ) }
        
        # extract best model
        tuned.svm.model <- tuned.svm
      }
      
      # prepare counts for sample prediction
      normalized.counts.prediction <- normalized.counts[selected.transcripts, evaluation.samples]
      
      # store data
      dgeTraining$biomarker.transcripts <- selected.transcripts
      dgeTraining$tuned.svm.model <- tuned.svm.model
      
      # this file is stored in the outputPSO-folder and can be employed for evaluation and validation 
      # if this particle is selected as the 'best option'
      save(dgeTraining, 
           file = paste("outputPSO/", paste(x, collapse = "-"), ".RData", sep = ""))
      
      # predict evaluation series
      prediction.class <- predict(tuned.svm.model,
                                  newdata = t(normalized.counts.prediction), 
                                  probability = TRUE)
      confusion.matrix <- table(prediction.class, real.groups.evaluation)
      confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
      
      # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
      if (nlevels(dgeTraining$samples$group) == 2){
        # create classification overview
        svm.summary <- data.frame(
          sampleName = attributes(prediction.class)$names,
          predicted = as.character((prediction.class)[1:length(prediction.class)]),
          real = real.groups.evaluation
        )
        svm.summary <- cbind(svm.summary, 
                             data.frame(attributes(prediction.class)$probabilities[,
                               which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                             data.frame(attributes(prediction.class)$probabilities[,
                               which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
        )
        colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
        # ROC
        rocra <- prediction(as.numeric(as.character(svm.summary[, 4])), ## adjusted 5 --> 4
                            svm.summary[, 3], 
                            label.ordering = rev(levels(dgeTraining$samples$group)) ## adjusted:added rev (it appeared cancer was on specificity side)
                            )
        perfa <- performance(rocra, "tpr", "fpr")
        if (verbose == TRUE){
          print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                      sep = ""))
        }
        roc.summary <- data.frame(
          cutOffs = unlist(attributes(rocra)$cutoffs),
          tp = unlist(attributes(rocra)$tp),
          tn = unlist(attributes(rocra)$tn),
          fp = unlist(attributes(rocra)$fp),
          fn = unlist(attributes(rocra)$fn),
          accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
            (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
          xValues = unlist(attributes(perfa)$x.values),
          yValues = unlist(attributes(perfa)$y.values)
        )
        roc.optimal.accuracy <- max(roc.summary$accuracy)
        
        if (rule.in.optimization == TRUE){
          # tune towards 99% specificity
          rule.in.out.setting = 99
          # get threshold at which the best specificity is most near by the user-provided specificity
          rule.threshold <- roc.summary[which(roc.summary$xValues == roc.summary$xValues[which.min(abs(roc.summary$xValues - (1 - (rule.in.out.setting / 100))))] & 
                                                roc.summary$yValues == max(roc.summary$yValues[which(
                                                  roc.summary$xValues == roc.summary$xValues[which.min(abs(roc.summary$xValues - (1 - (rule.in.out.setting / 100))))])])), ]
          # in case multiple results are returned, first take those with highest specificity,
          # it then multiple results remain, randomly select one option
          if (nrow(rule.threshold) > 1) {
            rule.threshold <- rule.threshold[which(rule.threshold$xValues == min(rule.threshold$xValues)), ]
          }
          if (nrow(rule.threshold) > 1) {
            set.seed(1000)
            rule.threshold <- rule.threshold[sample(1 : nrow(rule.threshold), size = 1), ]
          }
          
          svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
            if(svm.summary[i, 4] > rule.threshold$cutOffs - 1e-11) {
              svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
            } else {
              svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
            }
          }
          
          # prepare adjusted confusion matrix
          svm.summary$predicted.group <- factor(svm.summary$predicted.group, levels = levels(dgeTraining$samples$group)) ## updated 19th June 2022, MB (covers relevant groups now as levels)
          confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real)
          sensitivity <- confusion.matrix.updated[1] / (confusion.matrix.updated[1] + confusion.matrix.updated[2])
          # metrics for readout; report sensitivity
          if (verbose == TRUE){
            print(paste("Sensitivity: ", sensitivity, sep = ""))
          }
          AUC <- 1 - sensitivity
          return(AUC)
        } else {
          # metrics for readout
          AUC <- 1 - attributes(performance(rocra, 'auc'))$y.values[[1]]
          return(AUC)
        }
      } else {
        AUC <- 1 - confusion.matrix.evaluated$diag
        return(AUC)
      } 
    } else {
      AUC <- Inf
      return(AUC)
    }
  } else {
    AUC <- Inf
    return(AUC)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### svmrfeFeatureRanking ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Functions to identify the relative contribution of features (transcripts)
# to the SVM model for binary (svmrfeFeatureRanking-function) and multiclass
# (svmrfeFeatureRankingForMulticlass-function) algorithm development.
# the svm.weights function is employed in the other two individual functions
# adapted from http://www.bdmg.com.ar/wp-content/uploads/2011/11/SVM_RFE_R_implementation.pdf
svmrfeFeatureRanking <- function(x, 
                                 y,
                                 bestGamma,
                                 bestCost){
  n = ncol(x)
  survivingFeaturesIndexes = seq(1 : n)
  featureRankedList = vector(length = n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes) > 0){
    # train the support vector machine
    svmModel <- svm(x[, survivingFeaturesIndexes], y, cost = bestCost, gamma = bestGamma,
                    cachesize = 1000, scale = T, type = "C-classification")
    # compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    #compute ranking criteria
    rankingCriteria = w * w
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    # update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
  }
  return (featureRankedList)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### svmrfeFeatureRankingForMulticlass ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

svmrfeFeatureRankingForMulticlass <- function(x,
                                              y,
                                              bestCost,
                                              bestGamma){ 
  n = ncol(x) 
  survivingFeaturesIndexes = seq(1 : n) 
  featureRankedList = vector(length = n) 
  rankedFeatureIndex = n 
  while(length(survivingFeaturesIndexes) > 0){ 
    #train the support vector machine 
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = bestCost, gamma = bestGamma, cachesize = 1000,  
                   scale = T, type = "C-classification")
    #compute the weight vector 
    multiclassWeights = svm.weights(svmModel) 
    #compute ranking criteria 
    multiclassWeights = multiclassWeights * multiclassWeights 
    rankingCriteria = 0 
    for(i in 1:ncol(multiclassWeights)) 
      rankingCriteria[i] = mean(multiclassWeights[, i]) 
    #rank the features 
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix) 
    # update feature ranked list 
    (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]) 
    rankedFeatureIndex = rankedFeatureIndex - 1
    # eliminate the feature with smallest ranking criterion 
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]]) 
  } 
  return (featureRankedList)
} 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### svm.weights ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

svm.weights <- function(model){ 
  w = 0 
  if(model$nclasses == 2){ 
    w = t(model$coefs)%*%model$SV 
  } else {    
    #when we deal with OVO svm classification 
    ## compute start-index 
    start <- c(1, cumsum(model$nSV) + 1) 
    start <- start[-length(start)] 
    calcw <- function (i,j) { 
      ## ranges for class i and j: 
      ri <- start[i] : (start[i] + model$nSV[i] - 1) 
      rj <- start[j] : (start[j] + model$nSV[j] - 1) 
      ## coefs for (i,j): 
      coef1 <- model$coefs[ri, j-1] 
      coef2 <- model$coefs[rj, i] 
      ## return w values: 
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,] 
      return(w) 
    } 
    W=NULL 
    for (i in 1 : (model$nclasses - 1)){ 
      for (j in (i + 1) : model$nclasses){ 
        wi=calcw(i,j) 
        W=rbind(W,wi) 
      } 
    } 
    w=W 
  } 
  return(w) 
} 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### do ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Svm functions, binary and multiclass.
# Optimize towards AUC in binary classifier, and accuracy in multiclass classifier

# Function do, required by thrombo.svm.gamma.cost and thrombo.svm.gamma.cost.multiclass.
# determines whether the value outputted by do can be handled by R.
do <- function(x){
  2^x
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.svm.gamma.cost ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.svm.gamma.cost <- function(x) {
  # Function to improve gamma and cost settings in trained SVM model, binary comparison.
  # 
  # Args
  #   x: Vector with PSO-proposed gamma and cost values.
  #
  # Returns:
  #   Inverse AUC-value of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings.
  
  set.seed(1000) # lock randomness
  gammaRange <- x[1] # first value in x-vector is gamma value
  costRange <- x[2] # second value in x-vector is gamma value
  internalTraining <- colnames(normalized.counts) # training series
  internalTesting <- evaluation.samples # evaluation series
  realGroupsTesting <- droplevels(dgeTraining$samples[internalTesting, "group"])
  if(!do(costRange) == Inf){ # ensure the value for cost does not result in non-handible values
    if(!do(gammaRange) == Inf){ # ensure the value for gamma does not result in non-handible values
      # train SVM model with the proposed gamma and cost values
      m <- svm(t(normalized.counts)[internalTraining, selected.transcripts],
               droplevels(dgeTraining$samples[internalTraining, "group"]), 
               kernel = "radial", 
               cost = do(costRange), 
               gamma = do(gammaRange),
               type = "C-classification", 
               probability = TRUE, 
               cachesize = 500, 
               scale = T)
      # predict the evaluation series in the trained SVM model
      yscore <- attr(predict(m, t(normalized.counts)[internalTesting, selected.transcripts], probability= T), "probabilities")[, levels(realGroupsTesting)[1]]
      # construct AUC-value
      rocra <- prediction(as.numeric(as.character(yscore)), realGroupsTesting)
      perfa <- performance(rocra, "tpr", "fpr")
      AUC <- 1 - attributes(performance(rocra, 'auc'))$y.values[[1]]
      return(AUC)
    } else { 
      AUC <- Inf
      return(AUC)
    }
  } else { 
    AUC <- Inf
    return(AUC)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.svm.gamma.cost.multiclass ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.svm.gamma.cost.multiclass <- function(x) {
  # Function to improve gamma and cost settings in trained SVM model, multiclass comparison.
  # 
  # Args
  #   x: Vector with PSO-proposed gamma and cost values.
  #
  # Returns:
  #   Inverse accuracy of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings.
  
  set.seed(1000) # lock randomness
  gammaRange <- x[1] # first value in x-vector is gamma value
  costRange <- x[2] # second value in x-vector is gamma value
  internalTraining <- colnames(normalized.counts) # training series
  internalTesting <- evaluation.samples # evaluation series
  realGroupsTesting <- droplevels(dgeTraining$samples[internalTesting, "group"])
  if(!do(costRange) == Inf){# ensure the value for cost does not result in non-handible values
    if(!do(gammaRange) == Inf){ # ensure the value for gamma does not result in non-handible values
      # train SVM model with the proposed gamma and cost values
      m <- svm(t(normalizedCounts)[internalTraining, selected.transcripts],
               droplevels(dgeTraining$samples[internalTraining,"group"]), 
               kernel = "radial", 
               cost = do(costRange), 
               gamma = do(gammaRange),
               type = "C-classification", 
               probability = TRUE, 
               cachesize = 500, 
               scale = T)
      # predict the evaluation series in the trained SVM model and calculate inverse accuracy
      Acc <- 1 - classAgreement(table(predict(m, t(normalizedCounts)[
        internalTesting,selectedTranscripts_svmRanking], probability = T), realGroupsTesting), match.names = T)$diag
      return(Acc)
    } else { 
      Acc <- Inf
      return(Acc)
    }
  } else { 
    Acc <- Inf
    return(Acc)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thromboSeqPSO.readout ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thromboSeqPSO.readout <- function(dge = dgeIncludedSamples,
                                  pso.output.log.file = "ppso.log",
                                  pso.output.pro.file = "ppso.pro",
                                  readout.training = TRUE,
                                  readout.evaluation = TRUE,
                                  readout.validation = TRUE,
                                  replace.counts.validation = 0,
                                  filter.clinical.characteristics.validation = NA,
                                  filter.clinical.characteristics.group = NA,
                                  filter.clinical.characteristics.specified = NA,
                                  clinical.data.in.output = NA,
                                  apply.rule.in.readout = FALSE,
                                  rule.in.setting = NA,
                                  apply.rule.out.readout = FALSE,
                                  rule.out.setting = NA,
                                  additional.dge.to.predict = NA,
                                  R.snapshot = "Pre-PSO-snapshot.RData",
                                  figureDir = "figureOutputFolder",
                                  number.cores = 2,
                                  verbose = TRUE){
  # Perform readout of the PSO-enhanced thromboSeq classification algorithm. 
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   pso.output.log.file: String with log-file name of the PSO-optimization process.
  #   pso.output.pro.file: String with pro-file name of the PSO-optimization process.
  #   readout.training: Whether or not to perform LOOCV analysis of training series (TRUE/FALSE).
  #   readout.evaluation: Whether or not to classify samples of the evaluation series (TRUE/FALSE).
  #   readout.validation: Whether or not to classify samples of the validation series (TRUE/FALSE).
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value)
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   filter.clinical.characteristics.validation: String indicating a single clinical characteristics that may be
  #   limited for the validation series, e.g. only stage I-II cancer
  #   filter.clinical.characteristics.group: to which group in validation should filter be applied
  #   filter.clinical.characteristics.specified = String specifying the included clinical characteristic.
  #   clinical.data.in.output: String with additional clinical info to be presented in the svm.summary.
  #   apply.rule.in.readout: Whether or not to set algorithm readout towards rule-in setting (TRUE/FALSE)
  #   rule.in.setting: Numeric value indicating percentage of specificity in evaluation series at which 
  #   rule-in setting should be applied
  #   apply.rule.out.readout: Whether or not to set algorithm readout towards rule-out setting (TRUE/FALSE)
  #   rule.out.setting: Numeric value indicating percentage of sensitivity in evaluation series at which 
  #   rule-out setting should be applied
  #   additional.dge.to.predict: Provide a new DGE-object with samples to be classificated by trained algorithm
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Vector with best parameter settings. Also returns in working directory table with algorithm metrics. 
  #   In addition, outputs RData output files with specific data regarding the training, evaluation and 
  #   validation series, and a pdf-plot of the ROC-curve.
  
  if (all(!is.numeric(replace.counts.validation), !is.nan(replace.counts.validation))) {
    stop("provide numeric value or NaN for replace.counts.validation")
  }
  
  if (all(!clinical.data.in.output %in% c(colnames(dge$samples),NA))) {
    stop("Provide existing column names in clinical.data.in.output")
  }
  
  if (!rule.in.setting %in% c(NA, seq(1,100,by=1e-3))) {
    stop("rule.in.setting should be within 1 and 100")  
  }
  
  if (!rule.out.setting %in% c(NA, seq(1,100,by=1e-3))) {
    stop("rule.out.setting should be within 1 and 100")  
  }
  
  if (!filter.clinical.characteristics.validation %in% c(colnames(dge$samples), NA)) {
    stop("filter.clinical.characteristics.validation not available in dge$samples")
  }
  
  if (!filter.clinical.characteristics.group %in% c(levels(dge$samples$group), NA)) {
    stop("filter.clinical.characteristics.group not available in dge$samples$group")
  }
  
  if (!is.na(filter.clinical.characteristics.specified)) {
    if (!filter.clinical.characteristics.specified %in% c(levels(dge$samples[, filter.clinical.characteristics.validation]), NA)) {
      stop("filter.clinical.characteristics.specified not available in provided filter.clinical.characteristics.validation")
    }
  }
  
  if (!is.na(additional.dge.to.predict)) {
    if (verbose == TRUE) {
      print("Additional DGE for prediction provided.")
    }
  }
  
  if (missing(additional.dge.to.predict)) {
    stop("Provide DGElist object")
  }
  
  if (!is.na(additional.dge.to.predict)) {
    stopifnot(class(additional.dge.to.predict) == "DGEList")
  }
  
  # load required packages
  suppressMessages(library(foreach))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(e1071))
  suppressMessages(library(ppso))
  suppressMessages(library(ROCR))
  suppressMessages(library(pROC))
  suppressMessages(library(reshape))
  suppressMessages(library(ggplot2))
  
  # ensure all count tables in DGE are of same size and sample composition
  dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
  dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
  
  # Select the best parameter setting
  logged.PSO.distribution <- read.csv(pso.output.log.file, sep = "\t")
  logged.PSO.distribution <- logged.PSO.distribution[
    order(logged.PSO.distribution$objective_function), 
    ]
  logged.PSO.distribution.Index <- logged.PSO.distribution[
    which(logged.PSO.distribution$objective_function == min(logged.PSO.distribution$objective_function)), 
    ]
  
  set.seed(1000)
  # in case more particles have the same AUC output value, select randomly one as the particle for readout
  if (nrow(logged.PSO.distribution.Index) > 1){
    logged.PSO.distribution.Index <- logged.PSO.distribution.Index[
      sample(1 : nrow(logged.PSO.distribution.Index), size = 1),
      ]
  }
  best.selection <- paste(logged.PSO.distribution.Index[, c(seq(2,ncol(logged.PSO.distribution.Index)-2, by=1))], collapse = "-") # collapse data
  
  if (verbose == TRUE){
    print(paste("Best selection: ", best.selection, sep = ""))
  }
  
  if (readout.training == TRUE){
    results.classification.training <- thrombo.algo.classify.training.set(dge = dge,
                                                                          best.particle = best.selection,
                                                                          clinical.info.output = clinical.data.in.output,
                                                                          R.snapshot = "Pre-PSO-snapshot.RData",
                                                                          number.cores = 2)
    save(results.classification.training, file = "results.classification.training.RData")
  }
  
  if (readout.evaluation == TRUE){
    # If no rule-in or rule-out read-out is required, run default evaluation series read-out
    if (apply.rule.in.readout == FALSE & apply.rule.out.readout == FALSE) {
        results.classification.evaluation <- thrombo.algo.classify.evaluation.set(dge = dge,
                                                                                  best.particle = best.selection,
                                                                                  clinical.info.output = clinical.data.in.output,
                                                                                  rule.in = FALSE,
                                                                                  rule.out = FALSE,
                                                                                  R.snapshot = "Pre-PSO-snapshot.RData")
        save(results.classification.evaluation, file = "results.classification.evaluation.RData")
    } 
    if (apply.rule.in.readout == TRUE) {
      # apply rule-in setting
      results.classification.evaluation.rule.in <- thrombo.algo.classify.evaluation.set(dge = dge,
                                                                                best.particle = best.selection,
                                                                                clinical.info.output = clinical.data.in.output,
                                                                                rule.in = TRUE,
                                                                                rule.out = FALSE, 
                                                                                rule.in.out.setting = rule.in.setting,
                                                                                R.snapshot = "Pre-PSO-snapshot.RData")
      save(results.classification.evaluation.rule.in, file = "results.classification.evaluation.rule.in.RData")
      results.classification.evaluation <- results.classification.evaluation.rule.in
    } 
    if (apply.rule.out.readout == TRUE) {
      # apply rule-out setting
      results.classification.evaluation.rule.out <- thrombo.algo.classify.evaluation.set(dge = dge,
                                                                                best.particle = best.selection,
                                                                                clinical.info.output = clinical.data.in.output,
                                                                                rule.in = FALSE,
                                                                                rule.out = TRUE, 
                                                                                rule.in.out.setting = rule.in.setting,
                                                                                R.snapshot = "Pre-PSO-snapshot.RData")
      save(results.classification.evaluation.rule.out, file = "results.classification.evaluation.rule.out.RData")
      results.classification.evaluation <- results.classification.evaluation.rule.out
    }
  }  
  
  if (readout.validation == TRUE){
    # If no rule-in or rule-out read-out is required, run default evaluation series read-out
    if (apply.rule.in.readout == FALSE & apply.rule.out.readout == FALSE) {
      results.classification.validation <- thrombo.algo.classify.validation.set(dge.tool = dge,
                                                                                best.particle = best.selection,
                                                                                replace.counts.validation = replace.counts.validation,
                                                                                clin.characteristics.validation = filter.clinical.characteristics.validation,
                                                                                clin.characteristics.group = filter.clinical.characteristics.group,
                                                                                clin.characteristics.specified = filter.clinical.characteristics.specified,
                                                                                clinical.info.output = clinical.data.in.output,
                                                                                rule.in = FALSE,
                                                                                rule.out = FALSE,
                                                                                external.DGE = additional.dge.to.predict,
                                                                                rule.in.out.proposed.threshold = results.classification.evaluation$rule.threshold$cutOffs,
                                                                                R.snapshot = "Pre-PSO-snapshot.RData")
      if (!is.na(filter.clinical.characteristics.validation)){
        save(results.classification.validation, file = paste("results.classification.validation.",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".RData",sep=""))
      } else {
        save(results.classification.validation, file = "results.classification.validation.RData")  
      }
    } 
    if (apply.rule.in.readout == TRUE & readout.evaluation == TRUE) {
      # apply rule-in setting
      results.classification.validation.rule.in <- thrombo.algo.classify.validation.set(dge.tool = dge,
                                                                                        best.particle = best.selection,
                                                                                        clinical.info.output = clinical.data.in.output,
                                                                                        replace.counts.validation = replace.counts.validation,
                                                                                        clin.characteristics.validation = filter.clinical.characteristics.validation,
                                                                                        clin.characteristics.group = filter.clinical.characteristics.group,
                                                                                        clin.characteristics.specified = filter.clinical.characteristics.specified,
                                                                                        rule.in = TRUE,
                                                                                        rule.out = FALSE, 
                                                                                        external.DGE = additional.dge.to.predict,
                                                                                        rule.in.out.proposed.threshold = results.classification.evaluation.rule.in$rule.threshold$cutOffs,
                                                                                        R.snapshot = "Pre-PSO-snapshot.RData")
      if (!is.na(filter.clinical.characteristics.validation)){
        save(results.classification.validation.rule.in, file = paste("results.classification.validation.rule.in.",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".RData",sep=""))
      } else {
        save(results.classification.validation.rule.in, file = "results.classification.validation.rule.in.RData")  
      }
      results.classification.validation <- results.classification.validation.rule.in
    }
    if (apply.rule.out.readout == TRUE & readout.evaluation == TRUE) {
      # apply rule-out setting
      results.classification.validation.rule.out <- thrombo.algo.classify.validation.set(dge.tool = dge,
                                                                                         best.particle = best.selection,
                                                                                         clinical.info.output = clinical.data.in.output,
                                                                                         replace.counts.validation = replace.counts.validation,
                                                                                         clin.characteristics.validation = filter.clinical.characteristics.validation,
                                                                                         clin.characteristics.group = filter.clinical.characteristics.group,
                                                                                         clin.characteristics.specified = filter.clinical.characteristics.specified,
                                                                                         rule.in = FALSE,
                                                                                         rule.out = TRUE, 
                                                                                         external.DGE = additional.dge.to.predict,
                                                                                         rule.in.out.proposed.threshold = results.classification.evaluation.rule.out$rule.threshold$cutOffs,
                                                                                         R.snapshot = "Pre-PSO-snapshot.RData")
      if (!is.na(filter.clinical.characteristics.validation)){
        save(results.classification.validation.rule.out, file = paste("results.classification.validation.rule.out.",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".RData",sep=""))
      } else {
        save(results.classification.validation.rule.out, file = "results.classification.validation.rule.out.RData")  
      }
      results.classification.validation <- results.classification.validation.rule.out
    }
  }

  if (verbose == TRUE){
    print("summarize data")
  }
  
  # store biomarker panel
  # load particle-specific output file
  load(paste("outputPSO/", best.selection, ".RData", sep = ""))
  matrix <- as.matrix(dgeTraining$biomarker.transcripts)
  matrix <- cbind(matrix, as.character(dge$genes$hgnc_symbol[dgeTraining$biomarker.transcripts]),
                  as.character(dge$genes$description[dgeTraining$biomarker.transcripts]))
  write.csv(matrix, file = 'Swarm-enhanced-biomarker-panel.csv')
  
  # summarize data in a matrix
  matrix <- matrix(nrow = 3, ncol = 4)
  rownames(matrix) <- c("Training", "Evaluation", "Validation")
  colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy")
  if (readout.training == TRUE){
  matrix[1, ] <- c(results.classification.training$number.samples,
                   results.classification.training$AUCorDiagonal,
                  paste(round(results.classification.training$roc.95ci$ci[1], digits = 2), "-",
                        round(results.classification.training$roc.95ci$ci[3], digits = 2), sep = ""),
                  results.classification.training$roc.optimal.accuracy)
  }
  if (readout.evaluation == TRUE){
  matrix[2,] <- c(length(results.classification.evaluation$samples.for.evaluation),
                  results.classification.evaluation$AUCorDiagonal,
                  paste(round(results.classification.evaluation$roc.95ci$ci[1], digits = 2), "-",
                        round(results.classification.evaluation$roc.95ci$ci[3], digits = 2), sep = ""),
                  results.classification.evaluation$roc.optimal.accuracy)
  }
  if (readout.validation == TRUE){
  matrix[3,] <- c(length(results.classification.validation$samples.for.validation),
                  results.classification.validation$AUCorDiagonal,
                  paste(round(results.classification.validation$ci.roc$ci[1], digits = 2), "-",
                        round(results.classification.validation$ci.roc$ci[3], digits = 2), sep = ""),
                  results.classification.validation$roc.optimal.accuracy)
  }
  if (!is.na(filter.clinical.characteristics.validation)){
    write.csv(matrix, file = paste("ROCcurve_Metrics",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".csv",sep=""))
  } else {
    write.csv(matrix, file = "ROCcurve_Metrics.csv")  
  }
  
  # if binary classifier print ROC-curve. For multiclass comparisons print confusion matrices
  if (class(additional.dge.to.predict)=="DGEList") {
    n.levels.groups <- levels(additional.dge.to.predict$samples$group)
  } else {
    n.levels.groups <- NA
  }
  
  if (class(additional.dge.to.predict)=="DGEList" & n.levels.groups == 2 | is.na(additional.dge.to.predict)) {
    if (nlevels(dge$samples$group) == 2) {
      # Plot ROC-curve
      # check for correct output-directory for ROC Curve
      currentDir <- getwd()
      if (!is.na(filter.clinical.characteristics.validation)){
        if (file.exists(figureDir) == TRUE){
          pdf(paste(figureDir, "/ROCcurve-",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".pdf", sep = ""))  
        } else {
          setwd('..') # move to (previous) mother directory
          if (file.exists(figureDir) == TRUE){
            pdf(paste(figureDir, "/ROCcurve-",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".pdf", sep = ""))  
          } else {
            dir.create(figureDir)
            pdf(paste(figureDir, "/ROCcurve-",filter.clinical.characteristics.validation,"-",filter.clinical.characteristics.group,"-",paste(filter.clinical.characteristics.specified, collapse = "-"),".pdf", sep = ""))  
          }
        }
      } else {
        if (file.exists(figureDir) == TRUE){
          pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
        } else {
          setwd('..') # move to (previous) mother directory
          if (file.exists(figureDir) == TRUE){
            pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
          } else {
            dir.create(figureDir)
            pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
          }
        }
      }
      
      if (readout.training == TRUE){
        plot(results.classification.training$perfa, 
             lwd = 4, 
             # col = "grey")
             col = "#C0C0C0", lty=2)
        par(new=T)
      }
      if (readout.evaluation == TRUE){
        plot(results.classification.evaluation$perfa, 
             lwd = 6, 
             # col = "#B03B3D")
             col = "#969696")
        par(new=T)
      }
      if (readout.validation == TRUE){
        plot(results.classification.validation$perfa, 
             lwd = 6, 
             # col = "#3C66A6"
             col = "#B03B3D"
             # col = "#3C66A6"
               )
        par(new=T)
      }
      dev.off()
      setwd(currentDir)
    } else {
      # confusion matrices
      # check for correct output-directory for confusion matrices
      currentDir <- getwd()
      if (!file.exists(figureDir) == TRUE){
        setwd('..') # move to (previous) mother directory
        if (!file.exists(figureDir) == TRUE){
          dir.create(figureDir)
        }
      }
      if (readout.training == TRUE){
        confusion.matrix <- as.data.frame(
                                cast(results.classification.training$svm.summary, 
                                     predicted.group ~ real.group, 
                                     length, 
                                     value = "sample.ID"))
        rownames(confusion.matrix) <- confusion.matrix$predicted
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
        
        tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
          geom_tile(aes(fill = frequency)) +
          scale_fill_continuous(low = "white", high = "red") + 
          geom_text(aes(label = frequency)) +
          ggtitle("SVM classification results") + 
          labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
          theme_bw() + 
          theme(legend.position = "top")
        
        pdf(paste(figureDir, "/TrainingSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
        print(tiles)
        dev.off()
      }
      if (readout.evaluation == TRUE){
        confusion.matrix <- as.data.frame(
          cast(results.classification.evaluation$svm.summary, 
               predicted.group ~ real.group, 
               length, 
               value = "sampleName"))
        rownames(confusion.matrix) <- confusion.matrix$predicted
        confusion.matrix[order(colnames(confusion.matrix)[2:ncol(confusion.matrix)]),]
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
        tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
          geom_tile(aes(fill = frequency)) +
          scale_fill_continuous(low = "white", high = "red") + 
          geom_text(aes(label = frequency)) +
          ggtitle("SVM classification results") + 
          labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
          theme_bw() + 
          theme(legend.position = "top")
        pdf(paste(figureDir, "/EvaluationSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
        print(tiles)
        dev.off()
      }
      if (readout.validation == TRUE){
        confusion.matrix <- as.data.frame(
          cast(results.classification.validation$svm.summary, 
               predicted.group ~ real.group, 
               length, 
               value = "sampleName"))
        rownames(confusion.matrix) <- confusion.matrix$predicted
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
        
        tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
          geom_tile(aes(fill = frequency)) +
          scale_fill_continuous(low = "white", high = "red") + 
          geom_text(aes(label = frequency)) +
          ggtitle("SVM classification results") + 
          labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
          theme_bw() + 
          theme(legend.position = "top")
        pdf(paste(figureDir, "/ValidationSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
        print(tiles)
        dev.off()
      }
      
      setwd(currentDir)
    }
  } else {
    currentDir <- getwd()
    if (!file.exists(figureDir) == TRUE){
      setwd('..') # move to (previous) mother directory
      if (!file.exists(figureDir) == TRUE){
        dir.create(figureDir)
      }
    }
    write.csv(results.classification.validation$rule.in.out.confusion.matrix, 
              file = paste(figureDir, "/ValidationSeriesConfusionMatrix_additional.dge.to.predict.csv", sep = ""))
    setwd(currentDir)
  }
  return(best.selection)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thromboSeqPSO.controls ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thromboSeqPSO.controls <- function(dge = dgeIncludedSamples,
                                   filter.clinical.characteristics.validation = NA,
                                   filter.clinical.characteristics.group = NA,
                                   filter.clinical.characteristics.specified = NA,
                                   thromboSeqPSO.shuffled = TRUE,
                                   thromboSeqPSO.iterations = TRUE,
                                   n.shuffled = 1000,
                                   n.iterations = 1000,
                                   best.particle.input = thromboPSOreadout,
                                   replace.counts.validation = 0,
                                   n.particles.gamma.cost.optimization = 50,
                                   n.iterations.gamma.cost.optimization = 4,
                                   R.snapshot = "Pre-PSO-snapshot.RData",
                                   number.cores.controls = 2,
                                   verbose = TRUE) {
  # Perform shuffled labels and training series iterations for PSO-enhanced thromboSeq classifier.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   filter.clinical.characteristics.validation: String indicating a single clinical characteristics that may be
  #   limited for the validation series, e.g. only stage I-II cancer
  #   filter.clinical.characteristics.group: to which group in validation should filter be applied
  #   filter.clinical.characteristics.specified = String specifying the included clinical characteristic.
  #   thromboSeqPSO.shuffled: TRUE/FALSE whether or not to perform shuffled labels experiment in which the group 
  #                           labels of samples in the training series are randomly assigned. Indicates specificity 
  #                           of the developed classification algorithm.
  #   thromboSeqPSO.iterations: TRUE/FALSE whether or not to perform training series iterations experiment, in which 
  #                             the samplesassigned to the training series are randomly shuffled. Indicates the sensitivity 
  #                             of the developed classification algorithm
  #   n.shuffled: Vector with number of iterations for thromboSeqPSO.shuffled.
  #   n.iterations: Vector with number of iterations for thromboSeqPSO.iterations.
  #   best.particle.input: Vector with the merged best.selection values.
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value)
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores.controls: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Table with algorithm metrics.
  
  if (all(!is.numeric(replace.counts.validation), !is.nan(replace.counts.validation))) {
    stop("provide numeric value or NaN for replace.counts.validation")
  }
  
  if (!filter.clinical.characteristics.validation %in% c(colnames(dge$samples), NA)) {
    stop("filter.clinical.characteristics.validation not available in dge$samples")
  }
  
  if (!filter.clinical.characteristics.group %in% c(levels(dge$samples$group), NA)) {
    stop("filter.clinical.characteristics.group not available in dge$samples$group")
  }
  
  if (!is.na(filter.clinical.characteristics.specified)) {
    if (!filter.clinical.characteristics.specified %in% c(levels(dge$samples[, filter.clinical.characteristics.validation]), NA)) {
      stop("filter.clinical.characteristics.specified not available in provided filter.clinical.characteristics.validation")
    }
  }
  
  dge.tool <- dge
  
  if (thromboSeqPSO.shuffled == TRUE) {
    if (verbose == TRUE){
      print("start shuffled labels")
    }
    
    if (!file.exists("outputShuffled")) {
      dir.create("outputShuffled", recursive = T)
    }
    
    # loop the training process i times (total number of iterations) and perform in each iteration algorithm
    # training and validation. Results are stored in a container and summarized afterwards.
    suppressMessages(library(foreach))
    suppressMessages(library(doMC))
    registerDoMC(cores = number.cores.controls)
    shuffled.loop <- foreach(i = 1 : n.shuffled) %dopar% {
      
      results.classification.evaluation.shuffled <- thrombo.algo.classify.evaluation.shuffled(best.particle = best.particle.input,
                                                                                              shuffle.iteration = i)
      save(results.classification.evaluation.shuffled, 
           file = paste("outputShuffled/", i, "-evaluation.RData", sep = ""))
      results.classification.validation.shuffled <- thrombo.algo.classify.validation.shuffled(dge = dge.tool,
                                                                                              n.to.impute = replace.counts.validation,
                                                                                              best.particle = best.particle.input,
                                                                                              shuffle.iteration = i,
                                                                                              clin.characteristics.validation = filter.clinical.characteristics.validation,
                                                                                              clin.characteristics.group = filter.clinical.characteristics.group,
                                                                                              clin.characteristics.specified = filter.clinical.characteristics.specified)
      save(results.classification.validation.shuffled, 
           file = paste("outputShuffled/", i, "-validation.RData", sep = ""))
      
      results <- list()
      results[["i"]] <- i
      results[["evaluation.acc"]] <- results.classification.evaluation.shuffled$roc.optimal.accuracy
      results[["evaluation.AUC"]] <- results.classification.evaluation.shuffled$AUCorDiagonal
      results[["validation.acc"]] <- results.classification.validation.shuffled$roc.optimal.accuracy
      results[["validation.AUC"]] <- results.classification.validation.shuffled$AUCorDiagonal
      results
    }
    
    # summarize
    output.shuffled <- data.frame(
      i = unlist(lapply(shuffled.loop, function(x){x[["i"]]})),
      evaluation.acc = unlist(lapply(shuffled.loop, function(x){x[["evaluation.acc"]]})),
      evaluation.AUC = unlist(lapply(shuffled.loop, function(x){x[["evaluation.AUC"]]})),
      validation.acc = unlist(lapply(shuffled.loop, function(x){x[["validation.acc"]]})),
      validation.AUC = unlist(lapply(shuffled.loop, function(x){x[["validation.AUC"]]}))
    )
  }
  
  if (thromboSeqPSO.iterations == TRUE) {
    if (verbose == TRUE){
      print("start training set iterations")
    }
    
    if (!file.exists("outputIterations")) {
      dir.create("outputIterations", recursive = T)
    }
    
    n.iterations.datasets <- n.iterations
    
    # load R snapshot data
    load(R.snapshot)
   
     # loop the training process i times (total number of iterations) and perform in each iteration algorithm
    # training and validation. Results are stored in a container and summarized afterwards.
    suppressMessages(library(foreach))
    suppressMessages(library(doMC))
    registerDoMC(cores = number.cores.controls)
    iterations.loop <- foreach(iter = 1 : n.iterations.datasets) %dopar% {
      # reset seed for random selection of training-evaluation series
      rm(list=".Random.seed", envir = globalenv())
      
      dge.tool <- dge
      
      # load particle
      load(paste("outputPSO/", best.particle.input, ".RData", sep = ""))
      dgeParticle <- dgeTraining
      
      # randomly select samples for the training and evaluation series.
      # here it is assumed the group size and potential confounding factors 
      # (e.g. age of the individuals) are similar among both groups
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
      # store training series
      write.csv(training.samples, 
                paste("outputIterations/", iter, "-trainingSamples_subsampling.csv", sep = ""))
      
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
      # store evaluation series
      write.csv(evaluation.samples,
                paste("outputIterations/", iter, "-evaluationSamples_subsampling.csv", sep = ""))
      # save R environment
      save(list = ls(envir = environment(), all.names = TRUE), 
           file = paste("outputIterations/", iter, "-Pre-PSO-like.RData", sep = ""), 
           envir = environment()) 
      
      dge <- dge.tool
      
      # select series
      real.groups.training <- dge$samples[training.samples, "group"]
      real.groups.prediction <- dge$samples[evaluation.samples, "group"]
      
      dgeTraining <- perform.RUVg.correction(dge = dge[, c(training.samples, evaluation.samples)], 
                                             k.variables = k.variables, 
                                             variable.to.assess = variable.to.assess,
                                             variable.threshold = variable.threshold, 
                                             ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                             ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                             training.series.only = TRUE, 
                                             training.series.only.samples = training.samples,
                                             verbose = verbose
                                             )
      dgeTraining$counts <- dgeTraining$ruv.counts
      dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
      # TMM-normalisation
      dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.sample.readout = FALSE) # calculate normalization factors
      dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
      dgeTraining$samples <- droplevels(dgeTraining$samples)
      
      # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
      normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
      
      # first SVM-model, with grid search for gamma and cost
      tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                            y           = real.groups.training,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            probability = TRUE
      )
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      ## optimize gamma and cost by a second PSO-optimization algorithm
      # employ the training series for SVM algorithm training, and the evaluation series for optimization
      # make sure the values for cost and gamma are not infinite
      if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
        # if necessary - broaden the range of cost and gamma
        if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
          svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
        }
        if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
          svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
        }
        svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                             svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
        svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                            svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
        
        # input parameters:
        parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
        rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
        parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
        parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
        
        # PSO
        set.seed(2000)
        selected.transcripts <- dgeParticle$biomarker.transcripts
        result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
          thrombo.svm.gamma.cost
        } else {
          thrombo.svm.gamma.cost.multiclass
        },
        number_of_parameters      = nrow(parameter.bounds.gamma.cost),
        plot_progress             = FALSE,
        number_of_particles       = n.particles.gamma.cost.optimization,
        max_number_of_iterations  = n.iterations.gamma.cost.optimization,
        max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
        parameter_bounds          = parameter.bounds.gamma.cost,
        tryCall                   = TRUE, 
        verbose                   = FALSE,
        lhc_init                  = FALSE, 
        wait_complete_iteration   = TRUE,
        logfile                   = NULL,
        projectfile               = NULL,
        break_file                = "stopPPSO.txt"
        )
        
        # employ PSO proposed gamma and cost parameters
        tuned.svm <- svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                         y           = real.groups.training,
                         gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                         cost        = as.numeric(result.internal.gamma.cost$par[2]),
                         tunecontrol = tune.control(cross = number.cross.splits),
                         probability = TRUE
        )
        
        # extract best model
        tuned.svm.model <- tuned.svm
      }
      
      # prepare counts for sample prediction
      normalized.counts.prediction <- normalized.counts[dgeParticle$biomarker.transcripts, evaluation.samples]
      
      # store data
      dgeParticle$tuned.svm.model <- tuned.svm.model
      save(dgeParticle, 
           file = paste("outputIterations/", iter, "-dgeParticle.RData", sep = ""))
      
      # perform classification of evaluation and validation series and store output
      thrombo.algo.classify.evaluation.iterations <- thrombo.algo.classify.evaluation.set(dge = dgeTraining,
                                                                                          iterations = TRUE,
                                                                                          n_iter = iter,
                                                                                          clinical.info.output = NA,
                                                                                          rule.in = FALSE,
                                                                                          rule.out = FALSE,
                                                                                          R.snapshot = "Pre-PSO-snapshot.RData"
      )
      save(thrombo.algo.classify.evaluation.iterations, 
           file = paste("outputIterations/", iter, "-results.evaluation.iterations.RData", sep = ""))
      thrombo.algo.classify.validation.iterations <- thrombo.algo.classify.validation.set(dge.tool = dge,
                                                                                          iterations = TRUE,
                                                                                          n_iter = iter,
                                                                                          replace.counts.validation = replace.counts.validation,
                                                                                          clin.characteristics.validation = filter.clinical.characteristics.validation,
                                                                                          clin.characteristics.group = filter.clinical.characteristics.group,
                                                                                          clin.characteristics.specified = filter.clinical.characteristics.specified,
                                                                                          clinical.info.output = NA,
                                                                                          rule.in = FALSE,
                                                                                          rule.out = FALSE,
                                                                                          external.DGE = NA,
                                                                                          rule.in.out.proposed.threshold = 0.5,
                                                                                          R.snapshot = "Pre-PSO-snapshot.RData"
      )
      save(thrombo.algo.classify.validation.iterations, 
           file = paste("outputIterations/", iter, "-results.validation.iterations.RData", sep = ""))
      
      results <- list()
      results[["i"]] <- iter
      results[["evaluation.acc"]] <- thrombo.algo.classify.evaluation.iterations$roc.optimal.accuracy
      results[["evaluation.AUC"]] <- thrombo.algo.classify.evaluation.iterations$AUCorDiagonal
      results[["validation.acc"]] <- thrombo.algo.classify.validation.iterations$roc.optimal.accuracy
      results[["validation.AUC"]] <- thrombo.algo.classify.validation.iterations$AUCorDiagonal
      results
    }
    
    # summarize
    output.iterations <- data.frame(
      i = unlist(lapply(iterations.loop, function(x){x[["i"]]})),
      evaluation.acc = unlist(lapply(iterations.loop, function(x){x[["evaluation.acc"]]})),
      evaluation.AUC = unlist(lapply(iterations.loop, function(x){x[["evaluation.AUC"]]})),
      validation.acc = unlist(lapply(iterations.loop, function(x){x[["validation.acc"]]})),
      validation.AUC = unlist(lapply(iterations.loop, function(x){x[["validation.AUC"]]}))
    )
  }
  
  # summarize data
  # also check whether output files from the classifications are available
  matrix <- matrix(nrow = 6, ncol = 4)
  rownames(matrix) <- c("Training", "Evaluation", "Validation", "NA", "NA", "NA")
  colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy")
  if (file.exists("results.classification.training.RData") == TRUE){
    load("results.classification.training.RData")
    matrix[1, ] <- c(results.classification.training$number.samples,
                     results.classification.training$AUCorDiagonal,
                     paste(round(results.classification.training$roc.95ci$ci[1], digits = 2), "-",
                           round(results.classification.training$roc.95ci$ci[3], digits = 2), sep = ""),
                     results.classification.training$roc.optimal.accuracy)
  }
  
  if (file.exists("results.classification.evaluation.RData") == TRUE){
    load("results.classification.evaluation.RData")
    matrix[2,] <- c(length(results.classification.evaluation$samples.for.evaluation),
                    results.classification.evaluation$AUCorDiagonal,
                    paste(round(results.classification.evaluation$roc.95ci$ci[1], digits = 2), "-",
                          round(results.classification.evaluation$roc.95ci$ci[3], digits = 2), sep = ""),
                    results.classification.evaluation$roc.optimal.accuracy)
  }
  
  if (file.exists("results.classification.validation.RData") == TRUE){
    load("results.classification.validation.RData")
    matrix[3,] <- c(length(results.classification.validation$samples.for.evaluation),
                    results.classification.validation$AUCorDiagonal,
                    paste(round(results.classification.validation$ci.roc$ci[1], digits = 2), "-",
                          round(results.classification.validation$ci.roc$ci[3], digits = 2), sep = ""),
                    results.classification.validation$roc.optimal.accuracy)
  }
  if (thromboSeqPSO.shuffled == TRUE){
    matrix[5, 1] <- paste("Shuffled class labels (n=", nrow(output.shuffled), "): Median AUC Validation series: ",
                         round(median(output.shuffled$validation.AUC), digits = 2), ", and IQR: ",
                         round(IQR(output.shuffled$validation.AUC), digits = 2), sep = "")
    
  }
  if (thromboSeqPSO.iterations == TRUE){
    matrix[6, 1] <- paste("Iterations (n=", nrow(output.iterations), "): Median AUC Validation series: ",
                         round(median(output.iterations$validation.AUC), digits = 2), ", and IQR: ",
                         round(IQR(output.iterations$validation.AUC), digits = 2), sep = "")
    
  }
  write.csv(matrix, file = "ROCcurve_Metrics.csv")
  return(matrix)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### perform.RUVg.correction.validation ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

perform.RUVg.correction.validation <- function(dge = dge,
                                               readout.setting = NULL,
                                               training.samples.set = training.samples,
                                               evaluation.samples.set = evaluation.samples,
                                               validation.samples.set = validation.samples,
                                               output.particle = dgeParticle){
  # Performs the RUVSeq confounding variable correction specifically in a validation setting, i.e.
  # with known stable transcripts, confounding axes, and correction factors.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   readout.setting: Indicated whether training, evaluation, or validation series need to be RUV-corrected.
  #   thereby adjusts samples selected for RUVg correction.
  #   training.samples.set: provide vector with training samples.
  #   evaluation.samples.set: provide vector with evaluation samples.
  #   validation.samples.set: provide vector with validation samples.
  #   output.particle: DGEList compiled during the PSO process with the setting specific
  #                    stable transcripts, ruv-axes, and correction factors.
  #
  # Returns:
  #   DGEList including the corrected raw read counts.
  
  # remove the factors that were identified as potential confounding variables from the dataset
  # collect all output from the specific particle for correction of the counts of the to-be classified sample series
  axis.group <- output.particle$axis.group
  axis.na <- output.particle$axis.na
  axis.confounding <- output.particle$axis.confounding
  axis.all <- output.particle$axis.all
  axis.drop <- 0
  axis.removed <- FALSE
  k.variables <- length(output.particle$axis.all)
  
  # in case no RUV-correction is applied, skip this module
  if (all(!is.na(axis.confounding) & !is.na(axis.group))) {
    
    # steps:
    # for each validation sample loop the RUVg correction
    # previously all corrected at once has shown minor variation resulting in 
    # different predictive strength scores.
    # prepare output matrix
    
    dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
    dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
    
    # first correct all training samples (also include evaluation samples to reproduce 
    # RUVg correction in the PSO-optimised function)
    tmp.loop <- foreach(i = 1 : k.variables) %do% {
      if (i == 1) {
        if (i %in% axis.confounding) {
          RUVg.post.correction <- RUVg(as.matrix(dge[,c(training.samples.set, evaluation.samples.set)]$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = 0)
          axis.removed <- TRUE
        } else {
          axis.drop <- 1
        }
      } else {
        if (i %in% axis.confounding) {
          if(axis.removed == TRUE) {
            RUVg.post.correction <- RUVg(as.matrix(RUVg.post.correction$normalizedCounts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
          } else {
            RUVg.post.correction <- RUVg(as.matrix(dge[,c(training.samples.set, evaluation.samples.set)]$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
            axis.removed <- TRUE
          }
        } else {
          axis.drop <- axis.drop + 1
        }
      }
    }
    
    if (readout.setting %in% c('validation','validation.digitalSWARM')) { 
      if (!identical(axis.confounding, character(0))) {
        # loop all validation samples one-by-one (if only n=1 once)
        for (sample in validation.samples.set) {
          # show percentage processed
          # if (verbose == T) {
          if (which(sample == validation.samples.set) %in% seq(round(length(validation.samples.set) / 10), 
                                                               round(length(validation.samples.set) / 10) * 10, 
                                                               by = round(length(validation.samples.set) / 10))){
            print(paste(
              which(which(sample == validation.samples.set) == seq(round(length(validation.samples.set) / 10), 
                                                                   round(length(validation.samples.set) / 10) * 10, 
                                                                   by = round(length(validation.samples.set) / 10))) * 10, '% processed', sep = ""
            ))
            # }
          }
          
          axis.drop <- 0
          axis.removed <- FALSE
          tmp.loop <- foreach(i = 1 : k.variables) %do% {
            if (i == 1) {
              if (i %in% axis.confounding) {
                RUVg.post.correction.sample <- RUVg(as.matrix(dge[,c(training.samples.set, evaluation.samples.set, sample)]$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = 0)
                axis.removed <- TRUE
              } else {
                axis.drop <- 1
              }
            } else {
              if (i %in% axis.confounding) {
                if(axis.removed == TRUE) {
                  RUVg.post.correction.sample <- RUVg(as.matrix(RUVg.post.correction.sample$normalizedCounts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
                } else {
                  RUVg.post.correction.sample <- RUVg(as.matrix(dge[,c(training.samples.set, evaluation.samples.set, sample)]$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
                  axis.removed <- TRUE
                }
              } else {
                axis.drop <- axis.drop + 1
              }
            }
          }
          RUVg.post.correction$normalizedCounts <- cbind(RUVg.post.correction$normalizedCounts, RUVg.post.correction.sample$normalizedCounts[,sample])
        }
        colnames(RUVg.post.correction$normalizedCounts) <- c(training.samples.set, evaluation.samples.set, validation.samples.set)
      } else {
        RUVg.post.correction <- NULL
        RUVg.post.correction$normalizedCounts <- (dge$counts)
      }
    }
    
    # prepare a new corrected countmatrix, update the total library size
    # reduce dge to only training and validation series
    if (readout.setting %in% c('LOOCV','evaluation')) {
      dge$raw.counts <- dge$counts
      if (length(axis.confounding) > 0) {
        # if RUVg correction has been performed, add the updated countmatrix to the dge-object
        dge$ruv.counts <- as.data.frame(RUVg.post.correction$normalizedCounts)
      } else {
        dge$ruv.counts <- dge$raw.counts
      }
      dge$samples$ruv.lib.size <- as.numeric(colSums(dge$ruv.counts)) 
    } else if (readout.setting == 'validation') {
      dge <- dge[,c(training.samples.set, validation.samples.set)]
      dge$samples <- droplevels(dge$samples)
      dge$raw.counts <- dge[,c(training.samples.set, validation.samples.set)]$counts
      if (length(axis.confounding) > 0) {
        # if RUVg correction has been performed, add the updated countmatrix to the dge-object
        dge$ruv.counts <- as.data.frame(RUVg.post.correction$normalizedCounts[,c(training.samples.set, validation.samples.set)])
      } else {
        dge$ruv.counts <- dge$ruv.counts[,c(training.samples.set, validation.samples.set)]
      }
      dge$samples$ruv.lib.size <- as.numeric(colSums(dge$ruv.counts))
    } else if (readout.setting == 'validation.digitalSWARM') {
      dge <- dge[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]
      dge$samples <- droplevels(dge$samples)
      dge$raw.counts <- dge[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]$counts
      if (length(axis.confounding) > 0) {
        # if RUVg correction has been performed, add the updated countmatrix to the dge-object
        dge$ruv.counts <- as.data.frame(RUVg.post.correction$normalizedCounts[,c(training.samples.set, evaluation.samples.set, validation.samples.set)])
      } else {
        dge$ruv.counts <- dge$ruv.counts[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]
      }
      dge$samples$ruv.lib.size <- as.numeric(colSums(dge$ruv.counts))
    }
  } else {
    if (readout.setting %in% c('LOOCV','evaluation')){
      dge$counts <- dge$counts[,c(training.samples.set, evaluation.samples.set)]
      dge$samples <- dge$samples[c(training.samples.set, evaluation.samples.set),]
      dge$raw.counts <- dge$counts[,c(training.samples.set, evaluation.samples.set)]
      dge$ruv.counts <- dge$counts[,c(training.samples.set, evaluation.samples.set)]
      dge$samples$ruv.lib.size <- dge$samples$lib.size
      dge$samples <- droplevels(dge$samples)  
    } else if (readout.setting == 'validation') {
      dge$counts <- dge$counts[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]
      dge$samples <- dge$samples[c(training.samples.set, evaluation.samples.set, validation.samples.set),]
      dge$raw.counts <- dge$raw.counts[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]
      dge$ruv.counts <- dge$ruv.counts[,c(training.samples.set, evaluation.samples.set, validation.samples.set)]
      dge$samples <- droplevels(dge$samples)
    }
  }
  # return corrected DGEList
  return(dge)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.classify.training.set ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.classify.training.set <- function(dge = dge,
                                               n.particles.gamma.cost.optimization = 50,
                                               n.iterations.gamma.cost.optimization = 4,
                                               best.particle = best.particle.input,
                                               clinical.info.output = clinical.data.in.output,
                                               R.snapshot = "Pre-PSO-snapshot.RData",
                                               number.cores = number.cores,
                                               verbose = TRUE){
  # Performs leave-one-out cross-validation (LOOCV) analysis of the training samples series.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   best.particle: Vector with the merged best.selection values.
  #   clinical.info.output: Additional clinical info present in dge$samples to be included in svm.summary
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #   
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  
  suppressMessages(library(foreach))
  suppressMessages(library(doMC))
  registerDoMC(cores = number.cores)
  loocv.loop <- foreach(leave.out.sample = training.samples) %dopar% {
    # assign samples
    training.samples.loocv <- training.samples[leave.out.sample != training.samples]
    real.groups.training <- dge$samples[training.samples.loocv, "group"]
    real.groups.prediction <- dge$samples[leave.out.sample, "group"]
    # load particle data
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
   
    dgeParticle <- dgeTraining
    # perform RUVg correction
    dgeTraining <- perform.RUVg.correction.validation(dge = dge[,c(training.samples, evaluation.samples)], 
                                                      readout.setting = 'LOOCV',
                                                      training.samples.set = training.samples,
                                                      evaluation.samples.set = evaluation.samples,
                                                      output.particle = dgeParticle)
    dgeTraining$counts <- dgeTraining$ruv.counts
    dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
    
    # TMM-normalization
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = FALSE) # calculate normalization factors
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
    dgeTraining$samples <- droplevels(dgeTraining$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
    
    # train SVM-model, with grid search for gamma and cost.
    tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples.loocv]),
                          y           = real.groups.training,
                          gamma       = svm.gamma.range,
                          cost        = svm.cost.range,
                          tunecontrol = tune.control(cross = number.cross.splits),
                          probability = TRUE
    )
    # extract best model
    tuned.svm.model <- tuned.svm[["best.model"]]
    
    if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
      # if necessary - broaden the range of cost and gamma
      if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
        svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
      }
      if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
        svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
      }
      svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                           svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
      svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                          svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
      
      # input parameters:
      parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
      rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
      parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
      parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
      
      # PSO
      set.seed(2000)
      selected.transcripts <- dgeParticle$biomarker.transcripts
      result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
        thrombo.svm.gamma.cost
      } else {
        thrombo.svm.gamma.cost.multiclass
      },
      number_of_parameters      = nrow(parameter.bounds.gamma.cost),
      plot_progress             = FALSE,
      number_of_particles       = n.particles.gamma.cost.optimization,
      max_number_of_iterations  = n.iterations.gamma.cost.optimization,
      max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
      parameter_bounds          = parameter.bounds.gamma.cost,
      tryCall                   = TRUE,
      verbose                   = FALSE,
      lhc_init                  = FALSE,
      wait_complete_iteration   = TRUE,
      logfile                   = NULL,
      projectfile               = NULL,
      break_file                = "stopPPSO.txt"
      )
      
      # employ PSO proposed gamma and cost parameters
      tuned.svm <- svm(x           = t(normalized.counts[selected.transcripts, training.samples.loocv]),
                       y           = real.groups.training,
                       gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                       cost        = as.numeric(result.internal.gamma.cost$par[2]),
                       tunecontrol = tune.control(cross = number.cross.splits),
                       probability = TRUE
      )
      
      # extract best model
      tuned.svm.model <- tuned.svm
    }
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, leave.out.sample]
    # perform prediction
    prediction.class <- predict(tuned.svm.model,
                                newdata = t(normalized.counts.prediction), 
                                probability = TRUE)
    
    # summarize results
    result <- list()
    tryName <- paste("Try", leave.out.sample, sep = "")
    result[["training.samples"]] <- training.samples.loocv
    result[["prediction.sample"]] <- leave.out.sample
    result[["predicted.group"]] <- as.character(prediction.class)
    result[["real.group"]] <- as.character(real.groups.prediction)
    result[["biomarker.panel"]] <- length(dgeParticle$biomarker.transcripts)
    result[["predictive.strength"]] <- attributes(prediction.class)$probabilities
    result
  }
  save(loocv.loop, file = 'loocv-loop.RData')
  # load additional packages
  suppressMessages(library(reshape, warn.conflicts = F, quietly = T))
  suppressMessages(library(ROCR, warn.conflicts = F, quietly = T))
  suppressMessages(library(pROC, warn.conflicts = F, quietly = T))
  
  # summarize data into a data frame
  svm.summary <- data.frame(
    sample.ID = unlist(lapply(loocv.loop, function(x){x[["prediction.sample"]]})),
    predicted.group = unlist(lapply(loocv.loop, function(x){x[["predicted.group"]]})),
    real.group = unlist(lapply(loocv.loop, function(x){x[["real.group"]]})),
    biomarker.panel = unlist(lapply(loocv.loop, function(x){x[["biomarker.panel"]]}))
  )
  # add predictive values to the data frame
  svm.summary <- suppressWarnings(cbind(svm.summary, do.call("rbind", lapply(loocv.loop, function(x){x[["predictive.strength"]]}))))
  svm.summary <- cbind(svm.summary[,1:4], 
                       svm.summary[which(colnames(svm.summary)==levels(svm.summary$real.group)[1])],
                       svm.summary[which(colnames(svm.summary)==levels(svm.summary$real.group)[2])]
  )
  for.roc <- ncol(svm.summary)
  
  # add clinical info to svm.summary
  if (!is.na(clinical.info.output)) {
    col.names <- colnames(svm.summary)
    svm.summary <- cbind(svm.summary, dge$samples[training.samples, clinical.info.output])
    colnames(svm.summary) <- c(col.names, clinical.info.output)
  }
  
  # summarize in confusion matrix
  confusion.matrix <- cast(svm.summary, 
                           predicted.group ~ real.group, 
                           length, 
                           value = "sample.ID")
  confusion.matrix.evaluated <- classAgreement(as.matrix(confusion.matrix), match.names = T)
  
  if (nlevels(dgeTraining$samples$group) == 2){
    # in case two classes are included, generate ROC curve
    rocra <- prediction(as.numeric(as.character(svm.summary[, 6])), 
                        svm.summary[, 3], 
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Training Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, for.roc], 
                    ci = TRUE
    )
  } else {
    perfa <- NA
    AUC <- confusion.matrix.evaluated$diag
    roc.optimal.accuracy <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.summary <- NA
  }
  # summarize data
  result <- list()
  result[["number.samples"]] <- length(training.samples)
  result[["svm.summary"]] <- svm.summary
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["roc.95ci"]] <- roc.95ci
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  return(result)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.classify.evaluation.set ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.classify.evaluation.set <- function(dge = dge, 
                                                 best.particle = best.particle.input, 
                                                 iterations = FALSE,
                                                 n_iter = NULL,
                                                 clinical.info.output = clinical.data.in.output,
                                                 rule.in = NA,
                                                 rule.out = NA,
                                                 rule.in.out.setting = NA,
                                                 R.snapshot = "Pre-PSO-snapshot.RData",
                                                 verbose = TRUE){
  # Performs classification of samples included in the evaluation series in a PSO-optimized
  # classification algorithm.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   iterations: Whether or not (TRUE/FALSE) this function is called in the iterations control experiments.
  #   n_iter: Numeric value with the number of the iteration (only when iterations == TRUE).
  #   clinical.info.output: Additional clinical info present in dge$samples to be included in svm.summary
  #   rule.in: Whether or not rule-in setting should be applied (TRUE/FALSE)
  #   rule.out: Whether or not rule-out setting should be applied (TRUE/FALSE)
  #   rule.in.out.setting: Numeric value at which percentage rule-in or -out should be assessed
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  
  # ensure all count tables in DGE are of same size and sample composition
  dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
  dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
  
  # load data according to whether the iterations are enabled or not
  if (iterations == TRUE){
    load(paste("outputIterations/", n_iter, "-Pre-PSO-like.RData", sep = ""))
    load(paste("outputIterations/", n_iter, "-dgeParticle.RData", sep=""))
    
    training.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-trainingSamples_subsampling.csv", sep = ""))[, 2])
    evaluation.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-evaluationSamples_subsampling.csv", sep = ""))[, 2])
  } else {
    # load particle-specific output file
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
    dgeParticle <- dgeTraining
  }
  
  # assign samples to training and evaluation
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  # perform RUV correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge[, c(training.samples, evaluation.samples)], 
                                                    training.samples.set = training.samples,
                                                    evaluation.samples.set = evaluation.samples,
                                                    validation.samples.set = validation.samples,
                                                    readout.setting = 'evaluation',
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  # normalize using TMM normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE, 
                                           refColumn = dgeParticle$refSample) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, evaluation.samples]
  # perform classification
  prediction.class <- predict(dgeParticle$tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.evaluation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.evaluation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  for.roc <- ncol(svm.summary)
  
  # add clinical info to svm.summary
  if (!is.na(clinical.info.output)) {
  col.names <- colnames(svm.summary)
  svm.summary <- cbind(svm.summary, dgeTraining$samples[evaluation.samples, clinical.info.output])
  colnames(svm.summary) <- c(col.names, clinical.info.output)
  }
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC
    rocra <- prediction(predictions = as.numeric(as.character(svm.summary[, 4])), 
                        labels = svm.summary[, 3], 
                        label.ordering = rev(levels(dgeTraining$samples$group)))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, for.roc], 
                    ci = TRUE
    )
    
    # if rule-in or rule-out read-out is required, select predictive strength threshold at which
    # the provided rule-in/out setting should be applied
    if (is.numeric(rule.in.out.setting)) {
      if (rule.in == TRUE) {
        # get threshold at which the best specificity is most near by the user-provided sensitivity
        rule.threshold <- roc.summary[which(roc.summary$xValues == roc.summary$xValues[which.min(abs(roc.summary$xValues - (1 - (rule.in.out.setting / 100))))] & 
                                              roc.summary$yValues == max(roc.summary$yValues[which(
                                                roc.summary$xValues == roc.summary$xValues[which.min(abs(roc.summary$xValues - (1 - (rule.in.out.setting / 100))))])])), ]
        # in case multiple results are returned, first take those with highest specificity,
        # it then multiple results remain, randomly select one option
        if (nrow(rule.threshold) > 1) {
          rule.threshold <- rule.threshold[which(rule.threshold$xValues == min(rule.threshold$xValues)), ]
        }
        if (nrow(rule.threshold) > 1) {
          set.seed(1000)
          rule.threshold <- rule.threshold[sample(1 : nrow(rule.threshold), size = 1), ]
        }
      } else if (rule.out == TRUE) {
        # get threshold at which the best sensitivity is most near by the user-provided sensitivity 
        rule.threshold <- roc.summary[which(roc.summary$yValues == roc.summary$yValues[which.min(abs(roc.summary$yValues - (rule.in.out.setting / 100)))] & 
                                      roc.summary$xValues == min(roc.summary$xValues[which(
                                        roc.summary$yValues == roc.summary$yValues[which.min(abs(roc.summary$yValues - (rule.in.out.setting / 100)))])])), ]
        # in case multiple results are returned, first take those with highest sensitivity,
        # it then multiple results remain, randomly select one option
        if (nrow(rule.threshold) > 1) {
          rule.threshold <- rule.threshold[which(rule.threshold$yValues == max(rule.threshold$yValues)), ]
        } 
        if (nrow(rule.threshold) > 1) {
          set.seed(1000)
          rule.threshold <- rule.threshold[sample(1 : nrow(rule.threshold), size = 1), ]
        }
      }
      
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fifth column of svm.summary is over the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.threshold$cutOffs - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
    } else {
      # get threshold at which the best accuracy is achieved
      rule.threshold <- roc.summary[which(roc.summary$accuracy == max(roc.summary$accuracy)),]
      # in case multiple results are returned, randomly select one option
      if (nrow(rule.threshold) > 1) {
        set.seed(1000)
        rule.threshold <- rule.threshold[sample(1 : nrow(rule.threshold), size = 1), ]
      }
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fifth column of svm.summary is over the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.threshold$cutOffs - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
    }
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
    rule.threshold <- NA
    confusion.matrix.updated <- NA
  }
  
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.evaluation"]] <- evaluation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.95ci"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result[["rule.threshold"]] <- rule.threshold
  result[["rule.in.out.confusion.matrix"]] <- confusion.matrix.updated
  result
  return(result)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.classify.validation.set ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.classify.validation.set <- function(dge.tool = dge, 
                                                 best.particle = best.particle.input, 
                                                 iterations = FALSE,
                                                 n_iter = NULL,
                                                 replace.counts.validation = 0,
                                                 clin.characteristics.validation = NA,
                                                 clin.characteristics.group = NA,
                                                 clin.characteristics.specified = NA,
                                                 clinical.info.output = clinical.data.in.output,
                                                 rule.in = NA,
                                                 rule.out = NA,
                                                 rule.in.out.proposed.threshold = 0.5,
                                                 external.DGE = NA,
                                                 R.snapshot = "Pre-PSO-snapshot.RData",
                                                 verbose = TRUE){
  # Performs classification of samples included in the validation series in a PSO-optimized
  # classification algorithm.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   iterations: Whether or not (TRUE/FALSE) this function is called in the iterations control experiments.
  #   iter: Numeric value with the number of the iteration (only when iterations == TRUE).
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value) 
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   clinical.info.output: Additional clinical info present in dge$samples to be included in svm.summary
  #   rule.in: Whether or not rule-in setting read-out should be applied (TRUE/FALSE)
  #   rule.out: Whether or not rule-out setting read-out should be applied (TRUE/FALSE)
  #   rule.in.out.proposed.threshold: Numeric value derived from evaluation series at which
  #   predictive strength threshold the rule-in or -out should be assessed
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  
  # replace loaded dge by dge provided in function input arguments, prevents wrong or larger DGE previously compiled
  # to be employed for coming analyses
  dge <- dge.tool
  # ensure all count tables in DGE are of same size and sample composition
  dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
  dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
  
  # load data according to whether the iterations are enabled or not
  if (iterations == TRUE){
    load(paste("outputIterations/", n_iter, "-Pre-PSO-like.RData", sep = ""))
    load(paste("outputIterations/", n_iter, "-dgeParticle.RData", sep = ""))
    
    training.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-trainingSamples_subsampling.csv", sep = ""))[, 2])
    evaluation.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-evaluationSamples_subsampling.csv", sep = ""))[, 2])
  } else {
    # load particle-specific output file
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
    dgeParticle <- dgeTraining
  }
  
  # assign samples to training and evaluation
  if (!is.na(clin.characteristics.validation) & is.na(external.DGE)) {
    validation.samples <- c(
      colnames(dge)[dge$samples$group %in% clin.characteristics.group &
                      !colnames(dge) %in% c(training.samples, evaluation.samples) &
                      dge$samples[,clin.characteristics.validation] %in% clin.characteristics.specified],
      colnames(dge)[!dge$samples$group %in% clin.characteristics.group &
                      !colnames(dge) %in% c(training.samples, evaluation.samples)]) 
  } else if (!is.na(external.DGE)) {
    # prepare new DGE with both dataset on which has been trained and
    # via external.DGE the newest dataset to be validated
    dgeNew <- DGEList(counts = cbind(dge$counts, external.DGE$counts[rownames(dge$counts), ]),
                      group = c(as.character(dge$samples$group), as.character(external.DGE$samples$group)),
                      genes = dge$genes
    )
    validation.samples <- colnames(dgeNew)[which(!colnames(dgeNew) %in% colnames(dge))]
    dge <- dgeNew
  } else {
    validation.samples <- colnames(dge)[!colnames(dge) %in% c(training.samples, evaluation.samples)]    
  }
    
  # narrow the dge to those samples relevant
  # dge <- dge[, c(training.samples, evalvalidation.samples)]
  # ensure all count tables in DGE are of same size and sample composition
  dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
  dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
  dge$samples <- droplevels(dge$samples)
  
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.validation <- dge$samples[validation.samples, "group"]
  
  # perform RUV correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge, 
                                                    readout.setting = 'validation',
                                                    training.samples.set = training.samples,
                                                    evaluation.samples.set = evaluation.samples,
                                                    validation.samples.set = validation.samples,
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  
  # enable to replace counts with 0 to provided counts in the validation series
  # by the median of those in the training series.
  if (!replace.counts.validation %in% c("NaN", NaN)){ # omit this function when NaN is inputted
    for (assess.sample in validation.samples) { # for each sample in the validation series
      tmpA <- matrix(dgeTraining$counts[, assess.sample]) # select the counts
      sel <- which(tmpA %in% seq(0, replace.counts.validation)) # identify which counts have too little raw reads detected
      if (length(sel) > 1){ # if more than one selected, calculate median read counts of these genes in the training series and replace
        tmpB <- round(apply(dgeTraining$counts[sel, training.samples], 1, median))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, replace.counts.validation)), assess.sample] <- tmpB
      } else if (length(sel) == 1) { # if one selected, calculate median read counts of this gene in the training series and replace
        tmpB <- round(median(as.numeric(dgeTraining$counts[sel, ])))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, replace.counts.validation)), assess.sample] <- tmpB
      }
    }
  }
  
  # calculate newest total read counts (lib.size)
  dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = F, 
                                           refColumn = dgeParticle$refSample) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, validation.samples]
  # perform classification
  prediction.class <- predict(dgeParticle$tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.validation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.validation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  for.roc <- ncol(svm.summary)
  
  # add clinical info to svm.summary
  if (!is.na(clinical.info.output) & is.na(external.DGE)) {
    col.names <- colnames(svm.summary)
    svm.summary <- cbind(svm.summary, dgeTraining$samples[validation.samples, clinical.info.output])
    colnames(svm.summary) <- c(col.names, clinical.info.output)
  }
  if (!is.na(external.DGE) & !is.na(clinical.info.output)) {
    col.names <- colnames(svm.summary)
    svm.summary <- cbind(svm.summary, external.DGE$samples[validation.samples, clinical.info.output])
    colnames(svm.summary) <- c(col.names, clinical.info.output)
  }
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC
    rocra <- prediction(as.numeric(as.character(svm.summary[, 4])), 
                        svm.summary[, 3], 
                        label.ordering = rev(levels(dgeTraining$samples$group)))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Validation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, for.roc], 
                    ci = TRUE
    )
    
    # if rule-in or rule-out read-out is required, take predictive strength threshold as determined in 
    # the evaluation series, and apply to validation series
    if (any(c(rule.in, rule.out)) == TRUE){
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fourth column of svm.summary is more than the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.in.out.proposed.threshold - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
    } else {
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fourth column of svm.summary is more than the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.in.out.proposed.threshold - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
    }
  } else if (nlevels(dgeTraining$samples$group) > 2 & !is.na(external.DGE)) {
    # if rule-in or rule-out read-out is required, take predictive strength threshold as determined in 
    # the evaluation series, and apply to validation series
    if (any(c(rule.in, rule.out)) == TRUE) {
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fourth column of svm.summary is more than the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.in.out.proposed.threshold - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
      AUC <- NA
      roc.95ci <- list(NA)
      roc.95ci$ci <- NA
      roc.optimal.accuracy <- confusion.matrix.evaluated$diag
      perfa <- NA
      roc.summary <- NA
    } else {
      # re-assigned predicted class labels in svm.summary according to newly provided threshold
      # in case the observed predictive strength in the fifth column of svm.summary is over the proposed threshold,
      # assign this group name (encoded in the colnames). Otherwise, assign the other group name.
      svm.summary$predicted.group <- foreach(i = 1 : nrow(svm.summary)) %do% {
        if(svm.summary[i, 4] > rule.in.out.proposed.threshold - 1e-11) {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[4]
        } else {
          svm.summary$predicted.group[i] <- colnames(svm.summary)[5]
        }
      }
      
      # prepare adjusted confusion matrix
      confusion.matrix.updated <- table(unlist(svm.summary$predicted.group), svm.summary$real.group)
      AUC <- NA
      roc.95ci <- list(NA)
      roc.95ci$ci <- NA
      roc.optimal.accuracy <- confusion.matrix.evaluated$diag
      perfa <- NA
      roc.summary <- NA
    }
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
    confusion.matrix.updated <- NA
  }
  
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.validation"]] <- validation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["ci.roc"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result[["rule.in.out.confusion.matrix"]] <- confusion.matrix.updated
  result
  return(result)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.classify.evaluation.shuffled ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.classify.evaluation.shuffled <- function(dge = dge, 
                                                      best.particle = best.particle.input, 
                                                      svm.gamma.range = 2^(-20:0),
                                                      svm.cost.range = 2^(0:20),
                                                      number.cross.splits = 2,
                                                      n.particles.gamma.cost.optimization = 50,
                                                      n.iterations.gamma.cost.optimization = 4,
                                                      shuffle.iteration = i,
                                                      R.snapshot = "Pre-PSO-snapshot.RData",
                                                      number.cores = number.cores,
                                                      verbose = TRUE){
  # Performs development of a classification algorithm based on the biomarker panel in the best PSO 
  # particle but with shuffled group labels of the training series.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   input: Vector with the merged best.selection values.
  #   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
  #   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
  #   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   shuffle.iteration: Numeric value to be provided by thromboSeqPSO.controls-function.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Result-container with classification details and metrics.
  
  # load R environment
  load(R.snapshot)
  # reset seed for random selection of training-evaluation series
  rm(list=".Random.seed", envir = globalenv())
  
  # load particle data
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  dgeParticle <- dgeTraining
  
  # collect group conditions
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  # perform RUVg correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge[, c(training.samples, evaluation.samples)], 
                                                    readout.setting = 'evaluation',
                                                    training.samples.set = training.samples,
                                                    evaluation.samples.set = evaluation.samples,
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  
  # TMM-normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
  
  ### 
  ### shuffle here the group labels of the samples in the training series
  ### of note, the set.seed-function is not employed here to ensure each classification is random
  real.groups.training <- sample(dge$samples[training.samples, "group"])
  
  # first SVM-model, with grid search for gamma and cost
  tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                        y           = real.groups.training,
                        gamma       = svm.gamma.range,
                        cost        = svm.cost.range,
                        tunecontrol = tune.control(cross = number.cross.splits),
                        probability = TRUE
  )
  # extract best model
  tuned.svm.model <- tuned.svm[["best.model"]]
  
  ## optimize gamma and cost by a second PSO-optimization algorithm
  # employ the training series for SVM algorithm training, and the evaluation series for optimization
  # make sure the values for cost and gamma are not infinite
  if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
    # if necessary - broaden the range of cost and gamma
    if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
      svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
    }
    if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
      svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
    }
    svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                         svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
    svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                        svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
    
    # input parameters:
    parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
    rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
    parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
    parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
    
    # PSO
    set.seed(2000)
    selected.transcripts <- dgeParticle$biomarker.transcripts
    result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
                                                                    thrombo.svm.gamma.cost
                                                                  } else {
                                                                    thrombo.svm.gamma.cost.multiclass
                                                                  },
    number_of_parameters      = nrow(parameter.bounds.gamma.cost),
    plot_progress             = FALSE,
    number_of_particles       = n.particles.gamma.cost.optimization,
    max_number_of_iterations  = n.iterations.gamma.cost.optimization,
    max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
    parameter_bounds          = parameter.bounds.gamma.cost,
    tryCall                   = TRUE, 
    verbose                   = FALSE,
    lhc_init                  = FALSE, 
    wait_complete_iteration   = TRUE,
    logfile                   = NULL,
    projectfile               = NULL,
    break_file                = "stopPPSO.txt"
    )
    
    # employ PSO proposed gamma and cost parameters
    tuned.svm <- svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                     y           = real.groups.training,
                     gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                     cost        = as.numeric(result.internal.gamma.cost$par[2]),
                     tunecontrol = tune.control(cross = number.cross.splits),
                     probability = TRUE
    )
    
    # extract best model
    tuned.svm.model <- tuned.svm
  }
  
  # prepare counts for sample classification
  normalized.counts.prediction <- normalized.counts[dgeParticle$biomarker.transcripts, evaluation.samples]
  
  # store support vectors
  save(tuned.svm.model, 
       file = paste("outputShuffled/", shuffle.iteration, "-Testing-SupportVectors.RData", sep = "")
  )
  # perform classification
  prediction.class <- predict(tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  # summarize in confusion matrix
  confusion.matrix <- table(prediction.class, real.groups.evaluation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.evaluation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2) {
     # ROC
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3],
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  # summarize results
  result <- list()
  result[["training.samples"]] <- training.samples
  result[["evaluation.samples"]] <- evaluation.samples
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  
  return(result)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.classify.validation.shuffled ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.classify.validation.shuffled <- function(dge.tool = dge,
                                                      best.particle = best.particle.input,
                                                      shuffle.iteration = i,
                                                      n.to.impute = 0, 
                                                      clin.characteristics.validation = NA,
                                                      clin.characteristics.group = NA,
                                                      clin.characteristics.specified = NA,
                                                      R.snapshot = "Pre-PSO-snapshot.RData",
                                                      verbose = TRUE){
  # Performs validation of a classification algorithm trained using shuffled group labels.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   shuffle.iteration: Numeric value to be provided by thromboSeqPSO.controls-function.
  #   n.to.impute: Numeric-value indicating the number of reads counts (0 - value) 
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R environment
  load(R.snapshot)
  
  # replace loaded dge by dge provided in function input arguments, prevents wrong or larger DGE previously compiled
  # to be employed for coming analyses
  dge <- dge.tool
  # ensure all count tables in DGE are of same size and sample composition
  dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
  dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
  
  # load particle
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  dgeParticle <- dgeTraining
  
  if (!is.na(clin.characteristics.validation)) {
    validation.samples <- c(
      colnames(dge)[dge$samples$group %in% clin.characteristics.group &
                      !colnames(dge) %in% c(training.samples, evaluation.samples) &
                      dge$samples[,clin.characteristics.validation] %in% clin.characteristics.specified],
      colnames(dge)[!dge$samples$group %in% clin.characteristics.group &
                      !colnames(dge) %in% c(training.samples, evaluation.samples)]) 
  } else {
    validation.samples <- colnames(dge)[!colnames(dge) %in% c(training.samples, evaluation.samples)]    
  }
  
  # assign groups and select samples
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.validation <- dge$samples[validation.samples, "group"]
  
  # RUVg correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge.tool, 
                                                    readout.setting = 'validation',
                                                    training.samples.set = training.samples,
                                                    evaluation.samples.set = evaluation.samples,
                                                    validation.samples.set = validation.samples,
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  
  # enable to replace counts with 0 to provided counts in the validation series
  # by the median of those in the training series.
  if (!n.to.impute %in% c("NaN", NaN)){ # omit this function when NaN is inputted
    for (assess.sample in validation.samples) { # for each sample in the validation series
      tmpA <- matrix(dgeTraining$counts[, assess.sample]) # select the counts
      sel <- which(tmpA %in% seq(0, n.to.impute)) # identify which counts have too little raw reads detected
      if (length(sel) > 1) { # if more than one selected, calculate median read counts of these genes in the training series and replace
        tmpB <- round(apply(dgeTraining$counts[sel, training.samples], 1, median))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, n.to.impute)), assess.sample] <- tmpB
      } else if (length(sel) == 1) { # if one selected, calculate median read counts of this gene in the training series and replace
        tmpB <- round(median(as.numeric(dgeTraining$ruv.counts[sel, ])))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, n.to.impute)), assess.sample] <- tmpB
      }
    }
  }
  
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE, 
                                           refColumn = dgeParticle$refSample) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, validation.samples]
  # load SVM model
  load(paste("outputShuffled/", shuffle.iteration, "-Testing-SupportVectors.RData", sep = ""))
  
  # perform classification
  prediction.class <- predict(tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.validation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.validation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC 
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3],
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Validation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, ncol(svm.summary)], 
                    ci = TRUE
    )
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.evaluation"]] <- validation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["ci.roc"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  
  return(result)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### tune.svm ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# e1071 custom svm function
# tune function implemented from e1071 package and included in SVM training (tune.svm)
# function has been adjusted to lock randomness of the internal training series selection 
# (i.e. cross split). Adjustment highlighted by 'MB: ...'.
# additional functions were included in this file to ensure proper use while analyzing
tune.svm <- function(x, y = NULL, data = NULL, degree = NULL, gamma = NULL,
                     coef0 = NULL, cost = NULL, nu = NULL, class.weights = NULL,
                     epsilon = NULL, ...) {
  call <- match.call()
  call[[1]] <- as.symbol("best.svm")
  ranges <- list(degree = degree, gamma = gamma,
                 coef0 = coef0, cost = cost, nu = nu,
                 class.weights = class.weights, epsilon = epsilon)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    ranges = NULL
  modeltmp <- if (inherits(x, "formula"))
    tune("svm", train.x = x, data = data, ranges = ranges, ...)
  else
    tune("svm", train.x = x, train.y = y, ranges = ranges, ...)
  if (!is.null(modeltmp$best.model))
    modeltmp$best.model$call <- call
  modeltmp
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### tune.control ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

tune.control <- function(random = FALSE,
                         nrepeat = 1,
                         repeat.aggregate = mean,
                         sampling = c("cross", "fix", "bootstrap"),
                         sampling.aggregate = mean,
                         sampling.dispersion = sd,
                         cross = 10,
                         fix = 2 / 3,
                         nboot = 10,
                         boot.size = 9 / 10,
                         best.model = TRUE,
                         performances = TRUE,
                         error.fun = NULL) {
  structure(list(random = random,
                 nrepeat = nrepeat,
                 repeat.aggregate = repeat.aggregate,
                 sampling = match.arg(sampling),
                 sampling.aggregate = sampling.aggregate,
                 sampling.dispersion = sampling.dispersion,
                 cross = cross,
                 fix = fix,
                 nboot = nboot,
                 boot.size = boot.size,
                 best.model = best.model,
                 performances = performances,
                 error.fun = error.fun
  ),
  class = "tune.control"
  )
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

tune <- function(method, train.x, train.y = NULL, data = list(),
                 validation.x = NULL, validation.y = NULL,
                 ranges = NULL, predict.func = predict,
                 tunecontrol = tune.control(),
                 ...
) {
  call <- match.call()
  
  ## internal helper functions
  resp <- function(formula, data) {
    
    model.response(model.frame(formula, data))
  }
  
  classAgreement <- function (tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev])) / n
    } else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m])) / n
    }
    p0
  }
  
  ## parameter handling
  if (tunecontrol$sampling == "cross")
    validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0))
    data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x))
  if (is.data.frame(train.y))
    train.y <- as.matrix(train.y)
  
  ## prepare training indices
  if (!is.null(validation.x)) tunecontrol$fix <- 1
  n <- nrow(if (useFormula) data else train.x)
  
  # MB: fix set.seed here
  getOption("myseed")
  foo <- function() {
    if (!is.null(seed <- getOption("myseed")))
      set.seed(seed)
    sample(n)
  }
  options(myseed = 1000)
  perm.ind <- foo()
  
  # perm.ind <- sample(n) # sampling - this is random
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n)
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1)
      stop(sQuote("cross"), " must be greater than 1!")
  }
  train.ind <- if (tunecontrol$sampling == "cross")
    tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
  else if (tunecontrol$sampling == "fix")
    list(perm.ind[1:trunc(n * tunecontrol$fix)])
  else ## bootstrap
    lapply(1:tunecontrol$nboot,
           function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))
  
  ## find best model
  parameters <- if (is.null(ranges))
    data.frame(dummyparameter = 0)
  else
    expand.grid(ranges)
  p <- nrow(parameters)
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1)
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random),]
  }
  model.variances <- model.errors <- c()
  
  ## - loop over all models
  for (para.set in 1:p) {
    sampling.errors <- c()
    
    ## - loop over all training samples
    for (sample in 1:length(train.ind)) {
      repeat.errors <- c()
      
      ## - repeat training `nrepeat' times
      for (reps in 1:tunecontrol$nrepeat) {
        
        ## train one model
        pars <- if (is.null(ranges))
          NULL
        else
          lapply(parameters[para.set,,drop = FALSE], unlist)
        
        
        # set more decimals instead of round at 7 digits
        options(digits = 8)
        #if(para.set==1){
        # print(list(train.x[train.ind[[sample]],],
        #            y = train.y[train.ind[[sample]]]))
        # }
        model <- if (useFormula)
          do.call(method, c(list(train.x,
                                 data = data,
                                 subset = train.ind[[sample]]),
                            pars, list(...)
          )
          )
        else
          do.call(method, c(list(train.x[train.ind[[sample]],],
                                 y = train.y[train.ind[[sample]]]),
                            pars, list(...)
          )
          )
        
        ## predict validation set
        pred <- predict.func(model,
                             if (!is.null(validation.x))
                               validation.x
                             else if (useFormula)
                               data[-train.ind[[sample]],,drop = FALSE]
                             else if (inherits(train.x, "matrix.csr"))
                               train.x[-train.ind[[sample]],]
                             else
                               train.x[-train.ind[[sample]],,drop = FALSE]
        )
        # print(pred)
        ## compute performance measure
        true.y <- if (!is.null(validation.y))
          validation.y
        else if (useFormula) {
          if (!is.null(validation.x))
            resp(train.x, validation.x)
          else
            resp(train.x, data[-train.ind[[sample]],])
        } else
          train.y[-train.ind[[sample]]]
        # print(true.y)
        # print(validation.y)
        
        if (is.null(true.y)) true.y <- rep(TRUE, length(pred))
        
        repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun))
          tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) ## classification error
          1 - classAgreement(table(pred, true.y))
        else if (is.numeric(true.y) && is.numeric(pred)) ## mean squared error
          crossprod(pred - true.y) / length(pred)
        else
          stop("Dependent variable has wrong type!")
      }
      # print(repeat.errors[reps])
      sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
    }
    model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
    model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors)
  }
  
  ## return results
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))
    NULL
  else
    lapply(parameters[best,,drop = FALSE], unlist)
  structure(list(best.parameters  = parameters[best,,drop = FALSE],
                 best.performance = model.errors[best],
                 method           = if (!is.character(method))
                   deparse(substitute(method)) else method,
                 nparcomb         = nrow(parameters),
                 train.ind        = train.ind,
                 sampling         = switch(tunecontrol$sampling,
                                           fix = "fixed training/validation set",
                                           bootstrap = "bootstrapping",
                                           cross = if (tunecontrol$cross == n) "leave-one-out" else
                                             paste(tunecontrol$cross,"-fold cross validation", sep="")
                 ),
                 performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                 best.model       = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula)
                     do.call(method, c(list(train.x, data = data),
                                       pars, list(...)))
                   else
                     do.call(method, c(list(x = train.x,
                                            y = train.y),
                                       pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }
  ),
  class = "tune"
  )
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### best.tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

best.tune <- function(...) {
  call <- match.call()
  modeltmp <- tune(...)$best.model
  modeltmp$call <- call
  modeltmp
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### print.tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

print.tune <- function(x, ...) {
  if (x$nparcomb > 1) {
    cat("\nParameter tuning of ", sQuote(x$method), ":\n\n", sep="")
    cat("- sampling method:", x$sampling,"\n\n")
    cat("- best parameters:\n")
    tmp <- x$best.parameters
    rownames(tmp) <- ""
    print(tmp)
    cat("\n- best performance:", x$best.performance, "\n")
    cat("\n")
  } else {
    cat("\nError estimation of ", sQuote(x$method), " using ", x$sampling, ": ",
        x$best.performance, "\n\n", sep="")
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### summary.tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

summary.tune <- function(object, ...)
  structure(object, class = "summary.tune")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### print.summary.tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

print.summary.tune <- function(x, ...) {
  print.tune(x)
  if (!is.null(x$performances) && (x$nparcomb > 1)) {
    cat("- Detailed performance results:\n")
    print(x$performances)
    cat("\n")
  }
}

hsv_palette <- function(h = 2/3, from = 0.7, to = 0.2, v = 1)
  function(n) hsv(h = h, s = seq(from, to, length.out = n), v = v)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### plot.tune ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

plot.tune <- function(x,
                      type=c("contour","perspective"),
                      theta=60,
                      col="lightblue",
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      swapxy = FALSE,
                      transform.x = NULL,
                      transform.y = NULL,
                      transform.z = NULL,
                      color.palette = hsv_palette(),
                      nlevels = 20,
                      ...)
{
  if (is.null(x$performances))
    stop("Object does not contain detailed performance measures!")
  k <- ncol(x$performances)
  if (k > 4) stop("Cannot visualize more than 2 parameters")
  type = match.arg(type)
  
  if (is.null(main))
    main <- paste("Performance of `", x$method, "'", sep="")
  
  if (k == 3)
    plot(x$performances[,1:2], type = "b", main = main)
  else  {
    if (!is.null(transform.x))
      x$performances[,1] <- transform.x(x$performances[,1])
    if (!is.null(transform.y))
      x$performances[,2] <- transform.y(x$performances[,2])
    if (!is.null(transform.z))
      x$performances[,3] <- transform.z(x$performances[,3])
    if (swapxy)
      x$performances[,1:2] <- x$performances[,2:1]
    x <- xtabs(error~., data = x$performances[,-k])
    if (is.null(xlab)) xlab <- names(dimnames(x))[1 + swapxy]
    if (is.null(ylab)) ylab <- names(dimnames(x))[2 - swapxy]
    if (type == "perspective")
      persp(x=as.double(rownames(x)),
            y=as.double(colnames(x)),
            z=x,
            xlab=xlab,
            ylab=ylab,
            zlab="accuracy",
            theta=theta,
            col=col,
            ticktype="detailed",
            main = main,
            ...
      )
    else
      filled.contour(x=as.double(rownames(x)),
                     y=as.double(colnames(x)),
                     xlab=xlab,
                     ylab=ylab,
                     nlevels=nlevels,
                     color.palette = color.palette,
                     main = main,
                     x, ...)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### predict.svm ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

predict.svm <-
  function (object, newdata,
            decision.values = FALSE,
            probability = FALSE,
            ...,
            na.action = na.omit)
  {
    if (missing(newdata))
      return(fitted(object))
    
    if (object$tot.nSV < 1)
      stop("Model is empty!")
    
    
    if(inherits(newdata, "Matrix")) {
      loadNamespace("SparseM")
      loadNamespace("Matrix")
      newdata <- as(newdata, "matrix.csr")
    }
    if(inherits(newdata, "simple_triplet_matrix")) {
      loadNamespace("SparseM")
      ind <- order(newdata$i, newdata$j)
      newdata <- new("matrix.csr",
                     ra = newdata$v[ind],
                     ja = newdata$j[ind],
                     ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                     dimension = c(newdata$nrow, newdata$ncol))
    }
    
    sparse <- inherits(newdata, "matrix.csr")
    if (object$sparse || sparse)
      loadNamespace("SparseM")
    
    act <- NULL
    if ((is.vector(newdata) && is.atomic(newdata)))
      newdata <- t(t(newdata))
    if (sparse)
      newdata <- SparseM::t(SparseM::t(newdata))
    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
      rownames(newdata)
    else
      1:nrow(newdata)
    if (!object$sparse) {
      if (inherits(object, "svm.formula")) {
        if(is.null(colnames(newdata)))
          colnames(newdata) <- colnames(object$SV)
        newdata <- model.matrix(delete.response(terms(object)),
                                as.data.frame(newdata))
        newdata <- na.action(newdata)
        act <- attr(newdata, "na.action")
      } else {
        newdata <- na.action(as.matrix(newdata))
        act <- attr(newdata, "na.action")
      }
    }
    
    if (!is.null(act) && !preprocessed)
      rowns <- rowns[-act]
    
    if (any(object$scaled))
      newdata[,object$scaled] <-
      scale(newdata[,object$scaled, drop = FALSE],
            center = object$x.scale$"scaled:center",
            scale  = object$x.scale$"scaled:scale"
      )
    
    if (ncol(object$SV) != ncol(newdata))
      stop ("test data does not match model !")
    
    ret <- .C ("svmpredict",
               as.integer (decision.values),
               as.integer (probability),
               
               ## model
               as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.integer (object$compprob),
               as.double  (if (object$compprob) object$probA else 0),
               as.double  (if (object$compprob) object$probB else 0),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),
               
               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),
               
               ## test matrix
               as.double  (if (sparse) newdata@ra else t(newdata)),
               as.integer (nrow(newdata)),
               as.integer (if (sparse) newdata@ia else 0),
               as.integer (if (sparse) newdata@ja else 0),
               as.integer (sparse),
               
               ## decision-values
               ret = double(nrow(newdata)),
               dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
               prob = double(nrow(newdata) * object$nclasses),
               
               PACKAGE = "e1071"
    )
    
    ret2 <- if (is.character(object$levels)) # classification: return factors
      factor (object$levels[ret$ret], levels = object$levels)
    else if (object$type == 2) # one-class-classification: return TRUE/FALSE
      ret$ret == 1
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
      ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    else
      ret$ret
    
    names(ret2) <- rowns
    ret2 <- napredict(act, ret2)
    
    if (decision.values) {
      colns = c()
      for (i in 1:(object$nclasses - 1))
        for (j in (i + 1):object$nclasses)
          colns <- c(colns,
                     paste(object$levels[object$labels[i]],
                           "/", object$levels[object$labels[j]],
                           sep = ""))
        attr(ret2, "decision.values") <-
          napredict(act,
                    matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE,
                           dimnames = list(rowns, colns)
                    )
          )
    }
    
    if (probability && object$type < 2) {
      if (!object$compprob)
        warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
      else
        attr(ret2, "probabilities") <-
          napredict(act,
                    matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE,
                           dimnames = list(rowns, object$levels[object$labels])
                    )
          )
    }
    
    ret2
  }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thromboSeq.LOOCV ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thromboSeq.LOOCV <- function(dge = dgeIncludedSamples,
                             k.variables = 3,
                             variable.to.assess = c("lib.size"),
                             variable.threshold = c(0.8),
                             ruvg.pvalue.threshold.group = 1e-2,
                             ruvg.pvalue.threshold.strongest.variable = 1e-2,
                             select.biomarker.FDR = TRUE,
                             significance.threshold = 0.05,
                             minimum.n.transcripts.biomarkerpanel = 2,
                             svm.gamma.range = 2^(-20:0),
                             svm.cost.range = 2^(0:20),
                             number.cross.splits = 2,
                             clinical.info.output = NA,
                             figureDir = "figureOutputFolder",
                             number.cores = 8,
                             verbose = TRUE) {
  # Performs leave-one-out cross-validation (LOOCV) analysis of the provided dataset.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   k.variables: Number of (expected) variables/axes in the dataset.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
  #                                be reached by the correlation between the counts and  any variable
  #                                in order to bypass the wanted variable 'group' to be selected.
  #   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
  #                                be reached by the correlation between the counts and the specific 
  #                                variable in order to this variable to be assigned to the RUVg axis.
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
  #                                         biomarker panel.
  #   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
  #   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
  #   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
  #   clinical.data.in.output: String with additional clinical info to be presented in the svm.summary.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #   
  # Returns:
  #  Result-container with classification details and metrics. Prints ROC curve and confusion matrix
  
  # load required packages
  if (verbose == TRUE){
    print("Load required packages ppso, edgeR, e1071, and RUVSeq")
  } 
  suppressMessages(library(ppso))
  suppressMessages(library(edgeR))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(e1071))
  suppressMessages(library(foreach))
  suppressMessages(library(doMC))
  
  if (verbose == TRUE){
    print("Perform LOOCV analysis")
  } 
  
  # loop all samples one by one with foreach loop.
  # select leave-out sample, perform RUV correction employing only the training samples
  # continue with ANOVA differential expression analysis of splice junctions and select
  # biomarker gene panel based on provided significance threshold.
  # start training of SVM, and predict remaining left-out sample.
  registerDoMC(cores = number.cores)
  loocv.loop <- foreach(leave.out.sample = colnames(dge)) %dopar% {
    if (verbose == TRUE) {
      print(leave.out.sample)
    }
    # assign samples
    training.samples.loocv <- colnames(dge)[leave.out.sample != colnames(dge)]
    real.groups.training <- dge$samples[training.samples.loocv, "group"]
    real.groups.prediction <- dge$samples[leave.out.sample, "group"]
    
    # ensure all count tables in DGE are of same size and sample composition
    dge$raw.counts <- dge$raw.counts[,which(colnames(dge$raw.counts) %in% colnames(dge$counts))]
    dge$ruv.counts <- dge$ruv.counts[,which(colnames(dge$ruv.counts) %in% colnames(dge$counts))]
    
    dgeTraining <- perform.RUVg.correction(dge = dge, 
                                           k.variables = k.variables, 
                                           variable.to.assess = variable.to.assess,
                                           variable.threshold = variable.threshold, 
                                           ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                           ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                           training.series.only = T,
                                           training.series.only.samples = training.samples.loocv,
                                           verbose = verbose)
    dgeTraining$counts <- dgeTraining$ruv.counts
    dgeTraining$samples$raw.lib.size <- dgeTraining$samples$lib.size
    dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
    
    # perform TMM normalization
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = FALSE) # calculate normalization factors
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
    dgeTraining$samples <- droplevels(dgeTraining$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
    
    # Likelihood-ratio test modified for thromboSeq (ANOVA)
    dgeTraining$ruv.counts <- dgeTraining$counts
    dgeTraining$counts <- dgeTraining$raw.counts
    dgeTraining$samples$lib.size <- dgeTraining$samples$raw.lib.size
    thromboSeq.anova <- anovaLikeEdgeRthromboSeq(dgeTraining[, training.samples.loocv], 
                                                 method = "TMM",
                                                 normalize.on.training.series = TRUE, 
                                                 samples.for.training = training.samples.loocv,
                                                 ref.column = dgeTraining$refSample)
    
    # select biomarker panel using either FDR or p-value statistics
    if (select.biomarker.FDR == TRUE){
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$FDR < significance.threshold & 
          thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22, "X")
        ]  
    } else {
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$PValue < significance.threshold & 
          thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22, "X")
        ]  
    }
    
    # only continue when biomarker panel size is more than provided threshold
    if (length(selected.transcripts) > minimum.n.transcripts.biomarkerpanel) {
      # second SVM-model, re-train with adjusted biomarker panel and employ a grid search for gamma and cost
      tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples.loocv]),
                            y           = real.groups.training,size,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            probability = TRUE
      )
      
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      # prepare counts for sample prediction
      normalized.counts.prediction <- normalized.counts[selected.transcripts, leave.out.sample]
      
      # perform prediction
      prediction.class <- predict(tuned.svm.model,
                                  newdata = t(normalized.counts.prediction), 
                                  probability = TRUE)
      
    } else {
      predicted.group <- NA
      biomarker.panel <- NA
      prediction.class <- c(NA, NA)
    }
    
    # summarize results
    result <- list()
    tryName <- paste("Try", leave.out.sample, sep = "")
    result[["training.samples"]] <- training.samples.loocv
    result[["prediction.sample"]] <- leave.out.sample
    result[["predicted.group"]] <- as.character(prediction.class)
    result[["real.group"]] <- as.character(real.groups.prediction)
    result[["biomarker.panel"]] <- length(selected.transcripts)
    result[["predictive.strength"]] <- attributes(prediction.class)$probabilities
    result
  } 
  
  # load additional packages
  suppressMessages(library(reshape, warn.conflicts = F, quietly = T))
  suppressMessages(library(ROCR, warn.conflicts = F, quietly = T))
  suppressMessages(library(pROC, warn.conflicts = F, quietly = T))
  suppressMessages(library(ggplot2, warn.conflicts = F, quietly = T))
  
  # summarize data into a data frame
  svm.summary <- data.frame(
    sample.ID = unlist(lapply(loocv.loop, function(x){x[["prediction.sample"]]})),
    predicted.group = unlist(lapply(loocv.loop, function(x){x[["predicted.group"]]})),
    real.group = unlist(lapply(loocv.loop, function(x){x[["real.group"]]})),
    biomarker.panel = unlist(lapply(loocv.loop, function(x){x[["biomarker.panel"]]}))
  )
  # add predictive values to the data frame
  svm.summary <- suppressWarnings(cbind(svm.summary, do.call("rbind", lapply(loocv.loop, function(x){x[["predictive.strength"]]}))))
  # warning suppressed, indicates: duplicated rownames, hence not used.
  for.roc <- ncol(svm.summary)
  
  # add clinical info to svm.summary
  if (!is.na(clinical.info.output)) {
    col.names <- colnames(svm.summary)
    svm.summary <- cbind(svm.summary, dge$samples[training.samples, clinical.info.output])
    colnames(svm.summary) <- c(col.names, clinical.info.output)
  }
  
  # summarize in confusion matrix
  confusion.matrix <- cast(svm.summary, 
                           predicted.group ~ real.group, 
                           length, 
                           value = "sample.ID")
  confusion.matrix.evaluated <- classAgreement(as.matrix(confusion.matrix), match.names = T)
  
  if (nlevels(dge$samples$group) == 2){
    # in case two classes are included, generate ROC curve
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3], 
                        label.ordering = levels(dge$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Training Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, for.roc-1], 
                    ci = TRUE
    )
  } else {
    perfa <- NA
    AUC <- confusion.matrix.evaluated$diag
    roc.optimal.accuracy <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.summary <- NA
  }
  
  # prepare and print output figures
  
  # Plot ROC-curve
  # check for correct output-directory for ROC Curve
  if( nlevels(dge$samples$group) == 2){
    currentDir <- getwd()
    if (file.exists(figureDir) == TRUE){
      pdf(paste(figureDir, "/ROCcurve_LOOCV.pdf", sep = ""))  
    } else {
      setwd('..') # move to (previous) mother directory
      if (file.exists(figureDir) == TRUE){
        pdf(paste(figureDir, "/ROCcurve_LOOCV.pdf", sep = ""))  
      } else {
        dir.create(figureDir)
        pdf(paste(figureDir, "/ROCcurve_LOOCV.pdf", sep = ""))  
      }
    }
    plot(perfa, 
         lwd = 4, 
         col = "grey")
    
    dev.off()
    setwd(currentDir)
  }
  # prepare confusion matrix to print
  confusion.matrix <- cast(svm.summary, 
                           predicted.group ~ real.group, 
                           length, 
                           value = "sample.ID")
  rownames(confusion.matrix) <- confusion.matrix$predicted.group
  confusion.matrix <- confusion.matrix[, -1]
  lev <- sort(unique(c(colnames(confusion.matrix), rownames(confusion.matrix))))
  confusion.matrix <- confusion.matrix[lev, lev]
  colnames(confusion.matrix) <- lev
  rownames(confusion.matrix) <- lev
  # confusion.matrix <- as.matrix(confusion.matrix)
  # confusion.matrix[is.na(confusion.matrix)] <- 0
  melted.confusion.matrix <- data.frame(confusion.matrix)
  melted.confusion.matrix$predicted <- rownames(melted.confusion.matrix)
  melted.confusion.matrix <- melt(melted.confusion.matrix, id.vars = "predicted")
  colnames(melted.confusion.matrix) <- c("predicted", "real", "frequency")
  
  tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
    geom_tile(aes(fill = frequency)) +
    scale_fill_continuous(low = "white", high = "red") + 
    geom_text(aes(label = frequency)) +
    ggtitle("SVM classification results") + 
    labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
    theme_bw() + 
    theme(legend.position = "top")
  
  currentDir <- getwd()
  if (file.exists(figureDir) == TRUE){
    pdf(paste(figureDir, "/TrainingSeriesConfusionMatrix.pdf", sep = ""))  
  } else {
    setwd('..') # move to (previous) mother directory
    if (file.exists(figureDir) == TRUE){
      pdf(paste(figureDir, "/TrainingSeriesConfusionMatrix.pdf", sep = ""))  
    } else {
      dir.create(figureDir)
      pdf(paste(figureDir, "/TrainingSeriesConfusionMatrix.pdf", sep = ""))  
    }
  }
  print(tiles)
  dev.off()
  
  setwd(currentDir)
  
  # summarize output
  results.classification.loocv <- list()
  results.classification.loocv$svm.summary <- svm.summary
  results.classification.loocv$confusion.matrix <- confusion.matrix
  results.classification.loocv$confusion.matrix.evaluated <- confusion.matrix.evaluated
  results.classification.loocv$AUCorDiagonal <- AUC
  results.classification.loocv$roc.95ci <- roc.95ci
  results.classification.loocv$roc.optimal.accuracy <- roc.optimal.accuracy
  results.classification.loocv$perfa <- perfa
  results.classification.loocv$roc.summary <- roc.summary
  save(results.classification.loocv, file = "results.classification.loocv.RData")
  return(results.classification.loocv)
}

# End



#### thrombo.algo.anova.up ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.anova.up <- function(x){
  # Function to perform digitalSWARM together with ppso-functions. Selects from the input variable 'x'
  # the values > or < than 0.5. Function specifically for the signature increased biomarkers 
  #
  # Args:
  #   x: vector provided by ppso-function with proposed biomarker values.
  #   
  # Returns:
  #  P-value of t-test comparing both groups
  
  x <- x[-length(x)]
  selected.transcripts <- signature_up[which(x > 0.5)]
  # ensure >1 transcript is selected, otherwise return p-value = 1
  if(length(selected.transcripts) > 1){
    tTest <- t.test(as.numeric(apply(normalized.counts[selected.transcripts, colnames(dge)[dge$samples$group == levels(dge$samples$group)[1] &
                                                                                             colnames(dge) %in% evaluation.set]],
                                     2, median)),
                    as.numeric(apply(normalized.counts[selected.transcripts, colnames(dge)[dge$samples$group == levels(dge$samples$group)[2] &
                                                                                             colnames(dge) %in% evaluation.set]],
                                     2, median)))
    return(tTest$p.value)
  } else {
    return(1)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### thrombo.algo.anova.down ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

thrombo.algo.anova.down <- function(x){
  # Function to perform digitalSWARM together with ppso-functions. Selects from the input variable 'x'
  # the values > or < than 0.5. Function specifically for the signature decreased biomarkers 
  #
  # Args:
  #   x: vector provided by ppso-function with proposed biomarker values.
  #   
  # Returns:
  #  P-value of t-test comparing both groups
  
  x <- x[-length(x)]
  selected.transcripts <- signature_down[which(x > 0.5)]
  # ensure >1 transcript is selected, otherwise return p-value = 1
  if(length(selected.transcripts) > 1){
    tTest <- t.test(as.numeric(apply(normalized.counts[selected.transcripts, colnames(dge)[dge$samples$group == levels(dge$samples$group)[1] &
                                                                                             colnames(dge) %in% evaluation.set]],
                                     2, median)),
                    as.numeric(apply(normalized.counts[selected.transcripts, colnames(dge)[dge$samples$group == levels(dge$samples$group)[1] &
                                                                                             colnames(dge) %in% evaluation.set]],
                                     2, median)))
    return(tTest$p.value)
  } else {
    return(1)
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### TEPscore ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

TEPscore <- function(signature.up = signature.up.swarm,
                     signature.down = signature.down.swarm,
                     normalized.count = normalized.counts,
                     dgeIncludedSample = dge,
                     series = NA){
  # Function calculates the TEPscore following digitalSWARM biomarker selection. First selects either groups included
  # (at maximum only two groups allowed), then calculates median expression values in case only the increased RNAs are included,
  # in case only the decreased RNAs are included, in case both increased RNAs plus the inverted decreased RNAs are included, or
  # in case both the decreased RNAs plus the inverted increased RNAs are included.
  # 
  # Args:
  #   signature.up: list containing the ppso-selected increased RNA biomarkers
  #   signature.down: list containing the ppso-selected decreased RNA biomarkers
  #   normalized.count = matrix containing the normalized counts
  #   dgeIncludedSample = DGE-object containing the employed sample counts and sample info
  #   series: either 'training.set', 'evaluation.set', 'verification.set', or 'validation.set' indicating the series
  #           employed for TEPscore calculation
  #
  # Returns:
  #  Container list containing per calculated score the scores per sample and calculated AUC-values
  # 
  
  # select either groups
  group.a <- levels(dgeIncludedSample$samples$group)[1]
  group.b <- levels(dgeIncludedSample$samples$group)[2]
  
  # test on only up
  group.a.up <- as.matrix(as.numeric(apply(normalized.count[signature.up, colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & 
                                                                                                        colnames(dgeIncludedSample) %in% series]], 2, median)))
  rownames(group.a.up) = colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & colnames(dgeIncludedSample) %in% series]
  group.b.up <- as.matrix(as.numeric(apply(normalized.count[signature.up, colnames(dgeIncludedSample)[dgeIncludedSample$samples$group== group.b & 
                                                                                                        colnames(dgeIncludedSample) %in% series]], 2, median)))
  rownames(group.b.up) = colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.b & colnames(dgeIncludedSample) %in% series]
  
  t.test.only.up <- t.test(as.numeric(group.a.up[, 1]), as.numeric(group.b.up[, 1]))
  
  veri.rocr.up <- prediction(predictions = c(as.numeric(as.character(group.a.up[, 1])), as.numeric(as.character(group.b.up[, 1]))),
                             labels = c(rep(group.a, times = nrow(group.a.up)), rep(group.b, times = nrow(group.b.up)))
  )
  perf.veri.rocr.up <- performance(veri.rocr.up, "tpr", "fpr")
  auc.veri.rocr.up <- attributes(performance(veri.rocr.up, 'auc'))$y.values[[1]]
  
  # test on only down
  group.a.down <- as.matrix(as.numeric(apply(normalized.count[signature.down, colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & 
                                                                                                            colnames(dgeIncludedSample) %in% series]], 2, median)))
  rownames(group.a.down) = colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & colnames(dgeIncludedSample) %in% series]
  group.b.down <- as.matrix(as.numeric(apply(normalized.count[signature.down, colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.b & 
                                                                                                            colnames(dgeIncludedSample) %in% series]], 2, median)))
  rownames(group.b.down) = colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.b & colnames(dgeIncludedSample) %in% series]
  
  t.test.only.down <- t.test(as.numeric(group.a.down[, 1]), as.numeric(group.b.down[, 1]))
  
  veri.rocr.down <- prediction(predictions = c(as.numeric(as.character(group.a.down[, 1])), as.numeric(as.character(group.b.down[, 1]))),
                               labels = c(rep(group.a, times = nrow(group.a.down)), rep(group.b, times = nrow(group.b.down)))
  )
  perf.veri.rocr.down <- performance(veri.rocr.down, "tpr", "fpr")
  auc.veri.rocr.down <- attributes(performance(veri.rocr.down, 'auc'))$y.values[[1]]
  
  # summarize; down will be inverted
  # calculate median RNAs with reduced spliced RNA levels
  medianDown <- median(as.numeric(apply(normalized.count[signature.down,
                                                         colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)))
  # calculate correction factor by calculating median expression value by adding twice the delta median to true value
  delta.down <- as.matrix(as.numeric(apply(normalized.count[signature.down, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)) +
                            ((medianDown - as.numeric(apply(normalized.count[signature.down, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median))) +
                               (medianDown - as.numeric(apply(normalized.count[signature.down, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)))
                            ))
  rownames(delta.down) <- colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]
  
  group.a.delta.down <- as.matrix(delta.down[(colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & colnames(dgeIncludedSample) %in% series]), ])
  group.b.delta.down <- as.matrix(delta.down[(colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.b & colnames(dgeIncludedSample) %in% series]), ])
  
  group.a.delta.down.combi <- group.a.up + group.a.delta.down
  group.a.delta.down.combi <- cbind(group.a.delta.down.combi, rep(group.a, times = nrow(group.a.delta.down.combi)))
  group.b.delta.down.combi <- group.b.up + group.b.delta.down
  group.b.delta.down.combi <- cbind(group.b.delta.down.combi, rep(group.b, times = nrow(group.b.delta.down.combi)))
  
  t.test.combi.down.inv <- t.test(as.numeric(group.a.delta.down.combi[, 1]), as.numeric(group.b.delta.down.combi[, 1]))
  
  veri.rocr.d.down <- prediction(predictions = c(as.numeric(as.character(group.a.delta.down.combi[, 1])), as.numeric(as.character(group.b.delta.down.combi[, 1]))),
                                 labels = c(rep(group.a, times = nrow(group.a.delta.down.combi)), rep(group.b, times = nrow(group.b.delta.down.combi)))
  )
  perf.veri.rocr.d.down <- performance(veri.rocr.d.down, "tpr", "fpr")
  auc.veri.rocr.d.down <- attributes(performance(veri.rocr.d.down, 'auc'))$y.values[[1]]
  
  # summarize; up will be inverted
  # calculate median RNAs with reduced spliced RNA levels
  medianUp <- median(as.numeric(apply(normalized.count[signature.up,
                                                       colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)))
  # calculate correction factor by calculating median expression value by adding twice the delta median to true value
  delta.up <- as.matrix(as.numeric(apply(normalized.count[signature.up, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)) +
                          ((medianUp - as.numeric(apply(normalized.count[signature.up, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median))) +
                             (medianUp - as.numeric(apply(normalized.count[signature.up, colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]], 2, median)))
                          ))
  rownames(delta.up) <- colnames(dgeIncludedSample)[colnames(dgeIncludedSample) %in% series]
  group.a.delta.up <- as.matrix(delta.up[(colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.a & colnames(dgeIncludedSample) %in% series]),])
  group.b.delta.down <- as.matrix(delta.up[(colnames(dgeIncludedSample)[dgeIncludedSample$samples$group == group.b & colnames(dgeIncludedSample) %in% series]),])
  
  group.a.delta.up.combi <- group.a.down + group.a.delta.up
  group.a.delta.up.combi <- cbind(group.a.delta.up.combi, rep(group.a, times = nrow(group.a.delta.up.combi)))
  group.b.delta.up.combi <- group.b.down + group.b.delta.down
  group.b.delta.up.combi <- cbind(group.b.delta.up.combi, rep(group.b, times = nrow(group.b.delta.up.combi)))
  
  t.test.combi.up.inv <- t.test(as.numeric(group.a.delta.up.combi[, 1]), as.numeric(group.b.delta.up.combi[, 1]))
  
  veri.rocr.d.up <- prediction(predictions = c(as.numeric(as.character(group.a.delta.up.combi[, 1])), as.numeric(as.character(group.b.delta.up.combi[, 1]))),
                               labels = c(rep(group.a, times = nrow(group.a.delta.up.combi)), rep(group.b, times = nrow(group.b.delta.up.combi)))
  )
  perf.veri.rocr.d.up <- performance(veri.rocr.d.up, "tpr", "fpr")
  auc.veri.rocr.d.up <- attributes(performance(veri.rocr.d.up, 'auc'))$y.values[[1]]
  
  output <- list()
  output$group.a.up <- group.a.up
  output$group.b.up <- group.b.up
  output$t.test.only.up <- t.test.only.up
  output$rocr.only.up <- veri.rocr.up
  output$perf.only.up <- perf.veri.rocr.up
  output$auc.only.up <- auc.veri.rocr.up
  output$group.a.down <- group.a.down
  output$group.b.down <- group.b.down
  output$t.test.only.down <- t.test.only.down
  output$rocr.only.down <- veri.rocr.down
  output$perf.only.down <- perf.veri.rocr.down
  output$auc.only.down <- auc.veri.rocr.down
  output$group.a.delta.up.combi <- group.a.delta.up.combi
  output$group.b.delta.up.combi <- group.b.delta.up.combi
  output$t.test.combi.up.inv <- t.test.combi.up.inv
  output$rocr.d.up <- veri.rocr.d.up
  output$perf.d.up <- perf.veri.rocr.d.up
  output$auc.d.up <- auc.veri.rocr.d.up
  output$group.a.delta.down.combi <- group.a.delta.down.combi
  output$group.b.delta.down.combi <- group.b.delta.down.combi
  output$t.test.combi.down.inv <- t.test.combi.down.inv
  output$rocr.d.down <- veri.rocr.d.down
  output$perf.d.down <- perf.veri.rocr.d.down
  output$auc.d.down <- auc.veri.rocr.d.down
  return(output)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### digitalSWARM ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

digitalSWARM <- function(dge = dgeIncludedSamples,
                         percentage.for.training = 30,
                         percentage.for.evaluation = 30,
                         percentage.for.verification = 20,
                         percentage.for.validation = 20,
                         training.samples.provided = NULL,
                         evaluation.samples.provided = NULL,
                         verification.samples.provided = NULL,
                         validation.samples.provided = NULL,
                         n.training.iterations = 200,
                         variable.to.assess = c("ageatbiopsy","lib.size"),
                         variable.to.assess.thresholds = c(0.9, 1.0, 0.7, 0.8),
                         select.biomarker.FDR = FALSE,
                         FDR.value.threshold = 0.05,
                         minimum.n.transcripts.biomarkerpanel = 2,
                         n.particles = 100,
                         n.iterations = 20,
                         number.cores = 2, 
                         verbose = TRUE){
  # Perform digitalSWARM thromboSeq classification algorithm development. 
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   percentage.for.training: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the training series.
  #   percentage.for.evaluation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the evaluation series.
  #   percentage.for.verification: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the verification series.
  #   percentage.for.validation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the validation series.
  #   training.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           training series.
  #   evaluation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           evaluation series.
  #   verification.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           verification series.
  #   validation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           validation series.
  #   n.training.iterations: Numeric value indicating number of iterations (split training-evaluation series) that 
  #                           need to be performed.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
  #                                         biomarker panel.
  #   n.particles: Numeric-value with number of PSO particles per iteration for classifier development.
  #   n.iterations: Numeric-value with number of PSO iterations in total for classifier development.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  # Outputs a matrix containing the digitalSWARM classifiers performance of training/evaluation, verification
  # and validation series, including which readout setting was selected. Also provides a printed ROC-curve and boxplots
  # of the TEPscores among the mentioned data series.
  
  if (missing(dge)) {
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (nlevels(dge$samples$group) > 2) {
    stop("Only two-group comparisons are allowed for digitalSWARM")
  }
  
  if (!percentage.for.training %in% seq(1, 100, by =1e-3)) {
    stop("percentage.for.training should be within 1 and 100")
  }
  
  if (!percentage.for.evaluation %in% seq(1, 100, by = 1e-3)) {
    stop("percentage.for.evaluation should be within 1 and 100")
  }
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) > 1) {
    print("Both percentage for training and a separate training list provided. The provided training list will be used.")
  }
  
  if (length(training.samples.provided) > 1 & length(evaluation.samples.provided) == 1) {
    stop("list of training samples provided but no evaluation samples specified. Please specify evaluation samples.")
  }
  
  if (length(variable.to.assess) > 4) {
    stop("Only at maximum four variables.to.assess are allowed")
  }
  
  if (all(variable.to.assess %in% colnames(dge$samples)) != TRUE) {
    stop("Inputted variables do not match column names of the sample info table.")
  }
  
  if (!is.numeric(minimum.n.transcripts.biomarkerpanel)) {
    stop("Provide numeric value for minimum.n.transcripts.biomarkerpanel")
  }
  
  if (!is.numeric(n.particles)) {
    stop("Provide numeric value for n.particles")
  }
  
  if (!is.numeric(n.iterations)) {
    stop("Provide numeric value for n.iterations")
  }
  
  if (verbose == TRUE) {
    print("Load required packages ppso, edgeR, and RUVSeq")
  } 
  
  # load required packages
  suppressMessages(library(ppso))
  suppressMessages(library(edgeR))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(ROCR))
  suppressMessages(library(pROC))
  suppressMessages(library(doMC))
  suppressMessages(library(foreach))
 
  # create digitalSWARM subdirectory
  if (!file.exists("digitalSWARM")) {
    dir.create("digitalSWARM", recursive = T)
  }
  # store current directory and subdirectory
  workDir_main <- getwd()
  setwd("digitalSWARM/") # transfer to subdirectory
  workDir <- getwd()
  
  # either loop multiple training-evaluation sets or run the provided sample series
  # in case loop; use n.training.iterations, in case run provided sample series 
  # set n.training.iterations to 1.
  if (length(training.samples.provided) != 0 &
      length(evaluation.samples.provided) != 0 &
      length(verification.samples.provided) != 0 &
      length(validation.samples.provided) != 0){
    n.training.iterations <- 1
  }
  
  # set verification and validation series that should be fixed,
  # if not already provided
  if (is.numeric(percentage.for.verification) & length(verification.samples.provided) == 0 &
      is.numeric(percentage.for.validation) & length(validation.samples.provided) == 0) {
    # randomly select samples for the verification and validation series
    set.seed(1000) # lock random number generator
    series.verification <- foreach(p = 1 : nlevels(dge$samples$group)) %do% {
      n.samples.verification <- round(length(which(
        dge$samples$group == levels(dge$samples$group)[p])) * (percentage.for.verification / 100)
      ) 
      
      verification.samples.subset <- sample(
        colnames(dge)[dge$samples$group == levels(dge$samples$group)[p]],
        size = n.samples.verification,
        replace = F
      )
      
      # container
      series <- list()
      series[["verification.samples.subset"]] <- verification.samples.subset
      series
    }
    verification.set <- unlist(lapply(series.verification, function(x){x[["verification.samples.subset"]]}))  
    write.csv(verification.set, "verificationSamples.csv")
    
    series.validation <- foreach(p = 1 : nlevels(dge$samples$group)) %do% {
      n.samples.validation <- round(length(which(
        dge$samples$group == levels(dge$samples$group)[p])) * (percentage.for.validation / 100)
      ) 
      
      validation.samples.subset <- sample(
        colnames(dge[, dge$samples$group == levels(dge$samples$group)[p] & !colnames(dge) %in% verification.set]),
        size = n.samples.validation,
        replace = F
      )
      
      # container
      series <- list()
      series[["validation.samples.subset"]] <- validation.samples.subset
      series
    }
    validation.set <- unlist(lapply(series.validation,  function(x){x[["validation.samples.subset"]]}))  
    write.csv(validation.set, "validationSamples.csv")
  } else {
    verification.set <- verification.samples.provided
    validation.set <- validation.samples.provided
  }
  
  if (verbose == TRUE){
    print('Start training process')
  }
  
  # start loop in which the training and evaluaton series are employed as input for ANOVA-RUV-correction,
  # ANOVA-analysis, particle-swarm optimization biomarker selection employing the digitalSWARM approach, 
  # and calculation of the TEPscore area-under-the-curve.
  registerDoMC(cores = number.cores)
  tmp <- foreach(i = 1 : n.training.iterations) %dopar% {
    # load required packages
    suppressMessages(library(ppso))
    suppressMessages(library(edgeR))
    suppressMessages(library(ROCR))
    suppressMessages(library(pROC))
    suppressMessages(library(doMC))
    suppressMessages(library(foreach))
    suppressMessages(library(RUVSeq))
    suppressMessages(require(utils))
    
    # keep second DGE-object for validation
    dgeTmp <- dge
    # randomly select sample series if not already provided
    if (is.numeric(percentage.for.training) & length(training.samples.provided) == 0) {
      # randomly select samples for the training and evaluation series
      set.seed(i) # lock random number generator
      series.training <- foreach(p = 1 : nlevels(dge$samples$group)) %do% {
        n.samples.training <- round(length(which(
          dge$samples$group == levels(dge$samples$group)[p])) * (percentage.for.training / 100)
        ) 
        
        training.samples.subset <- sample(
          colnames(dge)[dge$samples$group == levels(dge$samples$group)[p] & 
                          !colnames(dge) %in% c(verification.set, validation.set)],
          size = n.samples.training,
          replace = F
        )
        
        # container
        series <- list()
        series[["training.samples.subset"]] <- training.samples.subset
        series
      }
      training.set <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  
      # all samples not yet assigned are included in the evaluation series
      evaluation.set <- colnames(dge)[which(!colnames(dge) %in% c(training.set, verification.set, validation.set))]
    } else {
      training.set <- training.samples.provided
      evaluation.set <- evaluation.samples.provided
    }
    
    # create matrix for grid search of best value for potential confounding factors
    if (length(variable.to.assess) == 1) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1))  
    } else if (length(variable.to.assess) == 2) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1))  
    } else if (length(variable.to.assess) == 3) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1),
                                  variable3 = seq(variable.to.assess.thresholds[5],
                                                  variable.to.assess.thresholds[6], 
                                                  by = 0.1))
    } else if (length(variable.to.assess) == 4) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1),
                                  variable3 = seq(variable.to.assess.thresholds[5],
                                                  variable.to.assess.thresholds[6], 
                                                  by = 0.1),
                                  variable4 = seq(variable.to.assess.thresholds[7],
                                                  variable.to.assess.thresholds[8], 
                                                  by = 0.1))
    }
    matrix_anova <- cbind(matrix_anova, matrix(NA, ncol = 1, nrow = 1))
    colnames(matrix_anova) <- c(variable.to.assess, "lowestFDR")
    
    # run for each potential setting an ANOVA comparison and observe with which settings the lowest FDR can be achieved
    for (row in seq(1, nrow(matrix_anova), by = 1)) {
      thromboSeq.anova <- thromboSeqANOVA(dge[,training.set],
                                          variable.to.assess = variable.to.assess,
                                          variable.threshold = as.numeric(matrix_anova[row, c(1 : length(variable.to.assess))]),
                                          plot = F, 
                                          iteration = i,
                                          swarm.optimization = F,
                                          verbose = F
      )
      
      matrix_anova[row,"lowestFDR"] <- as.numeric(thromboSeq.anova$FDR[1])
    }
    
    # select setting resulting in the lowest FDR and employ for final ANOVA
    if (nrow(data.frame(matrix_anova[which(matrix_anova[, 'lowestFDR'] == min(na.omit(matrix_anova[, 'lowestFDR']))),])) > 1) {
      set.seed(1)
      top_anova <- matrix_anova[sample(which(matrix_anova[, 'lowestFDR'] == min(na.omit(matrix_anova[, 'lowestFDR']))), size = 1), ]
    } else {
      top_anova <- matrix_anova[which(matrix_anova[, 'lowestFDR'] == min(na.omit(as.numeric(matrix_anova[, 'lowestFDR'])))), ]
    }
    
    # redo the ANOVA comparison to continue with that particular RUV-corrected dataset
    thromboSeq.anova <- thromboSeqANOVA(dge[, training.set],
                                        variable.to.assess = variable.to.assess,
                                        variable.threshold = as.numeric(top_anova[, c(1 : length(variable.to.assess))]),
                                        plot = F, 
                                        iteration = i,
                                        swarm.optimization = F,
                                        verbose = F
    )
    load(paste(i, "-RUVcorrectionSettings.RData", sep = ""))
    
    # perform correction of the verification series based on the correction factors calculated in training and evaluation series
    # omit in case no stable.transcript were identified with selected thresholds.
    if (!is.na(dgeRUV$stable.transcripts)) {
      dge <- perform.RUVg.correction.validation(dge = dge[, c(training.set, evaluation.set, verification.set)], 
                                                output.particle = dgeRUV, 
                                                readout.setting = 'validation.digitalSWARM',
                                                training.samples.set = training.set, 
                                                evaluation.samples.set = evaluation.set,
                                                validation.samples.set = verification.set
      )
      dge$counts <- dge$ruv.counts
      dge$samples$lib.size <- colSums(dge$counts)
    }
    
    # apply TMM-correction factors
    dge <- calcNormFactorsThromboseq(dge,
                                     normalize.on.training.series = TRUE, 
                                     samples.for.training = training.set,
                                     ref.sample.readout = FALSE)
    dge$samples <- droplevels(dge$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dge, log = T, normalized.lib.sizes = T)
    
    # select increased and decreased biomarkers
    if (select.biomarker.FDR == FALSE){
      signature_up <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < FDR.value.threshold  &
                                                   thromboSeq.anova$logCPM > 3 &
                                                   thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                   thromboSeq.anova$logFC > 0
                                                 ]
      signature_down <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < FDR.value.threshold &
                                                     thromboSeq.anova$logCPM > 3 &
                                                     thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                     thromboSeq.anova$logFC < 0
                                                   ]
    } else {
      signature_up <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < FDR.value.threshold  &
                                                   thromboSeq.anova$logCPM > 3 &
                                                   thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                   thromboSeq.anova$logFC > 0
                                                 ]
      signature_down <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < FDR.value.threshold &
                                                     thromboSeq.anova$logCPM > 3 &
                                                     thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                     thromboSeq.anova$logFC < 0
                                                   ]
    }
    
    # in case no biomarkers are selected, skip this part
    if (length(c(signature_up, signature_down)) > minimum.n.transcripts.biomarkerpanel) {
      
      # prepare particle swarm optimization using digital selection of particle parameters (between zero and one).
      pso.parameter.bounds <- matrix(ncol = 2, nrow = length(signature_up))
      rownames(pso.parameter.bounds) <- signature_up
      pso.parameter.bounds[,1] <- c(0)
      pso.parameter.bounds[,2] <- c(1.0)
      
      set.seed(1000) # lock random number generator, required for data reproducibility
      result.up.with.signature <- optim_pso(objective_function        = thrombo.algo.anova.up,
                                            number_of_parameters      = nrow(pso.parameter.bounds),
                                            plot_progress             = FALSE,
                                            number_of_particles       = n.particles, 
                                            max_number_of_iterations  = n.iterations, 
                                            max_number_function_calls = n.particles * n.iterations,
                                            parameter_bounds          = pso.parameter.bounds,
                                            tryCall                   = TRUE,
                                            lhc_init                  = TRUE, 
                                            wait_complete_iteration   = TRUE,
                                            logfile                   = paste(i, "-ppso-up.log", sep = ""),
                                            projectfile               = paste(i, "-ppso-up.pro", sep = ""),
                                            load_projectfile          = "no"
      )
      
      pso.parameter.bounds <- matrix(ncol = 2, nrow = length(signature_down))
      rownames(pso.parameter.bounds) <- signature_down
      pso.parameter.bounds[,1] <- c(0)
      pso.parameter.bounds[,2] <- c(1.0)
      
      set.seed(1000) # lock random number generator, required for data reproducibility
      result.down.with.signature <- optim_pso(objective_function        = thrombo.algo.anova.down,
                                              number_of_parameters      = nrow(pso.parameter.bounds),
                                              plot_progress             = FALSE,
                                              number_of_particles       = n.particles, 
                                              max_number_of_iterations  = n.iterations, 
                                              max_number_function_calls = n.particles * n.iterations,
                                              parameter_bounds          = pso.parameter.bounds,
                                              tryCall                   = TRUE,
                                              lhc_init                  = TRUE, 
                                              wait_complete_iteration   = TRUE,
                                              logfile                   = paste(i, "-ppso-down.log", sep = ""),
                                              projectfile               = paste(i, "-ppso-down.pro", sep = ""),
                                              load_projectfile          = "no"
      )
      # remove temporary PSO files
      file.remove(paste(i, "-ppso-up.log", sep = ""),
                  paste(i, "-ppso-up.pro", sep = ""),
                  paste(i, "-ppso-down.log", sep = ""),
                  paste(i, "-ppso-down.pro", sep = "")
      )
      
      # select the PSO-proposed biomarker panel for both increased and decrease spliced junction reads
      signature.up.swarm <- signature_up[result.up.with.signature$par[1 : length(result.up.with.signature$par) - 1] > 0.5]
      signature.up.swarm <- signature.up.swarm[!is.na(signature.up.swarm)]
      signature.down.swarm <- signature_down[result.down.with.signature$par[1 : length(result.down.with.signature$par) - 1] > 0.5]
      signature.down.swarm <- signature.down.swarm[!is.na(signature.down.swarm)]
      
      # run TEPscore function to calculate the median RNA expression TEP-score for the individual and combined biomarker panels.
      # outputs t-test values and AUC values
      TEPscore <- TEPscore(dgeIncludedSample = dge, 
                           signature.up = signature.up.swarm,
                           signature.down = signature.down.swarm,
                           normalized.count = normalized.counts,
                           series = verification.set)
      
      # summarize and store container
      cont <- list()
      cont[["i"]] <- i
      cont[["samples.for.training"]] <- training.set
      cont[["samples.for.evaluation"]] <- evaluation.set
      cont[["samples.for.verification"]] <- verification.set
      cont[["samples.for.validation"]] <- validation.set
      cont[["ANOVAsettings"]] <- matrix_anova
      cont[["initialSignature"]] <- length(signature_up) + length(signature_down)
      cont[["signUp"]] <- signature.up.swarm
      cont[["signDown"]] <- signature.down.swarm
      cont[["finalSignature"]] <- length(signature.up.swarm) + length(signature.down.swarm)
      cont[["lowestFDR"]] <- as.numeric(thromboSeq.anova$FDR[1])
      cont[["pValueSwarmUp"]] <- as.numeric(result.up.with.signature$value)
      cont[["pValueSwarmDown"]] <- as.numeric(result.down.with.signature$value)
      cont[["t.test.only.up.p"]] <- as.numeric(TEPscore$t.test.only.up$p.value)
      cont[["t.test.only.up.delta"]] <- as.numeric(TEPscore$t.test.only.up$estimate)[1] - as.numeric(TEPscore$t.test.only.up$estimate)[2]
      cont[["t.test.only.down.p"]] <- as.numeric(TEPscore$t.test.only.down$p.value)
      cont[["t.test.only.down.delta"]] <- as.numeric(TEPscore$t.test.only.down$estimate)[1] - as.numeric(TEPscore$t.test.only.down$estimate)[2]
      cont[["t.test.combi.up.p"]] <- as.numeric(TEPscore$t.test.combi.up.inv$p.value)
      cont[["t.test.combi.up.delta"]] <- as.numeric(TEPscore$t.test.combi.up.inv$estimate)[1] - as.numeric(TEPscore$t.test.combi.up.inv$estimate)[2]
      cont[["t.test.combi.down.p"]] <- as.numeric(TEPscore$t.test.combi.down.inv$p.value)
      cont[["t.test.combi.down.delta"]] <- as.numeric(TEPscore$t.test.combi.down.inv$estimate)[1] - as.numeric(TEPscore$t.test.combi.down.inv$estimate)[2]
      cont[["auc.only.up"]] <- TEPscore$auc.only.up
      cont[["auc.only.down"]] <- TEPscore$auc.only.down
      cont[["auc.combi.up"]] <- TEPscore$auc.d.up
      cont[["auc.combi.down"]] <- TEPscore$auc.d.down
      cont
    } else {
      # summarize and store container
      cont <- list()
      cont[["i"]] <- i
      cont[["samples.for.training"]] <- training.set
      cont[["samples.for.evaluation"]] <- evaluation.set
      cont[["samples.for.verification"]] <- verification.set
      cont[["samples.for.validation"]] <- validation.set
      cont[["ANOVAsettings"]] <- matrix_anova
      cont[["initialSignature"]] <- length(signature_up) + length(signature_down)
      cont[["signUp"]] <- NA
      cont[["signDown"]] <- NA
      cont[["finalSignature"]] <- NA
      cont[["lowestFDR"]] <- as.numeric(thromboSeq.anova$FDR[1])
      cont[["pValueSwarmUp"]] <- NA
      cont[["pValueSwarmDown"]] <- NA
      cont[["t.test.only.up.p"]] <- 1.00
      cont[["t.test.only.up.delta"]] <- 1.00
      cont[["t.test.only.down.p"]] <- 1.00
      cont[["t.test.only.down.delta"]] <- 1.00
      cont[["t.test.combi.up.p"]] <- 1.00
      cont[["t.test.combi.up.delta"]] <- 1.00
      cont[["t.test.combi.down.p"]] <- 1.00
      cont[["t.test.combi.down.delta"]] <- 1.00
      cont[["auc.only.up"]] <- 1.00
      cont[["auc.only.down"]] <- 1.00
      cont[["auc.combi.up"]] <- 1.00
      cont[["auc.combi.down"]] <- 1.00
      cont
    }
    save(cont, file = paste(i, "-cont.RData", sep = ""))
    return(cont)
  }
  
  # add summarization and readout part
  output.digital.swarm <- data.frame(
    i = unlist(lapply(tmp, function(x){x[["i"]]})),
    lowest.FDR = unlist(lapply(tmp, function(x){x[["lowestFDR"]]})),
    auc.only.up = unlist(lapply(tmp, function(x){x[["auc.only.up"]]})),
    auc.only.down = unlist(lapply(tmp, function(x){x[["auc.only.down"]]})),
    auc.combi.up = unlist(lapply(tmp, function(x){x[["auc.combi.up"]]})),
    auc.combi.down = unlist(lapply(tmp, function(x){x[["auc.combi.down"]]}))
  )
  
  if (verbose == TRUE){
    print('Start validation process')
  }
  
  # from summarized loop, select the best model
  a <- output.digital.swarm[which(output.digital.swarm$auc.only.up == max(output.digital.swarm$auc.only.up)), ]
  if (nrow(a) > 1){ a <- a[1, ]}
  b <- output.digital.swarm[which(output.digital.swarm$auc.only.down == max(output.digital.swarm$auc.only.down)), ]
  if (nrow(b) > 1){ b <- b[1, ]}
  c <- output.digital.swarm[which(output.digital.swarm$auc.combi.up == max(output.digital.swarm$auc.combi.up)), ]
  if (nrow(c) > 1){ c <- c[1, ]}
  d <- output.digital.swarm[which(output.digital.swarm$auc.combi.down == max(output.digital.swarm$auc.combi.down)), ]
  if (nrow(d) > 1){ d <- d[1, ]}
  select <- c('a', 'b', 'c', 'd')[which(c(a$auc.only.up,
                                          b$auc.only.down,
                                          c$auc.combi.up,
                                          d$auc.combi.down) == max(c(a$auc.only.up,
                                                                     b$auc.only.down,
                                                                     c$auc.combi.up,
                                                                     d$auc.combi.down)))]
  if(length(select) > 1){select <- select[length(select)]}
  if(select == 'a'){
    i <- output.digital.swarm$i[which(output.digital.swarm$auc.only.up == max(output.digital.swarm$auc.only.up))]
    if (length(i) > 1){ i <- i[1]}
  } else if(select == 'b'){
    i <- output.digital.swarm[which(output.digital.swarm$auc.only.down == max(output.digital.swarm$auc.only.down))]
    if (length(i) > 1){ i <- i[1]}
  } else if(select == 'c'){
    i <- output.digital.swarm$i[which(output.digital.swarm$auc.combi.up == max(output.digital.swarm$auc.combi.up))]
    if (length(i) > 1){ i <- i[1]}
  } else if(select == 'd'){
    i <- output.digital.swarm$i[which(output.digital.swarm$auc.combi.down == max(output.digital.swarm$auc.combi.down))]
    if (length(i) > 1){ i <- i[1]}
  }
  
  load(paste(i,"-cont.RData",sep=""))
  load(paste(i,"-RUVcorrectionSettings.RData",sep=""))
  load(paste(i,"-normalized.counts.RData",sep=""))
  training.set = cont$samples.for.training
  write.csv(training.set, file = "trainingSamples.csv")
  evaluation.set = cont$samples.for.evaluation
  write.csv(evaluation.set, file = "evaluationSamples.csv")
  
  # remove those temporary iteration files that are not useful anymore
  file.remove(c(list.files(pattern = "-cont.RData"),
                list.files(pattern = "-normalized.counts.RData"),
                list.files(pattern = "RUVcorrectionSettings.RData"))[which(!c(list.files(pattern = "-cont.RData"),
                                                                              list.files(pattern = "-normalized.counts.RData"),
                                                                              list.files(pattern = "RUVcorrectionSettings.RData")) %in% 
                                                                             c(paste(i, "-cont.RData", sep = ""), paste(i, "-RUVcorrectionSettings.RData", sep = ""), 
                                                                               paste(i, "-normalized.counts.RData", sep = "")))])
  
  # start validation process. First correct the validation series according to the selected best dataset
  if (!is.na(dgeRUV$stable.transcripts)) {
    dge <- perform.RUVg.correction.validation(dge = dge[, c(training.set, evaluation.set, verification.set, validation.set)], 
                                              output.particle = dgeRUV, 
                                              readout.setting = 'validation.digitalSWARM',
                                              training.samples.set = training.set, 
                                              evaluation.samples.set = evaluation.set,
                                              validation.samples.set = c(verification.set, validation.set)
    )
    dge$counts <- dge$ruv.counts
    dge$samples$lib.size <- colSums(dge$counts)
  }
  
  # apply TMM-correction factors
  dge <- calcNormFactorsThromboseq(dge,
                                   normalize.on.training.series = TRUE, 
                                   samples.for.training = training.set,
                                   ref.sample.readout = FALSE)
  dge$samples <- droplevels(dge$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dge, 
                           log = T, 
                           normalized.lib.sizes = T)
  
  # performed readout of the TEPscore
  TEPscore.training.evaluation <- TEPscore(dgeIncludedSample = dge, 
                                           signature.up = cont$signUp, 
                                           signature.down = cont$signDown,
                                           normalized.count = normalized.counts, 
                                           series = c(training.set, evaluation.set))
  TEPscore.verification <- TEPscore(dgeIncludedSample = dge, 
                                    signature.up = cont$signUp, 
                                    signature.down = cont$signDown,
                                    normalized.count = normalized.counts, 
                                    series = verification.set)
  TEPscore.validation <- TEPscore(dgeIncludedSample = dge, 
                                  signature.up = cont$signUp, 
                                  signature.down = cont$signDown,
                                  normalized.count = normalized.counts, 
                                  series = validation.set)
  
  save(TEPscore.training.evaluation, file = "TEPscoreTrainingEvaluation.RData")
  save(TEPscore.verification, file = "TEPscoreVerification.RData")
  save(TEPscore.validation, file = "TEPscoreValidation.RData")
  
  # visualize the ROC output and boxplots
  if (select == "a"){
    pdf("digitalSWARM-ROCcurve.pdf")
    plot(TEPscore.training.evaluation$perf.only.up, lwd = 4, col = "#C0C0C0")
    par(new = T)
    plot(TEPscore.verification$perf.only.up, lwd = 4, col = "#3C66A6")
    par(new = T)
    plot(TEPscore.validation$perf.only.up, lwd = 4, col = "#B03B3D")
    dev.off()
    
    pdf("digitalSWARM-BoxplotsTEPscores.pdf")
    boxplot(as.numeric(as.character(TEPscore.training.evaluation$group.a.up[, 1])),
            as.numeric(as.character(TEPscore.training.evaluation$group.b.up[, 1])),
            as.numeric(as.character(TEPscore.verification$group.a.up[, 1])),
            as.numeric(as.character(TEPscore.verification$group.b.up[, 1])),
            as.numeric(as.character(TEPscore.validation$group.a.up[, 1])),
            as.numeric(as.character(TEPscore.validation$group.b.up[, 1])),
            col = c('#feb24c', "#f03b20"), 
            pch = 19, 
            xaxt = 'n')
    dev.off()
    
    print.matrix <- "Employed score only increased"
    # 95%-CI
    roc.95.traineval <- roc(c(rep('groupA', times = nrow(TEPscore.training.evaluation$group.a.up)), rep('groupB', times = nrow(TEPscore.training.evaluation$group.b.up))),
                            as.numeric(c(TEPscore.training.evaluation$group.a.up[,1], TEPscore.training.evaluation$group.b.up[,1])), 
                            ci = TRUE
    )
    roc.95.veri <- roc(c(rep('groupA', times = nrow(TEPscore.verification$group.a.up)), rep('groupB', times = nrow(TEPscore.verification$group.b.up))),
                       as.numeric(c(TEPscore.verification$group.a.up[,1], TEPscore.verification$group.b.up[,1])), 
                       ci = TRUE
    )
    roc.95.val <- roc(c(rep('groupA', times = nrow(TEPscore.validation$group.a.up)), rep('groupB', times = nrow(TEPscore.validation$group.b.up))),
                      as.numeric(c(TEPscore.validation$group.a.up[,1], TEPscore.validation$group.b.up[,1])), 
                      ci = TRUE
    )
    # max accuracy
    rocra <- TEPscore.training.evaluation$rocr.only.up
    max.accuracy.traineval <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                    (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.verification$rocr.only.up
    max.accuracy.verification <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                       (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.validation$rocr.only.up
    max.accuracy.validation <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                     (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    
    # output results in csv file
    matrix <- matrix(nrow = 5, ncol = 4)
    rownames(matrix) <- c("Training plus evaluation", "Verification", "Validation", "", "")
    colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy (%)")
    matrix[1, ] <- c(c(length(cont$samples.for.training) + length(cont$samples.for.evaluation)),
                     round(TEPscore.training.evaluation$auc.only.up, digits = 2),
                     paste(round(roc.95.traineval$ci[1], digits = 2), "-",
                           round(roc.95.traineval$ci[3], digits = 2), sep = ""),
                     round(max.accuracy.traineval * 100, digits = 1)
    )  
    matrix[2,] <- c(length(cont$samples.for.verification),
                    round(TEPscore.verification$auc.only.up, digits = 2),
                    paste(round(roc.95.veri$ci[1], digits = 2), "-",
                          round(roc.95.veri$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.verification * 100, digits = 1)
    )
    matrix[3,] <- c(length(cont$samples.for.validation),
                    round(TEPscore.validation$auc.only.up, digits = 2),
                    paste(round(roc.95.val$ci[1], digits = 2), "-",
                          round(roc.95.val$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.validation * 100, digits = 1)
    )
    matrix[4,] <- c(print.matrix, "", "", "")
    # store in csv file
    write.csv(matrix, file = "ROCcurve_digitalSWARM-Metrics.csv")
    
  } else if (select == "b") {
    pdf("digitalSWARM-ROCcurve.pdf")
    plot(TEPscore.training.evaluation$perf.only.down, lwd = 4, col = "#C0C0C0")
    par(new = T)
    plot(TEPscore.verification$perf.only.down, lwd = 4, col = "#3C66A6")
    par(new = T)
    plot(TEPscore.validation$perf.only.down, lwd = 4, col = "#B03B3D")
    dev.off()
    
    pdf("digitalSWARM-BoxplotsTEPscores.pdf")
    boxplot(as.numeric(as.character(TEPscore.training.evaluation$group.a.down[, 1])),
            as.numeric(as.character(TEPscore.training.evaluation$group.b.down[, 1])),
            as.numeric(as.character(TEPscore.verification$group.a.down[, 1])),
            as.numeric(as.character(TEPscore.verification$group.b.down[, 1])),
            as.numeric(as.character(TEPscore.validation$group.a.down[, 1])),
            as.numeric(as.character(TEPscore.validation$group.b.down[, 1])),
            col = c('#feb24c', "#f03b20"), 
            pch = 19, 
            xaxt = 'n')
    dev.off()
    
    print.matrix <- "Employed score only decreased"
    
    # 95%-CI
    roc.95.traineval <- roc(c(rep('groupA', times = nrow(TEPscore.training.evaluation$group.a.down)), rep('groupB', times = nrow(TEPscore.training.evaluation$group.b.down))),
                            as.numeric(c(TEPscore.training.evaluation$group.a.down[,1], TEPscore.training.evaluation$group.b.down[,1])), 
                            ci = TRUE
    )
    roc.95.veri <- roc(c(rep('groupA', times = nrow(TEPscore.verification$group.a.down)), rep('groupB', times = nrow(TEPscore.verification$group.b.down))),
                       as.numeric(c(TEPscore.verification$group.a.down[,1], TEPscore.verification$group.b.down[,1])), 
                       ci = TRUE
    )
    roc.95.val <- roc(c(rep('groupA', times = nrow(TEPscore.validation$group.a.down)), rep('groupB', times = nrow(TEPscore.validation$group.b.down))),
                      as.numeric(c(TEPscore.validation$group.a.down[,1], TEPscore.validation$group.b.down[,1])), 
                      ci = TRUE
    )
    # max accuracy
    rocra <- TEPscore.training.evaluation$rocr.only.down
    max.accuracy.traineval <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                    (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.verification$rocr.only.down
    max.accuracy.verification <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                       (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.validation$rocr.only.down
    max.accuracy.validation <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                     (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    
    # output results in csv file
    matrix <- matrix(nrow = 5, ncol = 4)
    rownames(matrix) <- c("Training plus evaluation", "Verification", "Validation", "", "")
    colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy (%)")
    matrix[1, ] <- c(c(length(cont$samples.for.training) + length(cont$samples.for.evaluation)),
                     round(TEPscore.training.evaluation$auc.only.down, digits = 2),
                     paste(round(roc.95.traineval$ci[1], digits = 2), "-",
                           round(roc.95.traineval$ci[3], digits = 2), sep = ""),
                     round(max.accuracy.traineval * 100, digits = 1)
    )  
    matrix[2,] <- c(length(cont$samples.for.verification),
                    round(TEPscore.verification$auc.only.down, digits = 2),
                    paste(round(roc.95.veri$ci[1], digits = 2), "-",
                          round(roc.95.veri$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.verification * 100, digits = 1)
    )
    matrix[3,] <- c(length(cont$samples.for.validation),
                    round(TEPscore.validation$auc.only.down, digits = 2),
                    paste(round(roc.95.val$ci[1], digits = 2), "-",
                          round(roc.95.val$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.validation * 100, digits = 1)
    )
    matrix[4,] <- c(print.matrix,"","","")
    # store in csv file
    write.csv(matrix, file = "ROCcurve_digitalSWARM-Metrics.csv")
    
  } else if (select == "c") {
    pdf("digitalSWARM-ROCcurve.pdf")
    plot(TEPscore.training.evaluation$perf.d.up, lwd = 4, col = "#C0C0C0")
    par(new=T)
    plot(TEPscore.verification$perf.d.up, lwd = 4, col = "#3C66A6")
    par(new=T)
    plot(TEPscore.validation$perf.d.up, lwd = 4, col = "#B03B3D")
    dev.off()
    
    pdf("digitalSWARM-BoxplotsTEPscores.pdf")
    boxplot(as.numeric(as.character(TEPscore.training.evaluation$group.a.delta.up.combi[, 1])),
            as.numeric(as.character(TEPscore.training.evaluation$group.b.delta.up.combi[, 1])),
            as.numeric(as.character(TEPscore.verification$group.a.delta.up.combi[, 1])),
            as.numeric(as.character(TEPscore.verification$group.b.delta.up.combi[, 1])),
            as.numeric(as.character(TEPscore.validation$group.a.delta.up.combi[, 1])),
            as.numeric(as.character(TEPscore.validation$group.b.delta.up.combi[, 1])),
            col = c('#feb24c', "#f03b20"), 
            pch = 19, 
            xaxt = 'n')
    dev.off()
    
    print.matrix <- "Employed score combined down with inverted up"
    
    # 95%-CI
    roc.95.traineval <- roc(c(rep('groupA', times = nrow(TEPscore.training.evaluation$group.a.delta.up.combi)), rep('groupB', times = nrow(TEPscore.training.evaluation$group.b.delta.up.combi))),
                            as.numeric(c(TEPscore.training.evaluation$group.a.delta.up.combi[,1], TEPscore.training.evaluation$group.b.delta.up.combi[,1])), 
                            ci = TRUE
    )
    roc.95.veri <- roc(c(rep('groupA', times = nrow(TEPscore.verification$group.a.delta.up.combi)), rep('groupB', times = nrow(TEPscore.verification$group.b.delta.up.combi))),
                       as.numeric(c(TEPscore.verification$group.a.delta.up.combi[,1], TEPscore.verification$group.b.delta.up.combi[,1])), 
                       ci = TRUE
    )
    roc.95.val <- roc(c(rep('groupA', times = nrow(TEPscore.validation$group.a.delta.up.combi)), rep('groupB', times = nrow(TEPscore.validation$group.b.delta.up.combi))),
                      as.numeric(c(TEPscore.validation$group.a.delta.up.combi[,1], TEPscore.validation$group.b.delta.up.combi[,1])), 
                      ci = TRUE
    )
    # max accuracy
    rocra <- TEPscore.training.evaluation$rocr.d.up
    max.accuracy.traineval <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                    (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.verification$rocr.d.up
    max.accuracy.verification <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                       (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.validation$rocr.d.up
    max.accuracy.validation <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                     (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    
    # output results in csv file
    matrix <- matrix(nrow = 5, ncol = 4)
    rownames(matrix) <- c("Training plus evaluation", "Verification", "Validation", "", "")
    colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy (%)")
    matrix[1, ] <- c(c(length(cont$samples.for.training) + length(cont$samples.for.evaluation)),
                     round(TEPscore.training.evaluation$auc.d.up, digits = 2),
                     paste(round(roc.95.traineval$ci[1], digits = 2), "-",
                           round(roc.95.traineval$ci[3], digits = 2), sep = ""),
                     round(max.accuracy.traineval * 100, digits = 1)
    )  
    matrix[2,] <- c(length(cont$samples.for.verification),
                    round(TEPscore.verification$auc.d.up, digits = 2),
                    paste(round(roc.95.veri$ci[1], digits = 2), "-",
                          round(roc.95.veri$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.verification * 100, digits = 1)
    )
    matrix[3,] <- c(length(cont$samples.for.validation),
                    round(TEPscore.validation$auc.d.up, digits = 2),
                    paste(round(roc.95.val$ci[1], digits = 2), "-",
                          round(roc.95.val$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.validation * 100, digits = 1)
    )
    matrix[4,] <- c(print.matrix,"","","")
    # store in csv file
    write.csv(matrix, file = "ROCcurve_digitalSWARM-Metrics.csv")
    
  } else if (select == "d") {
    pdf("digitalSWARM-ROCcurve.pdf")
    plot(TEPscore.training.evaluation$perf.d.down, lwd = 4, col = "#C0C0C0")
    par(new = T)
    plot(TEPscore.verification$perf.d.down, lwd = 4, col = "#3C66A6")
    par(new = T)
    plot(TEPscore.validation$perf.d.down, lwd = 4, col = "#B03B3D")
    dev.off()
    
    pdf("digitalSWARM-BoxplotsTEPscores.pdf")
    boxplot(as.numeric(as.character(TEPscore.training.evaluation$group.a.delta.down.combi[, 1])),
            as.numeric(as.character(TEPscore.training.evaluation$group.b.delta.down.combi[, 1])),
            as.numeric(as.character(TEPscore.verification$group.a.delta.down.combi[, 1])),
            as.numeric(as.character(TEPscore.verification$group.b.delta.down.combi[, 1])),
            as.numeric(as.character(TEPscore.validation$group.a.delta.down.combi[, 1])),
            as.numeric(as.character(TEPscore.validation$group.b.delta.down.combi[, 1])),
            col = c('#feb24c', "#f03b20"), 
            pch = 19, 
            xaxt = 'n')
    dev.off()
    
    print.matrix <- "Employed score combined up with inverted down"
    
    # 95%-CI
    roc.95.traineval <- roc(c(rep('groupA', times = nrow(TEPscore.training.evaluation$group.a.delta.down.combi)), rep('groupB', times = nrow(TEPscore.training.evaluation$group.b.delta.down.combi))),
                            as.numeric(c(TEPscore.training.evaluation$group.a.delta.down.combi[, 1], TEPscore.training.evaluation$group.b.delta.down.combi[, 1])), 
                            ci = TRUE
    )
    roc.95.veri <- roc(c(rep('groupA', times = nrow(TEPscore.verification$group.a.delta.down.combi)), rep('groupB', times = nrow(TEPscore.verification$group.b.delta.down.combi))),
                       as.numeric(c(TEPscore.verification$group.a.delta.down.combi[, 1], TEPscore.verification$group.b.delta.down.combi[, 1])), 
                       ci = TRUE
    )
    roc.95.val <- roc(c(rep('groupA', times = nrow(TEPscore.validation$group.a.delta.down.combi)), rep('groupB', times = nrow(TEPscore.validation$group.b.delta.down.combi))),
                      as.numeric(c(TEPscore.validation$group.a.delta.down.combi[, 1], TEPscore.validation$group.b.delta.down.combi[, 1])), 
                      ci = TRUE
    )
    # max accuracy
    rocra <- TEPscore.training.evaluation$rocr.d.down
    max.accuracy.traineval <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                    (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.verification$rocr.d.down
    max.accuracy.verification <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                       (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    rocra <- TEPscore.validation$rocr.d.down
    max.accuracy.validation <- max((unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
                                     (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)))
    
    # output results in csv file
    matrix <- matrix(nrow = 5, ncol = 4)
    rownames(matrix) <- c("Training plus evaluation", "Verification", "Validation", "", "")
    colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy (%)")
    matrix[1, ] <- c(c(length(cont$samples.for.training) + length(cont$samples.for.evaluation)),
                     round(TEPscore.training.evaluation$auc.d.down, digits = 2),
                     paste(round(roc.95.traineval$ci[1], digits = 2), "-",
                           round(roc.95.traineval$ci[3], digits = 2), sep = ""),
                     round(max.accuracy.traineval * 100, digits = 1)
    )  
    matrix[2,] <- c(length(cont$samples.for.verification),
                    round(TEPscore.verification$auc.d.down, digits = 2),
                    paste(round(roc.95.veri$ci[1], digits = 2), "-",
                          round(roc.95.veri$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.verification * 100, digits = 1)
    )
    matrix[3,] <- c(length(cont$samples.for.validation),
                    round(TEPscore.validation$auc.d.down, digits = 2),
                    paste(round(roc.95.val$ci[1], digits = 2), "-",
                          round(roc.95.val$ci[3], digits = 2), sep = ""),
                    round(max.accuracy.validation * 100, digits = 1)
    )
    matrix[4,] <- c(print.matrix,"","","")
    # store in csv file
    write.csv(matrix, file = "ROCcurve_digitalSWARM-Metrics.csv")
  }
  
  # store biomarker panel
  panel <- dgeIncludedSamples$genes[c(cont$signUp, cont$signDown),
                                    c('ensembl_gene_id', 'hgnc_symbol', 'description')]
  write.csv(panel, file = 'digital.swarm.biomarker.panel.csv')
  
  # store 'select'
  cont$select <- select
  contName <- list.files(pattern = 'cont.RData')
  save(cont, file = contName)
  
  return(matrix)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### digitalSWARM.shuffled ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

digitalSWARM.shuffled <- function( dge = dgeIncludedSamples,
                                   percentage.for.training = 30,
                                   percentage.for.evaluation = 30,
                                   percentage.for.verification = 20,
                                   percentage.for.validation = 20,
                                   training.samples.provided = NULL,
                                   evaluation.samples.provided = NULL,
                                   verification.samples.provided = NULL,
                                   validation.samples.provided = NULL,
                                   n.training.iterations = 200,
                                   variable.to.assess = c("ageatbiopsy","lib.size"),
                                   variable.to.assess.thresholds = c(0.9, 1.0, 0.7, 0.8),
                                   select.biomarker.FDR = FALSE,
                                   FDR.value.threshold = 0.05,
                                   minimum.n.transcripts.biomarkerpanel = 2,
                                   n.particles = 100,
                                   n.iterations = 20,
                                   n.shuffled = 1000,
                                   matrix = digitalSWARM.output,
                                   number.cores = 2, 
                                   verbose = TRUE){
  # Perform digitalSWARM thromboSeq classification algorithm development shuffled lables (specificity-analysis). 
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   percentage.for.training: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the training series.
  #   percentage.for.evaluation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the evaluation series.
  #   percentage.for.verification: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the verification series.
  #   percentage.for.validation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the validation series.
  #   training.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           training series.
  #   evaluation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           evaluation series.
  #   verification.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           verification series.
  #   validation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           validation series.
  #   n.training.iterations: Numeric value indicating number of iterations (split training-evaluation series) that 
  #                           need to be performed.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
  #                                         biomarker panel.
  #   n.particles: Numeric-value with number of PSO particles per iteration for classifier development.
  #   n.iterations: Numeric-value with number of PSO iterations in total for classifier development.
  #   n.shuffled: Numeric- value indicating number of iterations employed in shuffling mode.
  #   matrix: matrix object outputted by the digitalSWARM function
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  # Updates the matrix-object created in digitalSWARM-function by adding the shuffled labels output
  
  if (missing(dge)) {
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (nlevels(dge$samples$group) > 2) {
    stop("Only two-group comparisons are allowed for digitalSWARM")
  }
  
  if (!percentage.for.training %in% seq(1, 100, by =1e-3)) {
    stop("percentage.for.training should be within 1 and 100")
  }
  
  if (!percentage.for.evaluation %in% seq(1, 100, by = 1e-3)) {
    stop("percentage.for.evaluation should be within 1 and 100")
  }
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) > 1) {
    print("Both percentage for training and a separate training list provided. The provided training list will be used.")
  }
  
  if (length(training.samples.provided) > 1 & length(evaluation.samples.provided) == 1) {
    stop("list of training samples provided but no evaluation samples specified. Please specify evaluation samples.")
  }
  
  if (length(variable.to.assess) > 4) {
    stop("Only at maximum four variables.to.assess are allowed")
  }
  
  if (all(variable.to.assess %in% colnames(dge$samples)) != TRUE) {
    stop("Inputted variables do not match column names of the sample info table.")
  }
  
  if (!is.numeric(minimum.n.transcripts.biomarkerpanel)) {
    stop("Provide numeric value for minimum.n.transcripts.biomarkerpanel")
  }
  
  if (!is.numeric(n.particles)) {
    stop("Provide numeric value for n.particles")
  }
  
  if (!is.numeric(n.iterations)) {
    stop("Provide numeric value for n.iterations")
  }
  
  if (verbose == TRUE) {
    print("Load required packages ppso, edgeR, and RUVSeq")
  } 
  
  # load required packages
  suppressMessages(library(ppso))
  suppressMessages(library(edgeR))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(ROCR))
  suppressMessages(library(pROC))
  suppressMessages(library(doMC))
  suppressMessages(library(foreach))
  
  # find output files of digitalSWARM, search in current or next directory
  if (getwd() == "digitalSWARM/"){
    toLoad <- c(list.files(pattern = "counts.RData"),
                list.files(pattern = "Settings.RData"),
                list.files(pattern = "cont.RData"))
    load(toLoad[1])
    load(toLoad[2])
    load(toLoad[3])
  } else if (file.exists("digitalSWARM/")) {
    setwd('digitalSWARM/')
    toLoad <- c(list.files(pattern = "counts.RData"),
                list.files(pattern = "Settings.RData"),
                list.files(pattern = "cont.RData"))
    load(toLoad[1])
    load(toLoad[2])
    load(toLoad[3])
  } else if (!identical(list.files(pattern = "cont.RData"),character(0))) {
    toLoad <- c(list.files(pattern = "counts.RData"),
                list.files(pattern = "Settings.RData"),
                list.files(pattern = "cont.RData"))
    load(toLoad[1])
    load(toLoad[2])
    load(toLoad[3])
  } else {
    stop("Ensure this function is started in a folder with the output of digitalSWARM")
  }
  select = cont$select
  
  # create digitalSWARM subdirectory
  if (!file.exists("digitalSWARMshuffled")) {
    dir.create("digitalSWARMshuffled", recursive = T)
  }
  # store current directory and subdirectory
  workDir_main <- getwd()
  setwd("digitalSWARMshuffled/") # transfer to subdirectory
  workDir <- getwd()
  
  training.set = cont$samples.for.training
  evaluation.set = cont$samples.for.evaluation
  verification.set = cont$samples.for.verification
  validation.set = cont$samples.for.validation
  
  if (verbose == TRUE){
    print('Start shuffled mode')
  }
  
  registerDoMC(cores = number.cores)
  output.shuffled <- foreach(i = 1 : n.shuffled) %dopar% {
    dgeTmp <- dge
    
    # shuffle labels training and evaluation series
    set.seed(i)
    groups <- dge$samples[c(training.set, evaluation.set),]$group
    dge$samples[c(training.set, evaluation.set),]$group <- sample(dge$samples[c(training.set, evaluation.set),]$group)
    # ensure that a new shuffle is done when either training or evaluation set has only PD or PR selected by chance. 
    if (length(unique(dge$samples[c(training.set), ]$group)) < 2){
      dge$samples[c(training.set, evaluation.set), ]$group <- sample(dge$samples[c(training.set, evaluation.set),]$group)  
    } else if ((length(unique(dge$samples[c(evaluation.set),]$group)) < 2)) {
      dge$samples[c(training.set, evaluation.set), ]$group <- sample(dge$samples[c(training.set, evaluation.set),]$group)  
    }
    # make sure that when after shuffling >90% of the sample groups is the same, reshuffle.
    while (as.numeric(summary(groups == dge$samples[c(training.set, evaluation.set), ]$group)[3]) / length(groups) > 0.80) {
      groups <- dge$samples[c(training.set, evaluation.set), ]$group
      dge$samples[c(training.set, evaluation.set), ]$group <- sample(dge$samples[c(training.set, evaluation.set), ]$group)
      # ensure that a new shuffle is done when either training or evaluation set has only PD or PR selected by chance. 
      if (length(unique(dge$samples[c(training.set), ]$group)) < 2){
        dge$samples[c(training.set, evaluation.set), ]$group <- sample(dge$samples[c(training.set, evaluation.set), ]$group)  
      } else if ((length(unique(dge$samples[c(evaluation.set), ]$group)) < 2)) {
        dge$samples[c(training.set, evaluation.set), ]$group <- sample(dge$samples[c(training.set, evaluation.set), ]$group)  
      } 
    }
    
    # create matrix for grid search of best value for potential confounding factors
    if (length(variable.to.assess) == 1) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1))  
    } else if (length(variable.to.assess) == 2) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1))  
    } else if (length(variable.to.assess) == 3) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1),
                                  variable3 = seq(variable.to.assess.thresholds[5],
                                                  variable.to.assess.thresholds[6], 
                                                  by = 0.1))
    } else if (length(variable.to.assess) == 4) {
      matrix_anova <- expand.grid(variable1 = seq(variable.to.assess.thresholds[1],
                                                  variable.to.assess.thresholds[2], 
                                                  by = 0.1),
                                  variable2 = seq(variable.to.assess.thresholds[3],
                                                  variable.to.assess.thresholds[4], 
                                                  by = 0.1),
                                  variable3 = seq(variable.to.assess.thresholds[5],
                                                  variable.to.assess.thresholds[6], 
                                                  by = 0.1),
                                  variable4 = seq(variable.to.assess.thresholds[7],
                                                  variable.to.assess.thresholds[8], 
                                                  by = 0.1))
    }
    matrix_anova <- cbind(matrix_anova, matrix(NA, ncol = 1, nrow = 1))
    colnames(matrix_anova) <- c(variable.to.assess, "lowestFDR")
    
    # run for each potential setting an ANOVA comparison and observe with which settings the lowest FDR can be achieved
    for (row in seq(1, nrow(matrix_anova), by = 1)) {
      thromboSeq.anova <- thromboSeqANOVA(dge[, training.set],
                                          variable.to.assess = variable.to.assess,
                                          variable.threshold = as.numeric(matrix_anova[row, c(1 : length(variable.to.assess))]),
                                          plot = F, 
                                          iteration = i,
                                          swarm.optimization = F,
                                          verbose = F
      )
      
      matrix_anova[row,"lowestFDR"] <- as.numeric(thromboSeq.anova$FDR[1])
    }
    
    # check lowest FDR and use those lib-hosp settings for final ANOVA
    if (nrow(data.frame(matrix_anova[which(matrix_anova[, 'lowestFDR'] == min(na.omit(matrix_anova[, 'lowestFDR']))),])) > 1) {
      set.seed(1)
      top_anova <- matrix_anova[sample(which(matrix_anova[, 'lowestFDR'] == min(na.omit(matrix_anova[, 'lowestFDR']))), size = 1), ]
    } else {
      top_anova <- matrix_anova[which(matrix_anova[, 'lowestFDR'] == min(na.omit(as.numeric(matrix_anova[, 'lowestFDR'])))), ]
    }
    
    # redo the ANOVA comparison to continue with that particular corrected dataset
    suppressMessages(library(RUVSeq))
    thromboSeq.anova <- thromboSeqANOVA(dge[, training.set],
                                        variable.to.assess = variable.to.assess,
                                        variable.threshold = as.numeric(top_anova[, c(1 : length(variable.to.assess))]),
                                        plot = F, 
                                        iteration = i,
                                        swarm.optimization = F,
                                        verbose = F
    )
    load(paste(i, "-RUVcorrectionSettings.RData", sep = ""))
    
    # perform correction of the verification series based on the correction factors calculated in training and evaluation series
    # omit in case no stable.transcript were identified with selected thresholds.
    if (!is.na(dgeRUV$stable.transcripts)) {
      dge <- perform.RUVg.correction.validation(dge = dge[, c(training.set, evaluation.set, verification.set)], 
                                                output.particle = dgeRUV, 
                                                readout.setting = 'validation.digitalSWARM',
                                                training.samples.set = training.set, 
                                                evaluation.samples.set = evaluation.set,
                                                validation.samples.set = verification.set
      )
      dge$counts <- dge$ruv.counts
      dge$samples$lib.size <- colSums(dge$counts)
    }
    
    # apply TMM-correction factors
    dge <- calcNormFactorsThromboseq(dge,
                                     normalize.on.training.series = TRUE, 
                                     samples.for.training = training.set,
                                     ref.sample.readout = FALSE)
    dge$samples <- droplevels(dge$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dge, log = T, normalized.lib.sizes = T)
    
    # select increased and decreased biomarkers
    if (select.biomarker.FDR == FALSE){
      signature_up <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < FDR.value.threshold  &
                                                   thromboSeq.anova$logCPM > 3 &
                                                   thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                   thromboSeq.anova$logFC > 0
                                                 ]
      signature_down <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < FDR.value.threshold &
                                                     thromboSeq.anova$logCPM > 3 &
                                                     thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                     thromboSeq.anova$logFC < 0
                                                   ]
    } else {
      signature_up <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < FDR.value.threshold  &
                                                   thromboSeq.anova$logCPM > 3 &
                                                   thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                   thromboSeq.anova$logFC > 0
                                                 ]
      signature_down <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < FDR.value.threshold &
                                                     thromboSeq.anova$logCPM > 3 &
                                                     thromboSeq.anova$chromosome_name %in% c(1:22,"X") &
                                                     thromboSeq.anova$logFC < 0
                                                   ]
    }
    
    # prepare particle swarm optimization using digital selection of particle parameters (between zero and one).
    pso.parameter.bounds <- matrix(ncol = 2, nrow = length(signature_up))
    rownames(pso.parameter.bounds) <- signature_up
    pso.parameter.bounds[,1] <- c(0)
    pso.parameter.bounds[,2] <- c(1.0)
    pso.parameter.bounds <- rbind(pso.parameter.bounds, c(i - 0.005, i))
    rownames(pso.parameter.bounds) <- c(rownames(pso.parameter.bounds)[1 : nrow(pso.parameter.bounds) - 1], "iter")
    
    set.seed(1000) # lock random number generator, required for data reproducibility
    result.up.with.signature <- optim_pso(objective_function        = thrombo.algo.anova.up,
                                          number_of_parameters      = nrow(pso.parameter.bounds),
                                          plot_progress             = FALSE,
                                          number_of_particles       = n.particles, 
                                          max_number_of_iterations  = n.iterations, 
                                          max_number_function_calls = n.particles * n.iterations,
                                          parameter_bounds          = pso.parameter.bounds,
                                          tryCall                   = TRUE,
                                          lhc_init                  = TRUE, 
                                          wait_complete_iteration   = TRUE,
                                          logfile                   = paste(i, "-ppso-up.log", sep = ""),
                                          projectfile               = paste(i, "-ppso-up.pro", sep = ""),
                                          load_projectfile          = "no"
    )
    
    pso.parameter.bounds <- matrix(ncol = 2, nrow = length(signature_down))
    rownames(pso.parameter.bounds) <- signature_down
    pso.parameter.bounds[,1] <- c(0)
    pso.parameter.bounds[,2] <- c(1.0)
    pso.parameter.bounds <- rbind(pso.parameter.bounds, c(i - 0.005, i))
    rownames(pso.parameter.bounds) <- c(rownames(pso.parameter.bounds)[1 : nrow(pso.parameter.bounds)-1], "iter")
    
    set.seed(1000) # lock random number generator, required for data reproducibility
    result.down.with.signature <- optim_pso(objective_function        = thrombo.algo.anova.down,
                                            number_of_parameters      = nrow(pso.parameter.bounds),
                                            plot_progress             = FALSE,
                                            number_of_particles       = n.particles, 
                                            max_number_of_iterations  = n.iterations, 
                                            max_number_function_calls = n.particles * n.iterations,
                                            parameter_bounds          = pso.parameter.bounds,
                                            tryCall                   = TRUE,
                                            lhc_init                  = TRUE, 
                                            wait_complete_iteration   = TRUE,
                                            logfile                   = paste(i, "-ppso-down.log", sep = ""),
                                            projectfile               = paste(i, "-ppso-down.pro", sep = ""),
                                            load_projectfile          = "no"
    )
    # remove temporary PSO files
    file.remove(paste(i, "-ppso-up.log", sep = ""),
                paste(i, "-ppso-up.pro", sep = ""),
                paste(i, "-ppso-down.log", sep = ""),
                paste(i, "-ppso-down.pro", sep = "")
    )
    
    # select the PSO-proposed biomarker panel for both increased and decrease spliced junction reads
    signature.up.swarm <- signature_up[result.up.with.signature$par[1 : length(result.up.with.signature$par) - 1] > 0.5]
    signature.up.swarm <- signature.up.swarm[!is.na(signature.up.swarm)]
    signature.down.swarm <- signature_down[result.down.with.signature$par[1 : length(result.down.with.signature$par) - 1] > 0.5]
    signature.down.swarm <- signature.down.swarm[!is.na(signature.down.swarm)]
    
    # verification and validation
    dge <- dgeTmp
    if (!is.na(dgeRUV$stable.transcripts)) {
      dge <- perform.RUVg.correction.validation(dge = dge[, c(training.set, evaluation.set, verification.set, validation.set)], 
                                                output.particle = dgeRUV, 
                                                readout.setting = 'validation.digitalSWARM',
                                                training.samples.set = training.set, 
                                                evaluation.samples.set = evaluation.set,
                                                validation.samples.set = c(verification.set, validation.set)
      )
      dge$counts <- dge$ruv.counts
      dge$samples$lib.size <- colSums(dge$counts)
    }
    
    # apply TMM-correction factors
    dge <- calcNormFactorsThromboseq(dge,
                                     normalize.on.training.series = TRUE, 
                                     samples.for.training = training.set,
                                     ref.sample.readout = FALSE)
    dge$samples <- droplevels(dge$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dge, 
                             log = T, 
                             normalized.lib.sizes = T)
    
    TEPscore.training.evaluation.shuffled <- TEPscore(dgeIncludedSample = dge, 
                                                      signature.up = signature.up.swarm,
                                                      signature.down = signature.down.swarm,
                                                      normalized.count = normalized.counts,
                                                      series = c(training.set, evaluation.set))
    TEPscore.verification.shuffled <- TEPscore(dgeIncludedSample = dge, 
                                               signature.up = signature.up.swarm,
                                               signature.down = signature.down.swarm,
                                               normalized.count = normalized.counts,
                                               series = verification.set)
    TEPscore.validation.shuffled <- TEPscore(dgeIncludedSample = dge, 
                                             signature.up = signature.up.swarm,
                                             signature.down = signature.down.swarm,
                                             normalized.count = normalized.counts,
                                             series = validation.set)
    if (select == 'a') {
      output.train.eval.auc <- TEPscore.training.evaluation.shuffled$auc.only.up
      output.verification.auc <- TEPscore.verification.shuffled$auc.only.up
      output.validation.auc <- TEPscore.validation.shuffled$auc.only.up
    } else if (select == 'b') {
      output.train.eval.auc <- TEPscore.training.evaluation.shuffled$auc.only.down
      output.verification.auc <- TEPscore.verification.shuffled$auc.only.down
      output.validation.auc   <- TEPscore.validation.shuffled$auc.only.down
    } else if (select == 'c') {
      output.train.eval.auc <- TEPscore.training.evaluation.shuffled$auc.d.up
      output.verification.auc <- TEPscore.verification.shuffled$auc.d.up
      output.validation.auc   <- TEPscore.validation.shuffled$auc.d.up
    } else if (select == 'd') {
      output.train.eval.auc <- TEPscore.training.evaluation.shuffled$auc.d.down
      output.verification.auc <- TEPscore.verification.shuffled$auc.d.down
      output.validation.auc <- TEPscore.validation.shuffled$auc.d.down
    }
    
    # cont
    conti <- list()
    conti[["i"]] <- i
    conti[["output.train.eval.auc"]] <- output.train.eval.auc
    conti[["output.verification.auc"]] <- output.verification.auc
    conti[["output.validation.auc"]] <- output.validation.auc
    conti
    
    return(conti)
  }
  setwd('..')
  
  # summarize output of shuffled loop
  output_shuffled <- data.frame(
    i = unlist(lapply(output.shuffled, function(x){x[["i"]]})),
    output.train.eval.auc = unlist(lapply(output.shuffled, function(x){x[["output.train.eval.auc"]]})),
    output.verification.auc = unlist(lapply(output.shuffled, function(x){x[["output.verification.auc"]]})),
    output.validation.auc = unlist(lapply(output.shuffled, function(x){x[["output.validation.auc"]]}))
  )
  shuffled.output <- max(output_shuffled$output.validation.auc)
  matrix[5,] <- c(paste("Shuffled maximum AUC validation: ", shuffled.output, sep = ""),"","","")
  write.csv(matrix, file = "ROCcurve_digitalSWARM-Metrics.csv")
  
  return(matrix)
}
