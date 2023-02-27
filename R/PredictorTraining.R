
#' Train GLM using H2O framework
#'
#' @param trainExpr Expression dataframe with genes in rows and samples in columns. Rownames must be gene symbols
#' @param trainAge Vector of chronological ages. The order must be that same as the samples in trainExpr.
#' @param trainSamples Vector of sample names or indices to select for training. If NULL, all samples are used.
#' @param legacy Fix lambda parameter for model training
#'
#' @return A trained H2O GLM model
#' @export
#'
trainModel <- function(trainExpr,
                       trainAge,
                       trainSamples = NULL,
                       legacy = FALSE
){
  all_procs <- setdiff(unique(MultiTIMER::mbotc$ProcessName),"Background genes")
  genes_to_procs <- lapply(all_procs,function(x){unique(MultiTIMER::mbotc$Symbol[MultiTIMER::mbotc$ProcessName == x])})
  names(genes_to_procs) <- all_procs
  genes_to_procs <- lapply(genes_to_procs,function(x){intersect(rownames(trainExpr),x)})
  genes_to_procs <- genes_to_procs[unname(which(sapply(genes_to_procs,length) >= 3))]
  level_to_procs <- sapply(names(genes_to_procs),function(x){unique(MultiTIMER::mbotc$ProcessLevel[MultiTIMER::mbotc$ProcessName == x])})
  levels_to_consider <- c(1)
  level_to_procs <- level_to_procs[level_to_procs %in% levels_to_consider]

  processes_to_consider <- names(level_to_procs)

  genesInLevel <- intersect(unname(do.call(c,genes_to_procs[processes_to_consider])),rownames(trainExpr))

  if(is.null(trainSamples)){
    trainSamples <- colnames(trainExpr)
  }

  trainData <- trainExpr[genesInLevel,]
  trainData <- trainData[,trainSamples]

  trainData <- t(trainData)
  trainData <- cbind(trainAge,trainData)
  colnames(trainData)[1] <- "Age"

  require(h2o)
  conn <- h2o.init(max_mem_size="8G")
  trainData_h2o <- as.h2o(trainData)

  if(legacy){
    lambda = 0.2599206
    lambda_search = FALSE
  }else{
    lambda = NULL
    lambda_search = TRUE
  }

  model <- h2o.glm(y = "Age",
                   training_frame = trainData_h2o,
                   nfolds = 10,
                   fold_assignment = "Random",
                   family = "AUTO",
                   link = "family_default",
                   lambda_search = lambda_search,
                   lambda = lambda,
                   standardize = T,
                   alpha = 0.5,
                   seed = 100,
                   max_active_predictors = ncol(trainData_h2o),
                   solver = "IRLSM"
  )
  return(model)
}

#' Predict age of samples
#'
#' @param model The H2O model used for prediction.
#' @param predExpr Dataframe of samples used for prediction.
#' Genes should be in rows and samples in columns. All genes used
#' in model should be present in the data.
#'
#' @return Vector of predicted ages
#' @export
#'
predictAge <- function(model,
                       predExpr
){
  predExpr <- predExpr[model@model$names,]
  predAges <- c()
  maxiter <- ceiling(ncol(predExpr)/1000)
  for(i in 1:maxiter){
    from <- ((i-1)*1000+1)
    to <- min(i*1000,ncol(predExpr))

    inp_h2o <- as.h2o(t(predExpr[,from:to]))

    pred <- h2o.predict(model,inp_h2o)
    pred <- as.data.frame(pred)$predict
    names(pred) <- colnames(predExpr)[from:to]

    predAges <- c(predAges,pred)
  }
  return(predAges)
}

#' Get activity matrix of processes
#'
#' @param model The H2O model used for prediction.
#' @param expr Dataframe of samples used for activity matrix calculation.
#' Genes should be in rows and samples in columns. All genes used
#' in model should be present in the data.
#'
#' @return A matrix of process activities
#' @export
#'
getActivityMatrix <- function(model,
                              expr
){
  all_procs <- setdiff(unique(MultiTIMER::mbotc$ProcessName),"Background genes")
  genes_to_procs <- lapply(all_procs,function(x){unique(MultiTIMER::mbotc$Symbol[MultiTIMER::mbotc$ProcessName == x])})
  names(genes_to_procs) <- all_procs
  genes_to_procs <- lapply(genes_to_procs,function(x){intersect(rownames(expr),x)})
  genes_to_procs <- genes_to_procs[unname(which(sapply(genes_to_procs,length) >= 3))]

  #Process activity
  coeff <- model@model$coefficients_table
  coeff_vec <- coeff$coefficients
  names(coeff_vec) <- coeff$names

  mat <- expr

  process_activity <- list()
  for(p in processes){
    print(p)
    genes_in_proc <- genes_to_procs[[p]]
    coeffs_model <- coeff_vec[genes_in_proc]
    coeffs_model <- coeffs_model[!is.na(coeffs_model)]


    tmp <- as.data.frame((t(mat[names(coeffs_model),]) %*% coeffs_model))
    colnames(tmp)[1] <- p

    process_activity[[p]] <- tmp
  }

  activitymat <- do.call("cbind",process_activity)

  return(activitymat)

}



