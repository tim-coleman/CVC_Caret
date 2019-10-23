######################
# The purpose of this document is to combine the
# caret package with the cross validation with confidence
# procedure proposed in Lei (2019).
######################

library(caret)
library(dplyr)


### Don't run the examples! ###
run <- FALSE

### First a function that runs the actual cross validation, and outputs the actual CV values.
caret_CV <- function(x, y, V = 5, models = c("rf"), ...){
  
  # Input: "x" - covariate data frame
  #        "y" - training responses
  #        "V" - number of folds to split the data into, for use in V-fold CV
  #         "models" - a list of caret model tags, to be used in the final result
  #         ... - additional arguments to be passed to caret::train
  
  
  ### Generating the Folds ###
  folds <- createFolds(y = y, k = V, returnTrain = T)
  
  ### Specifying our Train Control object ###
  ### Tunes with RMSE for regression, kappa for classification ###
  trC <- trainControl(index = folds,  savePredictions = "all")
  
  if(is.numeric(y)){
    metric <- "RMSE"
  }
  else{
    metric <- "Kappa"
  }

  fit_f <- function(model){
    obj <- train(x = x, y = y, method = model, preProcess = c("center", "scale"), trControl = trC, metric = metric, ...)
    out <- obj[["pred"]]
    out[["model"]] <- model
    
    parameters <- names(obj[["bestTune"]])
    
    for(p in parameters){
      if(is.numeric(out[[p]])){
        out[[p]] <- round(out[[p]], 3)
      }
      p_descriptive <- paste(p, out[[p]], sep = "=")
      out[[p]] <- p_descriptive
    }
    
    out <- out %>% tidyr::unite(col = model_tune_param, c("model", parameters), sep = "__")
    return(out)
  }
  
  model_out <- do.call(rbind, lapply(models, FUN = fit_f))
  
  #return(resamples(model_out, metric = metric))
  return(model_out)
}

### Quick example
if(run){
data(iris)
TrainData <- iris[,1:4]
TrainClasses <- iris[,5]

set.seed(1)
rf_svm_iris <- caret_CV(x = TrainData, y = rnorm(nrow(iris)), V = 5, models = c("ranger", "svmRadial"))
}

### Now a function that takes in a resampled table and conducts the wild bootstrap procedure from Lei (2019)
gauss_multiplier_boot <- function(preds, B = 250, alpha = 0.05){
  
  # Input: "preds" - data frame returned by caret_CV, containing out of sample predictions for each model/tuning parameter setting
  #        "B" - number of bootstraps
  #        "metric" - a string encoding which metric should be used to conduct the inference - not used right now
  
  ### Updating to include loss ###
  preds <- preds %>% dplyr::mutate(loss = (pred-obs)^2)
  
  ### Calculating the V-fold CV error ###
  model_CV <- preds %>% dplyr::group_by(model_tune_param) %>% summarise(VCV_err = mean(loss))
  
  ### Creating a difference matrix in loss for each observation ###
  
  get_model_diffs <- function(ind){
    preds_ind <- preds %>% dplyr::filter(rowIndex == ind)
    diff_loss <- outer(preds_ind[["loss"]], preds_ind[["loss"]], "-")
    #names(diff_loss) <- unique(preds[["model_tune_param"]])
    row.names(diff_loss) <- unique(preds[["model_tune_param"]])
    return(diff_loss)
  }
  
  n <- max(preds[["rowIndex"]]) # sample size
  model_diffs <- lapply(1:n, get_model_diffs)
  
  ### Now associating each observation with each fold ###
  
  fold_library <- preds %>% dplyr::select(c(rowIndex, Resample)) %>% unique() %>% 
    dplyr::arrange(rowIndex)
  V <- length(unique(fold_library[["Resample"]])) # number of folds
  

  ### Now a function for conducting the bootstrap for a single model m, returning the p-value ###
  
  boot_model_m <- function(m){
    model_m_diffs <- data.frame(do.call(rbind, lapply(model_diffs, FUN = function(x) return(x[which(row.names(x) == m),]))))
    names(model_m_diffs) <- unique(preds[["model_tune_param"]])
    
    model_m_diffs <- cbind(fold_library, model_m_diffs)
    #mu_mj_v <- model_m_diffs %>% dplyr::group_by(Resample) %>% dplyr::summarise_all(list(mean))
    
    ### Centered differences
    cntr <- function(x) x - (V/n)*sum(x)
    model_m_diffs_centered <- model_m_diffs %>% dplyr::group_by(Resample) %>% mutate_at(vars(-group_cols()),list(cntr))
    
    ### Mean differences
    mu_mj <- model_m_diffs %>% dplyr::select(-c(rowIndex, Resample)) %>% colMeans()
    
    ### sd of differences
    ### these differences are sad because the former differences were mean...
    sd_mj <- model_m_diffs_centered %>% dplyr::ungroup() %>% dplyr::select(-c(rowIndex, Resample)) %>%  dplyr::summarise_all(list(sd))
    
    
    ### Scaled & Centered differences, needed for the bootstrapping phase
    scl <- function(X) X/sd(X)
    model_m_diffs_sc <- model_m_diffs_centered %>% dplyr::mutate_at(vars(-dplyr::group_cols()),list(scl))
    
    ### Our test statistic - need to remove the comparisons of model m with itself here!
    T_m <- max(sqrt(n)*(mu_mj/sd_mj)[which(sd_mj > 0)])
    
    
    ### Now! The bootstrap phase...
    T_bm <- rep(NA, B)
    for(b in 1:B){
      zeta <- rnorm(n, mean = 0, sd = 1)
      model_m_zeta <- model_m_diffs_sc %>% dplyr::ungroup() %>%  dplyr::mutate_if(is.numeric, list(function(x) x*zeta))
      T_bj <- model_m_zeta %>% dplyr::summarise_if(is.numeric, list(function(x) n^(-0.5)*sum(x)))
      T_bm[b] <- max(T_bj[c(which(!is.na(T_bj)))][-1]) # removing the first entry, corresponding to the nonsensical (at this point) rowindex
    }
    
    ### Calculating the p-value
    p_out <- mean(T_bm > T_m)
    return(p_out)
  }
  # Getting the p-values
  out_all_models <- data.frame("model_tune_param" = unique(preds[["model_tune_param"]]), 
                               "VCV_Pval" = sapply(unique(preds[["model_tune_param"]]), boot_model_m))
  # Merging with the actual CV errors from earlier
  suppressWarnings(out <- inner_join(out_all_models, model_CV) %>% dplyr::mutate(in_ACV = VCV_Pval > alpha))
  return(out)
}


#### Finally, the main function that implements the enire procedure! ####
CVC_full <- function(x, y, V = 5, models = c("rf"), B = 250, alpha = 0.05, ...){
  
  ## First, training the models
  CV_obj <- caret_CV(x = x, y = y, V = V, models = models, ...)
  
  ## Now applying the VCV procedure
  VCV <- gauss_multiplier_boot(preds = CV_obj, B = B, alpha = alpha)
  
  ## Returning the results
  return(VCV)
}


#### Quick example
if(run){
set.seed(1)
CVC_full(x = TrainData, y = sin(pi*TrainData[,1]*TrainData[,3]) + rnorm(nrow(TrainData), sd = 15), V = 5, models = c("ranger", "svmRadial", "glmnet"))
}