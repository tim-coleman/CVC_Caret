######################
# The purpose of this document is to combine the
# caret package with the cross validation with confidence
# procedure proposed in Lei (2019).
######################

library(caret)
library(dplyr)


### Don't run the examples! ###
run <- F

### First a function that runs the actual cross validation, and outputs the actual CV values.
caret_CV <- function(x, y, V = 5, models = c("rf"), SS_data = NULL, SS_response = NULL, tune_params = "default", ...){
  
  # Input: "x" - covariate data.frame
  #        "y" - training responses
  #        "V" - number of folds to split the data into, for use in V-fold CV
  #         "models" - a list of caret model tags, to be used in the final result
  #         "tune_params"  - a list with names corresponding to models. Each named element contains the tune_grid parameters 
  #                          accepted by the corresponding model tag. Accepts the argument "default" which lets caret use its
  #                          default tuning grid.
  #         "SS_data" - Sample splitted data.frame to be used for validation. The "V" argument is then ignored. 
  #                     Currently, if this is non-null, then tune_params need to be provided.
  #         "SS_response" - character string indicating the name of the response variable in SS_data
  #         ... - additional arguments to be passed to caret::train
  
  
  ### Generating the Folds ###
  folds <- createFolds(y = y, k = V, returnTrain = T)
  
  ### Specifying our Train Control object ###
  ### Tunes with RMSE for regression, kappa for classification ###
  if(is.null(SS_data)){
    trC <- trainControl(index = folds,  savePredictions = "all")
  }
  else{
    trC <- trainControl(method = "none")
  }

  
  if(is.numeric(y)){
    metric <- "RMSE"
  }
  else{
    metric <- "Kappa"
  }
  

  fit_f <- function(model, tune_params = "default"){
    if(all(tune_params == "default")){
      obj <- train(x = x, y = y, method = model, preProcess = c("center", "scale"), 
                   trControl = trC, metric = metric, ...)
    }
    else{
      obj <- train(x = x, y = y, method = model, preProcess = c("center", "scale"), 
                   trControl = trC, metric = metric, tuneGrid = tune_params, ...)
    }
    
    if(is.null(SS_data)){
      out <- obj[["pred"]]
    }
    else{
      out <- data.frame(pred = predict(obj, newdata = SS_data), obs = SS_data[[SS_response]],
                        tune_params, rowIndex = 1:nrow(SS_data))
    }
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
  
  if(tune_params == "default"){
    model_out <- do.call(rbind, lapply(models, FUN = fit_f))
  }
  else{
    # Storing model name in each variable
    for(i in 1:length(tune_params)){
      tune_params[[i]][["model"]] <- names(tune_params)[i]
    }
    tune_params <- lapply(tune_params, FUN = function(x) split(x, seq(nrow(x))))
    tune_params <- do.call(c, tune_params)
    model_out <- do.call(rbind, lapply(tune_params, 
                                       FUN = function(x) fit_f(model = x$model, 
                                                               tune_params =  data.frame(as.list(x[-which(names(x) == "model")])))))
  }
  return(model_out)
}

### Quick example, where response is independent from features
if(run){
  data(iris)
  TrainData <- iris[,1:4]
  TrainClasses <- iris[,5]
  
  set.seed(1)
  rf_iris_SS <- caret_CV(x = TrainData, y = sample(unique(iris$Species), nrow(iris), T), V = 5, models = c("rf", "ranger"), SS_data = iris, SS_response = "Species",
                          tune_params = list(rf = expand.grid(mtry = 1:3), 
                                             ranger = expand.grid(mtry = 1:2, splitrule = "gini", min.node.size = 1:5)))
  rf_iris_noSS <- caret_CV(x = TrainData, y =  sample(unique(iris$Species), nrow(iris), T), V = 5, models = c("rf", "ranger"), 
                           tune_params = list(rf = expand.grid(mtry = 1:3)))
}

### Now a function that takes in a resampled table and conducts the wild bootstrap procedure from Lei (2019)
gauss_multiplier_boot <- function(preds, B = 250, alpha = 0.05, screen = T){
  
  # Input: "preds" - data frame returned by caret_CV, containing out of sample predictions for each model/tuning parameter setting
  #        "B" - number of bootstraps
  #        "alpha" - significance threshold - if all p-values from model m are above alpha, model m is included
  #        "screen" - Logical. Should the screening of "obviously inferior" models be performed, before the boostrap? 
  
  ### Updating to include loss ###
  preds <- preds %>% dplyr::mutate(loss = (pred-obs)^2)
  
  ### Calculating the V-fold CV error ###
  model_CV <- preds %>% dplyr::group_by(model_tune_param) %>% dplyr::summarise(VCV_err = mean(loss))
  
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
  if(!is.null(preds[["Resample"]])){
    fold_library <- preds %>% dplyr::select(c(rowIndex, Resample)) %>% unique() %>% 
      dplyr::arrange(rowIndex)
    V <- length(unique(fold_library[["Resample"]])) # number of folds
    M <- length(unique(preds[["model_tune_param"]]))
  }
  else{
    fold_library <- preds %>% dplyr::select(c(rowIndex)) %>% unique() %>% 
      dplyr::arrange(rowIndex)
    V <- 0 #0 because we haven't done CV
    M <- length(unique(preds[["model_tune_param"]]))
  }
  ### Now a function for conducting the bootstrap for a single model m, returning the p-value. Screening is implemented in here. ###
  
  boot_model_m <- function(m){
    model_m_diffs <- data.frame(do.call(rbind, lapply(model_diffs, FUN = function(x) return(x[which(row.names(x) == m),]))))
    names(model_m_diffs) <- unique(preds[["model_tune_param"]])
    model_m_diffs <- cbind(fold_library, model_m_diffs)

    if(V > 0){
      #mu_mj_v <- model_m_diffs %>% dplyr::group_by(Resample) %>% dplyr::summarise_all(list(mean))
      
      ### Centered differences
      cntr <- function(x) x - (V/n)*sum(x)
      model_m_diffs_centered <- model_m_diffs %>% dplyr::group_by(Resample) %>% mutate_at(vars(-group_cols()),list(cntr))
      
      ### Mean differences
      mu_mj <- model_m_diffs %>% dplyr::select(-c(rowIndex, Resample)) %>% colMeans()
      
      ### sd of differences
      sd_mj <- model_m_diffs_centered %>% dplyr::ungroup() %>% dplyr::select(-c(rowIndex, Resample)) %>%  dplyr::summarise_all(list(sd))
      
      
      ### Scaled & Centered differences, needed for the bootstrapping phase
      scl <- function(X) X/sd(X)
      model_m_diffs_sc <- model_m_diffs_centered %>% dplyr::mutate_at(vars(-dplyr::group_cols()),list(scl))
      
      ### Applying the screening, using a factor of 10 in keeping with recommendation of 
      if(screen){
        thresh <- -2 * (qnorm(1 - alpha/(5*M - 5)))/sqrt(1 - (qnorm(1-alpha/(5*M - 5))^2/n))
        #print(thresh)
        screened_in <- which(sqrt(n)*(mu_mj/sd_mj) > thresh)
        #print(min(sqrt(n)*(mu_mj/sd_mj)[!is.na(mu_mj/sd_mj)]))
        compare_inds <- intersect(screened_in, which(sd_mj >0))
        #print(length(screened_in))
      } else {
        compare_inds <- which(sd_mj > 0)
      }
      
      ### Now restricting our scaled differences matrix
      model_m_diffs_sc_screen <- model_m_diffs_sc[,c(1:2, 2+compare_inds)]
      
      ### Our test statistic - need to remove the comparisons of model m with itself here!
      T_m <- max(sqrt(n)*(mu_mj/sd_mj)[compare_inds])
      
      ### Now! The bootstrap phase...
      T_bm <- rep(NA, B)
      for(b in 1:B){
        zeta <- rnorm(n, mean = 0, sd = 1)
        model_m_zeta <- model_m_diffs_sc_screen %>% dplyr::ungroup() %>%  dplyr::mutate_if(is.numeric, list(function(x) x*zeta))
        T_bj <- model_m_zeta %>% dplyr::summarise_if(is.numeric, list(function(x) n^(-0.5)*sum(x)))
        T_bm[b] <- max(T_bj[c(which(!is.na(T_bj)))][-1]) # removing the first entry, corresponding to the nonsensical (at this point) rowindex
      }
      
      ### Calculating the p-value
      p_out <- mean(T_bm > T_m)
      cat("Model", m, "p-value:", p_out, "\n")
      return(p_out)
    }
    else{
      #mu_mj_v <- model_m_diffs %>% dplyr::group_by(Resample) %>% dplyr::summarise_all(list(mean))
      
      ### Centered differences
      cntr <- function(x) x - mean(x)
      model_m_diffs_centered <- model_m_diffs %>% dplyr::mutate_all(.funs = cntr)
      
      ### Mean differences
      mu_mj <- model_m_diffs %>% dplyr::select(-c(rowIndex)) %>% colMeans()
      
      ### sd of differences
      sd_mj <- model_m_diffs_centered %>% dplyr::select(-c(rowIndex)) %>%  dplyr::summarise_all(list(sd))
      
      
      ### Scaled & Centered differences, needed for the bootstrapping phase
      scl <- function(X) X/sd(X)
      model_m_diffs_sc <- model_m_diffs_centered %>% dplyr::mutate_at(vars(-dplyr::group_cols()),list(scl))
      
      ### Applying the screening, using a factor of 10 in keeping with recommendation of 
      print("made it to flag 0")
      
      if(screen){
        thresh <- -2 * (qnorm(1 - alpha/(5*M - 5)))/sqrt(1 - (qnorm(1-alpha/(5*M - 5))^2/n))
        #print(thresh)
        screened_in <- which(sqrt(n)*(mu_mj/sd_mj) > thresh)
        #print(min(sqrt(n)*(mu_mj/sd_mj)[!is.na(mu_mj/sd_mj)]))
        compare_inds <- intersect(screened_in, which(sd_mj >0))
        #print(length(screened_in))
      } else {
        compare_inds <- which(sd_mj > 0)
      }
      

      ### Now restricting our scaled differences matrix
      model_m_diffs_sc_screen <- model_m_diffs_sc[,c(1, 1+compare_inds)]
      

      ### Our test statistic - need to remove the comparisons of model m with itself here!
      if(length(compare_inds) == 0){
        p_out <- 0
        cat("Model", m, "p-value:", p_out, "\n")
        return(p_out)
      } else {
        T_m <- max(sqrt(n)*(mu_mj/sd_mj)[compare_inds])
        
        ### Now! The bootstrap phase...
        T_bm <- rep(NA, B)
        for(b in 1:B){
          zeta <- rnorm(n, mean = 0, sd = 1)
          model_m_zeta <- model_m_diffs_sc_screen %>% dplyr::ungroup() %>%  dplyr::mutate_if(is.numeric, list(function(x) x*zeta))
          T_bj <- model_m_zeta %>% dplyr::summarise_if(is.numeric, list(function(x) n^(-0.5)*sum(x)))
          T_bm[b] <- max(T_bj[c(which(!is.na(T_bj)))][-1]) # removing the first entry, corresponding to the nonsensical (at this point) rowindex
        }
        
        ### Calculating the p-value
        p_out <- mean(T_bm > T_m)
        cat("Model", m, "p-value:", p_out, "\n")
        return(p_out)
    }
    }
  }
  # Getting the p-values
  out_all_models <- data.frame("model_tune_param" = unique(preds[["model_tune_param"]]), 
                               "VCV_Pval" = sapply(unique(preds[["model_tune_param"]]), boot_model_m))
  cat("Made it to before the last line\n")
  # Merging with the actual CV errors from earlier
  suppressWarnings(out <- dplyr::inner_join(out_all_models, model_CV) %>% dplyr::mutate(in_ACV = VCV_Pval > alpha))
  return(out)
}


#### Finally, the main function that implements the entire procedure! ####
CVC_full <- function(x, y, V = 5, models = c("rf"), SS_data = NULL, tune_params = "default", B = 250, alpha = 0.05, screen = T, ...){
  
  ## First, training the models
  CV_obj <- caret_CV(x = x, y = y, V = V, models = models, SS_data = SS_data, tune_params = tune_params, ...)
  cat("Finished Model Training\n")
  
  ## Now applying the VCV procedure
  VCV <- gauss_multiplier_boot(preds = CV_obj, B = B, alpha = alpha, screen = screen)
  
  ## Returning the results
  return(VCV)
}


#### Quick example
if(run){
  set.seed(1)
  TrainData <- data.frame(replicate(15, runif(850)))
  TestData <- data.frame(replicate(15, runif(200))) %>% dplyr::mutate(y = sin(pi*sqrt(X1*X3)) + rnorm(200, sd = 0.5))
  CVC_noSS <- CVC_full(x = TrainData, y = sin(pi*sqrt(TrainData[,1]*TrainData[,3])) + rnorm(nrow(TrainData), sd = 0.5), 
           V = 5, models = c("ranger", "rf"))
  CVC_SS <-  CVC_full(x = TrainData, y = sin(pi*sqrt(TrainData[,1]*TrainData[,3])) + rnorm(nrow(TrainData), sd = 0.5), 
                      SS_data = TestData,  models = c("rf", "ranger"), SS_response = "y",
                      tune_params = list(rf = expand.grid(mtry = c(2, 8, 15)), 
                                         ranger = expand.grid(mtry = c(2, 8, 15), splitrule = "variance", min.node.size = c(1, 5))))
}
