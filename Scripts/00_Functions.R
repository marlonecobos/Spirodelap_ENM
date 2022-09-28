# Download worldclim variables version 2.1
#
# argument options
# period: "historical", "future"
# res: 10m, 5m, 2.5m, 30s (30s only for historical)
# time: 2021-2040, 2041-2060, 2061-2080, 2081-2100 (only for future) 
# SSP: 126, 245, 370, 585 (only for future) 
# GCM: BCC-CSM2-MR, CNRM-CM6-1, CNRM-ESM2-1, CanESM5, IPSL-CM6A-LR, MIROC-ES2L, 
#      MIROC6, MRI-ESM2-0 (only for future) 
# output_dir: the directory the user wants to put the results in 

get_NWC_bio <- function(period, res, time = NULL, SSP = NULL, GCM = NULL, 
                        output_dir = NULL) {
  if (is.null(output_dir)) {dir <- getwd()} else {dir <- output_dir}
  if (length(res) > 1) {stop("Argument 'res' must be of length 1.")}
  
  if (period[1] == "historical") {
    if (any(!res %in% c("10m", "5m", "2.5m", "30s"))) {
      stop("Argument 'res' must by any of: 10m, 5m, 2.5m, or 30s")
    }
    
    url <- paste0("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_", 
                  res, "_bio.zip")
    
    dir.create(dir)
    dfile <- paste0(dir, "/wc2.1_", res, "_bio.zip")
    
    download.file(url, destfile = dfile, method = "auto", 
                  quiet = FALSE, mode = "wb", cacheOK = TRUE)
    
    dfol <- paste0(dir, "/wc2.1_", res, "_bio")
    dir.create(dfol)
    unzip(zipfile = dfile, exdir = dfol)
    
    files <- list.files(dfol, pattern = ".tif$", full.names = TRUE)
    return(raster::stack(files))
    
  } else {
    if (any(!res %in% c("10m", "5m", "2.5m"))) {
      stop("Argument 'res' must by any of: 10m, 5m, or 2.5m")
    }
    
    lens <- c(length(time), length(SSP), length(GCM))
    if (any(lens > 1)) {
      stop("Arguments 'time', 'SSP', and 'GCM' must be of length 1.")
    }
    
    url <- paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/", 
                  res, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_",
                  time, ".zip")
    
    dir.create(dir)
    dfile <- paste0(dir, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_",
                    time, ".zip")
    
    download.file(url, destfile = dfile, method = "auto", 
                  quiet = FALSE, mode = "wb", cacheOK = TRUE)
    
    dfol <- paste0(dir, "/wc2.1_", res, "_bioc_", GCM, "_ssp", SSP, "_", time)
    dir.create(dfol)
    unzip(zipfile = dfile, exdir = dfol)
    
    files <- list.files(dfol, pattern = ".tif$", full.names = TRUE)
    return(raster::stack(files))
  }
}


citation_NWC_bio <- function(period) {
  if (period[1] == "historical") {
    cat("Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: New 1km spatial resolution climate surfaces\n\tfor global land areas. International Journal of Climatology 37 (12): 4302-4315.")
  } else {
    cat("For more information about climatic data derived from the CMIP6 project, please visit:\n\thttps://wcrp-cmip.github.io/CMIP6_CVs/docs/CMIP6_institution_id.html\n\thttp://bit.ly/2gBCuqM")
  }
}




plot_list_histMO <- function(lits_histMO, col_M = "gray65", col_O = "gray25", 
                             alpha = 0.4, par_mfrow, where_ylab, par_cex = 1, 
                             par_mar = c(4.3, 4.3, 0.5, 0.3), lines = TRUE) {
  mcol <- alpha(col_M, alpha); ocol <- alpha(col_O, alpha)
  vnames <- gsub("_", " ", names(lits_histMO))
  
  par(mfrow = par_mfrow, mar = par_mar)
  par(cex = par_cex)
  
  ps <- sapply(1:length(lits_histMO), function(y) {
    yla <- ifelse(y %in% where_ylab, "Frequency", "")
    plot(lits_histMO[[y]]$hgM, col = mcol, main = "", xlab = vnames[y],
         border = mcol, freq = T, ylab = yla) 
    plot(lits_histMO[[y]]$hgO, col = ocol, add = TRUE, border = ocol)
    if (lines == TRUE) {
      abline(v = lits_histMO[[y]]$meanM, col = col_M)
      abline(v = lits_histMO[[y]]$meanO, col = col_O)
      abline(v = lits_histMO[[y]]$clM, col = col_M, lty = 2)
      abline(v = lits_histMO[[y]]$clO, col = col_O, lty = 2)
    }
    box(bty = "l")
  })
}




# from dismo
window <- function(x)  { 
  lng <- length(x)
  x <- c(x,  x[1:3])
  m <- matrix(ncol=3, nrow=lng)
  for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
  apply(m, MARGIN=1, FUN=sum)
}



# functions to explore espace
f_1 <- function(data, mapping, subset = NULL) {
  ggplot(mapping = mapping) + geom_point(data = data[subset[[1]], 
  ], colour = "lightblue1", alpha = 0.3) + geom_point(data = data[subset[[2]], 
  ], colour = "grey25", shape = 16, alpha = 0.6)
}

f_2 <- function(data, mapping, subset = NULL) {
  ggplot(mapping = mapping) + stat_density_2d(data = data[subset[[1]], 
  ], aes(fill = (..level..)), geom = "polygon", show.legend = FALSE, 
  bins = 350) + scale_fill_continuous(low = "lightblue1", 
                                      high = "dodgerblue4") + stat_density_2d(data = data[subset[[2]], 
                                      ], colour = "black", lwd = 0.3, show.legend = FALSE, 
                                      bins = 10)
}


# GLM calibration
library(foreach)
glm_calibration <- function(all_data, train_data, test_data, background, 
                            variable_sets = "all_comb", min_number = 2,
                            response_types = "all", family = binomial(link = "logit"), 
                            weights_1_0 = c(1, 10000), allowed_omission = 0.05,
                            random_percent = 50, iterations = 500, 
                            selection = "OR_AIC", output_directory = NULL, 
                            parallel = FALSE, n_cores = NULL) {
  
  # threshold
  thr <- allowed_omission * 100
  
  # data for fitting models
  ## excluding data from background (just in case)
  exa <- which(paste(background[, 1], background[, 2]) %in% 
                 paste(all_data[, 1], all_data[, 2]))
  extr <- which(paste(background[, 1], background[, 2]) %in% 
                  paste(train_data[, 1], train_data[, 2]))
  
  ## response variable
  zerosa <- rep(0, nrow(background[-exa, ]))
  zerostr <- rep(0, nrow(background[-extr, ]))
  onesa <- rep(1, nrow(all_data))
  onestr <- rep(1, nrow(train_data))
  
  response_variable <- "Response"
  
  # weights
  weightsa <- c(rep(weights_1_0[1], nrow(all_data)), 
                rep(weights_1_0[2], nrow(background[-exa, ])))
  weightstr <- c(rep(weights_1_0[1], nrow(train_data)), 
                 rep(weights_1_0[2], nrow(background[-extr, ])))
  
  ## all data
  databa <- cbind(Response = c(onesa, zerosa), rbind(all_data, background[-exa, ]))
  databtr <- cbind(Response = c(onestr, zerostr), rbind(train_data, background[-extr, ]))
  
  # formulas according to parameterization
  ## preparing sets of predictors
  if (variable_sets == "all_comb") {
    predictor_sets <- kuenm::all_var_comb(var.names = colnames(all_data)[-(1:2)],
                                          min.number = min_number)
  } else {
    predictor_sets <- variable_sets
  }
  
  ## preparing formulas
  formulas <- get_formulas(response_variable, predictor_sets, response_types)
  
  message("A total of ", length(formulas), " glm models will be calibrated.")
  
  # model calibration
  if (parallel == TRUE) {
    ## preparing parallel running
    n_cores <- ifelse(is.null(n_cores), parallel::detectCores() - 1, n_cores)
    cl <- snow::makeSOCKcluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    
    ## progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(formulas), style = 3)
    progress <- function(n) {utils::setTxtProgressBar(pb, n)}
    opts <- list(progress = progress)
    
    ## processing
    all_res <- foreach::foreach(i = 1:length(formulas), .inorder = FALSE, 
                                .options.snow = opts, .combine = "rbind") %dopar% {
                                  ## fitting models
                                  mod <- glm(as.formula(formulas[i]), family = family, 
                                             data = databa, weights = weightsa)
                                  modtr <- glm(as.formula(formulas[i]), family = family, 
                                               data = databtr, weights = weightstr)
                                  
                                  ## predicting probabilities
                                  testres <- predict.glm(modtr, newdata = test_data, 
                                                         type = "response")
                                  trainres <- predict.glm(modtr, newdata = train_data, 
                                                          type = "response")
                                  allres <- predict.glm(modtr, newdata = databa, 
                                                        type = "response")
                                  
                                  ## evaluating model
                                  ### partial ROC
                                  procs <- kuenm::kuenm_proc(occ.test = testres * weights_1_0[2], 
                                                             model = allres * weights_1_0[2], 
                                                             threshold = thr, 
                                                             rand.percent = random_percent, 
                                                             iterations = iterations)
                                  
                                  ### omission rates
                                  nval <- ceiling(nrow(train_data) * allowed_omission)
                                  threshold <- sort(trainres)[nval]
                                  orate <- sum(testres < threshold) / length(testres)
                                  
                                  ### AIC
                                  aic <- data.frame(AIC = mod$aic, delta_AIC = NA_real_, 
                                                    weight_AIC = NA_real_)
                                  
                                  ## results of evaluation
                                  return(data.frame(formulas[i], 
                                                    t(as.data.frame(procs$pROC_summary)), 
                                                    orate, aic, row.names = NULL))
                                }
    
    snow::stopCluster(cl)
    
  } else {
    ## progress bar
    pb <- utils::txtProgressBar(min = 1, max = length(formulas), style = 3)
    
    ## processing
    all_res <- list()
    
    for (x in 1:length(formulas)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, x)

      ## fitting models
      mod <- glm(formula = as.formula(formulas[x]), family = family, data = databa, 
                 weights = weightsa)
      modtr <- glm(formula = as.formula(formulas[x]), family = family, data = databtr, 
                   weights = weightstr)
      
      ## predicting probabilities
      testres <- predict.glm(modtr, newdata = test_data[, -(1:2)], 
                             type = "response")
      trainres <- predict.glm(modtr, newdata = train_data[, -(1:2)], 
                              type = "response")
      allres <- predict.glm(modtr, newdata = databa[, -(1:3)], 
                            type = "response")
      
      ## maximum value to scale for partial ROC
      allmax <- max(allres)
      
      ## evaluating model
      ### partial ROC
      procs <- kuenm::kuenm_proc(occ.test = testres / allmax, 
                                 model = allres / allmax, 
                                 threshold = thr, rand.percent = random_percent, 
                                 iterations = iterations)
      
      ### omission rates
      nval <- ceiling(nrow(train_data) * allowed_omission)
      threshold <- sort(trainres)[nval]
      orate <- sum(testres < threshold) / length(testres)
      
      ### AIC
      aic <- data.frame(AIC = mod$aic, delta_AIC = NA_real_, weight_AIC = NA_real_)
      
      ## results of evaluation
      all_res[[x]] <- data.frame(formulas[x], t(as.data.frame(procs$pROC_summary)), 
                                 orate, aic, row.names = NULL)
    }
    
    all_res <- do.call(rbind, all_res)
  }
 
  colnames(all_res) <- c("Formula", paste0("Mean_AUC_ratio_at_", thr, "%"),
                         "pval_pROC", paste0("Omission_rate_at_", thr, "%"),
                         "AIC", "delta_AIC", "weight_AIC")
  
  # final calculations
  all_res$delta_AIC <- all_res$AIC - min(all_res$AIC, na.rm = TRUE)
  all_res$weight_AIC <- exp(-0.5 * all_res$delta_AIC) / 
    sum(exp(-0.5 * all_res$delta_AIC), na.rm = TRUE)
  
  # selected models
  selection <- gsub("AIC", "AICc", selection)
  list_res <- kuenm::summary_calibration(all_res, selection = selection)
  
  # writing if needed
  if (!is.null(output_directory)) {
    dir.create(output_directory)
    name <- paste0(output_directory, "/glm_calibration_results.csv")
    name0 <- paste0(output_directory, "/glm_calibration_stats.csv")
    name1 <- paste0(output_directory, "/glm_selected_models.csv")
    write.csv(list_res[[3]], file = name, row.names = FALSE)
    write.csv(list_res[[1]], file = name0, row.names = FALSE)
    write.csv(list_res[[2]], file = name1, row.names = FALSE)
  }
  
  return(list_res)
}



# helper to get all formulas for glms based on certain criteria
get_formulas <- function(response_variable, predictor_sets, response_types = "all") {
  if (response_types[1] == "all") {
    response_types <- c("l", "q", "p", "lq", "lp", "qp", "lqp")
  }
  
  predictors <- lapply(response_types, function(x) {
    if (x == "l") {
      rs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          paste(y, collapse = " + ")
        } else {
          y
        }
      })
    }
    
    if (x == "q") {
      rs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          pr <- vapply(y, FUN.VALUE = character(1), FUN = function(z) {
            paste0("I(", z, "^2)")
            
          })
          paste(pr, collapse = " + ")
        } else {
          paste0("I(", y, "^2)")
        }
      })
    }
    
    if (x == "p") {
      rs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          cb <- combn(y, m = 2)
          cb <- vapply(1:ncol(cb), FUN.VALUE = character(1), FUN = function(z) {
            paste(cb[, z], collapse = ":")
            
          })
          paste(cb, collapse = " + ")
        } else {
          NA_character_
        }
      })
      rs <- na.omit(rs)
    }
    
    if (x == "lq") {
      ls <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          paste(y, collapse = " + ")
        } else {
          y
        }
      })
      
      qs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          pr <- vapply(y, FUN.VALUE = character(1), FUN = function(z) {
            paste0("I(", z, "^2)")
            
          })
          paste(pr, collapse = " + ")
        } else {
          paste0("I(", y, "^2)")
        }
      })
      
      rs <- paste(ls, qs, sep = " + ")
    }

    if (x == "lp") {
      rs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          paste(y, collapse = " * ")
        } else {
          y
        }
      })
    }
    
    if (x == "qp") {
      ps <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          cb <- combn(y, m = 2)
          cb <- vapply(1:ncol(cb), FUN.VALUE = character(1), FUN = function(z) {
            paste(cb[, z], collapse = ":")
            
          })
          paste(cb, collapse = " + ")
        } else {
          "erase_this_oneSB"
        }
      })
      
      qs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          pr <- vapply(y, FUN.VALUE = character(1), FUN = function(z) {
            paste0("I(", z, "^2)")
            
          })
          paste(pr, collapse = " + ")
        } else {
          paste0("I(", y, "^2)")
        }
      })
      
      rs <- paste(qs, ps, sep = " + ")
      rs <- gsub(" + erase_this_oneSB", "", rs, fixed = T)
    }
    
    if (x == "lqp") {
      ls <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          paste(y, collapse = " + ")
        } else {
          y
        }
      })
      
      ps <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          cb <- combn(y, m = 2)
          cb <- vapply(1:ncol(cb), FUN.VALUE = character(1), FUN = function(z) {
            paste(cb[, z], collapse = ":")
            
          })
          paste(cb, collapse = " + ")
        } else {
          "erase_this_oneSB"
        }
      })
      
      qs <- vapply(predictor_sets, FUN.VALUE = character(1), FUN = function(y) {
        if (length(y) > 1) {
          pr <- vapply(y, FUN.VALUE = character(1), FUN = function(z) {
            paste0("I(", z, "^2)")
            
          })
          paste(pr, collapse = " + ")
        } else {
          paste0("I(", y, "^2)")
        }
      })
      
      rs <- paste(ls, qs, ps, sep = " + ")
      rs <- gsub(" + erase_this_oneSB", "", rs, fixed = T)
    }
    
    return(rs)
  })
  
  predictors <- unlist(predictors)
  predictors <- predictors[!duplicated(predictors)]
  
  
  return(paste(response_variable, "~", predictors))
}



# ENM glm
glm_enm <- function(formula, all_data, background, family = binomial(link = "logit"),
                    weights_1_0 = c(1, 10000), effect_evaluation = TRUE, 
                    return_prediction = TRUE, rescale_prediction = TRUE, 
                    return_data = FALSE) { 
  
  # data for fitting models
  ## excluding data from background (just in case)
  exa <- which(paste(background[, 1], background[, 2]) %in% 
                 paste(all_data[, 1], all_data[, 2]))

  ## response variable
  zerosa <- rep(0, nrow(background[-exa, ]))
  onesa <- rep(1, nrow(all_data))

  response_variable <- "Response"
  
  # weights
  weightsa <- c(rep(weights_1_0[1], nrow(all_data)), 
                rep(weights_1_0[2], nrow(background[-exa, ])))

  ## all data
  databa <- cbind(Response = c(onesa, zerosa), rbind(all_data, background[-exa, ]))

  # model 
  ## fitting models
  mod <- glm(formula = as.formula(formula), family = family, data = databa, 
             weights = weightsa)
  
  if (return_prediction == TRUE) {
    ## predicting probabilities
    allres <- predict.glm(mod, newdata = background[, -(1:2)], type = "response")
    
    if (rescale_prediction == TRUE) {
      ## maximum value to re-scale 
      allres <- allres / max(allres)
    }
    
    allres <- cbind(background[, 1:2], prediction = allres)
      
  } else {
    allres <- NULL
  }
  
  # effect evaluation
  if (effect_evaluation == TRUE) {
    eff_eval <- anova(mod, test = "Chisq")
  } else {
    eff_eval <- NULL
  }
  
  # data to return
  if (return_data == TRUE) {
    data_b <- list(all_data = all_data, background = background)
  } else {
    data_b <- NULL
  }
  
  # returning results
  return(list(model = mod, prediction = allres, effects = eff_eval, data = data_b))
}


# predict ENM glm
predict_glm_enm <- function(glm_enm, new_data, rescale_prediction = TRUE) {
  ## preparing data
  cl_nd <- class(new_data)[1]
  allowed_f <- c("RasterLayer", "RasterStack", "RasterBrick", 
                 "data.frame", "matrix")
  
  if (!cl_nd %in% allowed_f) {
    stop("Class of new data is not valid")
  } else {
    if (cl_nd %in% allowed_f[1:3]) {
      baser <- new_data[[1]]
      names(baser) <- "prediction"
      new_data <- as.data.frame(na.omit(new_data[]))
    } 
  }
  
  ## predicting probabilities
  allres <- predict.glm(glm_enm$model, newdata = new_data, type = "response")
  
  if (rescale_prediction == TRUE) {
    ## maximum value to re-scale 
    allres <- allres / max(allres)
  }
  
  ## returning results
  if (cl_nd %in% allowed_f[1:3]) {
    baser[!is.na(baser[])] <- allres
    return(baser)
  } else {
    return(allres)
  }
}
