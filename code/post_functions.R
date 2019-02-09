# Functions to process results of fitting



# Compute pvalues and other summary statistics using selectiveInference package
#    maindir = master directory  
#    resultsdir = directory where results are stored (prefit2.RData, unprocessed_results.RData)
#    writedir = directory to write outputs
#    net = Logical of whether network effects included (defaul TRUE, effects included)
#  Writes out results files of coefficients and their p-values   
compute_pvals <- function(maindir, resultsdir, writedir, net=TRUE, delete=FALSE)
{
  require("survival")
  require("intervals")
  nonet <- !as.logical(net)
  if(nonet){
    ps="_nonet"
  } else {
    ps=""
  }
  
  
  # resultsdir <- file.path(maindir, resultsdir)
  outdir <- writedir   # <- file.path(maindir, writedir)
  datadir <- file.path(maindir, "data")
  
  selectiveInfLib <- file.path(maindir, "code/selectiveInferenceFM")
  
  
  # source(file.path(maindir, "code/MLE_functions.R"))
  # library("RColorBrewer")
  
  
  if(!dir.exists(writedir)){dir.create(writedir)}
  
  
  load(file.path(datadir, "data.RData"))
  load(file.path(resultsdir, "Xreg.RData"))
  load(file.path(resultsdir, "Yreg.RData"))
  if(!nonet){
    lambda_min <- c(as.matrix(read.table(file.path(resultsdir, "lambda_min.txt"), header=FALSE)))
  } else {
    lambda_min <- c(as.matrix(read.table(file.path(resultsdir, "lambda_min_nonet.txt"), header=FALSE)))
  }
  S <- max(c(max(D$i), max(Y$i), max(X$i)))
  L <- max(c(max(D$j), max(Y$j), max(X$j)))
  tmax <- max(c(max(D$t), max(Y$t), max(X$t)))
  p <- ncol(X) - 3 - 1
  
  countries <- sort(unique(Y$cowcode))
  treaties <- sort(unique(Y$treaty))
  
  
  
  if(delete){
    # file.remove(file.path(resultsdir, "Xreg.RData"))
    # file.remove(file.path(resultsdir, "Yreg.RData"))
    f <- list.files(path=resultsdir,  pattern="lambda_min")
    if(length(f) > 0){
      for(i in 1:length(f)){
        file.remove(file.path(resultsdir, f[i]))  
      }
    }
  }
  
  # nonet <- FALSE
  
  if(!nonet){
    A <- as.matrix( read.table(file.path(resultsdir, "A.txt"), header=TRUE) )
    B <- as.matrix( read.table(file.path(resultsdir, "B.txt"), header=TRUE) )
    beta <- as.matrix( read.table(file.path(resultsdir, "beta.txt"), header=TRUE) )
    
    S <- ncol(A)
    L <- ncol(B)

  } else {
    beta <- as.matrix( read.table(file.path(resultsdir, "beta_nonet.txt"), header=TRUE) )

  }

  px <- length(c(as.matrix(beta))) - 1
  
  
  
  #### Standard errors using selectiveInference package
  install.packages(file.path(selectiveInfLib, "selectiveInference_1.2.2.tar.gz"), repos = NULL, type="source", lib = selectiveInfLib)
  library("selectiveInference", lib.loc = selectiveInfLib)  # after Frank's edits
  # library("selectiveInference", lib.loc = )  # after Frank's edits
  # Xreg <- results$Xreg
  # Lcheck <- lambda_min*nrow(Xreg) < max(fit$lambda) & lambda_min*nrow(Xreg) > min(fit$lambda)
  
  penbeta <- TRUE
  
  if(penbeta & !nonet){
    # theta <- coef(results$fit)[,imin]   # with intercept, but fixedLassoInf knows to remove first entry (in fact intercept required)
    intercept <- as.matrix(beta)[1]
    theta <- c(intercept, c(t(as.matrix(A))), c(as.matrix(B)), c(as.matrix(beta)[-1]))
    
    remove0 <- which(theta[-1] == 0)
    keep0 <- setdiff(1:(length(theta) - 1), remove0)
    keep1 <- intersect(keep, keep0)
    remove1 <- setdiff(1:(length(theta) - 1), keep1)
    
    
  
  } else { 
    theta <- c(as.matrix(beta))
    
    Xreg <- Xreg[, tail(1:ncol(Xreg), px)]
    remove1 <- which(theta[-1] == 0)
    keep1 <- setdiff(1:(length(theta) - 1), remove1)
  }
  
  Xreg1 <- Xreg[,keep1]
  theta1 <- c(theta[1], theta[-1][keep1])
  preds <- length(theta1) > 1
  
  
  
  if(preds){
    out <- fixedLassoInf(x=(Xreg1), y=Yreg, beta=theta1, lambda=lambda_min*nrow(Xreg1), alpha=.05, family = "binomial", tol.beta=1e-10)  
    
    # cat("Fraction of nonzero coefficients that are significant by p-value:  ", mean(out$pv <= .05), "\n")
    
    
    if(sum(range(Xreg[,keep1] - Xreg1)) != 0){stop("Something went wrong with indexing in terms of original indices")}
    
    se <- rep(NA, ncol(Xreg))
    se[keep1] <- out$sd
    
    if(!nonet){
      sea <- c(t(matrix(se[1:(S^2)], S, S)))   # transpose A results
      seb <- se[S^2 + 1:(L^2)]
      sebeta <- c(0, tail(se, px))
      
      # save(out,  remove1, keep1, se, sea, seb, sebeta, theta, file=file.path(outdir, "se_from_selectiveInference.RData"))
      
      r <- range(c(
        range(which(is.na(sea)) - which(is.na(c((A))) | c((A)) == 0)),
        range(which(is.na(seb)) - which(is.na(c(B)) | c(B) == 0)),
        range(which(is.na(sebeta)) - which(is.na(c(beta)) | c(beta) == 0))))   # check if NAs and zeros line up... they do!
      # cat("Match 0s and NA check, should be zero: ", r,"\n")
    } else {
      sebeta <- c(0, tail(se, px))
      # save(out,  remove1, keep1, se, sebeta, theta, file=file.path(outdir, "se_from_selectiveInference.RData"))
    }
  }
  
  # print(warnings())
  # assign("last.warning", NULL, envir = baseenv())  # clear warnings
  ####
  
  
  
  
  #### Shape into useful output
  if(preds){
    beta_names <- rownames(beta)
    alpha <- .05
    # px <- length(beta) - 1
    
    if(!nonet){
      ciA <- data.frame(matrix(NA, S^2, 13))   ;  names(ciA) <- c("i", "j", "cow1", "cow2", "estimate", "lower", "upper", "pval", "signif", "se", "alpha", "nonid_keep", "nonid_remove")
      ciB <- data.frame(matrix(NA, L^2, 13))   ;  names(ciB) <- c("i", "j", "treaty1", "treaty2", "estimate", "lower", "upper", "pval", "signif", "se", "alpha", "nonid_keep", "nonid_remove")
      cibeta <- data.frame(matrix(NA, px + 1, 8))   ;  names(cibeta) <- c("covariate", "estimate", "lower", "upper", "pval", "signif", "se", "alpha")
      
      ciA[,1:2] <- cbind(rep(1:(S), times=S), rep(1:(S), each=S))   # columnwise unfolding
      ciB[,1:2] <- cbind(rep(1:(L), times=L), rep(1:(L), each=L))
      cibeta$covariate <- beta_names
      ciA$cow1 <- countries[ciA$i]  ;  ciA$cow2 <- countries[ciA$j]
      ciB$treaty1 <- treaties[ciB$i]  ;  ciB$treaty2 <- treaties[ciB$j]
      ciA$alpha <- ciB$alpha <- cibeta$alpha <- alpha
      
      ciA$estimate <- c(A)
      ciB$estimate <- c(B)
      cibeta$estimate <- c(beta)
      
      ciA$se <- sea
      ciB$se <- seb
      cibeta$se <- sebeta 
      
      ci <- matrix(NA, ncol(Xreg), 2)
      ci[keep1,] <- out$ci
      ciA[,6:7] <- ci[c(t(matrix(1:(S^2), S, S))),]   # transposed!
      ciB[,6:7] <- ci[S^2 + 1:(L^2),]
      cibeta[,3:4] <- rbind(0, ci[tail(1:ncol(Xreg), px),])
      
      pv <- rep(NA, ncol(Xreg))
      pv[keep1] <- out$pv
      ciA$pval <- pv[c(t(matrix(1:(S^2), S, S)))]   # transposed!
      ciB$pval <- pv[S^2 + 1:(L^2)]
      cibeta$pval <- c(0, tail(pv, px))
      
      cisignif <-  1*!(ci[,1] < 0 & ci[,2] > 0 )   # significance based on confidence intervals
      # cat("Fraction of pvalue < alpha and non-0-containing confidence intervals agree is : ", mean(cisignif == 1*(pv < alpha), na.rm=T), "\n")
      cisignif2 <- 1*(pv < alpha)  # significance based on pvalue
      cisignif <- cisignif*cisignif2   # only significant if both significant
      
      ciA$signif <- cisignif[c(t(matrix(1:(S^2), S, S)))]   
      ciB$signif <- cisignif[S^2 + 1:(L^2)]
      cibeta$signif <- c(1, tail(cisignif, px))
      
      ciA$nonid_keep <- ciA$nonid_remove <- ciB$nonid_keep <- ciB$nonid_remove <- 0
      
      # if(length(nonid_keep) > 0 & !is.na(nonid_keep[1])){
      #   ia <- t(matrix(1:(S^2), S, S))
      #   ib <- matrix(1:(L^2) + S^2, L, L)
      #   nonid_remove2 <- unlist(nonid_remove)
      #   
      #   ka <- sapply(nonid_keep[nonid_keep <= S^2], function(z) which(z == ia))  # transposed (interpretable) A index
      #   if(length(ka) > 0){
      #     ciA$nonid_keep[ka] <- 1
      #     ra <- sapply(nonid_remove2[nonid_remove2 <= S^2], function(z) which(z == ia))  
      #     ciA$nonid_remove[ra] <- 1
      #   }
      #   
      #   kb <- sapply(nonid_keep[nonid_keep > S^2], function(z) which(z == ib))  
      #   if(length(kb) > 0){
      #     ciB$nonid_keep[kb] <- 1
      #     rb <- sapply(nonid_remove2[nonid_remove2 > S^2], function(z) which(z == ib))  
      #     ciB$nonid_remove[rb] <- 1
      #   }
      # } 
      
      write.table(ciA[, c("i", "j", "cow1", "cow2", "estimate", "pval")], file.path(writedir, paste0("ciA", ps, ".txt")), row.names=FALSE)
      write.table(ciB[, c("i", "j", "treaty1", "treaty2", "estimate", "pval")], file.path(writedir, paste0("ciB", ps, ".txt")), row.names=FALSE)
      write.table(cibeta[, c("covariate", "estimate", "pval")], file.path(writedir, paste0("cibeta", ps, ".txt")), row.names=FALSE)
      
    } else {
      cibeta <- data.frame(matrix(NA, px + 1, 8))   
      names(cibeta) <- c("covariate", "estimate", "lower", "upper", "pval", "signif", "se", "alpha")
      
      cibeta$covariate <- beta_names
      cibeta$alpha <- alpha
      cibeta$estimate <- c(beta)
      cibeta$se <- sebeta 
      
      ci <- matrix(NA, ncol(Xreg), 2)
      ci[keep1,] <- out$ci
      cibeta[,3:4] <- rbind(0, ci[tail(1:ncol(Xreg), px),])
      
      pv <- rep(NA, ncol(Xreg))
      pv[keep1] <- out$pv
      cibeta$pval <- c(0, tail(pv, px))
      
      cisignif <-  1*!(ci[,1] < 0 & ci[,2] > 0 )   # significance based on confidence intervals
      # cat("Fraction of pvalue < alpha and non-0-containing confidence intervals agree is : ", mean(cisignif == 1*(pv < alpha), na.rm=T), "\n")
      cisignif2 <- 1*(pv < alpha)  # significance based on pvalue
      cisignif <- cisignif*cisignif2   # only significant if both significant
      cibeta$signif <- c(1, tail(cisignif, px))
      
      write.table(cibeta[, c("covariate", "estimate", "pval")], file.path(writedir, paste0("cibeta", ps, ".txt")), row.names=FALSE)
    }
  }
  ####
  
}




# Run counterfactual study (includes call to "process_counterfactual(.) below)
#    maindir = master directory  
#    resultsdir = directory where results are stored (prefit2.RData, unprocessed_results.RData)
#    writedir = directory to write outputs
#    nsims = total number of simulations
#    ncores = number of parallel cores to use
#    remove = logical indicating whether output RData files should be removed (default TRUE, remove them)
#    cswap = country whose ratification was swapped; only included to save with results
#    tswap = treaty whose ratification was swapped; only included to save with results
#    year0 = Year of ratification that was swapped; only included to save with results
#   
#  Writes out "res_counterfactual.RData" containing desired posterior means
#
run_counterfactual <- function(maindir, resultsdir="results", writedir="results", nsims=1e2, ncores=1, cswap=2, tswap=40793, year0=1992, verbose=FALSE, delete=TRUE)
{
  resultsdir <- file.path(maindir, resultsdir)
  outdir <- writedir <- file.path(maindir, writedir)
  datadir <- file.path(maindir, "data")
  
  if(as.numeric(ncores == 1)){
    if(verbose){cat("Starting", 1, "process of n=", nsims, "simulations \n")}

    run_counterfactual_sub(maindir, resultsdir, writedir,
                           nsims=nsims, seed0=1,
                           cswap=cswap, tswap=tswap, year0=year0)
  } else {

    neach <- floor(as.numeric(nsims) / as.numeric(ncores))
    nprocess <- floor(as.numeric(nsims) / as.numeric(neach))
    if(verbose){cat("Starting", nprocess, "processes of n=", neach, "simulations \n")}

    registerDoMC(cores=nprocess)
    mcoptions <- list(preschedule=FALSE, set.seed=TRUE)
    fitlist <- foreach(i=1:nprocess, .options.multicore=mcoptions) %dopar% {   #, .packages=c("glmnet")
      seed0 <- (i- 1)*neach + 1   # seed based on task and number of sims

      run_counterfactual_sub(maindir, resultsdir, writedir,
                             nsims=neach, seed0=seed0,
                             cswap=cswap, tswap=tswap, year0=year0)
    }

    
  }

  # process all the results!
  process_counterfactual(maindir, resultsdir, writedir=resultsdir, remove=delete, cswap=cswap, tswap=tswap, year0=year0)
  
  if(delete){
    f <- list.files(path=resultsdir,  pattern="Yreg")
    if(length(f) > 0){
      for(i in 1:length(f)){
        file.remove(file.path(resultsdir, f[i]))  
      }
    }
    
    f <- list.files(path=resultsdir,  pattern="Xreg")
    if(length(f) > 0){
      for(i in 1:length(f)){
        file.remove(file.path(resultsdir, f[i]))  
      }
    }
    
    f <- list.files(path=resultsdir,  pattern="Yhat")
    if(length(f) > 0){
      for(i in 1:length(f)){
        file.remove(file.path(resultsdir, f[i]))  
      }
    }
  }
  
  
}



# Sub-function for counterfactual that does the heavy lifting
run_counterfactual_sub <- function(maindir, resultsdir="results", writedir="results", nsims=1e2, seed0=NULL, cswap=2, tswap=40793, year0=1992)
{
  lag = 3
  n = nsims
  
  datadir <- file.path(maindir, "data")
  if(!dir.exists(writedir)){dir.create(writedir)}
  
  region_income <- read.table(file.path(datadir, "region_income.txt"), header=TRUE)
  
  
  #### Load data
  setwd(resultsdir)
  A <- read.table("A.txt")
  B <- read.table("B.txt")
  beta <- read.table("beta.txt")

  load(file.path(datadir, "data.RData"))   # for Y,D,X
  load(file.path(resultsdir, "Yhat.RData"))   # for Yhat
  ####
  
  for(swap in c(TRUE, FALSE)){
  
    if(swap){
      ps <- "_swapped"
      t0 <- unique(Y$t[Y$year==year0])+1
      NApairs = matrix(c(unique(Y$i[Y$cowcode == cswap]), unique(Y$j[Y$treaty == tswap])), nrow=1)
      Din <- D  ;  Yin <- Y
      orig_rat <- NULL
      for(k in 1:nrow(NApairs)){
        orig_rat <- c(orig_rat, Yin$ratification_year[Yin$i == NApairs[k,1] & Yin$j == NApairs[k,2] & Yin$t == (t0-1)])
        Din$ratification_year[Din$i == NApairs[k,1] & Din$j == NApairs[k,2] & Din$t == (t0)] <- 0
        Yin$ratification_year[Yin$i == NApairs[k,1] & Yin$j == NApairs[k,2] & Yin$t == (t0-1)] <- 0
      }
    } else {
      ps <- "_unswapped"
      t0 <- unique(Y$t[Y$year==year0])+1
      Din <- D  ;  Yin <- Y
      orig_rat <- NULL
      NApairs <- NULL
    }
    
    outfile <- paste0("predict_lag", lag, "_", year0+1,"_2000_seed", seed0, "_nsims", n,  ps, ".RData")  # countries=NULL, treaties=NULL, 
    
    rswap <- roll_forward_predict(t0, n, A, B, beta, Yin, Din, X, lag, region_income, seed0=seed0, 
                                  verbose=FALSE, NApairs=NApairs, filename=outfile, 
                                  datadir=datadir, outdir=writedir)  # countries=NULL, treaties=NULL, 
    
  }
  
  # process_counterfactual(maindir, resultsdir, writedir)
  
}




# SUBFUNCTION: Process counterfactual study
#    maindir = master directory  
#    resultsdir = directory where results are stored 
#    writedir = directory to write outputs
#    remove = logical indicating whether output RData files should be removed (default TRUE, remove them)
#    cswap = country whose ratification was swapped; only included to save with results
#    tswap = treaty whose ratification was swapped; only included to save with results
#    year0 = Year of ratification that was swapped; only included to save with results
#    
#  Writes out "res_counterfactual.RData" containing required posterior means
#
process_counterfactual <- function(maindir, resultsdir, writedir, remove=TRUE, cswap=2, tswap=40793, year0=1992)
{
  
  delete <- as.logical(remove)
  # resultsdir <- file.path(maindir, resultsdir)
  outdir <- writedir  # <- file.path(maindir, writedir)
  # selectiveInfLib <- file.path(maindir, "code/selectiveInferenceFM")
  
  
  # source(file.path(maindir, "code/MLE_functions.R"))
  # library("RColorBrewer")
  
  
  if(!dir.exists(writedir)){dir.create(writedir)}
  
  
  #### Swapped
  # direc <- "swapped_results"
  fs <- list.files(path=resultsdir, pattern="_swapped.RData")
  
  if(length(fs) > 0){
    load(file.path(resultsdir, fs[1]))
    Ysum <- psum <- Ytot <- ptot <- matrix(0, S, L)
    
    for(i in 1:length(fs)){
      load(file.path(resultsdir, fs[i]))
      Ysum <- Ysum + apply(Yhat, 1:2, sum, na.rm=T)
      psum <- psum + apply(phat, 1:2, sum, na.rm=T)
      Ytot <- Ytot + apply(Yhat, 1:2, function(z) sum(!is.na(z)))
      ptot <- ptot + apply(phat, 1:2, function(z) sum(!is.na(z)))
      # rm(Yhat, phat)
    }
    denom <- max(Ytot, na.rm=T) / dim(Yhat)[3]
    Ymean_swap <- Ysum / denom
    pmean_swap <- psum / denom
    # is.infinite(Ymean_swap) <- is.infinite(pmean_swap) <- NA
    # save(Ysum, psum, Ytot, ptot, Ymean_swap, pmean_swap, cswap, tswap, year0, denom, file=file.path(writedir, "all_swapped.RData"))
    if(remove){
      for(i in 1:length(fs)){
        file.remove(file.path(resultsdir, fs[i]))
      }
    }
  } else {
    warning("No swapped results found")
    Ymean_swap <- pmean_swap <- denom <-  NA
    cswap <- tswap <- year0 <- NA
  }
  ####
  
  
  
  #### Unswapped
  # direc <- "unswapped_results"
  fs <- list.files(path=resultsdir, pattern="_unswapped.RData")
  
  if(length(fs) > 0){
    load(file.path(resultsdir, fs[1]))
    Ysum <- psum <- Ytot <- ptot <- matrix(0, S, L)
    
    for(i in 1:length(fs)){
      load(file.path(resultsdir, fs[i]))
      Ysum <- Ysum + apply(Yhat, 1:2, sum, na.rm=T)
      psum <- psum + apply(phat, 1:2, sum, na.rm=T)
      Ytot <- Ytot + apply(Yhat, 1:2, function(z) sum(!is.na(z)))
      ptot <- ptot + apply(phat, 1:2, function(z) sum(!is.na(z)))
      # rm(Yhat, phat)
    }
    Ymean0 <- Ysum / denom
    pmean0 <- psum / denom
    # is.infinite(Ymean0) <- is.infinite(pmean0) <- NA
    # save(Ysum, psum, Ytot, ptot, Ymean0, pmean0, cswap, tswap, year0, denom, file=file.path(writedir, "all_unswapped.RData"))
    if(remove){
      for(i in 1:length(fs)){
        file.remove(file.path(resultsdir, fs[i]))
      }
    }
  }else {
    warning("No unswapped results found")
    Ymean0 <- pmean0 <- denom <-  NA
  }
  ####
  
  
  save(Ymean_swap, pmean_swap, Ymean0, pmean0, cswap, tswap, year0, denom, file=file.path(writedir, "res_counterfactual.RData"))
  
}


