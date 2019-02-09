
# Estimate the BLIN model, specifically for sparse data settings 
fit_blin <- function(maindir, resultsfilename="results", ncores=4, verbose=FALSE, write=TRUE, delete=TRUE)
{
  
  cvtype <- "random"    
  seed <- 1   # random seed (i.e. for CV partitioning)
 
  
  writedir <- outdir <- file.path(maindir, resultsfilename)
  if(!dir.exists(writedir)){ dir.create(writedir) }
  datadir <- file.path(maindir, "data")
  cowsfile <- file.path(datadir, "COW_country_codes.csv")
  load(file.path(datadir, "data.RData"))
  
  
  #### Fits function
  results <- cv_pr_envtreat_final(Y, D, X, ncores, cowsfile, keep, remove, withnet=TRUE,
                                  seed=seed, write=write, outdir=outdir, verbose=verbose, delete=delete)
  if(!dir.exists(writedir)){dir.create(writedir)}
  write.table(results$A, file.path(writedir, "A.txt"), row.names=T, col.names = T, quote=F)
  write.table(results$B, file.path(writedir, "B.txt"), row.names=T, col.names = T, quote=F)
  write.table(results$beta, file.path(writedir, "beta.txt"), row.names=T, col.names = T, quote=F)
  write.table(results$lambda_min, file.path(writedir, "lambda_min.txt"), row.names=F, col.names = F, quote=F)
  compute_pvals(maindir, resultsdir=writedir, writedir, net=TRUE, delete=FALSE)
  
  
  results_nonet <- cv_pr_envtreat_final(Y, D, X, ncores, cowsfile, keep=NA, remove=NA, withnet=FALSE,
                                        seed=seed, write=write, outdir=outdir, verbose=verbose, delete=delete)  
  if(!dir.exists(writedir)){dir.create(writedir)}
  write.table(results_nonet$beta, file.path(writedir, "beta_nonet.txt"), row.names=T, col.names = T, quote=F)
  write.table(results_nonet$lambda_min, file.path(writedir, "lambda_min_nonet.txt"), row.names=F, col.names = F, quote=F)
  compute_pvals(maindir, resultsdir=writedir, writedir, net=FALSE, delete=TRUE)
  
}



# Cross-validation function wrapper for final env.treat results
cv_pr_envtreat_final <- function(Y, D, X, ncores, cowsfile, keep, remove, withnet=TRUE,
                              seed=1, verbose=FALSE, maxit=1e5, thresh=1e-7, ncv=10,  
                              write=FALSE, outdir=NA, delete=FALSE )
{
  cows <- read.csv(cowsfile, header=T, stringsAsFactors = F)   # read country names
  writedir <- outdir
  if(!is.na(writedir)){
    if(!dir.exists(as.character(writedir)) & write){ dir.create( as.character(writedir)) }
  }
  cvtype <- "random"    # cross-validation type, random vs. years
  penalty <- 1   # penalty for glmnet
  use_cov <- TRUE
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha") }
  
  
  S <- max(c(max(D$i), max(Y$i), max(X$i)))
  L <- max(c(max(D$j), max(Y$j), max(X$j)))
  tmax <- max(c(max(D$t), max(Y$t), max(X$t)))
  p <- ncol(X) - 3 - 1
  
  countries <- sort(unique(Y$cowcode))
  treaties <- sort(unique(Y$treaty))
  
  # Sort Y
  rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
  indices <- order(rows)
  Yreg <- Y$ratification_year[indices]  # 1's and zeros, in orde
  
  
  
  pf <- rep(1, length(keep))   # penalty factor, can manipulate to avoid penalizing betas if desired (omitted in this function version)
  
  if(withnet){
    Xreg <- build_design_additive(D, X, type="biten", use_cov=TRUE, sparsedata = TRUE, S=S, L=L, tmax=tmax)  # build design matrix
    rows <- rows[indices]  # sorted rows to keep
    Xreg <- Xreg[rows, ]  # subset
    Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
    Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
    
    # if(write){
      save(Xreg, file=file.path(outdir, "Xreg.RData"))
      save(Yreg, file=file.path(outdir, "Yreg.RData"))
    # }
    
  } else {
    Xreg <- X[,!(names(X) %in% c("i","j","t","intercept"))]
    Xreg[is.na(Xreg)] <- 0
    rows <- X$i + (X$j-1)*S + (X$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Xreg <- as.matrix(Xreg[indices, ])  # 1's and zeros, in order
    
    keep <- 1:ncol(Xreg)
    p <- ncol(Xreg)
    remove <- NA
    
  }
  
  
  results <- cv_pr_lasso(Yreg, Xreg[,keep], outdir, penalty, seed=seed, ncv=ncv, verbose=verbose, 
                         maxit=maxit, ncores=ncores, thresh=thresh, cvtype=cvtype, years=NA, 
                         write=write, delete=delete)  

  if(verbose){
    cat("Done with CV")
  }
  
  #### Form coefficient results
  allcoefs <- rep(NA, ncol(Xreg))
  allcoefs[keep] <- results$coefs[-1]
  allcoefs <- c(results$coefs[1], allcoefs)  # add in intercept
  
  if(withnet){
    A <- t(matrix(allcoefs[1:(S^2) + 1], S,S))   # transposed A for interpretability
    colnames(A) <- rownames(A) <- countries
    B <- matrix(allcoefs[1:(L^2) + 1 + S^2], L,L)   # A as in coding, not t(A)
    colnames(B) <- rownames(B) <- treaties
  } else { A <- B <- NA}
  
  if(use_cov){
    beta <- matrix(c(allcoefs[1], tail(allcoefs, p)), ncol=1)
    rownames(beta) <- c("intercept", names(X)[!(names(X) %in% c("i","j","t","intercept"))]) 
  } else {beta <- NA}
  
  Yhat <- results$Yhat
  fit <- results$fit
  cvms <- results$cvms
  cvm_full <- results$cvm_full
  lambda_min <- results$lambda_min
  imin <- results$imin
  nl_full <- results$nl_full
  nl_cv <- results$nl_cv
  prfull <- cvm_full[imin]
  ####
  
  
  if(write){
    if(withnet){
      write.table(A, file.path(writedir, "A.txt"), row.names=T, col.names = T, quote=F)   # NOTE TRANSPOSE
      write.table(B, file.path(writedir, "B.txt"), row.names=T, col.names = T, quote=F)
      write.table(beta, file.path(writedir, "beta.txt"), row.names=T, col.names = T, quote=F)
      
    } 
    write.table(beta, file.path(writedir, "beta_nonet.txt"), row.names=T, col.names = T, quote=F)
    
  }
  
  
  return(list(A=A, B=B, beta=beta, Yhat=results$Yhat, fit=results$fit, cvms=results$cvms, 
              cvm_full=results$cvm_full, lambda_min=results$lambda_min, imin=results$imin,
              nl_full=results$nl_full, nl_cv=results$nl_cv, prfull=prfull))
}
  
  
  

