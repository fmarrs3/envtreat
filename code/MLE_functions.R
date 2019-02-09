# BiTEN Project: MLE estimation functions
# Frank Marrs
# 09/13/16
#
# This file contains functions supporting MLE estimation of A and B matrices on trade data based on inputs
# 





##########################
### Fitting functions  ###
##########################




# P/R cross-validation with automatic lambda selection, standardized X
cv_pr_lasso <- function(Yreg, Xreg, outdir, penalty=1, seed=NA, ncv=10, verbose=F, maxit=1e5, ncores=1, 
                        thresh=1e-7, pf=rep(1, ncol(Xreg)), cvtype="random", years=NA, write=FALSE, delete=FALSE)  
{
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha for glmnet(.))") }
  
  if(write){
    dir.create(outdir, showWarnings = F)
  }
  
  
  #### Perform full fit and extract coefficients for lambda
  if(is.numeric(seed)){ set.seed(seed)}  # set seed
  fitfull <- fit <- glmnet(Xreg, y=as.factor(c(Yreg)), alpha=penalty, family='binomial', intercept=T, lambda.min.ratio=1e-12, standardize=T, maxit=maxit, nlambda=100, thresh=thresh, penalty.factor=pf)
  if(write){
    save(fit, file=file.path(outdir, paste0("fullfit.RData")))
  }
  workedfull <- length(fitfull$lambda) 
  if(verbose){
    cat("full fit done, nlambda", workedfull, "\n")
  }
  lambdas <- fitfull$lambda   # save selected lambda sequence
  ####
  
  
  #### Cross-validate
  lambdas <- sort(lambdas, decreasing = T)  # decreasing sequence
  
  if(cvtype == "random"){
    ones_check <- rep(0, ncv)
    if(is.numeric(seed)){ set.seed(seed)}  # set seed
    count <- 0
    while(any(ones_check == 0) & count < 1e4){   # if there are any partitions without 1's, and make loop finite
      count <- count + 1
      # cat("cv partition count", count, "\n")
      rowscramble <- sample(1:nrow(Xreg), size = nrow(Xreg), replace = F)
      cvs <- rowscramble %% ncv + 1
      ones_check <- sapply(1:ncv, function(z) sum(Yreg[cvs == z]))   # number of 1's in Yreg for each cv partition
    }
  } else if (strtrim(cvtype,4) == "year"){
    ncv <- length(unique(years))
    cvs <- match(years, sort(unique(years)))
    ones_check <- rowscramble <- NA
  } else {
    stop("invalid cross-validation type input cvtype")
  }
    
  cvms <- matrix(NA, length(lambdas), ncv)
  assign("last.warning", NULL, envir = baseenv())   # clear warnings
  nlambda <- rep(NA, ncv)  # number of lambda fit
  
  ####
  if(write){
    save(thresh,seed,ncv,maxit,ncores,penalty,pf,cvs,rowscramble,ones_check, file=file.path(outdir, paste0("cv_pr_lasso_inputs.RData")))
  }
  ####
  
  
  if(ncores==1){
    fitlist <- vector("list", ncv)
    for(i in 1:ncv){
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      set.seed(i*1984)
      fitlist[[i]] <- fit <- glmnet(Xreg[train,], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh, penalty.factor = pf)
      worked <- nlambda[i] <- length(fit$lambda)
      
      if(write){
        save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      }
      
      if(verbose){
        cat("\n************************************************\n")
        cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
        print(warnings())
        cat("\n************************************************")
        cat("\n")
        assign("last.warning", NULL, envir = baseenv())
      }
    }
    
  } else if (ncores > 1 ){
    if(write){
      writeLines(c(""), file.path(outdir, "log.txt"))
    }
    
    registerDoMC(cores=round(ncores))
    mcoptions <- list(preschedule=FALSE, set.seed=T)
    fitlist <- foreach(i=1:ncv, .options.multicore=mcoptions, .packages=c("glmnet") ) %dopar% {  # .combine =cbind
      set.seed(i*1984)
      
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      
      fit <- glmnet(Xreg[train,], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh, penalty.factor = pf)
      
      if(write){
        save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      }
      
      worked <- length(fit$lambda)
      
      if(verbose){
        sink(file.path(outdir, "log.txt"), append=TRUE)   # write out to log file
        cat("\n************************************************\n")
        cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
        print(warnings())
        cat("\n************************************************")
        cat("\n")
        fit
      }
    }
    
  } else {stop("ncores must be numeric >= 1")}
  
  for(i in 1:ncv){
    test <- which(cvs == i)
    train <- (1:nrow(Xreg))[-test]
    
    if(write & ncores > 1){
      load(file=file.path(outdir, paste0("pr_cv", i, ".RData")))
    } else {
      fit <- fitlist[[i]]
    }
    
    #glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit)
    worked <- nlambda[i] <- length(fit$lambda)
    
    if(worked > 1){   # if fit worked
      Yhats <- predict(fit, newx=Xreg[test,], type="response")
      
      pr_temp <- rep(NA, ncol(Yhats))
      for(j in 1:ncol(Yhats)){
        # pr_temp[j] <- pr_curve(Yhats[,j], Y[test], n=2000)$auc
        # pr_temp[j] <- roc_curve(Yhats[,j], Y[test], n=min(length(Y[test]), 2500))$auc
        # pr_temp[j] <- simple_roc_auc(Yhats[,j], Y[test])
        pr_temp[j] <- simple_pr_auc(Yhats[, j], Yreg[test])
        # pr_temp[j] <- auc(roc(Yhats[,j], Y[test]))
      }
      
      cvms[match(fit$lambda, lambdas), i] <- pr_temp   # save lambda values
      # save(fit, pr_temp, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
    }
    if(verbose){
      cat("processed cv", i, "of", ncv, ", nlambda", worked, "\n")
    }
  }
  
  # Calculate minimum lambda value
  mean_cvms <- apply(cvms, 1, mean, na.rm=T)
  imin <- which(mean_cvms == max(mean_cvms, na.rm=T))[1]
  lambda_min <- lambdas[imin]
  ####
  
  
  #### Calculate results from full fit
  fit <- fitfull
  if(workedfull > 1){
    Yhats <- predict(fitfull, newx=Xreg, type="response")
    cvm_full <- rep(NA, ncol(Yhats))
    for(j in 1:ncol(Yhats)){
      cvm_full[j] <- simple_pr_auc(Yhats[,j], Yreg) # pr_curve(Yhats[,j], Y, n=min(length(Y), 2500))$auc
    }
    Yhat <- predict(fit, newx=Xreg, type="response", s=lambda_min)
    coefs <- coef(fit, s=lambda_min)
    
  } else {
    fit <- Yhat <- coefs <- cvm_full <- NA
  }
  ####
  
  if(write & delete){
    f1 <- list.files(outdir, pattern="pr_cv")
    if(length(f1 > 0)){
      for(i in 1:length(f1)){
        file.remove(file.path(outdir, f1[i]))
      }
    }
    
    f1 <- list.files(outdir, pattern="fullfit")
    if(length(f1 > 0)){
      for(i in 1:length(f1)){
        file.remove(file.path(outdir, f1[i]))
      }
    }
    
    f1 <- list.files(outdir, pattern="log.txt")
    if(length(f1 > 0)){
      for(i in 1:length(f1)){
        file.remove(file.path(outdir, f1[i]))
      }
    }
    
    f1 <- list.files(outdir, pattern="cv_pr_lasso_inputs")
    if(length(f1 > 0)){
      for(i in 1:length(f1)){
        file.remove(file.path(outdir, f1[i]))
      }
    }
  }
  
  #### To remove??
  save(Yhat, file=file.path(outdir, "Yhat.RData"))
  ####
  
  
  return(list(coefs=coefs, Yhat=Yhat, fit=fitfull, Yreg=Yreg, Xreg=Xreg, cvms=cvms, cvm_full=cvm_full, lambda_min=lambda_min, imin=imin, nl_full=workedfull, nl_cv=nlambda))
}



##############################
###  Prediction functions  ###
##############################


# roll forward prediction from some initial time, including random realizations of signatures
#  t0 is first year of PREDICTION
#  stochastic prediction requires number of simulations
roll_forward_predict <- function(t0, nsims, A, B, beta, Y, D, X, lag, 
                                    region_income, tfinal=NULL, response="ratification_year", seed0=1, verbose=F, 
                                    NApairs=NULL, write_interval=ceiling(nsims/10), outdir=getwd(), 
                                    filename=NULL, datadir=NULL)  # countries=NULL, treaties=NULL, 
{
  if(is.null(tfinal)){ tfinal <- max(Y$t)}
  trange <- t0:tfinal
  year_range <- unique(Y$year[Y$t == t0] ) : unique(Y$year[Y$t == tfinal] )
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  beta <- as.matrix(beta)
  
  # Initialize arrays
  Ynew <- Y[Y$t >= t0 & Y$t <= tfinal,]
  Yold <- Y[Y$t < t0,]
  Dnew <- D[D$t >= t0 & D$t <= tfinal,]
  Xnew <- X[X$t >= t0 & X$t <= tfinal,]
  S <- max(Dnew$i)  ;   L <- max(Dnew$j)    
  Yhat <- phat <- array(0, c(S,L,length(trange),nsims))
  countries <- sapply(1:S, function(z) unique(Dnew$cowcode[Dnew$i == z]))   # all countries of interest in future
  treaties <- sapply(1:L, function(z) unique(Dnew$treaty[Dnew$j == z]))    # all treaties of interest in future
  dimnames(Yhat)[[1]] <- dimnames(phat)[[1]] <- countries
  dimnames(Yhat)[[2]] <- dimnames(phat)[[2]] <- treaties
  dimnames(Yhat)[[3]] <- dimnames(phat)[[3]] <- year_range
  region_income$i <- Y$i[match(region_income$cowcode, Y$cowcode)]   # add column for i in region_income indicators 
  region_income <- region_income[order(region_income$i),]    # reorder
  
  # Read in big X if X is not big
  if(nrow(X) <  S*L*length(trange)){
    X <- build_big_X(t0, X, S=S, L=L, tfinal=max(X$t), readfile=T)
  }
  
  
  load(file.path(datadir, "pre-prediction.RData"))
  
  
  # build possible signatures for all future times and simulations
  for(t in trange){
    k <- t - t0 + 1
    iposs <- cbind(possibles[rep(1:nrow(possibles), times=nsims), ], k, rep(1:nsims, each=nrow(possibles)))   # indices in first year that can be signed
    Yhat[,,k,] <- phat[,,k,] <- NA   # all NAs in time period
    Yhat[iposs] <- phat[iposs] <- 0   # possible ratifications

    jlate_entry <- which(t >= late_entry1[,3])   # AFTER late entry
    if(length(jlate_entry) > 0){
      # cat("late entry k=", k, "t=", t, "\n")
      ilate <- cbind(late_entry1[rep(jlate_entry, times=nsims), 1:2], k, rep(1:nsims, each=length(jlate_entry)) )
      Yhat[ilate] <- phat[ilate] <- 0   # 0s for i,j pairs that enter late
    }

    jearly_exit <- which(t > early_exit[,3])   # AFTER exit
    if(length(jearly_exit) > 0){
      # cat("early exit k=", k, "t=", t,"\n")
      iearly <- cbind(early_exit[rep(jearly_exit, times=nsims), 1:2], k, rep(1:nsims, each=length(jearly_exit)) )
      Yhat[iearly] <- phat[iearly] <- NA   # NAs for i,j pairs that leave early
    }
  }
  
  # Save all ratifications for future updates
  NAinit1 <- already_signed
  #unique(Yold[Yold$ratification_year == 1, c("i", "j")])   # save already signed treaty/country pairs
  if(!is.null(NApairs)){
    NAinit1 <- rbind(NAinit1, NApairs)
  }
  
  
  # Prediction
  allcoefs <- c(c(t(as.matrix(A))), c(as.matrix(B)), c(as.matrix(beta))[-1])   # coefs without intercept
  allcoefs[is.na(allcoefs)] <- 0
  keep <- which(allcoefs != 0)
  
  # Make dummy D and X variable arrays to avoid rebuilding
  #   shell is all pairs of i,j such that i<S and j<L for a single time period
  Dshell <- D[0,]
  Dshell[1:(S*L),] <- 0
  Xshell <- X[0,]
  Xshell[1:(S*L),] <- 0
  Dshell$t <- Xshell$t <- 1
  Dshell$i <- Xshell$i <- rep(1:S, times=L)
  Dshell$j <- Xshell$j <- rep(1:L, each=S)
  Dshell$cowcode <- Y$cowcode[match(Dshell$i, Y$i)]
  Dshell$treaty <- Y$cowcode[match(Dshell$j, Y$j)]
  Xshell$intercept <- 1
  Xshell[, -c(1:4)] <- NA   # NAs for all in Xshell
  
  # filename to save
  if(is.null(filename)){
    filename <-  paste0("predict_lag", lag, "_", min(year_range),"_", max(year_range), "_seed", seed0,  ".RData")
  }
  
  
  # Run loop to make predictions, write out results periodicially
  for(i in 1:nsims){
    set.seed(seed0+i-1)   # set seed for repeatability
    NAremove <- NAinit1   # initialize country/treaty pairs to remove from dataset for each simulation
    # Ytemp <- Y[Y$t < t0 & Y$ratification_year==1, c("i", "j", "t", "ratification_year")]   # save all ratifications
    
    for(t in t0:tfinal){    # t is the year of the prediction.
      k <- t - t0 + 1   # index in vavriables to save
      
      if(t == t0){    # initialize possibly lagged autoregressive array and covariate array
        Dpred <- Dnew[Dnew$t == t,] 
        Xpred <- Xnew[Xnew$t == t,]
      }   
      
      # Set impossible signatures to NA (should already be, but double-check)
      Yhat[cbind(NAremove, k, i)] <- phat[cbind(NAremove, k, i)] <- NA
      
      # Build design matrix for appropriate D, X
      Xpred$t <- Dpred$t <- 1
      Xreg <- build_design_additive(Dpred, Xpred, sparsedata=T, write=F, S=S, L=L, tmax=1)
      
      # Subset and pare design matrix to predict
      Xreg <- Xreg[,-(S^2 + L^2 + 1)]   # remove intercept
      keep_mat <- which(!is.na(Yhat[,,k,i]), arr.ind=T)   # i,j pairs to predict, in columnwise order 
      rows <- keep_mat[,1] + (keep_mat[,2]-1)*S   #  + (Y$t-1)*S*L   # unfolded indices
      Xreg <- Xreg[sort(rows), ]   # rows of Xreg that pertain to the entres in Yhat
      Xreg[which(is.na(Xreg))] <- 0   # Set NAs to zero
      
      # Calculate new prediction probabilities and predict
      keep_mat <- keep_mat[order(rows),]   # order entries as calculated
      keep_mat1 <- as.matrix(cbind(keep_mat, k, i))   # including time and sim indices
      # ptilde <- as.matrix(Xreg[, keep] %*% allcoefs[keep] + beta[1])
      ptilde <- as.matrix( as.matrix(Xreg[, keep] %*% allcoefs[keep]) + rep(beta[1], nrow(Xreg)) )
      phat_temp <- 1/(1 + exp( -c(ptilde) ))    # new probabilities
      phat[keep_mat1] <- phat_temp
      yhat_temp <- sapply(phat_temp, function(p) sample(c(0,1), 1, prob=c(1-p,p)))
      Yhat[keep_mat1] <- yhat_temp     # new matrix
      
      # Remove new ratifications from future entries 
      if(t < tfinal){
        NAnew <- which(Yhat[,,k,i] == 1, arr.ind=T)    # new signatures
        # temp <- as.matrix(cbind(NAnew, t, 1))
        # rownames(temp) <- 1:nrow(temp)
        # colnames(temp) <- c("i", "j", "t", "ratification_year")
        # Ytemp <- rbind(Ytemp, temp)   # save new signatures
        for(l in (k+1):dim(Yhat)[3]){
          NAentries <- as.matrix(cbind(NAnew, l, i))
          Yhat[NAentries] <- phat[NAentries] <- NA   # set all future i,j ratifications to NA
        }
        NAremove <- rbind(NAremove, NAnew)   
      }
      
      # Make new X array to predict from
      if(t < tfinal){
        nextyear <- unique(Y$year[Y$t == t+1])   # year for NEXT time period
        Xpred <- update_X_ijt(X, Yhat[,,,i], Yold, t+1, nextyear, region_income)  # covariates Xpred for NEXT time step
      }
      
      # Make new D array to predict from in NEXT time step
      if(t < tfinal){
        Dpred <- Dshell   # shell in which to save
        
        # save indices of which ratifications were signed in previous years
        signs <- NAnew   # signatures from last time period
        
        # If lag is greater than 1, need to augment Yhat with previous signatures
        if(lag > 1){
          if (t < t0+lag-1){   # if need to consult input Y for some signatures
            
            # split indices into those before prediction and after
            lagrange <- lag:1    
            keept <-  (t-lag+1):(t) < t0
            lold <- lagrange[keept]
            lhat <- lagrange[!keept]
            
            # augment with previous signatures from Y (i.e. not simulated)
            temp <- Y[Y$ratification_year==1 & Y$t %in% ((t-lag+1):(t))[keept], c("i", "j")]
            if(length(temp) > 0){
              signs <- rbind(signs, as.matrix(unique(temp)))  # old signatures
            }
            
            for(l in (lhat-1)){   # augment with previous signatures from Yhat, -1 accounts for fact that we are using this for next time period
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))   # updated signatures
            }
            
            
          } else {   # can use Yhat exclusively for signatures
            for(l in 1:(lag-1)){   # augment with previous signatures
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))
            }
          }
        }
        
        signs <- unique(signs)   # remove any duplicates
        Dpred[signs[,1] + (signs[,2] - 1)*S, response] <- 1   # save signatures regardless of lag
      } 
    }
    
    if(i %% write_interval == 0){
      save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
           file=file.path(outdir, filename))
      if(verbose){
        cat("done with sim", i, ";   ")
      }
    }
    
    
  }
  
  # Write out a final time
  save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
       file=file.path(outdir, filename))
  if(verbose){
    cat("\n*********************************\n")
    cat("DONE; saved to ", outdir, "\n")
  }
  
  return(list(Yhat=Yhat, phat=phat, S=S, L=L, countries=countries, treaties=treaties, year_range=year_range,trange=trange, response=response, seed0=seed0, NApairs=NApairs))
}







# update X covariates that depend on past ratifications, helper function for roll_forward_predict()
#   ijrats are i,j pairs that ratified at t
#   t is year for which X is updating to, and year is corresponding year
#   Y is an array from Yhat and Yold is a matrix of old data before t
update_X_ijt <- function(X, Y, Yold, t, year, region_income)   
{
  if(length(dim(Y)) != 3){ stop("Y must be a 3-mode array")}
  kold <- which(as.numeric(dimnames(Y)[[3]]) < year)
  
  Xnew <- X[X$t == t,]   # subset to t only to return
  Xold <- X[X$t == t-1,]   # previous step
  region_income_old <- region_income[region_income$year == year-1,]    # region and income of countries in t-1 year
  region_income_new <- region_income[region_income$year == year,]    # region and income of countries in t year
  
  Xnew$lagpercentincome <- Xnew$lagthreshold <- Xnew$lagpercentregion <- NA  # set all to NA to begin
  
  allregions <- unique(region_income$region)   # ALL regions of interest
  L <- dim(Y)[2]   # number of treaties
  S <- dim(Y)[1]   # number of countries
  allincomes <- 0:2    # unique incomes are same in every year
  
  treaties <- as.numeric(unique(dimnames(Y)[[2]]))   # unique treaties of interest, in order
  countries <- as.numeric(unique(dimnames(Y)[[1]]))   # unique countries of interest, in order
  cows2old <- lapply(1:L, function(z) unique(Yold$cowcode[Yold$j == z & Yold$t < t]))  # countries that have the opportunity sign each treaty in Yold 
  
  # compute lagthreshold based on ratifications
  ratsold <- sapply(1:L, function(z) sum(Yold$ratification_year[Yold$j == z & Yold$t < t], na.rm=T))  # ratifications for all years before t of each treaty
  ratsnew <- sapply(1:L, function(z) sum(Y[,z,kold], na.rm=T))  # ratifications for all years up to t-1
  ratsold[is.na(ratsold)] <- 0  # set NAs to zero
  ratsnew[is.na(ratsnew)] <- 0
  rats <- ratsold + ratsnew
  Xnew$lagthreshold <- rats[match(Xnew$j, 1:L)]    # updated ratifications, lagged by 1!
  
  cows1_byincome <- lapply(allincomes, function(z) unique(region_income_old$cowcode[region_income_old$income == z]))
  cows1_byregion <- lapply(allregions, function(z) unique(region_income_old$cowcode[region_income_old$region == z]))
  
  
  for(j in 1:L){
    cows2new_i <- unique( which(!is.na(Y[,j,kold, drop=F]), arr.ind=T)[,1] )   # i-country indices that had the opportunity to sign the given treaty
    cows2new <- countries[cows2new_i]   # countries that had the opportunity to ratify the given treaty in Y
    cows2 <- union(cows2old[[j]], cows2new)   # countries from any time < t that had the chance to sign that particular treaty
    cows2 <- union(cows2, c(265, 955))
    if(j == 14){   # account for treaties 40439 and 40434 signed in year 1950
      cows2 <- union(cows2, 390)
    }
    if(j == 16){   # account for treaties 40439 and 40434 signed in year 1950
      cows2 <- union(cows2, c(2, 94))
    }
    
    # Lag percent income updates
    for(i in 1:length(allincomes)){   # income code
      income <- as.numeric(allincomes[i])  # income in question
      cows1 <- cows1_byincome[[i]]
      
      cows <- intersect(cows1, cows2)   # countries in that particular income group in the previous year that had the chance to sign that particular the treaty
      cows <- cows[!is.na(cows)]   # remove any NAs
      cows_i <- (1:S)[match(cows, countries)]
      if(length(cows) > 0){  # save if there are any cows to update
        rowsold <- which(Yold$j == j & Yold$t < t & Yold$cowcode %in% cows)  # pertinent entries in Yold for previous years
        iold <- as.matrix( expand.grid(cows_i,j,kold) )   # pertinent entries in Yarray
        rowsnew <- which(Xnew$j == j & Xnew$t == t & Xnew$i %in% cows_i)  # pertinent rows for this year
        if(length(rowsnew) > 0){
          saveit <- length(rowsnew) > 0 | nrow(iold) > 0
          if(saveit){
            temp <- sum(Yold$ratification_year[rowsold], na.rm=T) + sum(Y[iold], na.rm=T)     # old and new ratifications
            Xnew$lagpercentincome[rowsnew] <- temp / length(cows)*100  # recalculate
          } else {
            Xnew$lagpercentincome[rowsnew] <- 0   # set to zero if no intersection
          }
        }
      }
    }
    
    
    # Lag percent region updates
    for(i in 1:length(allincomes)){   # income code
      cows1 <- cows1_byregion[[i]]
      
      cows <- intersect(cows1, cows2)   # countries in that particular income group in the previous year that had the chance to sign that particular the treaty
      cows <- cows[!is.na(cows)]   # remove any NAs
      cows_i <- (1:S)[match(cows, countries)]
      if(length(cows) > 0){  # save if there are any cows to update
        rowsold <- which(Yold$j == j & Yold$t < t & Yold$cowcode %in% cows)  # pertinent entries in Yold for previous years
        iold <- as.matrix( expand.grid(cows_i,j,kold) )   # pertinent entries in Yarray
        rowsnew <- which(Xnew$j == j & Xnew$t == t & Xnew$i %in% cows_i)  # pertinent rows for this year
        if(length(rowsnew) > 0){
          saveit <- length(rowsnew) > 0 | nrow(iold) > 0
          if(saveit){
            temp <- sum(Yold$ratification_year[rowsold], na.rm=T) + sum(Y[iold], na.rm=T)     # old and new ratifications
            Xnew$lagpercentregion[rowsnew] <- temp / length(cows)*100  # recalculate
          } else {
            Xnew$lagpercentregion[rowsnew] <- 0   # set to zero if no intersection
          }
        }
      }
    }
    
    
  }
  
  return(Xnew)
}




##########################
###  Helper Functions  ###
##########################

# Load unprocessed Y array to 
load_ytemp <- function(path)
{
  load(path)
  return(Ytemp)
}



# Function to build regression design matrix for additive models
# D is array to regress Y on for A and B, additive models only!
# X is array of covariates to regress Y upon for beta
# type is "biten" or "sadd"
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
# 
#
# Returns design matrix
#    includes intercept if there is one in X (i.e. all covariates in X are included)
#    returns S*L*tmax \times S^2 + L^2 + ncol(X) matrix, with columnwise-vectorized order of rows
build_design_additive <- function(D, X, type="biten", use_cov=T, sparsedata=F, write=T, response="ratification_year", S=NULL, L=NULL, tmax=NULL)
{
  
  if(!sparsedata){
    # Check size
    if(length(dim(D)) != 3 ){ stop("D is not a 3-mode array") }
    if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0 & use_cov){  stop("Dimenions of D and X don't match")}
    
    # Find sizes
    if(is.null(S) & is.null(L) & is.null(tmax)){
      S = nrow(D[,,1])
      L = ncol(D[,,1])
      tmax = dim(D)[3]
    }

    
    # Build X matrix
    if(use_cov){   # start with beta columns base on use_cov flag
      p <- dim(X)[4]
      Xreg <- matrix(0, S*L*tmax, S^2 + L^2 + p)  # initialize X
      Xreg[, S^2 + L^2 + 1:p] <- t(mat(X, 4))  # beta columns
    } else {  
      p <-  0
      Xreg <- matrix(0, S*L*tmax, S^2 + L^2 + p)  # initialize X
    } 
    
    if(strtrim(type, 3) == "sad"){
      Js <- matrix(1, S, S)
      Jl <- matrix(1, L, L)
    } else if (strtrim(type, 3) == "bit"){
      Js <- diag(S)
      Jl <- diag(L)
    } else { stop("Invalid model type") }
    
    for(t in 1:tmax){  # vec A then vec B
      Xreg[ 1:(S*L) + (t - 1)*S*L, 1:S^2] <- kronecker(t(D[,,t] %*% Jl), diag(S))   # A columns
      Xreg[ 1:(S*L) + (t - 1)*S*L, S^2 + 1:L^2] <- kronecker(diag(L), Js %*% D[,,t])    # B columns
    }
  
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(D)) | !("j" %in% names(D)) | !("t" %in% names(D))){stop("D must have column names i,j, and t")}
    if( !("i" %in% names(X)) | !("j" %in% names(X)) | !("t" %in% names(X))){stop("X must have column names i,j, and t")}
    
    if(strtrim(type, 3) == "sad"){stop("SADD sparse not implemented")} 
    else if (strtrim(type, 3) == "bit"){
      
      filename <- paste0("Xreg_", type, "")
      
      # Find sizes
      if(is.null(S) & is.null(L) & is.null(tmax)){
        S <- max(D$i)
        L <- max(D$j)
        tmax <- max(D$t)
      }
      
      numcols <- S^2 + L^2 
      if(use_cov){
        p <- ncol(X) - 3   # remove i,j,t values
        numcols <- numcols + p
        
        X <- X[X$i <= S & X$j <= L & X$t <= tmax,]
      }
      
      # Xreg <- sparseMatrix(i=1,j=1,x=0, dims=c(tmax*S*L, numcols))   # initialize
      onerows <- which(as.vector(D[,response]) == 1)    # rows of D that have 1s in response
      reg1s <- matrix(0, length(onerows)*(S+L), 2)   # indices in Xreg that are 1s
      count <- 0
      for(k in onerows){
        count <- count+1
        # Xreg[cbind((D$j[k]-1)*S + 1:S + (D$t[k]-1)*S*L, (D$i[k]-1)*S + 1:S)] <- rep(1, S)   # A portion
        # Xreg[cbind((0:(L-1))*S + D$i[k] + (D$t[k]-1)*S*L, S^2 + (0:(L-1))*L + D$j[k])] <- rep(1, L)   # B portion
        i <- D$i[k]    ;   j <- D$j[k]   ;   t <- D$t[k]
        reg1s[1:S + (S+L)*(count-1),] <- cbind(S*L*(t-1) + S*(j -1) + 1:S, S*(i-1) + 1:S)
        reg1s[S + 1:L + (S+L)*(count-1),] <- cbind(S*L*(t-1) + i + S*(0:(L-1)), S^2 + j + L*(0:(L-1)))
      }
      Xreg <- sparseMatrix(i=reg1s[,1],j=reg1s[,2], x=1, dims=c(tmax*S*L, numcols))   # initialize
      if(use_cov){
        keep <- which(!(names(X) %in% c("i","j","t")))  # columns to keep that aren't i,j, or t
        Xreg[X$i + (X$j-1)*S + (X$t-1)*S*L, S^2 + L^2 + 1:p] <- as.matrix(X[, keep])
      }
      
    } else { stop("Invalid model type") }
    
  } else { stop("sparsedata must be true/false")}
  
  return(Xreg)
}



# Build X for all ijt triples starting at t0 and ending at tfinal
#  contains rows for all ijt triples, not just thoes in the dataset
build_big_X <- function(t0, X, S=NULL, L=NULL, tfinal=NULL, readfile=F, wd=NULL)
{
  filename <- "Xbig.RData"
  if(readfile){
    if(!is.null(wd)){
      if(filename %in% list.files(wd)){
        load(filename)
        cat("Read-in big X \n")
        return(Xbig)
      } else {
        warning("Did not find Xbig.RData in wd... building manually")
      }
    } else {
      if(filename %in% list.files()){
        load(filename)
        cat("Read-in big X \n")
        return(Xbig)
      } else {
        warning("Did not find Xbig.RData in current working directory... building manually")
      }
    }
  }
  
  # Check data
  if(!all(c("i", "j", "t") %in% names(X))){stop("Need columns named i,j, and t")}
  if(!((t0-1) %in% unique(X$t))){stop("X must contain at least time period t-1")}
  
  # Build data size
  if(is.null(S)){S <- max(X$i)}
  if(is.null(L)){L <- max(X$j)}
  if(is.null(tfinal)){tfinal <- max(X$t)}
  Xsave <- X[X$t >= t0, ]   # only t greater than t0-1
  dt <- tfinal + 1 - t0   # number of time periods to save
  
  # Initialize new X
  Xnew <- X[0,]  # new X with same columns
  Xnew[1:(S*L*dt),] <- NA    # row size of new X
  Xnew[, c("i", "j", "t")] <- which(array(0, c(S,L,dt)) == 0, arr.ind=T) # ijt indices
  Xnew$t <- Xnew$t + t0 - 1   # increase t index to appropriate range, MUST DO THESE STEPS as relies on i,j,t presence
  Xnew$intercept <- 1   # intercept
  
  # Save existing entries in X
  rowsX <- Xsave$i + (Xsave$j-1)*S + (Xsave$t - t0)*S*L   # locations of X values in Xnew
  Xnew[rowsX,] <- Xsave   # save all old X values in larger data frame
  
  
  # entries in X that only vary by i,t
  irange <- c("polity", "lopen", "lnso2pc", "lgrgdpc", "memberships", "lnGDP")   
  # entries in X that only vary by j,t
  jrange <- c("hard_law2",  "global", "gdmix",  "ass_all", "ass_dev")   
  # entries in X that only vary by i,j,t
  ijtrange <- c("t_new", "t2_new", "t3_new")   # don't actually use these two, just notes for me
  ijonly <- c("leg_ap")
  
  # fill in i and j indices
  for(t in t0:tfinal){   # loop through time periods
    
    # i indices
    for(i in 1:S){    
      ii <- which(Xnew$i == i & Xnew$t == t)   # relevant rows
      itemp <- unique(Xnew[ii,irange])   # all rows should be the same except those that are all NA
      if(all(is.na(itemp))){
        itemp <- matrix(NA, 1, length(irange))   # if everything is NA, keep them that way
      } else {
        NArows <- apply(itemp, 1, function(z) all(is.na(z)))   # otherwise remove the NA row
        itemp <- itemp[!NArows,]
      }
      if(nrow(itemp) != 1){warning("There is some issue with filling in values that are constant among i,t")}   # if saving more than one row something went wrong
      
      NAii <- which(apply( Xnew[ii,irange], 1, function(z) all(is.na(z))))
      Xnew[ii[NAii],irange] <- itemp   # saving borrowed across all relevant pairs! only to those rows that are all NA
    }
    
    # j indices
    for(j in 1:L){    
      jj <- which(Xnew$j == j & Xnew$t == t)   # relevant rows
      jtemp <- unique(Xnew[jj,jrange])   # all rows should be the same except those that are all NA
      if(all(is.na(jtemp))){
        jtemp <- matrix(NA, 1, length(jrange))   # if everything is NA, keep them that way
      } else {
        NArows <- apply(jtemp, 1, function(z) all(is.na(z)))   # otherwise remove the NA row
        jtemp <- jtemp[!NArows,]
      }
      if(nrow(jtemp) != 1){warning("There is some issue with filling in values that are constant among j,t")}   # if saving more than one row something went wrong
      
      NAjj <- which(apply( Xnew[jj, jrange], 1, function(z) all(is.na(z))))
      Xnew[jj[NAjj], jrange] <- jtemp   # saving borrowed across all relevant pairs! only to those rows that are all NA    
    }
  }
  
  # fill in ij indices
  ijmin <- ijmax <- NULL
  for(i in 1:S){
    for(j in 1:L){
      ij <- which(Xnew$i == i & Xnew$j == j)
      NAij <- which(is.na(Xnew$t_new[ij]))
      ij0 <- which(Xnew$t_new[ij] == 0)
      
      # Work on zeros
      if(length(NAij) > 0){   # if there are possibly entries to replace
        if(length(ij0) > 0){  # if there is a zero 
          if(min(Xnew$t[ij[NAij]]) < Xnew$t[ij[ij0]]){   # if therre are NA entries before the zero
            ijreplace <- which( Xnew$t[ij[NAij]] < Xnew$t[ij[ij0]] )   # indices to replace
            Xnew$t_new[ij[ijreplace]] <- 0   # zeros before the zero
          }
        }
      }
      
      # Work on NAs beyond zeros
      NAij <- which(is.na(Xnew$t_new[ij]))
      if(length(NAij) > 0){   # if there are still entries to replace
        if(length(NAij) != length(ij)){   # if not all NAs to replace
          tm <- max(which(!is.na(Xnew$t_new[ij])))  # max index
          Xnew$t_new[ij[(tm+1):length(ij)]] <- Xnew$t_new[ij[tm]] + 1:length( (tm+1):length(ij) )
        } else {
          pasttemp <- X$t[X$i==i & X$j==j]  # go back to X for the ones greater then zero
          if(sum(is.na(pasttemp)) != length(pasttemp)){   # if not all past values are NAs
            t1 <- max(pasttemp, na.rm=T)   # reference time point
            maxt <- max(Xnew$t[ij], na.rm=T)
            ijreplace <- which(Xnew$t[ij] > maxt)   # replace those NAs for which t is larger
            Xnew$t_new[ij[ijreplace]] <- (X$t_new[X$i == i & X$j == j & X$t== t1] + Xnew$t[ij] - t1)[ijreplace]
          }
        }
      }
      
      
      # leg_ap
      latemp <- unique(X$leg_ap[X$i==i & X$j == j])  # possible legistlative approvals
      if(length(latemp) > 0){
        latemp1 <- latemp[!is.na(latemp)]  # non-NAs
        if(length(latemp1) > 1){
          warning("There exists at least one leg_ap variable that depends on t")
          cat("leg_ap: i",i,"j",j,"\n")
        }
        Xnew$leg_ap[ij] <- latemp1   # write regardless of whether there is more than 1
      }
      
      
    }
  }
  
  # Assign squared and cubed times
  Xnew$t2_new <- Xnew$t_new^2
  Xnew$t3_new <- Xnew$t_new^3
  
  
  return(Xnew)
}



# simple, direct auc
#  https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
simple_roc_auc <- function(p,l)
{
  x = p;  y = l
  x1 = x[y==1]; n1 = length(x1); 
  x2 = x[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
  auc
}

# simple, direct auc of precision/recall curve
simple_pr_auc <- function(p,l)
{
  if(sum(l)==0){stop("Labels are all zero. Need at least one 1. ")}
  
  l <- l[order(p, decreasing=T)]
  p <- p[order(p, decreasing=T)]
  tp <- cumsum(l)
  n <- length(l)
  n1 <- sum(l)
  n0 <- n - n1
  np <- 1:n
  
  prec <- tp/np
  # prec[is.infinite(prec)] <- 0
  rec <- tp/n1
  
  prec <- prec[order(rec)]
  rec <- rec[order(rec)]
  
  dx <- c(rec[1], rec[2:n] - rec[1:(n-1)])
  if(rec[1] !=0 ){
    prec_traps <- c(prec[1], .5*(prec[1:(n-1)] + prec[2:n]))
  } else {
    prec_traps <- c(prec[1]/2, .5*(prec[1:(n-1)] + prec[2:n]))
  }
  
  auc <- sum(prec_traps*dx)
  auc
}




#### Courtesy of Dr. Hoff below

#' Matricization
#' 
#' Matricize an array
#'
#' This functions matricizes an array along a given mode. 
#'
#' @param A an array
#' @param k a mode of \code{A}, along which to matricize
#' @keywords arrays matrices
#' @export
#' @examples
#' A<-rsan(4,3,2)
#' mat(A,2) 
mat<-function(A,k)
{
  Ak<-t(apply(A,k,"c"))
  if(nrow(Ak)!=dim(A)[k])  { Ak<-t(Ak) }
  Ak
}


#' Array-matrix product
#'
#' Multiply an array by a matrix along a given mode
#'
#' This function multiplies a matricized array by another 
#' matrix, and then reforms the result into a new array. 
#'
#' @param A a real valued array 
#' @param M a real matrix
#' @param k an integer, a mode of \code{A}
#' @author Peter Hoff
#' @keywords arrays
#' @export
#' @examples
#' A<-rsan(c(5,4,3))
#' B<-rsan(c(2,5))
#' amprod(A,B,1)
amprod<-function(A,M,k)
{
  K<-length(dim(A))
  AM<-M%*%mat(A,k)
  AMA<-array(AM, dim=c(dim(M)[1],dim(A)[-k]) )
  aperm(AMA,  match(1:K,c(k,(1:K)[-k]) ) )
}








##################
roll_forward_predict_v2 <- function(t0, nsims, A, B, beta, Y, D, X, lag, 
                                 region_income, tfinal=NULL, response="ratification_year", seed0=1, verbose=F, 
                                 NApairs=NULL, write_interval=ceiling(nsims/10), outdir=getwd(), 
                                 filename=NULL, datadir=NULL)  # countries=NULL, treaties=NULL, 
{
  if(is.null(tfinal)){ tfinal <- max(Y$t)}
  trange <- t0:tfinal
  year_range <- unique(Y$year[Y$t == t0] ) : unique(Y$year[Y$t == tfinal] )
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  beta <- as.matrix(beta)
  
  # Initialize arrays
  Ynew <- Y[Y$t >= t0 & Y$t <= tfinal,]
  Yold <- Y[Y$t < t0,]
  Dnew <- D[D$t >= t0 & D$t <= tfinal,]
  Xnew <- X[X$t >= t0 & X$t <= tfinal,]
  S <- max(Dnew$i)  ;   L <- max(Dnew$j)    
  Yhat <- phat <- array(0, c(S,L,length(trange),nsims))
  countries <- sapply(1:S, function(z) unique(Dnew$cowcode[Dnew$i == z]))   # all countries of interest in future
  treaties <- sapply(1:L, function(z) unique(Dnew$treaty[Dnew$j == z]))    # all treaties of interest in future
  dimnames(Yhat)[[1]] <- dimnames(phat)[[1]] <- countries
  dimnames(Yhat)[[2]] <- dimnames(phat)[[2]] <- treaties
  dimnames(Yhat)[[3]] <- dimnames(phat)[[3]] <- year_range
  region_income$i <- Y$i[match(region_income$cowcode, Y$cowcode)]   # add column for i in region_income indicators 
  region_income <- region_income[order(region_income$i),]    # reorder
  
  # Read in big X if X is not big
  if(nrow(X) <  S*L*length(trange)){
    X <- build_big_X(t0, X, S=S, L=L, tfinal=max(X$t), readfile=T)
  }
  
  
  
  
  load(file.path(datadir, "pre-prediction.RData"))
  Yold2 <- Yold
  
  
  # build possible signatures for all future times and simulations
  for(t in trange){
    k <- t - t0 + 1
    iposs <- cbind(possibles[rep(1:nrow(possibles), times=nsims), ], k, rep(1:nsims, each=nrow(possibles)))   # indices in first year that can be signed
    Yhat[,,k,] <- phat[,,k,] <- NA   # all NAs in time period
    Yhat[iposs] <- phat[iposs] <- 0   # possible ratifications
    
    jlate_entry <- which(t >= late_entry1[,3])   # AFTER late entry
    if(length(jlate_entry) > 0){
      # cat("late entry k=", k, "t=", t, "\n")
      ilate <- cbind(late_entry1[rep(jlate_entry, times=nsims), 1:2], k, rep(1:nsims, each=length(jlate_entry)) )
      Yhat[ilate] <- phat[ilate] <- 0   # 0s for i,j pairs that enter late
    }
    
    jearly_exit <- which(t > early_exit[,3])   # AFTER exit
    if(length(jearly_exit) > 0){
      # cat("early exit k=", k, "t=", t,"\n")
      iearly <- cbind(early_exit[rep(jearly_exit, times=nsims), 1:2], k, rep(1:nsims, each=length(jearly_exit)) )
      Yhat[iearly] <- phat[iearly] <- NA   # NAs for i,j pairs that leave early
    }
  }
  
  # Save all ratifications for future updates
  NAinit1 <- already_signed
  #unique(Yold[Yold$ratification_year == 1, c("i", "j")])   # save already signed treaty/country pairs
  if(!is.null(NApairs)){
    NAinit1 <- rbind(NAinit1, NApairs)
  }
  
  
  # Prediction
  allcoefs <- c(c(t(as.matrix(A))), c(as.matrix(B)), c(as.matrix(beta))[-1])   # coefs without intercept
  allcoefs[is.na(allcoefs)] <- 0
  keep <- which(allcoefs != 0)
  
  # Make dummy D and X variable arrays to avoid rebuilding
  #   shell is all pairs of i,j such that i<S and j<L for a single time period
  Dshell <- D[0,]
  Dshell[1:(S*L),] <- 0
  Xshell <- X[0,]
  Xshell[1:(S*L),] <- 0
  Dshell$t <- Xshell$t <- 1
  Dshell$i <- Xshell$i <- rep(1:S, times=L)
  Dshell$j <- Xshell$j <- rep(1:L, each=S)
  Dshell$cowcode <- Y$cowcode[match(Dshell$i, Y$i)]
  Dshell$treaty <- Y$cowcode[match(Dshell$j, Y$j)]
  Xshell$intercept <- 1
  Xshell[, -c(1:4)] <- NA   # NAs for all in Xshell
  
  # filename to save
  # if(is.null(filename)){
  filename <-  paste0("predict_lag", lag, "_", min(year_range),"_", max(year_range), "_seed", seed0, "_", ps,  ".RData")
  # }
  # 
  
  # Run loop to make predictions, write out results periodicially
  for(i in 1:nsims){
    set.seed(seed0+i-1)   # set seed for repeatability
    NAremove <- NAinit1   # initialize country/treaty pairs to remove from dataset for each simulation
    # Ytemp <- Y[Y$t < t0 & Y$ratification_year==1, c("i", "j", "t", "ratification_year")]   # save all ratifications
    
    for(t in t0:tfinal){    # t is the year of the prediction.
      k <- t - t0 + 1   # index in vavriables to save
      
      if(t == t0){    # initialize possibly lagged autoregressive array and covariate array
        Dpred <- Dnew[Dnew$t == t,] 
        Xpred <- Xnew[Xnew$t == t,]
      }   
      
      # Set impossible signatures to NA (should already be, but double-check)
      Yhat[cbind(NAremove, k, i)] <- phat[cbind(NAremove, k, i)] <- NA
      
      # Build design matrix for appropriate D, X
      Xpred$t <- Dpred$t <- 1
      Xreg <- build_design_additive(Dpred, Xpred, sparsedata=T, write=F, S=S, L=L, tmax=1)
      
      # Subset and pare design matrix to predict
      Xreg <- Xreg[,-(S^2 + L^2 + 1)]   # remove intercept
      keep_mat <- which(!is.na(Yhat[,,k,i]), arr.ind=T)   # i,j pairs to predict, in columnwise order 
      rows <- keep_mat[,1] + (keep_mat[,2]-1)*S   #  + (Y$t-1)*S*L   # unfolded indices
      Xreg <- Xreg[sort(rows), ]   # rows of Xreg that pertain to the entres in Yhat
      Xreg[which(is.na(Xreg))] <- 0   # Set NAs to zero
      
      # Calculate new prediction probabilities and predict
      keep_mat <- keep_mat[order(rows),]   # order entries as calculated
      keep_mat1 <- as.matrix(cbind(keep_mat, k, i))   # including time and sim indices
      # ptilde <- as.matrix(Xreg[, keep] %*% allcoefs[keep] + beta[1])
      ptilde <- as.matrix( as.matrix(Xreg[, keep] %*% allcoefs[keep]) + rep(beta[1], nrow(Xreg)) )
      phat_temp <- 1/(1 + exp( -c(ptilde) ))    # new probabilities
      phat[keep_mat1] <- phat_temp
      yhat_temp <- sapply(phat_temp, function(p) sample(c(0,1), 1, prob=c(1-p,p)))
      Yhat[keep_mat1] <- yhat_temp     # new matrix
      
      # Remove new ratifications from future entries 
      if(t < tfinal){
        NAnew <- which(Yhat[,,k,i] == 1, arr.ind=T)    # new signatures
        # temp <- as.matrix(cbind(NAnew, t, 1))
        # rownames(temp) <- 1:nrow(temp)
        # colnames(temp) <- c("i", "j", "t", "ratification_year")
        # Ytemp <- rbind(Ytemp, temp)   # save new signatures
        for(l in (k+1):dim(Yhat)[3]){
          NAentries <- as.matrix(cbind(NAnew, l, i))
          Yhat[NAentries] <- phat[NAentries] <- NA   # set all future i,j ratifications to NA
        }
        NAremove <- rbind(NAremove, NAnew)   
      }
      
      # Make new X array to predict from
      if(t < tfinal){
        nextyear <- unique(Y$year[Y$t == t+1])   # year for NEXT time period
        Xpred <- update_X_ijt(X, Yhat[,,,i], Yold2, t+1, nextyear, region_income)  # covariates Xpred for NEXT time step
        # Xpred <- update_X_ijt_2(X, Yhat[,,,i], t+1, nextyear, region_income, cows2old, ratsold, Yold2)   
      }
      
      
      
      # Make new D array to predict from in NEXT time step
      if(t < tfinal){
        Dpred <- Dshell   # shell in which to save
        
        # save indices of which ratifications were signed in previous years
        signs <- NAnew   # signatures from last time period
        
        # If lag is greater than 1, need to augment Yhat with previous signatures
        if(lag > 1){
          if (t < t0+lag-1){   # if need to consult input Y for some signatures
            
            # split indices into those before prediction and after
            lagrange <- lag:1    
            keept <-  (t-lag+1):(t) < t0
            lold <- lagrange[keept]
            lhat <- lagrange[!keept]
            
            # augment with previous signatures from Y (i.e. not simulated)
            temp <- Y[Y$ratification_year==1 & Y$t %in% ((t-lag+1):(t))[keept], c("i", "j")]
            if(length(temp) > 0){
              signs <- rbind(signs, as.matrix(unique(temp)))  # old signatures
            }
            
            for(l in (lhat-1)){   # augment with previous signatures from Yhat, -1 accounts for fact that we are using this for next time period
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))   # updated signatures
            }
            
            
          } else {   # can use Yhat exclusively for signatures
            for(l in 1:(lag-1)){   # augment with previous signatures
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))
            }
          }
        }
        
        signs <- unique(signs)   # remove any duplicates
        Dpred[signs[,1] + (signs[,2] - 1)*S, response] <- 1   # save signatures regardless of lag
      } 
    }
    
    if(i %% write_interval == 0){
      save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
           file=file.path(outdir, filename))
      if(verbose){
        cat("done with sim", i, ";   ")
      }
    }
    
    
  }
  
  
  outdir <- "~/Dropbox/BiTEN/reproduce_envtreat/results/process_counterfactual_check_3"
  
  # Write out a final time
  save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
       file=file.path(outdir, filename))
  if(verbose){
    cat("\n*********************************\n")
    cat("DONE; saved to ", outdir, "\n")
  }
  
  return(list(Yhat=Yhat, phat=phat, S=S, L=L, countries=countries, treaties=treaties, year_range=year_range,trange=trange, response=response, seed0=seed0, NApairs=NApairs))
}




