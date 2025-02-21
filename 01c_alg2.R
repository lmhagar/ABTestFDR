## These functions implement Algorithm 2 from the paper for an arbitrary 
## example with multiple binary outcomes; a function to construct bootstrap
## confidence intervals is also provided

## load required libraries
library(foreach)
library(doParallel)
library(doSNOW)

## this function outputs the optimal sample size and decision threshold
## along with some additional information used to later construct bootstrap
## confidence intervals and contour plots; these outputs are explained right 
## before the return function. This function uses a multinomial model to
## generate and analyze data

algorithm2Multi <- function(par.a, par.b, q, target.pwr, deltas, 
                            mat, hyper, n0, c = 1, seed = 1, 
                            setup.par = FALSE, contour = FALSE, boot = FALSE){
  
  ## the inputs are described as follows
  ## par.a: a matrix of multinomial probabilities for group A; each row corresponds to a simulation
  ##        repetition and the parameters from the same submodel should be in consecutive rows
  ## par.b: a matrix of multinomial probabilities for group B; each row corresponds to a simulation
  ##        repetition and the parameters from the same submodel should be in consecutive row
  ## q: an upper bound for the false discovery rate
  ## target.pwr: a lower bound for the average power across the true hypotheses
  ## deltas: a 2-column matrix that corresponds to the delta_L and delta_U interval endpoints
  ##         for each hypothesis; the first column corresponds to the lower interval endpoint 
  #          and the second one corresponds to the upper interval one
  ## mat: the matrix returned by createTree()
  ## hyper: a 4 x (# of categories -1) matrix of hyperparameters for the independent beta distributions that  
  ##        correspond to the multinomial model. We assume a generalized Dirichlet distribution as the prior: 
  ##        https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution. This prior consists of independent
  ##        beta distributions for conditional multinomial probabilities; the first row is the alpha parameters
  ##        for these distributions in group A, and the second row is the beta parameters in group A. The third
  ##        row is the alpha parameters for group B, and the fourth row is the beta parameters for group B.
  ## n0: the initial sample size
  ## c: the constant for imbalanced sample sizes (1 by default)
  ## seed: numeric seed for reproducibility
  ## setup.par: a binary indicator denoting whether the parallelization should be set up inside the function
  ##            if FALSE (the default), parallelization should be set up before the function call
  ## contour: a binary indicator denoting whether to return the linear approximations to construct contour
  ##          plots (FALSE by default)
  ## boot: a binary indicator denoting whether to return the estimated logits to construct bootstrap
  ##       confidence intervals (FALSE by default)
  
  ## check the inputs for errors
  if (ncol(par.a) <= 1){
    stop("par.a must contain probabilities for at least two categories")
  } else if (mean(is.numeric(par.a)) < 1){
    stop("par.a must be a numeric matrix")
  } else if (mean(par.a < 0) > 0){
    stop("probs in par.a must be positive")
  } else if (sum(round(rowSums(par.a), 3) != 1) > 0){
    stop("probs in each row of par.a must sum to 1")
  } else if (ncol(par.b) <= 1){
    stop("par.b must contain probabilities for at least two categories")
  } else if (mean(is.numeric(par.b)) < 1){
    stop("par.b must be a numeric matrix")
  } else if (mean(par.b < 0) > 0){
    stop("probs in par.b must be positive")
  } else if (sum(round(rowSums(par.b), 3) != 1) > 0){
    stop("probs in each row of par.b must sum to 1")
  } else if (ncol(par.a) != ncol(par.b) | nrow(par.a) != nrow(par.b)){
    stop("par.a and par.b must have the same dimensions")
  } else if (!is.numeric(q)){
    stop("q must be a number")
  } else if (q <= 0 | q >= 1){
    stop("q must be at between 0 and 1")
  } else if (!is.numeric(target.pwr)){
    stop("target.pwr must be a number")
  } else if (target.pwr <= 0 | target.pwr >= 1){
    stop("target.pwr must be at between 0 and 1")
  } else if (nrow(hyper) != 4){
    stop("hyper must have four rows")
  } else if (mean(is.numeric(hyper)) < 1){
    stop("hyper must be a numeric matrix")
  } else if (sum(hyper <= 0) > 0){
    stop("all elements of hyper must be positive")
  } else if (ncol(hyper) != (ncol(par.a) - 1)){
    stop("the dimensions of hyper and probs in par.a and par.b do not match")
  } else if (ncol(deltas) != 2){
    stop("deltas must have two columns")
  } else if (mean(is.numeric(deltas)) < 1){
    stop("deltas must be a numeric matrix")
  } else if (!is.numeric(n0)){
    stop("n0 must be a number")
  } else if (n0 < 1){
    stop("n0 must be at least 1")
  } else if (!is.numeric(c)){
    stop("c must be a number")
  } else if (c < 0){
    stop("c must be at least 0")
  } else if (!is.numeric(seed)){
    stop("seed must be a number")
  } else if (seed <= 0 | (seed%%1 != 0)){
    stop("seed must be a positive integer") 
  } else if (!(setup.par %in% c(TRUE,FALSE))){
    stop("setup.par must be TRUE or FALSE")
  } else if (!(contour %in% c(TRUE,FALSE))){
    stop("contour must be TRUE or FALSE")
  } else if (!(boot %in% c(TRUE,FALSE))){
    stop("boot must be TRUE or FALSE")
  }
  
  ## define functions used to approximate logits of posterior probabilities
  getLogitsSlopes <- function(post1, post2, seed, n1 = n1, n2 = n2,
                              outcomes = NULL, deltas = NULL, mm = 10000) {
    
    ## the inputs are described below
    ## post1: the 2 x (K-1) matrix of posterior parameters for group A returned by PostParams()
    ## post2: the 2 x (K-1) matrix of posterior parameters for group B returned by PostParams()
    ## seed: a numeric seed for reproducibility
    ## n1: sample size for group A
    ## n2: sample size for group B
    ## outcomes: a list where each item is a vector of the multinomial categories that correspond
    ##           to a given outcome; this list is created by our implementation of Algorithm 2
    ## deltas: a 2-column matrix that corresponds to the delta_L and delta_U interval endpoints
    ##         for each hypothesis on the *logarithmic* scale; the first column corresponds to 
    ##         the lower interval endpoint and the second one corresponds to the upper interval one
    ## mm: the number of simulation repetitions used to estimate the posterior probabilities
    
    ## extract the number of categories in each multinomial model (these should be the same)
    K1 <- nrow(post1) + 1
    K2 <- nrow(post2) + 1
    
    set.seed(seed)
    
    ## generate the z variables for group A; these are the conditional multinomial
    ## probabilities that follow beta distributions
    ## initialize the matrix for the unconditional multinomial probabilities
    z1 <- matrix(rbeta(mm*(K1-1), post1[,1], post1[,2]), nrow = K1 - 1, ncol = mm)
    p1 <- matrix(0, nrow = K1, ncol = mm)
    
    ## efficiently convert between the conditional and unconditional probabilities
    cprod_one_min_z1 <- apply(1 - z1, 2, cumprod)
    p1[1:(K1 - 1), ] <- z1 * rbind(1, cprod_one_min_z1[-(K1 - 1), , drop = FALSE])
    p1[K1, ] <- cprod_one_min_z1[K1 - 1, ]
    
    ## generate the z variables for group B and repeat the conversion process
    z2 <- matrix(rbeta((K2-1)*mm, post2[,1], post2[,2]), nrow = K2 - 1, ncol = mm)
    p2 <- matrix(0, nrow = K2, ncol = mm)
    
    cprod_one_min_z2 <- apply(1 - z2, 2, cumprod)
    p2[1:(K2 - 1), ] <- z2 * rbind(1, cprod_one_min_z2[-(K2 - 1), , drop = FALSE])
    p2[K2, ] <- cprod_one_min_z2[K2 - 1, ]

    ## approximate posteriors for the logarithm of lift
    out <- matrix(0, nrow = length(outcomes), ncol = mm)
    denoms <- NULL
    for (i in 1:length(outcomes)){
      inds <- outcomes[[i]]
      ## sum up the multinomial probabilities corresponding to each outcome
      temp2 <- log(pmax(.Machine$double.eps, colSums(p2[inds,])))
      temp1 <- log(pmax(.Machine$double.eps, colSums(p1[inds,])))
      ## calculate the limiting variance related to the Fisher information
      denoms[i] <- var(na.omit(temp1))*n1 + var(na.omit(temp2))*n1
      ## save the posterior draws for the outcome in the for loop
      out[i,] <- temp2 - temp1
    }
    
    if (n1 < 500 | n2 < 500){
      ## create kernel density estimates to compute logits to account for sample-sample
      ## deviations from normality
      logits <- NULL
      for (i in 1:length(outcomes)){
        ## get kernel density estimate for each outcome
        kd <- density(na.omit(out[i,]))
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          ## calculate the complementary posterior probability based on which interval endpoints are finite
          ## we calculated the complementary probabilities because rounding up to 1 is much more prevalent than
          ## underflow
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE)) +
            mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                            max(na.omit(out[i,])) + 1), kd$bw))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE))
        }
        else{
          logits[i] <- mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw))
        }
      }
    } else {
      ## using a limiting normal approximation that is more stable for large sample sizes
      ## with very large logits
      logits <- NULL
      for (i in 1:length(outcomes)){
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE) +
            pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE)
        }
        else{
          logits[i] <- pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
      }
    }
    
    ## calculate the logits based on the complementary probabilities calculated above
    logits <- log(1 - logits) - log(logits)
    
    ## return the logits and denominators for the limiting slopes that are 
    ## to be used later
    return(c(logits, denoms))
  }
  
  ## this function is the same as getLogitsSlopes() except the additional information
  ## to estimate the limiting slopes is not returned; this function is used by our
  ## implementation of Algorithm 2 during the second time we explore the sampling 
  ## distribution of posterior probabilities
  getLogits <- function(post1, post2, seed, n1 = n1, n2 = n2,
                        outcomes = NULL, deltas = NULL, mm = 10000) {
    
    K1 <- nrow(post1) + 1
    K2 <- nrow(post2) + 1
    
    set.seed(seed)
    
    z1 <- matrix(rbeta((K1-1)*mm, post1[,1], post1[,2]), nrow = K1 - 1, ncol = mm)
    p1 <- matrix(0, nrow = K1, ncol = mm)
    
    cprod_one_min_z1 <- apply(1 - z1, 2, cumprod)
    p1[1:(K1 - 1), ] <- z1 * rbind(1, cprod_one_min_z1[-(K1 - 1), , drop = FALSE])
    p1[K1, ] <- cprod_one_min_z1[K1 - 1, ]
    
    z2 <- matrix(rbeta((K2-1)*mm, post2[,1], post2[,2]), nrow = K2 - 1, ncol = mm)
    p2 <- matrix(0, nrow = K2, ncol = mm)
    
    cprod_one_min_z2 <- apply(1 - z2, 2, cumprod)
    p2[1:(K2 - 1), ] <- z2 * rbind(1, cprod_one_min_z2[-(K2 - 1), , drop = FALSE])
    p2[K2, ] <- cprod_one_min_z2[K2 - 1, ]
    
    out <- matrix(0, nrow = length(outcomes), ncol = mm)
    denoms <- NULL
    for (i in 1:length(outcomes)){
      inds <- outcomes[[i]]
      temp2 <- log(pmax(.Machine$double.eps, colSums(p2[inds,])))
      temp1 <- log(pmax(.Machine$double.eps, colSums(p1[inds,])))
      out[i,] <- temp2 - temp1
    }
    
    if (n1 < 500 | n2 < 500){
      logits <- NULL
      for (i in 1:length(outcomes)){
        kd <- density(na.omit(out[i,]))
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE)) +
            mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                            max(na.omit(out[i,])) + 1), kd$bw))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE))
        }
        else{
          logits[i] <- mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw))
        }
      }
    } else {
      logits <- NULL
      for (i in 1:length(outcomes)){
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE) +
            pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE)
        }
        else{
          logits[i] <- pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
      }
    }
    
    logits <- log(1 - logits) - log(logits)
    
    return(logits)
  }
  
  ## define logit function
  logit <- function(x){log(x)-log(1-x)}
  
  ## turn the matrix from createTree() to a list with the probabilities corresponding to 
  ## each outcome
  cat <- apply(mat, 2, function(x){which(x == 1)})
  
  ## extract number of categories, outcomes, and simulation repetitions
  n.cat <- ncol(par.a)
  n.out <- length(cat)
  m <- nrow(par.a)
  
  ## extract hyperparameters from matrix
  alphas1 <- hyper[1,]; betas1 <- hyper[2,]
  alphas2 <- hyper[3,]; betas2 <- hyper[4,]
  
  ## set up parallelization if setup.par == TRUE
  if (setup.par == TRUE){
    cores=detectCores()
    cl <- makeSOCKcluster(cores[1]-1)
    
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  
  ## identify the rows corresponding to each submodel
  ## determine which consecutive rows are not the same in par.a and par.b
  diff.a <- rowSums(apply(par.a, 2, function(x){abs(diff(x))}))
  diff.b <- rowSums(apply(par.b, 2, function(x){abs(diff(x))}))
  
  ## find the rows where either parameters change and create vector 
  ## that indexes the scenario number
  new.scen <- diff.a + diff.b > 0
  scen.vec <- as.numeric(c(1, cumsum(new.scen) + 1))
  
  ## explore the sampling distribution at the first sample size n0
  n <- n0
  samp0 <- foreach(i=1:m, .combine='rbind',
                   .options.snow=opts, .errorhandling = "remove") %dopar% {
                     ## get the numerators for the limiting slopes
                     eff <- NULL
                     numer <- NULL
                     for (j in 1:length(cat)){
                       eff[j] <- log(sum(par.b[i,cat[[j]]])) - log(sum(par.a[i,cat[[j]]]))
                       num.temp <- min((eff[j] - deltas[j,1])^2, (eff[j] - deltas[j,2])^2)
                       numer[j] <- 0.5*(eff[j] > deltas[j,1] & eff[j] < deltas[j,2])*num.temp
                     }
                     
                     ## generate data
                     set.seed(seed + i)
                     y1 <- as.numeric(rmultinom(1, n, par.a[i,]))
                     y2 <- as.numeric(rmultinom(1, round(c*n), par.b[i,]))
                     
                     ## get the posterior parameters for group A
                     hypa <- cbind(alphas1 + head(y1, length(y1) - 1),
                                   betas1 + n - cumsum(head(y1, length(y1) - 1)))
                     
                     ## get the posterior parameters for group B
                     hypb <- cbind(alphas2 + head(y2, length(y2) - 1),
                                   betas2 + round(c*n) - cumsum(head(y2, length(y2) - 1)))
                     
                     # get logits of posterior probabilities and denominators for slopes
                     out.temp <- getLogitsSlopes(hypa, hypb, seed = seed+ i + m, 
                                                 n1 = n, n2 = round(c*n),
                                                 outcomes = cat, deltas = deltas)
                     c(out.temp, numer)
                     
                   }
  
  ## correct any logits that might be +/- infinity (rare)
  logits0 <- samp0[, 1:n.out]
  for (j in 1:ncol(logits0)){
    logits0[,j] <- ifelse(logits0[,j] == -Inf, min(subset(logits0, is.finite(logits0[,j]))[,j]) - 1, logits0[,j])
    logits0[,j] <- ifelse(logits0[,j] == Inf, max(subset(logits0, is.finite(logits0[,j]))[,j]) + 1, logits0[,j])
  } 
  
  ## get limiting slopes to choose next sample size
  denom0 <- samp0[, (n.out+1):(2*n.out)]
  numer0 <- samp0[,(2*n.out+1):(3*n.out)]
  slopes0 <- numer0/denom0
  
  ## get optimal decision threshold for this sample size
  gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
  while (gamma_u - gamma_l > 0.0001){
    ## check if current gamma value controls the FDR
    gamma_m <- 0.5*gamma_u + 0.5*gamma_l
    ## check discoveries
    bin <- logits0 > logit(gamma_m)
    ## this avoids dividing by 0
    zeros <- rowMeans(bin) == 0
    ## only count false discoveries in the numerator; all discoveries in the denominator
    fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
    fdr.temp <- mean(fdrs)
    
    ## update gamma value based on whether FDR criterion is satisified
    if (fdr.temp <= q){
      gamma_u <- gamma_m
    }
    else {
      gamma_l <- gamma_m
    }
  }
  gamma0 <- gamma_m
  
  ## get average for first sample size and gamma combination
  bin <- logits0 > logit(gamma_m)
  zeros <- rowMeans(bin) == 0
  ## only count true discoveries in the numerator
  pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
  pwr.0 <- mean(pwrs)
  
  ## find next sample size based on the limiting slopes
  ## we find a smaller sample size if average power is already large enough
  if (pwr.0 >= target.pwr){
    n.upper <- n
    n.lower <- FALSE
    ## we first a sufficiently small lower bound
    while (n.lower == FALSE){
      n <- floor(n/2)
      ## new logits are based on linear approximations
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      ## we need to find the optimal gamma value for each sample size
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      ## if n is small enough that average power is not large enough, we stop
      ## otherwise, we need to check a smaller sample size
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to search between n.lower and n.upper to 
    ## get the next sample size n1
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    n1 <- n.lower
  } else {
    ## this clause is initiated if the first sample size n0 is not large enough
    ## and the criterion for average power is not satisfied
    ## we then need to find a larger sample size n.upper where average power is large enough
    n.lower <- n
    n.upper <- FALSE
    while (n.upper == FALSE){
      n <- ceiling(2*n)
      ## new logits are still based on linear approximations
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      ## we stop this while loop once average power and the current sample size is large enough
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to search between n.lower and n.upper to 
    ## get the next sample size n1
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    n1 <- n.upper
  }
  
  ## adjust to make sure the difference between sample sizes is not too small
  if (n1 > n0){
    if (n1/n0 < 1.1){
      n1 <- ceiling(1.1*n0)
    }
  } else {
    if (n1/n0 > 0.9){
      n1 <- floor(0.9*n0)
    }
  }
  
  ## explore the sampling distribution at the second sample size n1
  n <- n1
  samp1 <- foreach(i=1:m, .combine='rbind', .packages = c("qrng"),
                   .options.snow=opts, .errorhandling = "remove") %dopar% {
                     
                     ## generate data
                     set.seed(seed + i + 2*m)
                     y1 <- as.numeric(rmultinom(1, n, par.a[i,]))
                     y2 <- as.numeric(rmultinom(1, round(c*n), par.b[i,]))
                     
                     ## get the posterior parameters for group A
                     hypa <- cbind(alphas1 + head(y1, length(y1) - 1),
                                   betas1 + n - cumsum(head(y1, length(y1) - 1)))
                     
                     ## get the posterior parameters for group B
                     hypb <- cbind(alphas2 + head(y2, length(y2) - 1),
                                   betas2 + round(c*n) - cumsum(head(y2, length(y2) - 1)))
                     
                     # get logits of posterior probabilities and denominators for slopes
                     out.temp <- getLogitsSlopes(hypa, hypb, seed = seed + i + 3*m, 
                                                 n1 = n, n2 = round(c*n),
                                                 outcomes = cat, deltas = deltas)
                     out.temp
                     
                   }
  
  ## adjust for any logits that are +/- Inf
  logits1 <- samp1[, 1:n.out]
  for (j in 1:ncol(logits1)){
    logits1[,j] <- ifelse(logits1[,j] == -Inf, min(subset(logits1, is.finite(logits1[,j]))[,j]) - 1, logits1[,j])
    logits1[,j] <- ifelse(logits1[,j] == Inf, max(subset(logits1, is.finite(logits1[,j]))[,j]) + 1, logits1[,j])
  } 
  
  ## get optimal decision threshold for sample size n1
  ## same process as for n0
  gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
  while (gamma_u - gamma_l > 0.0001){
    gamma_m <- 0.5*gamma_u + 0.5*gamma_l
    bin <- logits1 > logit(gamma_m)
    zeros <- rowMeans(bin) == 0
    fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
    
    fdr.temp <- mean(fdrs)
    
    if (fdr.temp <= q){
      gamma_u <- gamma_m
      fdr.1 <- fdr.temp
    }
    else {
      gamma_l <- gamma_m
    }
  }
  gamma1 <- gamma_m
  
  ## get average_power
  bin <- logits1 > logit(gamma_m)
  zeros <- rowMeans(bin) == 0
  pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
  pwr.1 <- mean(pwrs)
  
  ## extract number of submodels
  num.scen <- max(scen.vec)
  
  ## construct linear approximations for each metric
  slopes.emp <- NULL
  ints.emp <- NULL
  for (k in 1:n.out){
    slopes.temp <- NULL
    ints.temp <- NULL
    for (j in 1:num.scen){
      ## extract logits from the submodel and index them with the final
      ## column to later rejoin the approximations across metrics
      ll0.temp <- cbind(logits0[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
      ll1.temp <- cbind(logits1[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
      
      ## sort the logits for each sample size
      ll0_s <- ll0.temp[order(ll0.temp[,k]),k]
      ll1_s <- ll1.temp[order(ll1.temp[,k]),k]
      
      ## get slopes and intercepts of linear approximations
      l_slope <- (ll1_s - ll0_s)/(n1-n0)
      l_int <- ll0_s - l_slope*n0
      
      ## reorder according to second sample size
      l_slope[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_slope 
      l_int[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_int
      
      ## save results for submodel
      slopes.temp <- c(slopes.temp, l_slope)
      ints.temp <- c(ints.temp, l_int)
    }
    # save results as metrics are iterated over
    slopes.emp <- cbind(slopes.emp, slopes.temp)
    ints.emp <- cbind(ints.emp, ints.temp)
  }
  
  ## get a final sample size recommendation n2 based on these slopes
  ## this process is exactly the same as that used for n1, but we now use the
  ## empirically estimated slopes instead of the limiting ones
  if (pwr.1 >= target.pwr){
    n.upper <- n
    # n0 <- n.upper
    n.lower <- FALSE
    while (n.lower == FALSE){
      n <- floor(n/2)
      ## the limiting slopes are no longer being used to explore new sample sizes
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to get final sample size
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
        ## save fdr, average power, and gamma from final sample size
        gamma_final <- gamma_m
        pwr_final <- pwr.temp
        fdr_final <- fdr.temp
      }
    }
    n2 <- n.upper
  } else {
    n.lower <- n
    # n0 <- n.lower
    n.upper <- FALSE
    while (n.upper == FALSE){
      n <- ceiling(2*n)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to get final sample size
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
        ## save fdr, average power, and gamma from final sample size
        gamma_final <- gamma_m
        pwr_final <- pwr.temp
        fdr_final <- fdr.temp
      }
    }
    n2 <- n.upper
  }
  
  ## take fdr, average power, and gamma from sample size n1 if this
  ## was actually the optimal sample size
  if (n1 == n2){
    gamma_final <- gamma1
    pwr_final <- pwr.1
    fdr_final <- fdr.1
  }
  
  ## return different outputs based on whether users wants to construct bootstrap
  ## confidence intervals or confidence plots
  if (!boot & !contour){
    ## if not constructing bootstrap confidence intervals or contour plots,
    ## return a vector with the final sample size recommendation n2, the 
    ## final decision threshold recommendation, the estimated average power 
    ## at the final sample size, the estimated FDR at the final sample size,
    ## and the initial and second sample size considered (n0 and n1) as well 
    ## as the criteria for the FDR and power
    return(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr))
  } else if (boot & contour){
    ## if constructing bootstrap confidence intervals or contour plots, return
    ## the summary vector above and the vector indexing the submodels (scen.vec) 
    ## (used for bootstrapping), the matrices of logits estimated at n0 and n1 
    ## (used for boostrapping), a matrix dictating which hypotheses having positive
    ## limiting slopes (i.e., which are true) for each simulation repetition 
    ## (used for contour plots and bootstrapping), and the vector of intercepts
    ## and slopes that were empirically estimated (used for contour plots) 
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                scen.vec, logits0, logits1, numer0, ints.emp, slopes.emp))
  } else if (contour){
    ## otherwise, return only the additional components used for the contour plots
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                numer0, ints.emp, slopes.emp))
  } else if (boot){
    ## otherwise, return only the additional components used for bootstrapping
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                scen.vec, logits0, logits1, numer0))
  }
}

## this function implements Algorithm 2 similar to the previous function when
## independent binomial models are used to generate and analyze data

algorithm2Bin <- function(par.a, par.b, q, target.pwr, deltas, 
                          hyper, n0, c = 1, seed = 1, 
                          setup.par = FALSE, contour = FALSE, boot = FALSE){
  
  ## the inputs are described as follows
  ## par.a: a matrix of marginal binomial probabilities for group A; each row corresponds to a simulation
  ##        repetition and the parameters from the same submodel should be in consecutive rows
  ## par.b: a matrix of marginal binomial probabilities for group B; each row corresponds to a simulation
  ##        repetition and the parameters from the same submodel should be in consecutive row
  ## q: an upper bound for the false discovery rate
  ## target.pwr: a lower bound for the average power across the true hypotheses
  ## deltas: a 2-column matrix that corresponds to the delta_L and delta_U interval endpoints
  ##         for each hypothesis; the first column corresponds to the lower interval endpoint 
  #          and the second one corresponds to the upper interval one
  ## hyper: a 4 x (# of outcomes) matrix of hyperparameters for the independent beta distributions that  
  ##        correspond to the independent binomial models.This prior consists of independent beta
  ##        distributions for binomial probabilities; the first row is the alpha parameters for these
  ##        distributions in group A, and the second row is the beta parameters in group A. The third
  ##        row is the alpha parameters in group B, and the fourth row is the beta parameters in group B.
  ## n0: the initial sample size
  ## c: the constant for imbalanced sample sizes (1 by default)
  ## seed: numeric seed for reproducibility
  ## setup.par: a binary indicator denoting whether the parallelization should be set up inside the function
  ##            if FALSE (the default), parallelization should be set up before the function call
  ## contour: a binary indicator denoting whether to return the linear approximations to construct contour
  ##          plots (FALSE by default)
  ## boot: a binary indicator denoting whether to return the estimated logits to construct bootstrap
  ##       confidence intervals (FALSE by default)
  
  ## check the inputs for errors
  if (ncol(par.a) <= 1){
    stop("par.a must contain probabilities for at least two categories")
  } else if (mean(is.numeric(par.a)) < 1){
    stop("par.a must be a numeric matrix")
  } else if (mean(par.a < 0) > 0){
    stop("probs in par.a must be positive")
  } else if (ncol(par.b) <= 1){
    stop("par.b must contain probabilities for at least two categories")
  } else if (mean(is.numeric(par.b)) < 1){
    stop("par.b must be a numeric matrix")
  } else if (mean(par.b < 0) > 0){
    stop("probs in par.b must be positive")
  } else if (ncol(par.a) != ncol(par.b) | nrow(par.a) != nrow(par.b)){
    stop("par.a and par.b must have the same dimensions")
  } else if (!is.numeric(q)){
    stop("q must be a number")
  } else if (q <= 0 | q >= 1){
    stop("q must be at between 0 and 1")
  } else if (!is.numeric(target.pwr)){
    stop("target.pwr must be a number")
  } else if (target.pwr <= 0 | target.pwr >= 1){
    stop("target.pwr must be at between 0 and 1")
  } else if (nrow(hyper) != 4){
    stop("hyper must have four rows")
  } else if (mean(is.numeric(hyper)) < 1){
    stop("hyper must be a numeric matrix")
  } else if (sum(hyper <= 0) > 0){
    stop("all elements of hyper must be positive")
  } else if (ncol(hyper) != ncol(par.a)){
    stop("the dimensions of hyper and probs in par.a and par.b do not match")
  } else if (ncol(deltas) != 2){
    stop("deltas must have two columns")
  } else if (mean(is.numeric(deltas)) < 1){
    stop("deltas must be a numeric matrix")
  } else if (!is.numeric(n0)){
    stop("n0 must be a number")
  } else if (n0 < 1){
    stop("n0 must be at least 1")
  } else if (!is.numeric(c)){
    stop("c must be a number")
  } else if (c < 0){
    stop("c must be at least 0")
  } else if (!is.numeric(seed)){
    stop("seed must be a number")
  } else if (seed <= 0 | (seed%%1 != 0)){
    stop("seed must be a positive integer") 
  } else if (!(setup.par %in% c(TRUE,FALSE))){
    stop("setup.par must be TRUE or FALSE")
  } else if (!(contour %in% c(TRUE,FALSE))){
    stop("contour must be TRUE or FALSE")
  } else if (!(boot %in% c(TRUE,FALSE))){
    stop("boot must be TRUE or FALSE")
  }
  
  ## define functions used to approximate logits of posterior probabilities
  getLogitsSlopes <- function(post1, post2, seed, n1 = n1, n2 = n2,
                              deltas = NULL, mm = 10000) {
    
    ## the inputs are described below
    ## post1: the 2 x (# of outcomes) matrix of posterior parameters for group A returned by PostParams()
    ## post2: the 2 x (# of outcomes) matrix of posterior parameters for group B returned by PostParams()
    ## seed: a numeric seed for reproducibility
    ## n1: sample size for group A
    ## n2: sample size for group B
    ## deltas: a 2-column matrix that corresponds to the delta_L and delta_U interval endpoints
    ##         for each hypothesis on the *logarithmic* scale; the first column corresponds to 
    ##         the lower interval endpoint and the second one corresponds to the upper interval one
    ## mm: the number of simulation repetitions used to estimate the posterior probabilities
    
    ## extract the number of binary outcomes for each group (these should be the same)
    K1 <- nrow(post1) 
    K2 <- nrow(post2) 
    
    set.seed(seed)
    
    ## generate the p variables for group A (unconditional probabilities)
    p1 <- matrix(rbeta(mm*K1, post1[,1], post1[,2]), nrow = K1, ncol = mm)
    
    ## generate the p variables for group B
    p2 <- matrix(rbeta(mm*K2, post2[,1], post2[,2]), nrow = K2, ncol = mm)
    
    ## approximate posteriors for the logarithm of lift
    out <- matrix(0, nrow = K1, ncol = mm)
    denoms <- NULL
    for (i in 1:K1){
      ## get the posterior for log of binomial probabilities
      temp2 <- log(pmax(.Machine$double.eps, p2[i,]))
      temp1 <- log(pmax(.Machine$double.eps, p1[i,]))
      ## calculate the limiting variance related to the Fisher information
      denoms[i] <- var(na.omit(temp1))*n1 + var(na.omit(temp2))*n1
      ## save the posterior draws for the outcome in the for loop
      out[i,] <- temp2 - temp1
    }
    
    if (n1 < 500 | n2 < 500){
      ## create kernel density estimates to compute logits to account for sample-sample
      ## deviations from normality
      logits <- NULL
      for (i in 1:K1){
        ## get kernel density estimate for each outcome
        kd <- density(na.omit(out[i,]))
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          ## calculate the complementary posterior probability based on which interval endpoints are finite
          ## we calculated the complementary probabilities because rounding up to 1 is much more prevalent than
          ## underflow
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE)) +
            mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                            max(na.omit(out[i,])) + 1), kd$bw))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE))
        }
        else{
          logits[i] <- mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw))
        }
      }
    } else {
      ## using a limiting normal approximation that is more stable for large sample sizes
      ## with very large logits
      logits <- NULL
      for (i in 1:K1){
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE) +
            pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE)
        }
        else{
          logits[i] <- pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
      }
    }
    
    ## calculate the logits based on the complementary probabilities calculated above
    logits <- log(1 - logits) - log(logits)
    
    ## return the logits and denominators for the limiting slopes that are 
    ## to be used later
    return(c(logits, denoms))
  }
  
  ## this function is the same as getLogitsSlopes() except the additional information
  ## to estimate the limiting slopes is not returned; this function is used by our
  ## implementation of Algorithm 2 during the second time we explore the sampling 
  ## distribution of posterior probabilities
  getLogits <- function(post1, post2, seed, n1 = n1, n2 = n2,
                        deltas = NULL, mm = 10000) {
    
    K1 <- nrow(post1) 
    K2 <- nrow(post2) 
    
    set.seed(seed)
    
    p1 <- matrix(rbeta(mm*K1, post1[,1], post1[,2]), nrow = K1, ncol = mm)
    p2 <- matrix(rbeta(mm*K2, post2[,1], post2[,2]), nrow = K2, ncol = mm)
    
    out <- matrix(0, nrow = K1, ncol = mm)
    denoms <- NULL
    for (i in 1:K1){
      temp2 <- log(pmax(.Machine$double.eps, p2[i,]))
      temp1 <- log(pmax(.Machine$double.eps, p1[i,]))
      out[i,] <- temp2 - temp1
    }
    
    if (n1 < 500 | n2 < 500){
      logits <- NULL
      for (i in 1:K1){
        kd <- density(na.omit(out[i,]))
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE)) +
            mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                            max(na.omit(out[i,])) + 1), kd$bw))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- mean(pnorm(deltas[i, 2], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw, lower.tail = FALSE))
        }
        else{
          logits[i] <- mean(pnorm(deltas[i, 1], ifelse(is.finite(out[i,]), out[i,],
                                                       max(na.omit(out[i,])) + 1), kd$bw))
        }
      }
    } else {
      logits <- NULL
      for (i in 1:K1){
        if (is.finite(deltas[i, 1]) & is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE) +
            pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
        else if (is.finite(deltas[i, 2])){
          logits[i] <- pnorm(deltas[i, 2], mean(na.omit(out[i,])), sd(na.omit(out[i,])), lower.tail = FALSE)
        }
        else{
          logits[i] <- pnorm(deltas[i, 1], mean(na.omit(out[i,])), sd(na.omit(out[i,])))
        }
      }
    }
    
    logits <- log(1 - logits) - log(logits)
    
    return(logits)
  }
  
  ## define logit function
  logit <- function(x){log(x)-log(1-x)}
  
  ## extract number of categories, outcomes, and simulation repetitions
  n.out <- ncol(par.a)
  m <- nrow(par.a)
  
  ## extract hyperparameters from matrix
  alphas1 <- hyper[1,]; betas1 <- hyper[2,]
  alphas2 <- hyper[3,]; betas2 <- hyper[4,]
  
  ## set up parallelization if setup.par == TRUE
  if (setup.par == TRUE){
    cores=detectCores()
    cl <- makeSOCKcluster(cores[1]-1)
    
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  
  ## identify the rows corresponding to each submodel
  ## determine which consecutive rows are not the same in par.a and par.b
  diff.a <- rowSums(apply(par.a, 2, function(x){abs(diff(x))}))
  diff.b <- rowSums(apply(par.b, 2, function(x){abs(diff(x))}))
  
  ## find the rows where either parameters change and create vector 
  ## that indexes the scenario number
  new.scen <- diff.a + diff.b > 0
  scen.vec <- as.numeric(c(1, cumsum(new.scen) + 1))
  
  ## explore the sampling distribution at the first sample size n0
  n <- n0
  samp0 <- foreach(i=1:m, .combine='rbind',
                   .options.snow=opts, .errorhandling = "remove") %dopar% {
                     ## get the numerators for the limiting slopes
                     eff <- NULL
                     numer <- NULL
                     for (j in 1:n.out){
                       eff[j] <- log(par.b[i,j]) - log(par.a[i,j])
                       num.temp <- min((eff[j] - deltas[j,1])^2, (eff[j] - deltas[j,2])^2)
                       numer[j] <- 0.5*(eff[j] > deltas[j,1] & eff[j] < deltas[j,2])*num.temp
                     }
                     
                     ## generate data
                     set.seed(seed + i)
                     y1 <- rbinom(n.out, n, par.a[i,])
                     y2 <- rbinom(n.out, round(c*n), par.b[i,])
                     
                     ## get posterior parameters for the two groups
                     hypa <- cbind(alphas1 + y1, betas1 + n - y1)
                     hypb <- cbind(alphas2 + y2, betas2 + round(c*n) - y2)
                     
                     # get logits of posterior probabilities and denominators for slopes
                     out.temp <- getLogitsSlopes(hypa, hypb, seed = seed + i + m, 
                                                 n1 = n, n2 = round(c*n),
                                                 deltas = deltas)
                     c(out.temp, numer)
                     
                   }
  
  ## correct any logits that might be +/- infinity (rare)
  logits0 <- samp0[, 1:n.out]
  for (j in 1:ncol(logits0)){
    logits0[,j] <- ifelse(logits0[,j] == -Inf, min(subset(logits0, is.finite(logits0[,j]))[,j]) - 1, logits0[,j])
    logits0[,j] <- ifelse(logits0[,j] == Inf, max(subset(logits0, is.finite(logits0[,j]))[,j]) + 1, logits0[,j])
  } 
  
  ## get limiting slopes to choose next sample size
  denom0 <- samp0[, (n.out+1):(2*n.out)]
  numer0 <- samp0[,(2*n.out+1):(3*n.out)]
  slopes0 <- numer0/denom0
  
  ## get optimal decision threshold for this sample size
  gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
  while (gamma_u - gamma_l > 0.0001){
    ## check if current gamma value controls the FDR
    gamma_m <- 0.5*gamma_u + 0.5*gamma_l
    ## check discoveries
    bin <- logits0 > logit(gamma_m)
    ## this avoids dividing by 0
    zeros <- rowMeans(bin) == 0
    ## only count false discoveries in the numerator; all discoveries in the denominator
    fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
    fdr.temp <- mean(fdrs)
    
    ## update gamma value based on whether FDR criterion is satisified
    if (fdr.temp <= q){
      gamma_u <- gamma_m
    }
    else {
      gamma_l <- gamma_m
    }
  }
  gamma0 <- gamma_m
  
  ## get average for first sample size and gamma combination
  bin <- logits0 > logit(gamma_m)
  zeros <- rowMeans(bin) == 0
  ## only count true discoveries in the numerator
  pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
  pwr.0 <- mean(pwrs)
  
  ## find next sample size based on the limiting slopes
  ## we find a smaller sample size if average power is already large enough
  if (pwr.0 >= target.pwr){
    n.upper <- n
    n.lower <- FALSE
    ## we first a sufficiently small lower bound
    while (n.lower == FALSE){
      n <- floor(n/2)
      ## new logits are based on linear approximations
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      ## we need to find the optimal gamma value for each sample size
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      ## if n is small enough that average power is not large enough, we stop
      ## otherwise, we need to check a smaller sample size
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to search between n.lower and n.upper to 
    ## get the next sample size n1
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    n1 <- n.lower
  } else {
    ## this clause is initiated if the first sample size n0 is not large enough
    ## and the criterion for average power is not satisfied
    ## we then need to find a larger sample size n.upper where average power is large enough
    n.lower <- n
    n.upper <- FALSE
    while (n.upper == FALSE){
      n <- ceiling(2*n)
      ## new logits are still based on linear approximations
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      ## we stop this while loop once average power and the current sample size is large enough
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to search between n.lower and n.upper to 
    ## get the next sample size n1
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- logits0 + slopes0*(n - n0)
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    n1 <- n.upper
  }
  
  ## adjust to make sure the difference between sample sizes is not too small
  if (n1 > n0){
    if (n1/n0 < 1.1){
      n1 <- ceiling(1.1*n0)
    }
  } else {
    if (n1/n0 > 0.9){
      n1 <- floor(0.9*n0)
    }
  }
  
  ## explore the sampling distribution at the second sample size n1
  n <- n1
  samp1 <- foreach(i=1:m, .combine='rbind', .packages = c("qrng"),
                   .options.snow=opts, .errorhandling = "remove") %dopar% {
                     
                     ## generate data
                     set.seed(seed + i + 2*m)
                     y1 <- rbinom(n.out, n, par.a[i,])
                     y2 <- rbinom(n.out, round(c*n), par.b[i,])
                     
                     ## get posterior parameters for the two groups
                     hypa <- cbind(alphas1 + y1, betas1 + n - y1)
                     hypb <- cbind(alphas2 + y2, betas2 + round(c*n) - y2)
                     
                     # get logits of posterior probabilities and denominators for slopes
                     out.temp <- getLogitsSlopes(hypa, hypb, seed = seed + i + 3*m, 
                                                 n1 = n, n2 = round(c*n),
                                                 deltas = deltas)
                     out.temp
                     
                   }
  
  ## adjust for any logits that are +/- Inf
  logits1 <- samp1[, 1:n.out]
  for (j in 1:ncol(logits1)){
    logits1[,j] <- ifelse(logits1[,j] == -Inf, min(subset(logits1, is.finite(logits1[,j]))[,j]) - 1, logits1[,j])
    logits1[,j] <- ifelse(logits1[,j] == Inf, max(subset(logits1, is.finite(logits1[,j]))[,j]) + 1, logits1[,j])
  } 
  
  ## get optimal decision threshold for sample size n1
  ## same process as for n0
  gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
  while (gamma_u - gamma_l > 0.0001){
    gamma_m <- 0.5*gamma_u + 0.5*gamma_l
    bin <- logits1 > logit(gamma_m)
    zeros <- rowMeans(bin) == 0
    fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
    
    fdr.temp <- mean(fdrs)
    
    if (fdr.temp <= q){
      gamma_u <- gamma_m
      fdr.1 <- fdr.temp
    }
    else {
      gamma_l <- gamma_m
    }
  }
  gamma1 <- gamma_m
  
  ## get average_power
  bin <- logits1 > logit(gamma_m)
  zeros <- rowMeans(bin) == 0
  pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
  pwr.1 <- mean(pwrs)
  
  ## extract number of submodels
  num.scen <- max(scen.vec)
  
  ## construct linear approximations for each metric
  slopes.emp <- NULL
  ints.emp <- NULL
  for (k in 1:n.out){
    slopes.temp <- NULL
    ints.temp <- NULL
    for (j in 1:num.scen){
      ## extract logits from the submodel and index them with the final
      ## column to later rejoin the approximations across metrics
      ll0.temp <- cbind(logits0[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
      ll1.temp <- cbind(logits1[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
      
      ## sort the logits for each sample size
      ll0_s <- ll0.temp[order(ll0.temp[,k]),k]
      ll1_s <- ll1.temp[order(ll1.temp[,k]),k]
      
      ## get slopes and intercepts of linear approximations
      l_slope <- (ll1_s - ll0_s)/(n1-n0)
      l_int <- ll0_s - l_slope*n0
      
      ## reorder according to second sample size
      l_slope[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_slope 
      l_int[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_int
      
      ## save results for submodel
      slopes.temp <- c(slopes.temp, l_slope)
      ints.temp <- c(ints.temp, l_int)
    }
    # save results as metrics are iterated over
    slopes.emp <- cbind(slopes.emp, slopes.temp)
    ints.emp <- cbind(ints.emp, ints.temp)
  }
  
  ## get a final sample size recommendation n2 based on these slopes
  ## this process is exactly the same as that used for n1, but we now use the
  ## empirically estimated slopes instead of the limiting ones
  if (pwr.1 >= target.pwr){
    n.upper <- n
    # n0 <- n.upper
    n.lower <- FALSE
    while (n.lower == FALSE){
      n <- floor(n/2)
      ## the limiting slopes are no longer being used to explore new sample sizes
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to get final sample size
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
        ## save fdr, average power, and gamma from final sample size
        gamma_final <- gamma_m
        pwr_final <- pwr.temp
        fdr_final <- fdr.temp
      }
    }
    n2 <- n.upper
  } else {
    n.lower <- n
    # n0 <- n.lower
    n.upper <- FALSE
    while (n.upper == FALSE){
      n <- ceiling(2*n)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
      }
    }
    ## now use binary search to get final sample size
    while (n.upper - n.lower > 1){
      n <- floor(0.5*n.upper + 0.5*n.lower)
      logits.temp <- ints.emp + slopes.emp*n
      gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
      while (gamma_u - gamma_l > 0.0001){
        gamma_m <- 0.5*gamma_u + 0.5*gamma_l
        bin <- logits.temp > logit(gamma_m)
        zeros <- rowMeans(bin) == 0
        fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
        
        fdr.temp <- mean(fdrs)
        
        if (fdr.temp <= q){
          gamma_u <- gamma_m
        }
        else {
          gamma_l <- gamma_m
        }
      }
      
      ## get average_power
      bin <- logits.temp > logit(gamma_m)
      zeros <- rowMeans(bin) == 0
      pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
      pwr.temp <- mean(pwrs)
      if (pwr.temp < target.pwr){
        n.lower <- n
      } else {
        n.upper <- n
        ## save fdr, average power, and gamma from final sample size
        gamma_final <- gamma_m
        pwr_final <- pwr.temp
        fdr_final <- fdr.temp
      }
    }
    n2 <- n.upper
  }
  
  ## take fdr, average power, and gamma from sample size n1 if this
  ## was actually the optimal sample size
  if (n1 == n2){
    gamma_final <- gamma1
    pwr_final <- pwr.1
    fdr_final <- fdr.1
  }
  
  ## return different outputs based on whether users wants to construct bootstrap
  ## confidence intervals or confidence plots
  if (!boot & !contour){
    ## if not constructing bootstrap confidence intervals or contour plots,
    ## return a vector with the final sample size recommendation n2, the 
    ## final decision threshold recommendation, the estimated average power 
    ## at the final sample size, the estimated FDR at the final sample size,
    ## and the initial and second sample size considered (n0 and n1) as well 
    ## as the criteria for the FDR and power
    return(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr))
  } else if (boot & contour){
    ## if constructing bootstrap confidence intervals or contour plots, return
    ## the summary vector above and the vector indexing the submodels (scen.vec) 
    ## (used for bootstrapping), the matrices of logits estimated at n0 and n1 
    ## (used for boostrapping), a matrix dictating which hypotheses having positive
    ## limiting slopes (i.e., which are true) for each simulation repetition 
    ## (used for contour plots and bootstrapping), and the vector of intercepts
    ## and slopes that were empirically estimated (used for contour plots) 
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                scen.vec, logits0, logits1, numer0, ints.emp, slopes.emp))
  } else if (contour){
    ## otherwise, return only the additional components used for the contour plots
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                numer0, ints.emp, slopes.emp))
  } else if (boot){
    ## otherwise, return only the additional components used for bootstrapping
    return(list(c(n2, gamma_final, pwr_final, fdr_final, n0, n1, q, target.pwr),
                scen.vec, logits0, logits1, numer0))
  }
}

## this function takes the output from algorithm2Multi() or algorithm2Bin()
## and creates bootstrap confidence intervals for the sample size recommendation
## and optimal decision threshold. The functions listed above must be run
## with "boot" = TRUE to return the computing results required to implement
## the bootstrap

bootCIs <- function(alg2, M, seed = 2025, setup.par = FALSE){
  
  ## the inputs are described as follows
  ## alg2: the output from algorithm 2 (saved to an object)
  ## M: the number of bootstrap resamples
  ## seed: numeric seed for reproducibility
  ## setup.par: a binary indicator denoting whether the parallelization should be set up inside the function
  ##            if FALSE (the default), parallelization should be set up before the function call
  
  if (!is.numeric(M)){
    stop("M must be a number")
  } else if (M < 1){
    stop("M must be at least 1")
  } else if (!is.numeric(seed)){
    stop("seed must be a number")
  } else if (seed <= 0 | (seed%%1 != 0)){
    stop("seed must be a positive integer") 
  } else if (!(setup.par %in% c(TRUE,FALSE))){
    stop("setup.par must be TRUE or FALSE")
  }
  
  ## set up parallelization if setup.par == TRUE
  if (setup.par == TRUE){
    cores=detectCores()
    cl <- makeSOCKcluster(cores[1]-1)
    
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = M, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  
  ## extract information from algorithm 2 output
  n0 <- alg2[[1]][5]
  n1 <- alg2[[1]][6]
  q <- alg2[[1]][7]
  target.pwr <- alg2[[1]][8]
  scen.vec <- alg2[[2]]
  logits0.orig <- alg2[[3]]
  logits1.orig <- alg2[[4]]
  numer0 <- alg2[[5]]
  ## extract number of submodels
  num.scen <- max(scen.vec)
  n.out <- ncol(logits0.orig)
  
  logit <- function(x){log(x) - log(1-x)}
  
  boot.summary <- foreach(i=1:M, .combine='rbind',
                          .options.snow=opts, .errorhandling = "remove") %dopar% {
                     
                     ## set seed
                     set.seed(seed + 1000*i)
                     logits0 <- NULL
                     logits1 <- NULL
                     for (j in 1:num.scen){
                       ## resample from each submodel separately
                       submod.size <- sum(scen.vec == j)
                       l0.temp <- logits0.orig[which(scen.vec == j),]
                       l1.temp <- logits1.orig[which(scen.vec == j),]
                       logits0 <- rbind(logits0, l0.temp[sample(1:submod.size, submod.size, replace = TRUE),])
                       logits1 <- rbind(logits1, l1.temp[sample(1:submod.size, submod.size, replace = TRUE),])
                     }
                     
                     ## construct linear approximations for each metric
                     slopes.emp <- NULL
                     ints.emp <- NULL
                     for (k in 1:n.out){
                       slopes.temp <- NULL
                       ints.temp <- NULL
                       for (j in 1:num.scen){
                         ## extract logits from the submodel and index them with the final
                         ## column to later rejoin the approximations across metrics
                         ll0.temp <- cbind(logits0[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
                         ll1.temp <- cbind(logits1[which(scen.vec == j),], seq(1, sum(scen.vec == j), 1))
                         
                         ## sort the logits for each sample size
                         ll0_s <- ll0.temp[order(ll0.temp[,k]),k]
                         ll1_s <- ll1.temp[order(ll1.temp[,k]),k]
                         
                         ## get slopes and intercepts of linear approximations
                         l_slope <- (ll1_s - ll0_s)/(n1-n0)
                         l_int <- ll0_s - l_slope*n0
                         
                         ## reorder according to second sample size
                         l_slope[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_slope 
                         l_int[ll1.temp[order(ll1.temp[,k]), n.out + 1]] <- l_int
                         
                         ## save results for submodel
                         slopes.temp <- c(slopes.temp, l_slope)
                         ints.temp <- c(ints.temp, l_int)
                       }
                       # save results as metrics are iterated over
                       slopes.emp <- cbind(slopes.emp, slopes.temp)
                       ints.emp <- cbind(ints.emp, ints.temp)
                     }
                     
                     rbind(slopes.emp, ints.emp)
                     
                     # get optimal decision threshold for sample size n1 based on resampled logits
                     gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
                     while (gamma_u - gamma_l > 0.0001){
                       gamma_m <- 0.5*gamma_u + 0.5*gamma_l
                       bin <- logits1 > logit(gamma_m)
                       zeros <- rowMeans(bin) == 0
                       fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))

                       fdr.temp <- mean(fdrs)
                       
                       if (fdr.temp <= q){
                         gamma_u <- gamma_m
                         fdr.1 <- fdr.temp
                       }
                       else {
                         gamma_l <- gamma_m
                       }
                     }
                     gamma1 <- gamma_m

                     ## get average_power
                     bin <- logits1 > logit(gamma_m)
                     zeros <- rowMeans(bin) == 0
                     pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
                     pwr.1 <- mean(pwrs)
                     n <- n1
                     
                     ## get a final sample size recommendation n2 based on these slopes
                     ## this process is exactly the same as that used for n1, but we now use the
                     ## empirically estimated slopes instead of the limiting ones
                     if (pwr.1 >= target.pwr){
                       n.upper <- n
                       # n0 <- n.upper
                       n.lower <- FALSE
                       while (n.lower == FALSE){
                         n <- floor(n/2)
                         ## the limiting slopes are no longer being used to explore new sample sizes
                         logits.temp <- ints.emp + slopes.emp*n
                         gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
                         while (gamma_u - gamma_l > 0.0001){
                           gamma_m <- 0.5*gamma_u + 0.5*gamma_l
                           bin <- logits.temp > logit(gamma_m)
                           zeros <- rowMeans(bin) == 0
                           fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))

                           fdr.temp <- mean(fdrs)

                           if (fdr.temp <= q){
                             gamma_u <- gamma_m
                           }
                           else {
                             gamma_l <- gamma_m
                           }
                         }

                         ## get average_power
                         bin <- logits.temp > logit(gamma_m)
                         zeros <- rowMeans(bin) == 0
                         pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
                         pwr.temp <- mean(pwrs)
                         if (pwr.temp < target.pwr){
                           n.lower <- n
                         } else {
                           n.upper <- n
                         }
                       }
                       ## now use binary search to get final sample size
                       while (n.upper - n.lower > 1){
                         n <- floor(0.5*n.upper + 0.5*n.lower)
                         logits.temp <- ints.emp + slopes.emp*n
                         gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
                         while (gamma_u - gamma_l > 0.0001){
                           gamma_m <- 0.5*gamma_u + 0.5*gamma_l
                           bin <- logits.temp > logit(gamma_m)
                           zeros <- rowMeans(bin) == 0
                           fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))

                           fdr.temp <- mean(fdrs)

                           if (fdr.temp <= q){
                             gamma_u <- gamma_m
                           }
                           else {
                             gamma_l <- gamma_m
                           }
                         }

                         ## get average_power
                         bin <- logits.temp > logit(gamma_m)
                         zeros <- rowMeans(bin) == 0
                         pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
                         pwr.temp <- mean(pwrs)
                         if (pwr.temp < target.pwr){
                           n.lower <- n
                         } else {
                           n.upper <- n
                           ## save fdr, average power, and gamma from final sample size
                           gamma_final <- gamma_m
                           pwr_final <- pwr.temp
                           fdr_final <- fdr.temp
                         }
                       }
                       n2 <- n.upper
                     } else {
                       n.lower <- n
                       # n0 <- n.lower
                       n.upper <- FALSE
                       while (n.upper == FALSE){
                         n <- ceiling(2*n)
                         logits.temp <- ints.emp + slopes.emp*n
                         gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
                         while (gamma_u - gamma_l > 0.0001){
                           gamma_m <- 0.5*gamma_u + 0.5*gamma_l
                           bin <- logits.temp > logit(gamma_m)
                           zeros <- rowMeans(bin) == 0
                           fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
                           fdr.temp <- mean(fdrs)

                           if (fdr.temp <= q){
                             gamma_u <- gamma_m
                           }
                           else {
                             gamma_l <- gamma_m
                           }
                         }

                         ## get average_power
                         bin <- logits.temp > logit(gamma_m)
                         zeros <- rowMeans(bin) == 0
                         pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
                         pwr.temp <- mean(pwrs)
                         if (pwr.temp < target.pwr){
                           n.lower <- n
                         } else {
                           n.upper <- n
                         }
                       }
                       ## now use binary search to get final sample size
                       while (n.upper - n.lower > 1){
                         n <- floor(0.5*n.upper + 0.5*n.lower)
                         logits.temp <- ints.emp + slopes.emp*n
                         gamma_u <- 1 - 10^(-7); gamma_l <- 0.5
                         while (gamma_u - gamma_l > 0.0001){
                           gamma_m <- 0.5*gamma_u + 0.5*gamma_l
                           bin <- logits.temp > logit(gamma_m)
                           zeros <- rowMeans(bin) == 0
                           fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))

                           fdr.temp <- mean(fdrs)

                           if (fdr.temp <= q){
                             gamma_u <- gamma_m
                           }
                           else {
                             gamma_l <- gamma_m
                           }
                         }

                         ## get average_power
                         bin <- logits.temp > logit(gamma_m)
                         zeros <- rowMeans(bin) == 0
                         pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
                         pwr.temp <- mean(pwrs)
                         if (pwr.temp < target.pwr){
                           n.lower <- n
                         } else {
                           n.upper <- n
                           ## save fdr, average power, and gamma from final sample size
                           gamma_final <- gamma_m
                           pwr_final <- pwr.temp
                           fdr_final <- fdr.temp
                         }
                       }
                       n2 <- n.upper
                     }

                     ## take fdr, average power, and gamma from sample size n1 if this
                     ## was actually the optimal sample size
                     if (n1 == n2){
                       gamma_final <- gamma1
                       pwr_final <- pwr.1
                       fdr_final <- fdr.1
                     }

                     ## return outputs for each bootstrap repetition
                     c(n2, gamma_final, pwr_final, fdr_final)
                   }

  return(boot.summary)
  
}