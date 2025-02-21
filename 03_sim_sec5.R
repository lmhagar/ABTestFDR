## this code file reproduces the numerical results in Section 5 of the paper
## please run the previous code files in this repository to ensure the 
## necessary functions and parameter vectors are loaded

require(foreach)
require(doParallel)
require(doSNOW)

## implement sample size calculation where data are modelled with a 
## multinomial distribution and there is uncertainty about the number
## of false hypotheses

## set up parallelization
m.ssd <- 30000
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement sample size calculation
tic <- Sys.time()
ssd.multi <- algorithm2Multi(par.a = params.full[rep(seq(1, 59, 2), each = m.ssd/30), ], 
                             par.b = params.full[rep(seq(2, 60, 2), each = m.ssd/30), ], 
                             q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                             mat = tree.mat, 
                             hyper = rbind(rep(1,12), seq(12, 1, by = -1), rep(1,12), seq(12, 1, by = -1)), 
                             n0 = 12000, c = 1, seed = 1, 
                             setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic
## save results
for (i in 1:length(ssd.multi)){
  write.csv(ssd.multi[[i]], paste0("ssd_multi", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.multi <- bootCIs(ssd.multi, M = M.boot, seed = 100, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.multi, "boot_multi.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.multi[,1], probs = c(0.025, 0.975))
quantile(boot.multi[,2], probs = c(0.025, 0.975))

## confirm the performance of the recommended sample size and decision threshold
## create a function to estimate posterior probabilities for the multinomial model
getProbsMulti <- function(post1, post2, seed, n1 = n1, n2 = n2,
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
  for (i in 1:length(outcomes)){
    inds <- outcomes[[i]]
    temp2 <- log(pmax(.Machine$double.eps, colSums(p2[inds,])))
    temp1 <- log(pmax(.Machine$double.eps, colSums(p1[inds,])))
    out[i,] <- temp2 - temp1
  }
  
  probs <- NULL
  for (i in 1:length(outcomes)){
    probs[i] <- mean(out[i,] > deltas[i, 1] & out[i,] < deltas[i,2])
  }
  
  return(probs)
}

m.con <- 99000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m.con, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## extract sample size recommendation
n.multi <- ssd.multi[[1]][1]
## define hyperparameters for multinomial distribution
alphas1 <- rep(1,12); betas1 <- seq(12, 1, by = -1)
alphas2 <- rep(1,12); betas2 <- seq(12, 1, by = -1)
## extract list of outcomes from matrix in createTree()
cats <- apply(tree.mat, 2, function(x){which(x == 1)})
tic <- Sys.time()
con.multi <- foreach(i=1:m.con, .combine='rbind',
                  .options.snow=opts, .errorhandling = "remove") %dopar% {
                  n <- n.multi
                  
                  ## generate multinomial data
                  set.seed(200 + i)
                  y1 <- as.numeric(rmultinom(1, n, params.full[(2*i-1)%%60,]))
                  y2 <- as.numeric(rmultinom(1, n, params.full[ifelse((2*i)%%60==0, 60, (2*i)%%60),]))
                  
                  ## get posterior parameters
                  alphas.1 <- alphas1 + head(y1, length(y1) - 1)
                  betas.1 <- betas1 + n - cumsum(head(y1, length(y1) - 1))
                  
                  alphas.2 <- alphas2 + head(y2, length(y2) - 1)
                  betas.2 <- betas2 + n - cumsum(head(y2, length(y2) - 1))
                  
                  ## get posterior probabilities
                  list.temp <- getProbsMulti(cbind(alphas.1, betas.1), cbind(alphas.2, betas.2),
                                            seed = i + m.con, 
                                            outcomes = cats, 
                                            deltas = cbind(rep(0, 5), rep(Inf, 5)), mm = 10000)
                  
                  as.numeric(list.temp)
                }
toc <- Sys.time()
toc - tic

## save results
write.csv(con.multi, "con_multi.csv", row.names = FALSE)

## confirm FDR and average power
## check which posterior probabilities are greater than the recommended threshold
bin <- con.multi > ssd.multi[[1]][2]
## determine which repetitions have no discoveries to avoid division by 0
zeros <- rowMeans(bin) == 0
## create a matrix that dictates which hypotheses are true (1) and false (0)
trues <- NULL
for (j in 1:4){
  comb.temp <- combn(5,j)
  for (k in 1:ncol(comb.temp)){
    trues <- rbind(trues, ifelse(1:5 %in% comb.temp[,k], 0, 1))
  }
}

## repeat for all simulation repetitions
trues.all <- trues[rep(1:30, m.con/30),]

## compute confirmatory power and FDR
fdr.multi <- mean(ifelse(zeros, 0, rowSums((1-trues.all)*bin)/rowSums(bin)))
pwr.multi <- mean(ifelse(zeros, 0, rowSums(trues.all*bin)/rowSums(trues.all)))

## implement sample size calculation where data are modelled with independent
## binomial distributions and there is uncertainty about the number
## of false hypotheses

## set up parallelization
m.ssd <- 30000

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## create functions to transform multinomial probabilities to 
## marginal binomial probabilities for each scenario
groupCat <- function(x, cats){
  temp <- NULL
  for (i in 1:length(cats)){
    temp[i] <- sum(x[cats[[i]]])
  }
  temp
}

## get marginal parameters
mars.full <- t(apply(params.full, 1, groupCat, cats = cats))

## implement sample size calculation
tic <- Sys.time()
ssd.bin <- algorithm2Bin(par.a = mars.full[rep(seq(1, 59, 2), each = m.ssd/30), ], 
                         par.b = mars.full[rep(seq(2, 60, 2), each = m.ssd/30), ], 
                         q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                         hyper = rbind(rep(1,5), rep(1,5), rep(1,5), rep(1,5)), 
                         n0 = 12000, c = 1, seed = 300, 
                         setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic
## save results
for (i in 1:length(ssd.bin)){
  write.csv(ssd.bin[[i]], paste0("ssd_bin", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.bin <- bootCIs(ssd.bin, M = M.boot, seed = 100, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.bin, "boot_bin.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.bin[,1], probs = c(0.025, 0.975))
quantile(boot.bin[,2], probs = c(0.025, 0.975))

## confirm the performance of the recommended sample size and decision threshold
## when the data are generated from independent binomial distributions
m.con <- 99000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m.con, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## extract sample size recommendation
n.bin <- ssd.bin[[1]][1]

## implement confirmation simulations with BETA(1,1) priors for each probability
tic <- Sys.time()
con.bin <- foreach(i=1:m.con, .combine='rbind',
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  n <- n.bin
                  
                  set.seed(i + 900)
                  y1 <- rbinom(5, n, mars.full[(2*i-1)%%60,])
                  y2 <- rbinom(5, n, mars.full[ifelse((2*i)%%60==0, 60, (2*i)%%60),])
                  
                  ## get posterior parameters for the two groups
                  alphas.1 <- 1 + y1; betas.1 <- 1 + n - y1
                  alphas.2 <- 1 + y2; betas.2 <- 1 + n - y2
                  
                  list.temp <- NULL
                  for (j in 1:5){
                    temp1 <- log(rbeta(100000, alphas.1[j], betas.1[j]))
                    temp2 <- log(rbeta(100000, alphas.2[j], betas.2[j]))
                    temp <- temp2 - temp1
                    list.temp[j] <- mean(temp > 0 & temp < Inf)
                  }
                  
                  list.temp
                }
toc <- Sys.time()
toc - tic

## save results
write.csv(con.bin, "con_bin.csv", row.names = FALSE)

## confirm FDR and average power
## check which posterior probabilities are greater than the recommended threshold
bin <- con.bin > ssd.bin[[1]][2]
## determine which repetitions have no discoveries to avoid division by 0
zeros <- rowMeans(bin) == 0
## use same trues.all matrix from before to compute confirmatory power and FDR
fdr.bin <- mean(ifelse(zeros, 0, rowSums((1-trues.all)*bin)/rowSums(bin)))
pwr.bin <- mean(ifelse(zeros, 0, rowSums(trues.all*bin)/rowSums(trues.all)))

## we now implement sample size calculations for the multinomial model
## where we separately consider the submodels where 1, 2, 3, and 4 
## hypotheses are false; this is done to consider the impact of 
## the assumed number of false hypotheses on the sample size calculation

## consider the models with one false hypothesis
m.ssd <- 30000

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement sample size calculation
tic <- Sys.time()
ssd.1 <- algorithm2Multi(par.a = params1[rep(seq(1, 9, 2), each = m.ssd/5), ], 
                         par.b = params1[rep(seq(2, 10, 2), each = m.ssd/5), ], 
                         q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                         mat = tree.mat, 
                         hyper = rbind(rep(1,12), seq(12, 1, by = -1), rep(1,12), seq(12, 1, by = -1)), 
                         n0 = 7000, c = 1, seed = 400, 
                         setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic

## save results
for (i in 1:length(ssd.1)){
  write.csv(ssd.1[[i]], paste0("ssd_1", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.1 <- bootCIs(ssd.1, M = M.boot, seed = 500, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.1, "boot_1.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.1[,1], probs = c(0.025, 0.975))
quantile(boot.1[,2], probs = c(0.025, 0.975))

## consider the models with two false hypotheses
m.ssd <- 30000

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement sample size calculation
tic <- Sys.time()
ssd.2 <- algorithm2Multi(par.a = params2[rep(seq(1, 19, 2), each = m.ssd/10), ], 
                         par.b = params2[rep(seq(2, 20, 2), each = m.ssd/10), ], 
                         q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                         mat = tree.mat, 
                         hyper = rbind(rep(1,12), seq(12, 1, by = -1), rep(1,12), seq(12, 1, by = -1)), 
                         n0 = 11000, c = 1, seed = 600, 
                         setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic

ssd.2[[1]]

## save results
for (i in 1:length(ssd.2)){
  write.csv(ssd.2[[i]], paste0("ssd_2", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.2 <- bootCIs(ssd.2, M = M.boot, seed = 700, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.2, "boot_2.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.2[,1], probs = c(0.025, 0.975))
quantile(boot.2[,2], probs = c(0.025, 0.975))

## consider the models with three false hypotheses
m.ssd <- 30000

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement sample size calculation
tic <- Sys.time()
ssd.3 <- algorithm2Multi(par.a = params3[rep(seq(1, 19, 2), each = m.ssd/10), ], 
                         par.b = params3[rep(seq(2, 20, 2), each = m.ssd/10), ], 
                         q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                         mat = tree.mat, 
                         hyper = rbind(rep(1,12), seq(12, 1, by = -1), rep(1,12), seq(12, 1, by = -1)), 
                         n0 = 15000, c = 1, seed = 800, 
                         setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic

ssd.3[[1]]

## save results
for (i in 1:length(ssd.3)){
  write.csv(ssd.3[[i]], paste0("ssd_3", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.3 <- bootCIs(ssd.3, M = M.boot, seed = 900, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.3, "boot_3.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.3[,1], probs = c(0.025, 0.975))
quantile(boot.3[,2], probs = c(0.025, 0.975))

## consider the models with four false hypotheses
m.ssd <- 30000

registerDoSNOW(cl)
pb <- txtProgressBar(max = m.ssd, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement sample size calculation
tic <- Sys.time()
ssd.4 <- algorithm2Multi(par.a = params4[rep(seq(1, 9, 2), each = m.ssd/5), ], 
                         par.b = params4[rep(seq(2, 10, 2), each = m.ssd/5), ], 
                         q = 0.05, target.pwr = 0.8, deltas = cbind(rep(0, 5), rep(Inf, 5)), 
                         mat = tree.mat, 
                         hyper = rbind(rep(1,12), seq(12, 1, by = -1), rep(1,12), seq(12, 1, by = -1)), 
                         n0 = 18000, c = 1, seed = 1000, 
                         setup.par = FALSE, contour = TRUE, boot = TRUE)
toc <- Sys.time()
toc - tic

ssd.4[[1]]

## save results
for (i in 1:length(ssd.4)){
  write.csv(ssd.4[[i]], paste0("ssd_4", i, ".csv"), row.names = FALSE)
}

## implement bootstrapping process
M.boot <- 5000
registerDoSNOW(cl)
pb <- txtProgressBar(max = M.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic <- Sys.time()
boot.4 <- bootCIs(ssd.4, M = M.boot, seed = 1100, setup.par = FALSE)
toc <- Sys.time()
toc - tic

## save results
write.csv(boot.4, "boot_4.csv", row.names = FALSE)

## get bootstrap confidence intervals
quantile(boot.4[,1], probs = c(0.025, 0.975))
quantile(boot.4[,2], probs = c(0.025, 0.975))

## get recommendations for frequentist sample size calculation with 
## the Bonferroni correction

## get the marginal Bernoulli parameters for each outcome
theta.a <- mars.full[1,]
theta.b <- c(mars.full[4,1], mars.full[2, 2:5])

## this function returns the sample size recommendation for a given
## alpha value (type I error), beta value (type II error), and marginal
## success probabilities for groups A and B
sampleSizeBon <- function(alpha, beta, thetas){
  ceiling((qnorm(1-alpha) - qnorm(beta))^2*(1/thetas[1] + 1/thetas[2] -2)/(log(thetas[2]) - log(thetas[1]))^2)
}

## get sample size recommendations
samps.bon <- as.numeric(apply(rbind(theta.a, theta.b), 2, sampleSizeBon, alpha = 0.01, beta = 0.2))

## return average sample size
mean(samps.bon)

## now get sample size recommendation based on average power
## (i.e., what is the sample size such that the total power of
## all the hypothesis tests exceeds 5*0.8 = 4)

## the inputs are the same for this function as sampleSizeBon except
## we now pass in the sample size instead of the type II error rate
totalPower <- function(alpha, n, thetas){
  total.pwr <- 0
  for (j in 1:ncol(thetas)){
    pwr.temp <- pnorm(qnorm(1 - alpha), mean = sqrt(n)*(log(theta.b[j])-log(theta.a[j]))/sqrt(1/theta.a[j] + 1/theta.b[j] - 2), 
                      sd = 1, lower.tail = FALSE)
    total.pwr <- total.pwr + pwr.temp
  }
  return(total.pwr)
}

## check total power at mean sample size recommendation
totalPower(0.01, 25349, rbind(theta.a, theta.b))

## we now create a function with a "target" for power
## to put into a root finding algorithm and find the
## sample size required
totalPowerRoot <- function(alpha, n, thetas, target){
  total.pwr <- 0
  for (j in 1:ncol(thetas)){
    pwr.temp <- pnorm(qnorm(1 - alpha), mean = sqrt(n)*(log(theta.b[j])-log(theta.a[j]))/sqrt(1/theta.a[j] + 1/theta.b[j] - 2), 
                      sd = 1, lower.tail = FALSE)
    total.pwr <- total.pwr + pwr.temp
  }
  return(total.pwr - target)
}

## find sample size to ensure average power is sufficient
ceiling(uniroot(totalPowerRoot, lower = 25349, upper = 50000, alpha = 0.01, 
        thetas = rbind(theta.a, theta.b), target = 4)$root)