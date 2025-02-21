## this code file reproduces the contour plots in Section 6 of the paper
## please run the previous code files in this repository to ensure the 
## necessary functions and parameter vectors are loaded

require(foreach)
require(doParallel)
require(doSNOW)

## this code chunk generates the confirmation samples for the right 
## contour plot

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m.con <- 99000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m.con, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## extract sample size recommendation
ns <- seq(10000, 15000, by = 500)
## define hyperparameters for multinomial distribution
alphas1 <- rep(1,12); betas1 <- seq(12, 1, by = -1)
alphas2 <- rep(1,12); betas2 <- seq(12, 1, by = -1)
## extract list of outcomes from matrix in createTree()
cats <- apply(tree.mat, 2, function(x){which(x == 1)})
tic <- Sys.time()
for (l in 1:length(ns)){
  con.multi.temp <- foreach(i=1:m.con, .combine='rbind',
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         n <- ns[l]
                         
                         ## generate multinomial data
                         set.seed(5000 + i)
                         y1 <- as.numeric(rmultinom(1, n, params.full[(2*i-1)%%60,]))
                         y2 <- as.numeric(rmultinom(1, n, params.full[ifelse((2*i)%%60==0, 60, (2*i)%%60),]))
                         
                         ## get posterior parameters
                         alphas.1 <- alphas1 + head(y1, length(y1) - 1)
                         betas.1 <- betas1 + n - cumsum(head(y1, length(y1) - 1))
                         
                         alphas.2 <- alphas2 + head(y2, length(y2) - 1)
                         betas.2 <- betas2 + n - cumsum(head(y2, length(y2) - 1))
                         
                         ## get posterior probabilities
                         list.temp <- getProbsMulti(cbind(alphas.1, betas.1), cbind(alphas.2, betas.2),
                                                    seed = 10000 + i + m.con, 
                                                    outcomes = cats, 
                                                    deltas = cbind(rep(0, 5), rep(Inf, 5)), mm = 10000)
                         
                         as.numeric(list.temp)
                       }
  toc <- Sys.time()
  toc - tic
  
  ## save results
  write.csv(con.multi.temp, paste0("con_multi_", ns[l] ,".csv"), row.names = FALSE)
}


## we now construct the contour plots for Section 6
## create matrices for left contour plot in Figure 2 (based on single
## sample size calculation)
res <- read.csv("ssd_multi1.csv")$x
n_low <- res[5]; n_high <- res[6]

## get matrix that dictates which hypotheses were true
numer0 <- read.csv("ssd_multi5.csv")

## read in the posterior probabilities corresponding to n0 and n1
ints <- read.csv("ssd_multi6.csv")
slopes <- read.csv("ssd_multi7.csv")

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(10000, 15000, 10)){
  assign(paste0("l_vec_", i), ints + slopes*i)
}

## create a vector of gamma values on the logit scale to compute power
## and type I error rate estimates
opt_gamma <- as.numeric(res[2])
opt_gamma <- log(opt_gamma) - log(1 - opt_gamma)

gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(10000, 15000, 10)
y <- 1/(1 + exp(-gammas))

## z matrix is for average power
## w matrix is for the FDR
z <- matrix(0, nrow = length(x), ncol = length(y))
w <- matrix(0, nrow = length(x), ncol = length(y))
for (i in 1:length(x)){
  print(i)
  for (j in 1:length(y)){
    bin <- get(paste0("l_vec_", x[i])) > gammas[j]
    zeros <- rowMeans(bin) == 0
    pwrs <- ifelse(zeros, 0, rowSums((numer0>0)*bin)/rowSums(numer0>0))
    fdrs <- ifelse(zeros, 0, rowSums((numer0<=0)*bin)/rowSums(bin))
    z[i,j] <- mean(pwrs)
    w[i,j] <- mean(fdrs)
  }
}

write.csv(w, "w_mat1.csv", row.names = FALSE)
write.csv(z, "z_mat1.csv", row.names = FALSE)

## create matrices for right contour plot in Figure 2 (based on
## simulating data and repeatedly estimating the sampling distribution)
z_full2 <- matrix(0, nrow = 11, ncol = 50)
w_full2 <- matrix(0, nrow = 11, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in ns){
  assign(paste0("ll_vec_", i), 
         read.csv(paste0("con_multi_", i ,".csv")))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)

x <- seq(10000, 15000, 500)
y <- 1/(1 + exp(-gammas))

## construct matrix to denote which hypotheses are true
trues <- NULL
for (j in 1:4){
  comb.temp <- combn(5,j)
  for (k in 1:ncol(comb.temp)){
    trues <- rbind(trues, ifelse(1:5 %in% comb.temp[,k], 0, 1))
  }
}

## repeat for all simulation repetitions
trues.all <- trues[rep(1:30, m.con/30),]

## z2 matrix is for average power
## w2 matrix is for the FDR
z2 <- matrix(0, nrow = length(x), ncol = length(y))
w2 <- matrix(0, nrow = length(x), ncol = length(y))
for (i in 1:length(x)){
  print(i)
  for (j in 1:length(y)){
    bin <- get(paste0("ll_vec_", x[i])) > y[j]
    zeros <- rowMeans(bin) == 0
    pwrs <- ifelse(zeros, 0, rowSums(trues.all*bin)/rowSums(trues.all))
    fdrs <- ifelse(zeros, 0, rowSums((1 - trues.all)*bin)/rowSums(bin))
    z2[i,j] <- mean(pwrs)
    w2[i,j] <- mean(fdrs)
  }
}

## write output to a .csv file  
write.csv(z2, "z_mat2.csv", row.names = FALSE)
write.csv(w2, "w_mat2.csv", row.names = FALSE)

## create the contour plots and output as .pdf file for the article
pdf(file = "Figure2.pdf",   # The directory you want to save the file in
    width = 6, 
    height = 6) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.35,1,0))

## read in matrices for left plot
z <- matrix(unlist(read.csv("z_mat1.csv")), nrow = 501, ncol = 51)
w <- matrix(unlist(read.csv("w_mat1.csv")), nrow =501, ncol = 51)
gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(10000, 15000, 10)
y <- 1/(1 + exp(-gammas))

contour(x, y, w, levels = c(seq(0.02, 0.035, 0.015), 0.065, seq(0.095, 0.14, 0.015)), 
        xlab = expression(italic("n")['A']~' (Thousands)'), xlim = c(10000,15000), ylim = c(0.8, 0.991),
        ylab = expression(gamma),  main = "False Discovery Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.08), col = "black", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
points(x = res[1], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
axis(side = 1, 
     at = seq(10000, 15000, 1000), 
     labels = seq(10, 15, 1),
     cex.axis = 1.15)
axis(side = 2, at = c(0.81,0.87, 0.93, 0.99), cex.axis = 1.15)
box()

## reset plotting parameters for plots based on repeated estimates of the 
## sampling distribution
gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)

x <- seq(10000, 15000, 500)
y <- 1/(1 + exp(-gammas))

contour(x, y, w2, levels = c(seq(0.02, 0.035, 0.015), 0.065, seq(0.095, 0.14, 0.015)), 
        xlab = expression(italic("n")['A']~' (Thousands)'), xlim = c(10000,15000),
        ylab = expression(gamma),  main = "False Discovery Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w2, levels = c(0.08), col = "black", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, z2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, 
     at = seq(10000, 15000, 1000), 
     labels = seq(10, 15, 1),
     cex.axis = 1.15)
axis(side = 2, at = c(0.81,0.87, 0.93, 0.99), cex.axis = 1.15)
box()

## reset plotting parameters for plots based on repeated estimates of the 
## sampling distribution
gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(10000, 15000, 10)
y <- 1/(1 + exp(-gammas))

contour(x, y, z, levels = c(seq(0.60, 0.76, 0.04), c(0.83, 0.85, 0.87, 0.89, 0.91)), 
        xlab = expression(italic("n")['A']~' (Thousands)'), xlim = c(10000,15000),
        ylab = expression(gamma),  main = "Average Power", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
points(x = res[1], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
axis(side = 1, 
     at = seq(10000, 15000, 1000), 
     labels = seq(10, 15, 1),
     cex.axis = 1.15)
axis(side = 2, at = c(0.81,0.87, 0.93, 0.99), cex.axis = 1.15)
box()

## reset plotting parameters for plots based on repeated estimates of the 
## sampling distribution
gammas <- seq(log(0.8) - log(0.2), log(0.99) - log(0.01), length.out = 50)

x <- seq(10000, 15000, 500)
y <- 1/(1 + exp(-gammas))

contour(x, y, z2, levels = c(seq(0.60, 0.76, 0.04), c(0.83, 0.85, 0.87, 0.89, 0.91)), 
        xlab = expression(italic("n")['A']~' (Thousands)'), xlim = c(10000,15000),
        ylab = expression(gamma),  main = "Average Power", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, z2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, 
     at = seq(10000, 15000, 1000), 
     labels = seq(10, 15, 1),
     cex.axis = 1.15)
axis(side = 2, at = c(0.81,0.87, 0.93, 0.99), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()