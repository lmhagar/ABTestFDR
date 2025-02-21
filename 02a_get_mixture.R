## this code file uses the functions in "01a_build_tree.R" and "01_lp.R"
## to get all 30 scenarios that define our mixture model Psi

## the createTree function is used to obtain the matrix that is input into the LP functions
tree.opt <- createTree("Optimizely", c("Engage", "Editor", "Price", "Dialog", "Create"))

## save this matrix
tree.mat <- tree.opt[[2]]

## specify the marginal probabilities for group A based on the Optimizely data
marginal.a <- c(0.489, 0.230, 0.156, 0.047, 0.032)

params.full <- NULL
## iterate over the number of false hypotheses
for (j in 1:4){
  ## get all combinations of false hypotheses
  combs <- combn(5, j)
  
  params.temp <- NULL
  ## iterate over each combination
  for (k in 1:ncol(combs)){
    ## get vector of lifts
    lifts.1 <- ifelse(1:5 %in% combs[,k], 0, 0.1)
    
    ## initialize LP
    lp.1 <- initializeLP(tree.mat, marginal.a, 
                           lifts.1, min.prob = 0.0025)
    
    ## check to ensure LP has a feasible solution
    params <- lpSolveFn(lp.1)
    params
    
    ## reduce the probability for looking at the price but not the editor
    ## change for group a
    lp.2 <- addConstraint(lp.1, group = "a", 
                            yes = c(1, 3), no = c(2, 4, 5), dir = "=", rhs = 0.05)
    
    params <- lpSolveFn(lp.2)
    params
    
    ## change for group b
    lp.3 <- addConstraint(lp.2, group = "b", 
                            yes = c(1, 3), no = c(2, 4, 5), dir = "=", rhs = 0.06)
    
    params <- lpSolveFn(lp.3)
    params
    
    ## increase probability of experiencing all five outcomes
    ## change for group a
    lp.4 <- addConstraint(lp.3, group = "a", 
                            yes = c(1, 2, 3, 4, 5), no = NULL, dir = "=", rhs = 0.01)
    
    params <- lpSolveFn(lp.4)
    params
    
    ## change for group b
    lp.5 <- addConstraint(lp.4, group = "b", 
                            yes = c(1, 2, 3, 4, 5), no = NULL, dir = "=", rhs = 0.0125)
    
    params <- lpSolveFn(lp.5)
    params
    
    params.temp <- rbind(params.temp, params)
  }
  assign(paste0("params", j), params.temp)
  params.full <- rbind(params.full, params.temp)
}

## save parameters for full mixture model, where each group of two
## rows corresponds to a scenario. The first row denotes the multinomial
## probabilities for group A and the second denotes the probabilities
## for group B
write.csv(params.full, "params_full.csv", row.names = FALSE)

## save submodels in groups based on the number of false hypotheses;
## for use later
write.csv(params1, "params1.csv", row.names = FALSE)
write.csv(params2, "params2.csv", row.names = FALSE)
write.csv(params3, "params3.csv", row.names = FALSE)
write.csv(params4, "params4.csv", row.names = FALSE)

