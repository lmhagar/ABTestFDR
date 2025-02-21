## These series of functions are used to obtain the multinomial models
## that define the mixture components of the data generating model Psi
## for the Optimizely example

## load required libraries
require(lpSolve)

## this function initializes the LP so that the following criteria are satisfied
## 1. the multinomial probabilities in each group should sum to 1
## 2. the marginal outcome probabilities in group A are attained
## 3. the multinomial probabilities in group B achieve targets for marginal lift
## 4. all multinomial probabilities in both models are greater than a minimum value

## the output from this function is described before the return statement, and this
## output can be input into the other functions in this file
initializeLP <- function(mat, as, lifts, min.prob){
  
  ## the inputs are described as follows
  ## mat: the matrix that is the second item in the list returned by createTree()
  ## as: a vector of marginal probabilities for group A
  ## lifts: a vector of target lifts expressed as percentage increases 
  ##        (e.g., 0, 0.1, etc.); there must be one lift for each outcome.
  ## min.prob: a lower bound for the probability of each category
  
  if (nrow(mat) <= 2){
    stop("mat must have at least two rows")
  } else if (ncol(mat) <= 2){
    stop("mat must have at least two columns")
  } else if (mean(mat == 0 | mat == 1)!=1){
    stop("mat must be a matrix output by createTree()")
  } else if (length(as) <= 2){
    stop("as must have at least two probabilitites")
  } else if (sum(as <= 0 | as >= 1)){
    stop("the probabilities must all be between 0 and 1")
  } else if (length(as) != ncol(mat)){
    stop("the number of probabilities in as must match the number of columns in mat")
  } else if (length(lifts) <= 2){
    stop("lifts must have at least two ratios")
  } else if (mean(ifelse(lifts > -1, 1, ifelse(lifts < 10))) != 1){
    stop("please choose valid percentage increases for lift")
  } else if (length(lifts) != ncol(mat)){
    stop("the number of lifts must match the number of columns in mat")
  } else if (length(min.prob) != 1){
    stop("please provide one minimum probability")
  } else if (!is.numeric(min.prob)){
    stop("please provide a numeric probability")
  } else if (min.prob <= 0 | min.prob > 1/length(as)){
    stop("please provide a valid probability")
  }
  
  ## start with creating an objective function of all 0s
  n.cat <- nrow(mat)
  objective <- rep(0, 2*n.cat)
  
  ## create constraints to ensure all multinomial probabilities sum to 1
  const.mat <- matrix(rep(c(1,0,0,1), each = n.cat), 
                        nrow = 2, byrow = TRUE)
  
  ## create constraints to enforce minimum probability in group A
  const.mat <- rbind(const.mat,
                     cbind(diag(n.cat), matrix(0, nrow = n.cat, ncol = n.cat)))
  
  ## create constraints to enforce minimum probability in group B
  const.mat <- rbind(const.mat,
                     cbind(matrix(0, nrow = n.cat, ncol = n.cat), diag(n.cat)))
  
  ## create constraints to attain target marginal probabilities in group A
  inds <- apply(mat, 2, function(x){which(x == 1)})
  
  for (j in 1:ncol(mat)){
    const.mat <- rbind(const.mat,
                       c(rep(0, n.cat), ifelse(1:n.cat %in% inds[[j]], 1, 0)))
  }
  
  ## create constraints to achieve target lifts 
  for (j in 1:ncol(mat)){
    b.temp <- ifelse(1:n.cat %in% inds[[j]], -1, 0)
    a.temp <- ifelse(1:n.cat %in% inds[[j]], 1 + lifts[j], 0)
    const.mat <- rbind(const.mat,
                       c(b.temp, a.temp))
  }
  
  ## get direction vector for constraints (same order as rows in const.mat above)
  const.dir <- c("=", "=", rep(">=", 2*n.cat), rep("=", 2*ncol(mat)))
  
  ## get right hand side for the linear (in)equalities
  const.rhs <- c(1, 1, rep(min.prob, 2*n.cat), as, rep(0, ncol(mat)))
  
  ## return a list with the following items:
  ## 1. the coefficients for the objective function of the LP
  ## 2. the constraint matrix for the LP
  ## 3. the direction vector for the LP
  ## 4. the right hand side vector for the LP
  ## 5. the matrix passed into this function as an input (used in later functions)
  return(list(objective, const.mat, const.dir, const.rhs, mat))
}


## this function adds constraints to the LP such that the constraints already 
## specified for the LP input into this function will be satisfied. We recommend
## adding additional constraints one at a time to ensure there is a feasible 
## solution for the LP. A feasible solution can be found using lpSolveFn() below.
## the LPs output by initializeLP() and addConstraint() can both be input into 
## this function

## this function can be used to specify probabilities of certain combinations
## of the binary outcomes; these combinations are such that the "yes" outcomes
## are experienced and the "no" outcomes are not. An outcome need not be specified
## as "yes" or "no". Refer to the "yes" and "no" outcomes using the column numbers 
## corresponding to the outcomes in the matrix returned by createTree()
addConstraint <- function(lp, group, yes = NULL, no = NULL, dir, rhs){
  
  ## the inputs are described as follows
  ## lp: an LP output by initializeLP() or addConstraint()
  ## group: "a" or "b" to denote which group the new constraint is for
  ## yes: a vector of numbers corresponding to the outcomes that are experienced;
  ##      leave blank if specifying an event only characterized by "no" 
  ## no: a vector of numbers corresponding to outcomes that are not experienced;
  ##     leave blank if specifying an event only characterized by "yes"
  ## dir: one of "=", ">=", or "<=" to reflect the direction of the new constraint
  ## rhs: a number denoting the right hand side of the linear constraint
  
  if (length(lp) != 5){
    stop("lp must be a list of length 5 returned by intializeLP() or addConstraint()")
  } else if (!(group %in% c("a","b"))){
    stop("group must be 'a' or 'b'")
  } else if (is.null(yes) & is.null(no)){
    stop("one of yes or no must be non-null")
  } else if (!is.null(yes)){
    if (mean(yes %in% seq(1, nrow(lp[[5]]))) < 1){
      stop("all items in yes must correspond to outcome numbers")
    }
  } else if (!is.null(no)){
    if (mean(no %in% seq(1, nrow(lp[[5]]))) < 1){
      stop("all items in no must correspond to outcome numbers")
    }
  } else if (!is.null(no) & !is.null(yes)){
    if (mean(no %in% yes) > 0){
      stop("the same outcome cannot appear in yes and no")
    }
  } else if (length(dir) != 1){
    stop("dir must have length 1")
  } else if (!(dir %in% c("=", ">=", "<="))){
    stop("dir must be '=', '>=', or '<='")
  } else if (length(rhs) != 1){
    stop("rhs must have length 1")
  } else if (!is.numeric(rhs)){
    stop("rhs must be numeric")
  }
  
  ## extract information from input LP
  const.mat <- lp[[2]]
  const.dir <- lp[[3]]
  const.rhs <- lp[[4]]
  mat <- lp[[5]]
  n.cat <- nrow(mat)
  
  ## extract a list of mutlinomial categories that correspond to each outcome
  lst <- apply(mat, 2, function(x){which(x == 1)})
  
  ## get indices for the multinomial categories corresponding to yes outcomes
  yes_numbers <- if (!is.null(yes)) Reduce(intersect, lst[yes]) else numeric()
  
  ## get indices for the multinomial categories corresponding to no outcomes
  no_numbers <- if (!is.null(no)) unique(unlist(lst[no])) else numeric()
  
  ## get indices that correspond to yes outcomes but not no outcomes
  inds <- setdiff(yes_numbers, no_numbers)
  
  if (group == "b"){
    ## add new constraint for group B
    new.row <- c(1:n.cat %in% inds, rep(0, n.cat))
  }
  else{
    ## add new constraint for group A
    new.row <- c(rep(0, n.cat), 1:n.cat %in% inds)
  }
  
  ## append new row to constraint matrix and new element to 
  ## the constraint direction vector and RHS
  const.mat <- rbind(const.mat, new.row)
  const.dir <- c(const.dir, dir)
  const.rhs <- c(const.rhs, rhs)
  
  ## return five element list in the same format as initializeLP()
  ## this output can be used in addConstraint() or lpSolveFn()
  return(list(lp[[1]], const.mat, const.dir, const.rhs, mat))
}

## this function outputs a solution to the LP as a matrix where each
## column corresponds to a multinomial category and each row corresponds 
## to a group. If there is no feasible solution, a warning message will be
## output. In that event, you can relax or reconsider some of the initial or
## added constraints. For this reason, we commend checking the feasibility
## of the LP after every constraint is added

lpSolveFn <- function(input.lp){
  
  ## the input is described as follows
  ## input.lp: an lp output by initializeLP() or addConstraint()
  
  ## extract matrix from createTree()
  mat <- input.lp[[5]]

  ## solve LP using lp() in LP solve package
  lp.solution <- lp("max", input.lp[[1]], input.lp[[2]], 
                    input.lp[[3]], input.lp[[4]], compute.sens=TRUE)

  ## format the output probabilities for both groups in a matrix
  params <- matrix(lp.solution$solution[c(seq(nrow(mat) + 1, 2*nrow(mat)), seq(1, nrow(mat)))], 
                   nrow = 2, byrow = TRUE)
  
  rownames(params) <- c("Group A", "Group B")
  colnames(params) <- sprintf("p%d", 1:nrow(mat))
  
  if (sum(as.numeric(params))==0){
    warning("this LP has no feasible solution; please reconsider your added or initial constraints")
  }
  
  ## output a matrix with the multinomial probabilities for each group
  params
}