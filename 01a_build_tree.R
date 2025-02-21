## This function creates a similar diagram to Figure 1 in the paper. It outputs
## a list with two items: the first item can be plotted in the graphics console to
## visualize the outcomes from the experiment, and the second item is a matrix that
## characterizes which multinomial categories correspond to which binary outcomes
## (i.e., if the ith row in the jth column is a 1, then the ith category corresponds 
## to the jth outcome). This matrix is used in later functions.

## load necessary libraries
require(DiagrammeR)
require(data.tree)

## this is the function to create the tree diagram and matrix
## for each node, you will be prompted to enter one of three options: "y", "n", or "b"
## if the new node represents a valid combination of outcomes, enter "y"
## if the new node does not represent a valid combination of outcomes, enter "n"
## if you need to go back to the previous node to correct a mistake, enter "b"

createTree <- function(name, outcomes){
  
  ## the inputs are described as follows 
  ## name: a single character string that names the experiment
  ## outcomes: a vector of character strings, where each string represents an outcome
  
  if (length(name) != 1){
    stop("name must be a single string")
  } else if (!is.character(name)){
    stop("name must be a single string")
  } else if (length(outcomes) <= 1){
    stop("there must be at least two outcomes")
  } 
  
  ## initialize tree with yes and no options for first outcome
  root <- Node$new(name)
  level1.no <- root$AddChild(paste0(outcomes[1], ": No"))
  level1.yes <- root$AddChild(paste0(outcomes[1], ": Yes"))
  
  ## set up numeric coding system that is used to create the matrix output
  level1.no$index <- 10^(length(outcomes))
  level1.yes$index <- 10^(length(outcomes)) + 10^(length(outcomes) - 1)
  
  ## set up visual formatting for the plot
  SetGraphStyle(root, rankdir = "TB")
  SetEdgeStyle(root, arrowhead = "vee", color = "grey35", penwidth = 1)
  SetNodeStyle(root, style = "filled,rounded", shape = "box", fillcolor = "White", 
               fontname = "helvetica", tooltip = GetDefaultTooltip, fontcolor = "Black")
  
  ## initialize the input flags from the console
  level1 <- c(level1.no, level1.yes)
  flag.yes <- "a"
  flag.no <- "a"
  
  ## loop through each of the remaining outcomes
  for (j in 2:length(outcomes)){
    assign(paste0("level", j), NULL)}
  
  ## initialize back indicator as false
  go.back <- FALSE
  for (j in 2:length(outcomes)){
    for (i in 1:length(get(paste0("level", j-1)))){
      ## start with adding a no child to the first available parent
      assign(paste0("level",j,".no", i), 
             get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": No")))
      
      ## use coloured formatting to denote which node is being considered
      SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "LightBlue", penwidth = 1)
      
      ## plot and solicit output from user
      print(plot(root))
      flag.no<-readline(prompt="Is this combination valid (y/n/b)? " )
      flag.no<-tolower(as.character(flag.no))
      
      ## make sure response is valid
      while(!(flag.no %in% c("y", "n", "b"))){
        flag.no<-readline(prompt="You must enter y, n, or b: " )
        flag.no<-tolower(as.character(flag.no))
      }
      
      ## removed the coloured formatting from the previous node
      SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "White", penwidth = 1)
      
      ## go back to the previous node if necessary
      if (flag.no == "b") {
        assign(paste0("level",j,".no", i), 
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
        go.back <- TRUE
        break
      }
      ## if yes response, retain the node and index using numeric coding to create the matrix
      else if(flag.no == "y"){
        assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".no", i))))
        obj_name <- paste0("level", j, ".no", i)
        obj <- get(obj_name)
        obj$index <- get(paste0("level", j - 1))[[i]]$index
        assign(obj_name, obj)
      }
      else{
        ## if no response, the remove the node from the tree
        assign(paste0("level",j,".no", i), 
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
      }
      
      ## repeat this process for the yes child of the current parent node
      assign(paste0("level",j,".yes", i), 
             get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": Yes")))
      
      SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "LightBlue", penwidth = 1)
      
      print(plot(root))
      
      flag.yes<-readline(prompt="Is this combination valid (y/n/b)? " )
      flag.yes<-tolower(as.character(flag.yes))
      
      while(!(flag.yes %in% c("y", "n", "b"))){
        flag.yes<-readline(prompt="You must enter y, n, or b: " )
        flag.yes<-tolower(as.character(flag.yes))
      }
      
      SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "White", penwidth = 1)
      
      if (flag.yes == "b") {
        assign(paste0("level",j,".yes", i), 
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
        go.back <- TRUE
        break
      }
      else if(flag.yes == "y"){
        assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".yes", i))))
        obj_name <- paste0("level", j, ".yes", i)
        obj <- get(obj_name)
        obj$index <- get(paste0("level", j - 1))[[i]]$index + 10^(length(outcomes) - j)
        assign(obj_name, obj)
      }
      else{
        assign(paste0("level",j,".yes", i), 
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
      }
      
      ## if no response was provided to both the yes and no children, prompt the user to 
      ## go back (each non-terminal node must have at least one child)
      if (length(get(paste0("level", j-1))[[i]]$children) == 0){
        cat("Error: The previous node needs a child. Please enter 'b' to go back.\n")
      }
      
    }
    if (go.back == TRUE){
      break
    }
  }
  
  ## if the user needed to go back, then we had to break the for loop above
  ## it is reinitialized here depending on whether the user went back
  first.node <- FALSE
  while (go.back){
    ## this indexing determines what the previous node was so that we can
    ## repopulate the correct plot. The while loop ensures that we can 
    ## continue to go back if necessary. The remainder of this process is
    ## similar to what was in the first for loop
    go.back <- FALSE
    prev.j <- j
    if (flag.yes == "b"){
      prev.i <- i
    }
    else {
      prev.i <- i - 1
    }
    if (prev.i == 0){
      prev.j <- j - 1
      if (prev.j == 1){
        prev.j <- 2
        prev.i <- 1
        if (flag.no == "b"){
          cat("Error: Cannot go back from the first node.\n")
          first.node <- TRUE
        }
      }
      else {
        prev.i <- length(get(paste0("level", j-2)))
      }
    }
    j <- prev.j
    if (first.node == FALSE){
      length.temp <- length(get(paste0("level", prev.j)))
      if (flag.no == "b"){
        parent.temp <- get(paste0("level",prev.j,".yes", prev.i))$parent
      }
      else{
        parent.temp <- get(paste0("level",prev.j,".no", prev.i))$parent
      }
      if (!is.null(parent.temp)){
        if (length.temp == 1){
          assign(paste0("level", prev.j), NULL)
        }
        else{
          assign(paste0("level", prev.j), get(paste0("level", j))[1:(length.temp - 1)])
        }
      }
    }
    for (i in prev.i:length(get(paste0("level", j-1)))){
      if (flag.yes == "b"){
        flag.yes <- "a"
      }
      if (!(i == prev.i & flag.no == "b")){
        assign(paste0("level",j,".no", i),
               get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": No")))

        SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "LightBlue", penwidth = 1)

        print(plot(root))

        flag.no<-readline(prompt="Is this combination valid (y/n/b)? " )
        flag.no<-tolower(as.character(flag.no))

        while(!(flag.no %in% c("y", "n", "b"))){
          flag.no<-readline(prompt="You must enter y, n, or b: " )
          flag.no<-tolower(as.character(flag.no))
        }

        SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "White", penwidth = 1)

        if (flag.no == "b") {
          assign(paste0("level",j,".no", i),
                 get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
          go.back <- TRUE
          break
        }
        else if(flag.no == "y"){
          assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".no", i))))
          obj_name <- paste0("level", j, ".no", i)
          obj <- get(obj_name)
          obj$index <- get(paste0("level", j - 1))[[i]]$index
          assign(obj_name, obj)
        }
        else{
          assign(paste0("level",j,".no", i),
                 get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
        }
      }
      else if (first.node == TRUE){
        first.node <- FALSE
        assign(paste0("level",j,".no", i),
               get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": No")))
        
        SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "LightBlue", penwidth = 1)
        
        print(plot(root))
        
        flag.no<-readline(prompt="Is this combination valid (y/n/b)? " )
        flag.no<-tolower(as.character(flag.no))
        
        while(!(flag.no %in% c("y", "n", "b"))){
          flag.no<-readline(prompt="You must enter y, n, or b: " )
          flag.no<-tolower(as.character(flag.no))
        }
        
        SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "White", penwidth = 1)
        
        if (flag.no == "b") {
          assign(paste0("level",j,".no", i),
                 get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
          go.back <- TRUE
          break
        }
        else if(flag.no == "y"){
          assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".no", i))))
          obj_name <- paste0("level", j, ".no", i)
          obj <- get(obj_name)
          obj$index <- get(paste0("level", j - 1))[[i]]$index
          assign(obj_name, obj)
        }
        else{
          assign(paste0("level",j,".no", i),
                 get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
        }
      }

      if (flag.no == "b"){
        flag.no <- "a"
      }

      assign(paste0("level",j,".yes", i),
             get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": Yes")))

      SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "LightBlue", penwidth = 1)

      print(plot(root))

      flag.yes<-readline(prompt="Is this combination valid (y/n/b)? " )
      flag.yes<-tolower(as.character(flag.yes))

      while(!(flag.yes %in% c("y", "n", "b"))){
        flag.yes<-readline(prompt="You must enter y, n, or b: " )
        flag.yes<-tolower(as.character(flag.yes))
      }

      SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "White", penwidth = 1)

      if (flag.yes == "b") {
        assign(paste0("level",j,".yes", i),
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
        go.back <- TRUE
        break
      }
      else if(flag.yes == "y"){
        assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".yes", i))))
        obj_name <- paste0("level", j, ".yes", i)
        obj <- get(obj_name)
        obj$index <- get(paste0("level", j - 1))[[i]]$index + 10^(length(outcomes) - j)
        assign(obj_name, obj)
      }
      else{
        assign(paste0("level",j,".yes", i),
               get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
      }
      
      if (length(get(paste0("level", j-1))[[i]]$children) == 0){
        cat("Error: The previous node needs a child. Please enter 'b' to go back.\n")
      }

    }
    if (go.back != TRUE & prev.j < length(outcomes)){
      for (j in (prev.j+1):length(outcomes)){
        for (i in 1:length(get(paste0("level", j-1)))){
            assign(paste0("level",j,".no", i),
                   get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": No")))

            SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "LightBlue", penwidth = 1)

            print(plot(root))

            flag.no<-readline(prompt="Is this combination valid (y/n/b)? " )
            flag.no<-tolower(as.character(flag.no))

            while(!(flag.no %in% c("y", "n", "b"))){
              flag.no<-readline(prompt="You must enter y, n, or b: " )
              flag.no<-tolower(as.character(flag.no))
            }

            SetNodeStyle(get(paste0("level",j,".no", i)), fillcolor = "White", penwidth = 1)

            if (flag.no == "b") {
              assign(paste0("level",j,".no", i),
                     get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
              go.back <- TRUE
              break
            }
            else if(flag.no == "y"){
              assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".no", i))))
              obj_name <- paste0("level", j, ".no", i)
              obj <- get(obj_name)
              obj$index <- get(paste0("level", j - 1))[[i]]$index
              assign(obj_name, obj)
            }
            else{
              assign(paste0("level",j,".no", i),
                     get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": No")))
            }

          if (flag.no == "b"){
            flag.no <- "a"
          }

          assign(paste0("level",j,".yes", i),
                 get(paste0("level", j-1))[[i]]$AddChild(paste0(outcomes[j], ": Yes")))

          SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "LightBlue", penwidth = 1)

          print(plot(root))

          flag.yes<-readline(prompt="Is this combination valid (y/n/b)? " )
          flag.yes<-tolower(as.character(flag.yes))

          while(!(flag.yes %in% c("y", "n", "b"))){
            flag.yes<-readline(prompt="You must enter y, n, or b: " )
            flag.yes<-tolower(as.character(flag.yes))
          }

          SetNodeStyle(get(paste0("level",j,".yes", i)), fillcolor = "White", penwidth = 1)

          if (flag.yes == "b") {
            assign(paste0("level",j,".yes", i),
                   get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
            go.back <- TRUE
            break
          }
          else if(flag.yes == "y"){
            assign(paste0("level", j), c(get(paste0("level", j)), get(paste0("level",j,".yes", i))))
            obj_name <- paste0("level", j, ".yes", i)
            obj <- get(obj_name)
            obj$index <- get(paste0("level", j - 1))[[i]]$index + 10^(length(outcomes) - j)
            assign(obj_name, obj)
          }
          else{
            assign(paste0("level",j,".yes", i),
                   get(paste0("level", j-1))[[i]]$RemoveChild(paste0(outcomes[j], ": Yes")))
          }
          
          if (length(get(paste0("level", j-1))[[i]]$children) == 0){
            cat("Error: The previous node needs a child. Please enter 'b' to go back.\n")
          }

        }
        if (go.back == TRUE){
          break
        }
      }
    }
  }
  
  ## process to create the binary matrix described earlier
  ## the first column is all ones for formatting purposes, so 
  ## we will remove it later
  m <- length(outcomes)
  indices <- matrix(0, nrow = length(get(paste0("level", m))), 
                    ncol = length(outcomes) + 1)
  for (i in 1:length(get(paste0("level", m)))){
    temp_index <- as.character(format(get(paste0("level", m))[[i]]$index, scientific = FALSE))
    indices[i,] <- as.numeric(unlist(strsplit(temp_index, "")))
  }
  
  ## plot final tree for experiment
  print(plot(root))
  
  ## return the tree object for future plotting and the matrix for use later
  return(list(root, indices[, -1]))
}