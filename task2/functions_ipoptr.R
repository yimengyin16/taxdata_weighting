


# Functions to compute constraint coefficients
nnz  <- function(var, weight) {(var != 0) * weight} # weighted number of returns with non-zero values
nz   <- function(var, weight) {(var != 0) * weight} # weighted number of returns with zero values
npos <- function(var, weight) {(var >= 0) * weight} # weighted number of returns with positive values}
nneg <- function(var, weight) {(var <= 0) * weight} # weighted number of returns with negative values

sumval <- function(var, weight) { weight * var} # wegithed sum of non zero value
sumpos <- function(var, weight) {(var > 0) * weight * var}
sumneg <- function(var, weight) {(var < 0) * weight * var}



# Objective function ----

eval_f_xtop <- function(x, inputs) {
  # .. objective function - evaluate to a single number ----
  # returns a single value
  
  #ipoptr requires that ALL functions receive the same arguments, so a list called inputs is passed to ALL functions
  
  # here are the objective function, the 1st deriv, and the 2nd deriv
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  p <- inputs$p
  w <- inputs$weight
  
  obj <- sum(w * (x^p + x^(-p) -2))
  
}


# Gradient
eval_grad_f_xtop <- function(x, inputs){
  #.. gradient of objective function - a vector length x ----
  # giving the partial derivatives of obj wrt each x[i]
  # returns one value per element of x
  
  # ipoptr requires that ALL functions receive the same arguments, so a list called inputs is passed to ALL functions
  
  # here are the objective function, the 1st deriv, and the 2nd deriv
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$weight
  
  gradf <- w * (p * x^(p-1) - p * x^(-p-1))
  
}


# Hessian

eval_h_xtop <- function(x, obj_factor, hessian_lambda, inputs){
  # The Hessian matrix ----
  # The Hessian matrix has many zero elements and so we set it up as a sparse matrix
  # We only keep the (potentially) non_zero values that run along the diagonal. 
  # the Hessian is returned as a long vector. Separately, we define which elements of this vector correspond to which cells of the Hessian matrix
  
  # obj_factor and hessian_lambda are required arguments of the function. They are created within ipoptr (?) - we do not create them
  
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # Make it easier to read
  p <- inputs$p
  w <- inputs$weight
  
  hess <- obj_factor * 
    ( p*w*x^(-p-2) * {(p-1)*x^(2*p) + p + 1} )
}

hess_structure <- function(n_variables) {lapply(1:n_variables, function(x) x)} # diagonal elements of our Hessian





# Constraints and Jacobian --  dense matrix ----


calc_constraints <- function(ccoef, constraint_vars, x = rep(1, nrow(ccoef))){
  # calculate weighted sums of constraint_vars that are in data, using an optional multiplier x
  # ccoef: df with 1 colum per constraint variable (plus possibly additional conlums); 1 row per person (tax return)
  # constraint_vars: character vector
  # x: numeric vector of adjustment factor
  # return: a named vector of constraint sums
  colSums(x * select(ccoef, !!constraint_vars) )
}


eval_g_dense <- function(x, inputs){
  #.. constraints that must hold in the solution ----
  # Just give the lLHS of the expression
  # return a vector where each element evaluates a constraint (i.e., sum of (x * a cc matrix column), for each column)
  
  # ipoptr requires that ALL functions receive the same arguments, so a list called inputs is passed to ALL functions
  
  # This is the dense version in that it uses a dense set of constraint coefficients, inputs$ccoef (passed to this 
  # function as a data frame but we could easily use a matrix instead). It has
  #   one column for every constraint
  #   one row for every variable
  # Many of the cells are likely to be zero.
  
  unname(calc_constraints(inputs$ccoef, inputs$constraint_vars, x)) # ipoptr doesn't seem to like a named vector
}

eval_jac_g <- function(x, inputs){
  # Jacobian of the constraints ----
  # The Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # This function evaluates the Jacobian at point x
  # It returns a vector with these partial derivatives 
  
  # THis is the dense version that returns a vector with one element for Every item in the Jacobian (including zeros)
  # Thus the vector has n_constraints * n_variables elements
  # First, all of the elements for constraint 1, then all for constraint 2, etc...
  
  # Because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requres that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  inputs$jac_vector
}


jac_flatten <- function(jac_df, nz = FALSE) {
  # Flatten the Jacobian matrix
  jac_vector <- c(as.matrix(jac_df))
  if(nz) jac_vector <- jac_vector[jac_vector != 0]
  return(jac_vector)
} 

define_jac_g_structure_dense <- function(n_constraints, n_variables){
  # .. function to define the structure of the Jacobian ---
  # list with n_constraints elements
  # each is 1:n_variables
  lapply(1:n_constraints, function(n_constraints, n_variables) 1:n_variables, n_variables)
}





make.sparse.structure <- function(ccoef) {
  # argument is a constraint coefficients data frame (or could be a matrix)
  #  -- columns are constraints
  #  -- rows are variables
  # return a list with indeces for nonzero constraint coefficients
  
  # this is much faster than the make.sparse function in ipoptr
  
  f <- function(x) which(x != 0, arr.ind = TRUE)
  
  rownames(ccoef) <- NULL # just to be safe
  indexes <- apply(ccoef, 2, f)
  #indexes.list <- as.list(indexes) # not necessary
}




get_ccoefs <- function(recipe, data){
  # recipe as passed to this function will have only one record - 
  #  that is, it has rules for a single constraint
  # Returns a data frame containing all non-zero coefficients of the constraint
  #data <- acssub
  
  recipe <- as_tibble(recipe) # force this - otherwise passed as list
  
  datause <- 
    data %>%  
    mutate(row_num = row_number()) %>% 
    filter(eval(parse(text = recipe$valgroup))) %>% 
    mutate(vname = recipe$vname) %>% 
    select(vname, row_num, pwgtp, value = !!recipe$vname)
  
  ccoef <- left_join(recipe, datause, by = "vname") %>% 
    mutate(ccoef = case_when(recipe$fn == "sumval" ~ value * pwgtp,
                             recipe$fn == "npos"   ~ (value > 0) * pwgtp,
                             recipe$fn == "nneg"   ~ (value < 0) * pwgtp
    )) %>% 
    # keep only the the nonzero coefficients
    filter(ccoef != 0)
  
  return(ccoef)
}


get_ccoefs2 <- function(recipe, data){
  # recipe as passed to this function will have only one record - 
  #  that is, it has rules for a single constraint
  # Returns a data frame containing all non-zero coefficients of the constraint
  #data <- acssub
  
  recipe <- as_tibble(recipe) # force this - otherwise passed as list
  
  datause <- 
    data %>%  
    mutate(row_num = row_number()) %>% 
    filter(eval(parse(text = recipe$valgroup))) %>% 
    mutate(vname = recipe$vname) %>% 
    select(stabbr, vname, row_num, pwgtp, value = !!recipe$vname)
  
  ccoef <- left_join(recipe, datause, by = "vname") %>% 
    mutate(ccoef = case_when(recipe$fn == "sumval" ~ value * pwgtp,
                             recipe$fn == "npos"   ~ (value > 0) * pwgtp,
                             recipe$fn == "nneg"   ~ (value < 0) * pwgtp
    )) %>% 
    # keep only the the nonzero coefficients
    filter(ccoef != 0)
  
  return(ccoef)
}


eval_g_sparse <- function(x, inputs){
  #.. constraints that must hold in the solution ----
  # just give the LHS of the expression
  # return a vector where each element evaluates a constraint (i.e., sum of (x * a cc matrix column), for each column)
  
  # ipoptr requires that ALL functions receive the same arguments, so a list called inputs is passed to ALL functions
  
  # This is the sparse version
  
  constraint_vals <- inputs$ccoef_sparse %>% 
    mutate(x = x[row_num]) %>% # extract adjustment factor from x corresponding to non-zero ccoefficients
    group_by(constraint_num) %>% 
    summarise(constraint_val = sum(x * ccoef))
  
  return(constraint_vals$constraint_val)
}



# now we need to flatten the constraint coefficients
jac_flatten_sparse <- function(ccoef_sparse){
  ccoef_sparse %>% 
    ungroup %>% 
    arrange(constraint_num, row_num) %>% 
    pull(ccoef)
}

define_jac_g_structure_sparse <- function(ccoef_sparse){
  #.. function to define the structure of the Jacobian ---
  # list with n_constraints elements
  # each is a vector indicating row indexes relevant to the constraint
  
  f <- function(cnum, ccoef_sparse){
    ccoef_sparse %>% 
      filter(constraint_num == cnum) %>% 
      pull(row_num)
  }
  
  jac_j_structure <- plyr::llply(unique(ccoef_sparse$constraint_num),
                                 f,
                                 ccoef_sparse = ccoef_sparse)
  return(jac_j_structure)
}













