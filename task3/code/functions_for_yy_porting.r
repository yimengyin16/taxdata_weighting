
#****************************************************************************************************
#                functions - custom ####
#****************************************************************************************************
ns <- function(df) {names(df) %>% sort}

ht <- function(df, nrows=6) {head(df, nrows); tail(df, nrows)}


#****************************************************************************************************
#                function to solve for x values for a single stub ####
#****************************************************************************************************
runstub <- function(AGI_STUB, tolerances, nzcc, pufbase_state, log_dir, interim_results_dir) {
  stub <- AGI_STUB
  # constraint_names <- tolerances %>% filter(AGI_STUB==stub) %>% .$table_desc %>% str_remove(., "Number of") %>% str_sub(., 1, 35)
  
  # scale the constraints so that the RHS is a given number
  contraint_value <- 1000 # each constraint will be this number when scaled
  constraints_unscaled <- tolerances %>% filter(AGI_STUB==stub) %>% .$target
  constraint_scales <- ifelse(constraints_unscaled==0, 1, abs(constraints_unscaled) / contraint_value)
  constraints <- constraints_unscaled / constraint_scales
  
  # create nzcc for the stub, and on each record create i to index constraints and j to index variables (the x elements)
  nzcc_stub <- nzcc %>%
    ungroup %>% # just to be sure
    filter(AGI_STUB==stub) %>%
    arrange(constraint_name, RECID) %>%
    # NOTE!!: create i and j, each of which will be consecutive integersm where
    #   i gives the index for constraints
    #   j gives the index for the RECID (for the variables)
    mutate(i=group_indices(., constraint_name),
           j=group_indices(., RECID)) %>% # file-wide indexes based on sort of RECID and constraint_name
    mutate(nzcc_unscaled=nzcc,
           nzcc=nzcc_unscaled / constraint_scales[i])
  
  inputs <- list()
  inputs$p <- 2
  inputs$wt <- pufbase_state %>% filter(AGI_STUB==stub) %>% arrange(RECID) %>% .$weight_initial
  inputs$RECID <- pufbase_state %>% filter(AGI_STUB==stub) %>% arrange(RECID) %>% .$RECID
  inputs$constraint_coefficients_sparse <- nzcc_stub
  inputs$n_variables <- length(inputs$wt)
  inputs$n_constraints <- length(constraints)
  inputs$objscale <- 1e6 # scaling constant used in the various objective function functions
  inputs$constraint_scales <- constraint_scales
  
  xlb <- rep(0, inputs$n_variables) # arbitrary
  xub <- rep(100, inputs$n_variables) # arbitrary
  x0 <- rep(1, inputs$n_variables)
  
  tol <- tolerances %>% filter(AGI_STUB==stub) %>% .$tol_default
  
  clb <- constraints - abs(constraints) * tol
  cub <- constraints + abs(constraints) * tol
  
  eval_jac_g_structure <- define_jac_g_structure_sparse2(inputs$constraint_coefficients_sparse, ivar="i", jvar="j")
  eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian
  
  logfile_name <- paste0(log_dir, "stub_", str_pad(stub, width=2, side="left", pad="0"), ".out")
  # ma86 was faster in one test I did
  opts <- list("print_level" = 0,
               "file_print_level" = 5, # integer
               "linear_solver" = "ma27", # mumps pardiso ma27 ma57 ma77 ma86 ma97
               "max_iter"= 100,
               "output_file" = logfile_name)
  
  result <- ipoptr(x0 = x0,
                   lb = xlb,
                   ub = xub,
                   eval_f = eval_f_xtop, # arguments: x, inputs
                   eval_grad_f = eval_grad_f_xtop,
                   eval_g = eval_g, # constraints LHS - a vector of values
                   eval_jac_g = eval_jac_g,
                   eval_jac_g_structure = eval_jac_g_structure,
                   eval_h = eval_h_xtop, # the hessian is essential for this problem
                   eval_h_structure = eval_h_structure,
                   constraint_lb = clb,
                   constraint_ub = cub,
                   opts = opts,
                   inputs = inputs)
  saveRDS(result, paste0(interim_results_dir, "ipopt_", str_pad(stub, width=2, side="left", pad="0"), ".rds"))
  print(result$message)
  
  df <- tibble(AGI_STUB=stub, RECID=inputs$RECID, wt_init=inputs$wt, x=result$solution)
  return(df)
}


#****************************************************************************************************
#                functions for ipoptr: constraint evaluation and coefficients ####
#****************************************************************************************************
eval_g <- function(x, inputs) {
  # constraints that must hold in the solution - just give the LHS of the expression
  # return a vector where each element evaluates a constraint (i.e., sum of (x * a ccmat column), for each column)
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  # inputs$constraint_coefficients_sparse has the fixed constraint coefficients in sparse form in a dataframe that has:
  #   i -- the constraint number
  #   j -- index into x (i.e., the variable number)
  #   nzcc -- the nonzero constraint coefficient
  
  constraints <- inputs$constraint_coefficients_sparse %>%
    group_by(i) %>%
    summarise(constraint_value=sum(nzcc * x[j]))
  
  return(constraints$constraint_value)
}


eval_jac_g <- function(x, inputs){
  # the Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # this function evaluates the Jacobian at point x
  
  # return: a vector where each element gives a NONZERO partial derivative of constraints wrt change in x
  # so that the first m items are the derivs with respect to each element of first column of ccmat
  # and next m items are derivs with respect to 2nd column of ccmat, and so on
  # so that it returns a vector with length=nrows x ncolumns in ccmat
  
  # because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  return(inputs$constraint_coefficients_sparse$nzcc)
}


define_jac_g_structure_sparse2 <- function(nzcc_df, ivar="i", jvar="j"){
  # the jacobian 
  # return a list that defines the non-zero structure of the "virtual" constraints coefficient matrix
  # the list has 1 element per constraint
  #   each element of the list has a vector of indexes indicating which x variables have nonzero constraint coefficents
  
  # nzcc_df is a nonzero constraints coefficients data frame
  # ivar gives the variable name for the integer index indicating each CONSTRAINT
  # jvar gives the variable name (character) for the integer index indicating the nonzero x variables for that constraint
  
  jac_sparse <- dlply(nzcc_df, ivar, function(x) return(x[[jvar]]))
  
  return(jac_sparse)
} 


#****************************************************************************************************
#                x^p + x^-p {xtop} -- functions for ipoptr ####
#****************************************************************************************************
eval_f_xtop <- function(x, inputs) {
  # objective function - evaluates to a single number
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  # here are the objective function, the 1st deriv, and the 2nd deriv
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  obj <- sum(w * {x^p + x^(-p) -2})
  
  return(obj)
}


eval_grad_f_xtop <- function(x, inputs){
  # gradient of objective function - a vector length x 
  # giving the partial derivatives of obj wrt each x[i]
  
  # ipoptr requires that ALL functions receive the same arguments, so I pass the inputs list to ALL functions
  
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  gradf <- w * (p * x^(p-1) - p * x^(-p-1))
  
  return(gradf)
}


eval_h_xtop <- function(x, obj_factor, hessian_lambda, inputs){
  # The Hessian matrix has many zero elements and so we set it up as a sparse matrix
  # We only keep the (potentially) non-zero values that run along the diagonal.
  
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  hess <- obj_factor * 
    { p*w*x^(-p-2) * ((p-1)*x^(2*p)+p+1) }
  
  return(hess)
}