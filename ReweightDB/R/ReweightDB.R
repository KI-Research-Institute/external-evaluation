# suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(CVXR)
library(glue)


#' Reweight an internal database to match the means of an external one.
#'
#' @description
#'
#' This function receives a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#' @param Z a data frame where every row stores a sample of the internal
#' databse.
#' @param mu a vector of means of the external dataset.
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy,
#' \eqn{-\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#'
#'
#' @param optimizationMethod primal or dual. Currently dual works only with entropy divergence
#'
#' @return
#' A vector of weights
#'
#' @export
reweightByMeans <- function(
  Z, mu, divergence = 'entropy', lambda=1e-6, minSd=1e-4, minW=1e-6, distance = 'l2', optimizationMethod = 'primal',
  verbose=FALSE) {
  # TODO in the regularized case, hyper param may depend on the number of features/samples
  # Find optimal weights

  if (optimizationMethod == 'primal')
    w_hat <- primalReweightByMeans(Z, mu, divergence, lambda, minSd, minW, distance, verbose)
  else
    w_hat <- dualReweightByMeans(Z, mu, lambda, minSd, minW, verbose)

  min_w <- min(w_hat)
  mean_w <- mean(w_hat)
  if (verbose)
    cat(sprintf('mean(w) = %.2f (sd = %.2f, min = %.2g, max = %.2g)\n',
                mean_w, sd(w_hat), min_w, max(w_hat)))
  if ((!is.na(min_w)) & (min_w < 0)) {
    warning("Trimming negative weights to zero")
    w_hat[w_hat<0] <- minW/n
  }
  return (w_hat)
}


#' Reweight an internal database to match the means of an external via solving a corresponding optimization problem.
#'
#' @description
#'
#' This function recieves a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#'
#' SEE THE PARENT FUNCTION
#'
#' @return
#' A vector of weights
#'
#' @export
primalReweightByMeans <- function(Z, mu, divergence,lambda, minSd, minW, distance, verbose) {
  normalized <- normalizeDataAndExpectations(Z, mu, minSd)
  n <- nrow(normalized$Z)
  w <- Variable(n)

  if (divergence == 'entropy')    fDivergence <- -mean(entr(w))
  else if (divergence == 'chi2')  fDivergence <- norm2(w-(1/n)) ** 2
  else                            stop(glue("unsuported divergence type {divergence}"))

  if (distance == 'l2')       expectationsDistance <- norm2(t(normalized$Z) %*% w - normalized$mu)
  else if (distance == 'l1')  expectationsDistance <- norm1(t(normalized$Z) %*% w - normalized$mu)
  else                        stop(glue("unsuported distance type {distance}"))

  if (lambda > 0) {
    if (verbose) cat(glue('Reweighting using {divergence}, {distance}, lambda = {lambda}, minW = {minW}'), '\n')
    objective <- Minimize(expectationsDistance + lambda*fDivergence)
    constr <- list(w >= minW, sum(w) == 1)
  } else {
    if (verbose) cat(glue('Reweighting using {divergence}, hard expectation constraints, minW = {minW}'), '\n')
    objective <- Minimize(fDivergence)
    constr <- list(w >= minW, sum(w) == 1, (t(normalized$Z) %*% w) == normalized$mu)
  }
  problem <- Problem(objective, constraints = constr)
  result <- solve(problem)
  # The status of the solution can be "optimal", "optimal_inaccurate", "infeasible", "infeasible_inaccurate",
  # "unbounded", "unbounded_inaccurate", or "solver_error"
  if (result$status != 'optimal') {
    warning(glue('non-optimal results, returning NaNs, data size = {n} * {length(mu)}'))
    w_hat <- rep(NaN, n)
  } else {
    w_hat <- result$getValue(w) * n
  }

  return (w_hat)
}


#' Reweight an internal database to match the means of an external one using
#' a dual formulation of a weight optimization problem.
#'
#' @description
#'
#' This function recieves a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#' @param Z a data frame where every row stores a sample of the internal
#' databse.
#' @param mu a vector of means of the internal dataset.
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#'
#' @return
#' A vector of weights
#'
#' @export
dualReweightByMeans <- function(X, mu, lambda, minSd, minW, verbose) {
  normalized <- normalizeDataAndExpectations(X, mu, minSd)
  m <- ncol(normalized$Z)
  n <- nrow(normalized$Z)

  nu <- Variable(m+1)
  C <- rbind(t(normalized$Z), matrix(1,1,n))
  d <- c(normalized$mu, 1)
  objective <- Minimize(t(d) %*% nu + exp(-1) * sum_entries(exp(- t(C) %*% nu)))
  if (lambda==0)
    problem <- Problem(objective)
  else {
    constr <- list(norm2(nu[1:m]) <= (1/lambda))
    problem <- Problem(objective, constr = constr)
  }
  result <- solve(problem)

  if (result$status != 'optimal') {
    warning(glue('non-optimal results, returning NaNs, data size = {n} * {length(mu)}'))
    w_hat <- rep(NaN, n)
  } else {
      nu_hat <- result$getValue(nu)
      cat('nu:', nu_hat[1], ',' , nu_hat[2],  ',' , nu_hat[3], ', ... ,', nu_hat[m],  ',',nu_hat[m+1], '\n')
      if (!any(is.na(nu_hat))) {
          w_hat <- exp(-1- t(C) %*% nu_hat) * n
      } else {
        warning(glue('Optimization resulted in NaNs, data size = {n} * {length(mu)}'))
        w_hat <- rep(NaN, n)
      }
  }
  return(w_hat)
}


maximizeWeightedObj <- function(X, b, loss, lambda=1, alpha=0, minSd=1e-4, minW=1e-6, verbose=FALSE) {
  # TODO - check correct usage of normalized expectations
  # TODO set a hyper param that depends on the number of features
  # Find optimal weights
  normalized <- normalizeDataAndExpectations(X, b, minSd)
  n <- nrow(normalized$Z)
  w <- Variable(n)
  if (alpha==0) {
    objective <- Maximize( (t(loss) %*% w) + lambda * sum(entr(w)) )
    constr <- list(w >= 0, sum(w) == 1, (t(normalized$Z) %*% w) == normalized$mu)
  } else {
    objective <- Maximize( (t(loss) %*% w) + lambda * sum(entr(w)) - alpha * norm2((t(normalized$Z) %*% w) - normalized$mu) )
    constr <- list(w >= 0, sum(w) == 1)
  }

  problem <- Problem(objective, constraints = constr)
  result <- solve(problem)
  if (result$status == 'solver_error') {
    warning(glue('Solver error, returning NaNs, data size = {n} * {length(b)}'))
    w_hat <- rep(NaN, n)
  } else {
    w_hat <- result$getValue(w) * n
  }
  if (verbose)
    cat(sprintf(
      'mean(w) = %.2f (sd = %.2f, min = %.2g, max = %.2g, n=%d)\n',
      mean(w_hat), sd(w_hat), min(w_hat), max(w_hat), length(w_hat)))
  min_w <- min(w_hat)
  if (is.na(min_w)) {
    warning(sprintf("Could not solve optimization problem, n = %d, q = %d", n, ncol(normalized$Z)))
  } else {
    if (min_w < 0) {
      warning("Trimming negative weights to zero")
      w_hat[w_hat<0] <- minW/n
    }
  }
  return (w_hat)
}


#' Normalize data and expectations
#'
#' @description
#'
#' Normalize data and expectations. 
#'
#' @param z data frame
#' @param mu vector of expectations
#' @param minSd float value of minimum standard deviation. Columns with a smaller sd will be removed
#'
#' @return
#' list
#'
#' @export
normalizeDataAndExpectations <- function(Z, mu, minSd) {
  # Preprocess and filter columns with low variance
  muZ <- colMeans(Z)
  sdZ <- apply(Z, 2, sd)
  useCols <- sdZ>minSd # remove columns with low variance
  if (sum(!useCols) > 0) {
    warning(
      sprintf("Removing %d columns with low variance (<%.2g). PLEASE EXAMINE THESE COLUMNS.", sum(!useCols), minSd))
    # print(sdX[!useCols])
    Z <- Z[,names(Z)[useCols]]  # remove low variance columns todo: flag a positivity issue
    sdZ <- sdZ[useCols]
    muZ <- muZ[useCols]
    mu <- mu[useCols]
  }
  # Normalize
  p <- ncol(Z)
  for (i in 1:p) Z[,i] <- (Z[,i]-muZ[i])/sdZ[i]
  mu <- (mu-muZ)/sdZ
  return(list(Z=Z, mu=mu, muZ=muZ, sdZ=sdZ))
}


#' Compute Table 1 like transformation
#'
#' @description
#'
#' Compute features and squared numeric feature multiplied by outcome and 1-outcome. 
#'
#' @param x data frame
#' @param outcomeBalance boolean specifying if to use outcome in balancing
#' @param outcomeCol string specifying the name of the outcome column
#'
#' @return
#' A data frame
#'
#' @export
computeTable1LikeTransformation <- function(X, outcomeBalance, outcomeCol='Y') {
  # Add squares of numeric features
  is_numeric <- sapply(X, function(x) length(unique(x))>2)  # TODO - consider a more elegant way
  for (s in names(is_numeric[is_numeric])) {
    sNew <- paste(s, '_2', sep='')
    X[[sNew]] <- (X[[s]]**2)
  }
  # Convert binary variables to numeric 0-1
  is_factor <- sapply(X, is.factor)
  X[is_factor] <- sapply(X[is_factor], as.numeric)
  X[is_factor] <- X[is_factor] - 1
  # Add interactions with outcome
  if (outcomeBalance) {
    Z = data.frame(row.names = row.names(X))
    Z[[outcomeCol]] <- X[[outcomeCol]]
    x_names <- names(X)
    ext_x <- x_names[x_names != outcomeCol]
    ext_x_y1 <- paste(ext_x, '_y1', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y1[i]]] <- X[[outcomeCol]] * X[[ext_x[i]]]
    }
    ext_x_y0 <- paste(ext_x, '_y0', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y0[i]]] <- (1-X[[outcomeCol]]) * X[[ext_x[i]]]
    }
  } else {
    Z <- X
  }
  return (Z)
}
