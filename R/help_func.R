softmax <- function(x){

  matrix(exp(c(0, x)) / sum(exp(c(0, x))), ncol = 1)

}

# Log likelihood function on a logistic scale


log_prev <- function(Q, obs, b){


  -t(obs) %*% log(Q %*% softmax(b))


}


# log-likelihood for multinomial regression model (on logistic scale)

logl_reg <- function(D,  P_trans, b ){

  beta_hat <- matrix(b, nrow = ncol(D))

  pr   <- t(apply(D %*% beta_hat, 1, softmax))

  logl <- log(rowSums(P_trans * pr))

  return(-sum(logl))

}


# extract predicted probabilities from the model object `rrmult` conditional on data and return a data frame.

fitted_probs <- function(designmat, coeff){

  pred      <- apply(designmat %*% coeff, 1, softmax)

  rownames(pred)   <- c("never", "former", "Last year")


  invisible(data.frame(t(pred)))
}


#help function to calculate marginal effect of numeric covariate. need to apply transformation on data

marginal_num <- function(model, data, variable, eps = 1e-7) {

  # set value of `h` based on `eps`

  h <- function(x) {

    x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
  }

  covs <- get_all_vars(model$formula, data = data)

  d0 <- d1 <- covs

  # calculate change in covariate from (x+h) to (x-h)

  d0[[variable]] <- d0[[variable]] - h(d0[[variable]])

  d1[[variable]] <- d1[[variable]] + h(d1[[variable]])


  # calculate fitted probabilities at (x+h) and (x-h)

  probs <- fitted_probs(model.matrix(model$formula, data.frame(rbind(d0, d1))), model$coefs)


  effect <- (probs[nrow(d0) + seq_len(nrow(d0)), ] - probs[seq_len(nrow(d0)), ]) / (unlist(d1[[variable]] - d0[[variable]]))

  return(structure(list(effect),
                   names = paste(variable),
                   class = c("data.frame"),
                   row.names = seq_len(nrow(covs))))
}


#help function to calculate marginal effect of factor covariate with 2 levels or more.

marginal_factor <- function(model, data, variable) {

  covs <- get_all_vars(model$formula, data = data)

  levs <- levels(as.factor(covs[[variable]]))

  base <- (levs[1L])

  main <- (levs[-1L])

  out <-  structure(rep(list(list()), (length(levs) - 1)),
                    class = "data.frame",
                    names = paste0(variable, main),
                    row.names = seq_len(nrow(covs)))

  outcolnames <- paste0(variable, main)

  d0 <- d1 <- covs

  # replace all values of covariate with base values

  d0[[variable]] <- factor(base, levels = levs)

  pred0 <- fitted_probs(model.matrix(model$formula, data.frame(d0)), model$coefs)


  # calculate difference between each factor level and base level

  for (i in seq_along(main)) {

    d1[[variable]] <- factor(main[i], levels = levs)

    pred1 <-  fitted_probs(model.matrix(model$formula, data.frame(d1)), model$coefs)

    out[[outcolnames[i]]] <- (pred1 - pred0)

  }

  return(out) # return data.frame with column of derivatives
}


marginal <- function(model, data, eps = 1e-7) {


  # identify classes of terms in `model`


  variables <- get_all_vars(model$formula, data = data)


  classes <- lapply(variables, function(x) paste(class(x), collapse = ','))

  varslist <- list(

    nnames = unique(names(classes)[!classes %in% c("factor", "ordered", "logical")]),

    fnames = unique(names(classes)[classes %in% c("factor", "ordered")])
  )

  # estimate numerical derivatives with respect to each variable

  out1 <- lapply(varslist$nnames, marginal_num, model = model, data = data, eps = eps)

  # add discrete differences for factor terms

  out2 <- list()

  for (i in seq_along(varslist$fnames)) {

    out2[[i]] <- marginal_factor(model = model, data = data, varslist$fnames[i])
  }

  out <- c(out1, out2)

  out <- do.call("cbind", out[vapply(out, function(x) length(x) > 0, FUN.VALUE = logical(1))])


  return(out)
}

# Jacobian function: derivatives of marginal effects of each level w.r.t each coeff
##  as a function of average marginal effect and coefficients


jac <- function(model, data, eps = 1e-7){

  coeff <- model$coefs

  M0 <- colMeans(marginal(model = model, data = data, eps = eps), na.rm = TRUE)

  mat <- matrix(NA, nrow = length(M0), ncol = length(coeff))

  for (i in seq_along(coeff)) {

    coeft <- coeff

    coeft[i] <- coeft[i] + eps

    model[["coefs"]] <- coeft

    mat[, i] <- (colMeans(marginal(model = model, data = data, eps = eps), na.rm = TRUE) - M0) / eps
  }

  mat

}


variance <- function(model, data,  method = c("delta", "simulation"), K = 500, eps = 1e-7){

  method <- match.arg(method)

  var_cov <- model$varmat

  if (method == "delta") {

    jac_fun <- jac(model, data = data, eps = eps)

    var_mat <- jac_fun %*% var_cov %*% t(jac_fun)

    se_mat <- sqrt(diag(var_mat))

  } else if (method == "simulation") {

    out <- matrix(NA, nrow = K, ncol = ncol(marginal(model, data = data, eps = eps)))

    for (i in 1:K) {

      coefficient <- MASS::mvrnorm(K, c(model$coefs), model$varmat)

      new <- model

      new[["coefs"]] <- matrix(coefficient[i, ], ncol = 2)

      means <- colMeans(marginal(model = new, data = data, eps = eps), na.rm = TRUE)

      out[i, ] <- means
    }
    se_mat <- sqrt(diag(var(out)))

  }

  return(se_mat)

}
