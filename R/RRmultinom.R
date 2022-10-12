#' Prevalence estimation for the bivariate model of the "ever"
#' and "last year" RR Questions
#'
#' @description Point and interval estimates of the
#' Maximum Likelihood Estimator (MLE) for prevalence
#' of never, former and recent categories through the bivariate model of the "ever"
#' and "last year" randomized response (RR) questions.
#' @param Q1 name of the question "Have you ever...?".
#' @param Q2 name of the question "In the last year, Did you ...?"
#' @param data data frame with individual records containing \code{Q1} and
#' \code{Q2} variables with values 1 denoting "yes" and 0 for "no" answer.
#' @param p00 The conditional probability of observed "no" while
#' true response is "no".
#' @param p11 The conditional probability of observed "yes" while
#' true response is "yes".
#' @param seed a single integer value. Default: \code{seed = NULL}. If you got \code{NAN}, you should specify
#' seeds.
#'
#' @return A list with the elements:
#' \describe{
#' \item{est}{data frame with point and interval estimates for the MLE.}
#' \item{gof}{Goodness-of-fit statistic for the bivariate model.}
#' \item{freqs}{data frame with the observed and estimated response frequencies.}
#' }
#'
#' @examples
#' prev_bivar(Q1 = Q9, Q2 = Q10, data = Doping, p00 = 5/6, p11 = 5/6)
#'
#' @importFrom stats optim pchisq xtabs
#' @importFrom dplyr arrange %>%
#' @importFrom numDeriv jacobian
#'
#' @export
prev_bivar <- function(Q1,  Q2, data, p00, p11, seed = NULL){

  ever   <- data[[substitute(Q1)]]
  past12 <- data[[substitute(Q2)]]

  if (!is.numeric(ever))  stop("Q1 variable is not numeric.")

  if (!is.numeric(past12)) stop("Q2 variable is not of class numeric.")

  if (sum(is.na(ever)) > 0 || sum(is.na(past12)) > 0)     stop("Q1 variable and|or Q2 have missing.")

  if (!all(ever %in% 0:1))    stop("Q1 variable has values other than 0 and 1.")

  if (!all(past12 %in% 0:1))  stop("Q2 variable has values other than 0 and 1.")


  R <- data.frame(xtabs(~ ever + past12)) %>%  arrange(ever, past12)

  obs  <- R$Freq

  P  <- matrix(c(p00, 1 - p11, 1 - p00, p11), 2, 2)

  Q0 <- (P %x% P)[, -2]

  if (!is.null(seed)) set.seed(seed)

  mle <- optim(rnorm(2, 0, .1), fn = log_prev, Q = Q0, obs = obs)

  mle <- optim(mle$par, fn = log_prev, Q = Q0,  obs = obs, method = "BFGS", hessian = T)

  estimates <- softmax(mle$par)

  ## derivatives using the existing function jacobian

  jaco  <- numDeriv::jacobian(softmax, x = mle$par)[-1, , drop = F] # derivatives w.r.t. parameters

  vpar <- t(jaco) %*% solve(mle$hessian) %*% jaco                    # variance probabilities former and recent

  se   <- sqrt(c(sum(vpar), diag(vpar)))                           # standard deviations including reference category

  prevalences  <- data.frame(est       = softmax(mle$par),

                             se        = se,

                             min95     = ifelse(estimates - 1.96 * se < 0, 0, estimates - 1.96 * se),

                             max95     = ifelse(estimates + 1.96 * se > 1, 1, estimates + 1.96 * se),

                             row.names = c("Never", "Former", "Recent"))

  fitted <- sum(obs) * Q0 %*% estimates

  gof <- data.frame(G2        =  2 * t(obs) %*% log(obs / fitted),
                    df        =  1,
                    pvalue    =  pchisq(2 * t(obs) %*% log(obs / fitted), df = 1, lower.tail = F),
                    row.names = "G2 statistic")

  freqs <- data.frame(response = c("nn", "ny", "yn", "yy"), obs, fitted)


  colnames(freqs) = c("response", "observed freq.", "fitted")

  cat("Log likelihood", paste(round(-mle$value, 1)),"\n", "\n")
  cat("prevalence estimates \n \n")
  print(prevalences %>% round(3))
  cat("\n")
  cat("Goodness-of-fit \n \n")
  print(gof %>% round(4))
  cat("\n")
  cat("observed and estimated frequencies \n \n")
  print(freqs)

  invisible(list(est         = round(prevalences, 3),
                 gof         = round(gof, 4),
                 freqs       = freqs,
                 Hessian     = mle$hessian))

}

#' Multinomial logistic RR regression model
#' @description Fits a multinomial logistic regression model to the bivariate analysis of the "ever"
#' and "last year" Questions to investigate potential covariates associated
#' with the sensitive characteristic.
#'
#' @param formula a "\code{\link{formula}}" in the form
#' \code{~ terms}, where  \code{terms} describes the linear predictors.
#' @param data data frame with individual records containing the
#' \code{response} with values 1 for "yes", 0 for "no", and the predictors.
#' @param response character vector with the names of the response variables of the
#' two questions "Have you ever...?" and "in the last year, Did you ...?", in that order.
#' @param p00 The conditional probability of observed "no" while
#' true response is "no".
#' @param p11 The conditional probability of observed "yes" while
#' true response is "yes".
#' @param iter The maximum number of iterations of optimization method "Nelder-Mead". Default is 500.
#'
#' @return A list with the following objects:
#' \describe{
#' \item{formula}{formula used in the model specification}
#' \item{Call}{the matched call}
#' \item{coefs}{matrix of the regression coefficients}
#' \item{parmsig}{data frame with the parameter estimates and its signficance.}
#' \item{gof}{fit statistics including AIC and a likelihood-ratio
#' test against the intercept-only model.}
#' \item{fitted}{data frame with the fitted probabilities.}
#' \item{varmat}{variance-covariance matrix of the parameter estimates.}
#' \item{modeldata}{data frame with the covariates used in \code{formula}}
#' \item{D}{design matrix of the fitted model.}
#' }
#'
#' @examples
#' regression <- rr_multinom(formula = ~BodyBuild+Age, data = Doping,
#' response = c("Q9", "Q10"), p00 = 5/6, p11 = 5/6)
#'
#' @importFrom stats model.matrix optim pchisq pnorm rnorm model.frame get_all_vars
#' @importFrom dplyr mutate case_when %>%
#'
#' @export

rr_multinom <- function(formula, data, response, p00, p11, iter = 500){

  if (!is.data.frame(data))   stop("The input data must be data frame")

  D       <- model.matrix(formula, data)

  depend  <- data[response]

  ever    <- depend[1]

  past12  <- depend[2]

  R       <- mutate(depend, case_when(ever == 0 & past12 == 0 ~ 1,
                                      ever == 0 & past12 == 1 ~ 2,
                                      ever == 1 & past12 == 0 ~ 3,
                                      ever == 1 & past12 == 1 ~ 4),
                    .keep = "none") %>%
    unlist() %>%
    unname()

  P  <- matrix(c(p00, 1 - p00, 1 - p11, p11), 2, 2)

  Q0 <- (P %x% P)[, -2]

  P_trans  <- Q0[R, ]

  reg <- optim(matrix(rep(0, 2 * ncol(D)), ncol = 2, nrow = ncol(D)), logl_reg, P_trans = P_trans, D = D, control = list(maxit = iter))

  reg <- optim(reg$par, logl_reg, P_trans = P_trans, D = D, method = "BFGS", hessian = TRUE)


  se <- tryCatch(c(sqrt(diag(solve(reg$hessian)))),

                 warning = function(w){
                   message("Unreliable solution due a non-invertible Hessian")
                   NA})

  est <- data.frame(Coef.    = c(reg$par),

                    SE.      = c(sqrt(diag(solve(reg$hessian)))),

                    z.value  = c(reg$par)/c(sqrt(diag(solve(reg$hessian)))),

                    p.value  = 2 * pnorm(-abs(c(reg$par)/c(sqrt(diag(solve(reg$hessian)))))),

                    row.names = c(paste0("Former: ", colnames(D)), paste0("Recent: ", colnames(D)))

  )


  prev <- optim(matrix(0, ncol = 2), logl_reg, P_trans = P_trans, D =  model.matrix(~1, as.data.frame(data)), control = list(maxit = iter))


  LR <- 2 * (prev$value - reg$value)

  D       <- model.matrix(formula, as.data.frame(data))

  gof <-  data.frame(AIC     = 2 * (reg$value + ncol(D) * 2) %>%  round(1),
                     LR      = LR %>%  round(2),
                     df      = (ncol(D) - 1) * 2,
                     p       = pchisq(LR, df = (ncol(D) - 1) * 2, lower.tail = FALSE) %>%  round(4),
                     row.names = "")


  probs           <-  apply(D %*% reg$par, 1, softmax)

  rownames(probs) <- c("Never", "Former", "Recent")

  cat(" ", "Call: ", paste(call("rr_multinom", formula)), "\n", "\n",

      "Log likelihood: ", paste(round(-reg$value, 2)),"\n", "\n")

  cat("Parameter estimates:", "\n")
  print(est %>% round(3))
  cat("\n")
  cat("\n")
  cat("Fit statistics")
  cat("\n")
  print(gof)
  cat("\n")

  if (reg$convergence != 0) warning("The model did not successfully converge")

  invisible(list(formula      = formula,
                 Call         = paste(call("rr_multinom", formula)),
                 modeldata     = get_all_vars(formula, data),
                 coefs         = reg$par,
                 gof          = gof,
                 fitted       = data.frame(t(probs)),
                 varmat       = solve(reg$hessian),
                 modelmatrix  = D

  )
  )


}


#' @title Average Marginal effect
#' @description calculates average marginal effect (AME) along with its standard error
#' for all numerical and categorical covariates of the model object \code{\link{rr_multinom}}.
#' @param model the model object of class \code{\link{rr_multinom}}.
#' @param data The data set on which to calculate the marginal effect.
#' @param method character specifying the method of estimation: delta or simulation.
#' @param K integer number indicating the number of iterations if method "simulation" is specified.
#' @param eps the value of the step used in calculating the marginal effect of numerical covariates.
#' @details For numerical covariates, marginal effect is calculated using two sided numerical differentiation
#' adopted from \url{https://github.com/leeper/margins}.
#'
#' @return data frame with each row is an observations and columns are the marginal effect for all covariates on the estimated probabilities
#' of never, former and recent categories.
#'
#' @references
#' Leeper TJ (2021). margins: Marginal Effects for Model Objects. R package version 0.3.26.
#'
#' @examples
#'regression <- rr_multinom(formula = ~BodyBuild+Age, data = Doping,
#'                          response = c("Q9", "Q10"), p00 = 5/6, p11 = 5/6)
#' ME <- AME(model = regression, data = subset(Doping, BodyBuild == 0),
#'           method = "delta")
#' df <- ME$Age
#' @importFrom stats pnorm var
#' @importFrom MASS mvrnorm
#' @export

AME <- function(model, data, method = c("delta", "simulation"), K = 500, eps = 1e-7){

  method <- match.arg(method)

  Std <- variance(model, data,method = method, K = K, eps = eps)

  ME <- data.frame(AME      =  c(colMeans(marginal(model, data,eps = eps), na.rm = TRUE)),
                   SE.      =  c(Std),
                   min95    =  c(colMeans(marginal(model, data, eps = eps), na.rm = TRUE)) - 1.96 * c(Std),
                   max95    =  c(colMeans(marginal(model, data, eps = eps), na.rm = TRUE)) + 1.96 * c(Std),
                   p.value  =  2 * pnorm(-abs(c(colMeans(marginal(model, data, eps = eps), na.rm = TRUE))/c(Std))),
                   row.names = names(colMeans(marginal(model, data, eps = eps))))
  print(ME %>% round(3))

  invisible(data.frame(marginal(model, data, eps = eps)))
}




