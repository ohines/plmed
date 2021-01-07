#' Baron-Kenny style Mediated Effects
#'
#' \code{OLS} is used to estimate mediated effects
#' (Baron and Kenny, 1986) in linear models using Ordinary Least Squares based methods.
#' Linear models are assumed for the conditional expectation of mediator given exposure and confounders, 
#' and for the conditional expectation of outcome given mediator, exposure and confounders. As with the \code{\link{plmed}} function,
#' the confounder set is the union terms in the
#' \code{mediator.formula}, and \code{outcome.formula}. Missing data behaviour is always \code{\link[stats]{na.action}=na.omit}.
#'
#' @param exposure.variable a character string naming the exposure variable of interest.
#' @param mediator.formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class) 
#' where the left hand side of the formula contains the mediator variable of interest.
#' @param outcome.formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class) 
#' where the left hand side of the formula contains the outcome variable of interest.
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{OLS} is called.
#' @return An object of class \code{plmed} with, Natural Direct Effect (NDE) and 
#' Natural Indirect Effect (NIDE) estimates, as well as the effect of exposure on mediator (X_on_M)
#' and the effect of mediator on outcome (M_on_Y). Estimated standard errors, and
#' Wald based test statistics are also returned, as is the Likelihood Ratio (LR) based test statitic.
#' @examples
#' #Example on Generated data
#' N <- 100
#' beta <- c(1,0,1) #Some true parameter values
#' #Generate data on Confounders (Z), Exposure (X)
#' #Mediator (M), Outcome (Y)
#' Z <- rnorm(N)
#' X <- rbinom(N,1,1/(exp(-Z)+1))
#' M <- beta[1]*X + Z +rnorm(N)
#' Y <- beta[2]*M + beta[3]*X + Z +rnorm(N)
#'
#' OLS("X",M~Z,Y~Z)
#' @export
OLS <- function(exposure.variable,mediator.formula,outcome.formula,data,weights){
  cl <- match.call(expand.dots = FALSE)
  mf <- cl
  m <- match(c("mediator.formula","outcome.formula","data","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$method <- "OLS"
  mf$exposure.formula <- reformulate("1",response=exposure.variable)
  mf[[1L]] <- quote(plmed)
  a <- eval(mf, parent.frame())
  a$call <- cl
  return(a)
}

#' @export
fit.ols <- function(Y,M,X,Z,weights=rep(1,N)){
  N = NROW(X)
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  
  M.mod <- glm.fit(cbind(X,Z),M,weights=weights)
  Y.mod <- glm.fit(cbind(M,X,Z),Y,weights=weights)
  beta <- c(M.mod$coefficients[1],Y.mod$coefficients[1:2])
  
  seM <- diag(chol2inv(M.mod$qr$qr))*sum((weights*M.mod$residuals^2)[weights>0])/M.mod$df.residual
  seY <- diag(chol2inv(Y.mod$qr$qr))*sum((weights*Y.mod$residuals^2)[weights>0])/Y.mod$df.residual
  
  T1_ols = beta[1]^2/seM[1]
  T2_ols = beta[2]^2/seY[1]
  T3_ols = beta[3]^2/seY[2]
  
  a <- list()
  Sob <- T1_ols*T2_ols/(T1_ols+T2_ols)
  a$Wald = c(T3_ols, Sob,T1_ols,T2_ols )
  a$coefs = c(beta[3], beta[1]*beta[2] , beta[1], beta[2])
  a$std.err = abs(a$coef/sqrt(a$Wald))
  
  a$score <-  min(T1_ols,T2_ols)
  names(a$score) <- "LR Mediation"
  names(a$Wald) <- names(a$coef) <- names(a$std.err) <- c('NDE','NIDE','X_on_M','M_on_Y')
  (a)
}
