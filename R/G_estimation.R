#' Partially Linear Mediation using G-estimation
#'
#' \code{plmed} is used to fit mediation effects for a continuous outcome
#' given a single continuous mediator, binary exposure, and set of
#' variables which adjust for confounding. It is assumed that the
#' outcome and mediator models are linear and the exposure model is
#' logisitic. It uses a G-estimation procedure to estimate indirect and
#' direct effect, with a Bias-Reduced strategy used to estimate
#' parameters in confounder models. Missing data behaviour is always \code{\link[stats]{na.action}=na.omit}
#'
#' @param exposure.formula an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to
#' that class) where the left hand side of the formula contains the binary exposure variable of interest.
#' @param mediator.formula an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to
#' that class) where the left hand side of the formula contains the continuous mediator variable of interest.
#' @param outcome.formula an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to
#' that class) where the left hand side of the formula contains the continuous outcome variable of interest.
#' @param exposure.family link function for the exposure model, can be can be a character string naming a family function,
#'  a family function or the result of a call to a family function. (See \link[stats]{family} for details of family functions.)
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{plmed} is called.
#' @param weights an optional vector of ‘prior weights’ to be used in 
#' the fitting process. Should be \code{NULL} or a numeric vector.
#'
#' @return An object of class \code{plmed} with unconstrained parameter estimates,
#'  estimated standard errors, Wald based and CUE score based test statistics.
#' @examples
#' #Example on Generated data
#' N <- 100
#' beta <- c(1,0,1) #Some true parameter values
#' #Generate data on Confounders (Z), Exposure (X)
#' #Mediator (M), Outcome (Y)
#' Z <- rnorm(N)
#' X <- rbinom(N,1,plogis(Z))
#' M <- beta[1]*X + Z +rnorm(N)
#' Y <- beta[2]*M + beta[3]*X + Z +rnorm(N)
#'
#' plmed(X~Z,M~Z,Y~Z)
#'
#'
#' #Example on JobsII data from the mediation package
#' jobs <- mediation::jobs
#'
#' Z.formula = c('econ_hard','sex','age','occp',
#'               'marital','nonwhite','educ','income')
#' plmed(reformulate(Z.formula,response='treat'),
#'       reformulate(Z.formula,response='job_seek'),
#'       reformulate(Z.formula,response='depress2'),
#'       data=jobs)
#' @export
G_estimation <- function(exposure.formula,mediator.formula,outcome.formula,exposure.family="gaussian",data,weights){
  cl <- match.call(expand.dots = FALSE)
  mf <- cl
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","exposure.family","data","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$method <- "G-estimation"
  mf[[1L]] <- quote(plmed)
  a <- eval(mf, parent.frame())
  a$call <- cl
  return(a)
}


#' @export
fit.G_estimation <- function(Y,M,X,Z,Xfam=quasibinomial(),compute_CUE=TRUE,weights=rep(1,N)){
  N <- NROW(Y)
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  
  #Get initial parameter estimates for Newton Raphson using MLE
  X.lm <- glm.fit(Z,X,family=Xfam,weights=weights)$coefficients
  M.lm <- glm.fit(cbind(X,Z),M,weights=weights)$coefficients
  Y.lm <- glm.fit(cbind(M,X,Z),Y,weights=weights)$coefficients
  
  beta  <- c(M.lm[1],Y.lm[1:2]) #target parameters
  gam.x <- X.lm #nuisance parameters
  gam.m <- M.lm[-1]
  gam.y <- Y.lm[-(1:2)]
  theta <- c(beta,gam.x,gam.x,gam.m,gam.m,gam.y,gam.y)
  
  #Do unconstrained fit
  fit.unconstr <-   newton_raph(CUE_vec_J_bin,theta,X=X,M=M,Y=Y,Z=Z,Xfam=Xfam,method='G',weights=weights)
  theta.unc = fit.unconstr$par
  
  a <- list()
  a$coef <- c(fit.unconstr$par[1:3],fit.unconstr$par[1]*fit.unconstr$par[2])
  a$std.err      <- sqrt(c(fit.unconstr$val$var,fit.unconstr$par[1]^2*fit.unconstr$val$var[2] +
                             fit.unconstr$par[2]^2*fit.unconstr$val$var[1])  )
  a$Wald         <- c(fit.unconstr$val$T_stats,fit.unconstr$val$RSTest)
  names(a$Wald) <- names(a$coef) <- names(a$std.err) <- c('X_on_M','M_on_Y','NDE','NIDE')
  a$Wald <- a$Wald[c(3,4,1,2)]
  a$coef <- a$coef[c(3,4,1,2)]
  a$std.err <- a$std.err[c(3,4,1,2)]
  
  if(compute_CUE){
    w = c(rep.int(1,length(theta.unc)),length(theta)) #upweight solving the lagrange multiplier = faster
    
    fit.H0.cue <- tryCatch({
      newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,Xfam=Xfam,method='CUE',weights=weights,med_prop=0,
                  LSearch = FALSE, w=w,Max.it=50)
    },error = function(e){
      tryCatch({
        newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,Xfam=Xfam,method='CUE',weights=weights,med_prop=0,
                    LSearch = TRUE, w=w)
      },error = function(e){
        stop(gettextf("Numerical Error in CUE Calculation."), domain = NA)})})
    
    fit.H1.cue <- tryCatch({
      newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,Xfam=Xfam,method='CUE',weights=weights,med_prop=1,
                  LSearch = FALSE, w=w,Max.it=50)
    },error = function(e){
      tryCatch({
        newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,Xfam=Xfam,method='CUE',weights=weights,med_prop=1,
                    LSearch = TRUE, w=w)
      },error = function(e){
        stop(gettextf("Numerical Error in CUE Calculation."), domain = NA)})})
    
    a$score = c(CUE_0 = fit.H0.cue$val$score,
                CUE_1 = fit.H1.cue$val$score)
    names(a$score) <- c("CUE Mediation", "CUE Direct-Effect")
  }
  
  return(a)
  
}