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
plmed <- function(exposure.formula,mediator.formula,outcome.formula,
                     data,weights){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","data","weights"),names(mf), 0L)
  
  vars <- sapply(1:3,function(i){
    x = as.list(mf[m[c(i)]])
    tf = do.call(terms,list(x=x[[1]]))
    if(attr(tf,'response')!=1){
      stop(gettextf("Missing reponse variable in %s",names(mf)[i]), domain = NA)
    }
    varnames = vapply(tf,function(x) {
      paste(deparse(x,width.cutoff = 500L, backtick = !is.symbol(x) && is.language(x)),collapse = " ")}, " ")[-1]
    (varnames)
  })

  PL_formula = reformulate(c(vars[1,],vars[2,]))
  dft = append(list(formula = PL_formula,
                    drop.unused.levels = TRUE,
                    na.action=na.omit),as.list(mf[m[4:5]]))
  df = do.call(stats::model.frame,dft,envir = parent.frame())
  
  weights <- model.weights(df)
  if( !is.null(weights) && any(weights < 0) ){
    stop("negative weights not allowed")
  }
    
  df <- model.matrix(attr(df, "terms"),df)
  
  X = as.numeric(df[,2])
  M = as.numeric(df[,3])
  Y = as.numeric(df[,4])
  Z = cbind(1,scale(df[,-(1:4)]))

  a <- fit.plmed(Y,M,X,Z,compute_CUE=TRUE,weights=weights) 
  a$call <- mf
  class(a) <- "plmed"
  a
}

#' @export
fit.plmed <- function(Y,M,X,Z,compute_CUE=TRUE,weights=rep(1,N)){
  N <- NROW(Y)
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  
  #Get initial parameter estimates for Newton Raphson using MLE
  X.lm <- glm.fit(Z,X,family=quasibinomial(),weights=weights)$coefficients
  M.lm <- glm.fit(cbind(X,Z),M,weights=weights)$coefficients
  Y.lm <- glm.fit(cbind(M,X,Z),Y,weights=weights)$coefficients
  
  beta  <- c(M.lm[1],Y.lm[1:2]) #target parameters
  gam.x <- X.lm #nuisance parameters
  gam.m <- M.lm[-1]
  gam.y <- Y.lm[-(1:2)]
  theta <- c(beta,gam.x,gam.x,gam.m,gam.m,gam.y,gam.y)
  
  #Do unconstrained fit
  fit.unconstr <-   newton_raph(CUE_vec_J_bin,theta,X=X,M=M,Y=Y,Z=Z,method='G',weights=weights)
  theta.unc = fit.unconstr$par
  
  a <- list()
  a$coef <- c(fit.unconstr$par[1:3],fit.unconstr$par[1]*fit.unconstr$par[2])
  a$std.err      <- sqrt(c(fit.unconstr$val$var,fit.unconstr$par[1]^2*fit.unconstr$val$var[2] +
                             fit.unconstr$par[2]^2*fit.unconstr$val$var[1])  )
  a$Wald         <- c(fit.unconstr$val$T_stats,fit.unconstr$val$RSTest)
  names(a$Wald) <- c('beta1','beta2','NDE','NIDE')
  
  if(compute_CUE){
    w = c(rep.int(1,length(theta.unc)),length(theta)) #upweight solving the lagrange multiplier = faster
    
    fit.H0.cue <- tryCatch({
      newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,method='CUE',weights=weights,med_prop=0,
                  LSearch = FALSE, w=w,Max.it=50)
    },error = function(e){
      tryCatch({
        newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,method='CUE',weights=weights,med_prop=0,
                    LSearch = TRUE, w=w)
      },error = function(e){
        stop(gettextf("Numerical Error in CUE Calculation."), domain = NA)})})
    
    fit.H1.cue <- tryCatch({
      newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,method='CUE',weights=weights,med_prop=1,
                  LSearch = FALSE, w=w,Max.it=50)
    },error = function(e){
      tryCatch({
        newton_raph(CUE_vec_J_bin,c(theta.unc,0),X=X,M=M,Y=Y,Z=Z,method='CUE',weights=weights,med_prop=1,
                    LSearch = TRUE, w=w)
      },error = function(e){
        stop(gettextf("Numerical Error in CUE Calculation."), domain = NA)})})
    
    a$Score.cue = c(fit.H1.cue$val$score,fit.H0.cue$val$score)
    
  }

  return(a)
  
}



#' @export
print.plmed <- function(object){
  a <- object

  cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Coefficients:\n")

  df <- data.frame(Estimate = a$coef, Std.Error = a$std.err ,
                   Wald.value = a$Wald,
                   Wald.pval = pchisq(a$Wald,df=1,lower.tail = F) )

  printCoefmat(df, digits = 6, signif.stars = T,na.print = "NA",
               tst.ind = 3,P.values=T,has.Pvalue=T)

  print_scores = !is.null(a$Score.cue)
  if(print_scores){
    cat('\nCUE Score Test No-Direct-Effect:',formatC(a$Score.cue[1], digits = 6),'with p-value:',
        format.pval(pchisq(a$Score.cue[1],df=1,lower.tail = F),
                    digits = 6))
    
    cat('\nCUE Score Test No-Mediation:\t',formatC(a$Score.cue[2], digits = 6),'with p-value:',
        format.pval(pchisq(a$Score.cue[2],df=1,lower.tail = F),
                    digits = 6))
    
  }

}

