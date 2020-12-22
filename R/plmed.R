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
#' a family function or the result of a call to a family function. (See \link[stats]{family} for details of family functions.)
#' Must be \code{"binomial"} when using \code{Method="TTS"}. Must be \code{"none"} when using \code{Method="OLS"}.
#' @param mediator.family link function for the mediator model, can be either \code{"gaussian"} or \code{"binomial"}. 
#' Must be \code{"gaussian"} when using \code{Method="G-estimation"}
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{plmed} is called.
#' @param weights an optional vector of ‘observation weights’ to be used in 
#' the fitting process. Should be \code{NULL} or a numeric vector.
#' @param method The mediation fitting method to be used. Can be either \code{"G-estimation","TTS"} or \code{"OLS"}
#'
#' @return An object of class \code{plmed} with unconstrained parameter estimates,
#'  estimated standard errors, Wald based and CUE score based test statistics (G-estimation only).
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
#' plmed(X~Z,M~Z,Y~Z,method="G-estimation")
#' plmed(X~Z,M~Z,Y~Z,method="TTS")
#' plmed(X~Z,M~Z,Y~Z,method="OLS")
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
                  exposure.family="binomial",mediator.family="gaussian",
                  data,weights,method="G-estimation"){
  cl <- match.call(expand.dots = FALSE)
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","data","weights","method"),names(cl), 0L)
  
  # Process formulas in a similar way to stats::glm
  
  xf <- cl[c(1,m[c(1,4,5)])]
  mf <- cl[c(1,m[c(2,4,5)])]
  yf <- cl[c(1,m[c(3,4,5)])]
  xf[[1L]] <- mf[[1L]]  <- yf[[1L]] <- quote(stats::model.frame)
  xf$drop.unused.levels <- mf$drop.unused.levels <- yf$drop.unused.levels <- TRUE
  xf$na.action <- mf$na.action <- yf$na.action <- na.pass
  names(xf)[2] <- names(mf)[2] <- names(yf)[2] <- "formula"
  
  xf <- eval(xf, parent.frame())
  mf <- eval(mf, parent.frame())
  yf <- eval(yf, parent.frame())
  
  X <- model.response(xf, "any")
  M <- model.response(mf, "any")
  Y <- model.response(yf, "any")
  
  if(is.null(X)|is.null(M)|is.null(Y)){
    stop("All formulas must have a response variable")
  }
  
  xt <- attr(xf, "terms") #allow model.frame to have updated it
  mt <- attr(mf, "terms")
  yt <- attr(yf, "terms") 
  
  
  vars <- lapply(list(xt,mt,yt),function(tf){
    varnames = vapply(tf,function(x) {
      paste(deparse(x,width.cutoff = 500L, backtick = !is.symbol(x) && is.language(x)),collapse = " ")}, " ")[-1]
    (varnames)
  })
  vars <- simplify2array(vars)
  
  dft = append(list(formula = reformulate(c(vars[2,])),
                    drop.unused.levels = TRUE,
                    na.action=na.pass),as.list(cl[m[4:5]]))
  df = do.call(stats::model.frame,dft,envir = parent.frame())
  weights <- model.weights(df)
  
  df <- model.matrix(attr(df, "terms"),df)
  Z = cbind(1,scale(df[,-1]))
  
  keeps <- (complete.cases(X)&
              complete.cases(Y)&
              complete.cases(M)&
              complete.cases(Z))
  X = X[keeps]
  M = M[keeps]
  Y = Y[keeps]
  Z = Z[keeps,]
  
  if( !is.null(weights) && any(weights < 0) ){
    stop("negative weights not allowed")
  }else{
    weights <- weights[keeps]
  }
  
  
  allowed_families <- as.environment(list(gaussian = stats::gaussian,binomial = stats::binomial))
  ## exposure family
  if (method=="G-estimation"){
    if(is.character(exposure.family))
      Xfam <- get(exposure.family, mode = "function", envir = parent.frame())
    if(is.function(exposure.family)) Xfam <- exposure.family()
    if(is.null(Xfam()$family)) {
      stop("'exposure.family' not recognized")
    }
    if(mediator.family != "gaussian"){
      warning("G-estimation methods require a binary exposure.\n  Using 'mediator.family' = 'gaussian'")
    }
    Mfam <- allowed_families$gaussian
    a <- fit.G_estimation(Y,M,X,Z,Xfam(),compute_CUE=TRUE,weights=weights) 
    
  }  else if (method=="TTS"){
    ## mediator family
    if(is.character(mediator.family)) {
      Mfam <- tryCatch({
        get(mediator.family,mode = "function",
            envir = allowed_families)
      },error=function(e){
        e$message = gettextf("'mediator.family' must be either 'gaussian' or 'binomial' not: %s",mediator.family)
        stop(e)
      })
    }else{
      stop("'mediator.family' must be either 'gaussian' or 'binomial'")
    }
    if(exposure.family != "binomial"){
      warning("TTS methods require a binary exposure.\n  Using 'exposure.family' = 'binomial'")
    }
    Xfam <- allowed_families$binomial
    a <- fit.TTS(Y,M,X,Z,Mfam(),weights=weights) 
  } else if(method == "OLS"){
    Xfam <- function() list(family="none")
    Mfam <- allowed_families$gaussian
    a <- fit.ols(Y,M,X,Z,weights=weights) 
  }else{
    stop("Method not recognized")
  }
  
  
  a$Method <- method
  a$call <- cl
  a$exposure.family = Xfam()$family
  a$mediator.family = Mfam()$family
  a$outcome.family =  "gaussian"
  class(a) <- "plmed"
  a
}


#' @export
print.plmed <- function(object){
  a <- object

  cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Fitting Method: \t",a$Method,
      "\nExposure GLM family:\t",a$exposure.family)
  if (!is.null(a$mediator.family)){
    cat("\nMediator GLM family:\t",a$mediator.family)
  }
  cat("\nOutcome  GLM family:\t",a$outcome.family)

  cat("\n\nCoefficients:\n")

  df <- data.frame(Estimate = a$coef, Std.Error = a$std.err ,
                   Wald.value = a$Wald,
                   Wald.pval = pchisq(a$Wald,df=1,lower.tail = F) )

  printCoefmat(df, digits = 6, signif.stars = T,na.print = "NA",
               tst.ind = 3,P.values=T,has.Pvalue=T)

  print_scores = !is.null(a$score)
  if(print_scores){
    for (sc in (1:length(a$score)) ){
      cat('\n',names(a$score)[sc],'Score:',formatC(a$score[sc], digits = 6),'with p-value:',
          format.pval(pchisq(a$score[sc],df=1,lower.tail = F),
                      digits = 6))
      
    }
  }
  
}

