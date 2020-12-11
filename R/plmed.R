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
#' @param exposure.family link function for the exposure model, can be either \code{"gaussian"} or \code{"binomial"}
#' Must be \code{"binomial"} when using \code{Method="TTS"}
#' @param mediator.family link function for the mediator model, can be either \code{"gaussian"} or \code{"binomial"}. 
#' Must be \code{"gaussian"} when using \code{Method="G-estimation"}
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{plmed} is called.
#' @param weights an optional vector of ‘observation weights’ to be used in 
#' the fitting process. Should be \code{NULL} or a numeric vector.
#' @param method The mediation fitting method to be used. Can be either \code{"G-estimation"} or \code{"TTS"}
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
                  exposure.family="binomial",mediator.family="gaussian",
                  data,weights,method="G-estimation"){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","data","weights","method"),names(mf), 0L)
  
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
  #Z = cbind(1,scale(df[,-(1:4)]))
  Z = cbind(1,df[,-(1:4)])
  
  allowed_families <- as.environment(list(gaussian = stats::gaussian,binomial = stats::binomial))
  ## exposure family
  if (method=="G-estimation"){
    if(is.character(exposure.family)) {
      Xfam <- tryCatch({
        get(exposure.family,mode = "function",
            envir = allowed_families)
      },error=function(e){
        e$message = gettextf("'exposure.family' must be either 'gaussian' or 'binomial' not: %s",mediator.family)
        stop(e)
      })
    }else{
      stop("'exposure.family' must be either 'gaussian' or 'binomial'")
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
  }else{
    stop("Method not recognized")
  }
  

  a$Method <- method
  a$call <- mf
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
      "\nExposure GLM family:\t",a$exposure.family,
      "\nMediator GLM family:\t",a$mediator.family,
      "\nOutcome  GLM family:\t",a$outcome.family)

  cat("\n\nCoefficients:\n")

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

