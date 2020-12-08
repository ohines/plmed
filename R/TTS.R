#' Estimating Mediated Effects
#'
#' \code{TTS} is used to estimate natural direct and natural indirect effects using the methods of 
#' Tchetgen Tchetgen Shpitser 2012. This implementation can be used with a binary exposure, binary or continuous
#' mediator, and continuous outcome. Linear models are used to model the conditional expectation of continuous quantities,
#' and logisitc regression models are used to model binary quantities.
#'
#' @param exposure.formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class) 
#' where the left hand side of the formula contains the binary exposure variable of interest.
#' @param mediator.formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class) 
#' where the left hand side of the formula contains the mediator variable of interest.
#' @param outcome.formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class) 
#' where the left hand side of the formula contains the outcome variable of interest.
#' @param mediator.family link function for the mediator model, can be either \code{"gaussian"} or \code{"binomial"}
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{TTS} is called.
#' @return An object of class \code{TTS} with Total Effect (TE), Natural Direct Effect (NDE) and 
#' Natural Indirect Effect (NIDE) estimates, with estimated standard errors, and
#' Wald based test statistics.
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
#' TTS(X~Z,M~Z,Y~Z,"gaussian")
#' 
#' #Example on Generated data
#' N <- 100
#' beta <- c(1,0,1) #Some true parameter values
#' #Generate data on Confounders (Z), Exposure (X)
#' #Mediator (M), Outcome (Y)
#' Z <- rnorm(N)
#' X <- rbinom(N,1,plogis(Z))
#' M <- rbinom(N,1,plogis(X+Z))
#' Y <- beta[2]*M + beta[3]*X + Z +rnorm(N)
#'
#' TTS(X~Z,M~Z,Y~Z,"binomial")
#' @export

TTS <- function(exposure.formula,mediator.formula,outcome.formula,mediator.family="gaussian",data){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","mediator.family","data"),names(mf), 0L)
  
  ## mediator family
  if(is.character(mediator.family)) {
    allowed_families <- as.environment(list(gaussian = stats::gaussian,binomial = stats::binomial))
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
    
  #if(is.function(mediator.family)) mediator.family <- mediator.family()
  #if(is.null(family$family)) {
  #  print(family)
  #  stop("'mediator.family' not recognized")
  #}
  
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
  dft = append(list(object = PL_formula,
                    drop.unused.levels = TRUE,
                    na.action=na.omit),as.list(mf[m[4]]))
  df = do.call(stats::model.matrix,dft,envir = parent.frame())
  
  X = as.numeric(df[,2])
  M = as.numeric(df[,3])
  Y = as.numeric(df[,4])
  Z = cbind(1,scale(df[,-(1:4)]))

  a <- fit.TTS(Y,M,X,Z,Mfam) 
  a$call <- mf
  class(a) <- "TTS"
  a
}



#' @export
fit.TTS <- function(Y,M,X,Z,Mfam){

  m_linkinv <- Mfam()$linkinv
  m_dev.res <- Mfam()$dev.resids
  binM <- Mfam()$family == "binomial"
  contM <- Mfam()$family == "gaussian"
  if(!(binM|contM)) stop('family not recognized')

  Zx=Z
  Zm = cbind(X,Z)
  Zy = cbind(M,X,Z)
  
  X.mod <- glm.fit(Z,X,family=binomial())
  M.mod <- glm.fit(Zm,M,family=Mfam())
  Y.mod <- glm.fit(Zy,Y)
  
  X.lm <- X.mod$coefficients
  M.lm <- M.mod$coefficients
  Y.lm <- Y.mod$coefficients
  beta  <- c(M.lm[1],Y.lm[1:2]) #target parameters
  N = length(X)
  M.var <- ifelse(binM,0,M.mod$deviance/N)  #only used for continuous mediators
  
  PS = X.mod$fitted.values #E(X|Z_i)
  fZ = M.mod$linear.predictors -beta[1]*X
  m0 = m_linkinv(fZ) #E(M|X=0,Z)
  m1 = m_linkinv(fZ+beta[1]) #E(M|X=1,Z)
  y0 = Y.mod$fitted.values-beta[3]*X  #E(Y|M_i,X=0,Z_i)  
  y1 = y0 + beta[3] #E(Y|M_i,X=1,Z_i) 
  M.res0 = M- m0 #M_i-E(M|X=0,Z_i)
  eta_10 = y1 - beta[2]*M.res0
  eta_00 = y0 - beta[2]*M.res0
  eta_11 = y1 - beta[2]*(M- m1)
  
  X_PS = X/PS #comes up a lot - calculate once for speed
  cX_PS = (1-X)/(1-PS)
  
  K    <-   ifelse(binM, m0/m1    , exp((beta[1]/2 -M.res0)*beta[1]/M.var)  )  #f(M_i|X=0,Z_i)/f(M_i|X=1,Z_i)
  D_m0 <-   ifelse(binM, m0*(1-m0) ,1)                      #derivative of m0  wrt. f(z)
  D_m1 <-   ifelse(binM, m1*(1-m1) ,1)                      #derivative of m1  wrt. f(z) (and also beta1)
  D_K  <- K*ifelse(binM, m1-m0    ,beta[1]/M.var)          #derivative of K   wrt. f(z)
  D_Kb1<- K*ifelse(binM, m1-1     ,(beta[1]-M.res0)/M.var) #derivative of K   wrt. beta1
  
  Po0 = cX_PS*(Y- eta_00 ) + eta_00
  Po1 = X_PS*(Y- eta_11) + eta_11
  Po1M0 = X_PS*K*(Y-y1) + cX_PS*(y1-eta_10) + eta_10
  
  h = 1 #these methods were buitl to accomodate other weighting strategies (e.g. TTS 2014)
  h.mean <- function(a) mean(h*a)
  '%h%' <- function(a,b) (a*h)%*%b
  
  Y0 = h.mean(Po0)
  Y1 = h.mean(Po1)
  Y1M0 = h.mean(Po1M0)
  
  ### Expected score function derivatives ###
  #gamma_x (Through the Propensity score)
  bar_gamx_Po0 = (cX_PS*(Y- eta_00 )*PS)%h%Z/N
  bar_gamx_Po1 = -(X_PS*(Y- eta_11)*(1-PS))%h%Z/N
  bar_gamx_Po1M0 = -(X_PS*K*(Y-y1)*(1-PS))%h%Z/N + (cX_PS*(y1-eta_10)*PS)%h%Z/N
  
  #gamma_m (Through eta10,00,11 and K)    
  bar_gamm_Po0 = ((1-cX_PS)*beta[2]*D_m0 )%h%Z/N
  bar_gamm_Po1 = ((1-X_PS)*beta[2]*D_m1)%h%Z/N
  bar_gamm_Po1M0 = (X_PS*(Y-y1)*D_K  + (1-cX_PS)*beta[2]*D_m0)%h%Z/N
  
  #gamma_y (Through eta10,00,11 and y1) 
  bar_gamy_Po0 = (1 - cX_PS)%h%Z/N
  bar_gamy_Po1 = (1 - X_PS)%h%Z/N
  bar_gamy_Po1M0 = (1 - X_PS*K)%h%Z/N
  
  #sigma_M (Through K) #only for continuous mediators
  bar_sig_Po1M0 = ifelse(binM,0, -h.mean(K*X_PS*(Y-y1)*(beta[1]/2 - M.res0))*beta[1]/M.var^2)
  
  #beta1 (Through K, eta11)  #These change for logit vs linear
  bar_beta1_Po1 = h.mean((1 - X_PS)*beta[2]*D_m1 )
  bar_beta1_Po1M0 = h.mean( X_PS*(Y-y1)*D_Kb1 )
  
  #beta2 (through eta10,00,11 and y1)
  bar_beta2_Po0 = h.mean( (1 - cX_PS)*m0 )
  bar_beta2_Po1 = h.mean( (1 - X_PS)*m1 )
  bar_beta2_Po1M0 = h.mean( -X_PS*K*M + cX_PS*M.res0 + m0 )
  
  #beta3 (through y1, eta11,10)
  bar_beta3_Po1 = h.mean(1 - X_PS)
  bar_beta3_Po1M0 = h.mean(1 - X_PS*K)
  
  #Parameter Influence functions
  IFx = IF.glm(X.mod,Zx)
  IFm = IF.glm(M.mod,Zm)
  IFy = IF.glm(Y.mod,Zy)
  IFsig = ifelse(binM, 0, (M.mod$residuals)^2 - M.var)
  
  IF_Po0 = h*Po0 - Y0 + 
    IFx%*%t(bar_gamx_Po0) +
    IFm[,-1]%*%t(bar_gamm_Po0) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po0) +
    bar_beta2_Po0 * IFy[,1] 
  
  IF_Po1 = h*Po1 - Y1 + 
    IFx%*%t(bar_gamx_Po1) +
    IFm[,-1]%*%t(bar_gamm_Po1) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po1) +
    bar_beta1_Po1 * IFm[,1] +
    bar_beta2_Po1 * IFy[,1] +
    bar_beta3_Po1 * IFy[,2]
  
  IF_Po1M0 = h*Po1M0 - Y1M0 + 
    IFx%*%t(bar_gamx_Po1M0) +
    IFm[,-1]%*%t(bar_gamm_Po1M0) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po1M0) +
    bar_beta1_Po1M0 * IFm[,1] +
    bar_beta2_Po1M0 * IFy[,1] +
    bar_beta3_Po1M0 * IFy[,2] +
    bar_sig_Po1M0 * IFsig
  
  
  TE.var = mean((IF_Po1-IF_Po0)^2)/N
  NIDE.var = mean((IF_Po1-IF_Po1M0)^2)/N
  NDE.var = mean((IF_Po1M0-IF_Po0)^2)/N
  coefs = c(Y1-Y0, Y1M0-Y0,Y1-Y1M0)
  vars = c(TE.var,NDE.var,NIDE.var)
  
  a <- list()
  a$coef <- coefs
  a$std.err      <- sqrt(vars )
  a$Wald         <- c(coefs^2/vars)
  names(a$Wald) <- c('TE','NDE','NIDE')
  return(a)
}


#' @export
print.TTS <- function(object){
  a <- object
  
  cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("Coefficients:\n")
  
  df <- data.frame(Estimate = a$coef, Std.Error = a$std.err ,
                   Wald.value = a$Wald,
                   Wald.pval = pchisq(a$Wald,df=1,lower.tail = F) )
  
  printCoefmat(df, digits = 6, signif.stars = T,na.print = "NA",
               tst.ind = 3,P.values=T,has.Pvalue=T)
  
}



