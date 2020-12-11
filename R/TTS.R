#' Natural Mediated Effects
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

TTS <- function(exposure.formula,mediator.formula,outcome.formula,mediator.family="gaussian",data,weights){
  cl <- match.call(expand.dots = FALSE)
  mf <- cl
  m <- match(c("exposure.formula","mediator.formula","outcome.formula","mediator.family","data","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$method <- "TTS"
  mf[[1L]] <- quote(plmed)
  a <- eval(mf, parent.frame())
  a$call <- cl
  return(a)
}



#' @export
fit.TTS <- function(Y,M,X,Z,Mfam,weights=rep(1,N)){
  m_linkinv <- Mfam$linkinv
  m_dev.res <- Mfam$dev.resids
  m_mu.eta  <- Mfam$mu.eta
  N = NROW(X)
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  wt = as.vector(weights)
  N.wt = sum(wt)

  binM <- Mfam$family == "binomial"
  contM <- Mfam$family == "gaussian"
  if(!(binM|contM)) stop('family not recognized')

  Zx=Z
  Zm = cbind(X,Z)
  Zy = cbind(M,X,Z)
  
  X.mod <- glm.fit(Z,X,family=binomial(),weights=wt)
  M.mod <- glm.fit(Zm,M,family=Mfam,weights=wt)
  Y.mod <- glm.fit(Zy,Y,weights=wt)
  
  X.lm <- X.mod$coefficients
  M.lm <- M.mod$coefficients
  Y.lm <- Y.mod$coefficients
  beta  <- c(M.lm[1],Y.lm[1:2]) #target parameters
  
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
  
  D_m0 <- m_mu.eta(fZ)          #derivative of m0  wrt. f(z)
  D_m1 <- m_mu.eta(fZ+beta[1])  #derivative of m1  wrt. f(z) (and also beta1)
  
  #density ratio, K, its derivatives and IF parts are different for different distributions
  if(binM){ 
    K    <- m0/m1       #f(M_i|X=0,Z_i)/f(M_i|X=1,Z_i)
    D_K  <- K*( m1-m0 ) #derivative of K   wrt. f(z)
    D_Kb1<- K*(m1-1)    #derivative of K   wrt. beta1
    density_comp = 0    #part of IF_Po1M0 which depends on other paramerters of f(m|x,z)
    
  } else if(contM) {
    M.var <- M.mod$deviance/N.wt       #Estimate of normal variance
    K    <-   exp((beta[1]/2 -M.res0)*beta[1]/M.var)  #f(M_i|X=0,Z_i)/f(M_i|X=1,Z_i)
    D_K  <- K*beta[1]/M.var                           #derivative of K   wrt. f(z)
    D_Kb1<- K*(beta[1]-M.res0)/M.var                  #derivative of K   wrt. beta1
    #sigma_M (Through K) #only for continuous mediators
    bar_sig_Po1M0 = -sum(K*X_PS*(Y-y1)*(beta[1]/2 - M.res0)*wt)*beta[1]/M.var^2/N.wt
    IFsig = ((M.mod$residuals)^2 - M.var)*wt
    density_comp = bar_sig_Po1M0 * IFsig     #part of IF_Po1M0 which depends on other paramerters of f(m|x,z)
  }

 
  Po0 = cX_PS*(Y- eta_00 ) + eta_00
  Po1 = X_PS*(Y- eta_11) + eta_11
  Po1M0 = X_PS*K*(Y-y1) + cX_PS*(y1-eta_10) + eta_10
  
  Y0 = sum(Po0*wt)/N.wt
  Y1 = sum(Po1*wt)/N.wt
  Y1M0 = sum(Po1M0*wt)/N.wt
  
  ### Expected score function derivatives ###
  #gamma_x (Through the Propensity score)
  bar_gamx_Po0 = (cX_PS*(Y- eta_00 )*PS*wt)%*%Z/N.wt
  bar_gamx_Po1 = -(X_PS*(Y- eta_11)*(1-PS)*wt)%*%Z/N.wt
  bar_gamx_Po1M0 = (-X_PS*K*(Y-y1)*(1-PS)*wt + cX_PS*(y1-eta_10)*PS*wt)%*%Z/N.wt
  
  #gamma_m (Through eta10,00,11 and K)    
  bar_gamm_Po0 = ((1-cX_PS)*beta[2]*D_m0*wt )%*%Z/N.wt
  bar_gamm_Po1 = ((1-X_PS)*beta[2]*D_m1*wt)%*%Z/N.wt
  bar_gamm_Po1M0 = (X_PS*(Y-y1)*D_K*wt  + (1-cX_PS)*beta[2]*D_m0*wt )%*%Z/N.wt
  
  #gamma_y (Through eta10,00,11 and y1) 
  bar_gamy_Po0 = ( (1- cX_PS)*wt)%*%Z/N.wt
  bar_gamy_Po1 = ( (1- X_PS)*wt)%*%Z/N.wt
  bar_gamy_Po1M0 = ( (1- X_PS*K)*wt)%*%Z/N.wt
  
  #beta1 (Through K, eta11)  #These change for logit vs linear
  bar_beta1_Po1 = sum((1 - X_PS)*beta[2]*D_m1*wt )/N.wt
  bar_beta1_Po1M0 = sum( X_PS*(Y-y1)*D_Kb1*wt )/N.wt
  
  #beta2 (through eta10,00,11 and y1)
  bar_beta2_Po0 = sum( (1 - cX_PS)*m0*wt )/N.wt
  bar_beta2_Po1 = sum( (1 - X_PS)*m1*wt )/N.wt
  bar_beta2_Po1M0 = sum( (-X_PS*K*M + cX_PS*M.res0 + m0)*wt )/N.wt
  
  #beta3 (through y1, eta11,10)
  bar_beta3_Po1 = sum((1 - X_PS)*wt)/N.wt
  bar_beta3_Po1M0 = sum((1 - X_PS*K)*wt)/N.wt
  
  #Parameter Influence functions
  IFx = IF.glm(X.mod,Zx)
  IFm = IF.glm(M.mod,Zm)
  IFy = IF.glm(Y.mod,Zy)
  
  IF_Po0 = (Po0 - Y0)*wt*N/N.wt +
    IFx%*%t(bar_gamx_Po0) +
    IFm[,-1]%*%t(bar_gamm_Po0) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po0) +
    bar_beta2_Po0 * IFy[,1]
  
  IF_Po1 = (Po1 - Y1)*wt*N/N.wt  +
    IFx%*%t(bar_gamx_Po1) +
    IFm[,-1]%*%t(bar_gamm_Po1) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po1) +
    bar_beta1_Po1 * IFm[,1] +
    bar_beta2_Po1 * IFy[,1] +
    bar_beta3_Po1 * IFy[,2]
  
  IF_Po1M0 = (Po1M0 - Y1M0)*wt*N/N.wt +
    IFx%*%t(bar_gamx_Po1M0) +
    IFm[,-1]%*%t(bar_gamm_Po1M0) +
    IFy[,-(1:2)]%*%t(bar_gamy_Po1M0) +
    bar_beta1_Po1M0 * IFm[,1] +
    bar_beta2_Po1M0 * IFy[,1] +
    bar_beta3_Po1M0 * IFy[,2] +
    density_comp

  TE.var = sum((IF_Po1-IF_Po0)^2)/(N^2)
  NIDE.var = sum((IF_Po1-IF_Po1M0)^2)/(N^2)
  NDE.var = sum((IF_Po1M0-IF_Po0)^2)/(N^2)
  
  coefs = c(Y1-Y0, Y1M0-Y0,Y1-Y1M0)
  vars = c(TE.var,NDE.var,NIDE.var)

  a <- list()
  a$coef <- coefs
  a$std.err      <- sqrt(vars )
  a$Wald         <- c(coefs^2/vars)
  names(a$Wald) <- c('TE','NDE','NIDE')
  return(a)
}





