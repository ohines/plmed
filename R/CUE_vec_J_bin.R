#' Newton Raphson Steps for plmed fitting
#'
#' \code{CUE_vec_J_bin} Generates newton raphson steps for the \code{\link[plmed]{plmed}}
#' fitting functions.
#'
#' @param par vector containing beta and augmented nuisance parameters, optionally with a
#' lagrange multiplier a the final componenet.
#' @param X Vector of exposures
#' @param M Vector of mediators
#' @param Y Vector of outcomes
#' @param Z.matrix Design matrix of confounder models
#' @param method The method to be used. See ‘Details’.
#' @param med_prop Hypothesised mediation proportion for constrained optimisation.
#' Requires method \code{'CUE'} or \code{'Two-step'}.
#' @param Sig Constant covariance matrix must be prespecified when using \code{'Two-step'} method.
#' @return A list containing the newton-raphson step for use with \code{\link[plmed]{newton_raph}}
#' as well as variance and score estimates.
#' @export
CUE_vec_J_bin <- function(par,X,M,Y,Z.matrix,method='G',med_prop=NULL,Sig=NULL){
  theta=par
  N = length(M)
  nt=1 #Assume X is binomial with nt trials
  
  cue   <- method=='CUE'
  G     <- method=='G'
  Tstep <- method=='Two-step'

  if (!any(cue,G,Tstep)){
    warning(gettextf("method = '%s' is not supported. Using 'G-estimation'", method),domain = NA)
    G = TRUE
  }
  
  tZ.matrix = t(Z.matrix) #compute once for speed
  #### Read in parameters and construct residuals ####
  d = length(theta)
  if(!(d%%6)%in%c(3,4)) { 
    warning(gettextf("Theta has length = '%d \n
                     Theta should cotain beta, 6 nuisance components and optionally a lagrange multiplier", d),domain = NA)
  }

  p = floor({d-3}/6) #length of nuisance parameter vector note #beta = theta[1:3]
  
  p1 = as.vector(plogis(Z.matrix%*%theta[4:{p+3}]))
  p2 = as.vector(plogis(Z.matrix%*%theta[{p+4}:{2*p+3}]))

  x.res1   = X-nt*p1
  x.res2   = X-nt*p2
  m.res1   = M - theta[1]*X - as.vector(Z.matrix%*%theta[{2*p+4}:{3*p+3}])
  m.res2   = M - theta[1]*X - as.vector(Z.matrix%*%theta[{3*p+4}:{4*p+3}])
  y.res1   = Y - theta[2]*M - theta[3]*X - as.vector(Z.matrix%*%theta[{4*p+4}:{5*p+3}])
  y.res2   = Y - theta[2]*M - theta[3]*X - as.vector(Z.matrix%*%theta[{5*p+4}:{6*p+3}])

  hdot1     = nt*p1*(1-p1)
  hdot2     = nt*p2*(1-p2)

  #### Estimating Equations and first derivatives ####

  U = c(sum(x.res1*m.res1),sum(m.res2*y.res1),sum(x.res2*y.res2))/N
  dUdbeta = -matrix(c(sum(X*x.res1),0,0, #dU1dbeta
                      sum(X*y.res1),sum(M*m.res2),sum(X*m.res2), #dU2dbeta
                      0,sum(M*x.res2),sum(X*x.res2)), #dU3dbeta
                    nrow=3,byrow = TRUE)/N

  if(cue|G){
    V1 = tZ.matrix%*%{m.res1*hdot1}/N
    V2 = tZ.matrix%*%x.res1/N
    V3 = tZ.matrix%*%y.res1/N
    V4 = tZ.matrix%*%m.res2/N
    V5 = tZ.matrix%*%{y.res2*hdot2}/N
    V6 = tZ.matrix%*%x.res2/N
    
    V = c(V1,V2,V3,V4,V5,V6)

    p0 = rep.int(0,p)

    dUdtheta = cbind(dUdbeta,matrix(c(-V1,p0,-V2,p0,p0,p0, #dU1dgamma
                                       p0,p0,p0,-V3,-V4,p0, #dU2dgamma
                                       p0,-V5,p0,p0,p0,-V6),#dU3dgamma
                                       nrow=3,byrow = TRUE))


    #dVdbeta = -rbind(tZ.matrix%*%cbind(X*hdot1,0,0), #dV1dbeta
    #                 matrix(0,p,3),                  #dV2dbeta
    #                 tZ.matrix%*%cbind(0,M,X),       #dV3dbeta
    #                 tZ.matrix%*%cbind(X,0,0),       #dV4dbeta
    #                 tZ.matrix%*%cbind(0,M*hdot2,X*hdot2), #dV5dbeta
    #                 matrix(0,p,3))/N                #dV6dbeta
    pp0 = matrix(0,nrow=p,ncol=p)
    
    Z.sq =    tZ.matrix%*%Z.matrix
    Z.hdot1 = tZ.matrix%*%{Z.matrix*hdot1} 
    Z.hdot2 = tZ.matrix%*%{Z.matrix*hdot2}

    Z.hdotdot1 = tZ.matrix%*%{Z.matrix*hdot1*{1-2*p1}*m.res1}
    Z.hdotdot2 = tZ.matrix%*%{Z.matrix*hdot2*{1-2*p2}*y.res2}

    #dV1dgamma = -cbind(Z.hdotdot1,pp0,Z.hdot1,pp0,pp0,pp0)/N
    #dV2dgamma = -cbind(Z.hdot1,pp0,pp0,pp0,pp0,pp0)/N
    #dV3dgamma = -cbind(pp0,pp0,pp0,pp0,Z.sq,pp0)/N
    #dV4dgamma = -cbind(pp0,pp0,pp0,Z.sq,pp0,pp0)/N
    #dV5dgamma = -cbind(pp0,Z.hdotdot2,pp0,pp0,pp0,Z.hdot2)/N
    #dV6dgamma = -cbind(pp0,Z.hdot2,pp0,pp0,pp0,pp0)/N
    #dVdtheta = cbind(dVdbeta,rbind(dV1dgamma,dV2dgamma,dV3dgamma,dV4dgamma,dV5dgamma,dV6dgamma))
    
    dV1dtheta = -cbind(tZ.matrix%*%{X*hdot1},0,0,Z.hdotdot1,pp0,Z.hdot1,pp0,pp0,pp0)/N
    dV2dtheta = -cbind(0,0,0                    ,Z.hdot1,pp0,pp0,pp0,pp0,pp0)/N
    dV3dtheta = -cbind(0,tZ.matrix%*%cbind(M,X) ,pp0,pp0,pp0,pp0,Z.sq,pp0)/N
    dV4dtheta = -cbind(tZ.matrix%*%X,0,0         ,pp0,pp0,pp0,Z.sq,pp0,pp0)/N
    dV5dtheta = -cbind(0,tZ.matrix%*%cbind(M*hdot2,X*hdot2),pp0,Z.hdotdot2,pp0,pp0,pp0,Z.hdot2)/N
    dV6dtheta = -cbind(0,0,0                    ,pp0,Z.hdot2,pp0,pp0,pp0,pp0)/N
    
    dVdtheta = rbind(dV1dtheta,dV2dtheta,dV3dtheta,dV4dtheta,dV5dtheta,dV6dtheta)
    
    sig11 = sum(x.res1^2*m.res1^2)/N
    sig22 = sum(m.res2^2*y.res1^2)/N
    sig33 = sum(x.res2^2*y.res2^2)/N

    sig12 = sum(x.res1*m.res1*m.res2*y.res1)/N
    sig13 = sum(x.res1*m.res1*x.res2*y.res2)/N
    sig23 = sum(m.res2*y.res1*x.res2*y.res2)/N

    Sig = matrix(c(sig11,sig12,sig13,
                   sig12,sig22,sig23,
                   sig13,sig23,sig33),ncol=3)

  } else if(is.null(Sig)){
    stop(gettextf("A constant covriance matrix must be pre-specified when using the Two-step"), domain = NA)
  }

  A = solve(Sig)
  AU = A%*%U
  f_out = N*crossprod(AU,U)

  if(cue){
    dsig11dbeta = -2*c(sum(X*x.res1^2*m.res1),0,0)/N
    dsig22dbeta = -2*c(sum(X*y.res1^2*m.res2),sum(M*m.res2^2*y.res1),sum(X*m.res2^2*y.res1))/N
    dsig33dbeta = -2*c(0,sum(M*x.res2^2*y.res2),sum(X*x.res2^2*y.res2))/N

    dsig12dbeta = c(sum(-X*x.res1*y.res1*{m.res1+m.res2}),
                    sum(-M*x.res1*m.res1*m.res2),
                    sum(-X*x.res1*m.res1*m.res2))/N
    dsig13dbeta = c(sum(-X*x.res1*x.res2*y.res2),
                    sum(-M*x.res1*x.res2*m.res1),
                    sum(-X*x.res1*x.res2*m.res1))/N
    dsig23dbeta = c(sum(-X*x.res2*y.res1*y.res2),
                    sum(-M*x.res2*m.res2*{y.res1+y.res2}),
                    sum(-X*x.res2*m.res2*{y.res1+y.res2}))/N

    dsig11dgamx1 = tZ.matrix%*%{m.res1^2*hdot1*x.res1}/N
    dsig11dgamm1 = tZ.matrix%*%{x.res1^2*m.res1}/N
    dsig22dgamm2 = tZ.matrix%*%{y.res1^2*m.res2}/N
    dsig22dgamy1 = tZ.matrix%*%{m.res2^2*y.res1}/N
    dsig33dgamx2 = tZ.matrix%*%{y.res2^2*x.res2*hdot2}/N
    dsig33dgamy2 = tZ.matrix%*%{x.res2^2*y.res2}/N

    dsig11dgamma = -2*c(dsig11dgamx1,p0,dsig11dgamm1,p0,p0,p0)
    dsig22dgamma=  -2*c(p0,p0,p0,dsig22dgamm2,dsig22dgamy1,p0)
    dsig33dgamma = -2*c(p0,dsig33dgamx2,p0,p0,p0,dsig33dgamy2)


    dsig12dgamma = -c(tZ.matrix%*%{hdot1*m.res1*m.res2*y.res1},
                      p0,
                      tZ.matrix%*%{x.res1*m.res2*y.res1},
                      tZ.matrix%*%{x.res1*m.res1*y.res1},
                      tZ.matrix%*%{x.res1*m.res1*m.res2},
                      p0)/N
    dsig13dgamma = -c(tZ.matrix%*%{hdot1*m.res1*y.res2*x.res2},
                      tZ.matrix%*%{hdot2*m.res1*y.res2*x.res1},
                      tZ.matrix%*%{x.res1*y.res2*x.res2},
                      p0,
                      p0,
                      tZ.matrix%*%{m.res1*x.res2*x.res1})/N
    dsig23dgamma = -c(p0,
                      tZ.matrix%*%{hdot2*y.res1*m.res2*y.res2},
                      p0,
                      tZ.matrix%*%{x.res2*y.res1*y.res2},
                      tZ.matrix%*%{m.res2*x.res2*y.res2},
                      tZ.matrix%*%{m.res2*x.res2*y.res1} )/N

    dsig11dtheta = c(dsig11dbeta,dsig11dgamma)
    dsig22dtheta = c(dsig22dbeta,dsig22dgamma)
    dsig33dtheta = c(dsig33dbeta,dsig33dgamma)
    dsig12dtheta = c(dsig12dbeta,dsig12dgamma)
    dsig13dtheta = c(dsig13dbeta,dsig13dgamma)
    dsig23dtheta = c(dsig23dbeta,dsig23dgamma)

    B = dUdtheta -
      AU[1]*rbind(dsig11dtheta,dsig12dtheta,dsig13dtheta) -
      AU[2]*rbind(dsig12dtheta,dsig22dtheta,dsig23dtheta) -
      AU[3]*rbind(dsig13dtheta,dsig23dtheta,dsig33dtheta)
  } else {B=dUdbeta}

  dfdbeta = N*crossprod(AU,dUdbeta + B[,1:3])
  MX = sum(M*X)
  XX = sum(X*X)
  dU2dbetadbeta = matrix(c(0,MX,XX,MX,0,0,XX,0,0),nrow=3)/N


  #### Lots of second derivatives ####
  #calculation of d2fdbetadgamma
  if (cue){
    XZ = tZ.matrix%*%X
    MZ = tZ.matrix%*%M

    p000 = rep.int(0,3+6*p)

    dU1dbetadtheta = matrix(c(0,0,0,tZ.matrix%*%{hdot1*X},p0,p0,p0,p0,p0,
                              p000,
                              p000),byrow=T,nrow=3)/N

    dU2dbetadtheta = cbind(dU2dbetadbeta,matrix(c(p0,p0,p0,p0,XZ ,p0,
                                                  p0,p0,p0,MZ ,p0,p0,
                                                  p0,p0,p0,XZ ,p0,p0),byrow=T,nrow=3)/N)

    dU3dbetadtheta = matrix(c( p000,
                               0,0,0,p0,tZ.matrix%*%{hdot2*M} ,p0,p0,p0,p0,
                               0,0,0,p0,tZ.matrix%*%{hdot2*X} ,p0,p0,p0,p0),byrow=T,nrow=3)/N

    
    dsig11db1db1 = 2*sum(X^2* x.res1^2)
    dsig22db1db1 = 2*sum(X^2* y.res1^2)
    dsig22db1db2 = 4*sum(X*M* y.res1*m.res2)
    dsig22db1db3 = 4*sum(X^2* y.res1*m.res2)
    dsig22db2db2 = 2*sum(M^2* m.res2^2)
    dsig22db2db3 = 2*sum(X*M* m.res2^2)
    dsig22db3db3 = 2*sum(X^2* m.res2^2)
    dsig33db2db2 = 2*sum(M^2* x.res2^2)
    dsig33db2db3 = 2*sum(X*M* x.res2^2)
    dsig33db3db3 = 2*sum(X^2* x.res2^2)

    dsig11db1dgamx1 = 4*tZ.matrix%*%{hdot1*X*x.res1*m.res1}
    dsig11db1dgamm1 = 2*tZ.matrix%*%{X*x.res1^2}
    dsig22db1dgamm2 = 2*tZ.matrix%*%{X*y.res1^2}
    dsig22db2dgamm2 = 4*tZ.matrix%*%{M*y.res1*m.res2}
    dsig22db3dgamm2 = 4*tZ.matrix%*%{X*y.res1*m.res2}
    dsig22db1dgamy1 = 4*tZ.matrix%*%{X*y.res1*m.res2}
    dsig22db2dgamy1 = 2*tZ.matrix%*%{M*m.res2^2}
    dsig22db3dgamy1 = 2*tZ.matrix%*%{X*m.res2^2}
    dsig33db2dgamx2 = 4*tZ.matrix%*%{M*y.res2*x.res2*hdot2}
    dsig33db3dgamx2 = 4*tZ.matrix%*%{X*y.res2*x.res2*hdot2}
    dsig33db2dgamy2 = 2*tZ.matrix%*%{M*x.res2^2}
    dsig33db3dgamy2 = 2*tZ.matrix%*%{X*x.res2^2}

    dsig11dbetadtheta = matrix(c( dsig11db1db1,0,0,dsig11db1dgamx1, p0,dsig11db1dgamm1,p0,p0,p0,
                                  p000,p000),byrow=T,nrow=3)/N

    dsig22dbetadtheta = matrix(c(dsig22db1db1,dsig22db1db2,dsig22db1db3,p0,p0,p0,dsig22db1dgamm2,dsig22db1dgamy1,p0,
                                 dsig22db1db2,dsig22db2db2,dsig22db2db3,p0,p0,p0,dsig22db2dgamm2,dsig22db2dgamy1,p0,
                                 dsig22db1db3,dsig22db2db3,dsig22db3db3,p0,p0,p0,dsig22db3dgamm2,dsig22db3dgamy1,p0),byrow=T,nrow=3)/N

    dsig33dbetadtheta = matrix(c(p000,
                                 0,dsig33db2db2,dsig33db2db3,p0,dsig33db2dgamx2,p0,p0,p0,dsig33db2dgamy2,
                                 0,dsig33db2db3,dsig33db3db3,p0,dsig33db3dgamx2,p0,p0,p0,dsig33db3dgamy2),byrow=T,nrow=3)/N

    dsig12db1db1 = 2*sum(X^2*   x.res1*y.res1)
    dsig12db1db2 =   sum(X*M*   x.res1*{m.res1+m.res2})
    dsig12db1db3 =   sum(X^2*   x.res1*{m.res1+m.res2})
    dsig13db1db2 =   sum(X*M*   x.res1*x.res2)
    dsig13db1db3 =   sum(X^2*   x.res1*x.res2)
    dsig23db1db2 =   sum(X*M*   x.res2*{y.res1+y.res2})
    dsig23db1db3 =   sum(X^2*   x.res2*{y.res1+y.res2})
    dsig23db2db2 = 2*sum(M^2*   x.res2*m.res2)
    dsig23db2db3 = 2*sum(X*M*   x.res2*m.res2)
    dsig23db3db3 = 2*sum(X^2*   x.res2*m.res2)

    dsig12db1dgamx1 = tZ.matrix%*%{hdot1*X*y.res1*{m.res1+m.res2}}
    dsig12db1dgamy1 = tZ.matrix%*%{X*x.res1*{m.res1+m.res2}}
    dsig12db1dgamm1 = tZ.matrix%*%{X*x.res1*y.res1}
    dsig12db1dgamm2 = dsig12db1dgamm1
    dsig12db2dgamx1 = tZ.matrix%*%{hdot1*M*m.res1*m.res2}
    dsig12db2dgamm1 = tZ.matrix%*%{M*x.res1*m.res2}
    dsig12db2dgamm2 = tZ.matrix%*%{M*x.res1*m.res1}
    dsig12db3dgamx1 = tZ.matrix%*%{hdot1*X*m.res1*m.res2}
    dsig12db3dgamm1 = tZ.matrix%*%{X*x.res1*m.res2}
    dsig12db3dgamm2 = tZ.matrix%*%{X*x.res1*m.res1}

    dsig13db1dgamx1 = tZ.matrix%*%{hdot1*X*y.res2*x.res2}
    dsig13db1dgamx2 = tZ.matrix%*%{X*hdot2*y.res2*x.res1}
    dsig13db1dgamy2 = tZ.matrix%*%{X*x.res1*x.res2}
    dsig13db2dgamx1 = tZ.matrix%*%{hdot1*M*x.res2*m.res1}
    dsig13db2dgamx2 = tZ.matrix%*%{M*hdot2*m.res1*x.res1}
    dsig13db2dgamm1 = tZ.matrix%*%{M*x.res2*x.res1}
    dsig13db3dgamx1 = tZ.matrix%*%{hdot1*X*x.res2*m.res1}
    dsig13db3dgamx2 = tZ.matrix%*%{X*hdot2*m.res1*x.res1}
    dsig13db3dgamm1 = tZ.matrix%*%{X*x.res2*x.res1}

    dsig23db1dgamx2 = tZ.matrix%*%{X*hdot2*y.res1*y.res2}
    dsig23db1dgamy1 = tZ.matrix%*%{X*x.res2*y.res2}
    dsig23db1dgamy2 = tZ.matrix%*%{X*x.res2*y.res1}
    dsig23db2dgamx2 = tZ.matrix%*%{M*hdot2*m.res2*{y.res1+y.res2}}
    dsig23db2dgamm2 = tZ.matrix%*%{M*x.res2*{y.res1+y.res2}}
    dsig23db2dgamy1 = tZ.matrix%*%{M*x.res2*m.res2}
    dsig23db2dgamy2 = dsig23db2dgamy1
    dsig23db3dgamx2 = tZ.matrix%*%{X*hdot2*m.res2*{y.res1+y.res2}}
    dsig23db3dgamm2 = tZ.matrix%*%{X*x.res2*{y.res1+y.res2}}
    dsig23db3dgamy1 = tZ.matrix%*%{X*x.res2*m.res2}
    dsig23db3dgamy2 = dsig23db3dgamy1



    dsig12dbetadtheta = matrix(c(dsig12db1db1,dsig12db1db2,dsig12db1db3,dsig12db1dgamx1,p0,dsig12db1dgamm1,dsig12db1dgamm2,dsig12db1dgamy1,p0,
                                 dsig12db1db2,0,0                      ,dsig12db2dgamx1,p0,dsig12db2dgamm1,dsig12db2dgamm2,p0,p0,
                                 dsig12db1db3,0,0                      ,dsig12db3dgamx1,p0,dsig12db3dgamm1,dsig12db3dgamm2,p0,p0),byrow=T,nrow=3)/N

    dsig13dbetadtheta = matrix(c(0,dsig13db1db2,dsig13db1db3,dsig13db1dgamx1 ,dsig13db1dgamx2,p0              ,p0,p0,dsig13db1dgamy2 ,
                                 dsig13db1db2,0,0           ,dsig13db2dgamx1 ,dsig13db2dgamx2,dsig13db2dgamm1 ,p0,p0,p0,
                                 dsig13db1db3,0,0           ,dsig13db3dgamx1 ,dsig13db3dgamx2,dsig13db3dgamm1 ,p0,p0,p0),byrow=T,nrow=3)/N

    dsig23dbetadtheta = matrix(c(0           ,dsig23db1db2,dsig23db1db3,p0,dsig23db1dgamx2,p0,p0              ,dsig23db1dgamy1,dsig23db1dgamy2,
                                 dsig23db1db2,dsig23db2db2,dsig23db2db3,p0,dsig23db2dgamx2,p0,dsig23db2dgamm2 ,dsig23db2dgamy1,dsig23db2dgamy2,
                                 dsig23db1db3,dsig23db2db3,dsig23db3db3,p0,dsig23db3dgamx2,p0,dsig23db3dgamm2 ,dsig23db3dgamy1,dsig23db3dgamy2),byrow=T,nrow=3)/N



    d2fdbetadgamma = N*(
      2*crossprod(B[,1:3],A%*%B) +
        2* AU[1]*dU1dbetadtheta +
        2* AU[2]*dU2dbetadtheta +
        2* AU[3]*dU3dbetadtheta -
        AU[1]^2 *dsig11dbetadtheta -
        AU[2]^2 *dsig22dbetadtheta -
        AU[3]^2 *dsig33dbetadtheta -
        2*AU[1]*AU[2]*dsig12dbetadtheta -
        2*AU[1]*AU[3]*dsig13dbetadtheta -
        2*AU[2]*AU[3]*dsig23dbetadtheta)
  }

  #### Return output ####
  #Use G-estimation equations or constrained lagrangians as required
  vec = c(U,V)
  J = rbind(dUdtheta,dVdtheta)

  if (cue){
    vec = c(dfdbeta,V)
    J = rbind(d2fdbetadgamma,dVdtheta)
  } else if (method == 'Two-step'){
    vec = dfdbeta  #dfdbeta with Sig constant
    J   = N*(2*crossprod(dUdbeta,A%*%dUdbeta) + 2*AU[2]*dU2dbetadbeta) #d2fdbetadbeta with Sig constant
    d = 4
  }

  if(!is.null(med_prop) & (cue|Tstep) ){ #adds the lagrangian constraint
    g = theta[1]*theta[2]*(1-med_prop)-med_prop*theta[3]
    dgdb = c(theta[2]*(1-med_prop),theta[1]*(1-med_prop),-med_prop)

    vec = c(vec,g) #extend vec and J
    J = cbind(rbind(J,0),0)

    vec[1:3] <- vec[1:3] +theta[d]*dgdb  #theta[d] is the lagrange multiplier (lambda)
    J[1,2] <- J[1,2] + theta[d]*(1-med_prop)
    J[2,1] <- J[2,1] + theta[d]*(1-med_prop)
    J[1:3,d] <- dgdb
    J[d,1:3] <- dgdb
  }

  step = solve(J,vec)

  if (Tstep){ #keeps nuisance parameters constant
    step = c(step[1:3],rep.int(0,6*p),step[4])
  }

  inv_dUdbeta = solve(dUdbeta)
  PhiPhiT = tcrossprod(inv_dUdbeta%*%Sig,inv_dUdbeta)

  T1 = N*(theta[1]^2)/PhiPhiT[1,1]
  T2 = N*(theta[2]^2)/PhiPhiT[2,2]
  T3 = N*(theta[3]^2)/PhiPhiT[3,3]

  rho = PhiPhiT[1,2]/sqrt(PhiPhiT[1,1]*PhiPhiT[2,2])
  RSTest = T1*T2/(T1+T2+2*rho*sqrt(T1*T2))

  return(list(step=step,
              var=diag(PhiPhiT)/N ,
              score=f_out,
              RSTest = RSTest,
              T_stats = c(T1,T2,T3),
              Sig=Sig
              ))
}
