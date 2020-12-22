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
