#' Multivariate Newton-Raphson fitting method
#'
#' \code{newton_raph} provides an environment for performing iterated
#' Newton Raphson root finding of multivariate functions with an optional Backtracking line search,
#' For more on line search see: http://cosmos.phy.tufts.edu/~danilo/AST16/Material/RootFinding.pdf page479
#'
#' @param vec_J a function which returns newton-raphson steps via \code{$step} and the functions via \code{$vec}
#' @param par0 Initial values for the parameters to be optimized over.
#' @param ... Further arguments to be passed to \code{vec_J}
#' @param LSearch Option to perform the Backtracking line search algorithm. Defaults to \code{FALSE}
#' @param w Weights for Backtracking algorithm. Must be non-negative with length equal to \code{length(par0)}
#' @param Max.it Maximum number of iterations (Function calls)
#' @param prec Parameter precision tolerance. Iteration will stop when the step is
#' smaller than the precision for all parameters.
#' @param step.scale Optional scaling step to perform damped Newton Raphson.
#'
#' @return A list containing the \code{par} estimates and
#' an evaluation of the \code{vec_J} function at these parameters.
#' @examples
#' #Optimisation of the function f(x,y) = (x*exp(y) - 1, -1-x^2 + y)
#' vec_J <- function(par){
#'     x=par[1]
#'     y=par[2]
#'     f = c(x*exp(y) - 1, -1-x^2 + y)
#'     fdash = matrix(c(exp(y),-2*x,x*exp(y),1),ncol=2)
#'     step = solve(fdash,f)
#'     return(list(step=step, vec=f))
#' }
#'
#' newton_raph(vec_J,par0=c(0,0))
#'
#' @export
newton_raph <- function(vec_J,par0,...,LSearch = FALSE, w=rep.int(1,length(par0)),Max.it=2000,prec=1e-5, step.scale = 1){
  # For more on Linesearch:
  # https://en.wikipedia.org/wiki/Backtracking_line_search
  # http://cosmos.phy.tufts.edu/~danilo/AST16/Material/RootFinding.pdf page479
  par = par0
  converged = FALSE
  i = 1 #records function calls
  if(!LSearch){ #Do standard newton raphson
    
    while(i<Max.it & !converged){
      step = vec_J(par,...)
      par = par - step.scale*step$step
      if(max(abs(step$step))<=prec){converged = TRUE}
      i = i+1
    }
  }else{#Do Newton raphson with additional Backtracking line search
    step = vec_J(par,...)
    f0  = sum(step$vec*step$vec*w) #scalar which we want to minimise
    p = -step$step
    while(i<Max.it & !converged){
      lam = 1 #step scale starts at normal Newton Raphson lenght
      Backtrack = TRUE
      while(i<Max.it & Backtrack){
        parK = par + lam*p
        stepK = vec_J(parK,...)
        fK = sum(stepK$vec*stepK$vec*w)
        
        accept = fK <= (1-2e-4*lam)*f0 #accept step with condition
        #cat('Trial =',fK,'previous = ',f0, 'accept = ', accept, 'lam was',lam,'\n')
        if(accept|lam==0.01){
          par  = parK
          f0 = fK
          p  = -stepK$step
          Backtrack  = FALSE
          accept=FALSE
        }
        lam = max(min(as.vector((f0/(fK+f0))),lam/2),0.01) #If not accepted, Try again with a smaller step (capped at 0.01 of a Newton Raphson step)
        i = i+1
      }
      if(max(abs(p))<=prec){converged = TRUE}
    }
  }
  if (!converged){
    if(LSearch){
      method = 'Line Search'
    }else{method = 'Standard Newton Raphson'}
    warning(gettextf("Newton Raphson step did not converge to desired tolerance level using: %s",method),domain = NA) }
  return(list(par=par,val = vec_J(par,...) ))
}
