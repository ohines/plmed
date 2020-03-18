#' Multivariate Newton-Raphson fitting method
#'
#' \code{newton_raph} provides an environment for performing iterated
#' Newton Raphson root finding of multivariate functions
#'
#' @param vec_J a function which returns newton-raphson steps via \code{$step}
#' @param Max.it Maximum number of iterations
#' @param prec Parameter precision tolerance. Iteration will stop when the step is
#' smaller than the precision for all parameters.
#' @param par0 Initial values for the parameters to be optimized over.
#' @param step.scale Scale the step size using a size on the interval (0,1]. This is sometimes
#' referred to as the relaxed or damped Newton's method.
#' @param ... Further arguments to be passed to \code{vec_J}
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
#'     return(list(step=step))
#' }
#'
#' newton_raph(vec_J,Max.it=1000,prec=1e-6,par0=c(0,0),step.scale=1)
#'
#' @export
newton_raph <- function(vec_J,Max.it,prec,par0,step.scale,...){
  # For usage about step.scale see:
  # https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
  par = par0
  converged = FALSE
  i = 1
  while(i<Max.it & !converged){
    step = vec_J(par,...)
    par = par - step.scale* step$step
    if(max(abs(step$step))<=prec){converged = TRUE}
    i = i+1
  }
  if (!converged){
    warning(gettextf("Newton Raphson step did not converge to desired tolerance level"),domain = NA)
    browser()}
  return(list(par=par,val = vec_J(par,...) ))
}
