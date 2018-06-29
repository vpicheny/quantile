#' Main function to find quantiles.
##' @title Main solver
##' @details Either the distribution of the input \code{Xdistrib} or a sample x must be given. 
##' If x is given, the problem is treated as discrete.
##' 
##' In the standard setting, a large number of points (\code{n.large}) is generated using \code{Xdistrib}, out of which 
##' \code{n} useful integration points are selected. The SUR criterion is then evaluated at the most promising \code{n.candidates} 
##' points. Finally, a local optimization is performed using \code{BFGS} from the best candidate.
##' 
##' Maximum recommended values are 5,000 for \code{n}, 1e6 for \code{n.large} and 1,000 for \code{n.candidates}.  
##'
##' @param model an object of class km
##' @param fun the function of interest
##' @param alpha the quantile level
##' @param n.ite number of iterations (points to add)
##' @param n the number of points used for integration
##' @param n.large (optional) a larger number of points used prior to integration
##' @param n.candidates (optional) a smaller number of candidate points on which the criterion is computed
##' @param Xdistrib,x \code{Xdistrib} is a function that returns a sample of x given an integer, and x is a matrix (see details)
##' @param n.cluster number of cores used (requires the libraries \code{forreach} and \code{doparallel})
##' @param seed the seed
##' @param cov.reestim Boolean; if TRUE, the GP parameters are re-estimated at each iteration
##' 
##' @return
##' A list with all.qn (all the quantiles estimated) and model (the last km model)
##'
##' @export
##' @importFrom stats pnorm qnorm rnorm dnorm runif
##' @importFrom MASS mvrnorm
##' @importFrom KrigInv precomputeUpdateData
##' @useDynLib Qlab, .registration = TRUE
##' @author Victor Picheny
##' @references
##' T. Labopin-Richard, V. Picheny, "Sequential design of experiments for estimating quantiles of black-box functions", 
##' Statistica Sinica, 2017, \emph{doi:10.5705/ss.202016.0160}
##' @examples
##' \dontrun{
##' library(DiceDesign)
##' #--------------------------------------------------------#
##' # Set problem parameters
##' fun <- branin
##' d <- 2
##' n.init <- 6
##' n.ite <- 24
##' seed <- 42
##' n <- 2e3
##' #--------------------------------------------------------#
##' # Define distribution over the input X
##' mu <- rep(.5, d)
##' Sigma <- matrix(rep(.05, d*d), d,d)
##' diag(Sigma) <- .1
##' Xdistrib <- function(n) return(mvrnorm(n=n, mu=mu, Sigma))
##' #--------------------------------------------------------#
##' # Initial set of observations (rescaled to fit Xdistrib)
##' x.init <- lhsDesign(n.init, d, seed=seed)$design
##' x.init <- mu + qnorm(x.init) %*% chol(Sigma)
##' y.init <- as.numeric(apply(x.init, 1, fun))
##' #--------------------------------------------------------#
##' # Initial kriging model
##' model <- km(~., design=data.frame(x=x.init), response=y.init, 
##'             lower=rep(.05,d), upper=rep(1,d), control=list(trace=FALSE))
##' #--------------------------------------------------------#
##' # Sequential design
##' res <- findQ(model=model, fun=fun, alpha=.95, n.ite=n.ite, n=n, n.cluster=1, seed=seed, n.large=1e5,
##'              n.candidates=100, Xdistrib=Xdistrib, cov.reestim = TRUE)
##' #--------------------------------------------------------#
##' # Plot DoE and quantile estimates
##' plot(res$model@X[,1], res$model@X[,2])
##' plot(res$all.qn)
##' }
##' 
findQ <- function(model, fun, alpha=.95, n.ite=45, n=1e3, n.large=NULL, n.cluster=1, seed=42,
                  n.candidates=NULL, Xdistrib=NULL, x=NULL, cov.reestim = TRUE){
  
  set.seed(seed)
  if (is.null(n.candidates)) n.candidates <- n
  
  if (is.null(Xdistrib) && is.null(x)) {
    cat("Either the distribution of x (Xdistrib) or a sample (x) should be provided")
    return(NA)
  }
  
  if (n.cluster>1) {
    library('foreach')
    library(doParallel)
    cl2 <- makeCluster(c(rep("localhost", n.cluster)), type = "SOCK")
    cluster = cl2
    registerDoParallel(n.cluster)
    getDoParWorkers()
  }

  #--------------------------------------------------------#
  if (!is.null(n.large))     x.large <- Xdistrib(n=n.large)
  
  all.qn <- rep(NA, n.ite+1)
  
  #--------------------------------------------------------#
  cat("--- Starting ", n.ite, " interations -------------------------------\n")
  cat("ite | qn | xnew | ynew \n")
  # Main loop
  for (ite in 1:n.ite) {
    # cat("--- ", ite, " -------------------------------\n")
    tstart <- Sys.time()
    
    # x sample
    if (!is.null(Xdistrib))   x <- Xdistrib(n=n)
    
    # Precalculations
    precalc.data.x <- precomputeUpdateData(model, integration.points=x)
    pred.x <- predict(model, data.frame(x=x), "UK", checkNames=FALSE, light.return = TRUE)
    data.x <- list(x=x, mean=pred.x$mean, sd=pred.x$sd, precalc.data=precalc.data.x)
    
    # Compute qn  
    if (!is.null(n.large)) {
      p.large  <- predict(model, data.frame(x=x.large), "UK", checkNames=FALSE, light.return = TRUE)
      all.qn[ite] <- qn <- p.large$mean[order(p.large$mean)[round(alpha*n.large)]]
    } else {
      all.qn[ite] <- qn <- pred.x$mean[order(pred.x$mean)[round(alpha*n)]]
    }
    
    # Reduce x to x.candidates if needed
    if (n == n.candidates) {
      x.candidates <- x
    } else {
      min.prob <- 0.001
      dens <- dnorm((qn - pred.x$mean)/pred.x$sd)
      dens.sum <- sum(dens)
      prob.n <- pmax(dens/dens.sum, min.prob/n)
      prob.n <- c(0, prob.n/sum(prob.n))
      prob.n.cum <- cumsum(prob.n)
      my.indices <- findInterval(runif(n.candidates), prob.n.cum, all.inside = TRUE)
      x.candidates <- x[my.indices,]
    }
    
    # Compute criterion
    if (n.cluster<2) {
      Vn <- apply(x.candidates, 1, crit_varQ, model=model, data.x=data.x, x=x, alpha=alpha, beta=4)
    } else {
      Vn <- foreach(i=1:n.candidates, .combine='rbind') %dopar% crit_varQ(x.candidates[i,], model=model, data.x=data.x, x=x, alpha=alpha, beta=4)
    }
    
    # Local descent based on the best x.candidates point
    Vn[is.na(Vn)] <- min(Vn, na.rm = TRUE)
    newX <- x.candidates[which.max(Vn),,drop=FALSE]
    res <- optim(par=newX, fn=crit_varQ, model=model, data.x=data.x, x=x, alpha=alpha, beta=4, 
                 lower=rep(min(x), d), upper=rep(max(x),d),control=list(fnscale=-1, maxit=5), method="L-BFGS-B")
    
    if (res$value > max(Vn)) newX <- res$par
    
    # New observation
    newy <- fun(newX)
    
    #--------------------------------------------------------#
    # Model update
    model <- update(object=model, newX=newX, newy=newy, cov.reestim=cov.reestim)

    tstop <- Sys.time()
    # cat(tstop - tstart, "\n")
    cat(c(ite, qn, newX, newy), "\n")
  }
  if (!is.null(n.large)) {
    p.large  <- predict(model, data.frame(x=x.large), "UK", checkNames=FALSE, light.return = TRUE)
    all.qn[ite+1] <- p.large$mean[order(p.large$mean)[round(alpha*n.large)]]
  } else {
    all.qn[ite+1] <- pred.x$mean[order(pred.x$mean)[round(alpha*n)]]
  }
  if (n.cluster>1)  stopCluster(cl2)
  return(list(all.qn=all.qn, model=model))
}