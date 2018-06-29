# ' Compute the future quantile variance given a potential new observation
##' @title Future quantile variance
##' @param xnew the new point
##' @param model object of class km
##' @param data.x data related to x 
##' @param x integration points
##' @param alpha quantile level
##' @param beta a threshold
##' @return the criterion
##' @importFrom KrigInv computeQuickKrigcov predict_nobias_km
##' @importFrom GPareto checkPredict
##' @useDynLib Qlab, .registration = TRUE
##' @export
##' @author Victor Picheny
##' @references
##' T. Labopin-Richard, V. Picheny, "Sequential design of experiments for estimating quantiles of black-box functions", 
##' Statistica Sinica, 2017, \emph{doi:10.5705/ss.202016.0160}
##' @examples
##' \dontrun{
##' compute_crit(xnew, model, data.x, x, alpha, beta=4)
##' }

crit_varQ <- function(xnew, model, data.x, x, alpha, beta=4){
  #-----------------------------------------------------------#
  # Cette fonction calcule soit l'espérance de la proba de dépasser 
  # le quantile, soit la variance du quantile. 
  # Entrees :
  # - xnew : nouveau point
  # - model : un modele de krigeage
  # - data.x : des precalculs aux points d'intégration
  # - x : les points d'intégration (vecteur, matrice ou data.frame)
  #-----------------------------------------------------------#
  if (is.null(dim(x)))    x    <- matrix(x, ncol=1)
  if (is.null(dim(xnew))) xnew <- matrix(xnew, nrow=1)
  
  if (checkPredict(x=xnew, model=list(model))){
    return(NA)
  } else {
    
    # Get precalculations
    m.x <- data.x$mean
    s.x <- data.x$sd
    precalc.data.x <- data.x$precalc.data
    n.x <- length(m.x)
    k <- round(alpha*n.x)
    
    # Precalculations for xnew
    p.xnew <- predict_nobias_km(model, data.frame(x=xnew), "UK", checkNames=FALSE)
    m.xnew <- p.xnew$mean
    s.xnew <- p.xnew$sd
    F.newdata <- p.xnew$F.newdata
    c.newdata <- p.xnew$c
    kn     <- computeQuickKrigcov(model=model, integration.points=(x), X.new=data.frame(x=xnew),
                                  precalc.data=precalc.data.x, F.newdata=F.newdata, c.newdata=c.newdata)
  
    # b et a comme sur le rapport
    b <- m.x
    a <- as.numeric(kn)
    
    # intervalle à 99% pour z
    ynew.range <- c(m.xnew - beta*s.xnew, m.xnew + beta*s.xnew)
    z.range    <- (ynew.range - m.xnew) / s.xnew^2
    
    i2 <- order(b + a * z.range[1])[k]
    .J2 <- i2
    .I3 <- i3 <- z.range[1]
    while (T) {
      v <- getIntervals(a, b, i2, i3)
      #print(v)
      if (v[2] > z.range[2]) {
        break;
      }
      i3 <- (b[i2] - b[v[1]]) / (a[v[1]] - a[i2])
      .J2 <- c(.J2, i2 <- v[1])
      .I3 <- c(.I3, i3)
    }
    .I3 <- c(.I3, z.range[2])
    
    # Points correspondants (reduits ou tous)
    xQ  <- x[.J2,,drop=FALSE]
    
    # Précalculations aux points quantile
    Kinv.F <- precalc.data.x$Kinv.F
    Kinv.c.olddata <- precalc.data.x$Kinv.c.olddata[,.J2,drop=FALSE]
    first.member   <- precalc.data.x$first.member[.J2,,drop=FALSE]
    precalc.data.xQ <- list(Kinv.F=Kinv.F, Kinv.c.olddata=Kinv.c.olddata, first.member=first.member)
    m.xQ <- m.x[.J2]
    s.xQ <- m.x[.J2]
    
    # Get precalculations
    I <- .I3
    n.xQ <- length(m.xQ)
    if ((length(I) - 1) != n.xQ) {
      cat("data.xQ and I are inconsistent")
      break;
    }
    
    # Covariance
    k.xQ.xnew <- computeQuickKrigcov(model=model, integration.points=data.frame(x=xQ), X.new=data.frame(x=xnew),
                                     precalc.data=precalc.data.xQ, F.newdata=F.newdata, c.newdata=c.newdata)
    # Main loop
    varX <- EX <- PA <- a <- b <- rep(0, n.xQ)
    a <- k.xQ.xnew / s.xnew^2
    b <- m.xQ - k.xQ.xnew / s.xnew^2*m.xnew
    
    all.d <- s.xnew^2*I + m.xnew
    all.u <- (all.d-m.xnew)/s.xnew
    phi.u <- dnorm(all.u)
    Phi.u <- pnorm(all.u)
    PA <- diff(Phi.u)
    I <- which(PA>0)
    term1 <- - diff(all.u*phi.u)
    term2 <- - diff(phi.u)
    EX[I]   <- m.xnew + s.xnew * term2[I] / PA[I]
    varX[I] <- s.xnew^2 *  ( 1 +  term1[I] / PA[I] - (term2[I]/PA[I])^2 )
    
    varQ <- sum(a^2*varX*PA) + sum((b + a*EX)^2*(1-PA)*PA)
    
    if (n.xQ >1) {
      baEXPA <- matrix((b + a*EX)*PA, nrow=1)
      Q <- 2*crossprod(baEXPA)
      varQ <- varQ - sum(Q[upper.tri(Q)]) 
    }
    
    return(varQ)
  }
}