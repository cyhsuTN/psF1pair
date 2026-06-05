
#' @title Confidence Interval for F.beta Difference
#' @description A function to calculate the CI of the F.beta difference between
#' two correlated classifiers.
#' @importFrom utils tail
#' @importFrom stats qnorm

#' @param data A data frame (columns: actual outcomes,
#' C1 predicted outcomes, C2 predicted outcomes).
#' Encoded as 0 (negative) and 1 (positive).
#' @param rho.N Correlation of two classifiers' predictions
#' among negative cases. rho.N will be estimated if rho.N = NULL.
#' @param rho.P Correlation of two classifiers' predictions
#' among positive cases. rho.P will be estimated if rho.P = NULL.
#' @param beta beta.
#' @param alpha Significance level.
#' @param normal.approx normal.approx = c("auto", "TRUE", "FALSE").

#' @examples
#' \dontrun{
#' set.seed(2028)
#' data <- data_Generation(N=100, s=40,
#'                         pp1=0.6, ps1=0.8,
#'                         pp2=0.7, ps2=0.9,
#'                         rho.N=0.2, rho.P=0.2,
#'                         beta=1)
#' ciFbTwo.rho(data)
#' ciFbTwo.rho(data, rho.N=0, rho.P=0) # wider
#' }

#' @export
ciFbTwo.rho <- function(data,
                        rho.N=NULL,
                        rho.P=NULL,
                        beta=1,
                        alpha=0.05,
                        normal.approx=c("auto", "TRUE", "FALSE")) {

  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  s <- sum(data[,1]==1)
  N <- nrow(data)

  d1 <- sum(data[,1]==1 & data[,2]==1)
  d2 <- sum(data[,1]==1 & data[,3]==1)

  b1 <- sum(data[,1]==0 & data[,2]==1)
  b2 <- sum(data[,1]==0 & data[,3]==1)

  if (b1==0) b1 <- 0.5
  if (b2==0) b2 <- 0.5

  pp1.hat <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
  ps1.hat <- (d1+0.5)/(s+1)

  pp2.hat <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
  ps2.hat <- (d2+0.5)/(s+1)

  rho.N.hat <- rho.N
  rho.P.hat <- rho.P

  if (is.null(rho.N.hat) | is.null(rho.P.hat)) {
    #### Pearson correlation
    #rho.N.hat <- cor(x=data[data[,1]==0,2], y=data[data[,1]==0,3])
    #rho.P.hat <- cor(x=data[data[,1]==1,2], y=data[data[,1]==1,3])

    #### Posterior Pearson correlation (Jeffrey’s non-informative prior)
    rho.N.hat <- bayes_phi(x=data[data[,1]==0,2],
                           y=data[data[,1]==0,3])$summary$phi_mean
    rho.P.hat <- bayes_phi(x=data[data[,1]==1,2],
                           y=data[data[,1]==1,3])$summary$phi_mean

    #p11 <- sum(x == 1 & y == 1)/length(x)  # p11
    #p10 <- sum(x == 1 & y == 0)/length(x)  # p10
    #p01 <- sum(x == 0 & y == 1)/length(x)  # p01
    #p00 <- sum(x == 0 & y == 0)/length(x)  # p00
    #px1 <- p11 + p10
    #px0 <- p01 + p00
    #py1 <- p11 + p01
    #py0 <- p10 + p00
    #num <- p11 * p00 - p10 * p01
    #den <- sqrt(px1 * px0 * py1 * py0)
    #(phi <- num / den)  # equal Pearson correlation
  }


  fb1.raw <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2.raw <- one.plus.beta2*d2/(d2 + b2 + s*beta2)
  fb1 <- one.plus.beta2*pp1.hat*ps1.hat/(ps1.hat + pp1.hat*beta2)
  fb2 <- one.plus.beta2*pp2.hat*ps2.hat/(ps2.hat + pp2.hat*beta2)

  dfb.raw <- fb1.raw - fb2.raw
  dfb <- fb1 - fb2

  pa1 <- paCalculation(N, s, pp1.hat, ps1.hat, beta)$pa
  pa2 <- paCalculation(N, s, pp2.hat, ps2.hat, beta)$pa

  if(normal.approx=="auto" & (N-s)*min(pa1, 1-pa1, pa2, 1-pa2) >= 5 &
     s*min(ps1.hat, 1-ps1.hat, ps2.hat, 1-ps2.hat) >= 5 ) normal.approx <- "TRUE"

  if(normal.approx=="TRUE") {

    out <- dFb_variance_approx(N=N, s=s,
                               pp1=pp1.hat, ps1=ps1.hat,
                               pp2=pp2.hat, ps2=ps2.hat,
                               rho.N=rho.N.hat, rho.P=rho.P.hat,
                               beta=beta)

    dmean <- out$Ew
    dvar <- out$Varw

    dFb.L <- dmean - qnorm(1-alpha/2)*sqrt(dvar)

    dFb.U <- dmean + qnorm(1-alpha/2)*sqrt(dvar)

    dFb.vec <- c(dfb, dFb.L, dFb.U, dfb.raw)

  } else {

    F1dist.diff <- F1.cond.two.b.rho(N, s,
                                     pp1=pp1.hat, ps1=ps1.hat,
                                     pp2=pp2.hat, ps2=ps2.hat,
                                     rho.N=rho.N.hat, rho.P=rho.P.hat,
                                     beta=beta)

    F1dist.diff <- F1dist.diff[F1dist.diff$pf1>1E-10,]
    pf1 <- F1dist.diff$pf1
    pf1 <- pf1/sum(pf1)
    f1s <- F1dist.diff$f1s

    C.obs <- alpha/2

    cumsum.pf1 <- cumsum(pf1)
    cumsum.rev.pf1 <- cumsum(rev(pf1))

    dFb.L <- ifelse(any(cumsum.pf1<=C.obs), tail(f1s[cumsum.pf1<=C.obs],1), 0)

    dFb.U <- ifelse(any(cumsum.rev.pf1<=C.obs), tail(rev(f1s)[cumsum.rev.pf1<=C.obs],1), 1)

    dFb.vec <- c(dfb, dFb.L, dFb.U, dfb.raw)

  }

  names(dFb.vec) <- c("dFb.est", paste0((1-alpha)*100,"%CI", c(".L", ".U")), "dFb.est.raw")

  return(c(dFb.vec, rho.N.hat=rho.N.hat, rho.P.hat=rho.P.hat))

}

