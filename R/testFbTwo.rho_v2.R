
#' @title Test for Difference in F.beta Between Two Classifiers
#' @description A function to test the difference in F.beta between two classifiers.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot

#' @param data A data frame (columns: actual outcomes,
#' C1 predicted outcomes, C2 predicted outcomes).
#' Encoded as 0 (negative) and 1 (positive).
#' @param rho.N Correlation of two classifiers' predictions
#' among negative cases. rho.N will be estimated if rho.N = NULL.
#' @param rho.P Correlation of two classifiers' predictions
#' among positive cases. rho.P will be estimated if rho.P = NULL.
#' @param Fb0 Null distribution is generated using the estimated F.beta
#' if Fb0 = "est"; Null distribution is generated using all possible F.beta
#' if Fb0 = NULL.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = c("auto", "TRUE", "FALSE").

#' @examples
#' \dontrun{
#' set.seed(2028)
#' data <- data_Generation(N=100, s=40,
#'                          pp1=0.6, ps1=0.8,
#'                          pp2=0.7, ps2=0.9,
#'                          rho.N=0.2, rho.P=0.2,
#'                          beta=1)
#' testFbTwo.rho(data)
#' testFbTwo.rho(data, rho.N=0, rho.P=0) # wider
#' }

#' @export
testFbTwo.rho <- function(data,
                          rho.N=NULL,
                          rho.P=NULL,
                          Fb0="est",
                          beta=1,
                          alternative=c("two.sided", "less", "greater"),
                          normal.approx=c("auto", "TRUE", "FALSE")) {

  alternative <- match.arg(alternative)
  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  s <- sum(data[,1]==1)
  N <- nrow(data)

  d1 <- sum(data[,1]==1 & data[,2]==1)
  d2 <- sum(data[,1]==1 & data[,3]==1)

  b1 <- sum(data[,1]==0 & data[,2]==1)
  b2 <- sum(data[,1]==0 & data[,3]==1)

  fb1 <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2 <- one.plus.beta2*d2/(d2 + b2 + s*beta2)

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

  if(!is.null(Fb0) & !is.numeric(Fb0)) {

    if(is.null(d1) | is.null(b1) | is.null(d2) | is.null(b2)) {
      stop("pp.bar and ps.bar cannot be estimated")
    }

    if (b1==0) b1 <- 0.5
    if (b2==0) b2 <- 0.5

    #pp1 <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
    #ps1 <- (d1+0.5)/(s+1)
    #pp2 <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
    #ps2 <- (d2+0.5)/(s+1)

    pp1 <- d1/(d1+b1); pp2 <- d2/(d2+b2)
    ps1 <- d1/s; ps2 <- d2/s
    pp.ave <- (pp1+pp2)/2; ps.ave <- (ps1+ps2)/2

    out <- Two.Classifier.Cond.Test.beta.rho(data,
                                             pp=pp.ave, ps=ps.ave,
                                             rho.N=rho.N.hat,
                                             rho.P=rho.P.hat,
                                             beta=beta,
                                             alternative=alternative,
                                             normal.approx=normal.approx)

    names(out) <- c("fb1", "fb2", "dfb", "P-value", "Fb0.est", "pp.bar", "ps.bar")

  } else {

    if(is.null(Fb0)) {

      if(is.null(d1) | is.null(b1) | is.null(d2) | is.null(b2)) {
        print("Search for (0, 1) instead of 95% intervals of hat{fb} due to d1, b1, d2, and b2 missing")

        fb.hat <- round(0.5 * (fb1 + fb2), 2)
        Fb0all <- seq( min(0.3, round(fb.hat,1)), 0.7, 0.1)

      } else {

        if (b1==0) b1 <- 0.5
        if (b2==0) b2 <- 0.5

        #pp1 <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
        #ps1 <- (d1+0.5)/(s+1)
        #pp2 <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
        #ps2 <- (d2+0.5)/(s+1)

        pp1 <- d1/(d1+b1); pp2 <- d2/(d2+b2)
        ps1 <- d1/s; ps2 <- d2/s
        pp.ave <- (pp1+pp2)/2; ps.ave <- (ps1+ps2)/2

        F1dist0 <- F1.cond.b(N=N, s=s, pp=pp.ave, ps=ps.ave, beta)

        cumsum.pf1 <- cumsum(F1dist0$pf1)
        cumsum.rev.pf1 <- cumsum(rev(F1dist0$pf1))

        id.left <- findInterval(0.025, cumsum.pf1)
        id.right <- findInterval(0.025, cumsum.rev.pf1)

        Fb0all <- seq( max(round(F1dist0$f1[id.left]-0.05, 1), 0),
                       min(round(rev(F1dist0$f1)[id.right]+0.05, 1), 1), 0.1)

      }

    } else {

      Fb0all <- Fb0

    }

    outall <- matrix(NA, nrow=length(Fb0all), ncol=7)
    for(i in 1:length(Fb0all)) {
      Fb0 <- Fb0all[i]

      Fb0.beta2 <- Fb0*beta2

      pp.low <- min(round(Fb0/(one.plus.beta2-Fb0.beta2), 1) + 0.05, 1)
      if(Fb0 <  0.5) {
        xx <- seq(pp.low, 0.5, 0.05)
      } else {
        xx <- seq(pp.low, 1, 0.05)
      }

      yy <- sapply(xx, function(x) {
        res <- try({
          ps <- Fb0.beta2*x/(one.plus.beta2*x - Fb0)

          Two.Classifier.Cond.Test.beta.rho(data,
                                            pp=x, ps=ps,
                                            rho.N=rho.N.hat,
                                            rho.P=rho.P.hat,
                                            beta=beta,
                                            alternative=alternative,
                                            normal.approx=normal.approx)
        }, silent = T)

        if (inherits(res, "try-error")) {
          return(rep(NA, 7))
        } else {
          return(res)
        }
      })
      sup.idx <- order(yy[4,], decreasing = T)[1]
      out1 <- yy[,sup.idx]
      outall[i,] <- out1
    }

    out <- outall[which.max(outall[,4]),]
    names(out) <- c("fb1", "fb2", "dfb", "P-value", "Fb0.argsup", "pp.argsup", "ps.argsup")


  }

  round(out, 8)

}


### Not export
Two.Classifier.Cond.Test.beta.rho <- function(data,
                                              pp=9/11, ps=0.75,
                                              rho.N=0,
                                              rho.P=0,
                                              beta=1,
                                              alternative=c("two.sided", "less", "greater"),
                                              normal.approx=c("auto", "TRUE", "FALSE")) {

  alternative <- match.arg(alternative)
  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  s <- sum(data[,1]==1)
  N <- nrow(data)

  d1 <- sum(data[,1]==1 & data[,2]==1)
  d2 <- sum(data[,1]==1 & data[,3]==1)

  b1 <- sum(data[,1]==0 & data[,2]==1)
  b2 <- sum(data[,1]==0 & data[,3]==1)

  fb1 <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2 <- one.plus.beta2*d2/(d2 + b2 + s*beta2)

  Fb0 <- one.plus.beta2*pp*ps/(beta2*pp + ps)

  fb <- fb1 - fb2

  pa <- paCalculation(N, s, pp, ps, beta)$pa


  #if(normal.approx=="auto" & ((N-s)*min(pa, 1-pa) >= 5 &
  #   s*min(ps, 1-ps) >= 5 | (N-s) > 500) ) normal.approx <- "TRUE"
  if(normal.approx=="auto" & ((N-s)*min(pa, 1-pa) >= 5 &
     s*min(ps, 1-ps) >= 5) ) normal.approx <- "TRUE"

  if (fb1==0 & fb2==0 & is.na(Fb0)) {

    pvalue <- 1

  } else {

    if(normal.approx=="TRUE") {

      out <- dFb_variance_approx(N=N, s=s,
                                 pp1=pp, ps1=ps,
                                 pp2=pp, ps2=ps,
                                 rho.N=rho.N, rho.P=rho.P,
                                 beta=beta)

      dmean <- out$Ew
      dvar <- out$Varw

      if(alternative=="greater") {
        pvalue <- pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=FALSE)
      } else if(alternative=="less") {
        pvalue <- pnorm(fb, mean=dmean, sd=sqrt(dvar))
      } else {
        ppp <- pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=TRUE)
        pvalue <- 2 * ifelse(ppp<0.5, ppp,
                             pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=FALSE))
      }

    } else {

      F1dist.diff <- F1.cond.two.b.rho(N, s,
                                       pp1=pp, ps1=ps,
                                       pp2=pp, ps2=ps,
                                       rho.N=rho.N, rho.P=rho.P,
                                       beta=beta)

      F1dist.diff <- F1dist.diff[F1dist.diff$pf1>1E-10,]
      pf1 <- F1dist.diff$pf1
      pf1 <- pf1/sum(pf1)
      f1s <- F1dist.diff$f1s

      if(alternative=="greater") {
        pvalue <- sum(pf1[f1s>fb]) + 0.5*sum(pf1[f1s=fb])
      } else if(alternative=="less") {
        pvalue <- sum(pf1[f1s<fb]) + 0.5*sum(pf1[f1s=fb])
      } else {
        C.obs <- min(sum(pf1[f1s<=fb]), sum(pf1[f1s>=fb]))

        cumsum.pf1 <- cumsum(pf1)
        cumsum.rev.pf1 <- cumsum(rev(pf1))

        pvalue <- as.numeric(
          ifelse(!any(cumsum.rev.pf1<=C.obs), 0,
                 tail(cumsum.rev.pf1[cumsum.rev.pf1<=C.obs],1)) +
            ifelse(!any(cumsum.pf1<=C.obs), 0,
                   tail(cumsum.pf1[cumsum.pf1<=C.obs],1)))
        pvalue <- min(pvalue, 1)
      }

    }

  }

  c(fb1=fb1, fb2=fb2, dfb=fb, 'P-value'=pvalue, Fb0=Fb0, pp=pp, ps=ps)

}
#Two.Classifier.Cond.Test.beta.rho(data)





