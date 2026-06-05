
#' @title Power Calculation for F.beta Comparison Between Two Correlated Classifiers
#' @description A function to calculate power for testing the difference in F.beta scores between two correlated classifiers.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot pnorm qnorm qt


#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param Fb1 F.beta score of 1st classifier.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param Fb2 F.beta score of 2nd classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param rho.N Correlation of two classifiers' predictions
#' among negative cases. rho.N will be estimated if rho.N = NULL.
#' @param rho.P Correlation of two classifiers' predictions
#' among positive cases. rho.P will be estimated if rho.P = NULL.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = c("auto", "TRUE", "FALSE").
#' @param alpha.level alpha level.

#' @examples
#' \dontrun{
#' powerFbTwo.rho(N=100, s=40,
#'                 Fb1=NULL, pp1=0.6, ps1=0.8,
#'                 Fb2=NULL, pp2=0.7, ps2=0.9,
#'                 rho.N=0.2, rho.P=0.2,
#'                 beta=1,
#'                 alternative=c("two.sided", "less", "greater"),
#'                 normal.approx=c("auto", "TRUE", "FALSE"),
#'                 alpha.level=0.05)
#' }

#' @export
powerFbTwo.rho <- function(N=100, s=40,
                           Fb1=NULL, pp1=0.6, ps1=0.8,
                           Fb2=NULL, pp2=0.7, ps2=0.9,
                           rho.N=0, rho.P=0,
                           beta=1,
                           alternative=c("two.sided", "less", "greater"),
                           normal.approx=c("auto", "TRUE", "FALSE"),
                           alpha.level=0.05) {

  alternative <- match.arg(alternative)
  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb1)) Fb1 <- one.plus.beta2*pp1*ps1/(beta2*pp1+ps1)
  if(is.null(Fb2)) Fb2 <- one.plus.beta2*pp2*ps2/(beta2*pp2+ps2)

  Fb1.beta2 <- Fb1*beta2
  Fb2.beta2 <- Fb2*beta2


  z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
  z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))

  if(is.null(pp2) & is.null(ps2)) {

    if(!is.null(pp1)) {
      xx0 <- pp1
    } else if(is.null(pp1) & !is.null(ps1)) {
      xx0 <- Fb1.beta2*ps1/(one.plus.beta2*ps1 - Fb1)
    } else {
      pp.low <- min(round(Fb1/(one.plus.beta2-Fb1.beta2), 1) + 0.05, 1)
      if(Fb1 <  0.5) {
        xx0 <- seq(pp.low, 0.5, 0.05)
      } else {
        xx0 <- seq(pp.low, 1, 0.05)
      }
    }

    allf1 <- lapply(xx0, function(x) {
      pp1 <- x
      ps1 <- Fb1.beta2*x/(one.plus.beta2*x - Fb1)
      list(pp1=pp1, ps1=ps1)
    })
    allf1 <- Filter(Negate(is.null), allf1)


    pp.low <- min(round(Fb2/(one.plus.beta2-Fb2.beta2), 1) + 0.05)
    if(Fb2 <  0.5) {
      xx <- seq(pp.low, 0.5, 0.05)
    } else {
      xx <- seq(pp.low, 1, 0.05)
    }

    yy <- sapply(xx, function(x) {
      pp2 <- x
      ps2 <- Fb2.beta2*x/(one.plus.beta2*x - Fb2)

      res <- try({

        pw.s <- sapply(allf1, function(u) {

          pp1 <- u$pp1
          ps1 <- u$ps1

          ## Null distribution
          pp.bar <- (pp1+pp2)/2
          ps.bar <- (ps1+ps2)/2

          pa.bar <- paCalculation(N, s, pp.bar, ps.bar, beta)$pa

          #if(normal.approx=="auto" & length(xf112n) > 500) normal.approx <- "TRUE"
          #if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
          #   s*min(ps.bar, 1-ps.bar) >= 5 | (N-s) > 500)  ) normal.approx <- "TRUE"
          if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
             s*min(ps.bar, 1-ps.bar) >= 5)  ) normal.approx <- "TRUE"

          if(normal.approx=="TRUE") {

            outn <- dFb_variance_approx(N=N, s=s,
                                        pp1=pp.bar, ps1=ps.bar,
                                        pp2=pp.bar, ps2=ps.bar,
                                        rho.N=rho.N, rho.P=rho.P,
                                        beta=beta)

            outa <- dFb_variance_approx(N=N, s=s,
                                       pp1=pp1, ps1=ps1,
                                       pp2=pp2, ps2=ps2,
                                       rho.N=rho.N, rho.P=rho.P,
                                       beta=beta)

            dmean.bar.00 <- 0
            dvar.bar.00 <- outn$Varw

            dmean.bar.01 <- outa$Ew
            dvar.bar.01 <- outa$Varw

            if(alternative=="greater") {
              cr <-   z.alpha * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
            } else if(alternative=="less") {
              cr <- - z.alpha * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
            } else {
              cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
              cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
                pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

            }

          } else {

            ## Alternative distribution
            F1dist12a <- F1.cond.two.b.rho(N, s,
                                           pp1=pp1, ps1=ps1,
                                           pp2=pp2, ps2=ps2,
                                           rho.N=rho.N, rho.P=rho.P,
                                           beta)

            F1dist12a <- F1dist12a[F1dist12a$pf1>1E-10,]
            yf112a <- F1dist12a$pf1
            yf112a <- yf112a/sum(yf112a)
            xf112a <- F1dist12a$f1s

            ## Null distribution
            #pp.bar <- (pp1+pp2)/2
            #ps.bar <- (ps1+ps2)/2

            F1dist12n <- F1.cond.two.b.rho(N, s,
                                           pp1=pp.bar, ps1=ps.bar,
                                           pp2=pp.bar, ps2=ps.bar,
                                           rho.N=rho.N, rho.P=rho.P,
                                           beta)

            F1dist12n <- F1dist12n[F1dist12n$pf1>1E-10,]
            yf112n <- F1dist12n$pf1
            yf112n <- yf112n/sum(yf112n)
            xf112n <- F1dist12n$f1s

            ## null distribution
            f1s.null <- xf112n
            pf1.null <- yf112n

            ## alternative distribution
            f1s.alt <- xf112a
            pf1.alt <- yf112a

            cumsum.pf1 <- cumsum(pf1.null)
            cumsum.rev.pf1 <- cumsum(rev(pf1.null))

            id.left <- findInterval(f1s.alt, f1s.null)
            id.right <- findInterval(-f1s.alt, -rev(f1s.null))

            if(alternative=="greater") {
              pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
            } else if(alternative=="less") {
              pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
            } else {
              prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
              prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

              C.obs.vec <- pmin(prob.left, prob.right)

              id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
              id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

              pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
              pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

              pvalue <- pv.left + pv.right
            }

            reject1 <- pvalue<alpha.level
            pw.v <- sum(pf1.alt[reject1])
          }

          c(pw.v, u$pp1, u$ps1)

        })

        min.pw.idx <- order(pw.s[1,])
        c(pw.s[,min.pw.idx[1]], Fb1, pp2, ps2)

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(rep(NA, 6))
      } else {
        return(res)
      }

    })

    min.pw.idx <- order(yy[1,])
    out <- yy[,min.pw.idx[1]]

  } else {
    if(is.null(pp2) & !is.null(ps2)) {
      pp2 <- Fb2.beta2*ps2/(one.plus.beta2*ps2 - Fb2)
    } else if(!is.null(pp2) & is.null(ps2)) {
      ps2 <- Fb2.beta2*pp2/(one.plus.beta2*pp2 - Fb2)
    }
    pw <- Power_Two_beta_rho_1(N=N, s=s, pp2=pp2, ps2=ps2,
                               Fb1=Fb1, pp1=pp1, ps1=ps1,
                               rho.N=rho.N, rho.P=rho.P,
                               beta=beta,
                               alternative=alternative,
                               normal.approx=normal.approx,
                               alpha.level=alpha.level)
    out <- pw
  }

  out <- c(out, Fb2, rho.N, rho.P)
  names(out) <- c("Power", "pp1", "ps1", "Fb1", "pp2", "ps2", "Fb2", "rho.N", "rho.P")
  return(round(out, 5))

}



#Power_Two_beta_rho_1(N=50, s=20,
#                     pp2=0.9, ps2=0.9,
#                     Fb1=0.783, pp1=NULL, ps1=NULL,
#                     rho.N=0.3, rho.P=0.3,
#                     beta=1,
#                     alternative=c("two.sided", "less", "greater"),
#                     normal.approx=c("auto", "TRUE", "FALSE"),
#                     alpha.level=0.05)

### Not export
Power_Two_beta_rho_1 <- function(N=50, s=20,
                                 pp2=0.9, ps2=0.9,
                                 Fb1=0.783, pp1=NULL, ps1=NULL,
                                 rho.N=0.1, rho.P=0.1,
                                 beta=1,
                                 alternative=c("two.sided", "less", "greater"),
                                 normal.approx=c("auto", "TRUE", "FALSE"),
                                 alpha.level=0.05) {

  alternative <- match.arg(alternative)
  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb1)) Fb1 <- one.plus.beta2*pp1*ps1/(beta2*pp1+ps1)

  Fb1.beta2 <- Fb1*beta2

  z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
  z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))

  if(is.null(pp1) & is.null(ps1)) {

    pp.low <- min(round(Fb1/(one.plus.beta2-Fb1.beta2), 1) + 0.05, 1)
    if(Fb1 <  0.5) {
      xx <- seq(pp.low, 0.5, 0.05)
    } else {
      xx <- seq(pp.low, 1, 0.05)
    }

    yy <- sapply(xx, function(x) {

      res <- try({
        pp1 <- x
        ps1 <- Fb1.beta2*x/(one.plus.beta2*x - Fb1)

        ## Null distribution
        pp.bar <- (pp1+pp2)/2
        ps.bar <- (ps1+ps2)/2

        pa.bar <- paCalculation(N, s, pp.bar, ps.bar, beta)$pa

        #if(normal.approx=="auto" & length(xf112n) > 500) normal.approx <- "TRUE"
        #if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
        #   s*min(ps.bar, 1-ps.bar) >= 5 | (N-s) > 500) ) normal.approx <- "TRUE"
        if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
           s*min(ps.bar, 1-ps.bar) >= 5) ) normal.approx <- "TRUE"

        if(normal.approx=="TRUE") {

          outn <- dFb_variance_approx(N=N, s=s,
                                      pp1=pp.bar, ps1=ps.bar,
                                      pp2=pp.bar, ps2=ps.bar,
                                      rho.N=rho.N, rho.P=rho.P,
                                      beta=beta)

          outa <- dFb_variance_approx(N=N, s=s,
                                      pp1=pp1, ps1=ps1,
                                      pp2=pp2, ps2=ps2,
                                      rho.N=rho.N, rho.P=rho.P,
                                      beta=beta)

          dmean.bar.00 <- 0
          dvar.bar.00 <- outn$Varw

          dmean.bar.01 <- outa$Ew
          dvar.bar.01 <- outa$Varw

          if(alternative=="greater") {
            cr <-   z.alpha * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
          } else if(alternative=="less") {
            cr <- - z.alpha * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
          } else {
            cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
            cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
              pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

          }

        } else {

          ## Alternative distribution
          F1dist12a <- F1.cond.two.b.rho(N, s,
                                         pp1=pp1, ps1=ps1,
                                         pp2=pp2, ps2=ps2,
                                         rho.N=rho.N, rho.P=rho.P,
                                         beta)

          F1dist12a <- F1dist12a[F1dist12a$pf1>1E-10,]
          yf112a <- F1dist12a$pf1
          yf112a <- yf112a/sum(yf112a)
          xf112a <- F1dist12a$f1s

          ## Null distribution
          #pp.bar <- (pp1+pp2)/2
          #ps.bar <- (ps1+ps2)/2

          F1dist12n <- F1.cond.two.b.rho(N, s,
                                         pp1=pp.bar, ps1=ps.bar,
                                         pp2=pp.bar, ps2=ps.bar,
                                         rho.N=rho.N, rho.P=rho.P,
                                         beta)

          F1dist12n <- F1dist12n[F1dist12n$pf1>1E-10,]
          yf112n <- F1dist12n$pf1
          yf112n <- yf112n/sum(yf112n)
          xf112n <- F1dist12n$f1s

          ## null distribution
          f1s.null <- xf112n
          pf1.null <- yf112n

          ## alternative distribution
          f1s.alt <- xf112a
          pf1.alt <- yf112a

          cumsum.pf1 <- cumsum(pf1.null)
          cumsum.rev.pf1 <- cumsum(rev(pf1.null))

          id.left <- findInterval(f1s.alt, f1s.null)
          id.right <- findInterval(-f1s.alt, -rev(f1s.null))

          if(alternative=="greater") {
            pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
          } else if(alternative=="less") {
            pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
          } else {
            prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
            prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

            C.obs.vec <- pmin(prob.left, prob.right)

            id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
            id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

            pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
            pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

            pvalue <- pv.left + pv.right
          }

          reject1 <- pvalue<alpha.level
          pw.v <- sum(pf1.alt[reject1])
        }

        c(pw.v, pp1, ps1)

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(rep(NA, 3))
      } else {
        return(res)
      }

    })
    min.pw.idx <- order(yy[1,])
    out <- yy[,min.pw.idx[1]]

  } else {

    if(is.null(pp1) & !is.null(ps1)) {
      pp1 <- Fb1.beta2*ps1/(one.plus.beta2*ps1 - Fb1)
    } else if(!is.null(pp1) & is.null(ps1)) {
      ps1 <- Fb1.beta2*pp1/(one.plus.beta2*pp1 - Fb1)
    }
    pw <- Power_Two_beta_rho(N=N, s=s,
                             pp1=pp1, ps1=ps1,
                             pp2=pp2, ps2=ps2,
                             rho.N=rho.N, rho.P=rho.P,
                             beta=beta,
                             alternative=alternative,
                             normal.approx=normal.approx,
                             alpha.level=alpha.level)
    out <- pw[1:3]
  }

  out <- c(out, Fb1, pp2, ps2)
  names(out) <- c("Power", "pp1", "ps1", "Fb1", "pp2", "ps2")
  return(out)

}


#Power_Two_beta_rho(N=100, s=40,
#                   pp1=9/11, ps1=0.75,
#                   pp2=0.9, ps2=0.9,
#                   rho.N=0.3, rho.P=0.3,
#                   beta=1,
#                   alternative=c("two.sided", "less", "greater"),
#                   normal.approx=c("auto", "TRUE", "FALSE"),
#                   alpha.level=0.05)

### Not export
Power_Two_beta_rho <- function(N=50, s=20,
                               pp1=9/11, ps1=0.75,
                               pp2=0.9, ps2=0.9,
                               rho.N=0.1, rho.P=0.1,
                               beta=1,
                               alternative=c("two.sided", "less", "greater"),
                               normal.approx=c("auto", "TRUE", "FALSE"),
                               alpha.level=0.05) {

  alternative <- match.arg(alternative)
  normal.approx <- match.arg(normal.approx)

  ## Null distribution
  pp.bar <- (pp1+pp2)/2
  ps.bar <- (ps1+ps2)/2

  pa.bar <- paCalculation(N, s, pp.bar, ps.bar, beta)$pa

  #if(normal.approx=="auto" & length(xf112n) > 500) normal.approx <- "TRUE"
  #if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
  #   s*min(ps.bar, 1-ps.bar) >= 5 | (N-s) > 500) ) normal.approx <- "TRUE"
  if(normal.approx=="auto" & ((N-s)*min(pa.bar, 1-pa.bar) >= 5 &
     s*min(ps.bar, 1-ps.bar) >= 5) ) normal.approx <- "TRUE"

  if(normal.approx=="TRUE") {

    outn <- dFb_variance_approx(N=N, s=s,
                                pp1=pp.bar, ps1=ps.bar,
                                pp2=pp.bar, ps2=ps.bar,
                                rho.N=rho.N, rho.P=rho.P,
                                beta=beta)

    outa <- dFb_variance_approx(N=N, s=s,
                                pp1=pp1, ps1=ps1,
                                pp2=pp2, ps2=ps2,
                                rho.N=rho.N, rho.P=rho.P,
                                beta=beta)

    dmean.bar.00 <- 0
    dvar.bar.00 <- outn$Varw

    dmean.bar.01 <- outa$Ew
    dvar.bar.01 <- outa$Varw

    z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
    z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))

    if(alternative=="greater") {
      cr <-   z.alpha * sqrt(dvar.bar.00)
      pw <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
    } else if(alternative=="less") {
      cr <- - z.alpha * sqrt(dvar.bar.00)
      pw <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
    } else {
      cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
      cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
      pw <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
        pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

    }

  } else {

    F1dist12n <- F1.cond.two.b.rho(N, s,
                                   pp1=pp.bar, ps1=ps.bar,
                                   pp2=pp.bar, ps2=ps.bar,
                                   rho.N=rho.N, rho.P=rho.P,
                                   beta)

    F1dist12n <- F1dist12n[F1dist12n$pf1>1E-10,]
    yf112n <- F1dist12n$pf1
    yf112n <- yf112n/sum(yf112n)
    xf112n <- F1dist12n$f1s


    ## Alternative distribution
    F1dist12a <- F1.cond.two.b.rho(N, s,
                                   pp1=pp1, ps1=ps1,
                                   pp2=pp2, ps2=ps2,
                                   rho.N=rho.N, rho.P=rho.P,
                                   beta)

    F1dist12a <- F1dist12a[F1dist12a$pf1>1E-10,]
    yf112a <- F1dist12a$pf1
    yf112a <- yf112a/sum(yf112a)
    xf112a <- F1dist12a$f1s

    ## null distribution
    f1s.null <- xf112n
    pf1.null <- yf112n

    ## alternative distribution
    f1s.alt <- xf112a
    pf1.alt <- yf112a

    cumsum.pf1 <- cumsum(pf1.null)
    cumsum.rev.pf1 <- cumsum(rev(pf1.null))

    id.left <- findInterval(f1s.alt, f1s.null)
    id.right <- findInterval(-f1s.alt, -rev(f1s.null))

    if(alternative=="greater") {
      pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
    } else if(alternative=="less") {
      pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
    } else {
      prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
      prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

      C.obs.vec <- pmin(prob.left, prob.right)

      id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
      id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

      pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
      pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

      pvalue <- pv.left + pv.right
    }

    reject1 <- pvalue<alpha.level
    pw <- sum(pf1.alt[reject1])

  }

  out <- c(pw, pp1, ps1, pp2, ps2)
  names(out) <- c("Power", "pp1", "ps1", "pp2", "ps2")
  return(out)

}

