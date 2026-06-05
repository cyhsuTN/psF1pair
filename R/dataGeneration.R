
#' @title Data Generation
#' @description A function to generate the predicted results of two classifiers.

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param rho.N Correlation of two classifiers' predictions
#' among negative cases.
#' @param rho.P Correlation of two classifiers' predictions
#' among positive cases.
#' @param beta beta.

#' @export
data_Generation <- function(N=100, s=40,
                            pp1=0.6, ps1=0.8,
                            pp2=0.7, ps2=0.9,
                            rho.N=0, rho.P=0,
                            beta=1) {

  pa1 <- paCalculation(N, s, pp1, ps1, beta)$pa
  pa2 <- paCalculation(N, s, pp2, ps2, beta)$pa

  idx.N <- sample.int(4, size = N-s, replace = TRUE, prob = pmulti(pa1, pa2, rho.N))
  idx.P <- sample.int(4, size =   s, replace = TRUE, prob = pmulti(ps1, ps2, rho.P))

  Y.N <- cbind(as.integer(idx.N <= 2),          # y1 = 1 for cats 1,2
               as.integer(idx.N %% 2 == 1))     # y2 = 1 for cats 1,3
  Y.P <- cbind(as.integer(idx.P <= 2),          # y1 = 1 for cats 1,2
               as.integer(idx.P %% 2 == 1))     # y2 = 1 for cats 1,3

  data <- data.frame(Actual=c(rep(0,N-s), rep(1,s)),
                     C=rbind(Y.N, Y.P))
  data
}

