#' GMEAN
#'
#' This function calculates Geometric Mean from a set of data
#' \deqn{\textrm{GMEAN}(x_i) = \exp{(\mu_l)}; \mu_l = \frac{1}{n}\sum_i^n \ln(x_i)}
#' @param x Data points, as a vector-like object
#' @keywords GMEAN
#' @return GMEAN
#' @export
#' @examples
#' GMEAN((9:20)*10) ## Returns 140.729
GMEAN <- function(x){
  if (length(x) > 1){
    gmean <- exp(mean(log(x)))
    return(gmean)
  } else {
    print( "Need 2 or more data points")
  }
}

#' GMEAN.MM
#'
#' This function calculates Geometric Mean using the Method of Moments, as well as
#' the Geometric SD. The only inputs are the mean and variance of the data. The data is assumed to be log-normal.
#' \deqn{\textrm{GMEAN}(m,V) = \exp{(\mu_l)}; ~\mu_l = 2 \ln(m) - \ln \left (\sqrt{V + m^2} \right )}
#' \deqn{\textrm{GSD}(m,V) = \exp(\sigma_l); ~\sigma_l = \sqrt{\ln(1.0 + V/m^2)}}
#' @param m The mean of the data, an estimator for \eqn{E[x]}.
#' @param V The variance of the data, an estimator for \eqn{E[x^2] - E[x]^2}.
#' @keywords GMEAN.MM
#' @return GMEAN.MM
#' @export
#' @examples
#' GMEAN.MM( mean( 9:20 * 10 ), var((9:20)*10 ) ) ## Returns $g.mean 140.715, $g.sd 1.277543
GMEAN.MM <- function(m,V){
  s2 <- log(1.0 + V/m^2)
  mu <- 2.0*log(m) - log(sqrt(V + m^2))
  return(list(g.mean = exp(mu), g.sd = exp(sqrt(s2))))
}

#' MOMENTS.MM
#'
#' This function calculates arithmetic Mean as well as the variance from the Geometric mean and sd using the Method of Moments. The only inputs are the geometric mean and variance of the data. The data is assumed to be log-normal.
#' \deqn{\textrm{GMEAN}(m,V) = \exp{(\mu_l)}; ~\mu_l = 2 \ln(m) - \ln \left (\sqrt{V + m^2} \right )}
#' \deqn{\textrm{GSD}(m,V) = \exp(\sigma_l); ~\sigma_l = \sqrt{\ln(1.0 + V/m^2)}}
#' @param gmean The geometric mean of the data.
#' @param gsd The geometric SD.
#' @keywords MOMENTS.MM
#' @return MOMENTS.MM
#' @export
#' @examples
#' MOMENTS.MM( 140.0, 1.3 ) ## Returns $mean 144.9023, $var 1496.212
MOMENTS.MM <- function(gmean,gsd){
  sd <- log(gsd)
  m <- gmean*exp(0.5*sd^2)
  V <- m^2 * (exp(sd^2) - 1)
  return(list(mean = m, var = V))
}

#' GSD
#'
#' This function calculates Geometric Standard Deviation from a set of data
#' \deqn{$$GSD = \exp(\sigma_l); \sigma_l = \sqrt{\frac{\sum_i^n \left ( \ln(X_i) - \mu_l \right )^2}{n-1} }$$}
#' @param x Data points, as a vector-like object
#' @keywords GSD
#' @return GSD
#' @export
#' @examples
#' GSD((9:20)*10) ## Returns 1.295986
GSD <- function(x){
  if (length(x) > 1){
    gsd <- exp(sd(log(x)))
    return(gsd)
  } else {
    print( "Need 2 or more data points")
  }
}

#' GSEM
#'
#' This function calculates Geometric SEM from a set of data
#' \deqn{\textrm{GSEM}(x_i) = \exp{(sem_l)};~ sem_l = \frac{\sigma_l}{\sqrt{n}}}
#' @param x Data points, as a vector-like object
#' @keywords GSEM
#' @return GSEM
#' @export
#' @examples
#' GSEM((9:20)*10) ## Returns 1.077717
GSEM <- function(x){
  n <- length(x)
  if (n > 1){
    gsem <- exp(sd(log(x))/sqrt(n))
    return(gsem)
  } else {
    print( "Need 2 or more data points")
  }
}

#' GCV
#'
#' This function calculates Geometric Coefficient of Variation
#' @param @param x Data points, as a vector-like object
#' @param percent Is the GCV to be expressed as a percent? Defaults to FALSE.
#' @keywords GCV
#' @return GCV
#' @export
#' @examples
#' GCV((9:20)*10) ### Returns 0.2636907
GCV <- function(x, percent = FALSE){
  n <- length(x)
  if (n > 1){
    f <- ifelse(percent, 100.0, 1.0)
    return(  f*sqrt(expm1(sd(log(x))^2)))
  } else {
    print( "Need 2 or more data points")
  }
}

#' kGCV
#'
#' This function calculates Geometric Coefficient of Variation
#' using the formula of Kirkwood, TBL (1979) Biometrics, Vol. 35, pp. 908-909
#' @param @param x Data points, as a vector-like object
#' @param percent Is the GCV to be expressed as a percent? Defaults to FALSE.
#' @keywords kGCV
#' @return kGCV
#' @export
#' @examples
#' kGCV((9:20)*10) ### Returns 0.295986
kGCV <- function(x, percent = FALSE){
  n <- length(x)
  if (n > 1){
    f <- ifelse(percent, 100.0, 1.0)
    return(  f* expm1( sd( log(x) ) )   )
  } else {
    print( "Need 2 or more data points")
  }
}

#' GCI
#'
#' This function calculates Geometric Confidence Interval
#' @param x Data points, as a vector-like object
#' @param alpha Type I error rate (default 0.05 for 95\% CI)
#' @keywords GCI
#' @return GCI (list of two elements $lower, $upper)
#' @export
#' @examples
#' GCI((9:20)*10) # Returns $lower 119.3551 $upper 165.9306
GCI <- function(x,alpha = 0.05){
  n <- length(x)
  if (n > 1){
    lgsem <- sd(log(x))/sqrt(n)
    lmn <- mean(log(x))
    qs <- abs(qt(alpha/2. , df=(n-1),lower.tail = TRUE))
    lCI <- lgsem * qs
    return(list('lower'=exp(lmn - lCI), 'upper'=exp(lmn + lCI)))
  } else {
    print( "Need 2 or more data points")
  }
}

#' GCV2GSD
#'
#' This function converts Geometric Coefficient of Variation to Geometric SD
#' @param GCV Geometric Coefficient of Variation
#' @param percent Is the GCV expressed as a percent? Defaults to FALSE.
#' @keywords GCV, GSD
#' @return GSD
#' @export
#' @examples
#' GCV2GSD()
GCV2GSD <- function(GCV,percent=FALSE){
  if(percent){GCV <- GCV/100.}
  return ( exp( sqrt( log1p(GCV^2)) ))
}
 
#' kGCV2GSD
#'
#' This function converts Geometric Coefficient of Variation to Geometric SD,
#' using the formula of Kirkwood, TBL (1979) Biometrics, Vol. 35, pp. 908-909
#' @param GCV Geometric Coefficient of Variation
#' @param percent Is the GCV expressed as a percent? Defaults to FALSE.
#' @keywords GCV, GSD
#' @return GSD
#' @export
#' @examples
#' kGCV2GSD()
kGCV2GSD <- function(GCV,percent=FALSE){ ## kirkwood version
  f <- ifelse(percent, 100.0, 1.0)
  return ( exp (  log1p(GCV/f) ))
}

#' GSD2GSEM
#'
#' This function converts Geometric Standard Deviation to Geometric SEM
#' @param GSD Geometric Standard Deviation
#' @param n degrees of freedom (or number of data points)
#' @keywords GSD, GSEM
#' @return GSEM
#' @export
#' @examples
#' GSD2GSEM()
GSD2GSEM <- function(GSD,n){
  return ( exp(log(GSD)/sqrt(n)))
}
#' GCV2GSEM
#'
#' This function converts Geometric Coefficient of Variation to Geometric SEM
#' @param GCV Geometric Coefficient of Variation
#' @param n degrees of freedom (or number of data points)
#' @param percent Is the GCV expressed as a percent? Defaults to FALSE.
#' @keywords GCV, GSEM
#' @return GSEM
#' @export
#' @examples
#' GCV2GSEM()
GCV2GSEM <- function(GCV,n,percent=FALSE){
  f <- ifelse(percent, 100.0, 1.0)
  GCV <- GCV/f
  return ( exp(sqrt(  log1p(GCV^2)) / sqrt(n) ))
}

#' GCV2GCI
#'
#' This function converts Geometric Coefficient of Variation to Geometric Confidence Interval
#' @param GCV Geometric Coefficient of Variation
#' @param n degrees of freedom (or number of data points)
#' @param gmn geometric mean of the samples
#' @param alpha Type I error rate (default 0.05 for 95\% CI)
#' @param percent Is the GCV expressed as a percent? Defaults to FALSE.
#' @keywords GCV, GCI
#' @return GSEM
#' @export
#' @examples
#' GCV2GCI(GCV=0.25,n=10,gmn=75.0) ## Returns 62.88782, 89.44499
GCV2GCI <- function(GCV,n,gmn,alpha = 0.05,percent = FALSE){
  f <- ifelse(percent, 100.0, 1.0)
  GCV <- GCV/f
  lgsem <- sqrt( log1p(GCV^2)) / sqrt(n)
  qs <- abs(qt(alpha/2.0, df=(n-1),lower.tail = TRUE))
  CI <- lgsem * c(-qs,qs)
  return(exp(log(gmn) + CI))
}

#' GCI2GCV
#'
#' This function converts Geometric Confidence Interval to Geometric Coefficient of Variation
#' @param GCI Geometric Confidence interval, as `c(lower,upper)`
#' @param alpha Type I error rate (0.05 for 95\% GCI)
#' @param n degrees of freedom (or number of data points)
#' @param percent Is the GCV to be expressed as a percent? Defaults to FALSE.
#' @keywords GCI, GCV
#' @return GCV (based on mean of lower, upper log SDs)
#' @export
#' @examples
#' GCI2GCV(GCI=c(45.6,66.4), n=15) ## returns 0.349297 
GCI2GCV <- function(GCI,n,alpha=0.05, percent = FALSE){
  f <- ifelse(percent, 100.0, 1.0)
  fac <- sqrt(max(GCI)/min(GCI))
  q <- qt(1.0 - alpha/2. , df=(n-1),lower.tail = TRUE)
  sd <- sqrt(n)*log(fac) / q
  return(f*sqrt(expm1(sd^2)))
}
