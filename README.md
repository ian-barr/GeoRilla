### Installation

This method depends on the `remotes` library:

`> remotes::install_github("ian-barr/geoDude", build_manuals=TRUE,INSTALL_opts = c("--html"))`


### Working with Geometric Statistics

Geometric means, Standard deviations, 95% Confidence intervals, and SEMs are defined based on the log-transformed data.
If  we have a set of measurements $X_1, X_2, \dots X_n$, 

The geometric mean is $\exp{(\mu_l)}$ where 

$$\mu_l = \frac{1}{n}\sum_i^n  \ln(X_i)$$

The geometric standard deviation (GSD) is $\exp{(\sigma_l)}$ where

$$\sigma_l = \sqrt{\frac{\sum_i^n \left ( \ln(X_i) - \mu_l \right )^2}{n-1} }$$

The geometric standard error of the mean (GSEM) is $\exp{(sem_l)}$ where 

$$sem_l = \frac{\sigma_l}{\sqrt{n}}$$

For a confidence level of 95% and a number of datapoints $n$, the confidence interval is 
given by multiplying the $sem_l$ by the 97.5% quantile of a $t$-distribution with $df= n-1$, 
then exponentiating to $\exp{\text{(lower, upper)}}$, where (lower, upper) are :

$$\text{(lower, upper)} = \mu_l \pm t(n-1, 0.975)\cdot sem_l$$
   
This means that the ends of the confidence intervals are
   
$$\exp(\mu_l - t(n-1, 0.975)\cdot sem_l), ~ \exp(\mu_l + t(n-1, 0.975)\cdot sem_l)$$
   
Dividing the upper limit by the lower limit gives
   
$$\frac{\exp(\mu_l + t(n-1, 0.975)\cdot sem_l)}{\exp(\mu_l - t(n-1, 0.975)\cdot sem_l)} = \exp(2 \cdot t(n-1, 0.975)\cdot sem_l)$$
   
Calling that result $R$, we will be able to get an estimate of the $sem_l$:

$$sem_l = \frac{\ln \left (\sqrt{R} \right )} {t(n-1, 0.975)}$$
   
Including the definition for $sem_l$ lets us solve for $\sigma_l$:
 
$$\sigma_l = \frac{\sqrt{n}\cdot \ln \left (\sqrt{R} \right )} {t(n-1, 0.975)}$$
  
And since the geometric coefficient of variation (GCV) is $\text{GCV} = \sqrt{\exp{(\sigma_l^2)} ~– 1}$
We are able to solve for GCV, and Convert between the different statistics.

### Using the Method of Moments

The Method of Moments is a technique for estimating the parameters of a distribution based on the moments of the distribution, typically the first and second (mean and variance). 

In order to estimate the geometric mean, SD, CIs, etc, a few assumptions will be made:

1. that the underlying distribution is log-normal

2. That the variance is sufficiently small so that the overall error is small.

The median of a lognormal distribution is $\exp(\mu)$, where $\mu$ is the log-mean of the distribution. This means that for a lognormal distribution, the geometric mean and the median are equal. The lognormal distribution has two parameters $\mu$ and $\sigma$, and is defined over the interval $(0,\infty)$. The mean of the distribution is 

$$   m = E[x] =  \exp \left (\mu + \frac{\sigma^2}{2} \right ) $$

and the Variance is 

$$   V = E \left [(x- E[x])^2 \right ] =  \left ( \exp \left ( \sigma^2 \right ) -1 \right ) \cdot \exp \left ( 2\mu + \sigma^2 \right ) = \left ( \exp \left ( \sigma^2 \right ) -1 \right ) \cdot m^2$$

By solving for $\mu$ and $\sigma$ in the above equations, we are able to get estimates for the geometric mean and standard deviation:

$$\textrm{GMEAN}(m,V) = \exp{(\mu)}; ~\mu = 2 \ln(m) - \ln \left (\sqrt{V + m^2} \right )$$

$$\textrm{GSD}(m,V) = \exp(\sigma); ~\sigma = \sqrt{\ln(1 + V/m^2)}$$

So we are able to convert reported arithmetic means and standard deviations into estimates of their geometric counterparts. There is, however, some associated error. 


### Functions Defined:

`GMEAN` : Calculates Geometric Mean of a set of data points.

`GSD` : Calculates Geometric Standard Deviation of a set of data points.

`GSEM` : Calculates Geometric Standard Error of the Mean of a set of data points.

`GCV` : Calculates Geometric Coefficient of Variation of a set of data points.

`GCI` : Calculates Geometric Confidence Interval of a set of data points.

`GCV2GCI` : Converts Geometric Coefficient of Variation (expressed as a percentage if `percent = TRUE`) to Confidence interval.

`kGCV2GSD` : Converts Geometric Coefficient of Variation (expressed as a percentage if `percent = TRUE`) to Confidence interval using the definition of Kirkwood (1979) *Biometrics*.

`GCV2GSD` : Converts Geometric Coefficient of Variation (expressed as a percentage if `percent = TRUE`) to Geometric Standard Deviation.

`GCI2GCV` : Converts Geometric Confidence Interval to Geometric Coefficient of Variation (expressed as a percentage if `percent = TRUE`).

`GCV2GSEM` : Converts Geometric Coefficient of Variation (expressed as a percentage if `percent = TRUE`) to Geometric Standard Error of the Mean

`GSD2GSEM` : Converts Geometric Standard Deviation to Geometric Standard Error of the Mean


