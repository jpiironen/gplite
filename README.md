

# gplite


An R package for fitting some of the most common Gaussian process (GP) models. Uses Laplace approximation for handling non-Gaussian observation models, performs hyperparameter optimization using maximum marginal likelihood (or posterior), and implements some common sparse approximations for handling larger datasets. Provides also interface for performing MCMC for the latent values using Stan.

The syntax has taken a lot of inspiration from that of [GPstuff](https://github.com/gpstuff-dev/gpstuff) but the intention of the package is *not* to be a GPstuff clone for R.


### Resources

* [Quickstart tutorial](https://htmlpreview.github.io/?https://github.com/jpiironen/gplite/blob/master/vignettes/quickstart.html) (notebook)
* [Open an issue / ask question](https://github.com/jpiironen/gplite/issues) (GitHub issues for bug reports, questions, and feature requests)


### Installation

* Install the latest release from CRAN:

```r
install.packages('gplite')
```

* Install latest development version from GitHub (requires [devtools](https://github.com/r-lib/devtools) package):

```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github('jpiironen/gplite', build_vignettes = TRUE)
```
    
### Example

```R
library(gplite)
library(ggplot2)

# create some toy 1d regression data
set.seed(32004)
n <- 200
sigma <- 0.1
x <- rnorm(n)
y <- sin(3*x)*exp(-abs(x)) +  rnorm(n)*sigma 

# set up the gp model, and optimize the hyperparameters
gp <- gp_init(cfs = cf_sexp(), lik = lik_gaussian())
gp <- gp_optim(gp, x, y)

# compute the predictive mean and variance in a grid of points
xt <- seq(-4, 4, len=300)
pred <- gp_pred(gp, xt, var=T)

# visualize
mu <- pred$mean
lb <- pred$mean - 2*sqrt(pred$var)
ub <- pred$mean + 2*sqrt(pred$var)
ggplot() + 
  geom_ribbon(aes(x=xt, ymin=lb, ymax=ub), fill='lightgray') +
  geom_line(aes(x=xt, y=mu), size=1) +
  geom_point(aes(x=x, y=y), size=0.5) +
  xlab('x') + ylab('y')
```


### References

Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. *MIT Press*. [Online](http://www.gaussianprocess.org/gpml/)



  [quickstart-vignette]: https://htmlpreview.github.io/?https://github.com/jpiironen/gplite/blob/master/vignettes/quickstart.html

