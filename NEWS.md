# News

## gplite 0.13.0

* Add argument `tol_param` in `gp_optim`.
* Warn about overparametrized magnitude in product kernels.
* Disable EP approximation for Gaussian likelihood.
* Improve some documentation.

## gplite 0.12.0

This release improves the hyperparameter optimization convergence, and increases the stability of the Laplace iteration in some cases. New features regarding functionality to the user:

* Added prior `prior_lognormal`.
* Added argument `restarts` to `gp_optim`.
* Added argument `ref` to `gp_loo`.
* Added argument `tol` to `approx_laplace`.

## gplite 0.11.0

Initial release.
