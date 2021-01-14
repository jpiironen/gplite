hermite_polynomials <- function(n) {
  # probabilists' hermite polynomials
  coeffs <- matrix(0, nrow = n + 1, ncol = n + 1)
  coeffs[1, 1] <- 1
  if (n == 0) {
    return(coeffs)
  }

  coeffs[2, 2] <- 1
  if (n > 1) {
    for (i in 2:n) {
      coeffs[i + 1, ] <- c(0, coeffs[i, 1:n]) - (i - 1) * coeffs[i - 1, ]
    }
  }
  coeffs
}

polyeval <- function(c, x) {
  # evaluate polynomial with coeffs c at x. the coeffs are so that
  # p(x) = c[1] + c[2]*x + c[3]*x^2 + ... + c[p+1]*x^p
  #
  p <- length(c) - 1
  v <- cbind(1, stats::poly(x, degree = p, raw = TRUE, simple = TRUE)) %*% c
  as.vector(v)
}


gauss_hermite_points <- function(p) {
  coeff <- hermite_polynomials(p)
  x <- polyroot(coeff[p + 1, ])
  suppressWarnings(x <- as.numeric(x))
  weights <- factorial(p) / (p^2 * polyeval(coeff[p, ], x)^2)
  list(x = x, weights = weights)
}


gauss_hermite_points_scaled <- function(mean, scale, order = 7, ...) {
  gh <- gauss_hermite_points(order)
  nrow <- length(mean)
  x <- t(t(matrix(1, nrow = nrow, ncol = length(gh$x))) * gh$x)

  x <- x * scale + mean
  list(x = x, weights = gh$weights)
}
