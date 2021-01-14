

inv_lemma_get <- function(K, U, D, V = NULL, logdet = TRUE) {
  # Return an object that can be used to compute inv(U * inv(K) * V + D)
  # using inversion lemma. D should be a vector giving only the diagonal.
  # If V is not given, then V = t(U).
  # The actual solution to the system of equations can be obtained by calling
  # inv_lemma_solve where the first argument is the object returned by this function.
  #
  # The inversion lemma solution is:
  # inv(D) - inv(D) * U * inv(K + V*inv(D)*U) * V * inv(D)

  # NOTE: inv_lemma_solve does not yet account for V != t(U) !

  if (is.null(V)) {
    V <- t(U)
    V_is_tU <- TRUE
  } else {
    V_is_tU <- FALSE
  }

  D <- as.vector(D)
  inv_D_U <- diag_times_dense(1 / D, U)
  V_inv_D_U <- V %*% inv_D_U
  if (V_is_tU) {
    V_inv_D_U <- 0.5 * (V_inv_D_U + t(V_inv_D_U))
  } # ensure symmetric
  L <- t(chol(K + V_inv_D_U))

  # log determinant of the matrix that is to be inverted
  if (logdet) {
    logdet <- sum(log(D)) - 2 * sum(log(diag(chol(K)))) + 2 * sum(log(diag(L)))
  } else {
    logdet <- NA
  }

  return(list(
    chol = L,
    D = D,
    inv_D_U = inv_D_U,
    logdet = logdet
  ))
}

inv_lemma_solve <- function(lemma_obj, rhs, lhs = NULL, rhs_diag = FALSE, lhs_diag = FALSE, diag_only = FALSE) {
  D <- lemma_obj$D
  inv_D_U <- lemma_obj$inv_D_U
  L <- lemma_obj$chol

  if (rhs_diag && !is.vector(rhs)) {
    stop("rhs must be vector if rhs_diag = TRUE")
  }
  if (lhs_diag && !is.vector(lhs)) {
    stop("lhs must be vector if lhs_diag = TRUE")
  }

  if (!is.null(lhs)) {
    if (rhs_diag) {
      aux_right <- diag_times_dense(t(inv_D_U), rhs, diag_right = TRUE)

      if (lhs_diag) {
        # both rhs and lhs diagonal
        if (diag_only) {
          solution <-
            lhs / D * rhs -
            colSums(t(diag_times_dense(lhs, inv_D_U)) * backsolve(t(L), forwardsolve(L, aux_right)))
        } else {
          solution <-
            diag(lhs / D * rhs) -
            diag_times_dense(lhs, inv_D_U) %*% backsolve(t(L), forwardsolve(L, aux_right))
        }
      } else {
        # rhs diagonal, lhs non-diagonal
        if (diag_only) {
          stop("not implmented yet.")
        }
        solution <-
          diag_times_dense(lhs, rhs / D, diag_right = TRUE) -
          (lhs %*% inv_D_U) %*% backsolve(t(L), forwardsolve(L, aux_right))
      }
    } else {
      aux_right <- t(inv_D_U) %*% rhs

      if (lhs_diag) {
        # rhs non-diagonal, lhs diagonal
        if (diag_only) {
          stop("not implmented yet.")
        }
        solution <-
          diag_times_dense(lhs / D, rhs) -
          diag_times_dense(lhs, inv_D_U) %*% backsolve(t(L), forwardsolve(L, aux_right))
      } else {
        # rhs and lhs both non-diagonal
        if (diag_only) {
          stop("not implmented yet.")
        }
        solution <-
          lhs %*% diag_times_dense(1 / D, rhs) -
          (lhs %*% inv_D_U) %*% backsolve(t(L), forwardsolve(L, aux_right))
      }
    }
  } else {
    if (rhs_diag) {
      # rhs diagonal
      if (diag_only) {
        stop("not implmented yet.")
      }
      aux_right <- diag_times_dense(t(inv_D_U), rhs, diag_right = TRUE)
      solution <- diag(rhs / D) - inv_D_U %*% backsolve(t(L), forwardsolve(L, aux_right))
    } else {
      # rhs either vector or matrix (but dense)
      if (diag_only) {
        stop("not implmented yet.")
      }
      aux_right <- t(inv_D_U) %*% rhs
      if (is.matrix(rhs)) {
        aux_left <- diag_times_dense(1 / D, rhs)
      } else {
        aux_left <- rhs / D
      }
      solution <- aux_left - inv_D_U %*% backsolve(t(L), forwardsolve(L, aux_right))
    }
  }

  return(solution)
}




diag_times_dense <- function(A, B, diag_right = FALSE) {
  if (diag_right) {
    if (!is.vector(B)) {
      stop("B should be a vector if diag_right=TRUE")
    }
    return(t(B * t(as.matrix(A))))
  } else {
    if (!is.vector(A)) {
      stop("A should be a vector if diag_right=FALSE")
    }
    return(A * as.matrix(B))
  }
}

is_lower_tri <- function(A, tol = 1e-12) {
  sum(abs(A[upper.tri(A)])) < tol
}

is_upper_tri <- function(A, tol = 1e-12) {
  sum(abs(A[lower.tri(A)])) < tol
}


chol_update <- function(L, x, downdate = FALSE) {
  # rank one update of lower Cholesky
  n <- length(x)
  sgn <- ifelse(downdate, -1, 1)
  for (k in 1:n) {
    if (downdate && L[k, k] < x[k]) {
      stop("The downdate leads to non-positive definite matrix.")
    }
    r <- sqrt(L[k, k]^2 + sgn * x[k]^2)
    c <- r / L[k, k]
    s <- x[k] / L[k, k]
    L[k, k] <- r
    if (k < n) {
      L[(k + 1):n, k] <- (L[(k + 1):n, k] + sgn * s * x[(k + 1):n]) / c
      x[(k + 1):n] <- c * x[(k + 1):n] - s * L[(k + 1):n, k]
    }
  }
  return(L)
}
