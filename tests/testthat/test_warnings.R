
## test warnings about multiple magnitude parameters in a product kernel

# kernels with non-fixed magnitude
cfs1 <- list(
  cf_const(),
  cf_lin(),
  cf_sexp(),
  cf_matern32(),
  cf_matern52(),
  cf_nn(),
  cf_periodic()
)

# kernels with fixed magnitude
cfs2 <- list(
  cf_const(prior_magn = prior_fixed()),
  cf_lin(prior_magn = prior_fixed()),
  cf_sexp(prior_magn = prior_fixed()),
  cf_matern32(prior_magn = prior_fixed()),
  cf_matern52(prior_magn = prior_fixed()),
  cf_nn(prior_magn = prior_fixed()),
  cf_periodic(cf_base=cf_sexp(prior_magn = prior_fixed()))
)

# all pairs of covariance functions with non-fixed magnitude
cf_comb <- combn(cfs1, 2)
for (i in 1:NCOL(cf_comb)) {
  testthat::expect_warning(
    prod <- cf_prod(cf_comb[[1,i]], cf_comb[[2,i]]),
    regexp = "more than one non-fixed magnitude parameters"
  )
  testthat::expect_warning(
    prod <- cf_prod(cf_comb[[2,i]], cf_comb[[1,i]]),
    regexp = "more than one non-fixed magnitude parameters"
  )
  testthat::expect_warning(
    prod <- cf_comb[[1,i]] * cf_comb[[2,i]],
    regexp = "more than one non-fixed magnitude parameters"
  )
  testthat::expect_warning(
    prod <- cf_comb[[2,i]] * cf_comb[[1,i]],
    regexp = "more than one non-fixed magnitude parameters"
  )
}

# all pairs, so that one of them has a fixed magnitude
for (i in seq_along(cfs1)) {
  for (j in seq_along(cfs2)) {
    
    
    testthat::expect_silent(
      prod <- cf_prod(cfs1[[i]], cfs2[[j]])
    )
    testthat::expect_silent(
      prod <- cf_prod(cfs2[[j]], cfs1[[i]])
    )
    testthat::expect_silent(
      prod <- cfs1[[i]] * cfs2[[j]]
    )
    testthat::expect_silent(
      prod <- cfs2[[j]] * cfs1[[i]]
    )
    
    # check combinations of three
    for (k in seq_along(cfs1)) {
      cf_temp <- cfs1[[i]] * cfs2[[j]]
      
      testthat::expect_silent(
        prod <- cf_prod(cf_temp, cfs2[[k]])
      )
      testthat::expect_silent(
        prod <- cf_prod(cfs2[[k]], cf_temp)
      )
      testthat::expect_silent(
        prod <- cf_temp * cfs2[[j]]
      )
      testthat::expect_silent(
        prod <- cfs2[[j]] * cf_temp 
      )
      
      testthat::expect_warning(
        prod <- cf_prod(cf_temp, cfs1[[k]]),
        regexp = "more than one non-fixed magnitude parameters"
      )
      testthat::expect_warning(
        prod <- cf_prod(cfs1[[k]], cf_temp),
        regexp = "more than one non-fixed magnitude parameters"
      )
      testthat::expect_warning(
        prod <- cfs1[[k]] * cf_temp,
        regexp = "more than one non-fixed magnitude parameters"
      )
      testthat::expect_warning(
        prod <- cf_temp * cfs1[[k]],
        regexp = "more than one non-fixed magnitude parameters"
      )
      
    }
  }
}

