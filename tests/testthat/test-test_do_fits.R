df <- data.frame(
  values = c(21127540944,10437524907,14063403068)
)


big_data <- readRDS("testdata/raw_data.RDS")

na_df <- data.frame(
  values = rep(NA, 3)
)

test_that("do_fits returns a data.frame", {
  expect_true("data.frame" %in% class(do_fits(df)))
})

test_that("do_fits contains only correct colnames", {
  expect_true(all(colnames(do_fits(df) %in% c("dist", "loglik", "aic"))))
  expect_length(unique(colnames(do_fits(df))), 3)
})


test_that("parse_fitdist has correct colnames", {
  fit <- fitdistrplus::fitdist(df[,1], "nbinom")
  fit <- parse_fitdist(fit)
  expect_true(all(colnames(fit) %in% c("dist", "loglik", "aic")))
})

test_that("parse_univariateML has correct colnames", {
  fit <- univariateML::mlgamma(df[,1], na.rm=TRUE)
  fit <- parse_univariateML(fit)
  expect_true(all(colnames(fit) %in% c("dist", "loglik", "aic")))
})

test_that("single_fit returns distname on NA inputs", {
  fit <- single_fit(na_df[,1], "gamma")
  expect_equal("gamma", fit$dist)
})

test_that("squash_fits returns tibbles for each input ", {
  fits_uml <- single_fit(df[,1], "gamma")
  sq_uml <- squash_fits(fits_uml)
  expect_true(tibble::is_tibble(sq_uml))
  fits_fdist <- fitdistrplus::fitdist(df[,1], "nbinom")
  sq_fdist <- squash_fits(fits_fdist)
  expect_true(tibble::is_tibble(sq_fdist))
})


test_that("find fits returns nested data frame", {
  fits <- find_fits(big_data[1:33,], id_col="pep_id", group_col = "group", value_col = "value")
  expect_type(fits$fits, "list")
})


test_that("plot_best_loglik stops on no column", {
  expect_error(plot_best_loglik(df))
})
