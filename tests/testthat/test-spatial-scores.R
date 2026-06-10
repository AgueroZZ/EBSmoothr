test_that("spatial_scores extracts vector, matrix, and flash-like list inputs", {
  set.seed(201)

  loc <- seq(0, 1, length.out = 10)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)
  X <- cbind(spatial = x, flat = mean(x) + rnorm(length(loc), sd = 0.08))

  vec_scores <- spatial_scores(
    x,
    locations = loc,
    spatial = list(family = "matern", link = "identity"),
    reference = list(family = "constant"),
    keep_fits = FALSE
  )
  expect_s3_class(vec_scores, "spatial_scores")
  expect_equal(nrow(vec_scores$scores), 1L)
  expect_equal(vec_scores$scores$status, "ok")
  expect_true(is.finite(vec_scores$scores$spatial_score))

  matrix_scores <- spatial_scores(
    X,
    locations = loc,
    components = "flat",
    spatial = list(family = "matern", link = "identity"),
    reference = list(family = "constant")
  )
  expect_equal(matrix_scores$scores$factor_k, 2L)
  expect_equal(names(matrix_scores$fits$spatial), "flat")
  expect_true(is.finite(matrix_scores$scores$reference_score))

  list_scores <- spatial_scores(
    list(L_pm = X),
    locations = loc,
    spatial = list(family = "matern", link = "identity"),
    reference = NULL,
    keep_fits = FALSE
  )
  expect_equal(nrow(list_scores$scores), 2L)
  expect_true(all(is.na(list_scores$scores$reference_score)))
})

test_that("spatial_scores validates constant rename and supports Matern log link", {
  set.seed(202)

  loc <- seq(0, 1, length.out = 8)
  x <- exp(0.1 + 0.2 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.02)

  expect_error(
    spatial_scores(
      x,
      locations = loc,
      spatial = list(family = "matern", link = "identity"),
      reference = list(family = "nonspatial")
    ),
    "renamed"
  )

  scores <- spatial_scores(
    x,
    locations = loc,
    spatial = list(family = "matern", link = "log"),
    reference = list(family = "point_exponential", profile_s_upper = 2),
    keep_fits = FALSE
  )
  expect_equal(scores$scores$status, "ok")
  expect_equal(scores$scores$spatial_link, "log")
  expect_equal(scores$scores$spatial_backend_requested, "auto")
  expect_equal(scores$scores$spatial_backend, "fisher_pql")
  expect_equal(scores$scores$spatial_laplace_curvature, "fisher")
  expect_equal(scores$scores$matern_backend, "fisher_pql")
  expect_equal(scores$scores$matern_laplace_curvature, "fisher")
  expect_match(scores$scores$spatial_score_semantics, "^fisher_laplace_at_fisher_pql_mode_")
  expect_equal(scores$scores$reference_family, "point_exponential")
  expect_true(is.finite(scores$scores$spatial_score))
  expect_true(is.finite(scores$scores$reference_score))
})

test_that("spatial_scores supports all reference families with profiled point-family noise", {
  set.seed(203)

  loc <- seq(0, 1, length.out = 10)
  x <- exp(0.1 + 0.15 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.03)

  for (family in c("constant", "point_exponential", "point_normal", "point_laplace")) {
    scores <- spatial_scores(
      x,
      locations = loc,
      spatial = list(family = "matern", link = "identity"),
      reference = list(family = family, profile_s_upper = 2),
      keep_fits = FALSE
    )
    expect_equal(scores$scores$status, "ok")
    expect_equal(scores$scores$reference_family, family)
    expect_true(is.finite(scores$scores$reference_score))
    expect_true(is.finite(scores$scores$reference_noise_sd))
  }
})

test_that("spatial_select works in score and permutation modes", {
  set.seed(204)

  loc <- seq(0, 1, length.out = 10)
  X <- cbind(
    smooth = 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08),
    noisy = rnorm(length(loc), sd = 0.3)
  )

  score_selection <- spatial_select(
    X,
    locations = loc,
    spatial = list(family = "matern", link = "identity"),
    reference = list(family = "constant"),
    method = "score"
  )
  expect_s3_class(score_selection, "spatial_select")
  expect_equal(nrow(score_selection$selection_summary), 2L)
  expect_true(all(score_selection$selection_summary$selection %in% c("spatial", "reference", "tie", NA)))
  expect_true(all(c("selection", "selected_family", "route_label") %in% names(score_selection$selection_summary)))

  permutation_selection <- spatial_select(
    X[, 1],
    locations = loc,
    spatial = list(family = "matern", link = "identity"),
    reference = list(family = "constant"),
    method = "permutation",
    B = 3,
    seed = 1,
    refit = FALSE
  )
  expect_s3_class(permutation_selection, "spatial_select")
  expect_true(is.finite(permutation_selection$selection_summary$p_value))
  expect_equal(permutation_selection$selection_summary$selection_method, "permutation")
  expect_true("matern_fits" %in% names(permutation_selection))
})

test_that("spatial_scores_permutation is deterministic and supports refit choices", {
  set.seed(205)

  loc <- seq(0, 1, length.out = 9)
  x <- 0.2 + sin(2 * pi * loc) + rnorm(length(loc), sd = 0.08)

  perm_refit_1 <- spatial_scores_permutation(
    x,
    locations = loc,
    B = 3,
    seed = 99,
    spatial = list(family = "matern", link = "identity"),
    refit = TRUE
  )
  perm_refit_2 <- spatial_scores_permutation(
    x,
    locations = loc,
    B = 3,
    seed = 99,
    spatial = list(family = "matern", link = "identity"),
    refit = TRUE
  )
  expect_equal(perm_refit_1$permutation_scores$score, perm_refit_2$permutation_scores$score)
  expect_true(is.finite(perm_refit_1$summary$p_value))
  expect_equal(perm_refit_1$summary$refit, TRUE)

  perm_fixed <- spatial_scores_permutation(
    x,
    locations = loc,
    B = 3,
    seed = 99,
    spatial = list(family = "matern", link = "identity")
  )
  expect_equal(nrow(perm_fixed$permutation_scores), 3L)
  expect_equal(perm_fixed$summary$refit, FALSE)
  expect_true(is.finite(perm_fixed$summary$p_value))
  expect_equal(
    perm_fixed$summary$p_value,
    mean(perm_fixed$permutation_scores$score >= perm_fixed$summary$observed_score)
  )

  x_log <- exp(0.1 + 0.15 * sin(2 * pi * loc)) + rnorm(length(loc), sd = 0.03)
  perm_log <- spatial_scores_permutation(
    x_log,
    locations = loc,
    B = 3,
    seed = 100,
    spatial = list(family = "matern", link = "log"),
    refit = FALSE
  )
  expect_equal(perm_log$summary$spatial_backend, "fisher_pql")
  expect_equal(perm_log$summary$spatial_laplace_curvature, "fisher")
  expect_match(perm_log$summary$observed_score_semantics, "^fisher_laplace_at_fisher_pql_mode_")
  expect_true(all(perm_log$permutation_scores$laplace_curvature == "fisher"))
  expect_equal(
    perm_log$summary$p_value,
    mean(perm_log$permutation_scores$score >= perm_log$summary$observed_score)
  )
})
