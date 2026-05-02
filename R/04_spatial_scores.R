.spatial_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.spatial_family_name <- function(spec, default_family, allowed, nm) {
  if (is.null(spec)) spec <- list()
  if (!is.list(spec)) {
    stop("`", nm, "` must be a named list of `eb_smoother()` arguments.")
  }
  family <- .spatial_null_coalesce(spec$family, default_family)
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("`", nm, "$family` must be a single character string.")
  }
  if (identical(family, "nonspatial")) {
    stop("`family = \"nonspatial\"` has been renamed to `family = \"constant\"`.")
  }
  if (!(family %in% allowed)) {
    stop(
      "`", nm, "$family` must be one of ",
      paste(sprintf('"%s"', allowed), collapse = ", "),
      "."
    )
  }
  spec$family <- family
  spec
}

.spatial_validate_spatial_spec <- function(spatial) {
  .spatial_family_name(
    spatial,
    default_family = "matern",
    allowed = "matern",
    nm = "spatial"
  )
}

.spatial_validate_reference_spec <- function(reference) {
  if (is.null(reference)) return(NULL)
  .spatial_family_name(
    reference,
    default_family = "constant",
    allowed = c("constant", .point_family_names),
    nm = "reference"
  )
}

.spatial_extract_matrix <- function(object) {
  if (is.numeric(object) && is.null(dim(object))) {
    return(matrix(as.numeric(object), ncol = 1L))
  }

  if (is.matrix(object) || is.data.frame(object)) {
    out <- as.matrix(object)
    storage.mode(out) <- "double"
    return(out)
  }

  if (is.list(object)) {
    candidate_names <- c("L_pm", "L", "loadings", "EL", "F_pm", "F")
    for (nm in candidate_names) {
      if (!is.null(object[[nm]])) {
        out <- as.matrix(object[[nm]])
        storage.mode(out) <- "double"
        return(out)
      }
    }
  }

  if (requireNamespace("flashier", quietly = TRUE)) {
    flash_get_ldf <- getExportedValue("flashier", "flash_get_ldf")
    ldf <- tryCatch(
      flash_get_ldf(object, type = "pm"),
      error = function(e) NULL
    )
    if (!is.null(ldf) && !is.null(ldf$L)) {
      out <- as.matrix(ldf$L)
      storage.mode(out) <- "double"
      return(out)
    }
  }

  stop(
    "`object` must be a numeric vector, a matrix/data frame, a list with ",
    "`L_pm`, or a flashier object with extractable posterior mean loadings."
  )
}

.spatial_component_indices <- function(x, components = NULL) {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Internal error: component data must be a matrix.")
  }
  if (ncol(x) < 1L) stop("`object` must contain at least one component.")

  if (is.null(components)) {
    idx <- seq_len(ncol(x))
  } else if (is.character(components)) {
    if (is.null(colnames(x))) {
      stop("Character `components` require column names on the extracted matrix.")
    }
    idx <- match(components, colnames(x))
    if (anyNA(idx)) {
      stop("Unknown component name(s): ", paste(components[is.na(idx)], collapse = ", "))
    }
  } else {
    idx <- as.integer(components)
    if (length(idx) != length(components) || anyNA(idx) ||
        any(idx < 1L) || any(idx > ncol(x))) {
      stop("Numeric `components` must be valid column indices.")
    }
  }

  list(
    index = idx,
    label = if (is.null(colnames(x))) as.character(idx) else colnames(x)[idx]
  )
}

.spatial_merge_args <- function(base, defaults) {
  if (!length(defaults)) return(base)
  base_names <- names(base)
  default_names <- names(defaults)
  if (is.null(base_names)) base_names <- rep("", length(base))
  if (is.null(default_names)) default_names <- rep("", length(defaults))
  add <- default_names != "" & !default_names %in% base_names
  c(base, defaults[add])
}

.spatial_fit_component <- function(x, spec, locations = NULL, common_args = list()) {
  args <- .spatial_merge_args(
    base = spec,
    defaults = common_args
  )
  args$x <- as.numeric(x)
  if (identical(args$family, "matern") &&
      is.null(args$locations) && is.null(args$setup)) {
    args$locations <- locations
  }
  do.call(eb_smoother, args)
}

.spatial_fit_score <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$log_likelihood)) return(NA_real_)
  as.numeric(fit$log_likelihood)
}

.spatial_fit_semantics <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$log_likelihood_semantics)) {
    return(NA_character_)
  }
  as.character(fit$log_likelihood_semantics)
}

.spatial_fit_link <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$link)) return(NA_character_)
  as.character(fit$link)
}

.spatial_fit_backend <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$backend)) return(NA_character_)
  as.character(fit$backend)
}

.spatial_fit_laplace_curvature <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$laplace_curvature)) return(NA_character_)
  as.character(fit$laplace_curvature)
}

.spatial_fit_beta <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$fitted_beta) || !length(fit$fitted_beta)) {
    return(NA_real_)
  }
  as.numeric(fit$fitted_beta)[1]
}

.spatial_fit_noise <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$fitted_noise_sd)) return(NA_real_)
  as.numeric(fit$fitted_noise_sd)[1]
}

.spatial_fit_matern_range <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$fitted_g$theta)) return(NA_real_)
  as.numeric(exp(fit$fitted_g$theta))[1]
}

.spatial_fit_matern_sigma <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$fitted_g$sigma)) return(NA_real_)
  as.numeric(fit$fitted_g$sigma)[1]
}

.spatial_fit_laplace_implementation <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$raw_fit$laplace_implementation)) {
    return(NA_character_)
  }
  as.character(fit$raw_fit$laplace_implementation)
}

.spatial_fit_point_mass <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$point_mass_probability)) {
    return(NA_real_)
  }
  as.numeric(fit$point_mass_probability)[1]
}

.spatial_fit_point_scale <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$point_nonzero_scale)) {
    return(NA_real_)
  }
  as.numeric(fit$point_nonzero_scale)[1]
}

.spatial_row_status <- function(spatial_fit, reference_fit = NULL) {
  errors <- character()
  if (inherits(spatial_fit, "error")) {
    errors <- c(errors, paste0("spatial_error: ", conditionMessage(spatial_fit)))
  }
  if (!is.null(reference_fit) && inherits(reference_fit, "error")) {
    errors <- c(errors, paste0("reference_error: ", conditionMessage(reference_fit)))
  }
  list(
    status = if (length(errors)) "error" else "ok",
    error_message = if (length(errors)) paste(errors, collapse = " | ") else NA_character_
  )
}

.spatial_score_row <- function(factor_k,
                               component,
                               spatial_fit,
                               reference_fit,
                               spatial_spec,
                               reference_spec) {
  status <- .spatial_row_status(spatial_fit, reference_fit)
  spatial_score <- .spatial_fit_score(spatial_fit)
  reference_score <- .spatial_fit_score(reference_fit)
  delta <- spatial_score - reference_score
  if (is.null(reference_fit)) {
    reference_score <- NA_real_
    delta <- NA_real_
  }

  spatial_semantics <- .spatial_fit_semantics(spatial_fit)
  reference_semantics <- .spatial_fit_semantics(reference_fit)
  score_semantics <- if (!is.na(spatial_semantics) && !is.na(reference_semantics)) {
    paste(spatial_semantics, reference_semantics, sep = "_vs_")
  } else {
    NA_character_
  }

  reference_family <- if (is.null(reference_spec)) NA_character_ else reference_spec$family
  point_family <- if (!is.na(reference_family) && reference_family %in% .point_family_names) {
    reference_family
  } else {
    NA_character_
  }

  data.frame(
    factor_k = factor_k,
    component = component,
    status = status$status,
    error_message = status$error_message,
    spatial_family = spatial_spec$family,
    spatial_link = .spatial_fit_link(spatial_fit),
    spatial_score = spatial_score,
    spatial_score_semantics = spatial_semantics,
    spatial_backend_requested = .spatial_null_coalesce(spatial_spec$backend, "auto"),
    spatial_backend = .spatial_fit_backend(spatial_fit),
    spatial_laplace_curvature = .spatial_fit_laplace_curvature(spatial_fit),
    spatial_laplace_implementation = .spatial_fit_laplace_implementation(spatial_fit),
    reference_family = reference_family,
    reference_link = .spatial_fit_link(reference_fit),
    reference_score = reference_score,
    reference_score_semantics = reference_semantics,
    reference_backend = .spatial_fit_backend(reference_fit),
    score_semantics = score_semantics,
    delta = delta,
    spatial_beta = .spatial_fit_beta(spatial_fit),
    reference_beta = .spatial_fit_beta(reference_fit),
    spatial_noise_sd = .spatial_fit_noise(spatial_fit),
    reference_noise_sd = .spatial_fit_noise(reference_fit),
    matern_range = if (identical(spatial_spec$family, "matern")) .spatial_fit_matern_range(spatial_fit) else NA_real_,
    matern_sigma = if (identical(spatial_spec$family, "matern")) .spatial_fit_matern_sigma(spatial_fit) else NA_real_,
    point_reference_family = point_family,
    point_reference_noise_sd = if (!is.na(point_family)) .spatial_fit_noise(reference_fit) else NA_real_,
    point_reference_mass = if (!is.na(point_family)) .spatial_fit_point_mass(reference_fit) else NA_real_,
    point_reference_scale = if (!is.na(point_family)) .spatial_fit_point_scale(reference_fit) else NA_real_,
    loglik_matern = if (identical(spatial_spec$family, "matern")) spatial_score else NA_real_,
    loglik_nonspatial = reference_score,
    loglik_point_exponential = if (identical(reference_family, "point_exponential")) reference_score else NA_real_,
    matern_backend_requested = .spatial_null_coalesce(spatial_spec$backend, "auto"),
    matern_backend = if (identical(spatial_spec$family, "matern")) .spatial_fit_backend(spatial_fit) else NA_character_,
    matern_laplace_curvature = if (identical(spatial_spec$family, "matern")) .spatial_fit_laplace_curvature(spatial_fit) else NA_character_,
    matern_laplace_implementation = if (identical(spatial_spec$family, "matern")) .spatial_fit_laplace_implementation(spatial_fit) else NA_character_,
    matern_beta = if (identical(spatial_spec$family, "matern")) .spatial_fit_beta(spatial_fit) else NA_real_,
    nonspatial_beta = if (identical(reference_family, "constant")) .spatial_fit_beta(reference_fit) else NA_real_,
    matern_noise_sd = if (identical(spatial_spec$family, "matern")) .spatial_fit_noise(spatial_fit) else NA_real_,
    nonspatial_noise_sd = .spatial_fit_noise(reference_fit),
    point_exponential_noise_sd = if (identical(reference_family, "point_exponential")) .spatial_fit_noise(reference_fit) else NA_real_,
    point_exponential_pi0 = if (identical(reference_family, "point_exponential")) .spatial_fit_point_mass(reference_fit) else NA_real_,
    point_exponential_scale = if (identical(reference_family, "point_exponential")) .spatial_fit_point_scale(reference_fit) else NA_real_,
    stringsAsFactors = FALSE
  )
}

.spatial_apply_score_selection <- function(scores, tol, reference_family = "constant") {
  out <- scores
  delta <- as.numeric(out$delta)
  selection <- ifelse(
    out$status != "ok" | !is.finite(delta),
    NA_character_,
    ifelse(delta > tol, "spatial", ifelse(delta < -tol, "nonspatial", "tie"))
  )
  route_label <- ifelse(selection == "spatial", "spatial", "nonspatial")
  route_label[is.na(selection)] <- "nonspatial"
  selected_family <- ifelse(selection == "spatial", out$spatial_family, out$reference_family)
  tie_idx <- !is.na(selection) & selection == "tie"
  selected_family[tie_idx] <- "tie"
  selected_family[is.na(selection)] <- reference_family

  out$selection <- selection
  out$route_label <- route_label
  out$selected_family <- selected_family
  out$selection_method <- "score"
  out
}

.spatial_fit_matern_g_init <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$fitted_g)) return(NULL)
  fit$fitted_g
}

#' Score components under spatial and reference smoothers
#'
#' @description
#' Fits each selected component under a spatial smoother and, optionally, a
#' reference smoother, then returns the fitted scores in a component-level table.
#'
#' @details
#' `spatial_scores()` is a public wrapper around [eb_smoother()] for the common
#' workflow of scoring factor/loadings columns from a fitted factorization. The
#' input `object` may be a numeric vector, a matrix/data frame whose columns are
#' components, a flash-like list with `L_pm`, or a flashier object whose
#' posterior mean loadings can be extracted.
#'
#' The spatial score is `as.numeric(fit$log_likelihood)` from the spatial
#' [eb_smoother()] fit. For Matern fits this is the marginal log likelihood, or
#' the PC-penalized marginal log likelihood when `spatial$pc.penalty` is
#' supplied. The score semantics for each fit are reported in
#' `spatial_score_semantics` and `reference_score_semantics`.
#'
#' The `spatial` list is passed to [eb_smoother()] after adding `family =
#' "matern"` by default. Use `spatial = list(family = "matern", link = "log")`
#' for positive smooth signals on the log link; with `backend = "auto"` this
#' uses Fisher Laplace scoring and reports `laplace_fisher_<beta_mode>` score
#' semantics. The `reference` list, when supplied, supports
#' `family = "constant"`, `"point_exponential"`,
#' `"point_normal"`, and `"point_laplace"`. The old name
#' `family = "nonspatial"` is intentionally rejected; use `family = "constant"`
#' for the constant Gaussian reference.
#'
#' @param object A numeric vector, matrix/data frame, flash-like list, or
#'   flashier object containing component loadings. Matrix columns are treated
#'   as components.
#' @param locations Spatial coordinates used for Matern fits unless
#'   `spatial$locations` or `spatial$setup` is supplied.
#' @param components Optional component indices or names to score.
#' @param spatial Named list of spatial [eb_smoother()] arguments. Defaults to
#'   `list(family = "matern", link = "identity")`.
#' @param reference Optional named list of reference [eb_smoother()] arguments.
#'   If supplied without `family`, `family = "constant"` is used.
#' @param keep_fits If `TRUE`, return fitted `eb_smoother_fit` objects.
#' @param ... Additional [eb_smoother()] arguments used as defaults for both
#'   `spatial` and `reference`; model-specific options should usually be placed
#'   directly inside the corresponding list.
#'
#' @return A list of class `"spatial_scores"` with fields:
#' \describe{
#'   \item{scores}{A data frame with one row per component. Important columns
#'   include `spatial_score`, `reference_score`, `delta`,
#'   `spatial_score_semantics`, `reference_score_semantics`,
#'   `spatial_laplace_curvature`, `matern_range`, `matern_sigma`, and
#'   point-reference summaries when applicable.}
#'   \item{fits}{When `keep_fits = TRUE`, a list with `spatial` and
#'   `reference` fit lists.}
#'   \item{screening_summary}{Alias of `scores` for simulation workflows.}
#'   \item{matern_fits,nonspatial_fits}{Compatibility aliases for spatial and
#'   reference fits when `keep_fits = TRUE`.}
#' }
#'
#' @examples
#' loc <- seq(0, 1, length.out = 12)
#' X <- cbind(sin(2 * pi * loc), rnorm(length(loc)))
#' scores <- spatial_scores(
#'   X,
#'   locations = loc,
#'   spatial = list(family = "matern", link = "identity"),
#'   reference = list(family = "constant")
#' )
#' scores$scores[, c("factor_k", "spatial_score", "reference_score", "delta")]
#'
#' positive_scores <- spatial_scores(
#'   exp(X[, 1]),
#'   locations = loc,
#'   spatial = list(family = "matern", link = "log"),
#'   reference = list(family = "point_exponential")
#' )
#' positive_scores$scores$spatial_score_semantics
#'
#' @export
spatial_scores <- function(object,
                           locations = NULL,
                           components = NULL,
                           spatial = list(family = "matern", link = "identity"),
                           reference = NULL,
                           keep_fits = TRUE,
                           ...) {
  X <- .spatial_extract_matrix(object)
  component_info <- .spatial_component_indices(X, components)
  spatial_spec <- .spatial_validate_spatial_spec(spatial)
  reference_spec <- .spatial_validate_reference_spec(reference)
  common_args <- list(...)

  spatial_fits <- vector("list", length(component_info$index))
  reference_fits <- if (is.null(reference_spec)) NULL else vector("list", length(component_info$index))

  rows <- lapply(seq_along(component_info$index), function(i) {
    k <- component_info$index[i]
    x <- X[, k]
    spatial_fits[[i]] <<- tryCatch(
      .spatial_fit_component(
        x = x,
        spec = spatial_spec,
        locations = locations,
        common_args = common_args
      ),
      error = function(e) e
    )
    reference_fit <- NULL
    if (!is.null(reference_spec)) {
      reference_fits[[i]] <<- tryCatch(
        .spatial_fit_component(
          x = x,
          spec = reference_spec,
          locations = NULL,
          common_args = common_args
        ),
        error = function(e) e
      )
      reference_fit <- reference_fits[[i]]
    }
    .spatial_score_row(
      factor_k = k,
      component = component_info$label[i],
      spatial_fit = spatial_fits[[i]],
      reference_fit = reference_fit,
      spatial_spec = spatial_spec,
      reference_spec = reference_spec
    )
  })

  score_df <- do.call(rbind, rows)
  out <- list(
    scores = score_df,
    screening_summary = score_df,
    components = component_info,
    spatial = spatial_spec,
    reference = reference_spec
  )
  if (isTRUE(keep_fits)) {
    names(spatial_fits) <- component_info$label
    if (!is.null(reference_fits)) names(reference_fits) <- component_info$label
    out$fits <- list(spatial = spatial_fits, reference = reference_fits)
    out$matern_fits <- spatial_fits
    out$nonspatial_fits <- reference_fits
    if (!is.null(reference_spec) && identical(reference_spec$family, "point_exponential")) {
      out$point_exponential_fits <- reference_fits
    }
  }
  class(out) <- c("spatial_scores", "list")
  out
}

#' Select spatial versus reference components
#'
#' @description
#' Classifies each component as spatial or reference using either a direct score
#' comparison or a permutation p-value for the spatial score.
#'
#' @details
#' With `method = "score"`, the function calls [spatial_scores()] and selects
#' the spatial model when `spatial_score - reference_score > tol`. It selects
#' the nonspatial/reference route when the difference is below `-tol`;
#' otherwise the row is marked as a tie.
#'
#' With `method = "permutation"`, the function calls
#' [spatial_scores_permutation()] and selects the spatial model when the
#' one-sided permutation p-value is at most `alpha`. The `reference` family is
#' still recorded as the fallback selected family for components that do not
#' pass the permutation screen.
#'
#' The `spatial` list may request a Matern log-link fit via
#' `spatial = list(family = "matern", link = "log")`; all other Matern options,
#' including PC-prior options such as `pc.penalty`, are passed through to
#' [eb_smoother()].
#'
#' @param object A numeric vector, matrix/data frame, flash-like list, or
#'   flashier object containing component loadings.
#' @param locations Spatial coordinates for Matern fits.
#' @param spatial Named list of spatial [eb_smoother()] arguments.
#' @param reference Named list of reference [eb_smoother()] arguments. Defaults
#'   to `list(family = "constant")`.
#' @param method Selection method: `"score"` or `"permutation"`.
#' @param tol Score-difference tolerance for `method = "score"`.
#' @param alpha One-sided p-value cutoff for `method = "permutation"`.
#' @param ... Additional arguments passed to [spatial_scores()] for
#'   `method = "score"` or to [spatial_scores_permutation()] for
#'   `method = "permutation"`; examples include `components`, `B`, `seed`,
#'   `refit`, and `ncores`.
#'
#' @return A list of class `"spatial_select"` with a `selection_summary` data
#'   frame. Score mode also returns `scores` and, when kept, fit lists.
#'   Permutation mode returns the `permutation` result from
#'   [spatial_scores_permutation()].
#'
#' @examples
#' loc <- seq(0, 1, length.out = 12)
#' X <- cbind(sin(2 * pi * loc), rnorm(length(loc)))
#' selection <- spatial_select(
#'   X,
#'   locations = loc,
#'   spatial = list(family = "matern", link = "identity"),
#'   reference = list(family = "constant"),
#'   method = "score"
#' )
#' selection$selection_summary[, c("factor_k", "selection", "selected_family")]
#'
#' @export
spatial_select <- function(object,
                           locations,
                           spatial = list(family = "matern", link = "identity"),
                           reference = list(family = "constant"),
                           method = c("score", "permutation"),
                           tol = 1e-6,
                           alpha = 0.05,
                           ...) {
  method <- match.arg(method)
  reference_spec <- .spatial_validate_reference_spec(reference)
  if (is.null(reference_spec)) {
    reference_spec <- list(family = "constant")
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("`tol` must be a single non-negative numeric.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("`alpha` must be between 0 and 1.")
  }

  if (identical(method, "score")) {
    scores <- spatial_scores(
      object = object,
      locations = locations,
      spatial = spatial,
      reference = reference_spec,
      keep_fits = TRUE,
      ...
    )
    selection_summary <- .spatial_apply_score_selection(
      scores = scores$scores,
      tol = tol,
      reference_family = reference_spec$family
    )
    scores$screening_summary <- selection_summary
    out <- c(
      list(
        selection_summary = selection_summary,
        screening_summary = selection_summary,
        method = "score",
        tol = tol,
        scores = scores
      ),
      scores[intersect(names(scores), c("fits", "matern_fits", "nonspatial_fits", "point_exponential_fits"))]
    )
    class(out) <- c("spatial_select", "list")
    return(out)
  }

  permutation <- spatial_scores_permutation(
    object = object,
    locations = locations,
    spatial = spatial,
    ...
  )
  summary <- permutation$summary
  selection <- ifelse(
    summary$status != "ok" | !is.finite(summary$p_value),
    NA_character_,
    ifelse(summary$p_value <= alpha, "spatial", "nonspatial")
  )
  summary$selection <- selection
  summary$route_label <- ifelse(!is.na(selection) & selection == "spatial", "spatial", "nonspatial")
  summary$route_label[is.na(selection)] <- "nonspatial"
  summary$reference_family <- reference_spec$family
  summary$selected_family <- ifelse(!is.na(selection) & selection == "spatial", summary$spatial_family, reference_spec$family)
  summary$selection_method <- "permutation"
  summary$alpha <- alpha

  out <- list(
    selection_summary = summary,
    screening_summary = summary,
    method = "permutation",
    alpha = alpha,
    permutation = permutation
  )
  if (!is.null(permutation$matern_fits)) {
    out$matern_fits <- permutation$matern_fits
  }
  if (!is.null(permutation$fits)) {
    out$fits <- permutation$fits
  }
  class(out) <- c("spatial_select", "list")
  out
}

#' Compute permutation scores for a spatial smoother
#'
#' @description
#' Computes observed spatial scores and permutation-null spatial scores by
#' randomly breaking the association between component values and locations.
#'
#' @details
#' `spatial_scores_permutation()` first fits each selected component under the
#' requested spatial model and records the observed spatial score. It then
#' permutes the component vector `B` times relative to the fixed locations and
#' re-scores with fixed Matern parameters by default. Set `refit = TRUE` to
#' refit the spatial parameters for every permuted data set.
#'
#' The current spatial model is Matern. Use `spatial = list(family = "matern",
#' link = "log")` for positive smooth signals. With the default
#' `refit = FALSE`, the
#' observed fit's Matern range and marginal standard deviation are supplied as
#' `g_init` with `fix_g = TRUE`; if the user also supplies `spatial$fix_params`,
#' those values are merged with the permutation fixed-parameter behavior. The
#' observation noise and intercept are still fitted according to the selected
#' [eb_smoother()] semantics.
#'
#' The p-value is the one-sided upper-tail permutation p-value
#' `number of permutation scores >= observed score / number of valid
#' permutation scores`. The observed score is not added to the permutation-null
#' scores. Larger scores indicate better support for the spatial model.
#'
#' @param object A numeric vector, matrix/data frame, flash-like list, or
#'   flashier object containing component loadings.
#' @param locations Spatial coordinates for Matern fits.
#' @param B Number of random permutations.
#' @param spatial Named list of spatial [eb_smoother()] arguments. Defaults to
#'   `list(family = "matern", link = "identity")`.
#' @param components Optional component indices or names to score.
#' @param refit If `FALSE`, reuse the observed fit's Matern parameters with
#'   `fix_g = TRUE`. If `TRUE`, refit Matern parameters for every permuted data
#'   set.
#' @param ncores Number of cores for permutation fits. On Unix-like platforms,
#'   values greater than one use `parallel::mclapply`; otherwise serial
#'   evaluation is used with a warning.
#' @param seed Optional random seed used to generate the permutations.
#' @param keep_fits If `TRUE`, return observed and permutation fit objects.
#' @param ... Additional [eb_smoother()] arguments used as defaults unless
#'   overridden by `spatial`.
#'
#' @return A list of class `"spatial_scores_permutation"` with fields:
#' \describe{
#'   \item{summary}{One row per component with observed score, p-value, number
#'   of valid permutations, fitted Matern parameters, and score semantics.}
#'   \item{permutation_scores}{One row per component and permutation.}
#'   \item{observed_scores}{Observed-score table before p-value calculation.}
#'   \item{matern_fits}{Observed spatial fits, useful for downstream smoothing
#'   of selected components.}
#'   \item{fits}{When `keep_fits = TRUE`, observed and permutation fit lists.}
#' }
#'
#' @examples
#' loc <- seq(0, 1, length.out = 10)
#' x <- sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)
#' perm <- spatial_scores_permutation(
#'   x,
#'   locations = loc,
#'   B = 3,
#'   seed = 1,
#'   spatial = list(family = "matern", link = "identity")
#' )
#' perm$summary[, c("factor_k", "observed_score", "p_value")]
#'
#' @export
spatial_scores_permutation <- function(object,
                                       locations,
                                       B = 100,
                                       spatial = list(family = "matern", link = "identity"),
                                       components = NULL,
                                       refit = FALSE,
                                       ncores = 1,
                                       seed = NULL,
                                       keep_fits = FALSE,
                                       ...) {
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1 || B != floor(B)) {
    stop("`B` must be a positive integer.")
  }
  B <- as.integer(B)
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit)) {
    stop("`refit` must be TRUE or FALSE.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1L || is.na(ncores) ||
      ncores < 1 || ncores != floor(ncores)) {
    stop("`ncores` must be a positive integer.")
  }
  ncores <- as.integer(ncores)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be NULL or a single numeric seed.")
    }
    set.seed(seed)
  }

  X <- .spatial_extract_matrix(object)
  component_info <- .spatial_component_indices(X, components)
  spatial_spec <- .spatial_validate_spatial_spec(spatial)
  common_args <- list(...)
  permutations <- replicate(B, sample.int(nrow(X)), simplify = FALSE)

  observed_fits <- vector("list", length(component_info$index))
  permutation_fits <- if (isTRUE(keep_fits)) vector("list", length(component_info$index)) else NULL
  permutation_rows <- list()

  observed_rows <- lapply(seq_along(component_info$index), function(i) {
    k <- component_info$index[i]
    x <- X[, k]
    observed_fits[[i]] <<- tryCatch(
      .spatial_fit_component(
        x = x,
        spec = spatial_spec,
        locations = locations,
        common_args = common_args
      ),
      error = function(e) e
    )
    status <- .spatial_row_status(observed_fits[[i]])
    data.frame(
      factor_k = k,
      component = component_info$label[i],
      status = status$status,
      error_message = status$error_message,
      spatial_family = spatial_spec$family,
      observed_score = .spatial_fit_score(observed_fits[[i]]),
      observed_score_semantics = .spatial_fit_semantics(observed_fits[[i]]),
      spatial_link = .spatial_fit_link(observed_fits[[i]]),
      spatial_backend = .spatial_fit_backend(observed_fits[[i]]),
      spatial_laplace_curvature = .spatial_fit_laplace_curvature(observed_fits[[i]]),
      spatial_laplace_implementation = .spatial_fit_laplace_implementation(observed_fits[[i]]),
      matern_range = .spatial_fit_matern_range(observed_fits[[i]]),
      matern_sigma = .spatial_fit_matern_sigma(observed_fits[[i]]),
      matern_beta = .spatial_fit_beta(observed_fits[[i]]),
      matern_noise_sd = .spatial_fit_noise(observed_fits[[i]]),
      stringsAsFactors = FALSE
    )
  })
  observed_df <- do.call(rbind, observed_rows)

  for (i in seq_along(component_info$index)) {
    k <- component_info$index[i]
    x <- X[, k]
    observed_fit <- observed_fits[[i]]
    perm_spec <- spatial_spec
    if (!isTRUE(refit)) {
      if (inherits(observed_fit, "error")) {
        perm_spec$g_init <- NULL
      } else {
        perm_spec$g_init <- .spatial_fit_matern_g_init(observed_fit)
        existing_fix_params <- perm_spec$fix_params
        if (is.null(existing_fix_params)) existing_fix_params <- character()
        perm_spec$fix_params <- unique(c(existing_fix_params, "range", "sigma"))
        perm_spec$fix_g <- TRUE
      }
    }

    fit_one <- function(b) {
      if (inherits(observed_fit, "error")) {
        fit <- observed_fit
      } else {
        fit <- tryCatch(
          .spatial_fit_component(
            x = x[permutations[[b]]],
            spec = perm_spec,
            locations = locations,
            common_args = common_args
          ),
          error = function(e) e
        )
      }
      status <- .spatial_row_status(fit)
      list(
        row = data.frame(
          factor_k = k,
          component = component_info$label[i],
          permutation = b,
          score = .spatial_fit_score(fit),
          score_semantics = .spatial_fit_semantics(fit),
          laplace_curvature = .spatial_fit_laplace_curvature(fit),
          status = status$status,
          error_message = status$error_message,
          stringsAsFactors = FALSE
        ),
        fit = fit
      )
    }

    perm_results <- if (ncores > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(seq_len(B), fit_one, mc.cores = ncores)
    } else {
      if (ncores > 1L) {
        warning("`ncores > 1` is only used on Unix-like platforms; using serial permutation fits.")
      }
      lapply(seq_len(B), fit_one)
    }
    permutation_rows[[i]] <- do.call(rbind, lapply(perm_results, `[[`, "row"))
    if (isTRUE(keep_fits)) {
      permutation_fits[[i]] <- lapply(perm_results, `[[`, "fit")
    }
  }

  permutation_df <- do.call(rbind, permutation_rows)
  summary_df <- observed_df
  summary_df$B <- B
  summary_df$refit <- refit
  summary_df$n_permutations_ok <- vapply(seq_len(nrow(summary_df)), function(i) {
    rows <- permutation_df$factor_k == summary_df$factor_k[i] & permutation_df$status == "ok"
    sum(rows)
  }, integer(1))
  summary_df$permutation_score_mean <- vapply(seq_len(nrow(summary_df)), function(i) {
    scores <- permutation_df$score[permutation_df$factor_k == summary_df$factor_k[i] & permutation_df$status == "ok"]
    if (!length(scores)) NA_real_ else mean(scores)
  }, numeric(1))
  summary_df$permutation_score_sd <- vapply(seq_len(nrow(summary_df)), function(i) {
    scores <- permutation_df$score[permutation_df$factor_k == summary_df$factor_k[i] & permutation_df$status == "ok"]
    if (length(scores) < 2L) NA_real_ else stats::sd(scores)
  }, numeric(1))
  summary_df$p_value <- vapply(seq_len(nrow(summary_df)), function(i) {
    obs <- summary_df$observed_score[i]
    scores <- permutation_df$score[permutation_df$factor_k == summary_df$factor_k[i] & permutation_df$status == "ok"]
    if (!is.finite(obs) || !length(scores)) return(NA_real_)
    sum(scores >= obs) / length(scores)
  }, numeric(1))

  out <- list(
    summary = summary_df,
    observed_scores = observed_df,
    permutation_scores = permutation_df,
    components = component_info,
    spatial = spatial_spec,
    B = B,
    refit = refit,
    seed = seed,
    matern_fits = stats::setNames(observed_fits, component_info$label)
  )
  if (isTRUE(keep_fits)) {
    names(observed_fits) <- component_info$label
    names(permutation_fits) <- component_info$label
    out$fits <- list(observed = observed_fits, permutation = permutation_fits)
  }
  class(out) <- c("spatial_scores_permutation", "list")
  out
}
