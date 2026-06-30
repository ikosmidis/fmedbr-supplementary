base_dir <- "."
results_file <- file.path(base_dir, "results/osm.rda")

pkgs <- c("ordinal", "brglm2", "focuson", "future.apply", "progressr", "tinytest")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 1)
plan(multisession, workers = n_cores)

acl_quantities <- function(theta, X, po = TRUE) {
    X <- as.matrix(X)
    p <- ncol(X)
    ncat <- infer_ncat_acl(theta, X, po = po)
    q <- ncat - 1L
    Tmat <- upper.tri(matrix(1, q, q), diag = TRUE) * 1
    if (po) {
        nslopes <- p - 1L
        if (length(theta) != q + nslopes) {
            stop("for `po = TRUE`, `length(theta)` must be `ncat - 1 + ncol(X) - 1`")
        }
        intercept_coefs <- drop(Tmat %*% theta[seq_len(q)])
        slope_coefs <- if (nslopes > 0L) {
            outer(q:1, theta[-seq_len(q)])
        } else {
            matrix(numeric(0), nrow = q, ncol = 0L)
        }
        coefs <- cbind(intercept_coefs, slope_coefs)
        map <- matrix(0, nrow = q * p, ncol = q + nslopes)
        map[seq_len(q), seq_len(q)] <- Tmat
        if (nslopes > 0L) {
            for (j in seq_len(nslopes)) {
                rows <- q * j + seq_len(q)
                map[rows, q + j] <- q:1
            }
        }
    } else {
        if (length(theta) != q * p) {
            stop("for `po = FALSE`, `length(theta)` must be `(ncat - 1) * ncol(X)`")
        }
        theta_mat <- matrix(theta, nrow = q)
        coefs <- Tmat %*% theta_mat
        map <- kronecker(diag(p), Tmat)
    }
    list(coefs = coefs, map = map, ncat = ncat)
}

acl_probs <- function(theta, X, ncat, po = TRUE) {
    X <- as.matrix(X)
    quantities <- acl_quantities(theta, X, po = po)
    if (!missing(ncat) && !identical(as.integer(ncat), as.integer(quantities$ncat))) {
        stop("`ncat` is incompatible with `theta`, `X`, and `po`")
    }
    ncat <- quantities$ncat
    coefs <- quantities$coefs
    eta <- X %*% t(coefs)
    fits <- cbind(eta, 0)
    exp_fits <- exp(fits - apply(fits, 1, max))
    exp_fits / rowSums(exp_fits)
}

infer_ncat_acl <- function(theta, X, po = TRUE) {
    p <- ncol(as.matrix(X))
    if (po) {
        ncat <- length(theta) - p + 2L
        if (ncat < 2L) {
            stop("`theta` is incompatible with `X` for `po = TRUE`")
        }
    } else {
        ncat <- length(theta) / p + 1L
        if (!isTRUE(all.equal(ncat, round(ncat)))) {
            stop("`theta` is incompatible with `X` for `po = FALSE`")
        }
        ncat <- as.integer(round(ncat))
    }
    ncat
}

make_X <- function(formula, data, contrasts = NULL, na.action = na.omit) {
    Terms <- delete.response(terms(formula))
    raw_data <- get_all_vars(Terms, data)
    mf <- model.frame(Terms, raw_data, na.action = na.action)
    model.matrix(Terms, mf, contrasts.arg = contrasts)
}

make_XY <- function(formula, data, contrasts = NULL, na.action = na.omit) {
    Terms <- terms(formula)
    raw_data <- get_all_vars(Terms, data)
    mf <- model.frame(Terms, raw_data, na.action = na.action)
    X <- model.matrix(delete.response(Terms), mf, contrasts.arg = contrasts)
    Y <- model.response(mf)
    list(Y = Y, X = X)
}

loglik <- function(theta, data, formula, po = TRUE, contrasts = NULL, ...) {
    XY <- make_XY(formula = formula, data = data, contrasts = contrasts)
    Y <- XY$Y
    X <- XY$X
    X <- as.matrix(X)
    if (is.factor(Y)) {
        y <- as.integer(Y)
        ncat <- nlevels(Y)
    } else {
        ylev <- sort(unique(Y))
        y <- match(Y, ylev)
        ncat <- length(ylev)
    }
    probs <- acl_probs(theta, X, ncat = ncat, po = po)
    sum(log(pmax(probs[cbind(seq_along(y), y)], .Machine$double.xmin)))
}

score <- function(theta, data, formula, po = TRUE,
                              contrasts = NULL, ...) {
    XY <- make_XY(formula = formula, data = data, contrasts = contrasts)
    Y <- XY$Y
    X <- as.matrix(XY$X)
    quantities <- acl_quantities(theta, X, po = po)
    ncat <- quantities$ncat
    probs <- acl_probs(theta, X, ncat = ncat, po = po)
    if (is.factor(Y)) {
        y <- as.integer(Y)
    } else {
        ylev <- sort(unique(Y))
        y <- match(Y, ylev)
    }
    Yind <- matrix(0, nrow = length(y), ncol = ncat)
    Yind[cbind(seq_along(y), y)] <- 1
    resid <- Yind[, seq_len(ncat - 1L), drop = FALSE] -
        probs[, seq_len(ncat - 1L), drop = FALSE]
    score_beta <- c(t(resid) %*% X)
    out <- drop(crossprod(quantities$map, score_beta))
    if (!is.null(names(theta))) {
        names(out) <- names(theta)
    }
    out
}


info <- function(theta, data, formula, po = TRUE,
                              contrasts = NULL, ...) {
    XY <- make_XY(formula = formula, data = data, contrasts = contrasts)
    X <- as.matrix(XY$X)
    quantities <- acl_quantities(theta, X, po = po)
    ncat <- quantities$ncat
    probs <- acl_probs(theta, X, ncat = ncat, po = po)
    q <- ncat - 1L
    out <- matrix(0, nrow = length(theta), ncol = length(theta))
    for (i in seq_len(nrow(X))) {
        pp <- probs[i, seq_len(q)]
        S <- diag(pp) - tcrossprod(pp)
        I_beta_i <- kronecker(crossprod(X[i, , drop = FALSE]), S)
        out <- out + crossprod(quantities$map, I_beta_i %*% quantities$map)
    }
    if (!is.null(names(theta))) {
        dimnames(out) <- list(names(theta), names(theta))
    }
    out
}

simulate_acl <- function(theta, simu_data, formula, po = TRUE,
                           ylev = NULL, contrasts = NULL, ...) {
    Terms <- terms(formula)
    raw_simu_data <- get_all_vars(Terms, simu_data)
    mf <- model.frame(Terms, raw_simu_data, na.action = na.omit)
    X <- model.matrix(delete.response(Terms), mf, contrasts.arg = contrasts)
    Y <- model.response(mf)
    X <- as.matrix(X)
    ncat <- infer_ncat_acl(theta, X, po = po)
    probs <- acl_probs(theta, X, ncat = ncat, po = po)
    if (is.null(ylev)) {
        ylev <- if (is.factor(Y)) levels(Y) else seq_len(ncat)
    }
    if (length(ylev) != ncat) {
        stop("`ylev` must have length equal to the number of response categories")
    }
    y <- apply(probs, 1, function(p) sample.int(ncat, 1L, prob = p))
    y <- if (is.factor(Y)) {
        if (is.ordered(Y)) {
            ordered(ylev[y], levels = ylev)
        } else {
            factor(ylev[y], levels = ylev)
        }
    } else {
        y
    }
    out <- as.data.frame(simu_data)[row.names(mf), , drop = FALSE]
    out[[as.character(formula[[2L]])]] <- y
    out
}

make_X0_X1 <- function(formula, data, group, contrasts = NULL) {
    if (inherits(group, "formula")) {
        group <- all.vars(group)
    }
    if (length(group) != 1L) {
        stop("`group` must identify exactly one binary variable")
    }
    Terms <- delete.response(terms(formula))
    raw_data <- get_all_vars(Terms, data)
    mf <- model.frame(Terms, raw_data, na.action = na.omit)
    xlevels <- .getXlevels(Terms, mf)
    X_template <- model.matrix(Terms, mf, contrasts.arg = contrasts)
    contrasts <- attr(X_template, "contrasts")
    if (!(group %in% names(raw_data))) {
        stop("`group` must be a variable in `data`")
    }
    raw_data <- raw_data[row.names(mf), , drop = FALSE]
    group_values <- if (is.factor(raw_data[[group]])) {
        levels(droplevels(raw_data[[group]]))
    } else {
        sort(unique(raw_data[[group]]))
    }
    if (length(group_values) != 2L) {
        stop("`group` must have exactly two levels")
    }
    base_data <- unique(raw_data[, setdiff(names(raw_data), group), drop = FALSE])
    X0data <- X1data <- base_data
    X0data[[group]] <- rep(group_values[1L], nrow(base_data))
    X1data[[group]] <- rep(group_values[2L], nrow(base_data))
    m0 <- model.frame(Terms, X0data, xlev = xlevels, na.action = na.pass)
    m1 <- model.frame(Terms, X1data, xlev = xlevels, na.action = na.pass)
    X0 <- model.matrix(Terms, m0, contrasts.arg = contrasts)
    X1 <- model.matrix(Terms, m1, contrasts.arg = contrasts)
    same_cols <- colSums(abs(X0 - X1)) == 0
    list(X0 = X0,
         X1 = X1,
         Xbase = X0[, same_cols, drop = FALSE],
         group = group,
         group_values = group_values)
}

ordsup <- function(theta, data, formula, group, po = TRUE,
                         measure = c("gamma", "Delta"),
                         contrasts = NULL) {
    measure <- match.arg(measure)
    XX <- make_X0_X1(formula = formula, data = data, group = group,
                     contrasts = contrasts)
    X0 <- XX$X0
    X1 <- XX$X1
    X0 <- as.matrix(X0)
    X1 <- as.matrix(X1)
    if (!identical(dim(X0), dim(X1))) {
        stop("`X0` and `X1` must have the same dimensions")
    }
    if (!identical(colnames(X0), colnames(X1))) {
        stop("`X0` and `X1` must have the same column names")
    }
    ncat <- infer_ncat_acl(theta, X0, po = po)
    probs0 <- acl_probs(theta, X0, ncat = ncat, po = po)
    probs1 <- acl_probs(theta, X1, ncat = ncat, po = po)
    gamma_fun <- function(p0, p1) {
        out <- outer(p0, p1, "*")
        sum(out[upper.tri(out)]) + sum(diag(out)) / 2
    }
    gamma <- vapply(seq_len(nrow(X0)),
                    function(i) gamma_fun(probs0[i, ], probs1[i, ]),
                    numeric(1))
    if (measure == "gamma") gamma else 2 * gamma - 1
}

on_est <- function(theta, on, components, estimator = "ML", ...) {
    out <- lapply(c("no", "median", "mean"), function(corr) {
        focus_engine(theta,
                     components = components,
                     on = on,
                     correction = corr,
                     estimator = estimator,
                     ...)
    })
    names(out) <- c("plugin", "medianBR", "meanBR")
    out
}

summarize_results <- function(object, truth) {
    data.frame(pu = rowMeans(sapply(object, function(x) x["estimate", ]) < truth),
               meanBias = rowMeans(sapply(object, function(x) x["estimate", ]) - truth),
               abs_meanBias = rowMeans(abs(sapply(object, function(x) x["estimate", ]) - truth)),
               rmse = sqrt(rowMeans((sapply(object, function(x) x["estimate", ]) - truth)^2)),
               coverage = rowMeans(sapply(object, function(x) (x["lower", ] < truth) & (x["upper", ]) > truth)))
}


data("wine", package = "ordinal")
c_fit <- bracl(rating ~ temp * contact,
               data = wine,
               parallel = FALSE,
               type = "AS_mean")

osm <- ordinal_superiority(c_fit, ~ contact, bc = FALSE)
theta <- coef(c_fit)
fm <- formula(c_fit)
po <- c_fit$parallel

## Some tests on the correctness of loglik, score, info
expect_equal(loglik(theta, wine, fm, po = po), unclass(logLik(c_fit)), check.attributes = FALSE)
expect_equal(score(theta, data = wine, formula = fm, po = po),
             numDeriv::grad(loglik, theta, data = wine, formula = fm, po = po),
             check.attributes = FALSE)
expect_equal(ordsup(theta, wine, fm, "contact", po = po), osm[, "gamma"], check.attributes = FALSE)
expect_equal(info(theta, data = wine, formula = fm, po = po),
             -numDeriv::hessian(loglik, theta, data = wine, formula = fm, po = po),
             check.attributes = FALSE, tol = 1e-05)
expect_equal(info(theta, data = wine, formula = fm, po = po) |> solve(),
             vcov(c_fit))

## FEF
R <- 10000
nsim <- 500
set.seed(111)
with_progress({
    pr <- progressor(R)
    res_fef <- lapply(1:R, function(i) {
        pr()
        library("brglm2")
        cdat <- simulate_acl(theta, wine, fm, po = po)
        c_coef <- rep(0, length(theta))
        names(c_coef) <- names(theta)
        coefs <- update(c_fit, data = cdat) |> coef()
        c_coef[names(coefs)] <- coefs
        comps <- estimate_focus_components_fef(c_coef,
                                               data = wine,
                                               loglik,
                                               score,
                                               info,
                                               simulate = simulate_acl,
                                               nsim = nsim,
                                               parallelize = TRUE,
                                               simu_data = wine,
                                               formula = fm,
                                               po = po)
        list(estimate = c_coef, components = comps)
    })
})

simu_res <- vector("list", 2)
for (wh in 1:2) {
    with_progress({
        progr <- progressor(length(res_fef))
        simu_res[[wh]] <- future_lapply(res_fef, function(x) {
            progr()
            co <- x$components
            th <- x$estimate
            sapply(on_est(th,
                          on = function(theta, ...) ordsup(theta, ...)[wh],
                          components = co,
                          estimator = "meanBR",
                          formula = fm,
                          group = "contact",
                          data = wine,
                          po = po), function(x) {
                c(estimate = coef(x), confint(x))
            })
        })
    })
}

summaries <- vector("list", 2)
for (wh in 1:2) {
    summaries[[wh]] <- summarize_results(simu_res[[wh]], osm[wh, "gamma"])
}

summaries[[1]] <- summaries[[1]] |>
    transform(temperature = "cold",
              method = rownames(summaries[[1]]))
summaries[[2]] <- summaries[[2]] |>
    transform(temperature = "warm",
              method = rownames(summaries[[1]]))

summaries <- do.call("rbind", summaries)
rownames(summaries) <- NULL

save(summaries, file = results_file)
