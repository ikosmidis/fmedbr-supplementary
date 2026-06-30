base_dir <- "."
results_file <- file.path(base_dir, "results/marginal-effects.rda")

pkgs <- c("enrichwith", "brglm2", "focuson", "future.apply", "tinytest", "ISLR2")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 1)
plan(multisession, workers = n_cores)

model_matrix_derivative <- function(model, data, variable,
                                    method = "Richardson") {
    stopifnot(nrow(data) == 1)
    tt <- delete.response(terms(model))
    f <- function(z) {
        data_z <- data
        data_z[[variable]] <- z
        drop(model.matrix(
            tt,
            data_z,
            contrasts.arg = model$contrasts,
            xlev = model$xlevels
        ))
    }
    z0 <- data[[variable]]
    J <- numDeriv::jacobian(f, z0, method = method)
    dx0 <- drop(J)
    names(dx0) <- names(coef(model))
    dx0
}

model_matrix_contrast <- function(model, data, variable, from, to) {
    stopifnot(nrow(data) == 1)
    tt <- terms(model)
    data_from <- data
    data_to <- data
    data_from[[variable]] <- from
    data_to[[variable]] <- to
    X_from <- model.matrix(tt, data_from, contrasts.arg = model$contrasts)
    X_to <- model.matrix(tt, data_to, contrasts.arg = model$contrasts)
    X_to - X_from
}

me <- function(theta, x0, dx0, link) {
    eta0 <- sum(theta * x0)
    d0 <- link$mu.eta(eta0)
    unname(sum(theta * dx0) * d0)
}

me_gradient <- function(theta, x0, dx0, link) {
    ## link <- enrich(link)
    eta0 <- sum(theta * x0)
    a <- sum(theta * dx0)
    d1 <- link$mu.eta(eta0)
    d2 <- link$d2mu.deta(eta0)
    unname(drop(dx0 * d1 + a * d2 * x0))
}

me_hessian <- function(theta, x0, dx0, link) {
    ## link <- enrich(link)
    x0 <- as.numeric(x0)
    dx0 <- as.numeric(dx0)
    eta0 <- sum(theta * x0)
    a <- sum(theta * dx0)
    d2 <- link$d2mu.deta(eta0)
    d3 <- link$d3mu.deta(eta0)
    H <- d2 * (outer(dx0, x0) + outer(x0, dx0)) +
        a * d3 * outer(x0, x0)
    unname(H)
}



run_experiment <- function(fm, variable, newdata, R = 1000, hulc_ci = FALSE, alpha = 0.05,
                           separation = TRUE, br_base = FALSE) {
    stopifnot(nrow(newdata) == 1)
    x0 <- model.matrix(delete.response(fm$terms), data = newdata)
    dx0 <- model_matrix_derivative(fm, data = newdata, variable = variable)
    simu_fun <- get_simulate_function(fm)
    link <- enrich(make.link(fm$family$link))
    psi <- me(coef(fm), x0 = x0, dx0 = dx0, link = link)
    hulc_ci <- isTRUE(hulc_ci)
    focus_fun <- function(corr, model) {
        focus(model,
              on = me,
              on_gradient = me_gradient,
              on_hessian = me_hessian,
              correction = corr,
              x0 = x0,
              dx0 = dx0,
              link = link)
    }
    simulated_responses <- simu_fun(nsim = R)
    res <- future_lapply(1:R, function(j) {
        cdat <- model.frame(fm)
        cdat$default <- simulated_responses[, j]
        cmod <- update(fm, data = cdat)
        cmod$call$data <- cdat
        method <- c("no", "mean", "median")
        if (isTRUE(br_base)) {
            brcmod <- update(cmod, method = "brglm_fit", type = "AS_mean")
            objects <- list("no" = focus_fun("no", cmod),
                            "mean" = focus_fun("mean", brcmod),
                            "median" = focus_fun("median", brcmod))
        } else {
            objects <- lapply(method,
                              focus_fun,
                              model = cmod)
        }

        se_comp <- sapply(objects, function(o) {
            sec <- try(focus_se(o, control = list(tol_opt = 1e-04))$se, silent = TRUE)
            if (inherits(sec, "try-error")) {
                NA
            } else {
                sec
            }
        })
        names(objects) <- method
        out <- list(estimate = sapply(objects, coef),
                    se = sapply(objects, "[[", "se"),
                    se_comp = se_comp)
        if (hulc_ci) {
            attr(out, "hulc_ci_no") <- confint(objects[["no"]], method = "hulc", level = 1 - alpha)
            attr(out, "hulc_ci_mean") <- confint(objects[["mean"]], method = "hulc", level = 1 - alpha)
            attr(out, "hulc_ci_median") <- confint(objects[["median"]], method = "hulc", level = 1 - alpha)
        }
        if (isTRUE(separation)) {
            attr(out, "separation") <- update(cmod,
                                              method = detectseparation::detect_separation)$outcome
        } else {
            attr(out, "separation") <- NA
        }
        out
    }, future.seed = TRUE, future.packages = c("brglm2", "focuson"))
    attr(res, "R") <- R
    attr(res, "alpha") <- alpha
    attr(res, "model") <- fm
    attr(res, "newdata") <- newdata
    attr(res, "hulc_ci") <- isTRUE(hulc_ci)
    attr(res, "variable") <- variable.names
    attr(res, "psi") <- psi
    attr(res, "separation") <- isTRUE(separation)
    res
}

summarize_results <- function(object) {
    psi <- attr(object, "psi")
    alpha <- attr(object, "alpha")
    estimates <- sapply(object, "[[", "estimate")
    ses <- sapply(object, "[[", "se")
    ses_comp <- sapply(object, "[[", "se_comp")
    mean_bias <- rowMeans(estimates - psi, na.rm = TRUE)
    mean_abs_bias <- rowMeans(abs(estimates - psi), na.rm = TRUE)
    pu <- rowMeans(estimates < psi, na.rm = TRUE)
    mse <- rowMeans((estimates - psi)^2, na.rm = TRUE)
    cover_wald <- rowMeans(abs(estimates - psi) / ses < qnorm(1 - alpha / 2), na.rm = TRUE)
    cover_wald_comp <- rowMeans(abs(estimates - psi) / ses_comp < qnorm(1 - alpha / 2), na.rm = TRUE)
    if (attr(object, "separation")) {
        separation <- sapply(object, attr, "separation")
    } else {
        separation <- FALSE
    }
    if (attr(object, "hulc_ci")) {
        check_sign <- function(x) all(x == c(-1, 1))
        hulc_no <- sapply(object, attr, "hulc_ci_no")
        hulc_mean <- sapply(object, attr, "hulc_ci_mean")
        hulc_median <- sapply(object, attr, "hulc_ci_median")
        ch_no <- mean(apply(sign(hulc_no - psi), 2, check_sign), na.rm = TRUE)
        ch_mean <- mean(apply(sign(hulc_mean - psi), 2, check_sign), na.rm = TRUE)
        ch_median <- mean(apply(sign(hulc_median - psi), 2, check_sign), na.rm = TRUE)
        cover_hulc <- c(ch_no, ch_mean, ch_median)
    } else {
        cover_hulc <- rep(NA, length(cover_wald))
    }
    ## HULC
    out <- data.frame(
        R = attr(object, "R"),
        psi = psi,
        mean_bias = mean_bias,
        mean_abs_bias = mean_abs_bias,
        pu = pu,
        mse = mse,
        cover_wald = cover_wald,
        cover_wald_comp = cover_wald_comp,
        cover_hulc = cover_hulc,
        separation = mean(separation))
    out
}


data("Default", package = "ISLR2")
Nsub <- nrow(Default) / 4
R <- 10000
def_data <- Default |>
    transform(balance = balance / 10000,
              income = income / 10000)
set.seed(123)
def_data <- def_data[sample(nrow(Default), Nsub), ]

## Get representative settings for individual marginal effects for a
## student and a non-student
summ_dat <- aggregate(cbind(income, balance) ~ student, data = def_data, mean)

## ML fit to simualate from
def_mod <- glm(default ~ student * (balance + income), family = binomial(probit),
               data = def_data)

## A bit inefficient as we refit the same models multiple times per
## row of summ_dat but OK here
set.seed(123)
results1_ml_base <- run_experiment(def_mod,
                                   variable = "balance",
                                   newdata = summ_dat[1, ],
                                   R = R,
                                   hulc_ci = TRUE,
                                   alpha = 0.05,
                                   separation = TRUE,
                                   br_base = FALSE)
set.seed(123)
results2_ml_base <- run_experiment(def_mod,
                                   variable = "balance",
                                   newdata = summ_dat[2, ],
                                   R = R,
                                   hulc_ci = TRUE,
                                   alpha = 0.05,
                                   separation = TRUE,
                                   br_base = FALSE)
set.seed(123)
results1_br_base <- run_experiment(def_mod,
                                   variable = "balance",
                                   newdata = summ_dat[1, ],
                                   R = R,
                                   hulc_ci = TRUE,
                                   alpha = 0.05,
                                   separation = TRUE,
                                   br_base = TRUE)
set.seed(123)
results2_br_base <- run_experiment(def_mod,
                                   variable = "balance",
                                   newdata = summ_dat[2, ],
                                   R = R,
                                   hulc_ci = TRUE,
                                   alpha = 0.05,
                                   separation = TRUE,
                                   br_base = TRUE)

summaries1_ml_base <- summarize_results(results1_ml_base) |> transform(student = "No", base = "ML")
summaries2_ml_base <- summarize_results(results2_ml_base) |> transform(student = "Yes", base = "ML")
summaries1_br_base <- summarize_results(results1_br_base) |> transform(student = "No", base = "mean BR")
summaries2_br_base <- summarize_results(results2_br_base) |> transform(student = "Yes", base = "mean BR")
summaries1_ml_base$method <- rownames(summaries1_ml_base)
summaries2_ml_base$method <- rownames(summaries2_ml_base)
summaries1_br_base$method <- rownames(summaries1_br_base)
summaries2_br_base$method <- rownames(summaries2_br_base)
summaries <- rbind(summaries1_ml_base, summaries2_ml_base, summaries1_br_base, summaries2_br_base)
rownames(summaries) <- NULL

save(summaries, file = results_file)


## Tests
if (FALSE) {
    library("marginaleffects")
    id <- 2
    vari <- "balance"
    ndat <- summ_dat[id, ]
    ll <- enrich(make.link(def_mod$family$link))
    mod_slopes <- slopes(def_mod, variables = vari, newdata = ndat)
    x0 <- model.matrix(delete.response(def_mod$terms), data = ndat)
    dx0 <- model_matrix_derivative(def_mod, data = ndat, variable = vari)
    expect_equal(mod_slopes[1, "estimate"], me(coef(def_mod), x0 = x0, dx0 = dx0, link = ll), tolerance = 1e-06)
    expect_equal(numDeriv::hessian(me, coef(def_mod), x0 = x0, dx0 = dx0, link = ll),
                 me_hessian(coef(def_mod), x0 = x0, dx0 = dx0, link = ll))
    expect_equal(numDeriv::grad(me, coef(def_mod), x0 = x0, dx0 = dx0, link = ll),
                 me_gradient(coef(def_mod), x0 = x0, dx0 = dx0, link = ll),
                 check.attributes = FALSE)
    a0 <- focus(def_mod, on = me, on_gradient = me_gradient, on_hessian = me_hessian, correction = "no", x0 = x0, dx0 = dx0, link = ll)
    a1 <- focus(def_mod, on = me, on_gradient = me_gradient, on_hessian = me_hessian, correction = "median", x0 = x0, dx0 = dx0, link = ll)
    a_mean <- focus(def_mod, on = me, on_gradient = me_gradient, on_hessian = me_hessian, correction = "median", x0 = x0, dx0 = dx0, link = ll, se_at = "corrected")
}
