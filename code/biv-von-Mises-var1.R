base_dir <- "."
results_file <- file.path(base_dir, "results/bvmsin-circular-variance.rda")

pkgs <- c("BAMBI", "focuson", "numDeriv", "future.apply")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 1)
plan(multisession, workers = n_cores)

## Focus: circular variance of first coordinate
circvar_on <- function(theta, ...) {
    BAMBI::circ_varcor_model(model = "vmsin",
                             kappa1 = exp(theta[3]),
                             kappa2 = exp(theta[4]),
                             kappa3 = theta[5],
                             mu1 = theta[1],
                             mu2 = theta[2])$var1
}

simu_n <- function(theta, n = 10) {
    BAMBI::rvmsin(n = n,
                  kappa1 = exp(theta[3]),
                  kappa2 = exp(theta[4]),
                  kappa3 = theta[5],
                  mu1 = theta[1],
                  mu2 = theta[2])
}

simu_one <- function(theta) {
    drop(simu_n(theta, n = 1))
}

mle_theta <- function(data) {
    fit <- BAMBI::vm2_mle(data, model = "vmsin")
    cf <- fit@coef
    c(mu1 = cf[["mu1"]],
      mu2 = cf[["mu2"]],
      log_kappa1 = log(cf[["kappa1"]]),
      log_kappa2 = log(cf[["kappa2"]]),
      lambda = cf[["kappa3"]])
}

loglik_one <- function(theta, data) {
    BAMBI::dvmsin(x = data,
                  kappa1 = exp(theta[3]),
                  kappa2 = exp(theta[4]),
                  kappa3 = theta[5],
                  mu1 = theta[1],
                  mu2 = theta[2],
                  log = TRUE)
}

estimate_components <- function(theta, n, nsim = 1000, seed = NULL,
                                parallelize = TRUE) {
    if (!is.null(seed)) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
        set.seed(seed)
    }
    focuson::estimate_focus_components_iid(
                 theta = theta,
                 n = n,
                 loglik = loglik_one,
                 simulate = simu_one,
                 nsim = nsim,
                 parallelize = parallelize,
                 diagnostics = TRUE)
}

#' ... further arguments to be passed to estimate_components
focus_object <- function(data, correction = "all", ...) {
    theta_hat <- mle_theta(data)
    if (any(!is.finite(theta_hat))) return(NULL)
    components <- estimate_components(theta_hat, n = nrow(data), ...)
    methods <- c("no", "mean", "median", "all")
    vmethods <- methods[methods != "all"]
    selected <- if (identical(correction, "all")) vmethods else correction
    focus_fun <- function(corr) {
        focuson::focus_engine(theta = theta_hat,
                              components = components,
                              on = circvar_on,
                              correction = corr,
                              estimator = "ML")
    }
    out <- lapply(selected, focus_fun)
    names(out) <- selected
    out
}

hulc <- function(data, correction ="no", level = 0.95, ...) {
    statistic <- function(df) {
        focus_object(df, correction = correction, ...)[[1]] |> coef()
    }
    focuson::hulc_ci(data.frame(data), statistic = statistic, level = level)
}

run_experiment <- function(theta, n = 50, R = 1000, level = 0.95, nsim = 1000) {
    psi <- circvar_on(theta)
    res <- future_lapply(1:R, function(j) {
        dat <- simu_n(theta, n = n)
        fos <- focus_object(dat, correction = "all", nsim = nsim, seed = 123,
                            parallelize = FALSE)
        estimates <- sapply(fos, coef)
        wald_cis <- sapply(fos, confint, level = level)
        for (i in 1:5) {
            hulc_cis <- sapply(c("no", "mean", "median"), function(co)
                hulc(dat, correction = co, level = level, nsim = nsim, seed = 123,
                     parallelize = FALSE))
            if (all(!is.na(hulc_cis))) break
        }
        list(estimates = estimates,
             wald_cis = wald_cis,
             hulc_cis = hulc_cis)
    }, future.seed = TRUE)
    attr(res, "psi") <- psi
    attr(res, "theta") <- theta
    attr(res, "n") <- n
    attr(res, "R") <- R
    attr(res, "level") <- level
    attr(res, "nsim") <- nsim
    res
}


summarize_results <- function(object) {
    psi <- attr(object, "psi")
    cover <- function(cis) {
        rowMeans(sapply(cis, function(z)
            apply(sign(z - psi) == c(-1, 1), 2, all)),
            na.rm = TRUE)
    }
    estimates <- sapply(object, "[[", "estimates")
    wald_cis <- lapply(object, "[[", "wald_cis")
    hulc_cis <- lapply(object, "[[", "hulc_cis")
    out <- data.frame(
        method = c("ML", "meanBR", "medianBR"),
        R = attr(object, "R"),
        n = attr(object, "n"),
        psi = psi,
        mean_bias = rowMeans(estimates - psi, na.rm = TRUE),
        mean_abs_bias = rowMeans(abs(estimates - psi), na.rm = TRUE),
        pu = rowMeans(estimates < psi, na.rm = TRUE),
        mse = rowMeans((estimates - psi)^2, na.rm = TRUE),
        cover_wald = cover(wald_cis),
        cover_hulc = cover(hulc_cis)
    )
    attr(out, "theta") <- attr(object, "theta")
    attr(out, "level") <- attr(object, "level")
    attr(out, "nsim") <- attr(object, "nsim")
    out
}

theta <- c(1, 2, 1, 1, 2)
set.seed(123)
dd <- simu_n(theta, 100)
focus_object(dd, seed = 123)

settings <- expand.grid(mu1 = pi/3,
                        mu2 = 2 * pi/3,
                        kappa1 = 2.5,
                        kappa2 = 3,
                        lambda = 1,
                        n = c(30, 60, 90, 120),
                        R = 5000,
                        level = 0.95,
                        nsim = 500)

summaries <- vector("list", nrow(settings))
set.seed(111)
for (s in seq.int(nrow(settings))) {
    theta <- c(mu1 = settings[s, "mu1"],
               mu2 = settings[s, "mu2"],
               log_kappa1 = log(settings[s, "kappa1"]),
               log_kappa2 = log(settings[s, "kappa2"]),
               lambda = settings[s, "lambda"])
    n <- settings[s, "n"]
    R <- settings[s, "R"]
    level <- settings[s, "level"]
    nsim <- settings[s, "nsim"]
    res <- run_experiment(theta = theta,
                          n = n,
                          R = R,
                          level = level,
                          nsim = nsim)
    summaries[[s]] <- summarize_results(res)
    cat("Setting:", "n =", n, "R =", R,
        "nsim =", nsim, "Done.\n")
}

summaries <- do.call("rbind", summaries)
rownames(summaries) <- NULL

## Save results
save(settings, summaries, file = results_file)




