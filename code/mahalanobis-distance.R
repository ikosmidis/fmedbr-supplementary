base_dir <- "."
results_file <- file.path(base_dir, "results/mahalanobis-distance.rda")

pkgs <- c("mvtnorm", "future.apply", "focuson")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 1)
plan(multisession, workers = n_cores)

simu_fun <- function(mu = c(0, 0), sigma = diag(2), n = 10) {
    data.frame(mvtnorm::rmvnorm(n = n, mean = mu, sigma = sigma))
}

maha <- function(data, mu0, what = "all") {
    what <- match.arg(what, c("all", "ML", "meanBR", "medianBR"))
    n <- nrow(data)
    p <- ncol(data)
    mu <- colMeans(data)
    sigma <- var(data) * (n - 1) / n
    phi_ml <- drop(t(mu - mu0) %*% solve(sigma) %*% (mu - mu0))
    if (isTRUE(what == "ML")) return(phi_ml)
    lam <- n * phi_ml
    m1 <- p + lam
    v <- n - p
    m2 <- m1^2 + 2 * (p + 2 * lam)
    m3 <- m1^3 + 6 * m1 * (p + 2 * lam) + 8 * (p + 3 * lam)
    b <- m1 / (v - 2) - phi_ml
    K2 <- m2 / ((v - 2) * (v - 4)) - m1^2 / (v - 2)^2
    K3 <- m3 / ((v - 2) * (v - 4) * (v - 6)) - 3 * m1 * m2 / ((v - 2)^2 * (v - 4)) + 2 * m1^3 / (v - 2)^3
    phi_meanBR <- phi_ml - b
    if (isTRUE(what == "meanBR")) return(phi_meanBR)
    phi_medianBR <- phi_meanBR + K3 / (6 * K2)
    if (isTRUE(what == "medianBR")) return(phi_medianBR)
    out <- c("ML" = phi_ml, "meanBR" = phi_meanBR, "medianBR" = phi_medianBR)
    out
}

maha_se <- function(phi, n, p) {
   lam <- n * phi
   m1 <- p + lam
   v <- n - p
   m2 <- m1^2 + 2 * (p + 2 * lam)
   K2 <- m2 / ((v - 2) * (v - 4)) - m1^2 / (v - 2)^2
   sqrt(K2)
}

AR1cor <- function(p, rho) {
    sapply(1:p, \(k) rho^abs(k - 1:p))
}

run_experiment <- function(mu0, mu, sigma, n = 100, R = 1000, hulc_ci = FALSE, wald_ci = FALSE, alpha = 0.05) {
    p <- length(mu)
    ## Checks for mu0, sigma
    out <- future_lapply(1:R, function(i) {
        dat <- simu_fun(mu = mu, sigma = sigma, n = n)
        res <- maha(dat, mu0)
        if (isTRUE(wald_ci)) {
            quant <- qnorm(1 - alpha / 2)
            ses <- maha_se(res, n, p)
            attr(res, "wald_ci_ML") <- res["ML"] + c(-1, 1) * quant * ses["ML"]
            attr(res, "wald_ci_meanBR") <- res["meanBR"] + c(-1, 1) * quant * ses["meanBR"]
            attr(res, "wald_ci_medianBR") <- res["medianBR"] + c(-1, 1) * quant * ses["medianBR"]
        }
        if (isTRUE(hulc_ci)) {
            attr(res, "hulc_ci_ML") <-
                focuson::hulc_ci(dat, maha, mu0 = mu0, what = "ML", level = 1 - alpha)
            attr(res, "hulc_ci_meanBR") <-
                focuson::hulc_ci(dat, maha, mu0 = mu0, what = "meanBR", level = 1 - alpha)
            attr(res, "hulc_ci_medianBR") <-
                focuson::hulc_ci(dat, maha, mu0 = mu0, what = "medianBR", level = 1 - alpha)
        }
        res
    }, future.seed = TRUE)
    attr(out, "phi") <- drop(t(mu - mu0) %*% solve(sigma) %*% (mu - mu0))
    attr(out, "n") <- n
    attr(out, "p") <- p
    attr(out, "R") <- R
    attr(out, "alpha") <- alpha
    attr(out, "mu") <- mu
    attr(out, "sigma") <- sigma
    attr(out, "mu0") <- mu0
    attr(out, "hulc_ci") <- isTRUE(hulc_ci)
    attr(out, "wald_ci") <- isTRUE(wald_ci)
    out
}

summarize_results <- function(object) {
    phi <- attr(object, "phi")
    n <- attr(object, "n")
    p <- attr(object, "p")
    alpha <- attr(object, "alpha")
    estimates <- do.call("cbind", object)
    mean_bias <- rowMeans(estimates - phi, na.rm = TRUE)
    mean_abs_bias <- rowMeans(abs(estimates - phi), na.rm = TRUE)
    pu <- rowMeans(estimates < phi, na.rm = TRUE)
    mse <- rowMeans((estimates - phi)^2, na.rm = TRUE)
    ## We compute Wald coverage always
    cover_wald <- rowMeans(abs(estimates - phi) / maha_se(estimates, n, p) < qnorm(1 - alpha / 2), na.rm = TRUE)
    if (attr(object, "hulc_ci")) {
        check_sign <- function(x) all(x == c(-1, 1))
        hulc_ml <- sapply(object, attr, "hulc_ci_ML")
        hulc_meanBR <- sapply(object, attr, "hulc_ci_meanBR")
        hulc_medianBR <- sapply(object, attr, "hulc_ci_medianBR")
        ch_ml <- mean(apply(sign(hulc_ml - phi), 2, check_sign), na.rm = TRUE)
        ch_meanBR <- mean(apply(sign(hulc_meanBR - phi), 2, check_sign), na.rm = TRUE)
        ch_medianBR <- mean(apply(sign(hulc_medianBR - phi), 2, check_sign), na.rm = TRUE)
        cover_hulc <- c(ch_ml, ch_meanBR, ch_medianBR)
    } else {
        cover_hulc <- rep(NA, length(cover_wald))
    }
    out <- data.frame(
        n = n,
        p = p,
        R = R,
        phi = phi,
        mean_bias = mean_bias,
        mean_abs_bias = mean_abs_bias,
        pu = pu,
        mse = mse,
        cover_wald = cover_wald,
        cover_hulc = cover_hulc)
    out$method <- rownames(out)
    attr(out, "phi") <- attr(object, "phi")
    attr(out, "n") <- attr(object, "n")
    attr(out, "p") <- attr(object, "p")
    attr(out, "R") <- attr(object, "R")
    attr(out, "mu") <- attr(object, "mu")
    attr(out, "sigma") <- attr(object, "sigma")
    attr(out, "mu0") <- attr(object, "mu0")
    attr(out, "alpha") <- alpha
    out
}

settings <- expand.grid(p = c(10, 20, 30),
                        n = 2^(7:10),
                        R = 1000000)

summaries <- vector("list", nrow(settings))
set.seed(111)
for (s in seq.int(nrow(settings))) {
    p <- settings[s, "p"]
    n <- settings[s, "n"]
    R <- settings[s, "R"]
    mu <- rep(0, p)
    sigma <- AR1cor(p, 0.5)
    mu0 <- (1:p) / p
    res <- run_experiment(mu0 = mu0,
                          mu = mu,
                          sigma = sigma,
                          n = n,
                          R = R)
    summ <- summarize_results(res)
    summaries[[s]] <- summ
    cat("Setting:", "n =", n, "p =", p, "R =", R, "Done.\n")
}

summaries <- do.call("rbind", summaries)
rownames(summaries) <- NULL

## Save results
save(settings, summaries, file = results_file)

