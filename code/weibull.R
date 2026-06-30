base_dir <- "."
results_file <- file.path(base_dir, "results/weibull.rda")

pkgs <- c("focuson", "pracma", "future.apply", "likelihoodAsy")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 2)
plan(multisession, workers = n_cores)

## P, Q, V
V <- function(theta, n) {
    tau <- unname(theta["tau"])
    g <- -digamma(1)
    Ainfo <- pi^2 / 6 + (1 - g)^2
    (6 / (n * pi^2)) * matrix(c(Ainfo / tau^2, 1 - g, 1 - g, tau^2), nrow = 2)
}

focus_components <- function(theta, n) {
    tau <- unname(theta["tau"])
    g <- -digamma(1)
    z3 <- pracma::zeta(3)
    Ainfo <- pi^2 / 6 + (1 - g)^2
    V <- (6 / (n * pi^2)) * matrix(c(Ainfo / tau^2, 1 - g, 1 - g, tau^2), nrow = 2)
    D  <- g^2 - 4 * g + 2 + pi^2 / 6
    E3 <- 2 * g^3 - 12 * g^2 + pi^2 * g + 12 * g -  2 * pi^2 - 2 + 4 * z3
    H  <- -g^3 + 5 * g^2 - (pi^2 / 2) * g - 4 * g - 2 * z3 + 5 * pi^2 / 6
    P_eta <- matrix(c(2 * n * tau^3, 2 * n * tau * (g - 2), 2 * n * tau * (g - 2), 2 * n * D / tau), nrow = 2)
    P_tau <- matrix(c(2 * n * tau * (g - 2), 2 * n * D / tau, 2 * n * D / tau, n * E3 / tau^3), nrow = 2)
    Q_eta <- matrix(c(-n * tau^3, n * tau * (3 - g), n * tau * (3 - g),  -n * D / tau), nrow = 2)
    Q_tau <- matrix(c(n * tau * (1 - g), -n * D / tau, -n * D / tau,  n * H / tau^3), nrow = 2)
    list(V = V, P = list(eta = P_eta, tau = P_tau), Q = list(eta = Q_eta, tau = Q_tau))
}


## Quantile focus
quant_phi <- function(eta, tau, alpha_q) {
    exp(eta + log(-log(alpha_q)) / tau)
}

quant_on <- function(theta, alpha_q = 0.01, ...) {
    quant_phi(theta[1], theta[2], alpha_q)
}

quant_grad <- function(theta, alpha_q = 0.01, ...) {
    eta <- theta[1]
    tau <- theta[2]
    c_alpha <- log(-log(alpha_q))
    phi <- exp(eta + c_alpha / tau)
    c(phi, -c_alpha * phi / tau^2)
}

quant_hess <- function(theta, alpha_q = 0.01, ...) {
    eta <- theta[1]
    tau <- theta[2]
    c_alpha <- log(-log(alpha_q))
    phi <- exp(eta + c_alpha / tau)
    H <- matrix(0, nrow = 2, ncol = 2)
    H[1, 1] <- phi
    H[1, 2] <- H[2, 1] <- -c_alpha * phi / tau^2
    H[2, 2] <- phi * (c_alpha^2 / tau^4 + 2 * c_alpha / tau^3)
    H
}


## MLE in eta = log(sigma), tau = alpha parameterization
mle_theta <- function(y, interval = c(0.05, 20), tol = 1e-10) {
  ly <- log(y)
  n <- length(y)
  nll_prof <- function(tau) {
    ## if (!is.finite(tau) || tau <= 0) return(Inf)
    yt <- y^tau
    ## if (any(!is.finite(yt))) return(Inf)
    m <- mean(yt)
    ## if (!is.finite(m) || m <= 0) return(Inf)
    eta_hat <- log(m) / tau
    - (n * log(tau) - n * tau * eta_hat + (tau - 1) * sum(ly) - n)
  }
  opt <- tryCatch(
      optimize(nll_prof, interval = interval, tol = tol),
      error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$minimum)) {
      return(c(eta = NA_real_, tau = NA_real_))
  }
  tau_hat <- opt$minimum
  eta_hat <- log(mean(y^tau_hat)) / tau_hat
  c(eta = eta_hat, tau = tau_hat)
}

focus_objects <- function(y, alpha_q = 0.01, correction = "all", ...) {
    theta_hat <- mle_theta(y, ...)
    N <- length(y)
    methods <- c("no", "mean", "median", "all")
    vmethods <- methods[methods != "all"]
    corrections <- correction[vmethods %in% correction]
    comps <- focus_components(theta_hat, N)
    focus_fun <- function(corr) {
        focus_engine(theta_hat,
                     components = comps,
                     on = quant_on,
                     on_gradient = quant_grad,
                     on_hessian = quant_hess,
                     alpha_q = alpha_q,
                     correction = corr,
                     estimator = "ML")
    }
    selected <- if (identical(correction, "all")) vmethods else correction
    out <- lapply(selected, focus_fun)
    names(out) <- selected
    out
}

hulc <- function(y, alpha_q = 0.01, correction = "no", level = 0.95, ...) {
    A <- function(df)
        coef(focus_objects(df$y, alpha_q = alpha_q, correction = correction, ...)[[1]])
    focuson::hulc_ci(data.frame(y = y), statistic = A, level = level)
}


## Need estimates and ci from rstar
mle_theta_b <- function(y, interval = c(0.05, 20), tol = 1e-10) {
    out <- mle_theta(y, interval = interval, tol = tol)
    out[2] <- log(out[2])
    names(out) <- c("eta", "xi")
    out
}

loglik_b <- function(theta, data) {
  y <- data$y
  n <- length(y)
  eta <- theta[1]
  xi  <- theta[2]
  tau <- exp(xi)
  ly <- log(y)
  z  <- ly - eta
  w  <- exp(tau * z)
  n * log(tau) - n * tau * eta + (tau - 1) * sum(ly) - sum(w)
}

score_b <- function(theta, data) {
  y <- data$y
  n <- length(y)
  eta <- theta[1]
  xi  <- theta[2]
  tau <- exp(xi)
  ly <- log(y)
  z  <- ly - eta
  w  <- exp(tau * z)
  score_eta <- tau * (sum(w) - n)
  score_tau <- n / tau + sum(z) - sum(w * z)
  c(score_eta, tau * score_tau)
}

gendat_b <- function(theta, data) {
  eta <- theta[1]
  tau <- exp(theta[2])
  data$y <- rweibull(length(data$y), shape = tau, scale = exp(eta))
  data
}

quant_b <- function(alpha_q = 0.01) {
    c_alpha <- log(-log(alpha_q))
    function(theta) {
        eta <- theta[1]
        xi  <- theta[2]
        eta + c_alpha * exp(-xi)
    }
}

compute_rstar <- function(data, alpha_q = 0.01, R = 200, seed = 123, ...) {
    out <- rstar.ci(data = data,
                    thetainit = mle_theta_b(data$y, ...),
                    floglik = loglik_b,
                    fscore = score_b,
                    fpsi = quant_b(alpha_q = alpha_q),
                    datagen = gendat_b,
                    R = R,
                    seed = seed,
                    trace = FALSE,
                    psidesc = "psi = log Weibull quantile")
    class(out) <- c("rs_res", class(out))
    out
}

confint.rs_res <- function(object, parm, level = 0.95) {
    alpha <- 1 - level
    idx <- match(round(alpha, 10), c(0.10, 0.05, 0.01))
    if (is.na(idx)) stop("level must be one of 0.90, 0.95, 0.99.")
    c(lower = unname(object$CIrs[idx, 1]),
      upper = unname(object$CIrs[idx, 2])) |> exp()
}

coef.rs_res <- function(object, ...) {
    spl <- smooth.spline(object$rsvals, object$psivals)
    exp(predict(spl, 0)$y)
}


####


run_experiment <- function(eta, tau, alpha_q = 0.01, n = 50, R = 1000, level = 0.95, ...) {
    quant <- unname(quant_on(c(eta = eta, tau = tau), alpha_q = alpha_q))
    res <- future_sapply(1:R, function(j) {
        y <- rweibull(n, shape = tau, scale = exp(eta))
        fos <- focus_objects(y, alpha_q = alpha_q, correction = "all", ...)
        rs <- compute_rstar(data.frame(y = y), alpha_q = alpha_q, R = 200, seed = NULL, ...)
        estimates <- c(sapply(fos, coef), "rs" = coef.rs_res(rs))
        cis_supp <- sapply(fos, confint, se_at = "supplied", level = level)
        cis_comp <- sapply(fos, confint, se_at = "compatible", level = level,
                           V_function = V, n = n)
        cis_rstar <- confint.rs_res(rs, level = level)
        cis_hulc <- sapply(c("no", "mean", "median"), function(co)
            hulc(y, alpha_q = alpha_q, correction = co, level = level))
        list(estimates = estimates, cis_supp = cis_supp, cis_comp = cis_comp,
             cis_hulc = cis_hulc,
             cis_rstar = cis_rstar)
    }, future.seed = TRUE, simplify = FALSE)
    attr(res, "quant") <- quant
    attr(res, "eta") <- eta
    attr(res, "tau") <- tau
    attr(res, "alpha_q") <- alpha_q
    attr(res, "n") <- n
    attr(res, "R") <- R
    attr(res, "level") <- level
    res
}

summarize_results <- function(object) {
    quant <- attr(object, "quant")
    co <- function(cis) {
        apply(sign(cis - quant) == c(-1, 1), 2, all)
    }
    cover <- function(cis) {
        rowMeans(sapply(cis, function(z) apply(sign(z - quant) == c(-1, 1), 2, all)))
    }
    estimates <- sapply(object, "[[", "estimates")
    cis_supp <- lapply(object, "[[", "cis_supp")
    cis_comp <- lapply(object, "[[", "cis_comp")
    cis_hulc <- lapply(object, "[[", "cis_hulc")
    cis_rstar <- sapply(object, "[[", "cis_rstar")
    data.frame(
        estimator = c("no", "mean", "median", "rs"),
        R = attr(object, "R"),
        n = attr(object, "n"),
        eta = attr(object, "eta"),
        tau = attr(object, "tau"),
        alpha_q = attr(object, "alpha_q"),
        level = attr(object, "level"),
        psi = quant,
        mean_bias = rowMeans(estimates - quant),
        mean_abs_bias = rowMeans(abs(estimates - quant)),
        mse = rowMeans((estimates - quant)^2),
        pu = rowMeans(estimates < quant),
        cover_wald_supp = c(cover(cis_supp), NA),
        cover_wald_comp = c(cover(cis_comp), NA),
        cover_hulc = c(cover(cis_hulc), NA),
        cover_rstar = c(NA, NA, NA, mean(co(cis_rstar)))
    )
}

settings <- expand.grid(eta = log(2),
                        tau = 1.5,
                        alpha_q = c(0.1, 0.05, 0.01),
                        R = 50000,
                        n = c(25, 50, 100, 200),
                        level = 0.95)

set.seed(123)
summ <- apply(settings, 1, function(set) {
    run_experiment(eta = set[["eta"]],
                   tau = set[["tau"]],
                   alpha_q = set[["alpha_q"]],
                   R = set[["R"]],
                   n = set[["n"]],
                   level = set[["level"]]) |>
        summarize_results()
}, simplify = FALSE)


summaries <- do.call("rbind", summ)
rownames(summaries) <- NULL
save(summaries, file = results_file)


