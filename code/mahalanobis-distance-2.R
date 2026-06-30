base_dir <- "."
results_file <- file.path(base_dir, "results/mahalanobis-distance-2sample.rda")

pkgs <- c("mvtnorm", "future.apply", "focuson")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 1)
plan(multisession, workers = n_cores)

simu_fun <- function(mu1 = c(0, 0), mu2 = c(0, 0), sigma = diag(2),
                     n1 = 10, n2 = 10) {
  list(
    Y1 = data.frame(mvtnorm::rmvnorm(n = n1, mean = mu1, sigma = sigma)),
    Y2 = data.frame(mvtnorm::rmvnorm(n = n2, mean = mu2, sigma = sigma))
  )
}

maha <- function(data, what = "all") {
  what <- match.arg(what, c("all", "ML", "meanBR", "medianBR"))
  Y1 <- data$Y1
  Y2 <- data$Y2
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  n <- n1 + n2
  p <- ncol(Y1)

  mu1 <- colMeans(Y1)
  mu2 <- colMeans(Y2)
  d <- mu1 - mu2

  Y1c <- sweep(Y1, 2, mu1, "-")
  Y2c <- sweep(Y2, 2, mu2, "-")
  W <- crossprod(as.matrix(Y1c)) + crossprod(as.matrix(Y2c))

  sigma <- W / n
  phi_ml <- drop(t(d) %*% solve(sigma) %*% d)
  if (isTRUE(what == "ML")) return(phi_ml)

  cst <- n^2 / (n1 * n2)
  lam <- (n1 * n2 / n) * phi_ml
  m1 <- p + lam
  v <- n - p - 1
  m2 <- m1^2 + 2 * (p + 2 * lam)
  m3 <- m1^3 + 6 * m1 * (p + 2 * lam) + 8 * (p + 3 * lam)

  b <- cst * m1 / (v - 2) - phi_ml

  K2 <- cst^2 * (
    m2 / ((v - 2) * (v - 4)) -
      m1^2 / (v - 2)^2
  )

  K3 <- cst^3 * (
    m3 / ((v - 2) * (v - 4) * (v - 6)) -
      3 * m1 * m2 / ((v - 2)^2 * (v - 4)) +
      2 * m1^3 / (v - 2)^3
  )

  phi_meanBR <- phi_ml - b
  if (isTRUE(what == "meanBR")) return(phi_meanBR)

  phi_medianBR <- phi_meanBR + K3 / (6 * K2)
  if (isTRUE(what == "medianBR")) return(phi_medianBR)

  out <- c("ML" = phi_ml, "meanBR" = phi_meanBR, "medianBR" = phi_medianBR)
  out
}

maha_se <- function(phi, n1, n2, p) {
  n <- n1 + n2
  cst <- n^2 / (n1 * n2)
  lam <- (n1 * n2 / n) * phi
  m1 <- p + lam
  v <- n - p - 1
  m2 <- m1^2 + 2 * (p + 2 * lam)

  K2 <- cst^2 * (
    m2 / ((v - 2) * (v - 4)) -
      m1^2 / (v - 2)^2
  )

  sqrt(K2)
}

AR1cor <- function(p, rho) {
  sapply(1:p, \(k) rho^abs(k - 1:p))
}

run_experiment <- function(mu1, mu2, sigma, n1 = 50, n2 = 50, R = 1000,
                           hulc_ci = FALSE, wald_ci = FALSE,
                           alpha = 0.05) {
  p <- length(mu1)

  out <- future_lapply(1:R, function(i) {
    dat <- simu_fun(mu1 = mu1, mu2 = mu2, sigma = sigma, n1 = n1, n2 = n2)
    res <- maha(dat)

    if (isTRUE(wald_ci)) {
      quant <- qnorm(1 - alpha / 2)
      ses <- maha_se(res, n1, n2, p)
      attr(res, "wald_ci_ML") <- res["ML"] + c(-1, 1) * quant * ses["ML"]
      attr(res, "wald_ci_meanBR") <- res["meanBR"] + c(-1, 1) * quant * ses["meanBR"]
      attr(res, "wald_ci_medianBR") <- res["medianBR"] + c(-1, 1) * quant * ses["medianBR"]
    }

    if (isTRUE(hulc_ci)) {
      attr(res, "hulc_ci_ML") <-
        focuson::hulc_ci(dat, maha, what = "ML", level = 1 - alpha)
      attr(res, "hulc_ci_meanBR") <-
        focuson::hulc_ci(dat, maha, what = "meanBR", level = 1 - alpha)
      attr(res, "hulc_ci_medianBR") <-
        focuson::hulc_ci(dat, maha, what = "medianBR", level = 1 - alpha)
    }

    res
  }, future.seed = TRUE)

  attr(out, "phi") <- drop(t(mu1 - mu2) %*% solve(sigma) %*% (mu1 - mu2))
  attr(out, "n1") <- n1
  attr(out, "n2") <- n2
  attr(out, "n") <- n1 + n2
  attr(out, "p") <- p
  attr(out, "R") <- R
  attr(out, "alpha") <- alpha
  attr(out, "mu1") <- mu1
  attr(out, "mu2") <- mu2
  attr(out, "sigma") <- sigma
  attr(out, "hulc_ci") <- isTRUE(hulc_ci)
  attr(out, "wald_ci") <- isTRUE(wald_ci)

  out
}

summarize_results <- function(object) {
  phi <- attr(object, "phi")
  n1 <- attr(object, "n1")
  n2 <- attr(object, "n2")
  n <- attr(object, "n")
  p <- attr(object, "p")
  R <- attr(object, "R")
  alpha <- attr(object, "alpha")

  estimates <- do.call("cbind", object)
  mean_bias <- rowMeans(estimates - phi, na.rm = TRUE)
  mean_abs_bias <- rowMeans(abs(estimates - phi), na.rm = TRUE)
  pu <- rowMeans(estimates < phi, na.rm = TRUE)
  mse <- rowMeans((estimates - phi)^2, na.rm = TRUE)

  cover_wald <- rowMeans(
    abs(estimates - phi) / maha_se(estimates, n1, n2, p) < qnorm(1 - alpha / 2),
    na.rm = TRUE
  )

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
    n1 = n1,
    n2 = n2,
    n = n,
    p = p,
    R = R,
    phi = phi,
    mean_bias = mean_bias,
    mean_abs_bias = mean_abs_bias,
    pu = pu,
    mse = mse,
    cover_wald = cover_wald,
    cover_hulc = cover_hulc
  )

  out$method <- rownames(out)
  attr(out, "phi") <- attr(object, "phi")
  attr(out, "n1") <- attr(object, "n1")
  attr(out, "n2") <- attr(object, "n2")
  attr(out, "n") <- attr(object, "n")
  attr(out, "p") <- attr(object, "p")
  attr(out, "R") <- attr(object, "R")
  attr(out, "mu1") <- attr(object, "mu1")
  attr(out, "mu2") <- attr(object, "mu2")
  attr(out, "sigma") <- attr(object, "sigma")
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

  n1 <- n / 2
  n2 <- n / 2

  mu1 <- rep(0, p)
  sigma <- AR1cor(p, 0.5)
  mu2 <- (1:p) / p

  res <- run_experiment(mu1 = mu1,
                        mu2 = mu2,
                        sigma = sigma,
                        n1 = n1,
                        n2 = n2,
                        R = R)

  summ <- summarize_results(res)
  summaries[[s]] <- summ

  cat("Setting:", "n1 =", n1, "n2 =", n2, "p =", p, "R =", R, "Done.\n")
}

summaries <- do.call("rbind", summaries)
rownames(summaries) <- NULL

## Save results
save(settings, summaries, file = results_file)



