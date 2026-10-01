base_dir <- "."

exp_set <- "c"
## 0.125 (a: n1 = n/8); 0.09375 (b: n1 = 3n/32); 0.0625 (c: n1 = n/16)

pkgs <- c("future", "future.apply")
for (pkg in pkgs) library(pkg, character.only = TRUE)

n_cores <- parallel::detectCores()
if (is.na(n_cores)) n_cores <- 1
n_cores <- max(1, n_cores - 2)
plan(multisession, workers = n_cores)

## Multiple-mediator causal model.
## Focus on the joint indirect effect psi = nu'beta.

skew_constants <- function(lambda) {
    delta_SN <- lambda / sqrt(1 + lambda^2)
    mu_lambda <- delta_SN * sqrt(2 / pi)
    c_lambda <- sqrt(1 - mu_lambda^2)
    expectation <- function(fun) {
        integrand <- function(t) {
            r_lambda <- exp(dnorm(lambda * t, log = TRUE) -
                            pnorm(lambda * t, log.p = TRUE))
            v_lambda <- t - lambda * r_lambda
            K_lambda <- 1 + lambda^2 * r_lambda *
                (lambda * t + r_lambda)
            fun(v_lambda, K_lambda) * 2 * dnorm(t) * pnorm(lambda * t)
        }
        integrate(integrand, -12, 0, rel.tol = 1e-10)$value +
                                                        integrate(integrand, 0, 12, rel.tol = 1e-10)$value
    }
    J <- c_lambda^2 *
        expectation(function(v_lambda, K_lambda) v_lambda^2)
    C_B <- -c_lambda^3 *
        expectation(function(v_lambda, K_lambda) {
            v_lambda^3 - v_lambda * K_lambda
        }) / (2 * J^2)
    C_3 <- c_lambda^3 *
        expectation(function(v_lambda, K_lambda) {
            -2 * v_lambda^3 + 3 * v_lambda * K_lambda
        }) / J^3
    if (lambda == 0) {
        J <- 1; C_B <- C_3 <- 0
    }
    list(lambda = lambda, delta_SN = delta_SN,
         mu_lambda = mu_lambda, c_lambda = c_lambda,
         J = J, C_B = C_B, C_3 = C_3)
}

rskew <- function(n, sn) {
    T <- sn$delta_SN * abs(rnorm(n)) +
        sqrt(1 - sn$delta_SN^2) * rnorm(n)
    (T - sn$mu_lambda) / sn$c_lambda
}

fit_skew <- function(Z, W_j, sn, maxit = 100, tolerance = 1e-8) {
    ## W_j contains the scalar responses W_ij for a fixed mediator j
    n <- length(W_j)
    initial <- lm.fit(Z, W_j)
    scale_j <- sqrt(sum(initial$residuals^2) / n)
    if (sn$lambda == 0) {
        return(list(xi = initial$coefficients, scale2 = scale_j^2))
    }
    scaled_precision_j <- sn$c_lambda / scale_j
    par <- c(scaled_precision_j * initial$coefficients, scaled_precision_j)
    optimizer_design <- cbind(-Z, W_j)
    last <- length(par)
    nll <- function(par) {
        if (par[last] <= 0) return(Inf)
        T_ij <- sn$mu_lambda + drop(optimizer_design %*% par)
        -n * log(par[last]) + sum(T_ij^2 / 2 -
                                  pnorm(sn$lambda * T_ij, log.p = TRUE))
    }
    converged <- FALSE
    for (iter in seq_len(maxit)) {
        T_ij <- sn$mu_lambda + drop(optimizer_design %*% par)
        lambda_T <- sn$lambda * T_ij
        r_lambda <- exp(dnorm(lambda_T, log = TRUE) -
                        pnorm(lambda_T, log.p = TRUE))
        v_lambda <- T_ij - sn$lambda * r_lambda
        K_lambda <- 1 + sn$lambda^2 * r_lambda * (lambda_T + r_lambda)
        gradient <- drop(crossprod(optimizer_design, v_lambda))
        gradient[last] <- gradient[last] - n / par[last]
        H <- crossprod(optimizer_design, K_lambda * optimizer_design)
        H[last, last] <- H[last, last] + n / par[last]^2
        step <- solve(H, gradient)
        decrement <- sum(gradient * step)
        if (decrement / n < tolerance^2) {
            converged <- TRUE
            break
        }
        old <- nll(par)
        accepted <- FALSE
        for (backtrack in 0:40) {
            amount <- 2^(-backtrack)
            proposal <- par - amount * step
            value <- nll(proposal)
            if (is.finite(value) &&
                value <= old - 1e-4 * amount * decrement + 1e-12) {
                accepted <- TRUE
                break
            }
        }
        if (!accepted) stop("Line search failed")
        par <- proposal
    }
    if (!converged) stop("Optimization algorithm did not converge")
    xi <- par[-last] / par[last]
    names(xi) <- colnames(Z)
    list(xi = xi, scale2 = (sn$c_lambda / par[last])^2)
}


AR1cor <- function(q, rho_AR) {
    rho_AR^abs(outer(seq_len(q), seq_len(q), "-"))
}

simu_fun <- function(theta, L, sn, n, prob_A = 0.1) {
    q <- length(theta$nu)
    X <- cbind(x1 = rbinom(n, 1, 0.5), x2 = rexp(n, rate = 1))
    n1 <- round(n * prob_A)
    A <- sample(rep(c(0, 1), c(n - n1, n1)))
    U <- matrix(rskew(n * q, sn), n, q)
    zeta_M <- sweep(U, 2, exp(theta$eta), "*") %*% t(L)
    M <- outer(rep(1, n), theta$mu) + outer(A, theta$nu) +
        X %*% t(theta$Xi) + zeta_M
    Y <- theta$alpha + theta$gamma * A + drop(M %*% theta$beta) +
        drop(X %*% theta$delta) + rnorm(n, sd = sqrt(theta$sigma2))
    colnames(M) <- paste0("m", seq_len(q))
    data.frame(Y = Y, A = A, X, M)
}

fit_model <- function(data, L, sn) {
    n <- nrow(data)
    q <- nrow(L)
    Z <- model.matrix(~ A + x1 + x2, data)
    M <- as.matrix(data[, paste0("m", seq_len(q)), drop = FALSE])
    ## Rows of W are W_i'; column j of xi_hat contains xi_j estimates
    W <- M %*% solve(t(L))
    xi_hat <- matrix(NA_real_, ncol(Z), q)
    eta <- numeric(q)
    for (j in seq_len(q)) {
        fit <- fit_skew(Z, W[, j], sn)
        xi_hat[, j] <- fit$xi
        eta[j] <- log(fit$scale2) / 2
    }
    nu <- drop(L %*% xi_hat[2, ])
    Sigma_M <- tcrossprod(sweep(L, 2, exp(eta), "*"))
    G <- solve(crossprod(Z))
    g <- drop(Z %*% G[, 2])
    h <- rowSums((Z %*% G) * Z)
    B_A <- sn$C_B * sum(g * h)
    T_A <- sn$C_3 * sum(g^3)
    V_nu <- G[2, 2] * Sigma_M / sn$J
    ## Gaussian outcome regression with all q mediators
    outcome_design <- cbind(Z, M)
    fit <- lm.fit(outcome_design, data$Y)
    beta <- unname(fit$coefficients[ncol(Z) + seq_len(q)])
    sigma2 <- sum(fit$residuals^2) / n
    list(nu = nu, beta = beta, eta = eta, sigma2 = sigma2,
         Sigma_M = Sigma_M, V_nu = V_nu,
         B_A = B_A, T_A = T_A, L = L, n = n)
}

focus_components <- function(object) {
    nu <- object$nu
    beta <- object$beta
    V_nu <- object$V_nu
    V_beta <- object$sigma2 * solve(object$Sigma_M) / object$n

    rho <- exp(object$eta) * drop(crossprod(object$L, beta))
    b_psi <- object$B_A * sum(rho)
    kappa2 <- drop(crossprod(beta, V_nu %*% beta) +
                   crossprod(nu, V_beta %*% nu))
    kappa3 <- object$T_A * sum(rho^3) +
        6 * drop(crossprod(beta, V_nu %*% V_beta %*% nu))
    list(b_psi = b_psi, kappa2 = kappa2, kappa3 = kappa3)
}

mediation <- function(data, L, sn, what = "all") {
    what <- match.arg(what, c("all", "ML", "meanBR", "medianBR"))
    fit <- fit_model(data, L, sn)
    comp <- focus_components(fit)
    psi_hat <- sum(fit$nu * fit$beta)
    psi_meanBR <- psi_hat - comp$b_psi
    psi_tilde <- psi_meanBR + comp$kappa3 / (6 * comp$kappa2)
    out <- c(ML = psi_hat, meanBR = psi_meanBR, medianBR = psi_tilde)
    if (what != "all") return(unname(out[what]))
    attr(out, "se") <- sqrt(comp$kappa2)
    out
}

hulc <- function(data, L, sn, level = 0.95) {
    alpha_level <- 1 - level
    B <- ceiling(1 - log2(alpha_level))
    prob <- (alpha_level - 2^(1 - B)) / (2^(2 - B) - 2^(1 - B))
    B <- B - as.integer(runif(1) < prob)
    counts <- tabulate(data$A + 1, nbins = 2)
    groups <- integer(nrow(data))
    for (arm in 0:1) {
        inds <- which(data$A == arm)
        groups[inds] <- sample(rep(seq_len(B), length.out = length(inds)))
    }
    estimates <- vapply(seq_len(B), function(i) {
        mediation(data[groups == i, ], L, sn)
    }, c(ML = 0, meanBR = 0, medianBR = 0))
    out <- t(apply(estimates, 1, range))
    colnames(out) <- c("lower", "upper")
    out
}

run_experiment <- function(theta, L, sn, n = 256, R = 5000,
                           prob_A = 0.1, level = 0.95) {
    q <- length(theta$nu)
    B <- ceiling(1 - log2(1 - level))
    counts <- c(round(n * prob_A), n - round(n * prob_A))
    if (min(counts) < B || sum(counts %/% B) <= q + 4)
        stop("Increase n or decrease q")
    out <- future_lapply(seq_len(R), function(i) {
        data <- simu_fun(theta, L, sn, n, prob_A)
        estimates <- mediation(data, L, sn)
        margin <- qnorm((1 + level) / 2) * attr(estimates, "se")
        ci_wald <- cbind(lower = estimates - margin, upper = estimates + margin)
        ci_hulc <- hulc(data, L, sn, level)
        if (i %% 1000 == 0) cat(i, "/", R, "\n", sep = "")
        list(estimates = estimates, ci_wald = ci_wald, ci_hulc = ci_hulc)
    }, future.seed = TRUE)
    attr(out, "psi") <- sum(theta$nu * theta$beta)
    attr(out, "n") <- n
    attr(out, "q") <- q
    attr(out, "R") <- R
    out
}

summarize_results <- function(object) {
    psi <- attr(object, "psi")
    estimates <- vapply(object, function(x) x$estimates,
                        c(ML = 0, meanBR = 0, medianBR = 0))
    cover <- function(name) {
        rowMeans(vapply(object, function(x) {
            ci <- x[[name]]
            ci[, "lower"] <= psi & psi <= ci[, "upper"]
        }, c(ML = FALSE, meanBR = FALSE, medianBR = FALSE)))
    }
    data.frame(n = attr(object, "n"), q = attr(object, "q"),
               R = attr(object, "R"), psi = psi,
               method = rownames(estimates),
               bias = rowMeans(estimates - psi),
               pu = rowMeans(estimates < psi),
               rmse = sqrt(rowMeans((estimates - psi)^2)),
               cover_wald = cover("ci_wald"),
               cover_hulc = cover("ci_hulc"), row.names = NULL)
}

### Settings
p <- 2; rho_AR <- 0.5; lambda <- 8
prob_A <- switch(exp_set,
                 "a" = 1 / 8,
                 "b" = 3 / 32,
                 "c" = 1 / 16)
level <- 0.95
sn <- skew_constants(lambda)

settings <- expand.grid(
  q = c(10, 20, 30), n = 2^(8:10), tau = c(1, 4), R = 100000
)

summaries <- vector("list", nrow(settings))
set.seed(111)
for (s in seq_len(nrow(settings))) {
  q <- settings$q[s]
  n <- settings$n[s]
  R <- settings$R[s]
  tau <- settings$tau[s]
  L <- t(chol(AR1cor(q, rho_AR)))

  stopifnot(q %% 10 == 0)
  r <- q / 2
  k <- q / 10

  theta <- list(
    mu = rep(0, q),
    nu = c(rep(0.15, r), rep(0, k),
           rep(0, 2 * k), rep(0.15, 2 * k)),
    beta = c(rep(0.125 / r, r), rep(0, k),
             rep(0.125 / r, 2 * k), rep(0, 2 * k)),
    Xi = cbind(rep(0.4, q), 0.3 * (-1)^(seq_len(q) + 1)),
    alpha = 0, gamma = 0.5, delta = c(0.5, -0.5),
    eta = rep(0, q), sigma2 = 0.25
  )

  ## The joint indirect effect is psi = 0.01875 * tau.
  theta$beta <- tau * theta$beta
  cat("Setting: n =", n, "q =", q, "tau =", tau, "R =", R, "\n")
  res <- run_experiment(theta, L, sn, n = n, R = R,
                        prob_A = prob_A, level = level)
  summaries[[s]] <- cbind(tau = tau, summarize_results(res))
}

summaries <- do.call("rbind", summaries)
rownames(summaries) <- NULL

## Save results
results_file <- file.path(base_dir, paste0("results/multiple-mediator-", exp_set, ".rda"))
save(settings, summaries, p, rho_AR, lambda, prob_A, level, file = results_file)


