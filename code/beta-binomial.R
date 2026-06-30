base_dir <- "."
results_file <- file.path(base_dir, "results/beta-binomial.rda")

## pak::pak("eulogepagui/brbetabinomial")
pkgs <- c("brbetabinomial", "focuson")
for (pkg in pkgs) library(pkg, character.only = TRUE)

bb_one_obs_derivatives <- function(y, m, mu, phi) {
  E <- function(j) (1 - phi) * mu + j * phi
  F <- function(j) (1 - mu) * (1 - phi) + j * phi
  G <- function(j) (1 - phi) + j * phi
  jy <- if (y > 0) 0:(y - 1) else integer(0)
  jf <- if ((m - y) > 0) 0:(m - y - 1) else integer(0)
  jg <- if (m > 0) 0:(m - 1) else integer(0)
  U_mu <- sum((1 - phi) / E(jy)) - sum((1 - phi) / F(jf))
  U_phi <- sum((jy - mu) / E(jy)) + sum((jf + mu - 1) / F(jf)) - sum((jg - 1) / G(jg))
  U_mumu <- - (1 - phi)^2 * (sum(1 / E(jy)^2) + sum(1 / F(jf)^2))
  U_muphi <- -sum(jy / E(jy)^2) + sum(jf / F(jf)^2)
  U_phiphi <- -sum((mu - jy)^2 / E(jy)^2) -sum((mu + jf - 1)^2 / F(jf)^2) + sum((jg - 1)^2 / G(jg)^2)
  c(U_mu = U_mu, U_phi = U_phi, U_mumu = U_mumu, U_muphi = U_muphi, U_phiphi = U_phiphi)
}

bb_expectations_one_obs <- function(m, mu, phi) {
  y <- 0:m
  kappa <- (1 - phi) / phi
  a <- mu * kappa
  b <- (1 - mu) * kappa
  pr <- exp(lchoose(m, y) + lbeta(y + a, m - y + b) - lbeta(a, b))
  pr <- pr / sum(pr)
  D <- t(vapply(y, bb_one_obs_derivatives, numeric(5), m = m, mu = mu, phi = phi))
  U_mu <- D[, "U_mu"]
  U_phi <- D[, "U_phi"]
  U_mumu <- D[, "U_mumu"]
  U_muphi <- D[, "U_muphi"]
  U_phiphi <- D[, "U_phiphi"]
  c(L1  = sum(pr * U_mumu),
    L2  = sum(pr * U_muphi),
    L3  = sum(pr * U_phiphi),
    L4  = sum(pr * U_mu^3),
    L5  = sum(pr * U_mu^2 * U_phi),
    L6  = sum(pr * U_mu * U_phi^2),
    L7  = sum(pr * U_phi^3),
    L8  = sum(pr * U_mumu * U_mu),
    L9  = sum(pr * U_mu^2),
    L10 = sum(pr * U_mu * U_muphi),
    L11 = sum(pr * U_mu * U_phiphi),
    L12 = sum(pr * U_mu * U_phi),
    L13 = sum(pr * U_phi * U_mumu),
    L14 = sum(pr * U_phi * U_muphi),
    L15 = sum(pr * U_phi * U_phiphi),
    L16 = sum(pr * U_phi^2))
}

components_bb <- function(theta, data, Z = NULL) {
    y <- data$y
    m <- data$n
    X <- data$X
    p <- ncol(X)
    if (is.null(Z)) {
        Z <- matrix(1, nrow = length(y), ncol = 1)
    }
    q <- ncol(Z)
    beta <- theta[seq_len(p)]
    gamma <- theta[p + seq_len(q)]
    eta <- drop(X %*% beta)
    zeta <- drop(Z %*% gamma)
    mu <- plogis(eta)
    phi <- plogis(zeta)
    d1_mu <- mu * (1 - mu)
    d2_mu <- d1_mu * (1 - 2 * mu)
    d1_phi <- phi * (1 - phi)
    d2_phi <- d1_phi * (1 - 2 * phi)
    L <- t(vapply(seq_along(y), function(i) {
        bb_expectations_one_obs(m = m[i], mu = mu[i], phi = phi[i])
    }, numeric(16)))
    w_cross <- function(A, w, B) {
        crossprod(A, B * as.numeric(w))
    }
    i_bb <- -w_cross(X, L[, "L1"] * d1_mu^2, X)
    i_bg <- -w_cross(X, L[, "L2"] * d1_mu * d1_phi, Z)
    i_gg <- -w_cross(Z, L[, "L3"] * d1_phi^2, Z)
    info <- rbind(
        cbind(i_bb, i_bg),
        cbind(t(i_bg), i_gg)
    )
    V <- solve(info)
    dimnames(V) <- list(names(theta), names(theta))
    P <- Q <- vector("list", p + q)
    make_block <- function(BB, BG, GG) {
        out <- rbind(
            cbind(BB, BG),
            cbind(t(BG), GG)
        )
        dimnames(out) <- list(names(theta), names(theta))
        out
    }
    for (s in seq_len(p)) {
        xs <- X[, s]
        P_bb <- w_cross(X, L[, "L4"] * d1_mu^3 * xs, X)
        P_bg <- w_cross(X, L[, "L5"] * d1_mu^2 * d1_phi * xs, Z)
        P_gg <- w_cross(Z, L[, "L6"] * d1_mu * d1_phi^2 * xs, Z)
        Q_bb <- w_cross(X, (L[, "L8"] * d1_mu^3 + L[, "L9"] * d1_mu * d2_mu) * xs, X)
        Q_bg <- w_cross(X, L[, "L10"] * d1_mu^2 * d1_phi * xs, Z)
        Q_gg <- w_cross(Z, (L[, "L11"] * d1_mu * d1_phi^2 +
                            L[, "L12"] * d1_mu * d2_phi) * xs, Z)
        P[[s]] <- make_block(P_bb, P_bg, P_gg)
        Q[[s]] <- make_block(Q_bb, Q_bg, Q_gg)
    }
    for (t in seq_len(q)) {
        zt <- Z[, t]
        k <- p + t
        P_bb <- w_cross(X, L[, "L5"] * d1_mu^2 * d1_phi * zt, X)
        P_bg <- w_cross(X, L[, "L6"] * d1_mu * d1_phi^2 * zt, Z)
        P_gg <- w_cross(Z, L[, "L7"] * d1_phi^3 * zt, Z)
        Q_bb <- w_cross(X, (L[, "L13"] * d1_mu^2 * d1_phi + L[, "L12"] * d1_phi * d2_mu) * zt, X)
        Q_bg <- w_cross(X, L[, "L14"] * d1_mu * d1_phi^2 * zt, Z)
        Q_gg <- w_cross(Z, (L[, "L15"] * d1_phi^3 + L[, "L16"] * d1_phi * d2_phi) * zt, Z)
        P[[k]] <- make_block(P_bb, P_bg, P_gg)
        Q[[k]] <- make_block(Q_bb, Q_bg, Q_gg)
    }
    names(P) <- names(Q) <- names(theta)
    list(V = V, P = P, Q = Q, info = info)
}


data("carrots", package = "robustbase")
carrots14 <- carrots[-14, ]
start_fit <- glm(success / total ~ logdose + block, weights = total,
                 family = binomial(), data = carrots14)
brbb_settings <- expand.grid(focus = c("logit", "identity"),
                        correction = c("no", "median"))

## Using brbetabinomial
brbb_results <- apply(brbb_settings, 1, function(set) {
    brbetabinomial(success ~ logdose + block | 1, weights = total,
                   data = carrots14,
                   link.mu = "logit",
                   link.phi = set[["focus"]],
                   type = ifelse(set[["correction"]] == "median", "AS_median", "AS_ml"),
                   start = c(coef(start_fit), 0),
                   maxit = 1000,
                   epsilon = 1e-10)
})
names(brbb_results) <- with(brbb_settings, paste(correction, focus, sep = ":"))
brbb_settings$estimate <- sapply(brbb_results, function(x) coef(x)$precision)


## Using focuson
focuson_settings <- expand.grid(focus = c("logit", "identity"),
                                correction = c("no", "median"),
                                base = c("ML", "meanBR"))

## Compute P, Q, V
carrots14_list <- list(y = carrots14$success, n = carrots14$total,
                       X = model.matrix(~ logdose + block, data = carrots14))
coefs_ml <- brbb_results[["no:logit"]] |>
    coef() |> unlist()
coefs_br <- update(brbb_results[["no:logit"]], type = "AS_mean", link.phi = "logit") |>
    coef() |> unlist()
comps_ml <- components_bb(coefs_ml, carrots14_list)
comps_br <- components_bb(coefs_br, carrots14_list)

focuson_results <- apply(focuson_settings, 1, function(set) {
    if (set[["focus"]] == "identity")
        on_fun <- function(theta) plogis(theta[length(theta)])
    else
        on_fun <- function(theta) theta[length(theta)]
    coefs <- if (set[["base"]] == "ML") coefs_ml else coefs_br
    comps <- if (set[["base"]] == "ML") comps_ml else comps_br
    focus_engine(theta = coefs,
                 components = comps,
                 on = on_fun,
                 correction = set[["correction"]],
                 estimator = set[["base"]])
})
names(focuson_results) <- with(focuson_settings, paste(correction, focus, base, sep = ":"))
focuson_settings$estimate <- sapply(focuson_results, coef)

save(brbb_settings, focuson_settings, file = results_file)
