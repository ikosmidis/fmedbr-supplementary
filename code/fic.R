base_dir <- "."
results_file <- file.path(base_dir, "results/fic.rda")

pkgs <- c("brglm2", "MuMIn", "numDeriv", "focuson", "enrichwith",
          "detectseparation", "dplyr", "progressr", "tidyr", "dplyr")
for (pkg in pkgs) library(pkg, character.only = TRUE)

get_fic <- function(wide, focus, focus_grad, candidates, narrow = NULL,
                    fit_candidates = 1, control = glm.control(), ...) {
    is_ML <- isTRUE(is.null(control$type))
    if (!is_ML) stopifnot(inherits(wide, brglm2::brglm_fit))
    candidates <- as.matrix(candidates) == 1
    if (is.null(narrow)) {
        inds_beta <- colMeans(candidates) == 1
    } else {
        inds_beta <- as.logical(unlist(narrow))
    }
    inds_gamma <- !inds_beta
    gammahat <- coef(wide)[inds_gamma]
    if (missing(focus_grad))
        dpsi <- numDeriv::grad(focus, coef(wide), ...)
    else
        dpsi <- focus_grad(coef(wide), ...)
    dpsi_dt <- dpsi[inds_beta]
    dpsi_dg <- dpsi[inds_gamma]
    V <- vcov(wide)
    i <- solve(V)
    i_bb <- i[inds_beta, inds_beta, drop = FALSE]
    i_gb <- i[inds_gamma, inds_beta, drop = FALSE]
    i_bb_inv <- solve(i_bb)
    V_gg <- V[inds_gamma, inds_gamma, drop = FALSE]
    V_gg_inv <- solve(V_gg)
    w <- drop(i_gb %*% i_bb_inv %*% dpsi_dt - dpsi_dg)
    t0sq <- drop(t(dpsi_dt) %*% i_bb_inv %*% dpsi_dt)
    A <- tcrossprod(gammahat) - V_gg
    q <- sum(inds_gamma)
    n_models <- nrow(candidates)
    res <- matrix(NA_real_, n_models, 8)
    colnames(res) <- c("r_V", "B", "r_B_sq", "r_B_sq_adj",
                       "npar", "focus", "AIC", "BIC")
    gamma_sets <- lapply(seq_len(n_models), function(i) {
        which(candidates[i, inds_gamma])
    })
    npars <- rowSums(candidates)
    for (i in seq_len(n_models)) {
        s <- gamma_sets[[i]]
        if (length(s) > 0) {
            V_gg_s <- solve(V_gg_inv[s, s, drop = FALSE])
            L_s <- V_gg_s %*% V_gg_inv[s, , drop = FALSE]
            Gtw <- drop(crossprod(L_s, w[s]))
            u <- w - Gtw
            B <- sum(u * gammahat)
            biassq <- B^2 - drop(crossprod(u, V_gg %*% u))
            Vs <- drop(t0sq + crossprod(w[s], V_gg_s %*% w[s]))
            res[i, "B"] <- B
        } else {
            B <- sum(w * gammahat)
            biassq <- B^2 - drop(crossprod(w, V_gg %*% w))
            Vs <- t0sq
            res[i, "B"] <- B
        }
        sgn <- sign(res[i, "B"])
        res[i, "r_B_sq"] <- sgn * ifelse(biassq < 0, NA_real_, sqrt(biassq))
        res[i, "r_B_sq_adj"] <- sgn * sqrt(max(0, biassq))
        res[i, "r_V"] <- sqrt(Vs)
        res[i, "npar"] <- npars[i]
    }
    res <- res |>
        as.data.frame() |>
        transform(rmse = sqrt(r_B_sq^2 + r_V^2),
                  rmse_adj = sqrt(r_B_sq_adj^2 + r_V^2))
    rownames(res) <- rownames(candidates)
    res$model <- apply(candidates, 1, function(x) paste0(as.numeric(x), collapse = ""))
    res <- res[order(res$rmse_adj), ]
    if (isTRUE(fit_candidates > 0)) {
        mf_wide <- model.frame(wide)
        y <- model.response(mf_wide)
        w <- model.weights(mf_wide)
        off <- model.offset(mf_wide)
        X_wide <- model.matrix(wide)
        zero_coefs <- coef(wide)
        zero_coefs[] <- 0
        for (i in rownames(res)[1:fit_candidates]) {
            wh <- candidates[i, ]
            XX <- X_wide[, wh, drop = FALSE]
            if (is_ML) {
                m <- glm(y ~ -1 + XX, weights = w, offset = off,
                         family = wide$family, control = control)
            } else {
                m <- glm(y ~ -1 + XX, weights = w, offset = off,
                         family = wide$family, method = brglm2::brglm_fit,
                         control = control)
            }
            c_coefs <- coef(m)
            c_coefs_wide <- zero_coefs
            c_coefs_wide[gsub("XX", "", names(c_coefs))] <- c_coefs
            res[i, "AIC"] <- AIC(m)
            res[i, "BIC"] <- BIC(m)
            res[i, "focus"] <- focus(c_coefs_wide, ...)
        }
    }
    res[, c("model", "focus", "r_V", "B", "rmse_adj",
            "npar", "AIC", "BIC",
            "r_B_sq_adj", "r_B_sq", "rmse")]
}

probdiff <- function(theta, x_setting1, x_setting2) {
    etas <- c(x_setting1 %*% theta,
              x_setting2 %*% theta)
    ps <- plogis(etas) #1 / (1 + exp(-etas))
    ps[2] - ps[1]
}

probdiff_grad <- function(theta, x_setting1, x_setting2) {
    eta1 <- drop(x_setting1 %*% theta)
    eta2 <- drop(x_setting2 %*% theta)
    p1 <- plogis(eta1)
    p2 <- plogis(eta2)
    out <- p2 * (1 - p2) * x_setting2 - p1 * (1 - p1) * x_setting1
    out
}


make_candidate_set <- function(m_wide, fixed = c("x1", "x2")) {
    all_calls <- dredge(m_wide, evaluate = FALSE, fixed = fixed)
    wide_vars <- colnames(model.matrix(m_wide))
    S_vars <- sapply(all_calls, function(x) {
        as.numeric(wide_vars %in% colnames(model.matrix(formula(x), data = eval(x$data))))
    }) |> t()
    colnames(S_vars) <- wide_vars
    as.data.frame(S_vars)
}

simulate_df <- function(object, theta) {
    dat <- model.frame(bwt_mod)
    simu_fun <- get_simulate_function(bwt_mod)
    dat[[as.character(formula(bwt_mod)[[2]])]] <- simu_fun(theta)[[1]]
    data.frame(dat)
}

get_model <- function(str) {
    strsplit(str, "")[[1]] == "1"
}



data("birthwt", package = "MASS")
## As in ?MASS::birthwt
bwt <- with(birthwt, {
    race <- factor(race, labels = c("white", "black", "other"))
    ptd <- ptl > 0
    ftv <- factor(ftv)
    levels(ftv)[-(1:2)] <- "2+"
    data.frame(low = factor(low), age, lwt, race, smoke = (smoke > 0),
               ptd, ht = (ht > 0), ui = (ui > 0), ftv)
})

fm <- low ~ age + lwt + race + smoke + ptd + ht + ui + ftv
bwt_mod <- glm(fm, family = binomial("logit"), data = bwt, na.action = na.fail)
## Protected covariate is `lwt`
prot <- c("race", "smoke")
theta <- coef(bwt_mod)
## Set of candidate models
S <- make_candidate_set(bwt_mod, fixed = prot)

## Get representative settings per smoking status
x_settings_df <- bwt |>
    summarize(age = mean(age),
              lwt = mean(lwt),
              ptd = FALSE,
              ht = FALSE,
              ui = FALSE,
              ftv = factor("1", levels(ftv)),
              smoke = FALSE,
              .by = c(race))
x_settings_non_sm <- x_settings_df |> mutate(smoke = FALSE) |> arrange(race)
x_settings_sm <- x_settings_df |> mutate(smoke = TRUE) |> arrange(race)
x_settings1 <- model.matrix( ~ ., data = x_settings_non_sm)[, names(theta)]
x_settings2 <- model.matrix( ~ ., data = x_settings_sm)[, names(theta)]


### Simulation results
contr <- brglm_control(maxit = 1000, check_aliasing = FALSE, slowit = 0.5)
X <- model.matrix(bwt_mod)
set.seed(111)
R <- 10000
res <- vector("list", R)
with_progress({
    pr <- progressor(R)
    for (r in 1:R) {
        pr()
        cdat <- simulate_df(bwt_mod, theta)
        mod_ml <- glm(formula = fm,
                      family = binomial("logit"),
                      data = cdat,
                      na.action = na.fail)
        separation <- update(mod_ml, method = detectseparation::detect_separation)$outcome
        mod_mean_br <- glm(formula = fm,
                           family = binomial("logit"),
                           data = cdat,
                           na.action = na.fail,
                           method = brglm2::brglm_fit,
                           type = "AS_mean",
                           control = contr)
        est_ml <- est_median_br <- est_fic_ml <- est_fic_ml_br <- rep(NA, 3)
        ci_ml <- ci_median_br <- ci_hulc <- ci_fic_ml <- ci_fic_ml_br <- matrix(NA, 3, 2)
        for (j in 1:nrow(x_settings1)) {
            ## Key objects
            obj_ml <- focus(mod_ml,
                            on = probdiff, on_gradient = probdiff_grad,
                            correction = "no",
                            x_setting1 = x_settings1[j, ],
                            x_setting2 = x_settings2[j, ])
            est_ml[j] <- coef(obj_ml)
            ci_ml[j, ] <- confint(obj_ml)
            obj_median_br <- focus(mod_mean_br,
                                   on = probdiff, on_gradient = probdiff_grad,
                                   correction = "median",
                                   x_setting1 = x_settings1[j, ],
                                   x_setting2 = x_settings2[j, ])
            fics <- get_fic(mod_ml,
                            focus = probdiff,
                            focus_grad = probdiff_grad,
                            candidates = S,
                            fit_candidates = FALSE,
                            x_setting1 = x_settings1[j, ],
                            x_setting2 = x_settings2[j, ])
            ## Wide model statistics
            ## est_ml[j] <- coef(obj_ml)
            est_median_br[j] <- coef(obj_median_br)
            ## ci_ml[j, ] <- confint(obj_ml)
            ci_median_br[j, ] <- confint(obj_median_br)
            for (i in 1:10) {
                cih <- try(confint(obj_median_br, method = "hulc"), silent = TRUE)
                if (inherits(cih, "try-error")) {
                    print("oops")
                } else {
                    if (all(!is.na(cih)))
                        break
                }
            }
            ci_hulc[j, ] <- cih
            ## Selected model statistics
            inds <- get_model(fics[1, "model"])
            best_ml <- glm(formula = low ~ -1 + X[, inds],
                           family = binomial("logit"),
                           data = cdat, na.action = na.fail)
            best_fit <- focus(best_ml,
                              on = probdiff, on_gradient = probdiff_grad,
                              correction = "no",
                              x_setting1 = x_settings1[j, inds],
                              x_setting2 = x_settings2[j, inds])
            est_fic_ml[j] <- coef(best_fit)
            est_fic_ml_br[j] <- est_fic_ml[j] - fics[1, "B"]
            ci_fic_ml[j, ] <- confint(best_fit)
            ci_fic_ml_br[j, ] <- est_fic_ml_br[j] + c(-1, 1) * qnorm(1 - 0.05/2) * obj_ml$se
        }
        res[[r]] <- list(est_ml = est_ml,
                         est_median_br = est_median_br,
                         est_fic_ml = est_fic_ml,
                         est_fic_ml_br = est_fic_ml_br,
                         ci_ml = ci_ml,
                         ci_median_br = ci_median_br,
                         ci_hulc = ci_hulc,
                         ci_fic_ml = ci_fic_ml,
                         ci_fic_ml_br = ci_fic_ml_br)
    }
})

## Summaries
truth <- predict(bwt_mod, type = "response", newdata = x_settings_sm) -
    predict(bwt_mod, type = "response", newdata = x_settings_non_sm)
bias <- rbind(ml = rowMeans(sapply(res, "[[", "est_ml") - truth),
              median_br = rowMeans(sapply(res, "[[", "est_median_br") - truth),
              fic_ml = rowMeans(sapply(res, "[[", "est_fic_ml") - truth),
              fic_ml_br = rowMeans(sapply(res, "[[", "est_fic_ml_br") - truth))
pu <- rbind(ml = rowMeans(sapply(res, "[[", "est_ml") < truth),
            median_br = rowMeans(sapply(res, "[[", "est_median_br") < truth),
            fic_ml = rowMeans(sapply(res, "[[", "est_fic_ml") < truth),
            fic_ml_br = rowMeans(sapply(res, "[[", "est_fic_ml_br") < truth))
mse <- rbind(ml = rowMeans((sapply(res, "[[", "est_ml") - truth)^2),
             median_br = rowMeans((sapply(res, "[[", "est_median_br") - truth)^2),
             fic_ml = rowMeans((sapply(res, "[[", "est_fic_ml") - truth)^2),
             fic_ml_br = rowMeans((sapply(res, "[[", "est_fic_ml_br") - truth)^2))
cover <- rbind(ml = rowMeans(sapply(res, function(x)
                   between(truth, x$ci_ml[, 1], x$ci_ml[, 2])), na.rm = TRUE),
               median_br = rowMeans(sapply(res, function(x)
                   between(truth, x$ci_median_br[, 1], x$ci_median_br[, 2])), na.rm = TRUE),
               fic_ml = rowMeans(sapply(res, function(x)
                   between(truth, x$ci_fic_ml[, 1], x$ci_fic_ml[, 2])), na.rm = TRUE),
               fic_ml_br = rowMeans(sapply(res, function(x)
                   between(truth, x$ci_fic_ml_br[, 1], x$ci_fic_ml_br[, 2])), na.rm = TRUE),
               hulc = rowMeans(sapply(res, function(x)
                   between(truth, x$ci_hulc[, 1], x$ci_hulc[, 2])), na.rm = TRUE))

nams <- levels(x_settings_sm$race)
bias <- bias |> data.frame() |>
    setNames(nams) |>
    transform(method = rownames(bias)) |>
    pivot_longer(!method, names_to = "race", values_to = "bias")
mse <- mse |> data.frame() |>
    setNames(nams) |>
    transform(method = rownames(mse)) |>
    pivot_longer(!method, names_to = "race", values_to = "mse")
pu <- pu |> data.frame() |>
    setNames(nams) |>
    transform(method = rownames(pu)) |>
    pivot_longer(!method, names_to = "race", values_to = "pu")
cover <- cover |> data.frame() |>
    setNames(nams) |>
    transform(method = rownames(cover)) |>
    pivot_longer(!method, names_to = "race", values_to = "cover")

summaries <- bias |>
    merge(mse, by = c("method", "race")) |>
    merge(pu, by = c("method", "race")) |>
    merge(cover, by = c("method", "race"), all.y = TRUE)
summaries$method <- factor(summaries$method, c("ml", "fic_ml", "fic_ml_br", "median_br", "hulc"),
                           order = TRUE)
summaries$race <- factor(summaries$race, c("white", "black", "other"), ordered = TRUE)
summaries <- summaries |> arrange(race, method)
summaries$psi <- rep(truth, each = 5)

save(summaries, file = results_file)


