pkgs <- c("ggplot2", "tinytable", "dplyr", "tidyr", "patchwork")
for (pkg in pkgs) library(pkg, character.only = TRUE)

base_dir <- "."
results_dir <- file.path(base_dir, "results/")
tables_dir <- file.path(base_dir, "tables/")
figures_dir <- file.path(base_dir, "figures/")

## Mahalanobis distance
#######################
img_file <- file.path(results_dir, "mahalanobis-distance.rda")
out_file <- file.path(tables_dir, "mahalanobis-distance.tex")
## Load results
load(img_file)
caption <- "Comparison of the ML estimator ($\\hat\\psi$) and the MBC focus estimator ($\\tilde\\psi$) of the squared Mahalanobis distance in terms of simulation-based estimates (see Example~\\ref{ex:mahalanobis-distance}) of mean bias (BIAS), mean absolute deviation (MAD), probability of underestimation (PU), root mean squared error (RMSE), and coverage of $95\\%$ Wald-type confidence intervals, where the standard error is estimated by $\\{K_2(\\hat\\psi)\\}^{1/2}$ and $\\{K_2(\\tilde\\psi)\\}^{1/2}$, respectively. All summaries are $\\times 100$."
## Create table
summaries$rmse <- sqrt(summaries$mse)
summaries$method <- factor(summaries$method, levels = c("ML", "medianBR", "meanBR"), ordered = TRUE)
keep <- c("p", "method", "n", "mean_bias", "mean_abs_bias", "pu", "rmse", "cover_wald")
new_names <- c("$p$", "Estimator", "$n$", "BIAS", "MAD", "PU", "RMSE", "Coverage")
summaries_keep <- summaries |>
    sort_by(~ p + method) |>
    subset(method %in% c("ML", "medianBR"),
           select = keep) |>
    mutate(method = recode_values(method,
                                  "ML" ~ "$\\hat\\psi$",
                                  "medianBR" ~ "$\\tilde\\psi$"))
summaries_keep |>
    transform(cover_wald = cover_wald * 100,
              mean_bias = mean_bias * 100,
              mean_abs_bias = mean_abs_bias * 100,
              pu = pu * 100,
              rmse = rmse * 100) |>
    tt(caption = caption) |>
    setNames(new_names) |>
    format_tt(j = c(1, 3), digits = 1,
              num_fmt = "decimal", num_zero = TRUE) |>
    format_tt(j = c(4:8), digits = 1,
              num_fmt = "decimal", num_zero = TRUE) |>
    style_tt(align = "rldddddd") |>
    style_tt(i = "colnames", align = "c") |>
    ## multirows
    style_tt(i = seq(1, 24, 8), j = 1, rowspan = 8, alignv = "t") |>
    style_tt(i = seq(1, 24, 4), j = 2, rowspan = 4, alignv = "t") |>
    ## Midrules
    style_tt(i = c(9, 17), line = "t", line_width = 0.05) |>
    style_tt(i = c(5, 13, 21), line = "t", line_width = 0.05) |>
    theme_latex(outer = "label={tab:mahalanobis-distance}",
                placement = "t!") |>
    save_tt(out_file, overwrite = TRUE)


## Individual marginal effects
##############################
img_file <- file.path(results_dir, "marginal-effects.rda")
out_file <- file.path(figures_dir, "marginal-effects.pdf")
load(img_file)
summaries <- summaries |>
    transform(cover_wald = cover_wald,
              cover_hulc = cover_hulc,
              bias2_o_var = mean_bias^2 / (mse - mean_bias^2),
              pu = pu,
              correction = factor(method, levels = c("median", "mean", "no"), ordered = TRUE),
              base = factor(base, levels = c("ML", "mean BR"), labels = c("ML", "meanBR"), ordered = TRUE)) |>
    transform(corbas = correction:base)
padd <- 0.1
ylims <- c(0, 0.8)
titl <- "Probability of underestimation"
fig_pu <- ggplot(summaries) +
    geom_col(aes(student, pu, fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0.50), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
ylims <- c(0, 0.025)
titl <- "Mean bias"
fig_bias <- ggplot(summaries) +
    geom_col(aes(student, mean_bias, fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
ylims <- c(0, 0.1)
titl <- "RMSE"
fig_rmse <- ggplot(summaries) +
    geom_col(aes(student, sqrt(mse), fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
fig <- ((fig_pu + fig_bias + fig_rmse + plot_layout(axes = "collect"))) + #/
        ## (fig_cover_wald + fig_cover_hulc + plot_layout(axes = "collect"))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal" ) &
    scale_fill_manual(name = "correction [base]",
                      values = hcl.colors(6),
                      labels = paste0(gsub(":", " [", levels(summaries$corbas)), "]")) &
    guides(colour = guide_legend(nrow = 1), fill   = guide_legend(nrow = 1))
pdf(out_file, width = 10, height = 3)
print(fig)
dev.off()

out_file <- file.path(figures_dir, "marginal-effects-cover.pdf")
ylims <- c(0.5, 1.0)
titl <- "Wald [S]"
fig_cover_wald <- ggplot(summaries) +
    geom_col(aes(student, cover_wald, fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0.95), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
ylims <- c(0.5, 1.0)
titl <- "Wald [C]"
fig_cover_wald_c <- ggplot(summaries) +
    geom_col(aes(student, cover_wald_comp, fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0.95), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
ylims <- c(0.5, 1.0)
titl <- "HulC"
fig_cover_hulc <- ggplot(summaries) +
    geom_col(aes(student, cover_hulc, fill = corbas), position = position_dodge2(padding = padd)) +
    geom_hline(aes(yintercept = 0.95), col = "grey") +
    labs(x = "Student", y = NULL, title = titl) +
    coord_flip(ylim = ylims) +
    theme_minimal()
fig <- (fig_cover_wald + fig_cover_wald_c + fig_cover_hulc + plot_layout(axes = "collect")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal" ) &
    scale_fill_manual(name = "correction [base]",
                      values = hcl.colors(6),
                      labels = paste0(gsub(":", " [", levels(summaries$corbas)), "]")) &
    guides(colour = guide_legend(nrow = 1), fill   = guide_legend(nrow = 1))
pdf(out_file, width = 10, height = 3)
print(fig)
dev.off()


## Beta binomial - carrots
##########################
img_file <- file.path(results_dir, "beta-binomial.rda")
out_file <- file.path(tables_dir, "beta-binomial.tex")
## Load results
load(img_file)
brbb_res <- brbb_settings |>
    subset(select = c("correction", "estimate", "focus")) |>
    pivot_wider(names_from = "correction", values_from = "estimate")
focuson_res <- focuson_settings |>
    subset(select = c("correction", "estimate", "focus", "base")) |>
    pivot_wider(names_from = c("correction", "base"), values_from = "estimate") |>
    subset(select = c("focus", "median_ML", "median_meanBR"))
res <- brbb_res |>  left_join(focuson_res, by = "focus")
levels(res$focus) <- c("$\\xi$", "$\\rho$")
caption <- "Estimates of the beta-binomial overdispersion parameter for the carrots data in Example~\\ref{ex:beta-binomial}. The table reports estimates on both the logit scale $\\xi=\\log\\{\\rho/(1 - \\rho)\\}$ and the intra-class correlation scale $\\rho$. The columns show the ML estimator $\\hat\\psi$, the reduced median-bias estimator $\\psi^*$ of \\citet{pagui+etal:2017}, and the MBC estimator $\\tilde\\psi$ in~\\eqref{eq:focused}, computed using either the ML estimator (ML) or the reduced mean-bias estimator (meanBR) of $\\btheta = (\\beta_1, \\ldots, \\beta_4, \\xi)^\\top$."
res |>
    rename(`$\\psi^*$` = median,
           `$\\hat\\psi$` = no,
           `ML` = median_ML,
           `meanBR` = median_meanBR,
           `$\\psi$` = focus) |>
    tt(caption = caption) |>
    group_tt(j = list(" " = 1:3, "$\\tilde\\psi$" = 4:5)) |>
    format_tt(j = 2:5, digits = 5, num_fmt = "decimal", num_zero = TRUE) |>
    style_tt(align = "lrrrr") |>
    style_tt(i = "colnames", align = "c") |>
    style_tt(i = -1, j = 4:5, align = "c") |>
    theme_latex(outer = "label={tab:beta-binomial}",
                placement = "t!") |>
    save_tt(out_file, overwrite = TRUE)



## Weibull quantiles
####################
img_file <- file.path(results_dir, "weibull.rda")
aqs <- c(0.01, 0.05, 0.10)
out_files <- file.path(tables_dir, paste0("weibull", aqs * 100, ".tex"))

## Load results
load(img_file)
summaries$rmse <- sqrt(summaries$mse)
keep <- c("n", "estimator", "pu", "rmse",
          "cover_rstar", "cover_hulc", "cover_wald_supp", "cover_wald_comp")
new_names <- c("$n$", "Estimator", "PU", "RMSE",
               "$r^*$", "HulC", "Wald [S]", "Wald [C]")
for (r in seq_along(aqs)) {
    caption <- paste0("Comparison of the ML estimator ($\\hat\\psi$), the MBC focus estimator ($\\tilde\\psi$), and the $r^*$-based estimator ($\\psi^*$) of the $", sprintf("%0.2f", 1 - aqs[r]), "$ quantile in terms of simulation-based estimates (see Example~\\ref{ex:weibull}) of probability of underestimation (PU), root mean squared error (RMSE), and coverage of nominally $95\\%$ $r^*$-based confidence intervals, HulC-type confidence intervals with $\\Delta = 0$, and Wald-type confidence intervals, based on supplied and compatible estimates of $\\btheta$ (``Wald [S]'' and ``Wald [C]'', respectively). All summaries are $\\times 100$.")
    summaries_keep <- summaries |>
        subset(estimator != "mean" & alpha_q == aqs[r],
               select = keep) |>
        transform(estimator = recode_values(estimator,
                                            "no" ~ "$\\hat\\psi$",
                                            "median" ~ "$\\tilde\\psi$",
                                            "rs" ~ "$\\psi^*$")) |>
        transform(pu = pu * 100,
                  rmse = rmse * 100,
                  cover_wald_supp = cover_wald_supp * 100,
                  cover_wald_comp = cover_wald_comp * 100,
                  cover_hulc = cover_hulc * 100,
                  cover_rstar = cover_rstar * 100) |>
        tt(caption = caption) |>
        setNames(new_names) |>
        format_tt(j = c(3:8), digits = 1, num_fmt = "decimal", num_zero = TRUE) |>
        format_tt(replace = TRUE) |>
        style_tt(align = "rlrrrrrr") |>
        style_tt(i = "colnames", align = "c") |>
        ## multirows
        style_tt(i = seq(1, 12, 3), j = 1, rowspan = 3, alignv = "t") |>
        ## Midrules
        style_tt(i = c(4, 7, 10), line = "t", line_width = 0.05) |>
        theme_latex(outer = paste0("label={tab:weibull", aqs[r] * 100, "}"),
                    placement = "t!") |>
        save_tt(out_files[r], overwrite = TRUE)
}


## Ordinal superiority measures
####################
img_file <- file.path(results_dir, "osm.rda")
out_file <- file.path(tables_dir, "osm.tex")
## Load results
load(img_file)
caption <- "Comparison of the plug-in estimator $\\gamma^\\dagger = h(\\btheta^\\dagger)$ of~(\\ref{eq:osm}) with the MBC focus estimator $\\tilde\\gamma$ based on Monte Carlo estimates $\\bP_1, \\ldots, \\bP_p$ at $\\btheta^\\dagger$ ($R = 500$), in terms of simulation-based estimates (see Example~\\ref{ex:ordinal-superiority}) of probability of underestimation (PU), mean bias (BIAS), root mean squared error (RMSE), and coverage of nominally $95\\%$ Wald-type confidence intervals based on the supplied estimates of $\\btheta$. All summaries are $\\times 100$."
new_names <- c("Temperature", "Estimator", "PU", "BIAS", "RMSE", "Wald")
summaries |>
    transform(pu = pu * 100,
              bias = meanBias * 100,
              rmse = rmse * 100,
              coverage = coverage * 100) |>
    subset(method != "meanBR", select = c("temperature", "method", "pu", "bias", "rmse", "coverage")) |>
    setNames(new_names) |>
    transform(Estimator = recode_values(Estimator,
                                        "plugin" ~ "$\\gamma^\\dagger$",
                                        "medianBR" ~ "$\\tilde\\gamma$")) |>
    tt(caption = caption) |>
    format_tt(j = c(3:6), digits = 1, num_fmt = "decimal", num_zero = TRUE) |>
    style_tt(align = "llrrrr") |>
    style_tt(i = "colnames", align = "c") |>
    ## multirows
    style_tt(i = c(1, 3), j = 1, rowspan = 2, alignv = "t") |>
    theme_latex(outer = "label={tab:osm}",
                placement = "t!") |>
    save_tt(out_file, overwrite = TRUE)


## FIC
####################
img_file <- file.path(results_dir, "fic.rda")
out_file <- file.path(tables_dir, "fic.tex")
caption <- paste0("Comparison of the ML focus estimator from the wide model ($\\hat\\psi_{(\\mA)}$), the ML focus estimator and its corrected version from the FIC-selected model ($\\hat\\psi_{(\\hat{\\mS})}$ and cor.~$\\hat\\psi_{(\\hat{\\mS})}$, respectively), and the MBC estimator from the wide model ($\\tilde\\psi_{(\\mA)}$) in terms of simulation-based estimates (see Example~\\ref{ex:fic}) of probability of underestimation (PU), mean bias (BIAS), and root mean squared error (RMSE). The focus parameters are the three race-specific risk differences for low birth weight due to smoking during pregnancy. The column ``Wald'' reports the coverage of nominally $95\\%$ Wald-type intervals. Those corresponding to the MBC estimator are based on supplied estimates of $\\btheta$, and those based on cor.~$\\hat\\psi_{(\\hat{\\mS})}$ use the estimated standard error from the wide model.  The column ``HulC'' reports the coverage of nominally $95\\%$ HulC-type intervals with $\\Delta = 0$. All summaries are $\\times 100$.")
## Load results
load(img_file)
summaries <- summaries |>
    transform(Wald = cover,
              HulC = cover)
summaries$Wald[summaries$method == "hulc"] <- NA
summaries$HulC[summaries$method == "median_br"] <- summaries$HulC[summaries$method == "hulc"]
summaries$HulC[summaries$method != "median_br"] <- NA
summaries <- subset(summaries, method != "hulc")
new_names <- c("Race", "Estimator", "PU", "BIAS", "RMSE", "Wald", "HulC")
summaries |>
    transform(
        bias = bias * 100,
        pu = pu * 100,
        rmse = sqrt(mse) * 100,
        Wald = Wald * 100,
        HulC = HulC * 100) |>
    select(race, method, pu, bias, rmse, Wald, HulC) |>
    mutate(method = recode_values(method,
                                  "ml" ~ "$\\hat\\psi_{(\\mA)}$",
                                  "fic_ml" ~ "$\\hat\\psi_{(\\hat{\\mS})}$",
                                  "fic_ml_br" ~ "cor.~$\\hat\\psi_{(\\hat{\\mS})}$",
                                  "median_br" ~ "$\\tilde\\psi_{(\\mA)}$")) |>
    setNames(new_names) |>
    tt(caption = caption) |>
    format_tt(j = c(3:7), digits = 1, num_fmt = "decimal", num_zero = TRUE) |>
    format_tt(replace = TRUE) |>
    style_tt(align = "llrrrrr") |>
    style_tt(i = "colnames", align = "c") |>
    ## multirows
    style_tt(i = seq(1, 12, 4), j = 1, rowspan = 4, alignv = "t") |>
    ## Midrules
    style_tt(i = c(5, 9), line = "t", line_width = 0.05) |>
    theme_latex(outer = "label={tab:fic}", placement = "t!") |>
    save_tt(out_file, overwrite = TRUE)



## Mahalanobis distance between 2 distributions
##################################################
img_file <- file.path(results_dir, "mahalanobis-distance-2sample.rda")
out_file <- file.path(tables_dir, "mahalanobis-distance-2sample.tex")
load(img_file)
caption <- "Comparison of the ML estimator ($\\hat\\psi$) and the MBC focus estimator ($\\tilde\\psi$) of the squared Mahalanobis distance between two multivariate normal distributions in terms of simulation-based estimates of mean bias (BIAS), mean absolute deviation (MAD), probability of underestimation (PU), root mean squared error (RMSE), and coverage of $95\\%$ Wald-type confidence intervals, where the standard error is estimated by $\\{K_2(\\hat\\psi)\\}^{1/2}$ and $\\{K_2(\\tilde\\psi)\\}^{1/2}$, respectively. All summaries are $\\times 100$."
summaries$rmse <- sqrt(summaries$mse)
summaries$method <- factor(summaries$method,
                           levels = c("ML", "medianBR", "meanBR"),
                           ordered = TRUE)
keep <- c("p", "method", "n", "mean_bias", "mean_abs_bias",
          "pu", "rmse", "cover_wald")
new_names <- c("$p$", "Estimator", "$n$", "BIAS", "MAD",
               "PU", "RMSE", "Coverage")
summaries_keep <- summaries |>
  sort_by(~ p + method) |>
  subset(method %in% c("ML", "medianBR"),
         select = keep) |>
  mutate(method = recode_values(method,
                                "ML" ~ "$\\hat\\psi$",
                                "medianBR" ~ "$\\tilde\\psi$"))
summaries_keep |>
  transform(cover_wald = cover_wald * 100,
            mean_bias = mean_bias * 100,
            mean_abs_bias = mean_abs_bias * 100,
            pu = pu * 100,
            rmse = rmse * 100) |>
  tt(caption = caption) |>
  setNames(new_names) |>
  format_tt(j = c(1, 3), digits = 0,
            num_fmt = "decimal", num_zero = TRUE) |>
  format_tt(j = c(4:8), digits = 1,
            num_fmt = "decimal", num_zero = TRUE) |>
  style_tt(align = "rldddddd") |>
  style_tt(i = "colnames", align = "c") |>
  style_tt(i = seq(1, 24, 8), j = 1, rowspan = 8, alignv = "t") |>
  style_tt(i = seq(1, 24, 4), j = 2, rowspan = 4, alignv = "t") |>
  style_tt(i = c(9, 17), line = "t", line_width = 0.05) |>
  style_tt(i = c(5, 13, 21), line = "t", line_width = 0.05) |>
  theme_latex(outer = "label={tab:mahalanobis-distance-two-sample}",
              placement = "t!") |>
  save_tt(out_file, overwrite = TRUE)



## Circular variance in the bivariate von Mises sine model
#########################################################
img_file <- file.path(results_dir, "bvmsin-circular-variance.rda")
out_file <- file.path(tables_dir, "bvmsin-circular-variance.tex")
load(img_file)
caption <- "Comparison of the plug-in estimator $\\hat\\psi = h(\\hat\\btheta)$ of the circular variance of the first angular component in the bivariate von Mises sine model with the mean bias-corrected focus estimator $\\hat\\psi - \\hatsabias$ and the MBC focus estimator $\\tilde\\psi$, both based on Monte Carlo estimates of $\\iinfo(\\btheta)$, $\\bP_1,\\ldots,\\bP_p$ and $\\bQ_1,\\ldots,\\bQ_p$ at $\\hat\\btheta$ ($R = 500$), in terms of simulation-based estimates of probability of underestimation (PU), mean bias (BIAS), mean absolute deviation (MAD), root mean squared error (RMSE), and coverage of nominally $95\\%$ Wald-type confidence intervals based on supplied and compatible estimates of $\\btheta$ (Wald [S] and Wald [C], respectively), and HulC-type confidence intervals. All summaries are $\\times 100$."
summaries$rmse <- sqrt(summaries$mse)
summaries$method <- factor(summaries$method,
                           levels = c("ML", "meanBR", "medianBR"),
                           ordered = TRUE)
keep <- c("n", "method", "mean_bias", "mean_abs_bias",
          "pu", "rmse", "cover_wald",
          "cover_hulc")
new_names <- c("$n$", "Estimator", "BIAS", "MAD",
               "PU", "RMSE", "Wald", "HulC")
summaries_keep <- summaries |>
  sort_by(~ n + method) |>
  subset(method %in% c("ML", "meanBR", "medianBR"),
         select = keep) |>
  mutate(method = recode_values(method,
                                "ML" ~ "$\\hat\\psi$",
                                "meanBR" ~ "$\\hat\\psi - \\hatsabias$",
                                "medianBR" ~ "$\\tilde\\psi$"))
summaries_keep |>
  transform(mean_bias = mean_bias * 100,
            mean_abs_bias = mean_abs_bias * 100,
            pu = pu * 100,
            rmse = rmse * 100,
            cover_wald = cover_wald * 100,
            cover_hulc = cover_hulc * 100) |>
  tt(caption = caption) |>
  setNames(new_names) |>
  format_tt(j = 1, digits = 0,
            num_fmt = "decimal", num_zero = TRUE) |>
  format_tt(j = c(3:8), digits = 1,
            num_fmt = "decimal", num_zero = TRUE) |>
  format_tt(replace = TRUE) |>
  style_tt(align = "rlrrrrrr") |>
  style_tt(i = "colnames", align = "c") |>
  style_tt(i = seq(1, 12, 3), j = 1, rowspan = 3, alignv = "t") |>
  style_tt(i = c(4, 7, 10), line = "t", line_width = 0.05) |>
  theme_latex(outer = "label={tab:bvmsin-circular-variance}",
              placement = "t!") |>
  save_tt(out_file, overwrite = TRUE)


## Multiple mediators
#########################################################
for (exp_set in c("a", "b", "c")) {
    img_file <- file.path(results_dir, paste0("multiple-mediator-", exp_set, ".rda"))
    load(img_file)
    methods <- c("ML", "medianBR")
    estimator_names <- c(ML = "$\\hat\\psi$", medianBR = "$\\tilde\\psi$")
    tables <- data.frame(tau = c(1, 4),
                         name = c(paste0("multiple-mediator-", exp_set, "1"),
                                  paste0("multiple-mediator-", exp_set, "4")))
    new_names <- c("$q$", "Estimator", "$n$", "PU", "BIAS", "RMSE", "Wald", "HulC")
    A <- switch(exp_set,
                 "a" = "n / 8",
                 "b" = "3n / 32",
                 "c" = "n / 16")
    for (s in seq_len(nrow(tables))) {
        tau_value <- tables$tau[s]
        out_file <- file.path(tables_dir, paste0(tables$name[s], ".tex"))
        caption <- paste0("Comparison of the ML focus estimator ",
                          "($\\hat\\psi$) and the MBC focus estimator ",
                          "($\\tilde\\psi$) of the joint indirect effect ",
                          "for $\\tau = ", tables$tau[s], "$ and ",
                          "$n_1 =", A, "$, in terms of simulation-based estimates ",
                          "(see Example~\\ref{ex:mediation-analysis}) of probability of ",
                          "underestimation (PU), mean bias (BIAS), and root mean squared ",
                          "error (RMSE). The column ``Wald'' reports the coverage of ",
                          "nominally $95\\%$ Wald-type intervals using the estimated ",
                          "standard error evaluated at the supplied maximum likelihood ",
                          "estimate of $\\btheta$. The column ``HulC'' reports the coverage ",
                          "of nominally $95\\%$ HulC-type intervals with $\\Delta = 0$, ",
                          "constructed using the corresponding estimator. All summaries ",
                          "are $\\times 100$.")
        summaries_keep <- summaries |>
            filter(tau == tau_value, method %in% methods) |>
            mutate(method = factor(method, levels = methods, ordered = TRUE)) |>
            arrange(q, method, n)
        stopifnot(nrow(summaries_keep) > 0)
        q_sizes <- rle(summaries_keep$q)$lengths
        q_rows <- cumsum(c(1L, head(q_sizes, -1L)))
        method_sizes <- rle(paste(summaries_keep$q,
                                  summaries_keep$method))$lengths
        method_rows <- cumsum(c(1L, head(method_sizes, -1L)))
        tab <- summaries_keep |>
            mutate(
                method = unname(estimator_names[as.character(method)]),
                pu = pu * 100,
                bias = bias * 100,
                rmse = rmse * 100,
                cover_wald = cover_wald * 100,
                cover_hulc = cover_hulc * 100
            ) |>
            select(q, method, n, pu, bias, rmse, cover_wald, cover_hulc) |>
            setNames(new_names) |>
            tt(caption = caption) |>
            format_tt(j = 4:8, digits = 1,
                      num_fmt = "decimal", num_zero = TRUE) |>
            format_tt(replace = TRUE) |>
            style_tt(align = "rlrrrrrr", bold = FALSE) |>
        style_tt(i = "colnames", align = "c", bold = FALSE)
        ## Multirows
        for (i in seq_along(q_rows)) {
            tab <- tab |>
                style_tt(i = q_rows[i], j = 1,
                         rowspan = q_sizes[i], alignv = "t")
        }
        for (i in seq_along(method_rows)) {
            tab <- tab |>
                style_tt(i = method_rows[i], j = 2,
                         rowspan = method_sizes[i], alignv = "t")
        }
        ## Midrules between estimators
        estimator_rows <- setdiff(method_rows, q_rows)
        if (length(estimator_rows) > 0) {
            tab <- tab |>
                style_tt(i = estimator_rows, j = 2:8,
                         line = "t", line_width = 0.05)
        }
        ## Midrules between q blocks
        if (length(q_rows) > 1) {
            tab <- tab |>
                style_tt(i = q_rows[-1], j = 1:8,
                         line = "t", line_width = 0.05)
        }
        tab |>
            theme_latex(
                outer = paste0("label={tab:", tables$name[s], "}"),
                placement = "t!"
            ) |>
            save_tt(out_file, overwrite = TRUE)
    }
}

