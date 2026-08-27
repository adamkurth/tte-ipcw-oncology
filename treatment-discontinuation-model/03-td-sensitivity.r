# 03-td-sensitivity.r
#
# Sensitivity analysis for the treatment-discontinuation TTE/IPCW estimator.
#
# QUESTION: the estimator recovers the right answer on ONE simulated world.
#           Does it still work when death is rare, or common, or when toxicity
#           competes heavily for follow-up time?
#
# To answer that you need two things that do not exist for free:
#   (1) a KNOB   -- set the event rate to a chosen value          -> BLOCK 2
#   (2) a TRUTH  -- what the answer actually is                   -> BLOCK 3
# The estimator itself is one line (BLOCK 4). Everything else is those two.
#
##############################################################################

library(dplyr); library(splines); library(survival); library(ggplot2); library(tidyr)

base.path <- rprojroot::find_root(rprojroot::is_git_root)
proj.dir  <- file.path(base.path, "treatment-discontinuation-model")
util.path <- file.path(base.path, "static", "util.r")     # repo root, not proj.dir
fig.path  <- file.path(base.path, "figures", "td-model")
res.path  <- file.path(proj.dir, "results")
stopifnot(file.exists(util.path))                          # fail now, not in 20 minutes
dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)
dir.create(res.path, recursive = TRUE, showWarnings = FALSE)

source(file.path(proj.dir, "td-fn", "td-sim-fn.r"))       # simulate.td.data()
source(file.path(proj.dir, "td-fn", "td-model-fn.r"))     # fit.td.gcomp()


##############################################################################
# BLOCK 1: CONFIG
##############################################################################
target.death <- 0.25    # desired 36-mo cumulative incidence of DEATH
target.irae  <- 0.10    # desired 36-mo cumulative incidence of irAE

eval.month <- 36        # WISH TO KNOW 36 MONTH FOLLOW-UP
n.calib    <- 4000      # patients per calibration run (sets calibration precision)
truth.n    <- 30000     # patients for oracle (sets truth precision)
n.rep.size <- 7837      # patients per replication/realistic cohort
n.reps     <- 20        # mc reps
calib.seed <- 1
truth.seed <- 999

strategies  <- c(paste0("disc.", seq(6, 42, by = 6), "mo"), "cont")
stop.months <- c(seq(6, 42, by = 6), Inf)                  # 'cont' never stops

# TRUE DGP coefficients, for checking the outcome models recover them.
# Read straight off lp.death / lp.irae in td-fn/td-sim-fn.r.
# Terms with truth = 0 are covariates in the model but NOT in the DGP -- they are
# the control: if the model inflates everything, these move too.
true.coef <- rbind(
  data.frame(outcome = "death", term = c(
    "age","ecog1","ecog2","ecog3","initial.stageIV","cancer.typeBladder",
    "cancer.typeKidney","cancer.typeMelanoma","on.ici","sexFemale","practice.typeCommunity"),
    truth = c(0.05, 0.55, 1.10, 1.50, 0.40, 0.20, 0.15, -0.10, -0.20, 0, 0)),
  data.frame(outcome = "irae", term = c(
    "on.ici","tx.combo","sexFemale","pd.l1High","pd.l1Low","cancer.typeMelanoma",
    "age","ecog1","ecog2","practice.typeCommunity"),
    truth = c(0.30, 0.25, 0.20, 0.25, 0, 0.15, 0, 0, 0, 0))
)


##############################################################################
# BLOCK 2: CALIBRATION
# find the logit intercept shifts that deliver the desired 36-mo cumulative incidence
# Note: Death and irAE COMPETE for follow-up, so the two shifts must be solved together, by alternating.
##############################################################################

# one sim run -> 36-mo cumulative incidence of `outcome`, at a GIVEN pair of shifts
ci36 <- function(ds, is, outcome) {
    o <- simulate.td.data(seed = calib.seed, n = n.calib,
                          death.intercept.shift = ds, irae.intercept.shift = is,
                          oracle.only = TRUE, util.path = util.path, verbose = FALSE)$X
    m <- if (outcome == "death") o$death.month else o$irae.month
    mean(!is.na(m) & m <= eval.month)
}

lo <- -3; hi <- 3
cat("== BLOCK 2: CALIBRATION ==\n")
cat(sprintf("  bracket: CI at %.0f is %.3f | at %+.0f is %.3f | death target %.3f\n",
            lo, ci36(lo, 0, "death"), hi, ci36(hi, 0, "death"), target.death))

death.shift <- 0; irae.shift <- 0
tol.ci <- 0.005; max.pass <- 8 # tolerance on the CI, not the shift; 8 passes is plenty
for (pass in seq_len(max.pass)) {
    death.shift <- uniroot(function(s) ci36(s, irae.shift, "death") - target.death,
                           lower = lo, upper = hi, tol = 0.002)$root
    irae.shift  <- uniroot(function(s) ci36(death.shift, s, "irae") - target.irae,
                           lower = lo, upper = hi, tol = 0.002)$root
    achieved.death <- ci36(death.shift, irae.shift, "death")
    achieved.irae  <- ci36(death.shift, irae.shift, "irae")
    cat(sprintf("  pass %d: death %+.4f -> CI %.4f | irae %+.4f -> CI %.4f\n",
                pass, death.shift, achieved.death, irae.shift, achieved.irae))
    if (abs(achieved.death - target.death) < tol.ci &&
        abs(achieved.irae  - target.irae)  < tol.ci) break
}
if (abs(achieved.death - target.death) >= tol.ci || abs(achieved.irae - target.irae) >= tol.ci)
    stop("joint calibration did not converge -- widen [lo,hi] or raise max.pass")

# Report ACHIEVED, never target: uniroot's tol is on the shift axis, and the
# objective is a step function with resolution 1/n.calib.
cat(sprintf("  CALIBRATED -> death %+.4f (CI %.4f) | irae %+.4f (CI %.4f)\n\n",
            death.shift, achieved.death, irae.shift, achieved.irae))
calib.df <- data.frame(target.death, achieved.death, death.shift,
                       target.irae,  achieved.irae,  irae.shift)


##############################################################################
# BLOCK 3: GROUND TRUTH
#
# force.on.ici.until makes EVERY patient follow the strategy exactly.
# this is the estimand the g-computation is trying to recover.
##############################################################################
cat("== BLOCK 3: GROUND TRUTH (forced adherence) ==\n")
truth.rows <- vector("list", length(strategies))
for (k in seq_along(strategies)) {
    Xo <- simulate.td.data(seed = truth.seed, n = truth.n,
                           death.intercept.shift = death.shift,
                           irae.intercept.shift  = irae.shift,
                           force.on.ici.until = stop.months[k],
                           oracle.only = TRUE, util.path = util.path, verbose = FALSE)$X
    # ESTIMAND = 1 - CIF. A patient removed by an irAE has death.month = NA and
    # counts as NOT having died. This matches what the estimator produces.
    # Polarity matters: is.na(...), never !is.na(...).
    truth.rows[[k]] <- data.frame(
        strategy = strategies[k], month = 1:48,
        true.surv.death = sapply(1:48, function(m) mean(is.na(Xo$death.month) | Xo$death.month > m)),
        true.surv.irae  = sapply(1:48, function(m) mean(is.na(Xo$irae.month)  | Xo$irae.month  > m)))
    cat(sprintf("  %-10s done\n", strategies[k]))
}
truth.curve <- do.call(rbind, truth.rows)
truth.df    <- truth.curve[truth.curve$month == eval.month,
                           c("strategy","true.surv.death","true.surv.irae")]
print(truth.df, row.names = FALSE)

# more ICI exposure -> better survival, but MORE toxicity: opposite directions
cat(sprintf("  OS increasing across arms?        %s\n", all(diff(truth.df$true.surv.death) >= 0)))
cat(sprintf("  irAE-free decreasing across arms? %s\n", all(diff(truth.df$true.surv.irae)  <= 0)))
if (anyNA(truth.df)) stop("NA in truth.df -- check is.na() polarity in BLOCK 3")


##############################################################################
# BLOCK 4: REPLICATE THE ESTIMATOR (on normal data (confounded, censored))
# 
# Accumulates THREE things per rep: survival estimates, model coefficients,
# and IPCW weight diagnostics.
##############################################################################
cat("\n== BLOCK 4: ESTIMATOR REPLICATIONS ==\n")
surv.rows <- list(); coef.rows <- list(); wt.rows <- list()
for (r in seq_len(n.reps)) {
    seed.r <- 1000 + r
    ok <- tryCatch({
        dl <- simulate.td.data(seed = seed.r, n = n.rep.size,
                               death.intercept.shift = death.shift,
                               irae.intercept.shift  = irae.shift,
                               oracle.only = FALSE,      # need X.lc for IPCW
                               util.path = util.path, verbose = FALSE)
        fit.td.gcomp(dl, verbose = FALSE)
    }, error = function(e) { cat(sprintf("  rep %2d FAILED: %s\n", r, conditionMessage(e))); NULL })
    if (is.null(ok)) next

    s <- ok$surv.results; s$rep <- r; s$seed <- seed.r
    surv.rows[[length(surv.rows)+1]] <- s
    cf <- ok$coef.df;     cf$rep <- r; cf$seed <- seed.r
    coef.rows[[length(coef.rows)+1]] <- cf
    wd <- ok$weight.diag; wd$rep <- r; wd$seed <- seed.r
    wt.rows[[length(wt.rows)+1]] <- wd
    cat(sprintf("  rep %2d/%d done\n", r, n.reps))
}
surv.all <- do.call(rbind, surv.rows)     # every arm x month x rep
coef.all <- do.call(rbind, coef.rows)     # all coefficients, all seeds
wt.all   <- do.call(rbind, wt.rows)
cat(sprintf("  completed %d/%d replications\n", length(surv.rows), n.reps))


##############################################################################
# BLOCK 5: DIAGNOSTIC 1. do the outcome models recover the DGP?
#
# These coefficients never enter the estimand (the g-computation does that)
# But if they are wrong, the estimand is being built on a misspecified model,
# and any 'bias' you measure is misspecification, not estimator performance.
##############################################################################
cat("\n== BLOCK 5: COEFFICIENT RECOVERY ==\n")
coef.summary <- coef.all %>%
    group_by(outcome, term) %>%
    summarise(n.reps = dplyr::n(),
              mean.est  = mean(estimate),
              sd.est    = sd(estimate),                     # spread across seeds
              mc.se     = sd(estimate)/sqrt(dplyr::n()),    # MC error in mean.est
              mean.se   = mean(std.error),                  # model-reported SE, avg
              .groups = "drop") %>%
    left_join(true.coef, by = c("outcome","term")) %>%
    filter(!is.na(truth)) %>%                                # only terms with a known truth
    mutate(diff = mean.est - truth,
           t    = diff / mc.se,
           verdict = ifelse(abs(t) > 3, "BIASED", ifelse(abs(t) > 2, "borderline", "ok")))
print(as.data.frame(coef.summary %>%
      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat(sprintf("\n  terms recovered (|t| <= 2): %d / %d\n",
            sum(coef.summary$verdict == "ok"), nrow(coef.summary)))




saveRDS(list(calib = calib.df, truth = truth.curve, surv = surv.all,
             coef = coef.all, coef.summary = NULL,
             summary = NULL, weights = wt.all),
        file.path(res.path, "sens_results.rds"))


##############################################################################
# BLOCK 6: DIAGNOSTIC 2. IPCW weight behavior
# Extreme weights are the usual reason an IPCW estimator degrades.
##############################################################################
cat("\n== BLOCK 6: IPCW WEIGHTS (across reps) ==\n")
print(as.data.frame(wt.all %>% summarise(
    mean.of.mean = mean(mean), sd.of.mean = sd(mean),
    mean.max = mean(max), worst.max = max(max),
    mean.pct.gt5 = mean(pct.gt5), any.nonfinite = sum(n.nonfinite))) , row.names = FALSE)
cat("  (stabilized weights should centre near 1; a mean far from 1 means the\n",
    "   numerator/denominator models disagree more than they should)\n")


##############################################################################
# BLOCK 7: HEADLINE RESULT. bias of the estimator against the truth.
##############################################################################
cat("\n== BLOCK 7: BIAS AT", eval.month, "MONTHS ==\n")
summary.df <- surv.all %>%
    filter(month == eval.month) %>%
    group_by(strategy) %>%
    summarise(n.reps    = dplyr::n(),
              est.death = mean(surv.death), sd.death = sd(surv.death),
              se.death  = sd(surv.death)/sqrt(dplyr::n()),
              est.irae  = mean(surv.irae),  sd.irae  = sd(surv.irae),
              se.irae   = sd(surv.irae)/sqrt(dplyr::n()), .groups = "drop") %>%
    left_join(truth.df, by = "strategy") %>%
    mutate(bias.death = est.death - true.surv.death,
           bias.irae  = est.irae  - true.surv.irae,
           # a bias is only real if it clears 2x its own Monte Carlo error
           real.death = abs(bias.death) > 2*se.death,
           real.irae  = abs(bias.irae)  > 2*se.irae)
print(as.data.frame(summary.df %>% mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)
cat(sprintf("\n  mean |bias| death: %.4f (%.2f pp) | irae: %.4f (%.2f pp)\n",
            mean(abs(summary.df$bias.death)), 100*mean(abs(summary.df$bias.death)),
            mean(abs(summary.df$bias.irae)),  100*mean(abs(summary.df$bias.irae))))
cat(sprintf("  arms with |bias| > 2*MCSE -- death: %d/%d | irae: %d/%d\n",
            sum(summary.df$real.death), nrow(summary.df),
            sum(summary.df$real.irae),  nrow(summary.df)))


##############################################################################
# BLOCK 8: PLOTS
##############################################################################
theme_set(theme_minimal(base_size = 12))

# (1) coefficient recovery: estimate +/- 2 SD vs the true DGP value
p1 <- ggplot(coef.summary, aes(x = reorder(term, mean.est))) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_pointrange(aes(y = mean.est, ymin = mean.est - 2*sd.est, ymax = mean.est + 2*sd.est,
                        colour = verdict), size = .4) +
    geom_point(aes(y = truth), shape = 4, size = 3.2, stroke = 1.2, colour = "firebrick") +
    coord_flip() + facet_wrap(~ outcome, scales = "free_y") +
    scale_colour_manual(values = c(ok = "#2c7fb8", borderline = "#e6a700", BIASED = "firebrick")) +
    labs(title = "Outcome model coefficient recovery",
         subtitle = "point = mean across reps (bars 2 SD); red X = true DGP value",
         x = NULL, y = "log-odds", colour = NULL)
ggsave(file.path(fig.path, "sens_coef_recovery.png"), p1, width = 11, height = 6, dpi = 300)

# (2) estimated vs true survival curves, by arm
est.curve <- surv.all %>% group_by(strategy, month) %>%
    summarise(est.death = mean(surv.death), est.irae = mean(surv.irae), .groups = "drop")
cmp <- est.curve %>% left_join(truth.curve, by = c("strategy","month")) %>%
    pivot_longer(c(est.death, est.irae, true.surv.death, true.surv.irae),
                 names_to = "k", values_to = "p") %>%
    mutate(outcome = ifelse(grepl("death", k), "Overall survival", "irAE-free survival"),
           source  = ifelse(grepl("^est", k), "estimate", "truth"))
p2 <- ggplot(cmp, aes(month, p, colour = strategy, linetype = source)) +
    geom_line(linewidth = .7) + facet_wrap(~ outcome, scales = "free_y") +
    scale_linetype_manual(values = c(truth = "solid", estimate = "22")) +
    scale_colour_viridis_d(option = "turbo") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Estimator vs ground truth, by strategy",
         subtitle = "solid = forced-adherence truth; dashed = IPCW + g-computation estimate",
         x = "Month", y = "1 - CIF", colour = "Strategy", linetype = NULL)
ggsave(file.path(fig.path, "sens_truth_vs_est.png"), p2, width = 12, height = 6, dpi = 300)

# (3) bias with Monte Carlo error bars
bias.long <- summary.df %>%
    select(strategy, bias.death, se.death, bias.irae, se.irae) %>%
    pivot_longer(-strategy, names_to = c(".value","outcome"), names_pattern = "(bias|se)\\.(.*)") %>%
    mutate(outcome = recode(outcome, death = "Overall survival", irae = "irAE-free survival"),
           strategy = factor(strategy, levels = strategies))
p3 <- ggplot(bias.long, aes(strategy, bias)) +
    geom_hline(yintercept = 0, colour = "firebrick") +
    geom_pointrange(aes(ymin = bias - 2*se, ymax = bias + 2*se)) +
    facet_wrap(~ outcome) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Estimator bias against forced-adherence truth",
         subtitle = "bars = +/- 2 Monte Carlo SE; a bar crossing zero is indistinguishable from unbiased",
         x = NULL, y = "estimate - truth") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig.path, "sens_bias.png"), p3, width = 11, height = 5, dpi = 300)

# (4) replication spread -- how much would ONE study bounce around?
p4 <- ggplot(surv.all %>% filter(month == eval.month) %>%
             mutate(strategy = factor(strategy, levels = strategies)),
             aes(strategy, surv.death)) +
    geom_boxplot(outlier.size = .8, fill = "#deebf7") +
    geom_point(data = truth.df %>% mutate(strategy = factor(strategy, levels = strategies)),
               aes(y = true.surv.death), colour = "firebrick", shape = 4, size = 3.2, stroke = 1.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = sprintf("Spread of the %d-month estimate across %d replications", eval.month, n.reps),
         subtitle = "box = replication-to-replication variability; red X = truth",
         x = NULL, y = "1 - CIF(death)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig.path, "sens_replication_spread.png"), p4, width = 10, height = 5, dpi = 300)


##############################################################################
# BLOCK 9 -- SAVE
##############################################################################
write.csv(coef.all,     file.path(res.path, "sens_coefficients_all_reps.csv"), row.names = FALSE)
write.csv(coef.summary, file.path(res.path, "sens_coefficient_recovery.csv"),  row.names = FALSE)
write.csv(summary.df,   file.path(res.path, "sens_bias_summary.csv"),          row.names = FALSE)
write.csv(truth.curve,  file.path(res.path, "sens_truth_curves.csv"),          row.names = FALSE)
write.csv(wt.all,       file.path(res.path, "sens_weight_diagnostics.csv"),    row.names = FALSE)
saveRDS(list(calib = calib.df, truth = truth.curve, surv = surv.all,
             coef = coef.all, coef.summary = coef.summary,
             summary = summary.df, weights = wt.all),
        file.path(res.path, "sens_results.rds"))
cat("\nSaved results to", res.path, "and figures to", fig.path, "\n")
