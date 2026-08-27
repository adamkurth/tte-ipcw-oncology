# td-model.r
rm(list=ls()); set.seed(202608)

library(ggplot2); library(dplyr) # visual
library(glmnet); library(splines); library(survival) # modeling

base.path <- rprojroot::find_root(rprojroot::is_git_root)
dir.path <- file.path(base.path, "treatment-discontinuation-model")
fig.path <- file.path(base.path, "figures", "td-model")
data.path <- file.path(dir.path, "td-data.rds")
data <- readRDS(data.path); X <- data$X; X.long <- data$X.long; X.lc <- data$X.lc
################################# 

# X: patient-level cohort data frame with baseline covariates, treatment types, dates of diagnosis and treatment, and event information (discontinuation, death, irAE)
#   covariates include: age, sex, race, ECOG performance status, smoking status, practice type, insurance type, initial cancer stage, cancer type, PD-L1 expression, histology, treatment type, IMDC risk score, and event types with corresponding months.
# X.long: long-format data frame representing the patient-level cohort in a person-month format, rows corresponding to specific months for each patient, indicators for being on immunotherapy, death events, irAE events, and any event occurrence, with baseline covariates merged from the patient-level cohort.
# X.lc: cloned longitudinal cohort data frame that includes 8 assigned strategies (one for every 6 month check-up) for each patient: discontinuation at each 6-month checkup and continuation through each checkup. Contains artificial censoring indicators based on the assigned strategy...

# -------------- calculate IPCW --------------
# same as grace-period-model, but now censoring is based on each 6-month checkup 
 
 
# numerator: predicts censoring probability given strategy and time only
# deliberately minimal -- this is the standard stabilized-weights numerator
#   baseline confounder adjustment is done in the denominator, not numerator.
m.num <- glm( censor.event ~ assigned.strategy + ns(month, df = 3), 
                  family = binomial(link = "logit"), data = X.lc)

# denominator: predicts censoring probability given strategy, time, and baseline cofounders + time-varying covariates (on.ici)
#   insurance.type and tx.combo are added here because they're true drivers of natural discontinuation timing in the simulated data, 
#   hence, true drivers of who gets artificially censored, differentially by arm. Omitting them leaves informative-censoring bias that IPCW
#   cannot correct for. (diagnosed empirically in td-crossarm-diagnostics.r: this turned out to be a smaller factor than the assigned.strategy/on.ici
#   collinearity issue fixed here, but still is a real factor in data instead of this being thrown away as a local variable)
m.denom <- glm( censor.event ~ assigned.strategy + ns(month, df = 3) +
                    age + sex + race + ecog + cancer.type + practice.type + 
                    insurance.type + tx.combo + on.ici, # time-varying/baseline
                    family = binomial(link = "logit"), data = X.lc)

X.lc$p.uncensored.num <- 1 - predict(m.num, type = "response")
X.lc$p.uncensored.denom <- 1 - predict(m.denom, type = "response")
# sort data to ensure chronological order cumulative calculation
X.lc <- X.lc[order(X.lc$patient.id, X.lc$month),]

# calculate cumulative product by clone.id using ave()
X.lc$cum.num <- ave(X.lc$p.uncensored.num, X.lc$clone.id, FUN = cumprod)
X.lc$cum.denom <- ave(X.lc$p.uncensored.denom, X.lc$clone.id, FUN = cumprod)

# calculate stabilized IPCW weights
X.lc$ipcw <- X.lc$cum.num / X.lc$cum.denom

# truncate weights
lb <- quantile(X.lc$ipcw, 0.01, na.rm = TRUE)
ub <- quantile(X.lc$ipcw, 0.99, na.rm = TRUE) 
if (ub > 10 ) { ub <- 10 } # truncate upper bound to 10 to avoid extreme weights

X.lc$ipcw <- pmax( lb, pmin(ub, X.lc$ipcw) )
 
p.weights <- ggplot(X.lc, aes(x = ipcw)) +
  geom_histogram(bins = 100, fill = "#3B7EA1", color = "white") +  
  geom_vline(xintercept = 1, linetype="dashed", color = "darkred", linewidth=1) + 
  theme_minimal() +
  labs(title = "Diagnostic: Truncated IPCW Distribution (forced ub <= 10)", x = "IPCW", y = "Person-Months")
print(p.weights)
ggsave(file.path(fig.path, "tte_ipcw_distribution.png"), width = 10, height = 6, dpi = 300)




# -------------- pooled logistic marginal structural model (quasibinomial) ---------------------
# assigned.strategy is INTENTIONALLY NOT in these formaulas (it stays in m.num/m.denom above, but not here since censoring is a deterministic function of it).
#   in the true GDP, assigned.strategy has no DIRECT effect on death/irae except through on.ici and clone-censoring rule that means among rows that survive
#   being uncensored, on.ici is almost fully determined by assigned.strategy and month. Fitting both models as seperate additive main effects is collinear and 
#   splits the signal across the collinear predictors (diagnosed in td-crossarm-diagnostics.r, where on.ici's coef  went from stat. null and wrong-sided to
#   to correctly signed and stat. significant after removing assigned.strategy from the outcome model).
#   IPCW (via the strategy-conditional weights above is what corrects for the arm-level confounding here; the outcome model's job is to estimate the on.ici effect
#   and adjust for genuine baseline confounders. g-comp in est.survival() below still recovers each strategy counterfactual curve by manipulating on.ici directly,
#   —— it doesn't need assigned.strategy in the outcome model to do that.
fit.msm.death <- glm(death.event ~ ns(month, df = 3) + 
                 age + sex + ecog + cancer.type + practice.type +
                 initial.stage +          # in the true death DGP (+0.40 stage IV)
                 on.ici,
               family = quasibinomial(link = "logit"), 
               weights = ipcw,
               data = X.lc)
fit.msm.irae <- glm(irae.event ~ ns(month, df = 3) +
                 age + sex + ecog + cancer.type + practice.type +
                 tx.combo + pd.l1 +       # in the true irAE DGP (+0.25 each)
                 on.ici,
               family = quasibinomial(link = "logit"), 
               weights = ipcw,
               data = X.lc)

print(summary(fit.msm.death))
print(summary(fit.msm.irae))



# –------------ DIAGNOSTICL aliased / non-identifiable coefficients ---------------------
# check for any coeffients in glm fit that are NA (rank deficiency, e.g. assigned.strategy x covariate stratum with ~0 events for rare irAE outcomes), 
# predict.glm() return NA for every row, not just rows in offending stratum. This is because R's NA propogates through the matrix multiplication and 
# cumprod() for that patient from the first affected month onward. the downstream through aggregate() which can make unrelated death observed survival 
# curve non-monotonic by silently changing how many patients are averaged at each month.
na.coef.death <- names(coef(fit.msm.death))[is.na(coef(fit.msm.death))]
na.coef.irae <- names(coef(fit.msm.irae))[is.na(coef(fit.msm.irae))]

if(length(na.coef.death) > 0) {
    cat("\n*** WARNING: fit.msm.death has aliased (NA) coefficients: ***\n")
    print(na.coef.death)
} 
if(length(na.coef.irae) > 0) {
    cat("\n*** WARNING: fit.msm.irae has aliased (NA) coefficients: ***\n")
    print(na.coef.irae)
    cat("This is consistent with the irAE model being unstable due to rare\n", 
        "events (see notes for survival bias) -– check event counts by\n",
        "assigned.strategy x covariate below:\n")
    print(with(X.lc, table(assigned.strategy, ecog, irae.event)))
}

# also check for non-finite/NaN IPCW weights that silently get dropped from 
# *training* data by glm()'s na.action=na.omit default (weights are part of model frame), 
# which can starve specific strata of events. 
cat("\nIPCW weight diagnostics (pre-truncation would have been dropped if NA):\n")
cat("    n rows with non-finite IPCW before truncation clamp:",
          sum(!is.finite(X.lc$cum.num / X.lc$cum.denom)), "\n")
cat("    n rows with non-finite IPCW after truncation clamp:",
          sum(!is.finite(X.lc$ipcw)), "\n")


# -------------- pooled logistic elastic net ---------------------
# X.mat <- model.matrix( any.event ~ assigned.strategy + ns(month, df = 3) + 
#                           age + sex + ecog + cancer.type + practice.type + on.ici - 1,
#                           data = X.lc)
# Y.vec <- X.lc$any.event
# W.vec <- X.lc$ipcw

# # dynamically isolate all assigned.strategy cols
# strat.cols <- grep("^assigned\\.strategy", colnames(X.mat), value = TRUE)
# # do not penalize "assigned.strategy" nor ns(month, df=3)
# pen.factors <- rep(1, ncol(X.mat)); names(pen.factors) <- colnames(X.mat)
# # identify exact colnames of variables to not penalize
# unpenalized.vars <- c(strat.cols, "on.ici", 
#                       "ns(month, df = 3)1", "ns(month, df = 3)2", "ns(month, df = 3)3")
# pen.factors[names(pen.factors) %in% unpenalized.vars] <- 0 


# # fit cross validated elastic net
# cv.fit <- cv.glmnet(x = X.mat, y = Y.vec, weights = W.vec, 
#                     family = "binomial", 
#                     alpha = 0.5, 
#                     nfolds = 10,
#                     penalty.factor = pen.factors)

# # penlization path plot
# png(file.path(fig.path, "tte_elastic_net_cv_curve.png"), width = 10, height = 6, units = "in", res = 300)
# plot(cv.fit); title("Diagnostic: Elastic Net CV Curve", line = 2.5, cex.main = 1.5, cex.lab = 1.2)
# dev.off()

# # extract nonzero coef at 1 SE optimal lambda
# cat("\nSelected Coefficients (Lambda 1SE with Penalty Factors):\n")
# # coef.1se <- coef(cv.fit, s = "lambda.1se")
# coef.1se <- coef(cv.fit, s = "lambda.min") # so clinical covariates are not shrunk to zero and eliminated
# print(coef.1se[coef.1se[, 1] != 0, , drop = FALSE])



# -------------- standardized marginal survival (g-computation) ---------------------
# create baseline cohort of unique patients to act as the reference population for g-computation
base.cohort <- X.lc[!duplicated(X.lc$patient.id), ]
 
# define prediciton horizon (entire follow-up period)
months.seq <- 1:48
 
# replicate baseline for each month and both strategies
pred.grid <- base.cohort[rep(seq_len(nrow(base.cohort)), each = length(months.seq)), ]
pred.grid$month <- rep(months.seq, times = nrow(base.cohort))

# estimate marginal survival under a specific treatment strategy 
est.survival <- function(strategy.name, stop.ici.month=NULL){
    tmp.data <- pred.grid
    # ensure levels match simulated target arms 
    all.strats <- c(paste0("disc.", seq(6, 48, by=6), "mo"), "cont")
    tmp.data$assigned.strategy <- factor(strategy.name, levels=all.strats)

    # enforce time-varying treatment status (on.ici) based on assigned.strategy
    if(is.null(stop.ici.month)){
      # continue ici 
      tmp.data$on.ici <- 1
    } else {
      # if discontinue, evaluate the month threshold safely
      tmp.data$on.ici <- ifelse(tmp.data$month >= stop.ici.month, 0, 1)
    }
    
    # SIMPLIFIED PREDICTION
    tmp.data$p.event.death <- as.numeric(predict(fit.msm.death, newdata = tmp.data, type="response"))
    tmp.data$p.event.irae <- as.numeric(predict(fit.msm.irae, newdata = tmp.data, type="response"))

    # diagnostic: how many predictions came back NA for each outcome, and at which months? 
    # (predict glm uses na.action=na.pass on newdata by default, so NAs are kept here only to get 
    # silently dropped at aggregate() step, which is exactly what breaks the monotonicity.
    n.na.death <- sum(is.na(tmp.data$p.event.death))
    n.na.irae <- sum(is.na(tmp.data$p.event.irae))
    if(n.na.death > 0 || n.na.irae > 0){
      cat(sprintf("   [%s] NA predictions -> death: %d rows, irae: %d rows (of %d)\n",
                  strategy.name, n.na.death, n.na.irae, nrow(tmp.data)))
    }

    # FAIL OPTION B: 
    # # monthly survival probs 
    # tmp.data$p.surv.death <- 1 - tmp.data$p.event.death
    # tmp.data$p.surv.irae  <- 1 - tmp.data$p.event.irae
    
    # # sort to ensure cumprod uses chronological per patient
    # tmp.data <- tmp.data[order(tmp.data$patient.id, tmp.data$month), ]
 
    # # cumulative survival provs per patient over time
    # tmp.data$surv.death <- ave(tmp.data$p.surv.death, tmp.data$patient.id, FUN=cumprod)
    # tmp.data$surv.irae  <- ave(tmp.data$p.surv.irae, tmp.data$patient.id, FUN=cumprod)

    
    # ---- competing-risks cumulative incidence (estimand: option B) ----------
    # the simulator checks BEFORE irAE within a month, so surviving month j, means surviving both death and irAE in month j.
    # miroring that order exactly:
    # 
    #   S(m)     = prod_{j<=m} (1-pd(j)) * (1-pr(j))        free of BOTH events
    #   CIF_d(m) = sum_{j<=m}  S(j-1) * pd(j)
    #   CIF_r(m) = sum_{j<=m}  S(j-1) * (1-pd(j)) * pr(j)
    #
    # cif.death now depends on irae model through S(), two outcomes are no longer seperable. 
    tmp.data <- tmp.data[order(tmp.data$patient.id, tmp.data$month), ]

    p.d <- tmp.data$p.event.death
    p.r  <- tmp.data$p.event.irae
    q <- (1-p.d) * (1-p.r) # prob of surviving both events in month j

    # S(m-1): within-patient lagged cumulative product
    tmp.data$S.lag <- ave(q, tmp.data$patient.id, FUN = function(x) c(1, head(cumprod(x), -1))) # lagged cumprod, S(0) = 1
    tmp.data$S.both <- ave(q, tmp.data$patient.id, FUN = cumprod) # cumulative product of surviving both events
    
    # CIFs: cumulative incidence functions for death and irAE
    tmp.data$cif.death <- ave(tmp.data$S.lag * p.d,                 tmp.data$patient.id, FUN = cumsum)
    tmp.data$cif.irae  <- ave(tmp.data$S.lag * (1-p.d) * p.r,       tmp.data$patient.id, FUN = cumsum)

    # plotted quantities: complements of subdistributions
    tmp.data$surv.death <- 1 - tmp.data$cif.death
    tmp.data$surv.irae  <- 1 - tmp.data$cif.irae

    # keep separate aggregate() calls from orignal fix. they still stop one outcome's missingness
    # from contaminating (changing) which patients are averaged into the other outcome's marginal curve at each month.
    marg.death <- aggregate(surv.death ~ month, data = tmp.data, FUN = mean, na.action = na.omit)
    marg.irae  <- aggregate(surv.irae  ~ month, data = tmp.data, FUN = mean, na.action = na.omit)
    marg.surv <- merge(marg.death, marg.irae, by="month", all=TRUE)
    marg.surv$strategy <- strategy.name

    # ------- guardrail 1: exact identity -------
    ident <- aggregate(cbind(cif.death, cif.irae, S.both) ~ month, data = tmp.data, FUN = mean)
    off <- abs(ident$cif.death + ident$cif.irae + ident$S.both - 1)
    if (any(off > 1e-8)) {
        cat(sprintf("  [%s] *** CIF IDENTITY VIOLATED: max |CIF_d + CIF_r + S - 1| = %.3e ***\n",
                    strategy.name, max(off)))
    }

    # ---- guardrail 2: monotonicity (1 - CIF is non-increasing by construction) ----
    for (nm in c("surv.death", "surv.irae")) {
        dd <- diff(marg.surv[[nm]])
        if (any(dd > 1e-8, na.rm = TRUE)) {
            cat(sprintf("  [%s] *** NON-MONOTONIC %s at month(s): %s ***\n", strategy.name, nm,
                        paste(marg.surv$month[-1][which(dd > 1e-8)], collapse = ", ")))
        }
    }

    # ---------
    # # build design matrix as specified by glmnet training
    # X.pred <- model.matrix( ~ assigned.strategy + ns(month, df=3) +
    #                           age + sex + ecog + cancer.type + practice.type + on.ici - 1,
    #                           data=tmp.data)
    # # ------- Add missing columns to prediction matrix to match training matrix
    # # model.matrix() drops other 7 columns not present in the prediction data,
    # missing.cols <- setdiff(colnames(X.mat), colnames(X.pred))

    # # add missing cols and fill with 0s
    # if (length(missing.cols) > 0) {

    #     zero.mat <- matrix(0, nrow = nrow(X.pred), ncol = length(missing.cols))
    #     colnames(zero.mat) <- missing.cols
    #     X.pred <- cbind(X.pred, zero.mat)

    # }

    # # ensure the same order as training matrix
    # X.pred <- X.pred[, colnames(X.mat), drop = FALSE]
    # -------
    # using 1SE lambda, predict monthly probability of event
    # pred.link <- predict(cv.fit, newx = X.pred, s = "lambda.1se", type="response")
    # tmp.data$p.event <- as.numeric(pred.link)
    # ---------

    # fix: marginalize each outcome seperately, not with a single
    #   cbind(surv.death, surv.irae) ~ month, which silently drops any rows with NA in either outcome

    # aggregate() formula uses na.omit by default combined with cbind() a single NA in either column frops the entire
    # (patient, month) row from both outcomes' averages. Since the irAE model is fit on rare outcomes across the 
    # 8 assigned.strategy strata (more prone to sparse/aliased coefficients -> NA predictions), 
    # an irAE-side NA was silently shrinking and changing the composition of patient set being averaged into death obs. curve
    # at some months but not others. Since the patient's own surv.death is guarenteed non-inc (cumprod of values [0,1]), 
    # the only was the marginal curve can go non-monotonic is if the population being averaged isnt the same fixed cohort at every month
    # which is what was happening here.

    # aggregate seperately (each with own na.action) means that NA in irAE predictions can no longer contaminate the obs. curve, vice versa
    # use: na.action = pass + explicit counting instead of omitting, so missingness is visible rather than swallowed.
    
    # marg.death <- aggregate(surv.death ~ month, data = tmp.data, FUN = mean, na.action = na.omit)
    # marg.irae  <- aggregate(surv.irae  ~ month, data = tmp.data, FUN = mean, na.action = na.omit)

    # # marg.surv <- aggregate(cbind(surv.death, surv.irae) ~ month, data=tmp.data, FUN=mean)
    
    # n.months.death <- nrow(marg.death); n.months.irae <- nrow(marg.irae)
    # if(n.months.death < length(months.seq) || n.months.irae < length(months.seq)){
    #   cat(sprintf("  [%s] WARNING: dropped months in aggregation -> death: %d/%d, irae: %d/%d\n",
    #               strategy.name, n.months.death, length(months.seq), n.months.irae, length(months.seq)))

    # }
    # marg.surv <- merge(marg.death, marg.irae, by="month", all=TRUE)
    # marg.surv$strategy <- strategy.name

    # # guard: each per-strategy marginal death curve must be non-increasing in month
    # # (up to floating point tolerance). If this trips something upstream is still wrong
    # d <- diff(marg.surv$surv.death)
    # if (any(d > 1e-6, na.rm=TRUE)){
    #     bad.months <- marg.surv$month[-1][which(d > 1e-8)]
    #     cat(sprintf("  [%s] *** NON-MONOTONIC OS CURVE at month(s): %s ***\n",
    #           strategy.name, paste(bad.months, collapse = ", ")))
    # }
    
    return(marg.surv)
}


# execute counterfactual predictions for each target arm 
target.arms <- seq(6, 48, by=6)
surv.list <- lapply(target.arms, function(mo){
    if (mo == 48){
        est.survival(strategy.name = "cont", stop.ici.month = NULL)
    } else {
        est.survival(strategy.name = paste0("disc.", mo, "mo"), stop.ici.month = mo)
    }
})

surv.results <- do.call(rbind, surv.list)
surv.results$strategy <- factor(surv.results$strategy, 
                                levels = c(paste0("disc.", seq(6, 42, by=6), "mo"), "cont"))

# final consolidated monotonicity check across all strategies (belt-and-braces
# on top of the per-strategy check inside est.survival())
mono.check <- surv.results %>%
    arrange(strategy, month) %>%
    group_by(strategy) %>%
    summarise(is.monotonic = all(diff(surv.death) <= 1e-8), .groups = "drop")
cat("\n--- Final monotonicity check for OS (surv.death) across all strategies ---\n")
print(mono.check)
if (!all(mono.check$is.monotonic)) {
    cat("\nIf any FALSE remain after the aggregate() fix above, re-run the\n",
        "coefficient diagnostics: the irAE model is the more likely source\n",
        "(rare outcome, 8-level assigned.strategy). Consider dropping\n",
        "assigned.strategy interaction complexity or pooling adjacent arms\n",
        "for irAE specifically -- this doesn't have to affect the OS model.\n")
}


# corss-arm ordering check, at end of followup. Unlike within-curve checks above, this is not a gard guarentee 
# — adjacent arms with small true differences can legitimately swap places due to estimation noise even in well-spec models.
# treat this as a sanity read: not pass/fail assertion, we want to see if survival is trending up as target month inc., roughly monotonoically
# not scrambled (e.g., what we say before, short-term exposure coming out on top)
cat("\n---48-mo survival by strategy, in target-month order (should be roughly monotonic)---\n")
end.of.fu <- surv.results[surv.results$month == max(months.seq), ] %>%
  arrange(strategy) %>%
  select(strategy, surv.death)
print(as.data.frame(end.of.fu))
cat(sprintf("spearman correlation (target month vs. 48-mo survival) : %.3f  (expect strongly positive)\n",
            suppressWarnings(cor(as.integer(end.of.fu$strategy), end.of.fu$surv.death, method = "spearman"))))

















library(tidyr)
plot.data <- surv.results %>%
  pivot_longer(
    cols = c(surv.death, surv.irae),
    names_to = "event.type",
    values_to = "probability"
  ) %>%
  mutate(
    event.type = ifelse(event.type == "surv.death", 
                        "Overall Survival (Mortality)", 
                        "irAE-Free Survival (Toxicity)")
  )
  p.survival.multi <- ggplot(plot.data, aes(x = month, y = probability, color = strategy)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ event.type, scales = "free_y") + # Side-by-side panels
  scale_color_viridis_d(option = "turbo", 
                        labels = c(paste(seq(6, 42, by = 6), "Mo Disc."), "Continuous")) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  theme_minimal(base_size = 14) +
  labs(title = "Clinical Trade-off: Mortality vs. Toxicity Risk",
       subtitle = "Marginal standardized survival evaluating discontinuation timing",
       x = "Follow-up Month from ICI Initiation",
       y = "Event-Free Probability",
       color = "Discontinuation Timing:") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 12))
 
print(p.survival.multi)
ggsave(file.path(fig.path, "td_marginal_competing_risks.png"), width = 14, height = 7, dpi = 300)

# ----------------- SUMMARY STATISTICS -----------------
cont.48 <- surv.results[surv.results$month == 48 & surv.results$strategy == "cont", ]
disc.42 <- surv.results[surv.results$month == 48 & surv.results$strategy == "disc.42mo", ]
 
rd.death <- cont.48$surv.death - disc.42$surv.death
rd.irae  <- cont.48$surv.irae - disc.42$surv.irae
 
cat("\n--- END OF FOLLOW-UP (48 MONTHS) EVALUATION ---\n")
cat("OVERALL SURVIVAL (MORTALITY)\n")
cat(sprintf("  Survival (Continue):             %.2f%%\n", cont.48$surv.death * 100))
cat(sprintf("  Survival (Discontinue at 42 Mo): %.2f%%\n", disc.42$surv.death * 100))
cat(sprintf("  Risk Difference:                 %.2f%% (Percentage Points)\n", rd.death * 100))
 
cat("\nIRAE-FREE SURVIVAL (TOXICITY)\n")
cat(sprintf("  irAE-Free (Continue):             %.2f%%\n", cont.48$surv.irae * 100))
cat(sprintf("  irAE-Free (Discontinue at 42 Mo): %.2f%%\n", disc.42$surv.irae * 100))
cat(sprintf("  Risk Difference:                 %.2f%% (Percentage Points)\n", rd.irae * 100))å