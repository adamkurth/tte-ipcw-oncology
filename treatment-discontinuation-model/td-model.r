# # Target Trial Emulation (TTE) using Inverse Probability Clone Censoring Weighting (IPCW) and Regularized Regression of Flatiron database for advanced NSCLC patients treated with immunotherapy

# treatment-discontinuation-model
# When is the optimal time to discontinue ICI therapy?
################################# 
rm(list=ls()); set.seed(20260717)

library(ggplot2); library(dplyr) # visual
library(glmnet); library(splines); library(survival) # modeling

base.path <- rprojroot::find_root(rprojroot::is_git_root)
dir.path <- file.path(base.path, "treatment-discontinuation-model")
fig.path <- file.path(base.path, "figures", "td-model")
data.path <- file.path(dir.path, "td-data.rds")
data <- readRDS(data.path); X <- data$X; X.long <- data$X.long; X.lc <- data$X.lc

# X: patient-level cohort data frame with baseline covariates, treatment types, dates of diagnosis and treatment, and event information (discontinuation, death, irAE)
#   covariates include: age, sex, race, ECOG performance status, smoking status, practice type, insurance type, initial cancer stage, cancer type, PD-L1 expression, histology, treatment type, IMDC risk score, and event types with corresponding months.
# X.long: long-format data frame representing the patient-level cohort in a person-month format, rows corresponding to specific months for each patient, indicators for being on immunotherapy, death events, irAE events, and any event occurrence, with baseline covariates merged from the patient-level cohort.
# X.lc: cloned longitudinal cohort data frame that includes 8 assigned strategies (one for every 6 month check-up) for each patient: discontinuation at each 6-month checkup and continuation through each checkup. Contains artificial censoring indicators based on the assigned strategy...

# -------------- calculate IPCW --------------
# same as grace-period-model, but now censoring is based on each 6-month checkup 


# numerator: predicts censoring probability given strategy and time
m.num <- glm( censor.event ~ assigned.strategy + ns(month, df = 3), 
                  family = binomial(link = "logit"), data = X.lc)
# denominator: predicts censoring probability given strategy, time, and baseline cofounders + time-varying covariates (on.ici)
m.denom <- glm( censor.event ~ assigned.strategy + ns(month, df = 3) +
                    age + sex + race + ecog + cancer.type + practice.type + 
                    on.ici, # time varying status included
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
fit.msm.death <- glm(death.event ~ assigned.strategy + ns(month, df = 3) + 
                 age + sex + ecog + cancer.type + practice.type + on.ici,
               family = quasibinomial(link = "logit"), 
               weights = ipcw,
               data = X.lc)
fit.msm.irae <- glm(irae.event ~ assigned.strategy + ns(month, df = 3) +
                 age + sex + ecog + cancer.type + practice.type + on.ici,
               family = quasibinomial(link = "logit"), 
               weights = ipcw,
               data = X.lc)

print(summary(fit.msm.death))
print(summary(fit.msm.irae))



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

    # monthly survival probs 
    tmp.data$p.surv.death <- 1 - tmp.data$p.event.death
    tmp.data$p.surv.irae  <- 1 - tmp.data$p.event.irae
    
    # sort to ensure cumprod uses chronological per patient
    tmp.data <- tmp.data[order(tmp.data$patient.id, tmp.data$month), ]

    # cumulative survival provs per patient over time
    tmp.data$cum.surv.death <- ave(tmp.data$p.surv.death, tmp.data$patient.id, FUN=cumprod)
    tmp.data$cum.surv.irae  <- ave(tmp.data$p.surv.irae, tmp.data$patient.id, FUN=cumprod)
    
    # marginalize/mean across all patients at each month
    marg.surv <- aggregate(cbind(cum.surv.death, cum.surv.irae) ~ month, data=tmp.data, FUN=mean)
    marg.surv$strategy <- strategy.name
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

library(tidyr)
plot.data <- surv.results %>%
  pivot_longer(
    cols = c(cum.surv.death, cum.surv.irae),
    names_to = "event.type",
    values_to = "probability"
  ) %>%
  mutate(
    event.type = ifelse(event.type == "cum.surv.death", 
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

rd.death <- cont.48$cum.surv.death - disc.42$cum.surv.death
rd.irae  <- cont.48$cum.surv.irae - disc.42$cum.surv.irae

cat("\n--- END OF FOLLOW-UP (48 MONTHS) EVALUATION ---\n")
cat("OVERALL SURVIVAL (MORTALITY)\n")
cat(sprintf("  Survival (Continue):             %.2f%%\n", cont.48$cum.surv.death * 100))
cat(sprintf("  Survival (Discontinue at 42 Mo): %.2f%%\n", disc.42$cum.surv.death * 100))
cat(sprintf("  Risk Difference:                 %.2f%% (Percentage Points)\n", rd.death * 100))

cat("\nIRAE-FREE SURVIVAL (TOXICITY)\n")
cat(sprintf("  irAE-Free (Continue):             %.2f%%\n", cont.48$cum.surv.irae * 100))
cat(sprintf("  irAE-Free (Discontinue at 42 Mo): %.2f%%\n", disc.42$cum.surv.irae * 100))
cat(sprintf("  Risk Difference:                 %.2f%% (Percentage Points)\n", rd.irae * 100))

