# # Target Trial Emulation (TTE) using Inverse Probability Clone Censoring Weighting (IPCW) and Regularized Regression of Flatiron database for advanced NSCLC patients treated with immunotherapy

################################# 
rm(list=ls()); set.seed(20260717)

library(ggplot2); library(dplyr) # visual
library(glmnet); library(splines); library(survival) # modeling

dir.path <- file.path("grace-period-model") # grace-period-model
fig.path <- file.path("figures", "gp-model")

data <- readRDS(file.path(dir.path, "gp-data.rds"))
X <- data$X; X.long <- data$X.long; X.lc <- data$X.lc

# X: patient-level cohort data frame with baseline covariates, treatment types, dates of diagnosis and treatment, and event information (discontinuation, death, irAE)
#   covariates include: age, sex, race, ECOG performance status, smoking status, practice type, insurance type, initial cancer stage, cancer type, PD-L1 expression, histology, treatment type, IMDC risk score, and event types with corresponding months.
# X.long: long-format data frame representing the patient-level cohort in a person-month format, rows corresponding to specific months for each patient, indicators for being on immunotherapy, death events, irAE events, and any event occurrence, with baseline covariates merged from the patient-level cohort.
# X.lc: cloned longitudinal cohort data frame that includes two assigned strategies for each patient: discontinuation by month 27 and continuation through month 27. Contains artificial censoring indicators based on the assigned strategy...

# -------------- calculate IPCW --------------
# model probability of artificial censoring (censor.event=1)
# use natural splines (ns) for 'month' to allow for non-linear effects of months 24-27 hazard on censoring probability

# numerator: predicts censoring prob given strategy and time
m.num   <- glm( censor.event ~ assigned.strategy + ns(month, df = 3), 
                  family = binomial(link = "logit"), data = X.lc)

# denominator: predicts censoring prob given strategy, time, and baseline cofounders + time-varying covariates (on.ici)
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
X.lc$ipcw <- pmax( lb, pmin(ub, X.lc$ipcw) )

p.weights <- ggplot(X.lc, aes(x = ipcw)) +
  geom_histogram(bins = 100, fill = "#3B7EA1", color = "white") +  
  geom_vline(xintercept = 1, linetype="dashed", color = "darkred", linewidth=1) + 
  theme_minimal() +
  labs(title = "Diagnostic: Truncated IPCW Distribution", x = "IPCW", y = "Person-Months")
print(p.weights)
ggsave(file.path(fig.path, "tte_ipcw_distribution.png"), width = 10, height = 6, dpi = 300)



# -------------- pooled logistic elastic net ---------------------
X.mat <- model.matrix( any.event ~ assigned.strategy + ns(month, df = 3) + 
                          age + sex + ecog + cancer.type + practice.type + on.ici - 1,
                          data = X.lc)

Y.vec <- X.lc$any.event
W.vec <- X.lc$ipcw

# do not penalize "assigned.strategy" nor ns(month, df=3)
pen.factors <- rep(1, ncol(X.mat)); names(pen.factors) <- colnames(X.mat)
# identify exact colnames of variables to not penalize
unpenalized.vars <- c("assigned.strategydisc.2yr", "on.ici", 
                      "ns(month, df = 3)1", "ns(month, df = 3)2", "ns(month, df = 3)3")
pen.factors[names(pen.factors) %in% unpenalized.vars] <- 0 


# fit cross validated elastic net
cv.fit <- cv.glmnet(x = X.mat, y = Y.vec, weights = W.vec, 
                    family = "binomial", alpha = 0.5, nfolds = 10,
                    penalty.factor = pen.factors)
# penlization path plot
png(file.path(fig.path, "tte_elastic_net_cv_curve.png"), width = 10, height = 6, units = "in", res = 300)
plot(cv.fit); title("Diagnostic: Elastic Net CV Curve", line = 2.5, cex.main = 1.5, cex.lab = 1.2)
dev.off()

# extract nonzero coef at 1 SE optimal lambda
cat("\nSelected Coefficients (Lambda 1SE with Penalty Factors):\n")
coef.1se <- coef(cv.fit, s = "lambda.1se")
print(coef.1se[coef.1se[, 1] != 0, , drop = FALSE])

# -------------- standardized marginal survival (g-computation) ---------------------
# create baseline cohort of unique patients to act as the reference population for g-computation
base.cohort <- X.lc[!duplicated(X.lc$patient.id), ]
# define prediciton horizon (months 21-48)
months.seq <- 21:48

# replicate baseline for each month and both strategies
pred.grid <- base.cohort[rep(seq_len(nrow(base.cohort)), each = length(months.seq)), ]
pred.grid$month <- rep(months.seq, times = nrow(base.cohort))

# estimate marginal survival under a specific treatment strategy 
est.survival <- function(strategy.name, stop.ici.month=NULL){
    tmp.data <- pred.grid
    tmp.data$assigned.strategy <- factor(strategy.name, levels=c("cont.2yr", "disc.2yr"))

    # enforce time-varying treatment status (on.ici) based on assigned.strategy
    if(is.null(stop.ici.month)){
      # continue ici 
      tmp.data$on.ici <- 1
    } else {
      # if discontinue, evaluate the month threshold safely
      tmp.data$on.ici <- ifelse(tmp.data$month >= stop.ici.month, 0, 1)
    }

    # build design matrix as specified by glmnet training
    X.pred <- model.matrix( ~ assigned.strategy + ns(month, df=3) +
                              age + sex + ecog + cancer.type + practice.type + on.ici - 1,
                              data=tmp.data)
    
    # using 1SE lambda, predict monthly probability of event
    pred.link <- predict(cv.fit, newx = X.pred, s = "lambda.1se", type="response")
    tmp.data$p.event <- as.numeric(pred.link)
    # monthly survival probs 
    tmp.data$p.surv <- 1 - tmp.data$p.event
    # sort to ensure cumprod uses chronological per patient
    tmp.data <- tmp.data[order(tmp.data$patient.id, tmp.data$month), ]
    # cumulative survival provs per patient over time
    tmp.data$cum.surv <- ave(tmp.data$p.surv, tmp.data$patient.id, FUN=cumprod)
    # marginalize/mean across all patients at each month
    marg.surv <- aggregate(cum.surv ~ month, data=tmp.data, FUN=mean)
    marg.surv$strategy <- strategy.name
    return(marg.surv)
}

# execute counterfactual predictions 
# world A: everybody continues ici treatment
surv.cont <- est.survival(
  strategy.name =  "cont.2yr", 
  stop.ici.month = NULL
)
# world B: everybody discontinues ici treatment
surv.disc <- est.survival(
  strategy.name =  "disc.2yr", 
  stop.ici.month = 24
)
surv.results <- rbind(surv.cont, surv.disc)

p.survival <- ggplot(surv.results, aes(x = month, y = cum.surv, color = strategy)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("cont.2yr" = "#1f78b4", "disc.2yr" = "#fb9a99"),
                     labels = c("Continue at 2 Yrs", "Discontinue at 2 Yrs")) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(title = "Estimated Marginal Event-Free Survival (Standardized)",
       x = "Follow-up Month",
       y = "Event-Free Survival Probability") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

print(p.survival)
ggsave(file.path(fig.path, "tte_marginal_survival.png"), width = 10, height = 6, dpi = 300)


# risk difference (causal contrast due to OR collapsibility) at end of followup
surv.48.cont <- surv.cont$cum.surv[surv.cont$month == 48]
surv.48.disc <- surv.disc$cum.surv[surv.disc$month == 48]
rd.48 <- surv.48.cont - surv.48.disc
cat(sprintf("Survival (Continue):    %.1f%%\n", surv.48.cont * 100))
cat(sprintf("Survival (Discontinue): %.1f%%\n", surv.48.disc * 100))
cat(sprintf("Risk Difference:        %.1f%% (Percentage Points)\n", rd.48 * 100))







