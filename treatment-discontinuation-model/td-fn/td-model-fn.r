# td-model-fn.r
#
# functionalized version of td-model.r, for use in simulation study replication loop. Takes the data list (as returned by td-sim-fn.r) 
# as input, and returns surv.results + diagnostics, with no plotting / file I/O, so its cheap and fast to call ~50 - 100 times.

library(dplyr); library(splines); library(survival)

fit.td.gcomp <- function(data.list, verbose=FALSE){
    # data.list is a list of data.frames, as returned by td-sim-fn.r. It contains:
    #   X, X.long, X.lc
    X <- data.list$X; X.long <- data.list$X.long; X.lc <- data.list$X.lc

    # ---- unordered factors, for INTERPRETABILITY (not a modelling change) ----
    # R gives ordered factors polynomial contrasts (ecog.L/.Q/.C). Those span the
    # same column space as treatment contrasts -- identical deviance, identical df,
    # identical fitted values (verified to 1e-15) -- but they cannot be read against
    # the DGP constants. Treatment contrasts give one coefficient per level, so
    # `ecog2` can be compared directly against its true +1.10.
    X.lc$ecog          <- factor(as.character(X.lc$ecog),          levels = c("0","1","2","3"))
    X.lc$pd.l1         <- factor(as.character(X.lc$pd.l1),         levels = c("Negative","Low","High"))
    X.lc$initial.stage <- factor(as.character(X.lc$initial.stage), levels = c("III","IV"))

    # -------------- calculate IPCW --------------
    # same as grace-period-model, but now censoring is based on each 6-month checkup 
    
    
    # numerator: predicts censoring probability given strategy and time only
    # (deliberately minimal -- standard stabilized-weights numerator spec)
    m.num <- glm( censor.event ~ assigned.strategy + ns(month, df = 3), 
                    family = binomial(link = "logit"), data = X.lc)
    # denominator: insurance.type + tx.combo added to account for potential confounding of censoring by these covariates
    m.denom <- glm( censor.event ~ assigned.strategy + ns(month, df = 3) +
                        age + sex + race + ecog + cancer.type + practice.type + 
                        insurance.type + tx.combo +
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
    
    # (IPCW distribution plot omitted here -- see td-model-fixed.r for the
    #  interactive/plotting version. Not needed for a repeated replication loop.)

    # -------------- pooled logistic marginal structural model (quasibinomial) ---------------------
    # assigned.strategy intentionally NOT here -- see td-model-fixed.r for the
    # full rationale (collinearity with on.ici among surviving clone-rows was
    # diluting/reversing the on.ici estimate; confirmed in
    # td-crossarm-diagnostics.r). It stays in m.num/m.denom above, where
    # censoring is genuinely a function of it.
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
    
    if (verbose) {
        print(summary(fit.msm.death))
        print(summary(fit.msm.irae))
    }

    # ---- tidy coefficient table (both outcomes, one data.frame) --------------
    # Returned so a replication loop can rbind across seeds with no extra plumbing.
    # These are DIAGNOSTIC: they never enter the estimand, which is produced by the
    # g-computation below. But they are how you check the outcome models are
    # recovering the DGP rather than absorbing an omitted variable.
    tidy.coef <- function(fit, outcome) {
        cc <- summary(fit)$coefficients
        data.frame(outcome   = outcome,
                   term      = rownames(cc),
                   estimate  = cc[, 1],
                   std.error = cc[, 2],
                   row.names = NULL, stringsAsFactors = FALSE)
    }
    coef.df <- rbind(tidy.coef(fit.msm.death, "death"),
                     tidy.coef(fit.msm.irae,  "irae"))

    # weight diagnostics, also returned -- extreme weights are the usual reason an
    # IPCW estimator degrades, and you want them per-replication, not just per-run.
    w <- X.lc$ipcw
    weight.diag <- data.frame(mean = mean(w), sd = sd(w), min = min(w), max = max(w),
                              p99 = unname(quantile(w, .99)),
                              pct.gt5 = mean(w > 5), n.nonfinite = sum(!is.finite(w)))

    # ------------- DIAGNOSTIC: aliased / non-identifiable coefficients ---------------------
    # check for any coeffients in glm fit that are NA (rank deficiency, e.g assigned.strategy x covariate stratum with ~0 events for rare irAE outcomes),
    # predict.glm() return NA for every row, not just rows in offending stratum. This is because R's NA propogates through the matrix multiplication and cumprod()
    na.coef.death <- names(coef(fit.msm.death))[is.na(coef(fit.msm.death))]
    na.coef.irae  <- names(coef(fit.msm.irae))[is.na(coef(fit.msm.irae))]
    if (length(na.coef.death) > 0) {
        cat("\n*** WARNING: fit.msm.death has aliased (NA) coefficients: ***\n")
        print(na.coef.death)
    }
    if (length(na.coef.irae) > 0) {
        cat("\n*** WARNING: fit.msm.irae has aliased (NA) coefficients: ***\n")
        print(na.coef.irae)
        cat("This is consistent with the irAE model being unstable due to rare\n",
            "events (see notes on survivor bias) -- check event counts by\n",
            "assigned.strategy x covariate below.\n")
        print(with(X.lc, table(assigned.strategy, irae.event)))
    }
    if (verbose) {
        cat("\nIPCW weight diagnostics (pre-truncation would have been dropped if NA):\n")
        cat("  n rows with non-finite ipcw before truncation clamp: ",
            sum(!is.finite(X.lc$cum.num / X.lc$cum.denom)), "\n")
    }

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
        # kept for bookkeeping; fit.msm.death/fit.msm.irae no longer use
        # assigned.strategy, so predict() ignores this -- on.ici below is what
        # actually sets the counterfactual.
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
    
        # diagnostic: how many predictions came back NA for each outcome, and at
        # which months? (predict.glm uses na.action = na.pass on newdata by
        # default, so NAs are *kept* here -- they only get silently dropped later,
        # at the aggregate() step, which is exactly what breaks monotonicity)
        n.na.death <- sum(is.na(tmp.data$p.event.death))
        n.na.irae  <- sum(is.na(tmp.data$p.event.irae))
        if (n.na.death > 0 || n.na.irae > 0) {
            cat(sprintf("  [%s] NA predictions -> death: %d rows, irae: %d rows (of %d)\n",
                        strategy.name, n.na.death, n.na.irae, nrow(tmp.data)))
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

    
        # guardrail: each per-strategy marginal death curve MUST be non-increasing
        # in month (up to floating point tolerance). If this trips, something
        # upstream is still wrong
        d <- diff(marg.surv$surv.death)
        if (any(d > 1e-8, na.rm = TRUE)) {
            bad.months <- marg.surv$month[-1][which(d > 1e-8)]
            cat(sprintf("  [%s] *** NON-MONOTONIC OS CURVE at month(s): %s ***\n",
                        strategy.name, paste(bad.months, collapse = ", ")))
        }
 
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
    if (verbose || !all(mono.check$is.monotonic)) {
        cat("\n--- OS curve monotonicity by strategy (should all be TRUE) ---\n")
        print(mono.check)
    }
    
    # (plots + point-estimate summary cat()s from td-model-fixed.r are omitted
    #  here -- this function is called many times in the replication loop; use
    #  td-model-fixed.r directly for a single, fully-annotated run.)
    
    return(list(
        surv.results   = surv.results,
        coef.df        = coef.df,       # tidy coefficients, both outcomes
        weight.diag    = weight.diag,   # IPCW distribution summary
        mono.check     = mono.check,
        na.coef.death  = na.coef.death,
        na.coef.irae   = na.coef.irae
    ))
 
} # end fit.td.gcomp()




 
# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
if (sys.nframe() == 0L) {
    
    cat("\n======================================================\n")
    cat(" Running Test Suite for fit.td.gcomp()\n")
    cat("======================================================\n\n")
    
    set.seed(42)
    
    # --------------------------------------------------------------------------
    # 1. Generate Synthetic Input Data (Matches simulate.td.data() output)
    # --------------------------------------------------------------------------
    cat("--> Generating synthetic cohort data for testing...\n")
    
    n_pts <- 80
    months <- 1:48
    strats <- c(paste0("disc.", seq(6, 42, by = 6), "mo"), "cont")
    
    # Build X (patient-level)
    X_mock <- data.frame(
        patient.id = sprintf("PT%05d", 1:n_pts),
        cancer.type = factor(sample(c("NSCLC", "Bladder", "Kidney", "Melanoma"), n_pts, replace = TRUE)),
        stringsAsFactors = FALSE
    )
    
    # Build X.lc (longitudinal cloned data)
    lc_grid <- expand.grid(
        patient.id = X_mock$patient.id,
        month = months,
        assigned.strategy = factor(strats, levels = strats),
        stringsAsFactors = FALSE
    )
    
    # Merge baseline variables into X.lc
    X.lc_mock <- merge(lc_grid, X_mock, by = "patient.id")
    
    # Add required covariates and event indicators
    X.lc_mock$clone.id <- paste0(X.lc_mock$patient.id, "_", X.lc_mock$assigned.strategy)
    X.lc_mock$age <- runif(nrow(X.lc_mock), 50, 80)
    X.lc_mock$sex <- factor(sample(c("Male", "Female"), nrow(X.lc_mock), replace = TRUE))
    X.lc_mock$race <- factor(sample(c("White", "Black", "Asian"), nrow(X.lc_mock), replace = TRUE))
    X.lc_mock$ecog <- factor(sample(0:2, nrow(X.lc_mock), replace = TRUE))
    X.lc_mock$practice.type <- factor(sample(c("Academic", "Community"), nrow(X.lc_mock), replace = TRUE))
    X.lc_mock$insurance.type <- factor(sample(c("Private", "Medicare"), nrow(X.lc_mock), replace = TRUE))
    X.lc_mock$tx.combo <- sample(0:1, nrow(X.lc_mock), replace = TRUE)
    X.lc_mock$on.ici <- sample(0:1, nrow(X.lc_mock), replace = TRUE)
    
    # Synthetic events (low probabilities to keep positive survival)
    X.lc_mock$censor.event <- rbinom(nrow(X.lc_mock), 1, 0.01)
    X.lc_mock$death.event <- rbinom(nrow(X.lc_mock), 1, 0.005)
    X.lc_mock$irae.event  <- rbinom(nrow(X.lc_mock), 1, 0.003)
    
    mock_data_list <- list(
        X = X_mock,
        X.long = NULL, # fit.td.gcomp relies on X and X.lc
        X.lc = X.lc_mock
    )
    
    # --------------------------------------------------------------------------
    # 2. Run fit.td.gcomp() Test Case
    # --------------------------------------------------------------------------
    cat("--> Fitting G-computation model on synthetic data...\n\n")
    res <- fit.td.gcomp(mock_data_list, verbose = TRUE)
    
    # --------------------------------------------------------------------------
    # 3. Structural Assertions
    # --------------------------------------------------------------------------
    cat("\n--> Validating model outputs...\n")
    
    stopifnot(
        "Result must be a list" = is.list(res),
        "Result must contain required keys" = all(c("surv.results", "mono.check", "na.coef.death", "na.coef.irae") %in% names(res)),
        "surv.results must be a data frame" = is.data.frame(res$surv.results),
        "surv.results must contain strategy, month, and cumulative survival columns" = 
            all(c("month", "surv.death", "surv.irae", "strategy") %in% names(res$surv.results)),
        "mono.check must evaluate monotonicity for all strategies" = nrow(res$mono.check) > 0
    )
    
    cat("✓ Test Passed: G-computation pipeline and output structures validated.\n\n")
    
    cat("======================================================\n")
    cat(" All tests passed successfully!\n")
    cat("======================================================\n\n")
}