# td-sensitivity-fn.r
# 
# Functionalized version of 03-td-sensitivity.r, for use in simulation study.
# Change cumulative incidence rate for each outcome at 36 months, to see what this does on the estimated survival curves. 
# per scenario in a grid of target 36 month cumulative incidence rates (one axis for death, one axis for irAE):
#   1. Calibrates the death/irae hazard intercept shifts (see death.intercept.shift / irae.intercept.shift in in td-model-fn.r) 
#       So the simulated 36 month cumulative incidence of each outcome hits the target, via uniroot.() (this is the "calibration" step)
#       on a fixed-seed oracle cohort (fixed seed makes the calibration objective a monotone, noise-free function of the intercept shift, see calibrate.intercept.shift()
#   2. Computes GROUND-TRUTH 36mo survival under each target trial strategy directly, by simulating a large cohort under FORCED adherence regime (force.ici.until.month)
#       rather than natural behavior —— this is the actual estimand est.survival() in td-model.r is trying to recover, so its the right truth to compare against.
#   3. Run n.reps full replications of the actual estimator (simulate -> IPCW -> pooled logistic MSM -> g-comp, via fit.td.gcomp() from td-model-fn.r, 
#       which includes the aggregate fix) at different seeds and summarizes the 36 month survival estimates per strategy per outcome mean, SD and SE = SD/sqrt(n.reps)
#       across replications, plus bias against the ground truth from step 2/
#   4. Writes tidy CSV + RDS results to disk (results/)
# 
# 
# NOTE: this is the full simulation study, which is slow. defaults are (3x3 scenario grid x 20 reps , n = 7837/ replication), and will take a while
# since the core simulator is uses nested per-patient, per-month for loops. Start with the quick test at the bottom of the CONFIG block to confirm everything runs,
# then scale up n.reps / scenario grid as desired.
########################################
library(dplyr)

proj.dir  <- file.path(rprojroot::find_root(rprojroot::is_git_root), "treatment-discontinuation-model")
util.path <- file.path(proj.dir, "static", "util.r")
source(file.path(proj.dir, "td-fn", "td-sim-fn.r"))    # defines simulate.td.data()
source(file.path(proj.dir, "td-fn",  "td-model-fn.r"))  # defines fit.td.gcomp()
 
results.dir <- file.path(proj.dir, "results")
if (!dir.exists(results.dir)) dir.create(results.dir, recursive = TRUE)


#-----------------------calibrate.intercept.shift()-----------------------
# finds the hazard intercept shift that makes the simulated 36month cumulative incidence of `outcome` equal to `target.ci` 
# holding other outcome's shifts fixed at `other.shift` (death/irAE compete for follow-up time, so calibrate death first, 
# then irAE conditional on the calibrated death shift—— see driver loop below)
calibrate.intercept.shift <- function( target.ci,
                                        outcome = c("death", "irae"),
                                        other.shift = 0,
                                        n.calib = 8000,
                                        calib.seed = 1,
                                        shift.range = c(-3,3),
                                        tol = 0.002, 
                                        util.path = util.path){

    outcome <- match.arg(outcome)
    
    # ci36() simulates a large cohort under the given intercept shift of linear predictor (offset), 
    # and returns the 36 month cumulative incidence of the specified outcome (death or irAE)
    ci36 <- function(shift){
        ds.shift <- if (outcome == "death") shift else other.shift
        is.shift <- if (outcome == "irae") shift else other.shift
        oracle <- simulate.td.data(
            seed = calib.seed,
            n=n.calib,
            death.intercept.shift = ds.shift,
            irae.intercept.shift = is.shift,
            oracle.only = TRUE,
            util.path = util.path, verbose = FALSE)
        X <- oracle$X
        if (outcome == "death") {
            mean(!is.na(X$death.time) & X$death.time <= 36)
        } else if (outcome == "irae") {
            mean(!is.na(X$irae.time) & X$irae.time <= 36)
        }
    }

    # this function is the root-finding objective: it returns the difference between the simulated 36 month cumulative 
    # incidence and the target cumulative incidence. The root of this function is the intercept shift that achieves the target.
    obj <- function(shift) ci36(shift) - target.ci
    lo <- shift.range[1]; hi <- shift.range[2]
    f.lo <- obj(lo); f.hi <- obj(hi) 
    if (sign(f.lo) == sign(f.hi)) {
            stop(sprintf(
                "target.ci = %.3f for '%s' unreachable in shift.range [%.1f, %.1f] (CuI at bounds: %.3f, %.3f). Widen shift.range.",
                target.ci, outcome, lo, hi, f.lo + target.ci, f.hi + target.ci))
    }

    #  it finds the root of the function obj(shift) = ci36(shift) - target.ci, i.e. the shift that makes ci36(shift) = target.ci
    root <- uniroot(obj, interval = c(lo, hi), tol = tol)
    return(list(shift = root$root, achieved.ci = ci36(root$root), target.ci = target.ci))

}


#-----------------------CONFIG-----------------------
# full scenario (adjust n.reps / scenario.grid as desired, but this is slow)

scenarios <- expand.grid(
    target.ci.death.36 = c(0.15, 0.25, 0.35),
    target.ci.irae.36 = c(0.05, 0.10, 0.15)
)

n.reps     <- 20     # MC reps of the estimator per scenario
n.rep.size <- 7837   # patients per replication (matches original design)
n.calib    <- 8000   # patients used to calibrate each intercept shift
truth.n    <- 50000  # patients used for the ground-truth oracle per strategy
truth.seed <- 999
eval.month <- 36

# -- Quick smoke test (uncomment to sanity-check the pipeline runs first) --
scenarios <- expand.grid(target.ci.death.36 = 0.25, target.ci.irae.36 = 0.10)
n.reps <- 3; n.rep.size <- 1500; n.calib <- 1500; truth.n <- 5000


#-----------------------DRIVER-----------------------
strategies <- c(paste0("disc.", seq(6, 42, by=6), "mo"), "cont")
stop.months <- c(seq(6, 42, by=6), Inf)

scenario.calib.list <- list(); scenario.results.list <- list()

for (s in seq_len(nrow(scenarios))){

    target.death <- scenarios$target.ci.death.36[s]
    target.irae  <- scenarios$target.ci.irae.36[s]

    cat(sprintf("\n=== Scenario %d/%d: target 36-mo CI  death=%.0f%%  irae=%.0f%% ===\n",
                s, nrow(scenarios), target.death * 100, target.irae * 100))

    # --- calibrate intercept shifts for (death first, then irae conditional on death) --- 
    # this it finds the intercept shift that makes the simulated 36 month cumulative incidence of death equal to target.death,
    # then it finds the intercept shift that makes the simulated 36 month cumulative incidence of irAE equal to target.irae,
    # holding the death shift fixed at the calibrated value. This is because death and irAE compete for follow-up time, so we calibrate death first, then irAE conditional on the calibrated
    cal.death <- calibrate.intercept.shift(target.ci = target.death, outcome = "death", other.shift = 0,               n.calib = n.calib, util.path = util.path)
    cal.irae  <- calibrate.intercept.shift(target.ci = target.irae,  outcome = "irae",   other.shift = cal.death$shift, n.calib = n.calib, util.path = util.path)

    cat(sprintf("  calibrated shifts -> death: %+.3f (achieved %.1f%%)   irae: %+.3f (achieved %.1f%%)\n",
                cal.death$shift, cal.death$achieved.ci * 100,
                cal.irae$shift,  cal.irae$achieved.ci * 100))

    scenario.calib.list[[s]] <- data.frame(
        scenario = s,
        target.ci.death.36 = target.death,
        target.ci.irae.36  = target.irae,
        death.intercept.shift = cal.death$shift,
        irae.intercept.shift  = cal.irae$shift,
        achieved.ci.death.36 = cal.death$achieved.ci,
        achieved.ci.irae.36  = cal.irae$achieved.ci
    )


    # ---- ground truth per strategy, via forced-adherence oracle cohorts ---
    truth.rows <- vector("list", length(strategies))
    for (k in seq_along(strategies)){
        oracle <- simulate.td.data(
            seed = truth.seed,
            n = truth.n,
            death.intercept.shift = cal.death$shift,
            irae.intercept.shift  = cal.irae$shift,
            force.ici.until.month = stop.months[k],
            oracle.only = TRUE,
            util.path = util.path, verbose = FALSE)
        Xo <- oracle$X
        truth.rows[[k]] <- data.frame(
            scenario = s,
            strategy = strategies[k],
            true.surv.death.36 = mean(is.na(Xo$death.month) | Xo$death.month > eval.month),
            true.surv.irae.36  = mean(is.na(Xo$irae.month)  | Xo$irae.month  > eval.month)
        )
    }
    truth.df <- do.call(rbind, truth.rows)
    
    # --- MC reps for full estimator ---
    rep.rows <- vector("list", n.reps)
    for (r in seq_len(n.reps)){
        rep.seed <- 100000L * s + r
        est.out <- tryCatch({
            data.list <- simulate.td.data(
                seed = rep.seed,
                n = n.rep.size,
                death.intercept.shift = cal.death$shift,
                irae.intercept.shift  = cal.irae$shift,
                oracle.only = FALSE,
                util.path = util.path, verbose = FALSE)
        fit.td.gcomp(data.list, verbose=FALSE)
    }, error = function(e){
        cat(sprintf("  [rep %d/%d] FAILED: %s\n", r, n.reps, conditionMessage(e)))
        NULL
    })
    if (is.null(est.out)) next
    at36 <- est.out$surv.results[est.out$surv.results$month == eval.month, ]
    at36$scenario <- s
    at36$rep <- r
    at36$seed <- rep.seed
    rep.rows[[r]] <- at36
}

    rep.df <- do.call(rbind, rep.rows)
 
    if (is.null(rep.df) || nrow(rep.df) == 0) {
        cat(sprintf("  *** scenario %d: every replication failed, skipping summary ***\n", s))
        next
    }

    summary.df <- rep.df %>%
        group_by(scenario, strategy) %>%
        summarise(
            n.reps.completed    = dplyr::n(),
            mean.surv.death.36  = mean(cum.surv.death, na.rm = TRUE),
            sd.surv.death.36    = sd(cum.surv.death, na.rm = TRUE),
            se.surv.death.36    = sd.surv.death.36 / sqrt(n.reps.completed),
            mean.surv.irae.36   = mean(cum.surv.irae, na.rm = TRUE),
            sd.surv.irae.36     = sd(cum.surv.irae, na.rm = TRUE),
            se.surv.irae.36     = sd.surv.irae.36 / sqrt(n.reps.completed),
            .groups = "drop"
        ) %>%
        left_join(truth.df, by = c("scenario", "strategy")) %>%
        mutate(
            bias.death.36 = mean.surv.death.36 - true.surv.death.36,
            bias.irae.36  = mean.surv.irae.36  - true.surv.irae.36
        )
 
    cat("  --- 36-month summary (this scenario) ---\n")
    print(as.data.frame(summary.df))

    scenario.results.list[[s]] <- list(rep.df = rep.df, summary.df = summary.df)
} # end scenario loop




# ============================ SAVE RESULTS ============================
summary.all <- do.call(rbind, lapply(scenario.results.list, `[[`, "summary.df"))
reps.all    <- do.call(rbind, lapply(scenario.results.list, `[[`, "rep.df"))
calib.all   <- do.call(rbind, scenario.calib.list)
 
summary.all <- merge(summary.all, calib.all,
                      by = c("scenario"), suffixes = c("", ".calib"))
 
write.csv(summary.all, file.path(results.dir, "td_36mo_sensitivity_summary.csv"), row.names = FALSE)
write.csv(reps.all,    file.path(results.dir, "td_36mo_sensitivity_replicates.csv"), row.names = FALSE)
write.csv(calib.all,   file.path(results.dir, "td_36mo_sensitivity_calibration.csv"), row.names = FALSE)
saveRDS(list(summary = summary.all, replicates = reps.all, calibration = calib.all),
        file.path(results.dir, "td_36mo_sensitivity_results.rds"))
 
cat("\nSaved sensitivity results to:", normalizePath(results.dir), "\n")
cat("  - td_36mo_sensitivity_summary.csv     (mean / SD / SE / bias per scenario x strategy)\n")
cat("  - td_36mo_sensitivity_replicates.csv  (every individual replicate's 36-mo estimate)\n")
cat("  - td_36mo_sensitivity_calibration.csv (calibrated intercept shifts per scenario)\n")
cat("  - td_36mo_sensitivity_results.rds     (all three, as R objects)\n")
 









