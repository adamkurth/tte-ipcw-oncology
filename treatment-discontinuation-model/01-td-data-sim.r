# td-data-sim.r
# 
# Treatment discontinuation model data simulation for multi-arm discontinuation decision point.
#################################
rm(list = ls()); set.seed(20260813)


source(file.path("static", "util.r")) # load helpers
library(dplyr)
n <- 7837 # total n initiating treatment
# Follow-up design aligned with treatment initiation
month.start <- 1L
month.admin.end <- 48L
# discontinuation arms (months)
target.arms <- seq(6, 48, by = 6) # 6, 12, 18, 24, 30, 36, 42, 48 months


# ---------- BASELINE COVARIATES ----------
patient.id <- sprintf("PT%05d", seq_len(n))

age.raw <- rnorm(n, mean = 65, sd = 10)             # N(65, 10) truncated to [18, 95]
age <- round(pmin(pmax(age.raw, 18), 95), 1)
sex <- sample_factor(levels = c("Male", "Female"), probs = c(0.55, 0.45), n)
race <- sample_factor(levels = c("White", "Black", "Asian", "Hispanic", "Other"), probs = c(0.6, 0.15, 0.1, 0.1, 0.05), n)
smoking.status <- sample_factor(levels = c("Never", "Former", "Current"), probs = c(0.3, 0.5, 0.2), n)
practice.type <- sample_factor(levels = c("Academic", "Community"), probs = c(0.3, 0.7), n)
insurance.type <- sample_factor(levels = c("Private", "Medicare", "Medicaid", "Other"), probs = c(0.4, 0.4, 0.1, 0.1), n)
initial.stage <- sample_factor(levels = c("III", "IV"), probs = c(0.3, 0.7), n, ordered = TRUE)

# ECOG worsens with age by construction.
ecog <- rep(NA_integer_, n)
idx.lt55 <- which(age < 55)                 # idx for age < 55
idx.55_69 <- which(age >= 55 & age < 70)    # idx for age 55-69
idx.70_79 <- which(age >= 70 & age < 80)    # idx for age 70-79
idx.ge80 <- which(age >= 80)                # idx for age >= 80
ecog[idx.lt55] <- as.integer(as.character(sample_factor(levels = c(0, 1, 2, 3), probs = c(0.72, 0.24, 0.04, 0.0), n = length(idx.lt55), ordered = TRUE))) 
ecog[idx.55_69] <- as.integer(as.character(sample_factor(levels = c(0, 1, 2, 3), probs = c(0.52, 0.38, 0.08, 0.02), n = length(idx.55_69), ordered = TRUE)))
ecog[idx.70_79] <- as.integer(as.character(sample_factor(levels = c(0, 1, 2, 3), probs = c(0.36, 0.45, 0.14, 0.05), n = length(idx.70_79), ordered = TRUE)))
ecog[idx.ge80] <- as.integer(as.character(sample_factor(levels = c(0, 1, 2, 3), probs = c(0.24, 0.50, 0.21, 0.05), n = length(idx.ge80), ordered = TRUE)))
ecog <- factor(ecog, levels = c(0, 1, 2, 3), ordered = TRUE) # total n = 7837, ECOG distribution: 0=0.46, 1=0.38, 2=0.16

# ---------- CANCER-SPECIFIC COVARIATES ----------
cancer.levels <- c("NSCLC", "Bladder", "Kidney", "Melanoma")
cancer.type.proportions <- c(5980, 624, 589, 644) / n           # proportions for NSCLC, Bladder, Kidney, Melanoma
cancer.type <- factor( draw_stratified_types(n = n, labels = cancer.levels, probs = cancer.type.proportions), levels = cancer.levels)

# placeholders vectors for cancer-specific fields.
# NSCLC histology: squamous vs. non-squamous
histology <- factor(rep(NA_character_, n), levels = c("Squamous", "Non-squamous"))
treatment.type <- factor(
    rep(NA_character_, n),
    levels = c("ICI alone", "ICI + Chemo", "ICI + Enfortumab Vedotin", "ICI doublet", "ICI + Targeted Therapy")
)
imdc.risk.score <- factor(rep(NA_character_, n), levels = c("Favorable", "Intermediate", "Poor"), ordered = TRUE)
pd.l1 <- factor(rep(NA_character_, n), levels = c("Negative", "Low", "High"), ordered = TRUE)

# indexes for each cancer type to populate
idx.nsclc <- which(cancer.type == "NSCLC") 
idx.bladder <- which(cancer.type == "Bladder")
idx.kidney <- which(cancer.type == "Kidney")
idx.melanoma <- which(cancer.type == "Melanoma")

# dependency: NSCLC histology depends on smoking (current smokers more likely squamous, never smokers more likely non-squamous) 
for (i in idx.nsclc) {
    p.squamous <- if (smoking.status[i] == "Current") 0.45 else if (smoking.status[i] == "Former") 0.34 else 0.18   # probabilities for squamous histology based on smoking status (current/former/never)
    histology[i] <- sample_factor(c("Squamous", "Non-squamous"), c(p.squamous, 1 - p.squamous), 1)                  # squamous histology more common in current smokers, non-squamous more common in never smokers  
}

# dependency: PD-L1 expression depends on cancer type, with NSCLC having higher prevalence of PD-L1 testing and expression levels.
# PD-L1 expression is more common in NSCLC, less common in bladder and kidney cancers, and least common in melanoma.
pd.l1[idx.nsclc] <- sample_factor(c("Negative", "Low", "High"), c(0.45, 0.33, 0.22), length(idx.nsclc), ordered = TRUE)
pd.l1[c(idx.bladder, idx.kidney, idx.melanoma)] <- sample_factor(c("Negative", "Low", "High"), c(0.72, 0.20, 0.08), length(c(idx.bladder, idx.kidney, idx.melanoma)), ordered = TRUE)
# dependency: treatment patterns associated with insurance and practice settings
# NSCLC: community practice and Medicaid insurance more likely to receive ICI + Chemo; older age and ECOG 2 less likely to receive combination therapy.
# Bladder: academic practice and private insurance more likely to receive ICI + Enfortumab Vedotin.
# Kidney: treatment choice influenced by IMDC risk score, practice type, and insurance; poor risk more likely to receive ICI + Targeted Therapy.
# Melanoma: academic practice and private insurance more likely to receive ICI + Targeted Therapy.

# probability of receiving ICI + Chemo among NSCLC patients
# among NSCLC patients, community practice and Medicaid insurance are more likely to receive ICI + Chemo, ECOG == 2 and older age are less likely to receive combination therapy.
for (i in idx.nsclc) {
    lp.combo <- (-0.45 +
        0.35 * (practice.type[i] == "Community") +          # community practice more likely to receive ICI + Chemo 
        0.28 * (insurance.type[i] == "Medicaid") +          # Medicaid 
        0.18 * (insurance.type[i] == "Private") -           # private insurance less likely to receive ICI + Chemo
        0.25 * (age[i] >= 75) -                             # older age less likely 
        0.20 * (as.integer(as.character(ecog[i])) == 2))     # ECOG 2 less likely
    
    p.combo <- expit(lp.combo) 
    
    treatment.type[i] <- sample_factor(levels = c("ICI alone", "ICI + Chemo"), probs = c(1 - p.combo, p.combo), n = 1)
}
 


# probability of receiving ICI + Enfortumab Vedotin among bladder cancer patients
for (i in idx.bladder) {
    lp.ev <- (-1.0 +
        0.45 * (practice.type[i] == "Academic") +           # academic practice more likely to receive ICI + Enfortumab Vedotin
        0.25 * (insurance.type[i] == "Private") -           # private insurance less likely
        0.25 * (insurance.type[i] == "Medicaid"))            # Medicaid insurance less likely
    
    p.ev <- expit(lp.ev)
 
    treatment.type[i] <- sample_factor(levels = c("ICI alone", "ICI + Enfortumab Vedotin"), probs = c(1 - p.ev, p.ev), n = 1) 
}
 
 
# probability of receiving ICI doublet or ICI + Targeted Therapy among kidney cancer patients
for (i in idx.kidney) {
    # IMDC influenced by age/ECOG.
    lp.poor <- (-1.1 + 
                0.35 * (as.integer(as.character(ecog[i])) == 2) +   # ECOG 2x more likely to be poor risk (= 2)
                0.25 * (age[i] >= 75))                               # older ages more likely to be poor risk
    lp.intermediate <- (0.0 + 
                0.20 * (as.integer(as.character(ecog[i])) >= 1)) # ECOG levels 1 or 2 more likely to be intermediate risk
    
    imdc.risk.score[i] <- sample_multicat(levels = c("Favorable", "Intermediate", "Poor"), logits = c(0, lp.intermediate, lp.poor))
 
    # Treatment choice influenced by IMDC risk score, practice type, and insurance.
    logits.tx <- c(
        "ICI alone" = 0,                                            # baseline 
        "ICI doublet" = -0.3 + 
                        0.45 * (practice.type[i] == "Academic") +   # academic prac more likely to receive doublet
                        0.2 * (insurance.type[i] == "Private"),     # private insurance more likely to receive doublet
        "ICI + Targeted Therapy" = -0.7 +
                         0.55 * (imdc.risk.score[i] == "Poor") +    # poor IMDC risk more likely to receive targeted therapy 
                         0.2 * (insurance.type[i] != "Medicaid")    # non-Medicaid insurance more likely to receive targeted therapy
    )
 
    treatment.type[i] <- sample_multicat(levels = names(logits.tx), logits = logits.tx)
}
 
 
# probability of receiving ICI + Targeted Therapy among melanoma patients
for (i in idx.melanoma) {
    lp.targeted <- (-1.2 + 
                    0.45 * (practice.type[i] == "Academic") +   # academic practice more likely to receive targeted therapy
                    0.30 * (insurance.type[i] == "Private"))     # private insurance more likely to receive targeted therapy
 
    p.targeted <- expit(lp.targeted)
 
    treatment.type[i] <- sample_factor(levels = c("ICI alone", "ICI + Targeted Therapy"), probs = c(1 - p.targeted, p.targeted), n = 1)
}
 
 
head(treatment.type[idx.kidney], 10) # check treatment distribution for kidney cancer patients
head(treatment.type[idx.melanoma], 10) # check treatment distribution for melanoma patients
head(treatment.type[idx.nsclc], 10) # check treatment distribution for lung cancer patients
head(treatment.type[idx.bladder], 10) # check treatment distribution for bladder cancer patients
 
check.table <- table(cancer.type, treatment.type)
check.table # check treatment distribution by cancer type


# ---------- DATES AND TRIAL TIMING ----------
# starts ICI withing 120 days of advanced diagnosis, follow-up for 48 months (4 years) after ICI initiation
advanced.dx.date <- as.Date("2015-01-01") + sample(0:3652, n, replace = TRUE)   # random dx date between 2015-01-01 and 2024-12-31
ici.start.date <- advanced.dx.date + sample(0:120, n, replace = TRUE)           # random ICI start date within 120 days of advanced dx date
admin.end.date <- add_months_approx(d = ici.start.date, m = month.admin.end) # administrative censoring at 48 months after ICI start


# ----------- simulate discontinuation / competing events --------------
discontinue.month <- rep(Inf, n) 
death.month <- rep(Inf, n)
irae.month <- rep(Inf, n) 


# tx.combo: combination-therapy indicator. Computed once, vectorized, over
# all n patients so it can be stored as a real column in X (previously this
# was only a loop-local variable) 
combo.treatment.types <- c("ICI + Chemo", "ICI + Enfortumab Vedotin", "ICI doublet", "ICI + Targeted Therapy")
tx.combo.vec <- as.integer(treatment.type %in% combo.treatment.types)

for (i in seq_len(n)) {
 
    ecogi <- as.integer(as.character(ecog[i]))  # convert ECOG factor to integer
    stage4 <- (initial.stage[i] == "IV")        # indicator
    tx.combo <- tx.combo.vec[i]                 # combination-therapy indicator for this patient (see vectorized def above)
 
 
    for (m in month.start:month.admin.end) {
 
        # natural discontinuation organically scales with time
        # + spikes at clinical decision/scanning intervals (3, 6, 9, 12 months)
        lp.disc <- (-4.5 
            + 0.03 * (age[i] - 65)                      # older age slightly inc disc risk
            + 0.35 * (ecogi == 1 ) + 0.75 * (ecogi == 2)# ECOG 1 or 2 inc disc risk 
            + 0.20 * (practice.type[i] == "Community")  # community practice slightly inc disc risk
            + 0.20 * (insurance.type[i] == "Medicaid")  # Medicaid insurance slightly inc disc risk
            - 0.10 * (insurance.type[i] == "Private")   # private insurance slightly dec disc risk
            + 0.15 * tx.combo                           # combo therapy slightly inc disc risk
            + 0.50 * (m %in% seq(6, 42, by = 6))        # clinical scan interval spikes every 6mo
            + 0.02 * m)                                  # gradual increase in risk over time
 
        u <- runif(1); p.disc <- expit(lp.disc) # probability of discontinuation at month m
        if (u < p.disc ) { 
            discontinue.month[i] <- m; break
        }
 
    }
 
    for (m in month.start:month.admin.end) {
    
 
        on.ici <- as.integer(m < discontinue.month[i]) # indicator for still on ICI therapy
 
        lp.death <- (-4.9 
            + 0.05 * (age[i] - 65)                          # older age slightly inc death risk
            + 0.55 * (ecogi == 1 ) + 1.10 * ( ecogi == 2)   # ECOG 1 or 2 inc death risk
            + 0.40 * stage4                                 # stage IV inc death risk
            + 0.20 * (cancer.type[i] == "Bladder")          # bladder slightly inc death risk
            + 0.15 * (cancer.type[i] == "Kidney")           # kidney slightly inc death risk
            - 0.10 * (cancer.type[i] == "Melanoma")         # melanoma slightly dec death risk
            + 0.20 * (1 - on.ici)                           # TIME-TREND: off ICI therapy slightly inc death risk as time progresses
            + 0.02 * m)                                     # TIME-TREND: gradual increase in risk over time
 
        lp.irae <- (-5.4 
            + 0.30 * on.ici                                 # on ICI therapy inc IRAE risk
            + 0.25 * tx.combo                               # combo therapy inc IRAE risk
            + 0.20 * (sex[i] == "Female")                   # female sex slightly inc IRAE risk
            + 0.25 * (pd.l1[i] == "High")                   # high PD-L1 expression inc IRAE risk
            + 0.15 * (cancer.type[i] == "Melanoma")         # melanoma slightly inc IRAE risk
            + 0.20 * (m <= 12))                             # TIME-TREND: first year of therapy inc IRAE risk
 
        u <- runif(1); p.death <- expit(lp.death) # probability of death at month m
        if (u < p.death ) { 
            death.month[i] <- m; break
        }
        u <- runif(1); p.irae <- expit(lp.irae) # probability of IRAE at month m
        if (u < p.irae ) { 
            irae.month[i] <- m; break
        }
    
    }       
    
}

first.event.month <- pmin( death.month, irae.month, month.admin.end) # either death or IRAE, whichever comes first

event.type.chr <- rep("Administrative censoring", n)
event.type.chr[death.month < irae.month & death.month <= month.admin.end] <- "Death"
event.type.chr[irae.month < death.month & irae.month <= month.admin.end] <- "irAE"
event.type <- factor(event.type.chr, levels = c("irAE", "Death", "Administrative censoring"))

# ---------- Build patient-level cohort X ----------
X <- data.frame(
    patient.id = patient.id,
    age = age,
    sex = sex,
    race = race,
    ecog = ecog,
    smoking.status = smoking.status,
    practice.type = practice.type,
    insurance.type = insurance.type,
    initial.stage = initial.stage,
    cancer.type = cancer.type,
    pd.l1 = pd.l1,
    histology = histology,
    treatment.type = treatment.type,
    tx.combo = tx.combo.vec,
    imdc.risk.score = imdc.risk.score,
    advanced.dx.date = advanced.dx.date,
    ici.start.date = ici.start.date,
    admin.end.date = admin.end.date,
    discontinue.month = ifelse(is.finite(discontinue.month), discontinue.month, NA_integer_), 
    death.month = ifelse(is.finite(death.month), death.month, NA_integer_),
    irae.month = ifelse(is.finite(irae.month), irae.month, NA_integer_),
    event.type = event.type,
    event.month = ifelse(event.type == "Administrative censoring", month.admin.end, first.event.month),
    stringsAsFactors = FALSE
)




# ---------- DISCRETE-TIME LONG FORMAT X.long ----------
person.month.list <- vector("list", n) # list to hold person-month data for each patient
 
for ( i in seq_len(n)) {
 
 
    end.month <- X$event.month[i] # end month for this patient
    months <- month.start:end.month 
 
    if (length(months) == 0)  next
        
    disc.month.i    <- ifelse(is.na(X$discontinue.month[i]), Inf, X$discontinue.month[i])
    death.month.i   <- ifelse(is.na(X$death.month[i]), Inf, X$death.month[i])
    irae.month.i    <- ifelse(is.na(X$irae.month[i]), Inf, X$irae.month[i])
 
    on.ici          <- as.integer(months < disc.month.i) # indicator for still on ICI therapy
    death.event     <- as.integer(months == death.month.i) # indicator for death event
    irae.event      <- as.integer(months == irae.month.i) # indicator for IRAE event
    
    person.month.list[[i]] <- data.frame(
        patient.id = X$patient.id[i],
        cancer.type = X$cancer.type[i],
        month = months,
        interval.start.date = add_months_approx(d = X$ici.start.date[i], m = months - 1L),
        interval.end.date   = add_months_approx(d = X$ici.start.date[i], m = months),
        on.ici = on.ici,
        death.event = death.event,
        irae.event = irae.event,
        any.event = as.integer(death.event | irae.event),
        stringsAsFactors = FALSE
    )
 
}


X.long <- do.call(rbind, person.month.list)
rownames(X.long) <- NULL

# data structure checks
# ------ check 1: follow-up time constraints 
# patients should start at month 1 and not exceed month.admin.end = 48
cat("Max follow-up month (should be <= 48):", max(X.long$month), "\n")
cat("Min follow-up month (should be 1):", min(X.long$month), "\n")

# ------ check 2: terminal event integrity
# death / irae (any.event = 1) should only occur on the very last obs. month, cannot have rows after 
term_check <- X.long %>%
    arrange(patient.id, month) %>%
    group_by(patient.id) %>%
    summarise(
        total_events = sum(any.event),
        # if have event, is it strictly on their final row?
        valid_terminal = if_else(total_events > 0, last(any.event) == 1, TRUE),
        # ensure they never have more than 1 total event
        single_event = total_events <= 1
    ) 
cat("Are all events strictly terminal? should be TRUE: ", all(term_check$valid_terminal), "\n")
cat("Do patients have max 1 event? ", all(term_check$single_event), "\n")

# ------ check 3: treatment monotonicity (no restarting)
# once `on.ici` switches from 1 to 0, it must NEVER switch back to 1.
# prevents immortal time bias or logic errors in the discontinuation hazard.
restarts <- X.long %>% 
    arrange(patient.id, month) %>%
    group_by(patient.id) %>%
    mutate(ici_diff = on.ici - lag(on.ici, default = 1)) %>%
    filter(ici_diff > 0)
cat("Number of patients who erroneously restarted ICI: ", nrow(restarts), "\n")
# ------ check 4: correct event mapping
# ensure that any death/irAE events in long format perfectly match the baseline cohort's recorded event month
event_match_check <- X.long %>%
  filter(any.event == 1) %>%
  left_join(X, by = "patient.id") %>%
  mutate(month_match = (month == event.month))
 
cat("Do all longitudinal events match baseline event months? ", all(event_match_check$month_match), "\n")
        
# ... so far nothing really has changed besides how we define treatment arms & start/stop times. 
#   next step is to implement the model-ready data structure for multi-arm treatment discontinuation


# --------- Multi-Arm Clone-based artificial censoring for TTE/IPCW ----------
# create clones for each target discontinuation strategy 
clone.list <- lapply(target.arms, function(target.mo){
    strat.name <- if(target.mo == 48) "cont" else paste0("disc.", target.mo, "mo")
    tmp <- transform(X.long,
        assigned.strategy = strat.name,
        target.disc.month = target.mo
    )
    return(tmp)
})
X.lc <- do.call(rbind, clone.list)

# artificial censoring logic 
X.lc$artificial.censor <- 0L # no censoring by default

#   defines censoring rules for each clone based on their assigned strategy and observed behavior.
#   if a patient stops ICI therapy before their assigned target month (based on the original X.long data), they are considered to have "stopped too early" and are censored at that point.
#   Likewise, if a patient continues ICI therapy beyond their assigned target month (by next follow-up month) and is still on therapy, they are considered to have "failed to stop" and are censored at that point.
 
# 1. stopped too early (before target month, patient is not on ICI)
stopped.early <- ( X.lc$month < X.lc$target.disc.month & X.lc$on.ici == 0 )
# 2. stopped too late / failed to stop (after target month + 2 mo grace period, patient is still on ICI)
failed.to.stop <- ( X.lc$month > (X.lc$target.disc.month + 2) & X.lc$on.ici == 1 & X.lc$target.disc.month != 48 )  # only applies to discontinuation arms, not the "continue" arm
X.lc$artificial.censor <- as.integer(stopped.early | failed.to.stop) 
 
# keeps rows up to first artificial censoring for each clone
X.lc$clone.id <- paste0(X.lc$patient.id, "_", X.lc$assigned.strategy) # unique clone identifier
first.censor <- tapply(
    ifelse(X.lc$artificial.censor == 1, X.lc$month, Inf), 
    X.lc$clone.id,
    min
) # first artificial censoring month for each clone
 
X.lc$first.censor.month <- first.censor[X.lc$clone.id] # map back to X.lc
X.lc <- X.lc[X.lc$month <= X.lc$first.censor.month, ] # keep rows up to first artificial censoring
X.lc$censor.event <- as.integer(X.lc$month == X.lc$first.censor.month & is.finite(X.lc$first.censor.month)) # indicator for artificial censoring event
X.lc$first.censor.month <- NULL # remove temporary column
X.lc$target.disc.month <- NULL # clean up temp vars 
 

# merge baseline covariates into long format for modeling
baseline.vars <- c(
    "patient.id", "age", "sex", "race", "ecog", "smoking.status", "practice.type", "insurance.type", 
    "initial.stage", "cancer.type", "pd.l1", "histology", "treatment.type", "tx.combo", "imdc.risk.score"
)
X.long <- merge(X.long, X[, baseline.vars], by = c("patient.id", "cancer.type"), all.x = TRUE, sort=FALSE)
X.lc <- merge(X.lc, X[, baseline.vars], by = c("patient.id", "cancer.type"), all.x = TRUE, sort=FALSE)
 

# diagnostics
cat("Patient-level cohort X:\n  n patients:", nrow(X), "\n")
cat("\nCloned longitudinal cohort X.lc (TD Model):\n  n rows:", nrow(X.lc), "\n")
cat("  artificial censor events by multiple strategies:\n")
print(with(X.lc, table(assigned.strategy, censor.event)))
 
invisible(list(X = X, X.long = X.long, X.lc = X.lc))
saveRDS(list(X = X, X.long = X.long, X.lc = X.lc), file = file.path("treatment-discontinuation-model", "td-data.rds"))









# visual 
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
plot.path <- file.path("figures", "td-model", "td-data")
# Visualization 1: Mechanics of Artificial Censoring
# Show not just WHEN clones are censored but WHY (mechanism of deviation)
# this dictates how the IPCW models must be structured
 
censor.data <- X.lc %>% 
    filter(censor.event == 1) %>%
    mutate(
        # extract numeric month from strat (e.g., "disc.6mo" -> 6, "cont" -> Inf) 
        target.mo = ifelse( 
            assigned.strategy == "cont", Inf, 
            as.integer(stringr::str_extract(assigned.strategy, "\\d+"))
        ),
        censor.reason = case_when(
            month < target.mo ~ "Stopped Too Early (Off ICI)",
            month >= target.mo ~ "Failed to Stop (Continued ICI)",
            TRUE ~ "Other"
        )
    )
 
p1 <- ggplot(censor.data, aes(x = month, fill = censor.reason)) +
  geom_histogram(binwidth = 1, position = "stack", alpha = 0.85, color = "black") +
  facet_wrap(~ assigned.strategy, scales = "free_y", ncol = 3) +
  theme_minimal() +
  scale_fill_manual(values = c("Stopped Too Early (Off ICI)" = "#E69F00", 
                               "Failed to Stop (Continued ICI)" = "#56B4E9")) +
  labs(
    title = "Mechanism and Timing of Artificial Censoring",
    subtitle = "Deconstructing clone deviation: Pre-target discontinuation vs. Post-target failure to stop",
    x = "Month of Follow-up",
    y = "Number of Clones Artificially Censored",
    fill = "Deviation Type:"
  ) +
  theme(
    legend.position = "bottom", 
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )
 
print(p1)
ggsave(filename = file.path(plot.path, "artificial_censoring_distribution.png"), plot = p1, width = 12, height = 8, dpi = 300)
 
# ---------------------------------------------------------
# Visualization 2: Denominator Attrition (At-Risk Pool)
# ---------------------------------------------------------
# highlight the structural spartsity by cloning and censoring.
# visually this explains why the IPCW weights might become unstable in later time-points
attrition.data <- X.lc %>%
  group_by(assigned.strategy, month) %>%
  summarise(patients.at.risk = n(), .groups = "drop")
 
p2 <- ggplot(attrition.data, aes(x = month, y = patients.at.risk, color = assigned.strategy, group = assigned.strategy)) +
    geom_step(linewidth = 1.2, alpha = 0.8) + # geom_step better reflects discrete monthly intervals
    theme_minimal() +
    scale_color_viridis_d(option = "turbo") +
    scale_y_continuous(labels = scales::comma) +
    labs(
        title = "Structural Attrition of the At-Risk Pool",
        subtitle = "Visualizing the shrinking denominator available for IPCW estimation per month",
        x = "Follow-up Month (Discrete)",
        y = "Active Compliant Clones (N)",
        color = "Assigned Target Trial Arm:"
    ) +
    theme(
        legend.position = "bottom",
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(face = "bold")
    )
 
print(p2)
ggsave(filename = file.path(plot.path, "attrition_over_time.png"), plot = p2, width = 10, height = 6, dpi = 300)
 
 
# ---------------------------------------------------------
# Visualization 3: Naive (Unweighted) Target Trial Emulation
# ---------------------------------------------------------
# establish biased baseline. This proves informative censoring exists in the data. 
# because patients who are sick are more likely to stop early, skew these crude curves (i.e. using g-formula without IPCW) to show the bias in naive estimates.
 
km.summ <- X.lc %>%
  group_by(clone.id, assigned.strategy) %>%
  summarise(
    # The last observed month for this clone
    exit.month = max(month),
    # 1 if they had an actual event (Death/irAE), 0 if artificially censored
    event.status = max(any.event), 
    .groups = "drop"
  )
 
# Fit the crude KM model
strat.levels <- c(paste0("disc.", seq(6, 48, by = 6), "mo"), "cont")
km.summ$assigned.strategy <- factor(km.summ$assigned.strategy, levels = intersect(strat.levels, unique(km.summ$assigned.strategy))) # ensure consistent ordering
km_fit <- survfit(Surv(exit.month, event.status) ~ assigned.strategy, data = km.summ)
 
p3 <- ggsurvplot(
  km_fit,
  data = km.summ,
  risk.table = TRUE,
  palette = "turbo",
  ggtheme = theme_classic(),
  censor = FALSE, # Hide censor ticks to prevent clutter from massive artificial censoring
  title = "Biased Baseline: Naive Event-Free Survival",
  subtitle = "Crude survival strictly subject to informative censoring bias (Pre-IPCW)",
  xlab = "Months from ICI Initiation",
  ylab = "Unweighted Event-Free Probability",
  legend.title = "Target Strategy:",
  risk.table.height = 0.35,
  risk.table.y.text = FALSE # Keeps the risk table clean
)
print(p3)
ggsave(filename = file.path(plot.path, "naive_km_event_free_survival.png"), plot = p3$plot, width = 12, height = 8, dpi = 300)