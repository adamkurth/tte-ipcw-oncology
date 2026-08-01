# gp-data-sim.r
# 
# Grace period model data simulation for fixed 2-year discontinuation decision point.
#################################
rm(list = ls()); set.seed(20260717)

source(file.path("static", "util.r")) # load helpers

n <- 7837 # total n at 2-year decision point
# Follow-up design aligned with protocol language.
month.start <- 21L
month.window.end <- 27L
month.admin.end <- 48L

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
    lp.combo <- -0.45 +
        0.35 * (practice.type[i] == "Community") +          # community practice more likely to receive ICI + Chemo 
        0.28 * (insurance.type[i] == "Medicaid") +          # Medicaid 
        0.18 * (insurance.type[i] == "Private") -           # private insurance less likely to receive ICI + Chemo
        0.25 * (age[i] >= 75) -                             # older age less likely 
        0.20 * (as.integer(as.character(ecog[i])) == 2)     # ECOG 2 less likely
    
    p.combo <- expit(lp.combo) 
    
    treatment.type[i] <- sample_factor(levels = c("ICI alone", "ICI + Chemo"), probs = c(1 - p.combo, p.combo), n = 1)
}



# probability of receiving ICI + Enfortumab Vedotin among bladder cancer patients
for (i in idx.bladder) {
    lp.ev <- -1.0 +
        0.45 * (practice.type[i] == "Academic") +           # academic practice more likely to receive ICI + Enfortumab Vedotin
        0.25 * (insurance.type[i] == "Private") -           # private insurance less likely
        0.25 * (insurance.type[i] == "Medicaid")            # Medicaid insurance less likely
    
    p.ev <- expit(lp.ev)

    treatment.type[i] <- sample_factor(levels = c("ICI alone", "ICI + Enfortumab Vedotin"), probs = c(1 - p.ev, p.ev), n = 1) 
}


# probability of receiving ICI doublet or ICI + Targeted Therapy among kidney cancer patients
for (i in idx.kidney) {
    # IMDC influenced by age/ECOG.
    lp.poor <- -1.1 + 
                0.35 * (as.integer(as.character(ecog[i])) == 2) +   # ECOG 2x more likely to be poor risk (= 2)
                0.25 * (age[i] >= 75)                               # older ages more likely to be poor risk
    lp.intermediate <- 0.0 + 
                0.20 * (as.integer(as.character(ecog[i])) >= 1) # ECOG levels 1 or 2 more likely to be intermediate risk
    
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
    lp.targeted <- -1.2 + 
                    0.45 * (practice.type[i] == "Academic") +   # academic practice more likely to receive targeted therapy
                    0.30 * (insurance.type[i] == "Private")     # private insurance more likely to receive targeted therapy

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
# window for discontinuation is months 24-27, with administrative censoring at month 48.

advanced.dx.date <- as.Date("2015-01-01") + sample(0:3652, n, replace = TRUE) 
ici.start.date <- advanced.dx.date + sample(0:120, n, replace = TRUE)
month21.date <- add_months_approx(d = ici.start.date, m = month.start - 1L)     # m = 20 months after ICI, month.start = 21 
month24.date <- add_months_approx(d = ici.start.date, m = 24L)                  # m = 24 after ICI 
month27.date <- add_months_approx(d = ici.start.date, m = month.window.end)     # m = 27 after ICI, month.window.end = 27
admin.end.date <- add_months_approx(d = ici.start.date, m = month.admin.end)    # m = 48 after ICI, month.admin.end = 48

# ---------- Simulate discontinuation and competing events ----------
discontinue.month <- rep(Inf, n)            # discontinuation month (Inf if no discontinuation)
death.month <- rep(Inf, n)                  # death month (Inf if no death)
irae.month <- rep(Inf, n)                   # irAE month (Inf if no irAE events)




for (i in seq_len(n)) {
    ecogi <- as.integer(as.character(ecog[i]))      # convert ecog factor to integer
    stage4 <- as.integer(initial.stage[i] == "IV") 
    tx.combo <- as.integer(
        treatment.type[i] %in% c("ICI + Chemo", "ICI + Enfortumab Vedotin", "ICI doublet", "ICI + Targeted Therapy")
    )

    # DISCONTINUATION PROCESS: after month 21, with a protocol-window spike (months 24-27).
    # - discontinuation process is discrete-time hazard
    # - hazard depends on age, ECOG, practice type, insurance type, and treatment type
    # - protocol window of disc hazard spike in months 24-27 (protocol-specified discontinuation window)
    for (m in month.start:month.admin.end) {
        lp.disc <- -4.0 +
            0.03 * (age[i] - 65) +                          # age contribution to discontinuation hazard
            0.35 * (ecogi == 1) + 0.75 * (ecogi == 2) +     # ECOG = 1 or 2 increases disc. hazard
            0.20 * (practice.type[i] == "Community") +      # community practice increases disc. hazard
            0.20 * (insurance.type[i] == "Medicaid") -      # Medicaid decreases disc. hazard
            0.10 * (insurance.type[i] == "Private") +       # private decreases disc. hazard
            0.15 * tx.combo +                               # combo regimens (chemo/enfortumab/doublet/targeted) increases disc. 
            0.85 * (m >= 24 & m <= 27)                      # protocol window spike for months 24-27 
        
        # disc occurs with probability expit(lp.disc) at each month
        if (runif(1) < expit(lp.disc)) {
            discontinue.month[i] <- m
            break
        }
    
    }

    # COMPETING EVENTS: death and irAE events 
    # Cause-specific monthly hazards for competing events.
    for (m in month.start:month.admin.end) {
        on.ici <- as.integer(m < discontinue.month[i])      # indicator for patients still on ICI therapy (not discontinued) at month m 

        # DEATH HAZARD: depends on age, ECOG, stage, cancer type, and on-treatment status (on ICI vs. discontinued)
        lp.death <- -4.9 +                              # baseline hazard
            0.05 * (age[i] - 65) +                      # age slightly increases
            0.55 * (ecogi == 1) + 1.10 * (ecogi == 2) + # ECOG 1 or 2 increases 
            0.40 * stage4 +                             # stage IV increases
            0.20 * (cancer.type[i] == "Bladder") +      # bladder increases
            0.15 * (cancer.type[i] == "Kidney") -       # kidney decreases
            0.10 * (cancer.type[i] == "Melanoma") +     # melanoma increases
            0.20 * (1 - on.ici) +                       # discontinuation increases death hazard !!
            0.02 * (m - month.start)                    # time trend: death hazard increases slightly over time
        
        # irAE HAZARD: depends on on-treatment status, treatment type, sex, PD-L1 expression, cancer type, and time (higher hazard in earlier months)
        lp.irae <- -5.4 +
            0.30 * on.ici +                             # on-treatment increases
            0.25 * tx.combo +                           # combo regimens increase
            0.20 * (sex[i] == "Female") +               # female increases
            0.25 * (pd.l1[i] == "High") +               # high PD-L1 increases
            0.15 * (cancer.type[i] == "Melanoma") +     # melanoma increases
            0.20 * (m <= 30)                            # time trend: higher hazard in earlier months (<= month 30)

        p.death <- expit(lp.death)
        p.irae <- expit(lp.irae)

        u <- runif(1) 
        if (u < p.death) {                              # death occurs
            death.month[i] <- m
            break
        }
        if (u < p.death + (1 - p.death) * p.irae)  {    # irAE occurs only if death did not occur
            irae.month[i] <- m
            break
        }
    }
}

first.event.month <- pmin(death.month, irae.month, month.admin.end)
event.type <- ifelse(
    death.month < irae.month & death.month <= month.admin.end,
    "Death",
    ifelse(irae.month < death.month & irae.month <= month.admin.end, "irAE", "Administrative censoring")
)
event.type <- factor(event.type, levels = c("irAE", "Death", "Administrative censoring"))

# visualize by event type and month of first event
# event.counts <- table(first.event.month, event.type)
# event.colors <- c("irAE" = "#3B7EA1","Death" = "#C44536","Administrative censoring" = "#9AA0A6")
# barplot(
#     t(event.counts),
#     col = event.colors[colnames(event.counts)],
#     border = NA,
#     main = "First Event Month by Event Type",
#     xlab = "Month",
#     ylab = "Frequency"
# )
# legend(
#     "topright",
#     legend = names(event.colors),
#     fill = event.colors,
#     bty = "n"
# )
table(event.type)

# ---------- Build patient-level cohort X ----------
X <- data.frame(
    # baseline covariates
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
    # cancer-specific covariates
    pd.l1 = pd.l1,
    histology = histology,
    treatment.type = treatment.type,
    imdc.risk.score = imdc.risk.score,
    # dates
    advanced.dx.date = advanced.dx.date,
    ici.start.date = ici.start.date,
    month21.date = month21.date,
    month24.date = month24.date,
    month27.date = month27.date,
    admin.end.date = admin.end.date,
    # events 
    discontinue.month = ifelse(is.finite(discontinue.month), discontinue.month, NA_integer_),
    death.month = ifelse(is.finite(death.month), death.month, NA_integer_),
    irae.month = ifelse(is.finite(irae.month), irae.month, NA_integer_),
    event.type = event.type,
    event.month = ifelse(event.type == "Administrative censoring", month.admin.end, first.event.month),
    stringsAsFactors = FALSE
)


# ADD MISSINGNESS (UNCOMMENT WHEN READY)
# # EHR-like missingness.
# missing.p <- c(ecog = 0.05, pd.l1 = 0.12, smoking.status = 0.03)
# for (v in names(missing.p)) {
#     idx <- sample.int(n, size = floor(n * missing.p[[v]]))
#     X[idx, v] <- NA
# }

head(X)






# ---------- DISCRETE-TIME LONGITUDINAL FORMAT X_LONG ----------
# turn into person-month format
person.month.list <- vector("list", n)

for (i in seq_len(n)) {
    end.month <- X$event.month[i]
    months <- month.start:end.month

    if (length(months) == 0) next 

    # if event (disc, death, irae) is NA, set to Inf for comparison
    disc.month.i    <- ifelse( is.na(X$discontinue.month[i]), Inf, X$discontinue.month[i] )
    death.month.i   <- ifelse( is.na(X$death.month[i]), Inf, X$death.month[i] )
    irae.month.i    <- ifelse( is.na(X$irae.month[i]), Inf, X$irae.month[i] )
    
    # indicators for each month: on ici therapy, death event, irAE event
    on.ici      <- as.integer(months < disc.month.i)    # on ici therapy before disc
    death.event <- as.integer(months == death.month.i) 
    irae.event  <- as.integer(months == irae.month.i)

    person.month.list[[i]] <- data.frame(
        patient.id = X$patient.id[i], 
        cancer.type = X$cancer.type[i], 
        month = months, 
        interval.start.date = add_months_approx(X$ici.start.date[i], months - 1L), # ici start date w/ realistic buffer
        interval.end.date   = add_months_approx(X$ici.start.date[i], months),       # ici end date using month indicator
        on.ici = on.ici, 
        death.event = death.event, 
        irae.event = irae.event, 
        any.event = as.integer(death.event == 1 | irae.event == 1),
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
        

# ----------- Clone-based artificial censoring for TTE/IPCW

X.lc <- rbind(
    transform(X.long, assigned.strategy = "disc.2yr"), # discontinue by month 27
    transform(X.long, assigned.strategy = "cont.2yr")  # continue through & beyond month 27 
)

X.lc$artificial.censor <- 0L
idx.disc.strategy   <- which( X.lc$assigned.strategy == "disc.2yr")
idx.cont.strategy   <- which( X.lc$assigned.strategy == "cont.2yr") 


# strategy: disc by end of month=27
X.lc$artificial.censor[idx.disc.strategy] <- as.integer(
    # discontinue strategy: censor if month > 27 and still on ici therapy
    X.lc$month[idx.disc.strategy] > month.window.end & X.lc$on.ici[idx.disc.strategy] == 1
)

# strategy: continue through month=27
X.lc$artificial.censor[idx.cont.strategy] <- as.integer(
    # continue strategy: censor if month <= 27 and not on ici therapy
    X.lc$month[idx.cont.strategy] <= month.window.end & X.lc$on.ici[idx.cont.strategy] == 0
)


# keep rows up to first artificial censoring for each clone. 
X.lc$clone.id <- paste0(X.lc$patient.id, "__", X.lc$assigned.strategy)
first.censor <- tapply(
    ifelse( X.lc$artificial.censor == 1, X.lc$month, Inf), # if artificial censoring, return month, else Inf
    X.lc$clone.id,
    min
)

X.lc$first.censor.month <- first.censor[X.lc$clone.id] # map first censor month to each row by clone.id
X.lc <- X.lc[
    X.lc$month <= X.lc$first.censor.month, # keep rows up to first artificial censoring month for each clone
]
X.lc$censor.event <- as.integer(
    X.lc$month == X.lc$first.censor.month & is.finite(X.lc$first.censor.month)
)
X.lc$first.censor.month <- NULL


# merge baseline covariates into long table format for modeling 
baseline.vars <- c(
    "patient.id", "age", "sex", "race", "ecog", "smoking.status", "practice.type", "insurance.type", 
    "initial.stage", "cancer.type", "pd.l1", "histology", "treatment.type", "imdc.risk.score"
)
X.long  <- merge(X.long, X[ , baseline.vars], by=c("patient.id", "cancer.type"), all.x = TRUE, sort = FALSE)
X.lc    <- merge(X.lc, X[ , baseline.vars], by=c("patient.id", "cancer.type"), all.x = TRUE, sort = FALSE)




# data quality checks 
stopifnot(nrow(X) == n, !anyDuplicated(X$patient.id))
stopifnot(all(ifelse(X$cancer.type == "NSCLC", !is.na(X$histology), is.na(X$histology))))
stopifnot(all(ifelse(X$cancer.type == "Kidney", !is.na(X$imdc.risk.score), is.na(X$imdc.risk.score))))
stopifnot(all(X.long$any.event %in% c(0, 1)))
stopifnot(all(X.lc$censor.event %in% c(0, 1)))




# ---------- Quick diagnostics ----------
cat("Patient-level cohort X:\n")
cat("  n patients:", nrow(X), "\n")
cat("  cancer counts:\n")
print(table(X$cancer.type))
cat("\nLongitudinal cohort X.long:\n")
cat("  n rows:", nrow(X.long), "\n")
cat("  first event type distribution:\n")
print(table(X$event.type))
cat("\nCloned longitudinal cohort X.lc:\n")
cat("  n rows:", nrow(X.lc), "\n")
cat("  artificial censor events by strategy:\n")
print(with(X.lc, table(assigned.strategy, censor.event)))


invisible(list(X = X, X.long = X.long, X.lc = X.lc))
saveRDS(list(X = X, X.long = X.long, X.lc = X.lc), file = file.path("grace-period-model", "gp-data.rds"))
