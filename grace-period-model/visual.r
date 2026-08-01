
################################# 
rm(list=ls()); set.seed(20260717)

library(ggplot2); library(dplyr) # visual
source("sim-data.r")

# ---------- visualization schematic for TTE -------------



# visualization of TTE
tte.viz.data <- data.frame(
    patient.id = rep(c("PT1", "PT2", "PT3", "PT4", "PT5"), each = 2),
    assigned.strategy = rep(c("Continue", "Discontinue"), 5),
    # everybody starts follow-up at month 21
    start.time = rep(21, 10),
    # define end times based on artificial censoring or actual events
    end.time = c(
        42, 27,     #PT1: continue past grace period. Died at 42m
        24, 48,     #PT2: discontinue at 24m. Survived to end of trial
        22, 30,     #PT3: discontinue early at month 22m. Died at 30m
        48, 27,     #PT4: continue past grace period. Survived to end of trial
        25, 25      #PT5: died at 25m (during grace period) BEFORE discontinuation decision
    ),
    event.type = c(
        "Death", "Artificial Censoring",
        "Artificial Censoring", "Administrative Censoring",
        "Artificial Censoring", "Death",
        "Administrative Censoring", "Artificial Censoring",
        "Death", "Death"
    ),
    y.position = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1) # space out on Y-axis
)
max.y <- max(tte.viz.data$y.position)

# swimmer plot
ggplot(tte.viz.data) +
  # Shaded region for "grace period" (Months 21 - 27)
  annotate("rect", xmin = 21, xmax = 27, ymin = 0.5, ymax = max.y + 0.5, 
           alpha = 0.15, fill = "black") +
  
  # Draw patient trajectory lines
  geom_segment(aes(x = start.time, xend = end.time, 
                   y = y.position, yend = y.position, 
                   color = assigned.strategy), 
               linewidth = 1.5) +
  
  # Event points at the end of the trajectories
  geom_point(aes(x = end.time, y = y.position, shape = event.type, fill = event.type), 
             size = 3.5, color = "black") +
  
  # Vertical lines for milestones
  geom_vline(xintercept = 21, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 48, linetype = "solid", color = "darkred", linewidth = 1) +
  
  # Annotations on timeline
  annotate("text", x = 21, y = max.y + 1, label = "Cloning\n(Month 21)", hjust = 0.5, size = 3.5, fontface = "bold") +
  annotate("text", x = 27, y = max.y + 1, label = "End Grace Period\n(Month 27)", hjust = 0.5, size = 3.5, fontface = "bold") +
  annotate("text", x = 48, y = max.y + 1, label = "End Trial\n(Month 48)", hjust = 1, size = 3.5, fontface = "bold") +
  
  # Formatting Y-axis to label patient archetypes clearly
  scale_y_continuous(
    breaks = c(1.5, 3.5, 5.5, 7.5, 9.5), 
    labels = c(
      "PT5: Died during\nGrace Period (25m)",
      "PT4: Continued tx,\nAdmin Censored (48m)",
      "PT3: Discontinued early (22m),\nDied (30m)",
      "PT2: Discontinued (24m),\nAdmin Censored (48m)", 
      "PT1: Continued tx,\nDied (42m)"
    )
  ) +
  
  # Manual scales for shapes, fills, and colors
  scale_shape_manual(values = c("Artificial Censoring" = 4, 
                                "Administrative Censoring" = 21, 
                                "Death" = 24)) +
  scale_fill_manual(values = c("Artificial Censoring" = "black", 
                               "Administrative Censoring" = "white", 
                               "Death" = "red")) +
  scale_color_manual(values = c("Continue" = "#1f78b4", "Discontinue" = "#fb9a99")) +
  
  # Theming
  theme_minimal(base_size = 12) +
  labs(
    title = "Target Trial Emulation: Clone-Censor Trajectories by Patient Archetype",
    subtitle = "Tracking clone outcomes through the 21-27 month grace period",
    x = "Months Since First-Line Treatment Initiation",
    y = "",
    color = "Assigned Clone Strategy",
    shape = "Event Type",
    fill = "Event Type"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "gray30")
  ) +
  coord_cartesian(xlim = c(20, 50), ylim = c(0.5, max.y + 1.5))

ggsave(file.path("figures", "tte_swimmer_plot.png"), width = 10, height = 6, dpi = 300)



# ---------- visual of num/denom models ----------

# run model.r until m.num and m.denom are created, then run the following code to visualize the models
est <- coef(m.denom)
se <- summary(m.denom)$coefficients[, "Std. Error"]

or.data <- data.frame(
  Variable = names(est),
  OR = exp(est),
  Lower = exp(est - 1.96 * se),
  Upper = exp(est + 1.96 * se),
  stringsAsFactors = FALSE
)
# filter out intercept and sline terms to focus on covariates
keep.idx <- !grepl("(Intercept)|ns\\(month", or.data$Variable)
or.data <- or.data[keep.idx, ]

# clean up names perserve order on plot
or.data$Variable <- factor(or.data$Variable, levels = rev(or.data$Variable))

p.forest <- 
ggplot(or.data, aes(x = OR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 1, color = "#3B7EA1") +
  geom_point(size = 3, color = "black") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 1.5, 2, 4)) + # Log scale for Odds Ratios
  theme_minimal(base_size = 12) +
  labs(
    title = "Forest Plot: Predictors of Artificial Censoring (Denominator Model)",
    subtitle = "Odds Ratios (Log Scale) with 95% Confidence Intervals",
    x = "Odds Ratio (OR)",
    y = ""
  ) +
  theme(panel.grid.minor = element_blank())

print(p.forest)
ggsave(file.path("figures", "tte_forest_plot.png"), width = 10, height = 6, dpi = 300)
# visual 




# ---------- visual of num/denom models ----------

# run model.r until m.num and m.denom are created, then run the following code to visualize the models
est <- coef(m.denom)
se <- summary(m.denom)$coefficients[, "Std. Error"]

or.data <- data.frame(
  Variable = names(est),
  OR = exp(est),
  Lower = exp(est - 1.96 * se),
  Upper = exp(est + 1.96 * se),
  stringsAsFactors = FALSE
)
# filter out intercept and sline terms to focus on covariates
keep.idx <- !grepl("(Intercept)|ns\\(month", or.data$Variable)
or.data <- or.data[keep.idx, ]

# clean up names perserve order on plot
or.data$Variable <- factor(or.data$Variable, levels = rev(or.data$Variable))

p.forest <- 
ggplot(or.data, aes(x = OR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 1, color = "#3B7EA1") +
  geom_point(size = 3, color = "black") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 1.5, 2, 4)) + # Log scale for Odds Ratios
  theme_minimal(base_size = 12) +
  labs(
    title = "Forest Plot: Predictors of Artificial Censoring (Denominator Model)",
    subtitle = "Odds Ratios (Log Scale) with 95% Confidence Intervals",
    x = "Odds Ratio (OR)",
    y = ""
  ) +
  theme(panel.grid.minor = element_blank())

print(p.forest)
ggsave(file.path("figures", "tte_forest_plot.png"), width = 10, height = 6, dpi = 300)
# ---------- visual of spline trajectories (censoring over time) ----------
# uses numerator model to isolate the effect of grace period

# dummy dataset spanning follow-up period (months 21-48)
months.seq <- 21:48
dummy.data <- data.frame(
  month = rep(months.seq, 2),
  assigned.strategy = factor(
    rep(c("cont.2yr", "disc.2yr"), each = length(months.seq)),
    levels = c("cont.2yr", "disc.2yr")
  )
)
# predict prob. of censoring for each month and both strategies
dummy.data$prob.censor <- predict(m.num, newdata = dummy.data, type = "response")

p.spline <- ggplot(dummy.data, aes(x = month, y = prob.censor, color = assigned.strategy)) +
  # Add shaded region for the decision grace period (Months 21 to 27)
  annotate("rect", xmin = 21, xmax = 27, ymin = 0, ymax = max(dummy.data$prob.censor) + 0.05, 
           alpha = 0.1, fill = "black") +
  
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1) +
  
  scale_color_manual(
    values = c("cont.2yr" = "#1f78b4", "disc.2yr" = "#fb9a99"),
    labels = c("Continue Clone", "Discontinue Clone")
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Probability of Artificial Censoring Over Time (Spline Fit)",
    subtitle = "Shaded area represents the 21-27 month grace period",
    x = "Month of Follow-up",
    y = "Predicted Probability of Censoring",
    color = "Assigned Clone Strategy"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p.spline)
ggsave(file.path("figures", "tte_spline_censoring.png"), width = 10, height = 6, dpi = 300)


# ----------------- visualization of IPCW weights -----------------
# create ipcw object in model.r and then run the following code
p.weights <- ggplot(X.lc, aes(x = ipcw)) +
  geom_histogram(bins = 40, fill = "#3B7EA1", color = "white") +  
  geom_vline(xintercept = 1, linetype="dashed", color = "darkred", linewidth=1) + 
  theme_minimal() +
  labs(title = "Diagnostic: Truncated IPCW Distribution", x = "IPCW", y = "Person-Months")
print(p.weights)

