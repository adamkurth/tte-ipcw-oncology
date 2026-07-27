# # Target Trial Emulation (TTE) using Inverse Probability Clone Censoring Weighting (IPCW) and Regularized Regression of Flatiron database for advanced NSCLC patients treated with immunotherapy

################################# 
rm(list=ls()); set.seed(20260717)

library(ggplot2); library(dplyr) # visual
library(glmnet); library(splines); library(survival) # modeling
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


# ---------- calculate IPCW ----------
# model probability of artificial censoring (censor.event=1)
# use natural splines (ns) for 'month' to allow for non-linear effects of months 24-27 hazard on censoring probability



