# Target Trial Emulation (TTE) using Inverse Probability Clone Censoring Weighting (IPCW) and Regularized Regression of Flatiron database for advanced NSCLC patients treated with immunotherapy


## Overview 
 - objective: basic risk prediction using regularized regression for outcomes 1. overall survival probability, 2. irAEs
 - covariates: patient demographics, dates of advanced cancer diagnosis, treatment dates, type of systemic therapy, ECOG performance status, comorbidities, smoking history, laboratory values (liver/renal function, tumor biomarkers), medications, and validated endpoints (mortality, real-world progression, select irAEs). 
 - adults of 18+, diagnosed with advanced bladder, kidney, melanoma, or non-small cell lung cancer (NSCLC) between 2015-2024. 
 - results stratified by: cancer type, age group, BMI, ECOG performance status, programmed death ligand (PD-L1) status, tumor mutational burden (TMB), and calendar period of ICI initiation.
  - sensitivity analysis to analyze the 2-year time-point decision of discontinuation/continuation of ICI therapy to earlier time points (e.g., 18 months).
- estimand: primary estimand is the per-protocol effect of discontinuation vs. continuation on overall survival (primary), and irAEs (secondary) at 2 years after ICI initiation.
- eligibility criteria: free from evidence of progression at start of 21st month after treatment initiation. 
- per-protocol TTE: discontinuation anywhere from 21-27th month after treatment initiation is considered as consistent with "discontinuation after 2 years" (per-protocol). 
- artificial censoring: patients assigned to "discontinue" strategy will be considered to deviate from protocol/strategy if they continue treatment (not discontinue) by 27 months. Likewise, patients assigned to "continue" strategy will be censored if discontinuation occurs prior to the end of 27 months. 
- power/sample size: estimate 7,837 patients (5,980 lung, 624 bladder, 589 kidney, and 644 melanoma) that were progression-free and initiated ICI therapy at 2-years. 2,405 discontinued at 2-year decision point, 5,432 continued ICI therapy beyond 2-years. 
- cloning: each patient is cloned at month 21, and given an assigned protocol of either discontinuation or continuation of ICI therapy. If they deviate from their assigned protocol, they are artificially censored at the time of deviation. This allows for estimation of the per-protocol effect of discontinuation vs. continuation on overall survival and irAEs at 2 years after ICI initiation.
- administrative censoring occurs exactly at 48 months after ICI initiation, clinical decision anchored at 2-year decision point, and follow-up is limited to 2-years after the 2-year decision point (i.e., 4-years after ICI initiation). In follow-up period, the study aims to estimate the difference in overall survival and irAEs over the 2-year follow-up period after the 2-year decision point. Thus, any patient who survives and remains event free at the end of the two-year follow-up period will be administratively censored at 48 months after ICI initiation.
- visual: Cloned at the start of follow-up period (21 months), and at cloning, progression-free patients is are cloned into two identical profiles: one "continue" and one "discontinue" ICI therapy. Deviating from "discontinue" protocol means that if a patient does not stop treatment by 27 months, their "discontinue" clone is censored at month 27. Likewise, if a patient deviates from "continue" protocol by discontinuing treatment prior to month 27, their "continue" clone is censored at the time of discontinuation. 

![Visualization of TTE/IPCW](figures/tte_swimmer_plot.png)
- PT1: (Continued, Died). Patient never stops treatment. Their "discontinue" clone is artificially censored exactly at month 27  (end of grace period), while their "continue" is tracked until death at month 42. 
- PT2: (Discontinued, Survived). Patient decides to stop treatment at month 24. Their "continue" clone is artificially censored at month 24. Their "discontinue" clone is tracked until end of trial (administrative censoring at month 48).
- PT3: (Discontinued Early, Died). Similar to PT2, but they stop at month 22. Their "continue" clone is artificially censored at month 22, and "discontinue" clone captures death at month 30.
- PT4: (Continued, Survived). Best-case scenario. Their "discontinue" clone is artificially censored at month 27, and their "continue" clone is successfully reaches the end of trial (48 months).
- PT5: (Died during Grace Period). A crucial edge case in TTE, since they died at month 25 BEFORE making any decision to continue or discontinue, both clones are still valid at the time of death. Both clones register the death event at month 25. 

## Aim 1 
### MODELING STRATEGY
- Aim 1, step 1: 1) descriptive analyses of patient characteristics, treatment patterns (freq and timing of treatment discontinuation), and length of available follow-up; 2) Baseline covariates summarized at start of 21st month decision window; 3) will examine distribution of discontinuation, progression, and irAEs. 4) data quality checks, missingness, internal consistency of key dates and events, plausiblility checks to ensure validity of analytic variables to inform model specification.
 - Aim 1, step 2: 1) estimate the per-protocol effect of discontinuation at 2-years vs continued ICI treatment until progression or irAE on primary outcome of overall survival using pooled logistic regression with IPCW to account for artificial censoring due to deviation from assigned strategy. Confounding eliminated due to all patients assigned both treatment strategies at baseline (cloning), and adjustment for confounding is incorperated into models for the censoring process rather than treatment assignment or outcome models.
 - fitted logistic model will estimate: 1) standardized survival curves under each treatment strategy and use these to calculate the primary estimand of interest: the difference in 2-year survival probability between the two treatment strategies. 2) estimate the restricted mean survival time (RMST) over the 2-year follow-up period. Noninferiority (margin of 5%) will be assessed by the difference in overall survival between continuation vs. discontinuation at 2-years. 
- EDA: 1) varying discontinuation timepoints to explore discontinuation at 18 months or 1 year; 2) stratified subgroup analyses will evaluate heterogeneity of treatment effects (HTE) across factors including cancer type, age group, BMI, ECOG performance status, PD-L1, and/or TMB status, and calendar period of ICI initiation by including interaction terms between treatment strategy and HTE factor of interest in the pooled logistic regression model. 

### COVARIATES & CONFOUNDERS 
- covariates and confounders: 1) common to all cancer types: age, sex, race, ethnicity, ECOG performance status, smoking status, practice type (academic vs. community), insurance type, and initial stage of diagnosis (all may influence both treatment and outcomes). Seperate IPCW models will be fit for each cancer type to allow for confounders unique to each cancer type, including cancer-specific confounders to account for differences in disease biology and treatment selection. 

Cancer-specific confounders:

 - NSCLC specific confounders: PD-L1 expression level, treatment type (ICI alone or ICI + chemotherapy), histology (squamous vs. non-squamous), key genetic alterations (e.g., KRAS, STK11, EGFR). 
 - bladder cancer specific confounders: treatment regimen (ICI alone or ICI + enfortumab vedotin). 
 - kidney cancer specific confounders: IMDC risk score  (poor, intermediate, favorable), treatment type (ICI alone, ICI doublet, or ICI + targeted therapy). 
 - other variables will be included based on EDA, including tumor biomarkers TMB, microsatellite instability status, or area-level socioeconomic status (Yost index), and site of metastasis (e.g., central nervous system involvement).

## Aim 2
- Aim 2, step 1: 1) examine composite outcome of time to first occurrence of any irAE under study (colitis, pneumonitis, hypothyroidism, adrenal insufficiency, or hypophysitis); 2) followed by analyses of time to first occurrence of specific irAE types. Since death preludes the occurrence of irAEs, it is treated as a competing event; 3) use descriptive analyses to summarize the frequency and timing of irAEs events, including variation by cancer types, and patient subgroups. 
- Aim 2, step 2: 1) estimate the effect of treatment strategy on irAE outcomes using pooled logistic with IPCW (structured in discrete time intervals). 2) Competing risks will be addressed by modeling the cause-specific hazard of irAE and death, using these models to estimate standardized cumulative incidence functions under each treatment strategy. 3) the primary estimand will be the difference in cumulative incidence of irAE at 2-years. 

## Aim 3
- develop individualized risk prediction models for each of the two outcomes of interest (overall survival and irAEs) using regularized regression (elastic net penalty to stabilize estimation with tuning parameters selected via cross-validation). Models will include: 1) binary indicator for discontinuation vs. continuation of ICI therapy, 2) patient and tumor characteristics available at that timepoint, interaction between continuation/discontinuation and these characteristics, 3) characteristics include " cancer type, age, sex, race, ethnicity, ECOG performance status, smoking status, practice setting (academic versus community), insurance type, initial stage at diagnosis, PD- L1 expression level, treatment type (ICI alone or in combination [another ICI, chemotherapy, or targeted therapy]), histology, genomic alterations (e.g., KRAS, STK11), IMDC risk score, TMB, microsatellite instability status, Yost index, and site of metastasis." 4) for risk factors only relevant for subset of cancer types (e.g., IMDC risk score, PD-L1 expression level) we will include interactions between the risk factors and relevant cancer types to allow these variables to contribute to risk estimation only for the relevant cancer types. 5) This will allow estimation of individualized predicted risks of overall survival and irAEs under each treatment strategy conditional on patient constellation of clinical and tumor characteristics. 6) Given the limited sample size for some characteristics, the elastic net penalty allows for stabilized parameter estimates borrowing information across patient subgroups. 7) Model performance will be based on optimism-corrected internal-validation using bootstrap resampling evaluating discrimination and calibration. Performance evaluated based on bootstrap-corrected time-varying c-statistic, expected-to-observed ratio, calibration slope, and calibration intercept — all assessed at 2-years after 2-year decision point.
- Aim 3, step 2: develop RShiny application. 1) adapting the final models into backend functions which will be assessed by RShiny to calculate overall survival probabilities and cumulative incidence of irAEs at 1-year and 2-years after the 2-year decision point. 
- Aim 3, step 3: Missingness. 1) we expect some confounder missingness including ECOG performance status. Missingness will be addressed using multiple imputation under a missingness at random assumption (MAR), with imputation models including all covariates used in the analysis as well as outcomes and indicators for artificial censoring. 2) Sensitivity analyses using tipping point approach to assess robustness of departures from MAR assumption. 


## Modeling (Grace Period Model): Regularized IPCW 

In the clone-censor-weighting approach, clinical decisions naturally introduce a risk of selection bias (and immortal time bias). For instance, a patient's clinical status may worsen (developing severe irAE at month 23), their clinician will most likely stop treatment, causing artificial censoring for their "continue" clone by stopping treatment. Consequently, the patient remaining uncensored in the "continue" arm are artificially healthier than average baseline population because those prone to irAE or progression have been censored and selectively removed. If we model the survival data directly on the remaining patients, our results will be hopelessly biased in favor of the "continue" arm—a form of selection bias under informative censoring. 


To fix this, we use inverse probability of censoring weights (IPCW) to create a pseudo-population. The intuition is: if an unhealthy patient is artificially censored, we find another patient who looks similar to them (same age, sex, ECOG, etc.) but was not censored and we up-weight that uncensored patient, it appears as though nobody was censored because the remaining patients are now representative of the original cohort. 


$$
  \tilde{W}_{i,t} = \prod_{k=0}^t \frac{ P(C_{i,k} = 0 \mid \bar{V}_{i}, A_i, \bar{C}_{i,k-1} = 0) }{ P(C_{i,k} = 0 \mid \bar{X}_{i,k}, A_i, \bar{C}_{i,k-1} = 0) } \qquad \text{stabilized IPCW for patient $i$ at time $t$}
$$

Where $i=1, \ldots, n$ indexes patients, $t=0, \ldots, T = 48$ indexes discrete time intervals (months of follow-up), $\bar{C}_{i,t}$ is the full patient censoring history, $\bar{V}_{i}$ is the vector of baseline covariates, and $\bar{X}_{i,t}$ is the covariate history encompassing the entire longitudinal trajectory of their baseline traits and time-varying clinical status (irAEs labs, performance status) up to month $k$. 


The denominator must contain all information related to censoring, namely time-varying information in $\bar{X}_{i,k}$ that is predictive of censoring. Specifically, $\bar{X}_{i,k}$ provides information about the patient's evolving clinical status (prognostic variables), and thus their likelihood of censorship from the study. The denominator model includes $\bar{X}_{i,k}$, in addition to treatment assignment $A_i$ and the censoring history $\bar{C}_{i,k-1}$.  This captures the true probability of censorship given the evolving clinical history. If we fail to include treatment assignment in the denominator, we would fail to account for the fact that treatment assignment is predictive of censorship. Additionally, to make the weights stable, we include baseline covariates $\bar{V}_i$ in the numerator model to capture a coarse summary of the patient's baseline characteristics. Taking the ratio of the two probabilities gives us the stabilized weight that controls variance while adjusting for time-dependent confounding.


Once every patient has their stabilized IPCW weights $\tilde{W}_{i,t}$ calculated, we can estimate the causal effect of discontinuing versus continuing treatment. TTE handles this via a discrete person-period data structure, where each patient contributes a distinct row for every month they remain at risk. Because target events (e.g., death or incidence of a specific irAE) are relatively rare within any single month, a pooled logistic regression model closely approximates a time-dependent Cox-PH model. Rather than estimating a continuous hazard ratio, we predict the discrete-time-hazard (binary probability of event occurrence) in a given month. We model time flexibility using restricted cubic splines for the month variable because the hazard is unlikely to be linear over time. The pooled logistic model is fit using the stabilized IPCW weights, $\tilde{W}_{i,t}$, to account for the artificial censoring introduced by the TTE design. 

In oncology cohorts, the covariate history $\bar{X}_{i,k}$ is often likely extremely high-dimensional. To avoid the multicollinearity and overfitting from including too many covariates, the elastic net penalty stabilizes the estimation by shrinkage and variable selection. The penalty term is a convex combination of the LASSO and ridge penalties:
$$
  \lambda \left( \frac{ 1-\alpha }{2} \underbrace{ \| \beta \|_2^2 }_{\text{Ridge}} + \alpha \underbrace{ \| \beta \|_1 }_{\text{LASSO}} \right)
$$
The $L_1$ norm (LASSO) performs variable selection by shrinking some coefficients to zero, while the $L_2$ norm (ridge) handles multicollinearity by shrinking correlated coefficients together. 


### Why do we not penalize treatment or time variables? 
[Hahn et. al. (2018)](https://projecteuclid.org/journals/bayesian-analysis/volume-13/issue-1/Regularization-and-Confounding-in-Linear-Regression-for-Treatment-Effect-Estimation/10.1214/16-BA1044.full) states that naive regularization can corrupt causal estimation and yield exceptionally poor estimators, even using a perfectly specified model. 

$$Y_i = \alpha Z_i + X_i \beta + \nu_i$$ 

Where in our framework, $Y_i$ is the discrete-time hazard of the event, $Z_i$ represents our primary structural variables (`assigned.strategy`, `on.ici`, and time spline month variables), $\alpha$ is the causal effect of interest, $X_i$ is the high dimensional vector of nuisance covariates, and $\beta$ represents the control effects of these nuisance covariates. The error term $\nu_i$ is assumed to be independent of the covariates $X_i$ and $Z_i$. Under exogeneity conditions for valid causal interpretations (for $\alpha$) the model must satisfy:

$$\mathrm{Cov}(Z_i, \nu_i \mid X_i) = 0$$

Ensuring that, conditional on the controls, treatment assignment is independent of the error term. Regularization optimizes for prediction, not causal estimation, which violates this condition in two ways: 
  
  1)  Direct shrinkage of causal effect: If we allow the elastic net to penalize the treatment variable $Z_i$ (previously denoted $A_i$) or the baseline hazard, the model applies a shrinkage prior directly to $\alpha$. This shrinks the causal effect estimate towards zero in exchange for a reduction in predictive variance.
  2)  Regularization-induced confounding: Even if we attempt to isolate $\alpha$, applying shrinkage to the high-dimensional controls $X_i$ corrupts the estimation of $\beta$. This is because the prior "prefers" small coefficients for the controls ($\beta$), the model will often achieve a similar in-sample fit by over-shrinking the control coefficients and shift explanatory weight to adversely bias the treatment effect $\alpha$. Hahn et. al formally demonstrate this by defining the bias for naively regularized regression ($\hat{\alpha}_{rr}$) as: $$\mathrm{bias}(\hat{\alpha}_{rr}) = - (z^\top z)^{-1} z^\top X \left( I_p + X^\top (X - \hat{X}_Z) \right)^{-1} \beta$$ where $\hat{X}_Z$ is the projection of $X$ onto $Z$ (fitted values). This formula here shows that the bias of treatment effect is a function of every element in the true unknown vector $\beta$, and the stronger the confounding in the data the worse the bias on treatment effect estimation $\hat{\alpha}$ will be. In other words, the bias is a function of the correlation between $Z$ and $X$, and the magnitude of $\beta$.

To avoid regularization-induced confounding Hahn et. al. propose reparameterizing the model into two-equation systems that isolate the treatment *selection* mechansim from the *response* mechanism. In our design we achieve this separation by using IPCW to resolve the selection mechanism by balancing the pseudo-population, which allows us to safely partition the response model. 

To resolve the shrinkage bias, we partition the model using penalty factors. Assigning penalty factors of 0 to structural variables ($Z_i$: strategy, month time, and time-varying treatment), we force the model to estimate their causal relationships $\alpha$ without shrinkage. Simultaneously, the model focuses its entire $L_1 / L_2$ penalty on the high-dimensional nuisance covariates $X_i$ to adjust for residual confounding without overfitting. 


## Diagnostics
### Numerator and Denominator Models

![Diagnostic: Numerator and Denominator Forest Plot](figures/gp-model/tte_forest_plot.png)
This plot shows that all the baseline covariates sit squarely on the null line. However, the `on.ici` is pushed far to the left, **indicating that patients who are on ICI therapy are much less likely to be censored** (i.e., more likely to continue treatment). In other words, this confirms that actively taking the drug significantly protects a patient from being artificially censored in the TTE design. 

Additionally, the `assigned.strategy` variable (discontinue 2yr) is pushed far to the right, **indicating that patients assigned to "discontinue at 2-years" are more likely to deviate from assigned strategy and be artificially censored, and continue treatment**. This is also consistent, since patients assigned to discontinue at 2-years are more likely to continue treatment (and thus be censored) than those assigned to continue treatment.


### Spline Visualization

![Diagnostic: Spline Visual](figures/gp-model/tte_spline_censoring.png)

In the gray shaded region (grace period), the predicted probability of censoring begins to climb and peaks at month 27. This indicates that the TTE design enforces month 27 as the end of the grace period, triggering lots of artificial censoring to occur for "discontinue" clones that fail to stop treatment by month 27. After this, the probability of censoring drops off near zero by month 30, indicating that the decision window has closed. 


### IPCW Distribution

![Diagnostic: Truncated IPCW Distribution](figures/gp-model/tte_ipcw_distribution.png)


The stabilized IPCW distribution is centered around 1 representing the null value of no censoring. While there is a slight right tail, the maximum truncated weights are approximately 1.5. The narrow range from ~ 0.9 to 1.5 indicates that the IPCW weights are highly stable with no extreme outliers that would threaten estimation / require aggressive truncation. This is a good sign that the IPCW model is well-specified and that the stabilized weights are not overly influenced by a small number of patients with extreme weights.




## Elastic Net CV Curve & Coefficient Selection

![Diagnostic: Elastic Net CV Curve](figures/gp-model/tte_elastic_net_cv_curve.png)

The cross-validation plot visualizes the tradeoff between complexity and performance. The y-axis is the binomial deviance (lower is better), while the x-axis shows the $-\log(\lambda)$ penalty, as you move from left to right the penalty decreases allowing more variables in the model. The numbers on the top indicate the actual amount of nonzero coefficients in the model at that penalty. 

- What is the "1SE" rule? Using CV, $CV(\lambda) = 1/K \sum_k E_k(\lambda)$ we find the mean squared error across the K folds, and identify the minimum error $\lambda_{min} = \arg\min_\lambda CV(\lambda)$, with the error denoted $CV_{min} = CV(\lambda_{min})$. This is the best performing model, but may be overfit (vertical dashed line on right side). To avoid overfitting, we first compute the standard error from sample variance ($s^2(\lambda_{min}) = \frac{1}{K-1} \sum_k (E_k(\lambda_{min}) - CV_{min})^2$), then scale accordingly to get the standard error $SE(\lambda_{min}) = \sqrt{\frac{ s^2(\lambda_{min})} { K }}$. The 1-SE rule, balancing performance and parsimony, creates a threshold for the simplest model within 1 standard error of the minimum error, $\lambda_{1se} = \max_\lambda \{ \lambda : CV(\lambda) \leq CV_{min} + SE(\lambda_{min}) \}$ (vertical dashed line on left side). This is the model we will use for prediction, as it is more parsimonious and less likely to overfit.

By extracting the coefficients at $\lambda_{1se}$ while utilizing our non-zero penalization factors, we observe how the model behaves. The unpenalized strucutural variables (`assigned.strategy`, `on.ici`, and month spline variables) are retained in the model by force. **Notably, `on.ici` has strong negative coefficient (-.162, OR = 0.85), indicating that active treatment is protective against the target event (death or irAE)**. 

From the high-dimensional covariate pool, the model selected only `age` and geometric contrasts of ECOG performance (`.L` and `.Q`). Age shows a mild positive association with the event risk (OR = 1.016). For ECOG, the linear contrast (`.L` = 0.096) indicating that performance status worsens, event risk increases. But the quadratic contrast (`.Q` = -0.102)  reveals a concave down risk curve, meaning that risk increases at a decreasing rate as ECOG increases.


## Aim 3: Standardized Marginal Survival Curves (G-Computation)

With the weighted, regularized, and fitted model, we translate this into a counterfactual survival curve for each treatment strategy using g-computation. We take the baseline cohort and clone it into two counterfactual worlds: A) everyone perfectly adheres to the "continue" strategy, and B) everyone perfectly adheres to the "discontinue" strategy. WE predict the month probability for every patient in both worlds, compute their cumulative event-free survival over time, and take the marginal average across the entire cohort. 

![Standardized Marginal Survival Curves](figures/gp-model/tte_marginal_survival.png)

The resulting standardized survival curves shows a clear divergence after the 24 month mark. At teh end of the 48 month follow-up window, the estimated survival for "continue" is 46.2% vs. 40.6% for "discontinue". This yields a risk difference (RD) (recall, non-collapsible) of 5.6%. In clinical terms, continuing ICI therapy for an additional 2 years yields an absolute survival benefit of roughly 5.6% over discontinuation, properly adjusted for both baseline and time-varying confounders. 





--- 
# Change TTE Design: When is the optimal time to discontinue ICI therapy?

Instead of "continue" or "discontinue" at 2-years, now we can explore when the best time to discontinue ICI therapy is. This extends the two-arm target trial into a multi-arm design by expanding the treatment strategy to include multiple timepoints for discontinuation (e.g., 3, 6, 9 and 12 months). 




##  Modeling (Treatment Discontinuation Model): 










