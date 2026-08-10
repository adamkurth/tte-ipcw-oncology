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


---
# Modeling ("Grace Period Model"): Regularized IPCW 

In the clone-censor-weighting approach, clinical decisions naturally introduce a risk of selection bias (and immortal time bias). For instance, a patient's clinical status may worsen (developing severe irAE at month 23), their clinician will most likely stop treatment, causing artificial censoring for their "continue" clone by stopping treatment. Consequently, the patient remaining uncensored in the "continue" arm are artificially healthier than average baseline population because those prone to irAE or progression have been censored and selectively removed. If we model the survival data directly on the remaining patients, our results will be hopelessly biased in favor of the "continue" arm—a form of selection bias under informative censoring. 

### 1. Clones, Artificial Censoring and IPCW
First, we must create the "clones". These are copies of the patients in the cohort, where each clone is assigned to a specific treatment strategy, here for one treatment continuation decision at month 21, $a \in \{ \text{discontinue, continue}\}$. A clone is artificially censored the moment its observed data becomes inconsistent with the assigned treatment strategy. Since this censoring is deterministic, it's also informative to the censoring process —— and must be corrected through stabilized inverse probability of censoring weights (IPCW). 

The IPCW approach uses two pooled logistic models to estimate the per-month hazard of artifical censoring $C_{i,t}$: 
$$
  \text{numerator: } \quad \mathrm{logit} [P(C_{i,t} = 0 \mid \bar{V}_{i}, A_i, \bar{C}_{i,t-1} = 0)]
      = \beta_0 + \beta_1 A_i + f_{\text{ns}}(t)
$$
$$
  \text{denominator: } \quad  \mathrm{logit} [P(C_{i,t} = 0 \mid \bar{X}_{i,t}, A_i, \bar{C}_{i,t-1} = 0)]
      = \gamma_0 + \gamma_1 A_i + f_{\text{ns}}(t) + \gamma^\top \bar{X}_{i,t}
$$
for $i=1, \ldots, n$ (patients) and $t = 21, \ldots, T = 48$ indexes the discrete-time calendar months since ICI initiation, where $t=21$ is the cloning/decision timepoint and $T=48$ is the administrative censoring month. $C_{i,t}$ is the month-specific indicator of whether a patient was censored, and $\bar{C}_{i,t}$ is the patient's censoring history through month $t$ (e.g., $\bar{C}_{i,T} = ( C_{i,21}=0, \dots, C_{i,t} = 0, \dots,  C_{i,T} = 0)$ for somebody who remains uncensored throughout). The baseline confounders $\bar{V}_{i} = (\text{age}, \text{sex}, \text{race}, \text{ECOG}, \text{cancer type}, \text{practice type}, \dots)$ are time-invariant, as well as cancer-specific baseline covariates.

Two distinct treatment-related trajectories must be distinguished: the **observed** treatment history $\bar{A}_i = (A_{i,21}, \dots, A_{i,t})$, where $A_{i,k} \in \{0,1\}$ ("`on.ici`") indicates whether patient $i$ was actively on ICI therapy at month $k$ — this is the actual, realized trajectory taken by the patient/clone; and the **assigned (target) regime** $\bar{A}^\ast_i = (A^\ast_{i,21}, \dots, A^\ast_{i,T})$, the deterministic treatment-status trajectory implied by the assigned strategy $a$ as a function of time (e.g., $A^\ast_{i,k}=1$ for all $k$ under "continue"; $A^\ast_{i,k} = \mathbb{1}\{k < 27\}$ under "discontinue"). A clone is artificially censored the first month its observed history $\bar{A}_i$ deviates from its assigned regime $\bar{A}^\ast_i$. Time-varying indicators used in the censoring models are collected in $\bar{X}_{i,t} = (\bar{A}_{i,t}, \bar{V}_i)$, combining the observed on-treatment history with the baseline confounders. $A_{i,t}$ (the current entry of $\bar{A}_i$) and the $f_{\text{ns}}(t)$ natural spline (df=3) for months of follow-up are included in both the numerator and denominator models to account for nonlinear baseline hazards in time.

Note on time-varying confounders in the denominator model: The denominator must contain all information related to censoring, namely time-varying information in $\bar{X}_{i,t}$ that is predictive of censoring. Specifically, $\bar{X}_{i,t}$ provides information about the patient's evolving clinical status (prognostic variables), and thus their likelihood of censorship from the study. The denominator model includes $\bar{X}_{i,t}$, in addition to treatment assignment $A_i$ and the censoring history $\bar{C}_{i,t-1}$. This captures the true probability of censorship given the evolving clinical history. If we fail to include treatment assignment in the denominator, we would fail to account for the fact that treatment assignment is predictive of censorship. Additionally, to make the weights stable, we include baseline covariates $\bar{V}_i$ in the numerator model to capture a coarse summary of the patient's baseline characteristics. Taking the ratio of the two probabilities gives us the stabilized weight that controls variance while adjusting for time-dependent confounding.

### 2. Cumulative Probability of Remaining Uncensored
Next, we take the cumulative discrete-survival product of those remaining uncensored as the "trial" continues. For each clone, ordered chronologically by month, we calculate the cumulative probability of remaining uncensored up to month $t$ for each patient $i$:
$$
  \widehat{\text{cum.num}}_{i,t} = \prod_{k \leq t} \left[ 1 - \hat{P}_{\text{num}}(C_{i,k} = 1)\right]
    = \prod_{k = 21}^{t} \left[ 1 - \hat{P}(C_{i,k} = 1 \mid \bar{V}_{i}, A_i, \bar{C}_{i,k-1} = 0)\right]
$$
$$
  \widehat{\text{cum.denom}}_{i,t} = \prod_{k \leq t} \left[ 1 - \hat{P}_{\text{denom}}(C_{i,k} = 1)\right]
    = \prod_{k \leq t} \left[ 1 - \hat{P}(C_{i,k} = 1 \mid \bar{X}_{i,k}, A_i, \bar{C}_{i,k-1} = 0)\right]
$$
Since artificial censoring is monotonic (absorbing; once a patient is censored, they remain censored), this cumulative product is equivalent to the discrete-time analogue of Kaplan-Meier estimation of calculating the probability of remaining uncensored through month $t$. 

### 3. Stabilized IPCW
Next, we take the ratio of these products and obtain the stabilized IPCW weights:
$$
  \tilde{W}_{i,t} = \frac{ \widehat{\text{cum.num}}_{i,t} }{ \widehat{\text{cum.denom}}_{i,t} }
    = \prod_{k=0}^t \frac{ P(C_{i,k} = 0 \mid \bar{V}_{i}, A_i, \bar{C}_{i,k-1} = 0) }{ P(C_{i,k} = 0 \mid \bar{X}_{i,k}, A_i, \bar{C}_{i,k-1} = 0) }
$$
The intuition is: if an unhealthy patient is artificially censored, we find another patient who looks similar to them (same age, sex, ECOG, etc.) but was not censored and we up-weight that uncensored patient, it appears as though nobody was censored because the remaining patients are now representative of the original cohort.

This Robins/Hernan stabilized time-varying IPCW weight is: the denominator removes the confounding by adjusting for $\bar{X}_{i,k}$, while the numerator stabilizes the variance by keeping the marginal (strategy/time-only) censoring distribution in the numerator. Weights are truncated at 1st and 99th percentiles to bound the variance from near-positivity violations, at the cost of a very small truncation bias. The final stabilized IPCW weights are:
$$
  \tilde{W}_{i,t} = \prod_{k=0}^t \frac{ P(C_{i,k} = 0 \mid \bar{V}_{i}, A_i, \bar{C}_{i,k-1} = 0) }{ P(C_{i,k} = 0 \mid \bar{X}_{i,k}, A_i, \bar{C}_{i,k-1} = 0) } \qquad \text{stabilized IPCW for patient $i$ at time $t$}
$$

### 4. Weighted Penalized Outcome Model
Once every patient has their stabilized IPCW weights $\tilde{W}_{i,t}$ calculated, we can estimate the causal effect of discontinuing versus continuing treatment. TTE handles this via a discrete person-period data structure, where each patient contributes a distinct row for every month they remain at risk. Because target events (e.g., death or incidence of a specific irAE) are relatively rare within any single month, a pooled logistic regression model closely approximates a time-dependent Cox-PH model: rather than estimating a continuous hazard ratio, we predict the discrete-time hazard (the binary probability of event occurrence) in a given month, using the observed treatment history $A_{i,t}$ (the current entry of $\bar{A}_i$) rather than the assigned strategy alone. We model time flexibility using restricted cubic splines $f_{\text{ns}}(t)$ for the month variable, since the hazard is unlikely to be linear over time. The pooled logistic model is fit using the stabilized IPCW weights $\tilde{W}_{i,t}$ to account for the artificial censoring introduced by the TTE design, and is specified as:
$$
  \mathrm{logit} P(Y_{i,t} = 1 \mid A_i, A_{i,t}, \bar{V}_i, f_{\text{ns}}(t)) = \alpha_0 + \alpha_1 A_i + \alpha_2 A_{i,t} + f_{\text{ns}}(t) + \alpha^\top \bar{V}_i
$$ 
fit with weighted elastic net pooled logistic: 
$$  
  \min_{\alpha} \left\{ - \sum_{i=1}^n \tilde{W}_{i,t} \ell_{i}(\alpha) + \lambda \left( \frac{ 1-\alpha_{\text{mix}} }{2} \| \alpha \|_2^2 + \alpha_{\text{mix}} \| \alpha \|_1 \right) \right\}, \qquad \alpha_{\text{mix}} = 0.5
$$
and $\ell_{i}(\alpha) = Y_{i,t} \log P(Y_{i,t} = 1 \mid A_i, A_{i,t}, \bar{V}_i, f_{\text{ns}}(t)) + (1 - Y_{i,t}) \log P(Y_{i,t} = 0 \mid A_i, A_{i,t}, \bar{V}_i, f_{\text{ns}}(t))$, this $\ell_i(\alpha)$ represents the log-likelihood contribution of patient $i$ at month $t$ under the pooled logistic (dependent on the month-specific outcome $Y_{i,t}$). In oncology cohorts, the covariate history $\bar{X}_{i,k}$ is often likely extremely high-dimensional. To avoid the multicollinearity and overfitting from including too many covariates, the elastic net penalty stabilizes the estimation by shrinkage and variable selection. The penalty term is a convex combination of the LASSO and ridge penalties:
$$
  \lambda \left( \frac{ 1-\alpha_{\text{mix}} }{2} \underbrace{ \| \alpha \|_2^2 }_{\text{Ridge}} + \alpha_{\text{mix}} \underbrace{ \| \alpha \|_1 }_{\text{LASSO}} \right)
$$
The $L_1$ norm (LASSO) performs variable selection by shrinking some coefficients to zero, while the $L_2$ norm (ridge) handles multicollinearity by shrinking correlated coefficients together. 

The causal exposure terms (`assigned.strategy` $A_i$, `on.ici` $A_{i,t}$, and time spline month variables) are exempt from penalization (set `penalty.factor = 0` in glmnet) so shrinkage cannot bias the effect of primary interest or distort baseline time-trend — only nuisance confounders are regularized (see [Hahn et. al. (2018)](https://projecteuclid.org/journals/bayesian-analysis/volume-13/issue-2/Regularization-and-confounding-in-causal-inference/10.1214/17-BA1091.full), and next section for more). $\lambda$ is chosen using 1-SE rule from 10-fold cross-validation, favoring a more parsimonious and stable model. 

### 5. G-Computation (standardization) of marginal survival curves

For a fixed hypothetical treatment strategy $a$, predict each patient's monthly event probability under the counterfactual exposure trajectory (holding baseline covariates $\bar{V}_{i}$ at $t=21$ fixed, forcing the observed on-treatment history to follow the assigned regime $\bar{A}^\ast_i$ implied by strategy $a$, with time-$k$ entry $A^\ast_{i,k}$), then compute the cumulative survival curve for each patient and average across all patients to obtain the marginal survival curve under that treatment strategy:
$$
  \hat{S}_a(t) = \frac{1}{n} \sum_{i=1}^n \prod_{k = 21}^{t} 
    \left[ 
      1 - \hat{P}(Y_{i,k} = 1 \mid A_i = a, A^\ast_{i,k}, \bar{V}_{i}) 
    \right]
$$
This discrete-time g-formula estimator of marginal (population average) survival curves under each static treatment strategy — averaging over the empirical distribution of baseline covariates $\bar{V}_{i}$, rather than reporting conditional model coefficients. The difference in survival probability at 2-years is the primary estimand of interest, and the restricted mean survival time (RMST) over the 2-year follow-up period is also estimated.

### 6. Risk Difference and Non-collapsibility
$$ RD(48) = \hat{S}_{\text{cont}}(48) - \hat{S}_{\text{disc}}(48) $$
Recall that the odds ratio is a non-collapsible measure of effect: the marginal (population-averaged) odds ratio implied by the model does not equal the conditional (patient-specific) odds ratio $\exp(\alpha_1)$ for a given covariate pattern $\bar{V}_i$, even absent confounding. Reporting $\exp(\alpha_1)$ directly would therefore not be a valid population-level causal contrast. Standardizing via g-computation converts the non-collapsible conditional model into the collapsible, marginal risk difference $RD(48)$, which is the estimand we report as our primary causal contrast rather than the conditional model coefficients.


This modeling approach relies on: 

- Exchangeability: no unmeasured confounding of the observed treatment continuation decision (`on.ici`/discontinuation, $A_{i,t}$, the current entry of the observed history $\bar{A}_i$) and outcome, given $\bar{V}_{i}$ and $\bar{X}_{i,t}$. i.e., 
$$
  Y_{i,t}(a) \perp A_{i,t} \mid \bar{V}_{i}, \bar{X}_{i,t} \quad \text{for all } a, t
$$
- Positivity: adequately supported by weight truncation (1st and 99th percentiles) to avoid extreme weights. 
$$
  0 < P(A_{i,t} = 1 \mid \bar{V}_{i}, \bar{X}_{i,t}) < 1 \quad \text{for all } t
$$
- Consistency and correct model specification: relies on spline flexibility and elastic-net covariate selection to avoid model misspecification. 

--- 

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

Instead of a binary "continue" or "discontinue" at a fixed 2-year timepoint, we can now explore the more dynamic question of *when* the best time to discontinue ICI therapy is. This extends the two-arm target trial into a multi-arm design by expanding the treatment strategy to include multiple timepoints for discontinuation. I used every 6 months up to 48 months to mirror realistic clinical decision points and standard oncology follow-ups.

This fundamentally shifts the research question. We are no longer asking "*Given a patient has survived 2 years, should they stop?*" Instead, we are asking for a more powerful prospective question, "*From the moment of ICI initiation, what is the optimal planned duration of therapy to maximize overall survival while minimizing irAEs?* (irAEs not the main focus so far, but will be explored in future work).

## Data Structure & Cloning
The data structure must accommodate the fact that at baseline (month 1) we do not know which treatment strategy a patient will comply with. Since this is an observational study rather than a randomized trial, a patient's true trajectory will organically align with some strategies and violate others over time. To resolve this, the clone-censor approach: 

- Each patient's discrete-time longitudinal data is replicated 8 times, once for each treatment/target strategy. 
- Intent-to-Treat Emulation: At month 1, all 8 clones for a given patient are identical and are perfectly compliant with their respective arms. 
- Divergence: As follow-up progresses, clones are evaluated monthly against the specific rules for their assigned strategy. If a clone deviates from its assigned strategy, it is artificially censored at the time of deviation (as before). This allows us to estimate the per-protocol effect of each treatment strategy on overall survival and irAEs at 48 months after ICI initiation.


By restructuring the data in this way, we can shift from static patient-level baseline matrix to a massive person-period format allowing us to evaluate survival outcomes across multiple overlapping treatment regimens simultaneously. 

## Artificial Censoring & Trail Enforcement 

Becuase trials were not randomized to these 6-month intervals, we must artificially enforce the rules of our target trial. Clones are systematically removed (artificially censored) the moment their observed data deviates from their assigned strategy. 

In two mechanisms, we can enforce the trial rules:

1. **Stopped too early (pre-target deviation)**: if a clone is assigned to 24 month discontinuation arm, but the patient stops ICI therapy at month 9, the clone is compliant with up until month 9. At month 9 the clone violates the trial protocol and is censored.
2. **Failed to stop (post-target deviation)**: If a clone is assigned to 6-month discontinuation arm, but the patient continues therapy through month 8 (exceeding the 2-month grace period), the clone has failed to stop and is censored at month 8.

The rigid enforcement ensures that the remaining at-risk population in any given arm at any given month consists strictly with those who have followed their specific treatment strategy up to that point. 


## (Pre-model) Structural Diagnostics

Before fitting any outcome models, I wanted to visualize the mechanics of the cloned censoring cohort to understand the magnitude of the bias we're introducing through artificial censoring. 

### Mechanism and Timing of Artificial Censoring
![Mechanism and Timing of Artificial Censoring](figures/td-model/td-data/artificial_censoring_distribution.png)
The first visual maps out *why* and *when* artificial censoring occurs. The massive systematic spikes in "Fail to Stop" (blue) censoring exactly at target discontinuation months (with grace-period of 2 months). Conversely, "Stopped too early" (yellow) censoring is distributed more evenly across the earlier months. This confirms the trail rules are implemented correctly. 

### Structural Attrition of At-Risk Pool
![Structural Attrition of At-Risk Pool](figures/td-model/td-data/attrition_over_time.png)

The second visual shows the rapid depletion of our denominator. As the clones are censored due to non-compliance, the sample size in earlier discontinuation arms drop precipitously in months following their target stop date. This is expected, as the majority of patients will naturally continue therapy beyond 6 months, and thus the 6-month discontinuation arm is quickly depleted. This demonstrates why the downstream models require robust sample sizes and careful weight stabilization, since the data becomes sparse by design.


### Biased Baseline (Naive Event-Free Survival)
![Biased Baseline (Naive Event-Free Survival)](figures/td-model/td-data/naive_km_event_free_survival.png)
The third visual shows the "crude" (think g-computation without IPCW) discrete Kaplan-Meier survival curves for each arm. These curves are inherently biased due to our informed censoring. In the real world, patients who stop therapy early often do so because they're expecting toxicity or disease progression (and they are sicker). Because these sick patients naturally comply with early-stopping arms (and are artificially censored in continuous arms), the early-stopping arms will falsely appear to have worse survival than the continuous arms. 

This informative censoring is the exact reason we cannot naively compare the survival curves without IPCW.


##  Modeling (Treatment Discontinuation Model): 

In answering the question of when the optimal time to discontinue ICI therapy is, we must define what "optimal" means in this setting. Since there are two primary outcomes (overall survival and irAEs), we can define true optimality as being a composite of the two (but will discuss in terms of survival-optimality, and irae-optimality separately). But first the modeling approach will only be described in terms of overall survival (then extended later to irAEs).


## Pitfall of Composite Endpoint (`any.event` = death or irAE)
When defining "optimal" treatment duration, it is easiest to define a composite endpoint of "any event" = death (disease progression/mortality) or irAE. This composite outcome for patient $i$ at month $t$ is:
$$ Y^{\text{any}}_{i,t} = \max(Y^{\text{death}}_{i,t}, Y^{\text{irae}}_{i,t}) $$

However, this approach is problematic in the context of immunotherapy discontinuation. The clinical dynamics of ICI therapy create a zero-sum (meaning that the two outcomes are inversely related) tradeoff between the outcomes. From the DGP, we know that remaining on ICI ($A_{i,t} = 1$) provides some protective effects against mortality by keeping the tumor growth at bay, but simultaneously increases the risk of irAEs. Conversely, discontinuing ICI therapy ($A_{i,t} = 0$) eliminates the hazard of irAEs (toxicity), but increases the hazard of tumor progression and death. 

If we model $Y_{i,t}^{\text{any}}$ directly, the opposing directional effects of $A_{i,t}$ on the two outcomes will cancel each other out. The composite hazard remains functionally flat across all treatment strategies in this setting, yielding a risk difference of nearly 0% and a null conclusion that "treatment duration does not matter". This phenomenon is a form of "collider bias" in the causal inference literature, where the two outcomes are competing events and the composite outcome is a collider. Conditioning on the collider (the composite outcome) induces a spurious correlation between the two outcomes, which biases the estimated effect of treatment duration.


## Dual MSMs for Competing Risks
An MSM is a marginal structural model, which is a causal model that estimates the marginal (population-level) effect of a treatment on an outcome while accounting for time-varying confounding. In this case, we are interested in estimating the effect of ICI therapy continuation/discontinuation on two competing outcomes: mortality and irAEs. To isolate the effects of the clinical trade-off, we partition the outcome space into two seperate pooled logistic MSMs:
$$
\begin{align*}  
  \text{Mortality: }& \quad \mathrm{logit} [P(Y^{\text{death}}_{i,t} = 1 \mid \bar{X}_{i,t}, A_{i,t})] &= \delta_0 + \delta_1 A_{i,t} + f_{\text{ns}}(t) + \delta^\top \bar{V}_{i,t} \\
  \text{Toxicity: }& \quad \mathrm{logit} [P(Y^{\text{irae}}_{i,t} = 1 \mid \bar{X}_{i,t}, A_{i,t})] &= \theta_0 + \theta_1 A_{i,t} + f_{\text{ns}}(t) + \theta^\top \bar{V}_{i,t} \\
\end{align*}
$$
Both models utilize the same stabilized IPCW weight distribution $\tilde{W}_{i,t}$  derived from the clone-censoring design, ensuring that the artificial censoring and time-varying confounding is properly accounted for in both models. 

By running g-computation seperately on both models, we can generate two counterfactual survival curves for every treatment strategy. Standardizing these outcomes reveals a true tradeoff between the two outcomes:
1. irAE-free survival: strategies with early discontinuation (e.g., 6, 12 months) will exhibit significantly higher probabilities of remaining toxicity-free (since patients are no longer exposed to therapy). 
2. Overall survival: strategies with prolonged or continuous therapy will exhibit higher probabilities of overall survival (since patients are protected from tumor progression).

## Quasibinomial GLM over Penalized Elastic Net?
Earlier models used the elastic-net algorithm in `glmnet` with penalty factors to handle the high-dimensional covariate space while protecting structural variables from shrinkage. However, for the current simulation. However, we transitioned to a quasibinomial GLM:
$$ 
  \hat{\beta} = \arg\max_\beta \sum_{i,t} \tilde{W}_{i,t} \left[ Y_{i,t} \log \pi_{i,t} + (1 - Y_{i,t}) \log (1 - \pi_{i,t}) \right], \quad \pi_{i,t} = \mathrm{expit}(X_{i,t}^\top \beta) 
$$
because of the following reasons:
1. Standard `glmnet` routinely assumes iid data and standard binomial variance. When we introduce the stabilized IPCW weights, we inherenely inflate the variance of our sample (psueo-population). Using a quasibinomial family allows for overdispersion to be estimated from the data, providing more accurate standard errors that account for a weighted pseudo-population.
2. Parsimonious nuisance space: elastic net is far superior for high-dimensional covariate spaces, bu once the covariate space is reduced to a small set of clinically relevant confounders, (e.g., age, sex, ECOG, baseline labs), the regularization is no longer necessary. 
3. Compete elimination of shrinkage bias: while seetting `penalty.factor = 0` in `glmnet` protects the structural variables from shrinkage, Hahn et. al. demonstrate that residual regularization-induced confounding can still leak into structural estimates if unpenalized covariates share covariance with heavily penalized covariates. By using unpenalized quasibinomial GLM on a small pre-selected set of confounders, we entirely eliminate the risk of the prior bias away from out causal effect of interest. 


![Estimated Marginal Competing Risk Survival Curves](figures/td-model/td_marginal_competing_risks.png)


### Validation of Simulation

The surface level conclusion of this modeling approach concludes that continuing therapy (or discontinuing as late as 42 months) shows the highest overall survival (early discontinuation shows lowest survival), which confirms the simulation result. Continuing therapy for as long as possible yields the highest overall survival rates, with the shortest discontinuation arms (6, 12 months) showing the lowest survival. Staying on ICI thearpy provides a protective effect against mortality

Suprisingly, continuing therapy for toxicity also yielded higher irAE-free survival. Strategies that stop at 18 or 24 dip the lowest in irAE-free survival. This tells us to never stop treatment...

However, in the underlying DGP I explicitly programmed `on.ici == 1` to increase the log-odds of irAE by `+0.30`. But in the estimated coefficients for irAE, we see: `on.ici -0.0979118`, meaning the negative coefficient for active treatment, concluding that the drug protects patients from toxicity -— which is not correct—— and projects highest irAE-free survival for continuous therapy. This is called survivor bias. 

Patients who are suseptible to irAEs will develop them early in treatment and subsequently drop out or discontinue/censored. Patients who survive the continuous arm through 24, 36, 48 are healthier and posess unmeasured resiliency to toxicity. Because the IPCW weights do not perfectly adjust the hidden frailty, the outcome model looks at the long term survivors, sees that they're on treatment, and conclude falsely that the drug is protective against toxicity. This is a classic example of survivor bias, and is a limitation of the current modeling approach. **The survival model works well but the toxicity model needs to be revised.**



