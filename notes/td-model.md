
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





<!-- 
The surface level conclusion of this modeling approach concludes that continuing therapy (or discontinuing as late as 42 months) shows the highest overall survival (early discontinuation shows lowest survival), which confirms the simulation result. Continuing therapy for as long as possible yields the highest overall survival rates, with the shortest discontinuation arms (6, 12 months) showing the lowest survival. Staying on ICI thearpy provides a protective effect against mortality

Suprisingly, continuing therapy for toxicity also yielded higher irAE-free survival. Strategies that stop at 18 or 24 dip the lowest in irAE-free survival. This tells us to never stop treatment...

However, in the underlying DGP I explicitly programmed `on.ici == 1` to increase the log-odds of irAE by `+0.30`. But in the estimated coefficients for irAE, we see: `on.ici -0.0979118`, meaning the negative coefficient for active treatment, concluding that the drug protects patients from toxicity -— which is not correct—— and projects highest irAE-free survival for continuous therapy. This is called survivor bias. 

Patients who are suseptible to irAEs will develop them early in treatment and subsequently drop out or discontinue/censored. Patients who survive the continuous arm through 24, 36, 48 are healthier and posess unmeasured resiliency to toxicity. Because the IPCW weights do not perfectly adjust the hidden frailty, the outcome model looks at the long term survivors, sees that they're on treatment, and conclude falsely that the drug is protective against toxicity. This is a classic example of survivor bias, and is a limitation of the current modeling approach. **The survival model works well but the toxicity model needs to be revised.** -->



NONMONOTONICITY IS WEIRD
RUN WITH DIFFERENT SEED
CHANGE CUMULATIVE INCIDENCE EACH OUTCOME AT 36 MONTHS 
 -> SUMMARIZE EACH RESULT IN TERMS OF STD ERROR 

MIMIC ROUGH SIMULATION ATTRITION RATES
READ LUNG CANCER PAPER ON
READ NUMBER OF PEOPLE WHO ARE TREATED (REBECCA TABLES / PRELIMINARY DATA)


sanity check: 
calculate true survival probabilites over all conditional probabilities, and calculate (not estimate) true counterfactual survival probabilities
validate with estimated survival probabilities from model, and compare to true survival probabilities.




(later) erroneous censoring
