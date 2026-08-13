# 03-td-sensitivity.r
# 
# Change cumulative incidence rate for each outcome at 36 months, to see what this does on the estimated survival curves. 
# per scenario in a grid of target 36 month cumulative incidence rates (one axis for death, one axis for irAE):
#   1. Calibrates the death/irae hazard intercept shifts (see death.intercept.shift / irae.intercept.shift in in td-model-fn.r) 
#      So the simulated 36 month 