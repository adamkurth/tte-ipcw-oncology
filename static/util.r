# ---------- HELPERS ----------
# helper: sample factor with probabilities
sample_factor <- function(levels, probs, n, ordered = FALSE) {
    stopifnot(length(levels) == length(probs), all(probs >= 0), n >= 0)
    probs <- probs / sum(probs)
    if (n == 0) {
        return(factor(character(0), levels = levels, ordered = ordered))
    }
    x <- sample(levels, n, replace = TRUE, prob = probs)
    factor(x, levels = levels, ordered = ordered)
}

# helper: draw stratified types given exact counts (e.g., cancer types)
draw_stratified_types <- function(n, labels, probs) {
    stopifnot(length(labels) == length(probs), n > 0)
    probs <- probs / sum(probs)
    base_counts <- floor(n * probs)
    remainder <- n - sum(base_counts)
    if (remainder > 0) {
        extra <- sample.int(length(labels), size = remainder, replace = TRUE, prob = probs)
        base_counts <- base_counts + tabulate(extra, nbins = length(labels))
    }
    out <- rep(labels, times = base_counts)
    sample(out, size = n, replace = FALSE)
}

# helper: expit
expit <- function(x) 1 / (1 + exp(-x))

# helper: approximate month addition (30.4375 days per month)
add_months_approx <- function(d, m) d + as.integer(round(30.4375 * m))

# helper: multinomial sampling function for ordered factors
sample_multicat <- function(levels, logits) { 
    logits <- as.numeric(logits)
    p <- exp(logits - max(logits))
    p <- p / sum(p)
    sample(levels, size = 1, prob = p)
}
