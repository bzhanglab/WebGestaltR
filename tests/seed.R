# The seed must reach the Rust permutations.
#
# This guards a failure that is silent by construction. The permutations are generated in
# Rust, and the compiled code reaching a user is not this repository's `webgestalt_rust`
# dependency but the frozen copy in `src/rust/vendor.tar.xz`. Re-point the dependency and
# forget to re-run `src/rust/vendor.sh` and everything still builds, every example still
# runs, and `seed` is silently ignored. Only comparing two seeded runs catches it.
#
# Deliberately no tolerance: same seed must give the same bits, not merely similar numbers.
# Kept small and offline so it costs little under R CMD check.

library(WebGestaltR)

set.seed(11) # shapes the synthetic input only; it cannot reach the Rust RNG
n <- 200
genes <- sprintf("G%03d", seq_len(n))
input_df <- data.frame(gene = genes, score = rnorm(n), stringsAsFactors = FALSE)
for (i in 1:6) {
    col <- rep(0, n)
    col[sample(n, 20)] <- 1
    input_df[[paste0("S", i)]] <- col
}

run <- function(seed) {
    r <- swGsea(input_df,
        thresh_type = "val", thresh = 0.5, p = 1, q = 1,
        perms = 100, min_set_size = 5, max_set_size = 500,
        nThreads = 1, rng_seed = seed
    )$Enrichment_Results
    r[order(rownames(r)), c("p_val", "fdr", "ES")]
}

a1 <- run(42)
a2 <- run(42)
b <- run(99)

# Same seed, identical to the bit.
stopifnot(isTRUE(all.equal(a1, a2, tolerance = 0)))

# Different seeds must actually differ, or the seed is being ignored in the other
# direction — a constant would pass the check above on its own.
stopifnot(!isTRUE(all.equal(a1, b, tolerance = 0)))

# The enrichment score is computed from the ranked list and the gene set with no
# randomness, so the seed must not move it. If this fails, the seed is reaching further
# than the permutations.
stopifnot(isTRUE(all.equal(a1$ES, b$ES, tolerance = 0)))
