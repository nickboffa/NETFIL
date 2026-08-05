library(jsonlite)

# Equation (1) from Additional File 1 (May 1977):
# p_mf = 2 - p_A - 2^(1+k0) * (1-p_A) * (1 + (1-p_A)^(1/k0))^(-k0)
# Solves for k0 given observed p_A and p_mf.
#
# For ASM: p_mf was not directly measured, so k0 was solved in Mathematica
# using p_mf inferred from the Xu et al. p_A/p_mf relationship (equation 2).
# The resulting table lives in data/ASM/initaggs.csv and is kept as-is.
#
# For WSM: p_mf IS directly observed, so equation (2) is skipped entirely.
# This script generates data/WSM/initaggs.csv from the survey values in
# country.json. Only one row is needed because initaggs is only ever looked
# up at init_ant_prev, which is fixed.

mf_from_k0 <- function(k0, p_A) {
  2 - p_A - 2^(1 + k0) * (1 - p_A) * (1 + (1 - p_A)^(1/k0))^(-k0)
}

k0_from_pA_pmf <- function(p_A, p_mf) {
  uniroot(
    f        = function(k0) mf_from_k0(k0, p_A) - p_mf,
    interval = c(1e-8, 1e3)
  )$root
}

cfg     <- fromJSON("data/WSM/country.json")
p_A     <- cfg$seeding$init_ant_prev
p_mf    <- (cfg$seeding$init_mf_prev_min + cfg$seeding$init_mf_prev_max) / 2

k0 <- k0_from_pA_pmf(p_A, p_mf)

cat("p_A  =", p_A,  "\n")
cat("p_mf =", p_mf, "(midpoint of acceptance bounds)\n")
cat("k0   =", k0,   "\n")
cat("Check: mf_from_k0 =", mf_from_k0(k0, p_A), "\n")

write.table(
  data.frame(p_A = p_A, k0 = k0),
  file      = "data/WSM/initaggs.csv",
  sep       = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)
cat("Written: data/WSM/initaggs.csv\n")
