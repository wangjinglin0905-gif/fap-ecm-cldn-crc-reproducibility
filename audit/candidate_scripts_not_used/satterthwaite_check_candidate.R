# Satterthwaite degrees-of-freedom sensitivity check for the mixed-effects
# models in Table 3 (manuscript Section 2.5 / 3.4)
# Run in R >= 4.4 with: install.packages(c("lme4", "lmerTest"))
#
# Input: a data.frame `cells` with one row per fibroblast-lineage cell and
# columns: patient (factor), FAP (numeric, z-scored), matrix4 (numeric,
# z-scored), receptor2 (numeric, z-scored) -- exactly as used for the
# lme4::lmer fits reported in the manuscript.

library(lmerTest)   # lmerTest::lmer = lme4::lmer + Satterthwaite df

m_matrix4   <- lmer(matrix4   ~ FAP + (1 | patient), data = cells)
m_receptor2 <- lmer(receptor2 ~ FAP + (1 | patient), data = cells)

s1 <- summary(m_matrix4)$coefficients
s2 <- summary(m_receptor2)$coefficients

cat("matrix4 ~ FAP   : slope =", round(s1["FAP","Estimate"], 3),
    " t =", round(s1["FAP","t value"], 2),
    " df(Satterthwaite) =", round(s1["FAP","df"], 1),
    " P =", format.pval(s1["FAP","Pr(>|t|)"], digits = 2), "\n")
cat("receptor2 ~ FAP : slope =", round(s2["FAP","Estimate"], 3),
    " t =", round(s2["FAP","t value"], 2),
    " df(Satterthwaite) =", round(s2["FAP","df"], 1),
    " P =", format.pval(s2["FAP","Pr(>|t|)"], digits = 2), "\n")

# Compare with the t approximation (df = n - 2) reported in the manuscript
# (slope = 0.112, t = 7.06, P < 0.001; slope = 0.051, t = 2.38, P = 0.017).
# If the Satterthwaite P values differ at the reported precision, update
# Table 3 and Section 3.4 accordingly.
