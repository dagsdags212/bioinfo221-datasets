# install.packages("qtl")
library(qtl)

## DATA IMPORT
# Load dataset as a crossfile
filepath <- "listeria.csv"
df <- read.cross(
  format = "csv", dir = "data", file = filepath,
  genotypes = c("A", "H", "B"),
  alleles = c("A", "B")
)

# View cross file
summary(df)

## DATA PREPROCESSING
## CONVERT AND ADJUST GENETIC DISTANCES TO GENERATE MARKER MAP

# Convert physical distance into genetic distance
listeria_cross <- est.rf(df)

# Create a new marker map containing estimated recombination fractions
newmap <- est.map(listeria_cross, error.prob = 0.01)

# Visually compare adjusted distances
plotMap(listeria_cross, newmap)

listeria <- replace.map(listeria_cross, newmap)

## QTL MAPPING: Assumes only one QTL

# Provide positions of candidate QTLs and their genotypes
listeria <- calc.genoprob(listeria, step = 1, error.prob = 0.01)
plotGeno(listeria)

# Simple QTL model: Interval Mapping
out_em <- scanone(listeria, method = "em")
em_targets <- summary(out_em, threshold = 3) # return only QTLs with LOD > 3
save(em_targets, file = "output/em_out.RData")
plot(out_em)

# Another QTL model: Haley-Knott
out_hk <- scanone(listeria, method = "hk")
hk_targets <- summary(out_hk, threshold = 3)
save(hk_targets, file = "output/hk_out.RData")

# Fit conflicting QTLs to a linear model
## Create a qtl object
target_qtl <- makeqtl(listeria, chr = 3, pos = 19, what = "prob")

## Initialize linear model
linear_formula <- y ~ Q1
out_fq <- fitqtl(listeria, qtl = target_qtl, formula = linear_formula)

## Display and save regression summary
fq_model <- summary(out_fq)
save(fq_model, file = "output/fq_model.RData")

# Fit another linear model with a different method
out_fq_2 <- fitqtl(
  listeria,
  qtl = target_qtl,
  formula = linear_formula, method = "imp"
)
summary(out_fq_2)

## QTL MAPPING: Assumes multiple QTLs
out_step <- stepwiseqtl(listeria, max.qtl = 5)
stepwise_model <- summary(out_step)
save(stepwise_model, file = "output/stepwise_model.RData")
