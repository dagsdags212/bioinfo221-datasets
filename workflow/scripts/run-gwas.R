# install.packages("qtl")
library(qtl)

## DATA IMPORT
## Load dataset as a crossfile

filepath <- snakemake@params[["dataset"]]
# Load dataset as a crossfile
df <- read.cross(
  format = "csv", dir = "data", file = filepath, 
  genotypes = c("A", "H", "B"), 
  alleles = c("A", "B"))

# View cross file
summary(df)

## DATA PREPROCESSING
## CONVERT AND ADJUST GENETIC DISTANCES TO GENERATE MARKER MAP

# Convert physical distance into genetic distance
listeria.cross <- est.rf(df)

# Create a new marker map containing estimated recombination fractions
newmap <- est.map(listeria.cross, error.prob = 0.01)

# Visually compare adjusted distances
plotMap(listeria.cross, newmap)

listeria <- replace.map(listeria.cross, newmap)

## QTL MAPPING: Assumes only one QTL

# Provide positions of candidate QTLs and their genotypes
listeria <- calc.genoprob(listeria, step = 1, error.prob = 0.01)
# plotGeno(listeria)

# Simple QTL model: Interval Mapping
out.em <- scanone(listeria, method = "em")
em.targets <- summary(out.em, threshold = 3) # return only QTLs with LOD > 3
save(em.targets, file=snakemake@output[["em_model"]])
plot(out.em)

# Another QTL model: Haley-Knott
out.hk = scanone(listeria, method = "hk")
hk.targets <- summary(out.hk, threshold = 3)
save(hk.targets, file = snakemake@output[["hk_model"]])

# Fit conflicting QTLs to a linear model
## Create a qtl object
target.qtl <- makeqtl(listeria, chr = 3, pos = 19, what = "prob")

## Initialize linear model
linear.formula <- y ~ Q1
out.fq <- fitqtl(listeria, qtl = target.qtl, formula = linear.formula)

## Display and save regression summary
fq.model <- summary(out.fq)
save(fq.model, file = snakemake@output[["fq_model"]])

# Fit another linear model with a different method
out.fq.2 <- fitqtl(
  listeria, qtl = target.qtl, 
  formula = linear.formula, method = "imp")
summary(out.fq.2)

## QTL MAPPING: Assumes multiple QTLs
out.step <- stepwiseqtl(listeria, max.qtl = 5)
stepwise.model <- summary(out.step)
save(stepwise.model, file = snakemake@output[["stepwise_model"]])