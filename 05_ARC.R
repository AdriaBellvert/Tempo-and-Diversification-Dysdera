
library(raster)      
library(plyr)        
library(dplyr)       
library(dismo)       
library(phyloclim)   
library(jSDM)        

# Load binary data file
load("data_occurence_predictors_GM.bin")

# Run jSDM model
mod_frogs_jSDM_probit <- jSDM_binomial_probit(
  burnin = 5000,
  mcmc = 10000,
  thin = 20,
  presence_data = data_jSDM$occurence,
  site_formula = ~.,   
  site_data = data_jSDM$predictors,
  trait_data = data_jSDM$traits,
  n_latent = 2,
  site_effect = "random",
  alpha_start = 0, beta_start = 0,
  lambda_start = 0, W_start = 0,
  V_alpha = 1,
  shape = 0.1, rate = 0.1,
  mu_beta = 0, V_beta = 1,
  mu_lambda = 0, V_lambda = 1,
  seed = 1234, verbose = 1
)

# Get residual correlation
jSDM_res <- get_residual_cor(mod_frogs_jSDM_probit, prob = 0.95)
# Get environmental correlation
jSDM_env <- get_enviro_cor(mod_frogs_jSDM_probit, type = "mean", prob = 0.95)

# Extract species names
sps <- colnames(data_jSDM$occurence)

# Remove 'lancerotensis' from results
total.res <- jSDM_res$cor.mean
total.res <- total.res[!(row.names(total.res) %in% "lancerotensis"), !(colnames(total.res) %in% "lancerotensis")]

# Calculate pairwise residual correlations
res.d <- vector()
sps <- sps[!sps %in% "lancerotensis"]

for (i in 1:length(sps)) {
  for (j in (i + 1):length(sps)) {
    va <- total.res[i, j]
    res.d <- append(res.d, va)
  }
}

# Extract environmental correlations
total.env <- jSDM_env$cor
total.env <- total.env[!(row.names(total.env) %in% "lancerotensis"), !(colnames(total.env) %in% "lancerotensis")]
env.d <- vector()

for (i in 1:length(sps)) {
  for (j in (i + 1):length(sps)) {
    va <- total.env[i, j]
    env.d <- append(env.d, va)
  }
}

# Read species/islands data
isl <- read.table("islands.txt", header = FALSE, sep = "\t")
isl <- isl[!isl$V1 %in% "lancerotensis",]


# Load necessary libraries
library(phytools)   
library(geiger)    


# Read the tree file
DysderaTree <- read.tree("Dys_tr_names.tre")

# Change tip labels to complete names
name_mapping <- c(
  "ale_l" = "alegranzaensis",
  "amb_t" = "ambulotenta",
  "and_c" = "andamanae",
  "aneris" = "",
  "ara_c" = "arabisenen",
  "ban_c" = "bandamae",
  "bre_t" = "brevisetae",
  "brv_t" = "brevispina",
  "cal_g" = "calderensis_G",
  "cal_p" = "calderensis_P",
  "chi_t" = "chioensis",
  "cri_t" = "cribellata",
  "cur_t" = "curvisetae",
  "esq_t" = "esquiveli",
  "eng_g" = "enghoffi",
  "fla_t" = "banot",
  "gae_h" = "garoe",
  "gai_g" = "",
  "gib_t" = "gibbifera",
  "gol_t" = "",
  "gom_g" = "gomerensis_G",
  "gom_h" = "gomerensis_H",
  "gua_g" = "guayota",
  "guy_t" = "",
  "hei_c" = "herii",
  "her_t" = "",
  "hir_g" = "hirguan",
  "igu_t" = "iguanensis",
  "ins_t" = "insulana",
  "lab_t" = "labradaensis",
  "lancerotensis" = "lancerotensis",
  "lea_t" = "levipes",
  "lea_t" = "",
  "lev_g" = "levipes_G",
  "lio_c" = "liostethus",
  "lon_f" = "longa",
  "mac_t" = "macra",
  "mad_t" = "madai",
  "mah_f" = "mahan",
  "mah_l" = "mahan",
  "mar_t" = "marmoratta",
  "mat_t" = "aniepa",
  "min_t" = "minutissima",
  "mon_t" = "montanetensis",
  "nes_l" = "nesiotes",
  "nin_t" = "new_insulana",
  "ora_h" = "orahan",
  "pau_c" = "paucispinosa",
  "ram_g" = "ramblae",
  "rat_p" = "ratonensis",
  "rug_c" = "rugichelis",
  "san_f" = "sanborondon",
  "Dyssib_JL7" = "sibyllina",
  "sil_g" = "silvatica_G",
  "sil_h" = "silvatica_H",
  "sil_p" = "silvatica_P",
  "sim_l" = "simbeque",
  "spi_f" = "spinidorsum",
  "til_c" = "tilosensis",
  "und_g" = "undupe",
  "ung_t" = "unguimmanis",
  "ver_t" = "verneaui",
  "volcania" = "volcania",
  "ygu_c" = "yguanirae"
)


DysderaTree$tip.label <- revalue(DysderaTree$tip.label, name_mapping)

# Read islands data
isl <- read.table("islands.txt", header = FALSE, sep = "\t")

# remove D. lancerotensis
isl2 <- isl[!isl$V1 %in% "lancerotensis", ]
rownames(isl2) <- isl2$V1

# Drop tips not found in the data
obj <- name.check(DysderaTree, isl2)
DysderaTreeTip <- drop.tip(DysderaTree, obj$tree_not_data)
obj <- name.check(DysderaTreeTip, isl2)


###ARC
# Calculate age-range correlation
x <- age.range.correlation(phy = DysderaTreeTip, overlap = jSDM_res$cor.mean, n = 999)

# Plot age-range correlation
plot(x$age.range.correlation)
