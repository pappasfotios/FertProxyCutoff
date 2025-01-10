library(rstan)
library(dplyr)
library(data.table)
library(ggplot2)
library(foreach)
library(doParallel)

setwd("/FertProxyCutoff")

options(mc.cores = parallel::detectCores()-1)
rstan_options(auto_write = TRUE)
options(future.globals.maxSize = 50 * 1024^3)  # 50 GB

sperm1 <- readRDS("CASA2024.rds")
sperm1 <- sperm1[,c(8,4,14)]
sperm2 <- readRDS("spermPheno2023full.rds")
sperm2 <- sperm2[,c(7,2,9)]

names(sperm1) <- names(sperm2)

sperm <- rbind(sperm1, sperm2)

sperm <- sperm[complete.cases(sperm),]

sperm$Pop <- "2023"
sperm$Pop[grep("2021", sperm$sample)] <- "2024"

stan_model_dep <- stan_model("./BLCMparam_dep.stan")

init_function <- function() {
  list(
    Se1 = 0.8,
    Sp1 = 0.8,
    Se2 = 0.8,
    Sp2 = 0.8, 
    pi = rep(0.35, length(unique(stan_data$Pop)))
  )
}

# List to store results
results_list <- list()

# Iterate over quantiles from 0.1 to 0.5 by 0.1
quantiles <- seq(0.1, 0.5, by = 0.1)

cl <- makeCluster(floor(52/4))
registerDoParallel(cl)


results_list <- foreach(vcl_cut = quantiles, .combine = "c", .packages = c("rstan", "data.table", "dplyr")) %:% 
  foreach(conc_cut = quantiles, .combine = "c") %dopar% {
    
    print(paste(vcl_cut, conc_cut, sep = "/"))
    
    sperm$bin_tr <- ifelse(sperm$VCL < quantile(sperm$VCL, vcl_cut, na.rm = TRUE), "1", "0")
    sperm$bin_den <- ifelse(sperm$Sperm_den < quantile(sperm$Sperm_den, conc_cut, na.rm = TRUE), "1", "0")
    
    tests <- sperm[, 4:6]
    tests_prepared <- as.data.table(tests)[, `:=`(
      Pop = as.numeric(as.factor(Pop)),
      bin_tr = as.numeric(bin_tr),
      bin_den = as.numeric(bin_den)
    )]
    
    stan_data <- list(
      N = nrow(tests_prepared),
      P = length(unique(tests_prepared$Pop)),
      n = as.numeric(table(tests_prepared$Pop)),
      Pop = tests_prepared$Pop,
      Test1 = tests_prepared$bin_tr,
      Test2 = tests_prepared$bin_den
    )
    
    fit_dep <- sampling(
      stan_model_dep,
      init = init_function,
      data = stan_data,
      iter = 200000,
      warmup = 50000,
      thin = 40,
      chains = 4,
      cores = 4,
      seed = 123
    )
    
    summary_stats <- summary(fit_dep, pars = c("Se1", "Sp1", "Se2", "Sp2", "pi", "cdn", "cdp"))$summary
    result_dt <- as.data.table(summary_stats, keep.rownames = "parameter")
    result_dt[, `:=`(vcl_cut = vcl_cut, conc_cut = conc_cut)]
    result_dt <- result_dt %>%
      select(parameter, vcl_cut, conc_cut, mean, sd, n_eff, Rhat)
    
    setNames(list(result_dt), paste0("VCL_", vcl_cut, "_Conc_", conc_cut))
  }


stopCluster(cl)


# Combine all results into a single data frame
final_results_df <- bind_rows(results_list)

# Save results for further exploration
write.csv(final_results_df, "BLCMdep_results.csv", row.names = FALSE)
