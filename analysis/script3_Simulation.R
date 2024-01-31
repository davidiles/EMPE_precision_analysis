# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2',
              'scales','tidyverse',
              'ggrepel','ggthemes','ggpubr')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/EMPE_Global/analysis")

rm(list=ls())

# --------------------------------------
# Load data package and fitted model
# --------------------------------------
load("output/EMPE_data_formatted.RData") # Data
load(file = "output/fitted_model.RData") # Fitted model

# --------------------------------------
# # JAGS script used to simulate data, based on jags.data.sim
# --------------------------------------

sink("EMPE_model_simulate_data.jags")
cat("
    model {
  
  # ------------------------------------
  # Population process model
  # ------------------------------------
  
  # Simulate initial abundances for each colony (X[s,1])
  logX1_tau <- pow(logX1_sigma,-2) 
  for (s in 1:n_sites){
    log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
    N[s,1] <- exp(log_X[s,1]) * z_occ[s,1]
  }
  
  # Simulate occupancy status in each year, at each colony
  for (s in 1:n_sites){
    for (t in 1:n_years){
    
      # True occupancy state at colony in each year
      z_occ[s,t] ~ dbern(prob_occ)
      
    }
  }
  
  # Simulate temporal dynamics at each colony
  r_mean_grandmean_tau <- pow(r_mean_grandmean_sigma,-2)
  r_tau <- pow(r_sigma,-2)

  for (s in 1:n_sites){
    
    r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau) 

    for (t in 1:(n_years-1)){
      log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s] - 1/(2*r_tau), r_tau)
      N[s,t+1] <- exp(log_X[s,t+1]) * z_occ[s,t+1]
    }
  } 
  
  # ------------------------------------
  # Aerial observation model
  # ------------------------------------
  
  aerial_tau <- pow(aerial_sigma,-2)
  
  for (i in 1:n_obs_aerial){
    
    # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
    lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] + aerial_DoS[i]*DoS_slope - 0.5*pow(aerial_sigma,2), aerial_tau)
    adult_expected[i] <- z_occ[aerial_site[i],aerial_year[i]] * lambda[i]
    adult_count[i] ~ dpois(adult_expected[i])
    
  }
  
  # ------------------------------------
  # Satellite observation model
  # ------------------------------------
  
  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_expected[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]] * exp(DoS_slope * satellite_DoS[i])
    sat_sigma[i] <- sat_expected[i] * sat_CV[img_qual[i]] + 0.00001 # Tiny constant ensures tau is define when var is zero
    sat_tau[i] <- pow(sat_sigma[i],-2)
    satellite[i] ~ dnorm(sat_expected[i], sat_tau[i])

  }
  
  # ------------------------------------
  # Derived quantities
  # ------------------------------------
  
  # Global abundance each year
  for (t in 1:n_years){
    N_global[t] <- sum(N[1:n_sites,t])
  }
  
  # Global trend (measured as least squares regression slope on log scale)
  global_trend <- inprod(log(N_global[1:n_years]),regression_weights[1,1:n_years])
  
}
    ",fill = TRUE)
sink()

# --------------------------------------
# Prepare an empty file to store results of simulations
# --------------------------------------

# Create empty output file to store simulation results
if (!file.exists("./output/simulation/sim_results.RData")){
  simulation_results <- data.frame()
  save(simulation_results, file = "./output/simulation/sim_results.RData") 
} else{
  load(file = "./output/simulation/sim_results.RData")
}

# --------------------------------------
# Conduct repeated simulations in which new data is simulated,
# and statistical model is applied to estimate global abundance and trend
# --------------------------------------

for (sim_run in seq(1,1000,1)){
  
  set.seed(sim_run)
  
  print(sim_run)
  if (file.exists("./output/simulation/sim_results.RData")) load(file = "./output/simulation/sim_results.RData")
  if (nrow(simulation_results) > 0 & sim_run %in% simulation_results$seed) next
  
  # --------------------------------------
  # - Prepare a data package to use for simulating new data
  # - Use sample sizes (e.g., n_sites, n_years, n_obs, etc) from observed data to
  #   ensure simulations have "realistic" data quantity and balance
  # --------------------------------------
  
  jags.data.sim <- jags.data[c("n_years","n_sites",
                               "n_obs_aerial","aerial_site","aerial_year","aerial_DoS",
                               "n_obs_satellite","satellite_site","satellite_year","satellite_DoS","img_qual",
                               "regression_weights")]
  
  # --------------------------------------
  # Select a sample from the (non-joint) posterior of the fitted model to specify
  # 'true' parameter values for this simulation
  # --------------------------------------
  
  # Non-joint samples, to represent a wider range of population dynamic scenarios
  jags.data.sim$prob_occ <- out$sims.list$prob_occ[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$r_mean_grandmean_mu <- out$sims.list$r_mean_grandmean_mu[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$r_mean_grandmean_sigma <- out$sims.list$r_mean_grandmean_sigma[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$logX1_mean <- out$sims.list$logX1_mean[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$logX1_sigma <- out$sims.list$logX1_sigma[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$r_sigma <- out$sims.list$r_sigma[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$aerial_sigma <- out$sims.list$aerial_sigma[sample(1:out$mcmc.info$n.samples,1)]
  jags.data.sim$sat_slope <- out$sims.list$sat_slope[sample(1:out$mcmc.info$n.samples,1),1:3]
  jags.data.sim$sat_CV <- out$sims.list$sat_CV[sample(1:out$mcmc.info$n.samples,1),1:3]
  jags.data.sim$DoS_slope <- out$sims.list$DoS_slope[sample(1:out$mcmc.info$n.samples,1)]
  
  out_sim <- jags(data=jags.data.sim,
                  model.file="EMPE_model_simulate_data.jags",
                  parameters.to.save=c(
                    
                    "N",
                    "adult_count",
                    "satellite",
                    "N_global",
                    "global_trend"
                    
                  ),
                  inits = NULL,
                  n.chains=1,
                  n.thin = 1,
                  n.iter= 2,
                  n.burnin= 1,
                  parallel = TRUE)
  
  # ------------------------------------
  # Extract simulated data and plot
  # ------------------------------------
  sim_data <- out_sim$sims.list
  
  # True (simulated) values of annual N at each colony, in this simulation run
  N_df <- sim_data$N[1,,] %>% melt() %>% rename(site = Var1, year = Var2, N = value)
  
  # Aerial observations
  aerial_df_sim <- data.frame(y = sim_data$adult_count[1,],
                              year = jags.data.sim$aerial_year,
                              site = jags.data.sim$aerial_site,
                              type = "Aerial count (adult)")
  
  # Satellite observations
  satellite_df_sim <- data.frame(y = sim_data$satellite[1,],
                                 year = jags.data.sim$satellite_year,
                                 site = jags.data.sim$satellite_site,
                                 type = "Satellite")
  # Combine into a dataframe
  obs_sim <- rbind(aerial_df_sim,satellite_df_sim)
  
  # Plot observed data at each colony (for visual confirmation that the simulation is reasonable)
  sim_data_plot <- ggplot(obs_sim) +
    geom_line(data = N_df, aes(x = year, y = N))+
    geom_point(aes(x = year, y = y, shape = type))+
    scale_shape_manual(name = 'Obs type',
                       values =c(19,4))+
    
    facet_wrap(site~., scales = "free")+
    theme_bw()
  
  # Global abundance each year
  N_global_df <- data.frame(year = 1:jags.data.sim$n_years,N = out_sim$sims.list$N_global[1,] )
  
  # Plot global abundance across the simulation
  global_N_plot <- ggplot(N_global_df) +
    geom_line(data = N_global_df, aes(x = year, y = N))+
    geom_hline(yintercept = 0, col = "transparent")+
    theme_bw()
  
  global_N_plot 
  
  # --------------------------------------
  # Fit model to simulated data
  # --------------------------------------
  
  # Data to fit
  jags.data.refit = list( n_years = jags.data.sim$n_years,
                          n_sites = jags.data.sim$n_sites,
                          
                          # aerial counts of adults
                          n_obs_aerial = length(sim_data$adult_count[1,]),
                          adult_count = sim_data$adult_count[1,],
                          aerial_site = jags.data.sim$aerial_site,
                          aerial_year = jags.data.sim$aerial_year,
                          aerial_DoS = jags.data.sim$aerial_DoS,
                          
                          # satellite counts
                          n_obs_satellite = length(sim_data$satellite[1,]),
                          satellite = sim_data$satellite[1,],
                          satellite_site = jags.data.sim$satellite_site,
                          satellite_year = jags.data.sim$satellite_year,
                          satellite_DoS = jags.data.sim$satellite_DoS,
                          img_qual = jags.data.sim$img_qual,
                          
                          regression_weights = jags.data.sim$regression_weights
  )
  
  z_init <- matrix(1,ncol = jags.data.refit$n_years, nrow = jags.data.refit$n_sites)
  inits <- function()list(z_occ = z_init,
                          prob_occ = jags.data.sim$prob_occ,
                          sat_slope = jags.data.sim$sat_slope,
                          sat_CV = jags.data.sim$sat_CV,
                          DoS_slope = jags.data.sim$DoS_slope
  )
  
  out_refit <- jags(data=jags.data.refit,
                    model.file="EMPE_model_empirical.jags",
                    parameters.to.save=c(
                      
                      # ------------------------
                      # Hyper-parameters
                      # ------------------------
                      "prob_occ",                # Probability colonies are "occupied"
                      "r_mean_grandmean_mu",     # Hierarchical grand mean of colony-level annual growth rates
                      "r_mean_grandmean_sigma",  # Hierarchical sd of colony-level growth rates
                      "logX1_mean",              # Hierarchical grand mean of colony-level initial abundances
                      "logX1_sigma",             # Hierarchical sd of colony-level initial abundances
                      "r_sigma",                 # Temporal variance of colony-level annual growth rates
                      
                      "DoS_slope",
                      "aerial_sigma",            # SD of aerial observations (on log scale)
                      "sat_slope",               # Bias in satellite observations
                      "sat_CV",                  # Coefficient of variation in satellite observations
                      
                      # ------------------------
                      # Colony-specific quantities
                      # ------------------------
                      
                      # Colony-specific mean annual differences (could be used as a measure of colony-specific "trend")
                      "r_mean",
                      
                      # Colony-specific index of abundance each year
                      "N",
                      
                      # ------------------------
                      # Global index and trend estimates
                      # ------------------------
                      
                      # Global index of abundance each year
                      "N_global",
                      
                      # Log-linear OLS slope across the study
                      "global_trend",
                      
                      # ------------------------
                      # Observation-specific quantities
                      # ------------------------
                      
                      # Fitted values
                      "adult_expected",
                      "sat_expected",
                      "sat_sigma",
                      
                      # Discrepancy measures for posterior predictive checks
                      "RMSE_adult_count_actual",
                      "RMSE_satellite_actual",
                      "RMSE_adult_count_sim",
                      "RMSE_satellite_sim",
                      
                      # Simulated datasets from fitted model (also for goodness-of-fit testing)
                      "sim_adult_count",
                      "sim_satellite"
                      
                    ),
                    
                    inits = inits,
                    n.chains=3,
                    n.thin = 10,
                    n.iter= 30000,
                    n.burnin= 10000,
                    parallel = TRUE)
  
  #----------------------------------------------------------
  # Examine model convergence
  #----------------------------------------------------------
  
  hyper_parameters <- c(
    
    # ------------------------
    # Hyper-parameters
    # ------------------------
    "prob_occ",               # Probability colonies are "occupied"
    "r_mean_grandmean_mu",    # Hierarchical grand mean of colony-level annual growth rates
    "r_mean_grandmean_sigma", # Hierarchical sd of colony-level growth rates
    "logX1_mean",             # Hierarchical grand mean of colony-level initial abundances
    "logX1_sigma",            # Hierarchical sd of colony-level initial abundances
    "r_sigma",                # Temporal variance of colony-level annual growth rates
    "DoS_slope",              # Effect of "day of season"
    "aerial_sigma",           # SD of aerial observations (on log scale)
    "sat_slope",              # Bias in satellite observations
    "sat_CV"                 # Coefficient of variation in satellite observations
    
  )
  
  Rhat_hyper <- unlist(out_refit$Rhat[which(names(out_refit$Rhat) %in% hyper_parameters)])
  mean(Rhat_hyper > 1.1, na.rm = TRUE) # Proportion of parameters with Rhat > 1.1

  indices_and_trends <- c(
    
    # ------------------------
    # Global estimates
    # ------------------------
    
    # Global abundance each year
    "N_global",
    
    # Log-linear OLS slope across the study
    "global_trend",
    
    # ------------------------
    # Colony-level estimates
    # ------------------------
    
    # Colony-specific mean annual differences
    "r_mean",
    
    # Colony-specific abundance each year
    "N"
    
  )
  
  Rhat_indices <- unlist(out_refit$Rhat[which(names(out_refit$Rhat) %in% indices_and_trends)])
  mean(Rhat_indices > 1.1, na.rm = TRUE) # Proportion of parameters with Rhat > 1.1
  
  # -----------------------------------------
  # Extract estimates of global trend and change estimates
  # -----------------------------------------
  
  # Percent change estimates
  percent_change_est = 100*(out_refit$sims.list$N_global[,jags.data$n_years] - out_refit$sims.list$N_global[,1])/out_refit$sims.list$N_global[,1]
  percent_change_true = 100*(out_sim$sims.list$N_global[,jags.data$n_years] - out_sim$sims.list$N_global[,1])/out_sim$sims.list$N_global[,1]
  
  # -----------------------------------------
  # Conduct goodness-of-fit evaluation (Bayesian p-values)
  # -----------------------------------------
  Bayesian_pval_adult <- mean(out_refit$sims.list$RMSE_adult_count_actual > out_refit$sims.list$RMSE_adult_count_sim)
  Bayesian_pval_satellite <- mean(out_refit$sims.list$RMSE_satellite_actual > out_refit$sims.list$RMSE_satellite_sim)
  
  # -----------------------------------------
  # Save results for this simulation
  # -----------------------------------------
  
  simulation_summary <- data.frame(seed = sim_run,
                                   Rhat_hyper = mean(Rhat_hyper > 1.1, na.rm = TRUE),      # Proportion of non-converged hyper parameters
                                   Rhat_indices = mean(Rhat_indices > 1.1, na.rm = TRUE),  # Proportion of non-converged hyper parameters
                                   Bayesian_pval_adult = Bayesian_pval_adult,
                                   Bayesian_pval_satellite = Bayesian_pval_satellite,
                                   
                                   global_trend_est_mean = mean(out_refit$sims.list$global_trend),
                                   global_trend_est_q025 = quantile(out_refit$sims.list$global_trend,0.025),
                                   global_trend_est_q500 = quantile(out_refit$sims.list$global_trend,0.500),
                                   global_trend_est_q975 = quantile(out_refit$sims.list$global_trend,0.975),
                                   global_trend_true = out_sim$sims.list$global_trend,
                                   
                                   percent_change_est_mean = mean(percent_change_est),
                                   percent_change_est_q025 = quantile(percent_change_est,0.025),
                                   percent_change_est_q500 = quantile(percent_change_est,0.500),
                                   percent_change_est_q975 = quantile(percent_change_est,0.975),
                                   percent_change_true = percent_change_true)
  
  
  # -----------------------------------------
  # Save output
  # -----------------------------------------
  
  load(file = "./output/simulation/sim_results.RData")
  simulation_results = rbind(simulation_results, simulation_summary)
  save(simulation_results,file = "./output/simulation/sim_results.RData")
  
}

# ***********************************************************************
# Summarize simulation results 
#     - (bias / coverage of trend estimates)
#     - distribution of Bayesian p-values under a correctly specified model
# ***********************************************************************

load(file = "./output/simulation/sim_results.RData")

simulation_results <- subset(simulation_results, Rhat_hyper == 0 & Rhat_indices == 0)

# ---------------------------------------
# Trend (average annual percent change from 2009 to 2018)
# ---------------------------------------

# Convert to percent change per year using 100*(exp(OLS_regression_slope)-1)
simulation_results$global_trend_est_mean = 100*(exp(simulation_results$global_trend_est_mean)-1)
simulation_results$global_trend_est_q025 = 100*(exp(simulation_results$global_trend_est_q025)-1)
simulation_results$global_trend_est_q500 = 100*(exp(simulation_results$global_trend_est_q500)-1)
simulation_results$global_trend_est_q975 = 100*(exp(simulation_results$global_trend_est_q975)-1)
simulation_results$global_trend_true = 100*(exp(simulation_results$global_trend_true)-1)

# Credible interval coverage
simulation_results$trend_cov <- simulation_results$global_trend_est_q025 <= simulation_results$global_trend_true & simulation_results$global_trend_est_q975 >= simulation_results$global_trend_true
trend_coverage <- mean(simulation_results$trend_cov) %>% round(2)
med_trend_bias <- median(simulation_results$global_trend_est_mean - simulation_results$global_trend_true) %>% round(1)

ylim = range(simulation_results[,c("global_trend_true","global_trend_est_q025","global_trend_est_q975")])
trend_plot <- ggplot(simulation_results,aes(x = global_trend_true, 
                                                   y = global_trend_est_mean, 
                                                   ymin = global_trend_est_q025,
                                                   ymax = global_trend_est_q975, 
                                                   col = trend_cov))+
  geom_abline(slope = 1, col = "gray90")+
  geom_errorbar(width = 0, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("orangered","dodgerblue"), name = "95% CRI overlaps\ntrue trend")+
  coord_cartesian(ylim = c(ylim[1],ylim[2]), xlim = c(ylim[1],ylim[2]))+
  xlab("True (simulated) global trend")+
  ylab("Estimated global trend")+
  labs(title = paste0("Simulation results"),
       subtitle = paste0("Median bias of trend estimate = ",med_trend_bias,"%\n95% credible interval coverage = ",trend_coverage*100,"%"))+
  theme_bw()
#print(trend_plot)

png("./output/simulation/trend_sim_results.png", width = 5, height = 4, units = "in",res=500)
print(trend_plot)
dev.off()

# ---------------------------------------
# Total percent change (2018 vs 2009)
# ---------------------------------------

# Credible interval coverage
simulation_results$percent_change_cov <- simulation_results$percent_change_est_q025 <= simulation_results$percent_change_true & simulation_results$percent_change_est_q975 >= simulation_results$percent_change_true
percent_change_coverage <- mean(simulation_results$percent_change_cov) %>% round(2)
med_percent_change_bias <- median(simulation_results$percent_change_est_mean - simulation_results$percent_change_true) %>% round(1)

# Just to help identify a useful scale
ymax = 400
ymax_scale = log(ymax/100+1)
ymin_scale = -ymax_scale
ymin = 100*(exp(ymin_scale)-1)

ylim = 400
ylim = log(ylim/100 + 1)
y_labels = data.frame(pchange = c(-80,-50,0,100,400),labels = c("-80%","-50%","0","+100%","+400%"))
y_labels$r = log(y_labels$pchange/100 + 1)

percent_change_plot <- ggplot(simulation_results,
                              aes(x = log(percent_change_true/100 + 1), 
                                  y = log(percent_change_est_mean/100 + 1), 
                                  ymin = log(percent_change_est_q025/100 + 1),
                                  ymax = log(percent_change_est_q975/100 + 1), 
                                  col = percent_change_cov))+
  geom_abline(slope = 1, col = "gray90")+
  geom_errorbar(width = 0, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("orangered","dodgerblue"), name = "95% CRI overlaps\ntrue trend")+
  scale_y_continuous(breaks = y_labels$r, labels = y_labels$labels)+
  scale_x_continuous(breaks = y_labels$r, labels = y_labels$labels)+
  
  coord_cartesian(ylim=c(-ylim,ylim),xlim = c(-ylim,ylim))+
  xlab("True (simulated) percent\nchange from 2009 to 2018")+
  ylab("Estimated percent\nchange from 2009 to 2018")+
  labs(title = paste0("Simulation results"),
       subtitle = paste0("Median bias of change estimate = ",med_percent_change_bias,"%\n95% credible interval coverage = ",percent_change_coverage*100,"%"))+
  theme_bw()
#print(percent_change_plot)

png("./output/simulation/change_sim_results.png", width = 5, height = 4, units = "in",res=500)
print(percent_change_plot)
dev.off()

# ---------------------------------------
# Distribution of Bayesian p-values under correctly specified model
# ---------------------------------------

pvals_adult <- ggplot(simulation_results, aes(x = Bayesian_pval_adult))+
  geom_histogram(fill = "dodgerblue",binwidth = 0.05)+
  theme_few()+
  xlab("Bayesian p-value")+
  ylab("Number of simulations")+
  ggtitle("Bayesian p-value (adult counts)")+
  scale_x_continuous(breaks = seq(-0.2,1.2,0.2), limits = c(0,1))

pvals_satellite <- ggplot(simulation_results, aes(x = Bayesian_pval_satellite))+
  geom_histogram(fill = "dodgerblue",binwidth = 0.05)+
  theme_few()+
  xlab("Bayesian p-value")+
  ylab("Number of simulations")+
  ggtitle("Bayesian p-value (satellite surveys)")+
  scale_x_continuous(breaks = seq(-0.2,1.2,0.2), limits = c(0,1))

pval_plot <- ggarrange(pvals_adult,pvals_satellite,nrow=2)

png("./output/simulation/Bayesian_pval_plot.png", width = 5, height = 7, units = "in",res=500)
print(pval_plot)
dev.off()