# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2',
              'scales','tidyverse',
              'ggrepel','ggthemes','ggpubr')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/EMPE_precision_analysis/analysis")

rm(list=ls())

# --------------------------------------
# Load data package and fitted model
# --------------------------------------
load("output/EMPE_data_formatted.RData") # Data
load(file = "output/fitted_model.RData") # Fitted model

# --------------------------------------
# JAGS script used to simulate data
# --------------------------------------

sink("EMPE_model_simulate_observations.jags")
cat("
    model {
  
  # Simulate occupancy status in each year, at each colony
  for (s in 1:n_sites){
    for (t in 1:n_years){
    
      # True occupancy state at colony in each year
      z_occ[s,t] ~ dbern(prob_occ)
      N[s,t] <- exp(log_X[s,t]) * z_occ[s,t]
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
if (!file.exists("./output/simulation/precision_results.RData")){
  precision_results <- data.frame()
  save(precision_results, file = "./output/simulation/precision_results.RData") 
} else{
  load(file = "./output/simulation/precision_results.RData")
}

# --------------------------------------
# Conduct repeated simulations in which new data is simulated,
# and statistical model is applied to estimate global abundance and trend
# --------------------------------------

for (sim_run in seq(1,1000,1)){
  
  # Scenarios that include additional sampling
  for (n_additional_years in seq(0,40,5)){
    
    n_years_total <- jags.data$n_years + n_additional_years
    
    set.seed(sim_run)
    
    print(sim_run)
    
    if (file.exists("./output/simulation/precision_results.RData")) load(file = "./output/simulation/precision_results.RData")
    if (nrow(precision_results) > 0 & sum(precision_results$sim_run == sim_run & precision_results$n_additional_years == n_additional_years)>0) next
    
    # --------------------------------------
    # - Prepare a data package to use for simulating new data
    # - Use sample sizes (e.g., n_sites, n_years, n_obs, etc) from observed data to
    #   ensure simulations have "realistic" data quantity and balance
    # --------------------------------------
    
    jags.data.sim <- jags.data[c("n_years","n_sites",
                                 "n_obs_aerial","aerial_site","aerial_year","aerial_DoS",
                                 "n_obs_satellite","satellite_site","satellite_year","satellite_DoS","img_qual")]
    
    # --------------------------------------
    # Select a sample from the posterior of the fitted model to specify 'true' parameter values for this simulation
    # --------------------------------------
    
    # Non-joint samples, to represent a wider range of population dynamic scenarios
    jags.data.sim$prob_occ <- out$sims.list$prob_occ[sample(1:out$mcmc.info$n.samples,1)]
    jags.data.sim$logX1_mean <- out$sims.list$logX1_mean[sample(1:out$mcmc.info$n.samples,1)]
    jags.data.sim$logX1_sigma <- out$sims.list$logX1_sigma[sample(1:out$mcmc.info$n.samples,1)]
    
    jags.data.sim$aerial_sigma <- out$sims.list$aerial_sigma[sample(1:out$mcmc.info$n.samples,1)]
    jags.data.sim$sat_slope <- out$sims.list$sat_slope[sample(1:out$mcmc.info$n.samples,1),1:3]
    jags.data.sim$sat_CV <- out$sims.list$sat_CV[sample(1:out$mcmc.info$n.samples,1),1:3]
    jags.data.sim$DoS_slope <- out$sims.list$DoS_slope[sample(1:out$mcmc.info$n.samples,1)]
    
    # --------------------------------------
    # Simulate dynamics at each colony
    # --------------------------------------
    
    # Assume a worst-case scenario where all colonies are declining rapidly enough to reach 50% of initial population size in 3 generations
    r_mean = log(0.5)/(16*3) 
    r_sigma <- out$sims.list$r_sigma[sample(1:out$mcmc.info$n.samples,1)]
    
    # Initial abundances
    logN_matrix <- matrix(NA,nrow = jags.data$n_sites, ncol = 100)
    logN_matrix[,1] <- rnorm(jags.data$n_sites, jags.data.sim$logX1_mean ,jags.data.sim$logX1_sigma)
    
    # Simulate dynamics.  Assume annual population growth is negative if N[t-1] is larger than K
    for (s in 1:jags.data$n_sites){
      
      for (t in 2:ncol(logN_matrix)){
        
        logN_matrix[s,t] <- logN_matrix[s,1] + rnorm(1,r_mean*(t-1),r_sigma)
        
      }
      
    }
    
    logN_matrix <- logN_matrix[,1:n_years_total]
    
    # --------------------------------------
    # Simulation (simulates occupancy process, and observations)
    # --------------------------------------
    
    XX <-  cbind(rep(1,ncol(logN_matrix)),1:ncol(logN_matrix))
    jags.data.sim$regression_weights <- matrix(c(0,1),1,2) %*% solve(t(XX) %*% XX) %*% t(XX)
    jags.data.sim$log_X <- logN_matrix
    jags.data.sim$n_years <- ncol(logN_matrix)
    
    # ***************************************************************************************
    # ***************************************************************************************
    # Evaluate resulting precision of estimates, under different future sampling scenarios
    # ***************************************************************************************
    # ***************************************************************************************
    
    # --------------------------------------------
    # Simulate new aerial surveys in future years
    # --------------------------------------------
    
    new_aerial_obs <- data.frame()
    
    for (s in 1:jags.data$n_sites){
      
      # Calculate proportion of historical years with aerial surveys for this colony
      aer_site <- subset(aer, site_number == s)
      if (nrow(aer_site)==0) next
      
      prob_annual_survey <- length(unique(aer_site$year))/jags.data$n_years
      new_surveys <- expand.grid(aerial_year = 1:ncol(logN_matrix),aerial_site = s)
      new_surveys$surveyed <- rbinom(nrow(new_surveys),1,prob_annual_survey)
      new_surveys <- subset(new_surveys, 
                            aerial_year > jags.data$n_years & 
                              aerial_year <= (jags.data$n_years + n_additional_years) &
                              surveyed == 1)
      
      if (nrow(new_surveys)>0){
        new_surveys$aerial_DoS <- sample(aer_site$day_of_season,nrow(new_surveys),replace=TRUE)
        new_aerial_obs <- rbind(new_aerial_obs,new_surveys)
      }
      
    }
    
    jags.data.sim$aerial_site <- c(jags.data.sim$aerial_site,new_aerial_obs$aerial_site)
    jags.data.sim$aerial_year <- c(jags.data.sim$aerial_year,new_aerial_obs$aerial_year)
    jags.data.sim$aerial_DoS <- c(jags.data.sim$aerial_DoS,new_aerial_obs$aerial_DoS)
    jags.data.sim$n_obs_aerial <- length(jags.data.sim$aerial_DoS)
    
    # --------------------------------------------
    # Simulate new satellite surveys in future years
    # --------------------------------------------
    
    new_satellite_obs <- data.frame()
    
    for (s in 1:jags.data$n_sites){
      
      # Calculate proportion of historical years with satellite surveys for this colony
      sat_site <- subset(sat, site_number == s)
      if (nrow(sat_site)==0) next
      
      prob_annual_survey <- length(unique(sat_site$year))/jags.data$n_years
      new_surveys <- expand.grid(satellite_year = 1:ncol(logN_matrix),satellite_site = s)
      new_surveys$surveyed <- rbinom(nrow(new_surveys),1,prob_annual_survey)
      new_surveys <- subset(new_surveys, 
                            satellite_year > jags.data$n_years & 
                              satellite_year <= (jags.data$n_years + n_additional_years) &
                              surveyed == 1)
      
      if (nrow(new_surveys)>0){
        new_surveys$satellite_DoS <- sample(sat_site$day_of_season,nrow(new_surveys),replace=TRUE)
        new_surveys$img_qual <- sample(sat$img_qualit,nrow(new_surveys),replace=TRUE)
        new_satellite_obs <- rbind(new_satellite_obs,new_surveys)
      }
      
    }
    
    jags.data.sim$satellite_site <- c(jags.data.sim$satellite_site,new_satellite_obs$satellite_site)
    jags.data.sim$satellite_year <- c(jags.data.sim$satellite_year,new_satellite_obs$satellite_year)
    jags.data.sim$satellite_DoS <- c(jags.data.sim$satellite_DoS,new_satellite_obs$satellite_DoS)
    jags.data.sim$img_qual <- c(jags.data.sim$img_qual,new_satellite_obs$img_qual)
    jags.data.sim$n_obs_satellite <- length(jags.data.sim$satellite_DoS)
    
    # --------------------------------------------
    # Simulate new data
    # --------------------------------------------
    
    out_sim <- jags(data=jags.data.sim,
                    model.file="EMPE_model_simulate_observations.jags",
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
      geom_point(aes(x = year, y = y, shape = type))+
      geom_line(data = N_df, aes(x = year, y = N), col = "blue")+
      scale_shape_manual(name = 'Obs type',
                         values =c(19,4))+
      
      facet_wrap(site~., scales = "free")+
      theme_bw()
    sim_data_plot
    
    # Plot global abundance across the simulation
    N_global_df <- data.frame(year = 1:jags.data.sim$n_years,N = out_sim$sims.list$N_global[1,] )
    
    global_N_plot <- ggplot(N_global_df) +
      geom_line(data = N_global_df, aes(x = year, y = N))+
      geom_hline(yintercept = 0, col = "transparent")+
      theme_bw()
    
    global_N_plot 
    
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
      "sat_CV"                  # Coefficient of variation in satellite observations
      
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
    percent_change_est = 100*(out_refit$sims.list$N_global[,jags.data.sim$n_years] - out_refit$sims.list$N_global[,1])/out_refit$sims.list$N_global[,1]
    percent_change_true = 100*(out_sim$sims.list$N_global[,jags.data.sim$n_years] - out_sim$sims.list$N_global[,1])/out_sim$sims.list$N_global[,1]
    
    global_trend_true <- out_sim$sims.list$global_trend
    global_trend_est <- out_refit$sims.list$global_trend
    
    # -----------------------------------------
    # Conduct goodness-of-fit evaluation (Bayesian p-values)
    # -----------------------------------------
    
    Bayesian_pval_adult <- mean(out_refit$sims.list$RMSE_adult_count_actual > out_refit$sims.list$RMSE_adult_count_sim)
    Bayesian_pval_satellite <- mean(out_refit$sims.list$RMSE_satellite_actual > out_refit$sims.list$RMSE_satellite_sim)
    
    # -----------------------------------------
    # Save results for this simulation
    # -----------------------------------------
    
    precision_summary <- data.frame(sim_run = sim_run,
                                    simulation_label = paste0("Sim # ",sim_run),
                                    n_additional_years = n_additional_years,
                                    
                                    r = r_mean, # true slope
                                    
                                    Rhat_hyper = mean(Rhat_hyper > 1.1, na.rm = TRUE),      # Proportion of non-converged hyper parameters
                                    Rhat_indices = mean(Rhat_indices > 1.1, na.rm = TRUE),  # Proportion of non-converged hyper parameters
                                    Bayesian_pval_adult = Bayesian_pval_adult,
                                    Bayesian_pval_satellite = Bayesian_pval_satellite,
                                    
                                    # The actual long-term slope of the population
                                    global_trend_true = out_sim$sims.list$global_trend,
                                    global_trend_est_q025 = quantile(out_refit$sims.list$global_trend,0.025),
                                    global_trend_est_q500 = quantile(out_refit$sims.list$global_trend,0.500),
                                    global_trend_est_q975 = quantile(out_refit$sims.list$global_trend,0.975),
                                    
                                    prob_slope_less_than_zero = mean(out_refit$sims.list$global_trend<0),
                                    
                                    percent_change_est_q025 = quantile(percent_change_est,0.025),
                                    percent_change_est_q500 = quantile(percent_change_est,0.500),
                                    percent_change_est_q975 = quantile(percent_change_est,0.975),
                                    percent_change_true = percent_change_true)
    
    
    # -----------------------------------------
    # Save output
    # -----------------------------------------
    
    load(file = "./output/simulation/precision_results.RData")
    precision_results = rbind(precision_results, precision_summary)
    save(precision_results,file = "./output/simulation/precision_results.RData")
    
    tmp_plot <- ggplot(data = precision_results, aes(x = n_additional_years, y = global_trend_est_q500,ymin = global_trend_est_q025, ymax = global_trend_est_q975))+
      
      geom_hline(yintercept = 0, linewidth=3, col = "gray80", alpha = 0.5)+
      geom_errorbar(width = 0, col = "gray30")+
      geom_point(col = "gray30")+
      geom_point(aes(x = n_additional_years, y = global_trend_true), col = "red", size = 2)+
      
      geom_hline(yintercept = log(0.7)/(16*3) , linetype = 1, linewidth = 1, col = "red", alpha = 0.2)+
      geom_text(x = 0, y = log(0.7)/(16*3), 
                label = "Trend if global population declines by 30% over 3 generations", 
                col = "red",hjust=-0.1,vjust=2,
                fontface = "bold", alpha = 0.2)+
      
      geom_hline(yintercept = mean(precision_results$r), linetype = 1, linewidth = 2, col = "red", alpha = 0.5)+
      geom_text(x = 0, y = mean(precision_results$r), 
                label = "Trend if global population declines by 50% over 3 generations", 
                col = "red",hjust=-0.1,vjust=2,
                fontface = "bold", alpha = 0.5)+
      
      theme_bw()+
      facet_grid(simulation_label~.)+
      ggtitle("Precision analysis for estimates of global trend with additional years of monitoring")+
      xlab("Additional years of monitoring (beyond 2018)")+
      ylab("Estimate of global trend")
    print(tmp_plot)
    
  }
}

