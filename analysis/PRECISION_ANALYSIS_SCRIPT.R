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
# # JAGS script used to simulate data, based on jags.data.sim
# --------------------------------------

sink("EMPE_model_simulate_observations.jags")
cat("
    model {

  # ------------------------------------
  # Simulate aerial observations
  # ------------------------------------
  
  aerial_tau <- pow(aerial_sigma,-2)
  for (i in 1:n_obs_aerial){
    
    # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
    lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] + aerial_DoS[i]*DoS_slope - 0.5*pow(aerial_sigma,2), aerial_tau)
    adult_expected[i] <- z_occ[aerial_site[i],aerial_year[i]] * lambda[i]
    adult_count[i] ~ dpois(adult_expected[i])
    
  }
  
  # ------------------------------------
  # Simulate satellite observations
  # ------------------------------------
  
  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_expected[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]] * exp(DoS_slope * satellite_DoS[i])
    sat_sigma[i] <- sat_expected[i] * sat_CV[img_qual[i]] + 0.00001 # Tiny constant ensures tau is define when var is zero
    sat_tau[i] <- pow(sat_sigma[i],-2)
    satellite[i] ~ dnorm(sat_expected[i], sat_tau[i])

  }
  
}
    ",fill = TRUE)
sink()

# --------------------------------------
# Prepare an empty file to store results of simulations
# --------------------------------------

# Create empty output file to store simulation results
if (!file.exists("./output/simulation/simulation_results.RData")){
  simulation_results <- data.frame()
  save(simulation_results, file = "./output/simulation/simulation_results.RData") 
} else{
  load(file = "./output/simulation/simulation_results.RData")
}

# --------------------------------------
# Conduct repeated simulations in which new data is simulated,
# and the statistical model is applied to estimate global abundance and trend
# --------------------------------------

for (sim_run in seq(1,1000,1)){
  
  # Scenarios that include additional sampling
  for (n_additional_years in seq(0,20,5)){
    
    set.seed(sim_run)
    
    n_years_total <- jags.data$n_years + n_additional_years

    # Skip this iteration, if it has already been run
    if (file.exists("./output/simulation/simulation_results.RData")) load(file = "./output/simulation/simulation_results.RData")
    if (nrow(simulation_results) > 0 & sum(simulation_results$sim_run == sim_run & simulation_results$n_additional_years == n_additional_years)>0) next
    
    print(sim_run)
    
    # --------------------------------------
    # Mean parameter values that will be used in simulations
    # --------------------------------------
    
    # Mean parameters
    prob_occ <- out$mean$prob_occ
    logX1_mean <- out$mean$logX1_mean
    logX1_sigma <- out$mean$logX1_sigma
    r_sigma <- out$mean$r_sigma
    aerial_sigma <- out$mean$aerial_sigma
    sat_slope <- out$mean$sat_slope
    sat_CV <- out$mean$sat_CV
    DoS_slope <- out$mean$DoS_slope
    
    # --------------------------------------
    # Simulate population process that leads to 50% reduction over 3 generations
    # --------------------------------------
    
    target_trend <- log(0.5)/(16*3)  # 50% reduction in 3 generations
    
    # Global population has a tendency to decline
    r_mean_grandmean_mu <- runif(1,-0.2,0)
    r_mean_grandmean_sigma <- 0
    
    # Repeatedly simulate random dynamics, until global population achieves target trend
    
    sim_trend <- NA
    while(is.na(sim_trend)){
      
      # Parameters of simulation
      n_sites <- jags.data$n_sites
      n_years <- n_years_total
      
      jags.data.sim <- jags.data[c("n_obs_aerial","aerial_site","aerial_year","aerial_DoS",
                                   "n_obs_satellite","satellite_site","satellite_year","satellite_DoS","img_qual")]
      
      jags.data.sim$n_sites <- n_sites
      jags.data.sim$n_years <- n_years
      
      # Simulate occupancy states of each colony in each year
      z <- matrix(rbinom(n_years*n_sites,1,prob_occ), ncol = n_years, nrow = n_sites)
      
      # Empty matrices to store abundances
      log_X <- N <- matrix(NA, ncol = n_years, nrow = n_sites)
      
      # Simulate initial abundances for each colony (X[,1])
      log_X[,1] <- rnorm(n_sites, logX1_mean,logX1_sigma)
      
      # Simulate temporal dynamics at each colony
      for (s in 1:n_sites){
        
        r_mean <- rnorm(1,r_mean_grandmean_mu,r_mean_grandmean_sigma) 
        
        for (t in 1:(n_years-1)) log_X[s,t+1] <- rnorm(1,log_X[s,t] + r_mean, r_sigma)
        for (t in 1:n_years) N[s,t] <- exp(log_X[s,t]) * z[s,t]
        
      } 
      
      t_vec <- 1:n_years
      global_N <- colSums(N)
      
      sim_trend <- lm(log(global_N)~t_vec)$coefficients[2]
      if (sim_trend < (target_trend-0.0005) | sim_trend > (target_trend+0.0005)) sim_trend <- NA
    }
    
    # --------------------------------------------
    # Choose years/sites where new aerial observations will occur
    # --------------------------------------------
    
    new_aerial_obs <- data.frame()
    
    for (s in 1:jags.data$n_sites){
      
      # Calculate proportion of historical years with aerial surveys for this colony
      aer_site <- subset(aer, site_number == s)
      if (nrow(aer_site)==0) next
      
      prob_annual_survey <- length(unique(aer_site$year))/jags.data$n_years
      new_surveys <- expand.grid(aerial_year = 1:n_years_total,aerial_site = s)
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
    # Choose years/sites where new satellite observations will occur
    # --------------------------------------------
    
    new_satellite_obs <- data.frame()
    
    for (s in 1:jags.data$n_sites){
      
      # Calculate proportion of historical years with satellite surveys for this colony
      sat_site <- subset(sat, site_number == s)
      if (nrow(sat_site)==0) next
      
      prob_annual_survey <- length(unique(sat_site$year))/jags.data$n_years
      new_surveys <- expand.grid(satellite_year = 1:n_years_total,satellite_site = s)
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
    
    # ------------------------------------
    # Use JAGS to simulate aerial and satellite observations
    # ------------------------------------
    
    jags.data.sim$N <- N
    jags.data.sim$log_X <- log_X
    jags.data.sim$z_occ <- z
    jags.data.sim$aerial_sigma <- aerial_sigma
    jags.data.sim$sat_CV <- sat_CV
    jags.data.sim$sat_slope <- sat_slope
    jags.data.sim$DoS_slope <- DoS_slope
    
    # Simulate data
    out_sim <- jags(data=jags.data.sim,
                    model.file="EMPE_model_simulate_observations.jags",
                    parameters.to.save=c("adult_count","satellite"),
                    n.chains=1,n.thin = 1,n.iter= 2,n.burnin= 1)
    
    # --------------------------------------
    # Fit model to simulated data
    # --------------------------------------
    
    # Data to fit
    jags.data.refit = list( n_years = jags.data.sim$n_years,
                            n_sites = jags.data.sim$n_sites,
                            
                            # aerial counts of adults
                            n_obs_aerial = length(out_sim$sims.list$adult_count[1,]),
                            adult_count = out_sim$sims.list$adult_count[1,],
                            aerial_site = jags.data.sim$aerial_site,
                            aerial_year = jags.data.sim$aerial_year,
                            aerial_DoS = jags.data.sim$aerial_DoS,
                            
                            # satellite counts
                            n_obs_satellite = length(out_sim$sims.list$satellite[1,]),
                            satellite = out_sim$sims.list$satellite[1,],
                            satellite_site = jags.data.sim$satellite_site,
                            satellite_year = jags.data.sim$satellite_year,
                            satellite_DoS = jags.data.sim$satellite_DoS,
                            img_qual = jags.data.sim$img_qual
    )
    
    # For calculating global trend in JAGS
    XX <-  cbind(rep(1,n_years_total),1:n_years_total)
    jags.data.refit$regression_weights <- matrix(c(0,1),1,2) %*% solve(t(XX) %*% XX) %*% t(XX)
    
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
                      n.iter= 10000,
                      n.burnin= 5000,
                      parallel = TRUE)
    
    #----------------------------------------------------------
    # Examine model convergence
    #----------------------------------------------------------
    
    hyper_parameters <- c(
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
    
    indices_and_trends <- c(
      # Global abundance each year
      "N_global",
      # Log-linear OLS slope across the study
      "global_trend",
      # Colony-specific mean annual differences
      "r_mean",
      # Colony-specific abundance each year
      "N")
    Rhat_indices <- unlist(out_refit$Rhat[which(names(out_refit$Rhat) %in% indices_and_trends)])

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
    
    simulation_summary <- data.frame(sim_run = sim_run,
                                     simulation_label = paste0("Sim # ",sim_run),
                                     n_additional_years = n_additional_years,
                                     
                                     Rhat_hyper = mean(Rhat_hyper > 1.1, na.rm = TRUE),      # Proportion of non-converged hyper parameters
                                     Rhat_indices = mean(Rhat_indices > 1.1, na.rm = TRUE),  # Proportion of non-converged hyper parameters
                                     Bayesian_pval_adult = Bayesian_pval_adult,
                                     Bayesian_pval_satellite = Bayesian_pval_satellite,
                                     
                                     global_trend_est_mean = mean(out_refit$sims.list$global_trend),
                                     global_trend_est_q050 = quantile(out_refit$sims.list$global_trend,0.050),
                                     global_trend_est_q500 = quantile(out_refit$sims.list$global_trend,0.500),
                                     global_trend_est_q950 = quantile(out_refit$sims.list$global_trend,0.950),
                                     global_trend_true = sim_trend)
    
    
    # -----------------------------------------
    # Save output
    # -----------------------------------------
    
    load(file = "./output/simulation/simulation_results.RData")
    simulation_results = rbind(simulation_results, simulation_summary)
    save(simulation_results,file = "./output/simulation/simulation_results.RData")
    
    # ----------------------------------------------------
    # Plot mean across repeated runs
    # ----------------------------------------------------
    
    simulation_results$coverage <- simulation_results$global_trend_est_q050 <= simulation_results$global_trend_true & simulation_results$global_trend_est_q950 >= simulation_results$global_trend_true
    
    results_mean <- simulation_results %>%
      subset(Rhat_hyper <= 0.001 & Rhat_indices <= 0.001) %>%
      group_by(n_additional_years) %>%
      summarize_all(mean)
    
    plot <- ggplot(data = results_mean, aes(x = n_additional_years, y = global_trend_est_q500,ymin = global_trend_est_q050, ymax = global_trend_est_q950))+
      
      geom_hline(yintercept = 0, linewidth=3, col = "gray80", alpha = 0.5)+
      geom_errorbar(width = 0, col = "gray30")+
      geom_point(col = "gray30")+
      geom_point(aes(x = n_additional_years, y = global_trend_true), col = "red", size = 2)+
      
      geom_hline(yintercept = log(0.5)/(16*3) , linetype = 1, linewidth = 2, col = "red", alpha = 0.5)+
      geom_text(x = 0, y = log(0.5)/(16*3) , 
                label = "Trend if global population declines by 50% over 3 generations", 
                col = "red",hjust=-0.1,vjust=2,
                fontface = "bold", alpha = 0.5)+
      
      theme_bw()+
      ggtitle("Precision analysis for estimates of global trend with additional years of monitoring")+
      xlab("Additional years of monitoring (beyond 2018)")+
      ylab("Global trend estimate")
    
    print(plot)
  }
}

# ----------------------------------------------------
# final plot
# ----------------------------------------------------

simulation_results$coverage <- simulation_results$global_trend_est_q050 <= simulation_results$global_trend_true & simulation_results$global_trend_est_q950 >= simulation_results$global_trend_true

results_mean <- simulation_results %>%
  subset(Rhat_hyper <= 0.001 & Rhat_indices <= 0.001) %>%
  group_by(n_additional_years) %>%
  summarize_all(mean)

plot <- ggplot(data = results_mean, aes(x = n_additional_years, y = global_trend_est_q500,ymin = global_trend_est_q050, ymax = global_trend_est_q950))+
  
  geom_hline(yintercept = 0, linewidth=3, col = "gray80", alpha = 0.5)+
  geom_errorbar(width = 0, col = "gray30")+
  geom_point(col = "gray30")+
  #geom_point(aes(x = n_additional_years, y = global_trend_true), col = "red", size = 2)+
  
  geom_hline(yintercept = log(0.7)/(16*3) , 
             linetype = 1, linewidth = 2, col = "orange", alpha = 0.5)+
  geom_text(x = 0, y = log(0.7)/(16*3) , 
            label = "Trend if global population declines by 30% over 3 generations", 
            col = "orange",hjust=-0.01,vjust=2,
            fontface = "bold", alpha = 0.5)+
  
  geom_hline(yintercept = log(0.5)/(16*3) , linetype = 1, linewidth = 2, col = "red", alpha = 0.5)+
  geom_text(x = 0, y = log(0.5)/(16*3) , 
            label = "Trend if global population declines by 50% over 3 generations", 
            col = "red",hjust=-0.01,vjust=2,
            fontface = "bold", alpha = 0.5)+
  
  theme_bw()+
  ggtitle("Precision analysis for global trend with additional years of monitoring\n
          (simulations of 50% decrease over 3 generations)")+
  xlab("Additional years of monitoring (beyond 2018)")+
  ylab("Global trend estimate\n(log-linear slope)\n")

print(plot)
