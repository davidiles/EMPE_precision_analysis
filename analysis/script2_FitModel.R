# install/load necessary packages
my.packs <- c(
  
  # Bayesian analysis
  'jagsUI','MCMCvis','DHARMa',
  
  # Data manipulation
  'tidyverse','reshape2','lubridate','readxl',
  
  # Spatial analysis
  'rgeos','raster','sp','sf',
  
  # Plotting
  'ggrepel','ggthemes','ggpubr','scales','viridis')

if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/EMPE_Global/analysis")

rm(list=ls())

# **************************************************************************************************
# **************************************************************************************************
# FIT MODEL
# **************************************************************************************************
# **************************************************************************************************

load("output/EMPE_data_formatted.RData")

# The jags script to fit the model
sink("EMPE_model_empirical.jags")
cat("
    model {
  
  # ------------------------------------
  # Population process model
  # ------------------------------------
  
  # Spatial random effect for X[s,1]
  X1_mean ~ dunif(0,50000)
  logX1_mean <- log(X1_mean)
  logX1_tau ~ dgamma(0.001,0.001)
  logX1_sigma <- 1/pow(logX1_tau,0.5)
  
  # Spatial random effect for r_mean[s]
  r_mean_grandmean_mu ~ dnorm(0,25)
  r_mean_grandmean_tau ~ dgamma(0.001,0.001)
  r_mean_grandmean_sigma <- 1/pow(r_mean_grandmean_tau,0.5)
  
  # Temporal variance in population growth
  r_tau ~ dgamma(0.001,0.001)
  r_sigma <- 1/pow(r_tau,0.5)
  
  # Growth dynamics at each colony
  for (s in 1:n_sites){
    
    # Initial abundance at colony 
    log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
    
    # Median annual change at colony
    r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau) 
    
    # Population growth from year to year at colony
    for (t in 1:(n_years-1)){
      log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s], r_tau)
    }
  } 
  
  # Probability each colony is occupied each year
  prob_occ ~ dunif(0,1)
  
  # Occupancy dynamics at each colony (not a Markov process)
  for (s in 1:n_sites){
    for (t in 1:n_years){
    
      # True occupancy state at colony in each year
      z_occ[s,t] ~ dbern(prob_occ)
      
    }
  }
  
  # ------------------------------------
  # Annual population indices at each colony
  # ------------------------------------
  
  for (s in 1:n_sites){
    for (t in 1:n_years){
      N[s,t] <- exp(log_X[s,t]) * z_occ[s,t]
    }
  }
  
  # ------------------------------------
  # Shared 'day of season' effect to account for phenology
  # ------------------------------------
  DoS_slope ~ dnorm(0,25)

  # ------------------------------------
  # Aerial observation model
  # ------------------------------------
  
  aerial_tau ~ dgamma(0.001,0.001) 
  aerial_sigma <- 1/pow(aerial_tau,0.5)
  
  for (i in 1:n_obs_aerial){
    
    lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] + aerial_DoS[i]*DoS_slope - 0.5*pow(aerial_sigma,2), aerial_tau)
    adult_expected[i] <- z_occ[aerial_site[i],aerial_year[i]] * lambda[i]
    adult_count[i] ~ dpois(adult_expected[i])
    
  }
  
  # ------------------------------------
  # Satellite observation model
  # ------------------------------------
  
  # Describes proportional bias (slope) and variance in satellite counts for each level of image quality (1, 2, or 3)
  for (i in 1:3){
    sat_slope[i] ~ dnorm(1,25)
    sat_CV[i] ~ dunif(0,2)
  }

  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_expected[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]] * exp(DoS_slope * satellite_DoS[i])
    sat_sigma[i] <- sat_expected[i] * sat_CV[img_qual[i]] + 0.00001 # Tiny constant ensures tau is define when var is zero
    sat_tau[i] <- pow(sat_sigma[i],-2)
    
    satellite[i] ~ dnorm(sat_expected[i], sat_tau[i])
  }
  
  # ------------------------------------
  # Global change and trend
  # ------------------------------------
  
  # Global abundance each year
  for (t in 1:n_years){
    N_global[t] <- sum(N[1:n_sites,t])
  }
  
  # Global trend (measured as least squares regression slope on log scale)
  global_trend <- inprod(log(N_global[1:n_years]),regression_weights[1,1:n_years])
  
  #------------------------------------
  # Posterior predictive checks
  #------------------------------------
  
  # Posterior predictive check for aerial count data
  for (i in 1:n_obs_aerial){
  
    # Simulate aerial observations
    sim_lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] + aerial_DoS[i]*DoS_slope - 0.5*pow(aerial_sigma,2), aerial_tau)
    sim_adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * sim_lambda[i])
    
    # Calculate discrepancy measures of actual and simulated data
    sqE_adult_count_actual[i] <- pow(adult_count[i] - N[aerial_site[i],aerial_year[i]] ,2)
    sqE_adult_count_sim[i] <- pow(sim_adult_count[i] - N[aerial_site[i],aerial_year[i]],2)
    
  }
  
  # Posterior predictive check for satellite data
  for (i in 1:n_obs_satellite){
  
    # Simulate satellite observations
    sim_satellite[i] ~ dnorm(sat_expected[i],sat_tau[i])
    
    # Calculate discrepancy measures of actual and simulated data
    sqE_satellite_actual[i] <- pow(satellite[i] - N[satellite_site[i],satellite_year[i]],2)
    sqE_satellite_sim[i] <- pow(sim_satellite[i] - N[satellite_site[i],satellite_year[i]],2)
    
  }
  
  # Root mean squared error of empirical data
  RMSE_adult_count_actual <- sqrt(mean(sqE_adult_count_actual[]))
  RMSE_satellite_actual <- sqrt(mean(sqE_satellite_actual[]))
  
  # Root mean squared error of simulated data
  RMSE_adult_count_sim <- sqrt(mean(sqE_adult_count_sim[]))
  RMSE_satellite_sim <- sqrt(mean(sqE_satellite_sim[]))
  
}
    ",fill = TRUE)
sink()

# Generate initial values
inits <- function(){list(r_mean = rnorm(jags.data$n_sites,0,0.03),
                         z_occ = matrix(1,ncol = n_years, nrow = n_sites)
)}

# out <- jags(data=jags.data,
#             model.file="EMPE_model_empirical.jags",
#             parameters.to.save=c(
#               
#               # ------------------------
#               # Hyper-parameters
#               # ------------------------
#               "prob_occ",                # Probability colonies are "occupied"
#               "r_mean_grandmean_mu",     # Hierarchical grand mean of colony-level annual growth rates
#               "r_mean_grandmean_sigma",  # Hierarchical sd of colony-level growth rates
#               "logX1_mean",              # Hierarchical grand mean of colony-level initial abundances
#               "logX1_sigma",             # Hierarchical sd of colony-level initial abundances
#               "r_sigma",                 # Temporal variance of colony-level annual growth rates
#               
#               "DoS_slope",
#               "aerial_sigma",            # SD of aerial observations (on log scale)
#               "sat_slope",               # Bias in satellite observations
#               "sat_CV",                  # Coefficient of variation in satellite observations
#               
#               # ------------------------
#               # Colony-specific quantities
#               # ------------------------
#               
#               # Colony-specific mean annual differences (could be used as a measure of colony-specific "trend")
#               "r_mean",
#               
#               # Colony-specific index of abundance each year
#               "N",
#               
#               # ------------------------
#               # Global index and trend estimates
#               # ------------------------
#               
#               # Global index of abundance each year
#               "N_global",
#               
#               # Log-linear OLS slope across the study
#               "global_trend",
#               
#               # ------------------------
#               # Observation-specific quantities
#               # ------------------------
#               
#               # Fitted values
#               "adult_expected",
#               "sat_expected",
#               "sat_sigma",
#               
#               # Discrepancy measures for posterior predictive checks
#               "RMSE_adult_count_actual",
#               "RMSE_satellite_actual",
#               "RMSE_adult_count_sim",
#               "RMSE_satellite_sim",
#               
#               # Simulated datasets from fitted model (also for goodness-of-fit testing)
#               "sim_adult_count",
#               "sim_satellite"
#               
#             ),
#             
#             inits = inits,
#             n.chains=3,
#             n.thin = 10*5,
#             n.iter= 110000*5,
#             n.burnin= 10000*5,
#             parallel = TRUE)
# 
# save(out, file = "output/fitted_model.RData")

# **************************************************************************************************
# **************************************************************************************************
# MODEL CONVERGENCE AND GOODNESS-OF-FIT
# **************************************************************************************************
# **************************************************************************************************

load(file = "output/fitted_model.RData")

# ----------------------------------------------------------
# Examine model convergence
# ----------------------------------------------------------

# Convergence of relevant hyper-parameters in model
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

Rhats <- unlist(out$Rhat[which(names(out$Rhat) %in% hyper_parameters)])
mean(Rhats > 1.1, na.rm = TRUE) # Proportion of parameters with Rhat > 1.1
max(Rhats, na.rm=TRUE) # max Rhat
Rhats[which(Rhats > 1.1)]
MCMCtrace(out, params = hyper_parameters, Rhat = TRUE, filename = "output/model_checks/traceplot_hyperparams.pdf")

# Convergence of relevant population indices and derived parameters
indices_and_trends <- c(
  
  # ------------------------
  # Global estimates
  # ------------------------
  
  "N_global",     # Global abundance each year
  "global_trend", # Log-linear OLS slope across the study
  
  # ------------------------
  # Colony-level estimates
  # ------------------------

  "r_mean", # Colony-specific mean annual differences
  "N"       # Colony-specific abundance each year
  
)

Rhats <- unlist(out$Rhat[which(names(out$Rhat) %in% indices_and_trends)])
mean(Rhats > 1.1, na.rm = TRUE) # Proportion of parameters with Rhat > 1.1
max(Rhats, na.rm=TRUE) # max Rhat
Rhats[which(Rhats > 1.1)] # N125     N353 
MCMCtrace(out, params = indices_and_trends, Rhat = TRUE, filename = "output/model_checks/traceplot_indices.pdf")

# Effective sample sizes
n.eff <- unlist(out$n.eff)
n.eff[n.eff > 1 & n.eff <= 2500] # Parameters with fewer than 2500 samples

# ----------------------------------------------------------
# Posterior predictive checks
# ----------------------------------------------------------

pval_adult <- mean(out$sims.list$RMSE_adult_count_actual > out$sims.list$RMSE_adult_count_sim) %>% round(2)
lim <- range(c(out$sims.list$RMSE_adult_count_actual,out$sims.list$RMSE_adult_count_sim))
df_aerial = data.frame(actual = out$sims.list$RMSE_adult_count_actual,
                       simulated = out$sims.list$RMSE_adult_count_sim)
plot1 <- ggplot(data = df_aerial,
                aes(x = actual, y = simulated )) +
  geom_hex(binwidth = diff(lim)/50) +
  scale_fill_gradientn(colors = c("gray95","darkblue"), name = "Number of\nsimulated\ndatasets") +
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim,ylim=lim)+
  ggtitle(paste0("Posterior predictive check: Adult counts\n\nBayesian p-value = ",pval_adult))+
  xlab("RMSE adult counts (actual)")+
  ylab("RMSE adult counts (simulated)")+
  theme_bw()

pval_satellite <- mean(out$sims.list$RMSE_satellite_actual > out$sims.list$RMSE_satellite_sim) %>% round(2)
lim <- range(c(out$sims.list$RMSE_satellite_actual,out$sims.list$RMSE_satellite_sim))
df_sat = data.frame(actual = out$sims.list$RMSE_satellite_actual,
                    simulated = out$sims.list$RMSE_satellite_sim)

plot2 <- ggplot(data = df_sat,
                aes(x = actual, y = simulated )) +
  geom_hex(binwidth = diff(lim)/50) +
  scale_fill_gradientn(colors = c("gray95","darkblue"), name = "Number of\nsimulated\ndatasets") +
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim,ylim=lim)+
  ggtitle(paste0("Posterior predictive check: Satellite surveys\n\nBayesian p-value = ",pval_satellite))+
  xlab("RMSE satellite counts (actual)")+
  ylab("RMSE satellite counts (simulated)")+
  theme_bw()

pval_plot <- ggarrange(plot1,plot2,nrow=2)

png("./output/model_checks/Bayesian_pval_plot.png", width = 5, height = 7, units = "in",res=500)
print(pval_plot)
dev.off()

# -----------------------------------------------------------------
# Use DHARMa for residual diagnostics
# -----------------------------------------------------------------
# NOTE: these were not especially useful because shrinkage of random effects in model
#       leads to pattern in residual vs predicted (see discussion of this on https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html)
#      - note some evidence of *under*dispersion in counts, likely caused by within-season autocorrelation
sim_aerial = createDHARMa(simulatedResponse = t(out$sims.list$sim_adult_count), # columns are posterior samples
                          observedResponse = jags.data$adult_count,
                          fittedPredictedResponse = apply(out$sims.list$adult_expected, 2, median),
                          integerResponse = T)

png("./output/model_checks/DHARMa_AdultCounts.png", width = 8, height = 5, units = "in",res=500)
plot(sim_aerial)
dev.off()

resid_aerial <- data.frame(site_id = aer$site_id,
                           resid = residuals(sim_aerial),
                           adult_count = aer$adult_count,
                           Date = aer$Date,
                           Year = aer$year,
                           yday = aer$yday)

sim_sat = createDHARMa(simulatedResponse = t(out$sims.list$sim_satellite),
                       observedResponse = jags.data$satellite,
                       fittedPredictedResponse = apply(out$sims.list$sat_expected, 2, median),
                       integerResponse = T)

png("./output/model_checks/DHARMa_Satellite.png", width = 8, height = 5, units = "in",res=500)
plot(sim_sat)
dev.off()

# **************************************************************************************************
# **************************************************************************************************
# MODEL INTERPRETATION
# **************************************************************************************************
# **************************************************************************************************

# ----------------------------------------------------------
# Output a table of parameter estimates
# ----------------------------------------------------------

parameter_estimates = out$summary[1:which(rownames(out$summary) == "sat_CV[3]"),] %>% as.data.frame()
write.csv(parameter_estimates, file = "output/model_results/parameter_estimates.csv", row.names = TRUE)

# ----------------------------------------------------------
# Function to calculate estimate of 'overall change' between endpoints (2009 and 2018) and
#   log-linear slope (a measure of trend) across study
# ----------------------------------------------------------

# Accepts an n_samp x n_year matrix of abundance estimates
change_trend_fn <- function(mat){
  
  n_samps = nrow(mat)
  n_years = ncol(mat)
  
  # Percent change between endpoints
  percent_change_samples <- 100*(mat[,n_years] - mat[,1])/mat[,1]
  percent_change_summary <- c(mean = mean(percent_change_samples),SE = sd(percent_change_samples), quantile(percent_change_samples,c(0.025,0.05,0.5,0.95,0.975)))
  
  prob_decline <- mean(percent_change_samples < 0)
  prob_30percent_decline <- mean(percent_change_samples < -30)
  prob_50percent_decline <- mean(percent_change_samples < -50)
  
  # Rate of log-linear change (calculated using OLS)
  XX = cbind(rep(1,n_years),1:n_years)
  regression_weights <- matrix(c(0,1),1,2) %*% solve(t(XX) %*% XX) %*% t(XX)
  OLS_regression_samples <- rep(NA,n_samps)
  for (i in 1:n_samps) OLS_regression_samples[i] <- regression_weights %*% log(mat[i,]) # OLS slope, on log scale
  
  OLS_regression_summary <- c(mean = mean(OLS_regression_samples),
                              SE = sd(OLS_regression_samples),
                              quantile(OLS_regression_samples,c(0.025,0.05,0.5,0.95,0.975)),
                              prob_negative = mean(OLS_regression_samples<0))
  
  
  return(list(prob_decline = prob_decline, 
              prob_30percent_decline = prob_30percent_decline,
              prob_50percent_decline = prob_50percent_decline,
              percent_change_samples = percent_change_samples,
              OLS_regression_samples = OLS_regression_samples,
              
              percent_change_summary = percent_change_summary,
              OLS_regression_summary = OLS_regression_summary))
}

# ----------------------------------------------------------
# Convert colony-level abundance mcmc samples to dataframe
# ----------------------------------------------------------

N_samples <- out$sims.list$N %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value)

N_samples$year <- year_range[N_samples$year_number]

#----------------------------------------------------------
# Summarize global dynamics
#----------------------------------------------------------

N_global <- out$sims.list$N_global 
colnames(N_global) <- year_range
N_global <- N_global %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, Year = Var2, N = value)

N_global_matrix = N_global %>% 
  spread(Year, N) %>%
  dplyr::select(-mcmc_sample) %>%
  as.matrix()

# Estimates of abundance/change/trend
global_abundance_summary <- N_global %>%
  group_by(Year) %>%
  summarize(
    N_mean = mean(N),
    N_se = sd(N),
    N_q025 = quantile(N,0.025),
    N_q05 = quantile(N,0.05),
    N_median = quantile(N,0.5),
    N_q95 = quantile(N,0.95),
    N_q975 = quantile(N,0.975))

global_change_trend <- change_trend_fn(N_global_matrix)

# Percent change from 2009 to 2018
global_change_summary <- data.frame(Region = "Global",
                                    Prob_Decline = global_change_trend$prob_decline,
                                    Prob_30percent_Decline = global_change_trend$prob_30percent_decline,
                                    Prob_50percent_Decline = global_change_trend$prob_50percent_decline,
                                    Quantile = names(global_change_trend$percent_change_summary),
                                    Estimate = global_change_trend$percent_change_summary) %>% 
  spread(Quantile, Estimate)

global_plot <- ggplot(global_abundance_summary, aes(x = Year, y = N_mean, ymin = N_q025, ymax = N_q975))+
  geom_ribbon(fill = "#0071fe", alpha = 0.3)+
  geom_line(col = "#0071fe")+
  ylab("Index of abundance")+
  xlab("Year")+
  #ggtitle("Global population")+
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.25, 0.25)))+
  theme_few()

tiff(filename = "output/model_results/1_Global_Level/global_trajectory.tif", width = 5, height = 4, units = "in", res = 300)
print(global_plot)
dev.off()

# ----------------------------------------------------------
# Contextualize global log-linear trend estimate
# ----------------------------------------------------------

# Log-linear slope from 2009 to 2018
global_trend_summary <- data.frame(Region = "Global",
                                   Quantile = names(global_change_trend$OLS_regression_summary),
                                   Estimate = global_change_trend$OLS_regression_summary) %>% 
  spread(Quantile, Estimate)

# Estimate of log-linear slope
100*(exp(global_trend_summary[,c("2.5%","50%","97.5%")])-1)

# Generation time = 16 years (Jenouvrier et al. 2014 Nature Climate Change); 3 generations = 48 years

# To achieve a 30% decline over 3 generations, log-linear trend must be at least: log(0.7) = log(1) + x*(GL*3)
# x = log(0.7)/(GL*3)
GL = 16
x = log(0.7)/(GL*3)

# Probability slope is more negative than x:
mean(global_change_trend$OLS_regression_samples <= x) # 0.69

# Probability global log-linear trend is negative
mean(global_change_trend$OLS_regression_samples<0) # 0.87

# Credible interval width
global_change_summary$`97.5%`-global_change_summary$`2.5%`  # width of 95% CI = 40.6
global_change_trend$OLS_regression_summary["97.5%"] - global_change_trend$OLS_regression_summary["2.5%"] # 0.044

# Percent change in population size between 2009 and 2018
global_change_summary$`50%` # -9.6
global_change_summary$Prob_Decline # 0.81

quantile(N_global_matrix[,10] - N_global_matrix[,1],c(0.025,0.5,0.975))

write.csv(global_abundance_summary, file = "output/model_results/1_Global_Level/GLOBAL_abundance.csv", row.names = FALSE)
write.csv(global_change_summary, file = "output/model_results/1_Global_Level/GLOBAL_change.csv", row.names = FALSE)
write.csv(global_trend_summary, file = "output/model_results/1_Global_Level/GLOBAL_trend.csv", row.names = FALSE)

# ----------------------------------------------------------
# Calculate change since 2009 at each colony
# ----------------------------------------------------------

delta_N_samples = data.frame()
N_2009 = out$sims.list$N[,,1]

for (t in 1:jags.data$n_years){
  
  # Abundance in this year
  N_t = out$sims.list$N[,,t]
  
  # Change in abundance relative to 2009
  delta_N_t = N_t - N_2009
  
  # Long format
  delta_N_t_samples <- delta_N_t %>% 
    reshape2::melt() %>% 
    rename(mcmc_sample = Var1, site_number = Var2, delta_N = value)
  delta_N_t_samples$year <- year_range[t]
  delta_N_samples = rbind(delta_N_samples, delta_N_t_samples)
  
}

# Join with N_samples dataframe
N_samples = full_join(N_samples, delta_N_samples)

# ----------------------------------------------------------
# Summarize dynamics at each colony
# ----------------------------------------------------------

colony_summary = N_samples %>%
  group_by(year, site_number) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975),
            
            change_since_2009_mean = mean(delta_N),
            change_since_2009_q025 = quantile(delta_N,0.025),
            change_since_2009_q05 = quantile(delta_N,0.05),
            change_since_2009_median = median(delta_N),
            change_since_2009_q95 = quantile(delta_N,0.95),
            change_since_2009_q975 = quantile(delta_N,0.975),
            
            prob_decline_since_2009 = mean(delta_N < 0)) %>%
  
  left_join(colony_attributes)

write.csv(colony_summary, file = "output/model_results/3_Colony_Level/colony_summary.csv", row.names = FALSE)

# ----------------------------------------------------------
# Plot dynamics within each colony
# ----------------------------------------------------------

aer$img_qualit = 3
p1 <- ggplot() +
  
  # Fitted dynamics
  geom_ribbon(data = colony_summary, aes(x = year, y = N_median, ymin = N_q05, ymax = N_q95),fill = "#0071fe", alpha = 0.3)+
  geom_line(data = colony_summary, aes(x = year, y = N_median),col = "#0071fe")+
  
  # Satellite observations
  geom_vline(data = subset(sat, area_m2 == 0), aes(xintercept = year), col = "gray80", size=2)+
  geom_point(data = sat, aes(x = year, y = area_m2, shape = "Satellite count (area in m2)", col = yday,
                             size = factor(img_qualit)), stroke = 0.5)+
  
  # Aerial observations
  geom_point(data = aer, aes(x = year, y = adult_count, shape = "Aerial count (adult)", col = yday,size = factor(img_qualit)))+
  
  # Scales
  scale_shape_manual(name = 'Survey Type', values =c('Satellite count (area in m2)'=10,'Aerial count (adult)'= 17))+
  scale_color_gradientn(name = 'Date of survey\n(day of year)', colors = magma(10)[1:9])+
  scale_size_manual(values = c(0.5,1,2), name = "Image Quality",
                    labels = c("1 - Poor","2 - Moderate","3 - Good"))+
  scale_x_continuous(limits = range(year_range))+
  
  ylab("Count")+
  xlab("Year")+
  facet_wrap(site_id_factor~., scales = "free",nrow=10,ncol=5)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, col = "transparent")

png("output/model_results/3_Colony_Level/colony_dynamics_fitted.png", units = "in", res = 300, width = 15, height = 20)
print(p1)
dev.off()

pdf("output/model_results/3_Colony_Level/colony_dynamics_fitted.pdf", width = 10, height = 30)
print(p1)
dev.off()


# ----------------------------------------------------------
# Summarize regional dynamics
# ----------------------------------------------------------

# ----------------------------------------------------------
# Function to calculate regional indices, and change/trend estimates for
#    any grouping of colonies (e.g., based on fast ice regions, ccamlr sectors,
#    colonies inside/outside protected areas, etc)
# ----------------------------------------------------------

# Accepts an n_samp x n_year matrix of abundance estimates
regional_estimate_fn <- function(region_names = NA, 
                                 region_colors = NA,
                                 N_samples = N_samples){
  
  tmp <- colony_attributes[,c("site_id","site_name","site_number")]
  tmp$region <- colony_attributes[,region_names]
  
  N_region <- N_samples %>%
    left_join(tmp)%>%
    group_by(region,mcmc_sample,year) %>%
    summarize(N = sum(N))
  
  regional_abundance_summary <- N_region %>%
    group_by(region,year) %>%
    summarize(N_mean = mean(N),
              N_se = sd(N),
              N_q025 = quantile(N,0.025),
              N_q05 = quantile(N,0.05),
              N_median = median(N),
              N_q95 = quantile(N,0.95),
              N_q975 = quantile(N,0.975))
  
  n_reg <- length(unique(N_region$region))
  regional_trend_summary <- data.frame()
  regional_change_summary <- data.frame()
  regional_trend_samples <- vector(mode = "list", length = 0)
  regional_change_samples <- vector(mode = "list", length = 0)
  
  for (i in 1:n_reg){
    
    reg <- unique(N_region$region)[i]
    N_reg_matrix = N_region %>%
      subset(region == reg) %>%
      spread(year, N) %>%
      ungroup() %>%
      dplyr::select(-region,-mcmc_sample) %>%
      as.matrix()
    
    # Estimates of change/trend
    reg_change_trend <- change_trend_fn(N_reg_matrix)
    
    regional_change_samples[[i]] <- reg_change_trend$percent_change_samples
    regional_trend_samples[[i]] <- reg_change_trend$OLS_regression_samples
    
    regional_change_summary <- rbind(regional_change_summary,
                                     data.frame(Region = reg,
                                                Prob_Decline = reg_change_trend$prob_decline,
                                                Prob_30percent_Decline = reg_change_trend$prob_30percent_decline,
                                                Prob_50percent_Decline = reg_change_trend$prob_50percent_decline,
                                                Quantile = names(reg_change_trend$percent_change_summary),
                                                Estimate = reg_change_trend$percent_change_summary) %>%
                                       spread(Quantile, Estimate))
    
    regional_trend_summary <- rbind(regional_trend_summary,
                                    data.frame(Region = reg,
                                               Quantile = names(reg_change_trend$OLS_regression_summary),
                                               Estimate = reg_change_trend$OLS_regression_summary) %>%
                                      spread(Quantile, Estimate))
    
    
  }
  
  # --------------------------------
  # Save summary tables
  # --------------------------------
  
  write.csv(regional_abundance_summary, file = paste0("output/model_results/2_Regional_Level/summary_REGIONAL_",region_names,"abundance.csv"), row.names = FALSE)
  write.csv(regional_change_summary, file = paste0("output/model_results/2_Regional_Level/summary_REGIONAL_",region_names,"_change.csv"), row.names = FALSE)
  write.csv(regional_trend_summary, file = paste0("output/model_results/2_Regional_Level/summary_REGIONAL_",region_names,"_trend.csv"), row.names = FALSE)
  
  # --------------------------------
  # Generate separate plots for each region
  # --------------------------------
  region_plots <- vector(mode = "list", length = 0)
  
  for (i in 1:n_reg){
    
    reg <- unique(N_region$region)[i]
    
    regional_abundance_summary$region <- factor(regional_abundance_summary$region, levels = names(region_colors))
    
    fill_col <- region_colors[reg]
    reg_plot1 <- ggplot()+
      
      geom_ribbon(data = subset(regional_abundance_summary, region == reg),
                  aes(x = year, y = N_mean, ymin = N_q025, ymax = N_q975),
                  fill = fill_col,alpha = 0.3)+
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_median),
                col = "black")+
      
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_q025),
                col = "black")+
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_q975),
                col = "black")+
      ylab("Index of abundance")+
      xlab("Year")+
      ggtitle(reg)+
      guides(fill="none", col = "none")+
      theme_few()
    region_plots[[i]] <- reg_plot1
    
    lim <- c(0,max(regional_abundance_summary$N_q975))
    reg_plot2 <- ggplot()+
      geom_ribbon(data = subset(regional_abundance_summary, region == reg),
                  aes(x = year, y = N_mean, ymin = N_q025, ymax = N_q975),
                  fill = fill_col,alpha = 0.3)+
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_median),
                col = "black")+
      
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_q025),
                col = "black")+
      geom_line(data = subset(regional_abundance_summary, region == reg),
                aes(x = year, y = N_q975),
                col = "black")+
      
      ylab("Index of abundance")+
      xlab("Year")+
      ggtitle(reg)+
      guides(fill="none", col = "none")+
      scale_y_continuous(limits = lim)+
      theme_few()
    
    # --------------------------------
    # Save figures
    # --------------------------------
    
    tiff(filename = paste0("output/model_results/2_Regional_Level/REGIONAL_",region_names,"_CustomYAxis_",reg,".tif"), width = 4, height = 3, units = "in", res = 300)
    print(reg_plot1)
    dev.off()
    
    tiff(filename = paste0("output/model_results/2_Regional_Level/REGIONAL_",region_names,"_FixedYAxis_",reg,".tif"), width = 4, height = 3, units = "in", res = 300)
    print(reg_plot2)
    dev.off()
    
  }
  
  names(regional_change_samples) <- names(regional_trend_samples) <- names(region_plots) <- unique(N_region$region)
  
  return(list(regional_abundance_summary = regional_abundance_summary,
              regional_trend_summary = regional_trend_summary,
              regional_change_summary = regional_change_summary,
              regional_trend_samples = regional_trend_samples,
              regional_change_samples = regional_change_samples,
              region_plots = region_plots))
}

# Regional trends based on fast ice regions
fast_ice_reg <- regional_estimate_fn(region_names = "ice_reg", 
                                     region_colors = c("Amundsen Sea" = "#686868",
                                                       "Australia" = "#dcdcdc",
                                                       "Bellingshausen Sea" = "#feaa02",
                                                       "Dronning Maud Land" = "#00a884",
                                                       "East Indian Ocean" = "#1558bf",
                                                       "Victoria Oates Land" = "#ff73de",
                                                       "Weddell Sea" = "#01c5ff",
                                                       "West Indian Ocean" = "#ffff00"), 
                                     N_samples = N_samples)

# Regional trends based on pack ice regions
pack_ice_reg <- regional_estimate_fn(region_names = "p_ice_reg",
                                     region_colors = c("Ross" = "#5a4dd8",
                                                       "Bell-Amundsen" = "#57abce",
                                                       "Weddell" = "#54edac",
                                                       "Indian" = "#ffff00",
                                                       "Pacific" = "#f6636d"), 
                                     N_samples = N_samples)

#----------------------------------------------------------
# Correlation between regional FAST ice trends and population trends
#----------------------------------------------------------

icetrend <- read.csv("../data/fast_ice_trends.csv")
popchange_fastice = fast_ice_reg$regional_change_summary %>% full_join(icetrend)

nameColor <- bquote(atop(Minimum~fast,
                         ice~extent~(km^2)))

sea_ice_plot <- ggplot(data = popchange_fastice,
                       aes(x = FastIceTrend,
                           y = Prob_Decline,
                           col = FastIceExtent_min*1000,
                           label = Region)) +
  geom_point(size = 4)+
  geom_text_repel(data = popchange_fastice,
                  aes(x = FastIceTrend,
                      y = Prob_Decline,
                      col = FastIceExtent_min*1000,
                      label = Region),
                  size = 3.5,
                  min.segment.length = 5,
                  hjust = 1,
                  nudge_x = -0.3,
                  alpha = 0.5)+
  #scale_color_gradientn(colors = c("gray90","blue"))+
  scale_color_gradientn(colors = rev(inferno(50)[12:40]))+
  xlab("Fast ice trend\n(% change per year)")+
  ylab("Probability of population decline")+
  coord_cartesian(xlim = c(-4,4))+
  theme_few()+
  labs(color = nameColor)+
  theme(legend.position = "right")
#print(sea_ice_plot)

# Save figure
tiff(filename = "output/model_results/2_Regional_Level/sea_ice_correlation.tif", width = 7, height = 3.5, units = "in", res = 300)
print(sea_ice_plot)
dev.off()

# Correlation between fast ice trend and regional population trends
rho = DescTools::SpearmanRho(popchange_fastice$Prob_Decline,popchange_fastice$FastIceTrend,conf.level = 0.95)
rho

# ***************************************************************************
# ***************************************************************************
# DEPRECATED CODE THAT COULD BE USEFUL FOR VISUALIZING/INTERPRETING DYNAMICS
# ***************************************************************************
# ***************************************************************************

# #----------------------------------------------------------
# # Plot magnitude of change at each colony on a map
# #----------------------------------------------------------
# 
# # Change categories (>100% decrease,50-100% decrease, 0-50% decrease, 0-50% increase, 50-100% increase)
# df_2009 = subset(colony_summary, year == 2009)
# df_2018 = subset(colony_summary, year == 2018)
# world <- map_data("world")
# lim <- max(abs(df_2018$change_since_2009_mean),na.rm = TRUE)
#
# trend_map <- ggplot(world, aes(x=long, y=lat, group=group)) +
#   geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
#   scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
#   scale_x_continuous(breaks=(-4:4) * 45) +
#   coord_map("ortho", orientation=c(-90, 0, 0)) +
#   theme_bw()+
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks=element_blank(),
#         panel.border = element_blank()) +
#   geom_point(data = df_2018,
#              aes(x=lon, y=lat,group=1,
#                  col = change_since_2009_mean#,size = change_since_2009_mean
#                  ))+
#   geom_label_repel(data = df_2018,aes(x=lon, y=lat,group=1,
#                                       label = site_id,
#                                       col = change_since_2009_mean
#                                       ))+
#   scale_color_gradientn(colors = c("red","gray90","blue"),limits = c(-lim,lim), name = "Change since 2009")+
#   scale_size_continuous(name = "Change since 2009")
# print(trend_map)
#
# pdf(file = "./output_empirical/FigX_trend_map.pdf", width = 8, height = 8)
# print(trend_map)
# dev.off()
# 
# # -------------------------------
# # Calculate sd of year-to-year change in log(global abundance)
# # ((Process variance of global abundance))
# # -------------------------------
# 
# global_sd = log(out$sims.list$N_global) %>% 
#   apply(.,1,diff) %>%
#   apply(.,2,sd)
# global_sd
# 
# hist(global_sd)
# 
# #----------------------------------------------------------
# # Evaluate magnitude of change between years at each colony
# #----------------------------------------------------------
# 
# change_annual <- data.frame()
# 
# # For each colony, in each year, for each mcmc sample, calculate percent change between years
# for (mc_sample in 1:dim(out$sims.list$N)[1]){
#   
#   N_matrix <- out$sims.list$N[mc_sample,,]
#   N_matrix[N_matrix == 0] <- NA
#   
#   tmp <- apply(log(N_matrix),1,function(x)diff(x)) %>%
#     reshape2::melt() %>%
#     dplyr::rename(year_number = Var1,
#                   site_number = Var2,
#                   log_change = value) %>%
#     mutate(mc_sample = mc_sample)
#   
#   
#   change_annual <- rbind(change_annual, tmp)
#   
#   print(mc_sample)
# }
# 
# # Calculate summaries of magnitude of change between consecutive years at each colony
# # test1 <- change_annual %>%
# #   group_by(mc_sample, site_number) %>%
# #   summarize(ratio_med = median(ratio, na.rm = TRUE),
# #             percent_change_med = median(percent_change, na.rm = TRUE))
# 
# test2 <- change_annual %>%
#   group_by(year_number, site_number) %>%
#   summarize(log_change_mean = mean(log_change, na.rm = TRUE),
#             log_change_q05 = quantile(log_change,0.05, na.rm = TRUE),
#             log_change_q95 = quantile(log_change,0.95, na.rm = TRUE)) %>%
#   full_join(colony_attributes[,c("site_id","site_number")],.)
# 
# scale_y <- data.frame(log_diff = c(log(0.5),log(0.75),log(1),-log(0.75),-log(0.5)))
# scale_y$percent_change <- (100*(exp(scale_y$log_diff)-1)) %>% round()
# scale_y$label <- paste0(scale_y$percent_change,"%")
# scale_y$label[4:5] <- paste0("+",scale_y$label[4:5])
# scale_y
# 
# ggplot(test2, aes(x = year_number, y = log_change_mean, ymin = log_change_q05, ymax = log_change_q95))+
#   geom_point()+
#   geom_errorbar(width=0)+
#   facet_wrap(site_id~.)+
#   scale_y_continuous(breaks = scale_y$log_diff, labels = scale_y$label)+
#   theme_bw()+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank())+
#   geom_hline(yintercept = c(log(0.5),-log(0.5)), col = "orangered", alpha = 0.5)+
#   geom_hline(yintercept = c(log(0.75),-log(0.75)), col = "orangered", alpha = 0.25)+
#   geom_hline(yintercept = 0, col = "orangered", alpha = 0.1)
# 
# #----------------------------------------------------------
# # Sequentially remove each colony and recalculate global trend
# #  - provides insight into the effect each colony has on the global trend
# #----------------------------------------------------------
# 
# colony_removal_results <- data.frame(site_removed = "All sites included",
#                                      
#                                      percent_change_q50 = global_change_trend$percent_change_summary["50%"],
#                                      percent_change_q05 = global_change_trend$percent_change_summary["5%"],
#                                      percent_change_q95 = global_change_trend$percent_change_summary["95%"],
#                                      
#                                      OLS_regression_q50 = global_change_trend$OLS_regression_summary["50%"],
#                                      OLS_regression_q05 = global_change_trend$OLS_regression_summary["5%"],
#                                      OLS_regression_q95 = global_change_trend$OLS_regression_summary["95%"])
# 
# 
# for (i in 1:nrow(colony_attributes)){
#   
#   # Remove focal colony from global sum
#   N_global_minus1 <- out$sims.list$N[,-i,] %>%
#     apply(.,c(1,3),sum)
#   
#   colnames(N_global_minus1) <- year_range
#   N_global_minus1 <- N_global_minus1 %>% 
#     reshape2::melt() %>% 
#     rename(mcmc_sample = Var1, Year = Var2, N = value)
#   
#   N_global_minus1_matrix = N_global_minus1 %>% 
#     spread(Year, N) %>%
#     dplyr::select(-mcmc_sample) %>%
#     as.matrix()
#   
#   global_change_trend_minus1 <- change_trend_fn(N_global_minus1_matrix)
#   
#   colony_removal_results <- rbind(colony_removal_results, 
#                                   data.frame(site_removed = colony_attributes$site_id[i],
#                                              
#                                              percent_change_q50 = global_change_trend_minus1$percent_change_summary["50%"],
#                                              percent_change_q05 = global_change_trend_minus1$percent_change_summary["5%"],
#                                              percent_change_q95 = global_change_trend_minus1$percent_change_summary["95%"],
#                                              
#                                              OLS_regression_q50 = global_change_trend_minus1$OLS_regression_summary["50%"],
#                                              OLS_regression_q05 = global_change_trend_minus1$OLS_regression_summary["5%"],
#                                              OLS_regression_q95 = global_change_trend_minus1$OLS_regression_summary["95%"]
#                                   ))
# }
# 
# # Merge with colony attributes table
# colony_removal_results <- colony_removal_results %>% full_join(colony_attributes, by = c("site_removed" = "site_id"))
# 
# colony_removal_results$p_ice_reg[1] <- ""
# colony_removal_results$ice_reg[1] <- ""
# 
# colony_removal_results$p_ice_reg <- factor(colony_removal_results$p_ice_reg,
#                                            levels = c(
#                                              "Ross",
#                                              "Bell-Amundsen",
#                                              "Weddell",
#                                              "Indian",
#                                              "Pacific",
#                                              ""))
# 
# colony_removal_results <- colony_removal_results %>%
#   arrange(p_ice_reg,img_lon_mean)
# 
# 
# colony_removal_results$site_removed <- factor(colony_removal_results$site_removed,
#                                               levels = colony_removal_results$site_removed)
# 
# 
# ggplot(colony_removal_results, aes(y = site_removed, 
#                                    x = OLS_regression_q50,
#                                    xmin = OLS_regression_q05,
#                                    xmax = OLS_regression_q95,
#                                    col = p_ice_reg))+
#   geom_vline(xintercept = colony_removal_results$OLS_regression_q50[which(colony_removal_results$site_removed == "All sites included")],
#              col = "gray80", size = 3)+
#   geom_vline(xintercept = colony_removal_results$OLS_regression_q05[which(colony_removal_results$site_removed == "All sites included")],
#              col = "gray80", size = 1.5)+
#   geom_vline(xintercept = colony_removal_results$OLS_regression_q95[which(colony_removal_results$site_removed == "All sites included")],
#              col = "gray80", size = 1.5)+
#   
#   geom_point(size=4)+
#   geom_errorbarh(height=0)+
#   theme_bw()+
#   #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   scale_color_manual(values=c(
#     "#5a4ed8",
#     "#57abcf",
#     "#54edac",
#     "#feff00",
#     "#f6626d",
#     
#     "black"),
#     guide = "none")+
#   xlab("Trend Estimate")+
#   ylab("Site Omitted")+
#   facet_grid(rows = vars(p_ice_reg), scales = "free_y", space = "free") +
#   theme(panel.spacing.x = unit(10, "cm"),
#         strip.text.x = element_blank())+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   ggtitle("Global trend estimate after omitting each colony\n")
# 
# #----------------------------------------------------------
# # Estimate change over 3 generations and IUCN thresholds
# # Generation time = 16 years (Jenouvrier et al. 2014 Nature Climate Change)
# # 3 generations = 48 years
# #----------------------------------------------------------
# 
# annual_trend = global_change_trend$OLS_regression_samples
# 
# # Estimated percent of population remaining after 3 generations
# IUCN_3generation = 100*exp(annual_trend*48)
# 
# # Probability the population will decline by more than 80% (Critically Endangered)
# mean(IUCN_3generation <= 20)
# 
# # Probability the population will decline by more than 50% (Endangered)
# mean(IUCN_3generation <= 50)
# 
# # Probability the population will decline by more than 30% (Vulnerable)
# mean(IUCN_3generation <= 70)