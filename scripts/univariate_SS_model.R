library(rjags)
library(tidybayes)
library(readr)
library(tidyverse)

beetles_target <- read_csv("./target_drivers/beetles-targets.csv.gz")

ss_univariate = ("ss_model.txt")
jagsscript = cat("
model {  
   # priors on parameters
   mu ~ dnorm(0, 0.01); 
   tau.pro ~ dgamma(0.001,0.001); 
   sd.q <- 1/sqrt(tau.pro);
   tau.obs ~ dgamma(0.001,0.001);
   sd.r <- 1/sqrt(tau.obs); 
   phi ~ dnorm(0,1);
   
   X[1] <- mu;
   predY[1] <- X[1];
   Y[1] ~ dnorm(X[1], tau.obs);

   for(i in 2:N) {
      predX[i] <- phi*X[i-1]; 
      X[i] ~ dnorm(predX[i],tau.pro); # Process variation
      predY[i] <- X[i];
      Y[i] ~ dnorm(X[i], tau.obs); # Observation variation
   }
}  
",file = ss_univariate)

# GET SITE NAMES
#site_names <- c(unique(beetles_target$siteID))

# FORECAST HORIZON (weekly for 20 months == 80 weeks (right?))
  
  # Select site
  site_data_var <- beetles_target %>%
    #filter(siteID == site_names[s])
    filter(siteID == "ORNL")
  
  observed_weeks <- as.data.frame(seq(as.Date("2013-07-01"), Sys.Date(), by = "weeks"))%>%
    rename(time = `seq(as.Date("2013-07-01"), Sys.Date(), by = "weeks")`)
    
  observed_weeks_44cast <- left_join(observed_weeks, site_data_var, by = "time")%>%
    #mutate(siteID = site_names[s])
    mutate(siteID == "ORNL")
  
  forecast_weeks <- as.data.frame(seq(max(observed_weeks_44cast$time), Sys.Date()+560, by = "weeks"))%>%
    rename(time = `seq(max(observed_weeks_44cast$time), Sys.Date() + 560, by = \"weeks\")`)%>%
    mutate(siteID = "ORNL",
           abundance = NA,
           richness = NA)
  
  forecast_timeline <- bind_rows(observed_weeks_44cast,forecast_weeks[-1,])
  
  jags.data = list(Y = forecast_timeline$richness, N = nrow(forecast_timeline))

  jags.params = c("predY", "mu")

  nchain = 5
  chain_seeds <- c(200,800,1400)
  #Initialize JAGS model
  j.model   <- jags.model(file = ss_univariate,
                          data = jags.data,
                          n.chains = 5)
  
  #Run JAGS model and sample from the posteriors
  jags.out   <- coda.samples(model = j.model,
                             variable.names = jags.params, 
                             n.iter = 50000, n.burnin = 10000)
  
  richness_forecast <- jags.out %>%
    spread_draws(predY[week]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = forecast_timeline$time[week]) %>%
    ungroup() %>%
    select(time, predY, ensemble)%>%
    group_by(time) %>% 
    summarise(mean = mean(predY),
              upper_95 = quantile(predY, 0.90),
              lower_95 = quantile(predY, 0.10),
              var = var(predY),
              sd = sd(predY),.groups = "drop")%>%
    mutate(lower_95 = ifelse(lower_95<=0,0.01,lower_95))
  
  current_forecast <- jags.out %>%
    spread_draws(predY[week]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = forecast_timeline$time[week]) %>%
    ungroup() %>%
    select(time, predY, ensemble)%>%
    group_by(time) %>% 
    summarise(mean = mean(predY),
              upper_95 = quantile(predY, 0.90),
              lower_95 = quantile(predY, 0.10),
              var = var(predY),
              sd = sd(predY),.groups = "drop")%>%
    mutate(lower_95 = ifelse(lower_95<=0,0.01,lower_95))%>%
    filter(time >= Sys.Date())
  
  ggplot(richness_forecast, aes(x = time, y = mean))+
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2, fill = "midnightblue") +
    geom_line(color = "black")+
    geom_point(data = forecast_timeline, aes(x = time, y = richness), color = "red") +
    labs(title = "ORNL")+
    theme_bw()+
    geom_vline(xintercept = Sys.Date())+
    labs(x = "Date", y = "Beetle Richness")
  
  ggsave(paste0("beetle_richness_ORNL_",Sys.Date(),"_figure.pdf"), device = "pdf")
  

