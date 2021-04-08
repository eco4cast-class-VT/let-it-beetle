library(readr)
library(tidyverse)



beetles_null_forecast <- read_csv("./forecast/beetles-2020-EFI_avg_null.csv.gz")
beetles_null_forecast_test <- beetles_null_forecast%>%
  group_by(time, siteID)%>%
  summarize(mean_rich = mean(richness),
            upper_rich = quantile(richness, 0.95, na.rm = T),
            lower_rich = quantile(richness, 0.05, na.rm = T),
            mean_abun = mean(abundance),
            upper_abun = quantile(abundance, 0.95, na.rm = T),
            lower_abun = quantile(abundance, 0.05, na.rm = T))%>%
  arrange(siteID, time)

beetles_target <- read_csv("./target_drivers/beetles-targets.csv.gz")


richness_null_forecasts <- beetles_null_forecast_test %>%
  ggplot(., aes(x = time, y = mean_rich)) +
  geom_ribbon(aes(ymin = lower_rich, ymax = upper_rich), alpha = 0.5, fill = "midnightblue") +
  geom_line(color = "black")+
  geom_point(data = beetles_target, aes(x = time, y = richness), inherit.aes = FALSE, pch = 21, color = "black", fill = "red", cex = 2) +
  theme_bw()+
  labs(title = "Beetle null")+
  ylab("Richness")+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2020-01-01"),as.Date("2021-01-01")))+
  theme(axis.text=element_text(size=8, color = "black"),
        axis.title=element_text(size=8, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 8),legend.position = "none",
        legend.text = element_text(size = 8, color = "black"))+
  facet_wrap(~siteID)

jpeg("./figures/BEETLE_NULL.jpg", width = 1600, height = 1000)
richness_null_forecasts
dev.off()
