library(readr)


beetles_target <- read_csv("./target_drivers/beetles-targets.csv.gz")

beetles_null_forecast <- read_csv("./forecast/beetles-2020-EFI_avg_null.csv.gz")


beetles_null_forecast <- beetles_null_forecast%>%
  group_by(siteID, time)%>%
  summarize(mean_rich = mean(richness),
            mean_abun = mean(abundance))
