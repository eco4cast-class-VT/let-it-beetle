
if (!require('pacman')) install.packages('pacman'); library('pacman')
pacman::p_load(tidyverse, ISOweek, neonstore)

source("./scripts/resolve_taxonomy.R")

# neonstore::neon_dir()
Sys.setenv("NEONSTORE_HOME" = "/Users/ryanmcclure/Documents/let-it-beetle/neonstore")
Sys.setenv("NEONSTORE_DB" = "/Users/ryanmcclure/Documents/let-it-beetle/neonstore")
# 
# neonstore::neon_download(product="DP1.10022.001")
# neonstore::neon_store(product="DP1.10022.001")


## Load data from raw files
sorting <- neon_table("bet_sorting")
para <- neon_table("bet_parataxonomistID")
expert <- neon_table("bet_expertTaxonomistIDProcessed")
field <- neon_table("bet_fielddata")


#### Generate derived richness table  ####################
beetles <- resolve_taxonomy(sorting, para, expert) %>% 
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1")))

richness <- beetles %>%  
  select(taxonID, siteID, collectDate, time) %>%
  distinct() %>%
  count(siteID, time) %>% 
  rename(richness = n)  %>%
  ungroup()



#### Generate derived abundance table ####################

effort <- field %>% 
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>% 
  group_by(siteID, time) %>% 
  summarise(trapnights = as.integer(sum(collectDate - setDate)),
            .groups = "drop")

counts <- sorting %>% 
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>%
  group_by(siteID, time) %>%
  summarise(count = sum(as.numeric(individualCount), na.rm = TRUE), 
            .groups = "drop")

abund <- counts %>% 
  left_join(effort) %>% 
  arrange(time) %>%
  mutate(abundance = count / trapnights) %>% 
  select(siteID, time, abundance) %>%
  ungroup()

targets <- full_join(abund, richness)



##  Write out the targets
write_csv(targets, "beetles-targets.csv.gz")