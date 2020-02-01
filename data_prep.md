Import and prepare data for analyses
================
Rebecca Batstone
2020-02-01

## load packages

``` r
library("tidyverse") ## includes ggplot2, dplyr, readr, stringr
library("reshape2") ## for melting datasheets
library("car") ## Anova function
library("cowplot") # paneled graphs
library("lme4") ## mixed effects models
library("knitr") ## knit analyses files
# library("multcomp") ## post-hoc comps
# library("multcompView") ## summarize key variables
library("fitdistrplus") ## determine best prob distributions
# library("logspline") ## visualize prob dists
library("DHARMa") ## residual diagnostics for GLMMs
library("emmeans") ## post-hoc contrasts, lsmeans
library("lattice") ## graphs

#model convergence packages
library("numDeriv")
# library("RCurl") ## to source() from Github
library("RColorBrewer") ## graph colors
```

# Plant-level data:

  - every row is an individual plant, unless melted
  - when melted, always include plant as a random factor (nested within
    tray and batch)

## load spreadsheet (plant\_info.csv)

  - compiled by Xue
  - columns 15 onward from nod-level data (Judith)
      - summarized using pivot table for plant\_ID
      - total count of plates with 1022 and 1021 on each root half
      - didn’t include any unsuccessfully scored or mixed plates
        (1021/1022)
      - 128 plants (out of 338) in which both roots have a propEm1022
        calculated

<!-- end list -->

``` r
# specify directory, column format
plants <- read_csv("./raw_data/plant_info1.csv", 
    col_types = cols(line = col_factor(levels = c("270", 
        "276", "267", "313", "279")), 
        treat = col_factor(levels = c("control", 
            "one", "two", "three"))))
plants$plant_ID <- as.factor(plants$plant_ID)

# rename treatments
levels(plants$treat)[levels(plants$treat)=="control"] <- "50:50"
levels(plants$treat)[levels(plants$treat)=="one"] <- "20:80"
levels(plants$treat)[levels(plants$treat)=="two"] <- "10:90"
levels(plants$treat)[levels(plants$treat)=="three"] <- "2:98"

# order treatment
plants$treat <- factor(plants$treat, levels=c("50:50", "20:80", "10:90", "2:98"),
                       ordered=T)

# calculate Em1022 proportions for each root half
plants$high_p1022 <- (plants$high_1022_plates)/(plants$high_scored)
plants$low_p1022 <- (plants$low_1022_plates)/(plants$low_scored)
# write.csv(plants, file="plants_prop.csv")

# calculate number of Em1022 and Em1021 on each side (use p1022 data)
plants$high_1022 <- (plants$high_p1022)*(plants$high_nod)
plants$low_1022 <- (plants$low_p1022)*(plants$low_nod)
plants$high_1021 <- (1-(plants$high_p1022))*(plants$high_nod)
plants$low_1021 <- (1-(plants$low_p1022))*(plants$low_nod)

# transform the data
plants$shoot.t <- plants$shoot*1000

# str(plants)
# 338 obs., 338 plants

# save data
save(plants, file = "./prepped_data/plants.Rdata")
```

## Half-root models

Prep spreadsheet:

``` r
# compare halves first
# creat long version of plants, melted for half (root)
plants_long_root <- melt(plants, measure.vars = c("low_root", "high_root"),
                         variable.name="half", value.name="root")
#rename levels
plants_long_root <- plants_long_root %>% mutate(half=recode_factor(half, 
                         low_root="Low-N",
                         high_root="High-N"))
# str(plants_long_root)
## 676 obs., 338 plants

# create long version of plants, melted for half (nod)
plants_long_nod <- melt(plants, measure.vars = c("low_nod", "high_nod"),
                         variable.name="half", value.name="nod")
# rename levels
plants_long_nod <- plants_long_nod %>% mutate(half=recode_factor(half, 
                         low_nod="Low-N",
                         high_nod="High-N"))
# str(plants_long_nod)
## 676 obs., 338 plants

# create long version of plants, melted for half (prop)
plants_long_prop <- melt(plants, measure.vars = c("low_p1022", "high_p1022"),
                         variable.name="half", value.name="prop")
#rename levels
plants_long_prop <- plants_long_prop %>% mutate(half=recode_factor(half, 
                         low_prop="Low-N",
                         high_prop="High-N"))
# str(plants_long_prop)
## 676 obs., 338 plants

# combine spreadsheets
plants_long1 <- cbind(plants_long_root, plants_long_nod$nod, plants_long_prop$prop)
plants_long <- plants_long1 %>%
  rename(nod = 'plants_long_nod$nod', prop = 'plants_long_prop$prop')
# transform root
plants_long$root.t <- plants_long$root*1000

# order half
plants_long$half <- factor(plants_long$half, levels=c("Low-N","High-N"), ordered=T)

# str(plants_long)
## 676 obs., 338 plants (2 obs per plant)

# save data
save(plants_long, file = "./prepped_data/plants_long.Rdata")
```

## Half-root models (accounting for strain)

Note: nodule numbers of each strain are based on proportions. In some
cases, even though plants formed nodules, a proportion wasn’t calculated
b/c no nodules were successfully plated. So, there are fewer plants w/
nod occupancy data (n = 128) compared to plants with nodule number
counts (n=338)

``` r
# use full dataset instead
# break down according to half

# creat long version of plants, melted for half (root)
plants_root <- melt(plants, measure.vars = c("low_root", "high_root"),
                         variable.name="half", value.name="root")
#rename levels
plants_root <- plants_root %>% mutate(half=recode_factor(half, 
                         low_root="Low-N",
                         high_root="High-N"))

# create long version of plants, melted for half (1021 nodules)
plants_1021 <- melt(plants, measure.vars = c("low_1021", "high_1021"),
                         variable.name="half", value.name="Em1021_nod")
#rename levels
plants_1021 <- plants_1021 %>% mutate(half=recode_factor(half, 
                         low_1021="Low-N",
                         high_1021="High-N"))

# create long version of plants, melted for half (1022 nodules)
plants_1022 <- melt(plants, measure.vars = c("low_1022", "high_1022"),
                         variable.name="half", value.name="Em1022_nod")
#rename levels
plants_1022 <- plants_1022 %>% mutate(half=recode_factor(half, 
                         low_1022="Low-N",
                         high_1022="High-N"))

#combine spreadsheets
plants_long2 <- cbind(plants_root, plants_1021$Em1021_nod,  plants_1022$Em1022_nod)
plants_long3 <- plants_long2 %>%
  rename(Em1021_nod = 'plants_1021$Em1021_nod', Em1022_nod = 'plants_1022$Em1022_nod')
plants_long3$Em1022_nod.r <- round(plants_long3$Em1022_nod)
plants_long3$Em1021_nod.r <- round(plants_long3$Em1021_nod)

# order half
plants_long3$half <- factor(plants_long3$half, levels=c("Low-N","High-N"), ordered=T)
# str(plants_long3)
## 676 obs., 338 plants (2 obs. per plant)

# Compare root responses between strains

# create long version of plants_long3, melted for strain
plants_long_strain <- melt(plants_long3, measure.vars = c("Em1021_nod.r", "Em1022_nod.r"),
                         variable.name="strain", value.name="nod")
# rename levels
plants_long_strain <- plants_long_strain %>% mutate(strain=recode_factor(strain, 
                         Em1021_nod.r="Em1021",
                         Em1022_nod.r="Em1022"))
plants_long_strain$nod <- as.integer(plants_long_strain$nod)
# str(plants_long_strain)
## 1352 obs., 338 plants (4 obs. per plant)

# save data
save(plants_long_strain, file = "./prepped_data/plants_long_strain.Rdata")
```
