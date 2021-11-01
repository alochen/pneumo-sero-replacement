############ Serotype Replacement Paper #################

library(tidyr)
library(dplyr)
library(ggplot2)
library(fmsb)
library(gridExtra)
library(purrr)
library(reshape)
library(lemon)
library(ggpubr)
library(ggstance)

##########################################################################################################################################################
########### DATA PREP ####################################################################################################################################
cols <- c("VT7" = "#35B779FF", "VT10"= "#31688EFF", "VT13" = "#E69F00", "NVT" = "#440154FF")

sero15bcagg <- function(df) { # fn that aggregates serotypes 15B and 15C
  sero15s <- df[startsWith(as.character(df$Serotype), prefix = "15"),]$Serotype
  sero15bc <- sero15s[!sero15s == "15A"]
  if (length(sero15bc) > 0) 
  {sero15bc.df <- df[df$Serotype %in% sero15bc,]
  sero15bc.aggreg <- aggregate(.~ Age.group, sero15bc.df, sum)
  sero15bc.aggreg$Serotype <- "15B/C"
  df <- df[!df$Serotype %in% sero15bc,]
  df <- rbind(df, sero15bc.aggreg) }
}

totdf <- function(df, country) { # fn that melts data and adds population of age group and age category
  colnames(df) <- gsub("^X", "", colnames(df))
  new.df <- melt(df)
  colnames(new.df) <- c("Serotype", "Age.group", "Year", "dis.cases")
  new.df$country <- country
  new.df$Year <- as.numeric(levels(new.df$Year))[new.df$Year]
  
  PCVyrs <- PCVintro[which(rownames(PCVintro) == country),]
  PCVyrs <- PCVyrs[!is.na(PCVyrs)]
  PCVera.spec <- PCVera[,which(colnames(PCVera) == country)]
  PCVera.spec <- PCVera.spec[!is.na(PCVera.spec)]
  
  prepcv_ind <- which(new.df$Year <= PCVyrs[1])
  prepcv_high <- which(new.df$Year > PCVyrs[1] & new.df$Year <= PCVyrs[2])
  postpcv_ind <- which(new.df$Year > PCVyrs[2])
  
  if(length(PCVyrs) < 2) {postpcv_ind <- which(new.df$Year > PCVyrs[1])}
  
  new.df$pcv_era <- NA
  new.df$pcv_era[prepcv_ind] <- as.character(PCVera.spec[1])
  new.df$pcv_era[prepcv_high] <- as.character(PCVera.spec[2])
  new.df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[3])
  
  if(length(PCVyrs) < 2) {new.df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[2])}
  
  new.df$agecat[new.df$Age.group %in% children] <- "children"
  new.df$agecat[new.df$Age.group %in% adults] <- "adults"
  
  new.df$population <- unlist(lapply(1:nrow(new.df), function(x) 
    pop.allages[which(pop.allages$Population == country & pop.allages$Age.group == new.df$agecat[x]), 
          which(colnames(pop.allages) == paste("X", new.df$Year[x], sep = ""))]))
  
  new.df$agecat[new.df$agecat == "children"] <- "Children"
  new.df$agecat[new.df$agecat == "adults"] <- "Adults"
  
  return(new.df)
}

freqdf <- function(df, country) { # fn that returns frequency of each serotype each year
  matrix <- prop.table(as.matrix(df[3:ncol(df)]), 2)
  new.df <- as.data.frame.matrix(matrix)
  freq.df <- data.frame(Serotype = df$Serotype, Age.group = df$Age.group, new.df)
  colnames(freq.df) <- gsub("^X", "", colnames(freq.df))
  freq.df <- melt(freq.df)
  colnames(freq.df) <- c("Serotype", "Age.group", "Year", "freq")
  freq.df$Year <- as.character(freq.df$Year)
  freq.df$country <- country
  
  PCVyrs <- PCVintro[which(rownames(PCVintro) == country),]
  PCVyrs <- PCVyrs[!is.na(PCVyrs)]
  PCVera.spec <- PCVera[,which(colnames(PCVera) == country)]
  PCVera.spec <- PCVera.spec[!is.na(PCVera.spec)]
  
  prepcv_ind <- which(freq.df$Year <= PCVyrs[1])
  prepcv_high <- which(freq.df$Year > PCVyrs[1] & freq.df$Year <= PCVyrs[2])
  postpcv_ind <- which(freq.df$Year > PCVyrs[2])
  
  if(length(PCVyrs) < 2) {postpcv_ind <- which(freq.df$Year > PCVyrs[1])}
  
  freq.df$pcv_era <- NA
  freq.df$pcv_era[prepcv_ind] <- as.character(PCVera.spec[1])
  freq.df$pcv_era[prepcv_high] <- as.character(PCVera.spec[2])
  freq.df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[3])
  
  if(length(PCVyrs) < 2) {freq.df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[2])}
  
  freq.df$Year <- as.numeric(freq.df$Year)
  
  freq.df$agecat[freq.df$Age.group %in% children] <- "Children"
  freq.df$agecat[freq.df$Age.group %in% adults] <- "Adults"

  return(freq.df)
}

PCVintro <- data.frame(PCV7intro = c(Australia = 2005, Finland = NA, #2009, #include Finland PCV7 date for Simpson's Diversity Index calculation 
                                     France = 2006, Italy = 2005, Norway = NA, USA = 2000), # 2006 for Norway, include Norway PCV7 date for SDI calculation
                       PCV10intro = c(NA, 2010, NA, NA, NA, NA), 
                       PCV13intro = c(2011, NA, 2010, 2010, 2011, 2010))

PCVera <- data.frame(Australia = c("Pre-PCV", "Pre-PCV13", "Post-PCV"), 
                     Canada = c("Pre-PCV", "Pre-PCV13", "Post-PCV"),
                     Finland = c(NA, "Pre-PCV10", "Post-PCV"),
                     France = c("Pre-PCV", "Pre-PCV13", "Post-PCV"),
                     Italy = c("Pre-PCV", "Pre-PCV13", "Post-PCV"),
                     Norway = c(NA, "Pre-PCV13", "Post-PCV"),# "Pre-PCV"
                     USA = c("Pre-PCV", "Pre-PCV13", "Post-PCV"))

#pop <- read.csv("Q:/Technical/R/SeroReplacement/Populationvsyr.csv")
pop <- read.csv("Q:/Technical/R/SeroReplacement/IPDincidence.csv")
pop.allages <- read.csv("Q:/Technical/R/SeroReplacement/Populationvsyr_allages.csv")

VT7 <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
VT10 <- c("1", "5", "7F")
VT13 <- c("3", "6A", "19A")

australia_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Australia age strat.csv")
australia_agestrat <- australia_agestrat[-which(australia_agestrat$Age.group == "Unknown"),] # remove Unknown age category
australia_agestrat$Age.group <- factor(australia_agestrat$Age.group, levels = c("0-1 yrs", "2-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs"))
australia_agestrat <- sero15bcagg(australia_agestrat)
# remove first two years of australia bc of small sample size
australia_agestrat <- australia_agestrat[,-c(3,4)]
#
australia_eachyr <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Australia age strat.csv")
australia_eachyr <- aggregate(. ~ Serotype, australia_eachyr, FUN = sum)
australia_eachyr$Age.group <- NA
australia_eachyr <- australia_eachyr[,-c(3,4)]
oz_lessthan18 <- australia_agestrat %>% filter(Age.group %in% c("0-1 yrs", "2-4 yrs", "5-17 yrs"))
oz_lessthan18 <- aggregate(. ~ Serotype, oz_lessthan18, FUN = sum)
oz_lessthan18$Age.group <- "< 18 yrs"
oz_morethaneq18 <- australia_agestrat %>% filter(Age.group %in% c("18-49 yrs", "50-64 yrs", "65+ yrs"))
oz_morethaneq18 <- aggregate(. ~ Serotype, oz_morethaneq18, FUN = sum)
oz_morethaneq18$Age.group <- "> 18 yrs"
oz_18 <- bind_rows(oz_lessthan18, oz_morethaneq18)

finland_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Finland age strat.csv")
finland_agestrat$Age.group <- factor(finland_agestrat$Age.group, levels = c("0 - 1 yrs", "2 - 4 yrs", "5 - 17 yrs", "18 - 49 yrs", "50 - 64 yrs", "65+ yrs"))
finland_agestrat <- sero15bcagg(finland_agestrat)
finland_eachyr <- aggregate(. ~ Serotype, finland_agestrat, FUN = sum)
finland_eachyr$Age.group <- NA
finland_lessthan18 <- finland_agestrat %>% filter(Age.group %in% c("0 - 1 yrs", "2 - 4 yrs", "5 - 17 yrs"))
finland_lessthan18 <- aggregate(. ~ Serotype, finland_lessthan18, FUN = sum)
finland_lessthan18$Age.group <- "< 18 yrs"
finland_morethaneq18 <- finland_agestrat %>% filter(Age.group %in% c("18 - 49 yrs", "50 - 64 yrs", "65+ yrs"))
finland_morethaneq18 <- aggregate(. ~ Serotype, finland_morethaneq18, FUN = sum)
finland_morethaneq18$Age.group <- "> 18 yrs"
finland_18 <- bind_rows(finland_lessthan18, finland_morethaneq18)

france_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/France age strat.csv")
france_agestrat[is.na(france_agestrat)] <- 0 # change empty cells to 0
france_agestrat$Age.group <- factor(france_agestrat$Age.group, levels = c("<16 yrs", ">16 yrs"))
france_agestrat <- sero15bcagg(france_agestrat)
france_eachyr <- aggregate(. ~ Serotype, france_agestrat, FUN = sum)
france_eachyr$Age.group <- NA
france_lessthan16 <- france_agestrat %>% filter(Age.group %in% "<16 yrs")
france_morethan16 <- france_agestrat %>% filter(Age.group %in% ">16 yrs")
france_16 <- france_agestrat
france_16$Age.group <- as.character(france_16$Age.group)

# italy_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Italy age strat.csv")
# italy_agestrat[is.na(italy_agestrat)] <- 0 # change empty cells to 0
# italy_agestrat$Age.group <- factor(italy_agestrat$Age.group, levels = c("0-4 yrs", "64+ yrs"))
# italy_agestrat <- sero15bcagg(italy_agestrat)
# italy_eachyr <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Italy.csv")
# italy_eachyr[is.na(italy_eachyr)] <- 0
# italy_eachyr$Age.group <- NA
# italy_eachyr <- italy_eachyr[,c(1,ncol(italy_eachyr), 2:(ncol(italy_eachyr)-1))]
# italy_inf <- italy_agestrat %>% filter(Age.group %in% "0-4 yrs")
# italy_eld <- italy_agestrat %>% filter(Age.group %in% "64+ yrs")
# italy_eld.colsums <- colSums(italy_eld[,3:ncol(italy_eld)])
# italy_eld <- italy_eld[,-(which(italy_eld.colsums == 0)+2)]
# italy_inf$Age.group <- as.character(italy_inf$Age.group)

norway_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/Norway age strat.csv")
norway_agestrat[is.na(norway_agestrat)] <- 0
norway_agestrat$Age.group <- factor(norway_agestrat$Age.group, levels = c("0-1 yrs", "2-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs"))
norway_agestrat <- sero15bcagg(norway_agestrat)
norway_eachyr <- aggregate(. ~ Serotype, norway_agestrat, FUN = sum)
norway_eachyr$Age.group <- NA
norway_lessthan18 <- norway_agestrat %>% filter(Age.group %in% c("0-1 yrs", "2-4 yrs", "5-17 yrs"))
norway_lessthan18 <- aggregate(. ~ Serotype, norway_lessthan18, FUN = sum)
norway_lessthan18$Age.group <- "< 18 yrs"
norway_morethaneq18 <- norway_agestrat %>% filter(Age.group %in% c("18-49 yrs", "50-64 yrs", "65+ yrs"))
norway_morethaneq18 <- aggregate(. ~ Serotype, norway_morethaneq18, FUN = sum)
norway_morethaneq18$Age.group <- "> 18 yrs"
norway_18 <- bind_rows(norway_lessthan18, norway_morethaneq18)

usa_agestrat <- read.csv("Q:/Technical/R/SeroReplacement/SeroRates/Raw Data/USA age strat.csv")
usa_agestrat$Age.group <- factor(usa_agestrat$Age.group, levels = c("< 5 yrs", "5 - 17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs"))
usa_agestrat <- sero15bcagg(usa_agestrat)
usa_eachyr <- aggregate(. ~ Serotype, usa_agestrat, FUN = sum)
usa_eachyr$Age.group <- NA
usa_lessthan18 <- usa_agestrat %>% filter(Age.group %in% c("< 5 yrs", "5 - 17 yrs"))
usa_lessthan18 <- aggregate(. ~ Serotype, usa_lessthan18, FUN = sum)
usa_lessthan18$Age.group <- "< 18 yrs"
usa_morethaneq18 <- usa_agestrat %>% filter(Age.group %in% c("18-49 yrs", "50-64 yrs", "65+ yrs"))
usa_morethaneq18 <- aggregate(. ~ Serotype, usa_morethaneq18, FUN = sum)
usa_morethaneq18$Age.group <- "> 18 yrs"
usa_18 <- bind_rows(usa_lessthan18, usa_morethaneq18)

ages <- unique(c(usa_18$Age.group, norway_18$Age.group, finland_18$Age.group, #italy_inf$Age.group, 
                 france_16$Age.group, oz_18$Age.group))
children <- ages[which(ages %in% c("< 18 yrs", "<16 yrs", "0-4 yrs"))]
adults <- ages[which(ages %in% c("> 18 yrs", ">16 yrs"))]

# disease cases over time df
melt.oz <- totdf(oz_18, "Australia")
melt.fin <- totdf(finland_18, "Finland")
melt.fra <- totdf(france_agestrat, "France")
# melt.ita <- totdf(italy_inf, "Italy")
melt.nor <- totdf(norway_18, "Norway")
melt.usa <- totdf(usa_18, "USA")

melt.tot <- rbind(melt.oz, melt.fin, melt.fra, #melt.ita, 
                  melt.nor, melt.usa)
melt.tot$pcv_era <- factor(melt.tot$pcv_era, levels = c("Pre-PCV", "Pre-PCV10", "Pre-PCV13", "Post-PCV"))
levels(melt.tot$pcv_era)[match("Pre-PCV10", levels(melt.tot$pcv_era))] <- "Pre-PCV10/13"
levels(melt.tot$pcv_era)[match("Pre-PCV13", levels(melt.tot$pcv_era))] <- "Pre-PCV10/13"

# percentage of each serotype over time df

# australia
aus_less18.freq <- freqdf(oz_lessthan18, "Australia")
aus_more18.freq <- freqdf(oz_morethaneq18, "Australia")

# finland
fin_less18.freq <- freqdf(finland_lessthan18, "Finland")
fin_more18.freq <- freqdf(finland_morethaneq18, "Finland")

# france
france_less16.freq <- freqdf(france_lessthan16, "France")
france_more16.freq <- freqdf(france_morethan16, "France")

# italy
# italy_inf.freq <- freqdf(italy_inf, "Italy")

# norway
norway_less18.freq <- freqdf(norway_lessthan18, "Norway")
norway_more18.freq <- freqdf(norway_morethaneq18, "Norway")

# usa
usa_less18.freq <- freqdf(usa_lessthan18, "USA")
usa_more18.freq <- freqdf(usa_morethaneq18, "USA")

tot.freq <- rbind(aus_less18.freq, aus_more18.freq, fin_less18.freq, fin_more18.freq, france_less16.freq, france_more16.freq, #italy_inf.freq, 
                  norway_less18.freq, norway_more18.freq, usa_less18.freq, usa_more18.freq)

tot.freq$pcv_era <- factor(tot.freq$pcv_era, levels = c("Pre-PCV", "Pre-PCV10", "Pre-PCV13", "Post-PCV"))
levels(tot.freq$pcv_era)[match("Pre-PCV10", levels(tot.freq$pcv_era))] <- "Pre-PCV10/13"
levels(tot.freq$pcv_era)[match("Pre-PCV13", levels(tot.freq$pcv_era))] <- "Pre-PCV10/13"

##########################################################################################################################################################
########### OVERALL IPD INCIDENCE PLOT ###################################################################################################################

colnames(pop) <- c("Country", gsub("^X", "", colnames(pop)[2:ncol(pop)]))
overall_inc_df <- melt(pop)
colnames(overall_inc_df) <- c("country", "year", "incidence")

assign_pcvera <- function(df) {
  PCVyrs <- PCVintro[which(rownames(PCVintro) == unique(df$country)),]
  PCVyrs <- as.integer(PCVyrs[!is.na(PCVyrs)])
  PCVera.spec <- PCVera[,which(colnames(PCVera) == unique(df$country))]
  PCVera.spec <- PCVera.spec[!is.na(PCVera.spec)]
  
  df$year <- as.numeric(levels(df$year))[df$year]
  
  prepcv_ind <- which(df$year <= PCVyrs[1])
  prepcv_high <- which(df$year > PCVyrs[1] & df$year <= PCVyrs[2])
  postpcv_ind <- which(df$year > PCVyrs[2])

  if(length(PCVyrs) < 2) {postpcv_ind <- which(df$year > PCVyrs[1])}
  
  df$pcv_era <- NA
  df$pcv_era[prepcv_ind] <- as.character(PCVera.spec[1])
  df$pcv_era[prepcv_high] <- as.character(PCVera.spec[2])
  df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[3])
  
  if(length(PCVyrs) < 2) {df$pcv_era[postpcv_ind] <- as.character(PCVera.spec[2])}
  return(df)
}

overall_inc_df <- lapply(unique(overall_inc_df$country), function(x) assign_pcvera(overall_inc_df %>% filter(country == x)))
overall_inc_df <- bind_rows(overall_inc_df)
overall_inc_df$pcv_era <- factor(overall_inc_df$pcv_era, levels = c("Pre-PCV", "Pre-PCV10", "Pre-PCV13", "Post-PCV"))
levels(overall_inc_df$pcv_era)[match("Pre-PCV10", levels(overall_inc_df$pcv_era))] <- "Pre-PCV10/13"
levels(overall_inc_df$pcv_era)[match("Pre-PCV13", levels(overall_inc_df$pcv_era))] <- "Pre-PCV10/13"

ggplot(overall_inc_df,aes(x = year, y = incidence, group = country, colour = country, shape = pcv_era)) + geom_point() + geom_line() +
  labs(x = "Year", y = "Overall IPD incidence / 100,000 people") + scale_shape_manual(values = c(1, 15, 2)) + labs(shape = "PCV era", colour = "Country") + 
  theme_bw() # pdf dim 4 x 8

# usa_pre13 <- overall_inc_df %>% filter(Country == "USA") %>% filter(variable %in% 2005:2010)
# usa_post13 <- overall_inc_df %>% filter(Country == "USA") %>% filter(variable %in% 2010:2013)
# usapre13.mod <- usa_pre13 %>% mutate(model = purrr::map(data, ~lm(value ~ variable, data = .)))
# rates <- usapre13.mod %>% filter(!grepl('(Intercept)', term))

# not used: IPD incidence calculation of datasets
# oz_totIPD <- colSums(australia_agestrat[,3:ncol(australia_agestrat)]) 
# pop[which(pop$Population == "Australia"), c('X1999', 'X2000')] <- NA
# finland_totIPD <- colSums(finland_agestrat[,3:ncol(finland_agestrat)])
# france_totIPD <- colSums(france_agestrat[,3:ncol(france_agestrat)])
# #italy_totIPD <- colSums(italy_eachyr[,3:ncol(italy_eachyr)])
# norway_totIPD <- colSums(norway_agestrat[,3:ncol(norway_agestrat)])
# usa_totIPD <- colSums(usa_agestrat[,3:ncol(usa_agestrat)])
# 
# Aus_inc <- (oz_totIPD/(pop[which(pop$Population == "Australia"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# Fin_inc <- (finland_totIPD/(pop[which(pop$Population == "Finland"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# Fran_inc <- (france_totIPD/(pop[which(pop$Population == "France"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# #Ita_inc <- (italy_totIPD/(pop[which(pop$Population == "Italy"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# Nor_inc <- (norway_totIPD/(pop[which(pop$Population == "Norway"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# usapop <- pop[which(pop$Population == "USA"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)
# usapop <- usapop[,2:10] # remove years that are not included in IPD data
# USA_inc <- (usa_totIPD/usapop)*100000
# #USA_inc <- (usa_totIPD/(pop[which(pop$Population == "USA"),] %>% select_if(function(x) is.na(x) == FALSE & map(x, is.factor) == FALSE)))*100000
# 
# overall_inc1 <- #full_join(Ita_inc, Aus_inc) 
#   Aus_inc %>% mutate(Country = c("Australia"))
# overall_inc2 <- full_join(Fran_inc, Fin_inc) %>% mutate(Country = c("France", "Finland"))
# overall_inc3 <- full_join(USA_inc, Nor_inc) %>% mutate(Country = c("USA", "Norway"))
# overall_inc <- bind_rows(overall_inc1, overall_inc2, overall_inc3)
# colnames(overall_inc) <- c(gsub("^X", "", colnames(overall_inc)[1:ncol(overall_inc)-1]), "Country")
# overall_inc_df <- melt(overall_inc)
# colnames(overall_inc_df) <- c("country", "year", "incidence")

##########################################################################################################################################################
########### ODDS RATIO CALCULATION AND PLOTTING ##########################################################################################################

setwd("Q:/Technical/R/SeroReplacement/OddsRatio")

cleanORdat.agestrat <- function(dat, df_nest) { # clean age-stratified OR data
  df <- data.frame(t(apply(rbind(dat), 2, FUN = unlist)))
  colnames(df) <- c("OR", "CI.low", "CI.high")
  df$Age.group <- df_nest$Age.group
  varnam <- deparse(substitute(dat))
  df$Category <-  strsplit(varnam,'_')[[1]][1]
  df$Vaccine <- strsplit(varnam,'_')[[1]][2]
  
  df$Category <- factor(df$Category, levels = c("VT7","VT10","VT13","NVT"))
  df$Vaccine <- factor(df$Vaccine, levels = c("PCV7","PCV10","PCV13"))
  return(df)
}

eachdf_OR <- function(df, country) { # fn to estimate OR for each age group in a country
  yrs <- gsub("^X", "", colnames(df)[3:ncol(df)])
  colnames(df) <- c("Serotype", "Age.group", yrs)
  
  pre7_yr <- PCVintro$PCV7intro[which(rownames(PCVintro) == country)]
  pre10_yr <- PCVintro$PCV10intro[which(rownames(PCVintro) == country)]
  pre13_yr <- PCVintro$PCV13intro[which(rownames(PCVintro) == country)]
  
  pre7_cols <- colnames(df)[2 + which(colnames(df[3:ncol(df)]) <= pre7_yr)]
  pre10_ind <- 2 + which(colnames(df[3:ncol(df)]) <= pre10_yr)
  pre10_cols <- colnames(df)[pre10_ind[!colnames(df)[pre10_ind] %in% pre7_cols]]
  pre13_ind <- 2 + which(colnames(df[3:ncol(df)]) <= pre13_yr)
  pre13_cols <- colnames(df)[pre13_ind[!colnames(df)[pre13_ind] %in% pre7_cols]]
  post_cols <- colnames(df)[(max(pre10_ind,pre13_ind)+1):ncol(df)]
  
  df_nest <- df %>% group_by(Age.group) %>% nest()
  n.agegrps <- nrow(df_nest)
  
  VT7_index <- lapply(df_nest$data, function(x) which(x$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F")))
  VT10_index <- lapply(df_nest$data, function(x) which(x$Serotype %in% c("1", "5", "7F")))
  VT13_index <- lapply(df_nest$data, function(x) which(x$Serotype %in% c("3", "6A", "19A")))
  NVT_index <- lapply(df_nest$data, function(x) which(!x$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F", "1", "5", "7F", "3", "6A", "19A")))
  
  # VT7
  
  VT7_pre7 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT7_index[[x]],colnames(df_nest$data[[x]]) %in% pre7_cols]),error = function(e) 0))
  VT7_pre10 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT7_index[[x]],colnames(df_nest$data[[x]]) %in% pre10_cols]),error = function(e) 0))
  VT7_pre13 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT7_index[[x]],colnames(df_nest$data[[x]]) %in% pre13_cols]),error = function(e) 0))
  VT7_post <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT7_index[[x]],colnames(df_nest$data[[x]]) %in% post_cols]),error = function(e) 0))
  
  # VT10
  
  VT10_pre7 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT10_index[[x]],colnames(df_nest$data[[x]]) %in% pre7_cols]),error = function(e) 0))
  VT10_pre10 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT10_index[[x]],colnames(df_nest$data[[x]]) %in% pre10_cols]),error = function(e) 0))
  VT10_pre13 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT10_index[[x]],colnames(df_nest$data[[x]]) %in% pre13_cols]),error = function(e) 0))
  VT10_post <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT10_index[[x]],colnames(df_nest$data[[x]]) %in% post_cols]),error = function(e) 0))
  
  # VT13
  
  VT13_pre7 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT13_index[[x]],colnames(df_nest$data[[x]]) %in% pre7_cols]),error = function(e) 0))
  VT13_pre10 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT13_index[[x]],colnames(df_nest$data[[x]]) %in% pre10_cols]),error = function(e) 0))
  VT13_pre13 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT13_index[[x]],colnames(df_nest$data[[x]]) %in% pre13_cols]),error = function(e) 0))
  VT13_post <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][VT13_index[[x]],colnames(df_nest$data[[x]]) %in% post_cols]),error = function(e) 0))
  
  # NVT
  
  NVT_pre7 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][NVT_index[[x]],colnames(df_nest$data[[x]]) %in% pre7_cols]),error = function(e) 0))
  NVT_pre10 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][NVT_index[[x]],colnames(df_nest$data[[x]]) %in% pre10_cols]),error = function(e) 0))
  NVT_pre13 <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][NVT_index[[x]],colnames(df_nest$data[[x]]) %in% pre13_cols]),error = function(e) 0))
  NVT_post <- lapply(1:n.agegrps, function(x) tryCatch(sum(df_nest$data[[x]][NVT_index[[x]],colnames(df_nest$data[[x]]) %in% post_cols]),error = function(e) 0))
  
  # average incidence for each period
  
  
  
  
  # OR VT7 vs all other serogroup categories for each vaccine
  VT7_PCV7_10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT7_pre10[[x]], VT7_pre7[[x]], sum(VT10_pre10[[x]], VT13_pre10[[x]], NVT_pre10[[x]]),
                                                                     sum(VT10_pre7[[x]], VT13_pre7[[x]], NVT_pre7[[x]]))), 
                           function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT7_PCV7_13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT7_pre13[[x]], VT7_pre7[[x]],  sum(VT10_pre13[[x]], VT13_pre13[[x]], NVT_pre13[[x]]),
                                                                     sum(VT10_pre7[[x]], VT13_pre7[[x]], NVT_pre7[[x]]))),
                           function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT7_PCV10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT7_post[[x]], VT7_pre10[[x]], sum(VT10_post[[x]], VT13_post[[x]], NVT_post[[x]]),
                                                                   sum(VT10_pre10[[x]], VT13_pre10[[x]], NVT_pre10[[x]]))),
                         function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT7_PCV13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT7_post[[x]], VT7_pre13[[x]], sum(VT10_post[[x]], VT13_post[[x]], NVT_post[[x]]),
                                                                   sum(VT10_pre13[[x]], VT13_pre13[[x]], NVT_pre13[[x]]))),
                         function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  
  VT7_PCV7_10 <- cleanORdat.agestrat(VT7_PCV7_10_OR, df_nest)
  VT7_PCV7_13 <- cleanORdat.agestrat(VT7_PCV7_13_OR, df_nest)
  VT7_PCV10 <- cleanORdat.agestrat(VT7_PCV10_OR, df_nest)
  VT7_PCV13 <- cleanORdat.agestrat(VT7_PCV13_OR, df_nest)
  
  # OR VT10 vs all other serogroup categories for each vaccine
  VT10_PCV7_10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT10_pre10[[x]], VT10_pre7[[x]], sum(VT7_pre10[[x]], VT13_pre10[[x]], NVT_pre10[[x]]), 
                                                                      sum(VT7_pre7[[x]], VT13_pre7[[x]], NVT_pre7[[x]]))),
                            function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT10_PCV7_13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT10_pre13[[x]], VT10_pre7[[x]], sum(VT7_pre13[[x]], VT13_pre13[[x]], NVT_pre13[[x]]), 
                                                                      sum(VT7_pre7[[x]], VT13_pre7[[x]], NVT_pre7[[x]]))),
                            function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT10_PCV10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT10_post[[x]], VT10_pre10[[x]], sum(VT7_post[[x]], VT13_post[[x]], NVT_post[[x]]), 
                                                                    sum(VT7_pre10[[x]], VT13_pre10[[x]], NVT_pre10[[x]]))),
                          function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT10_PCV13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT10_post[[x]], VT10_pre13[[x]], sum(VT7_post[[x]], VT13_post[[x]], NVT_post[[x]]), 
                                                                    sum(VT7_pre13[[x]], VT13_pre13[[x]], NVT_pre13[[x]]))),
                          function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  
  VT10_PCV7_10 <- cleanORdat.agestrat(VT10_PCV7_10_OR, df_nest)
  VT10_PCV7_13 <- cleanORdat.agestrat(VT10_PCV7_13_OR, df_nest)
  VT10_PCV10 <- cleanORdat.agestrat(VT10_PCV10_OR, df_nest)
  VT10_PCV13 <- cleanORdat.agestrat(VT10_PCV13_OR, df_nest)
  
  # OR VT13 vs all other serogroup categories for each vaccine
  VT13_PCV7_10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT13_pre10[[x]], VT13_pre7[[x]], sum(VT7_pre10[[x]], VT10_pre10[[x]], NVT_pre10[[x]]), 
                                                                      sum(VT7_pre7[[x]], VT10_pre7[[x]], NVT_pre7[[x]]))),
                            function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT13_PCV7_13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT13_pre13[[x]], VT13_pre7[[x]], sum(VT7_pre13[[x]], VT10_pre13[[x]], NVT_pre13[[x]]), 
                                                                      sum(VT7_pre7[[x]], VT10_pre7[[x]], NVT_pre7[[x]]))),
                            function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT13_PCV10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT13_post[[x]], VT13_pre10[[x]], sum(VT7_post[[x]], VT10_post[[x]], NVT_post[[x]]), 
                                                                    sum(VT7_pre10[[x]], VT10_pre10[[x]], NVT_pre10[[x]]))),
                          function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  VT13_PCV13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(VT13_post[[x]], VT13_pre13[[x]], sum(VT7_post[[x]], VT10_post[[x]], NVT_post[[x]]), 
                                                                    sum(VT7_pre13[[x]], VT10_pre13[[x]], NVT_pre13[[x]]))),
                          function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  
  VT13_PCV7_10 <- cleanORdat.agestrat(VT13_PCV7_10_OR, df_nest)
  VT13_PCV7_13 <- cleanORdat.agestrat(VT13_PCV7_13_OR, df_nest)
  VT13_PCV10 <- cleanORdat.agestrat(VT13_PCV10_OR, df_nest)
  VT13_PCV13 <- cleanORdat.agestrat(VT13_PCV13_OR, df_nest)
  
  # OR NVT vs all other serogroup categories for each vaccine
  NVT_PCV7_10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(NVT_pre10[[x]], NVT_pre7[[x]], sum(VT7_pre10[[x]], VT10_pre10[[x]], VT13_pre10[[x]]), 
                                                                     sum(VT7_pre7[[x]], VT10_pre7[[x]], VT13_pre7[[x]]))),
                           function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  NVT_PCV7_13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(NVT_pre13[[x]], NVT_pre7[[x]], sum(VT7_pre13[[x]], VT10_pre13[[x]], VT13_pre13[[x]]), 
                                                                     sum(VT7_pre7[[x]], VT10_pre7[[x]], VT13_pre7[[x]]))),
                           function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  NVT_PCV10_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(NVT_post[[x]], NVT_pre10[[x]], sum(VT7_post[[x]], VT10_post[[x]], VT13_post[[x]]), 
                                                                   sum(VT7_pre10[[x]], VT10_pre10[[x]], VT13_pre10[[x]]))),
                         function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  NVT_PCV13_OR <- lapply(lapply(1:n.agegrps, function(x) oddsratio(NVT_post[[x]], NVT_pre13[[x]], sum(VT7_post[[x]], VT10_post[[x]], VT13_post[[x]]), 
                                                                   sum(VT7_pre13[[x]], VT10_pre13[[x]], VT13_pre13[[x]]))),
                         function(r) data.frame(OR = r$estimate, CIlow = r$conf.int[1], CIhigh = r$conf.int[2]))
  
  NVT_PCV7_10 <- cleanORdat.agestrat(NVT_PCV7_10_OR, df_nest)
  NVT_PCV7_13 <- cleanORdat.agestrat(NVT_PCV7_13_OR, df_nest)
  NVT_PCV10 <- cleanORdat.agestrat(NVT_PCV10_OR, df_nest)
  NVT_PCV13 <- cleanORdat.agestrat(NVT_PCV13_OR, df_nest)
  
  OR_df <- rbind(VT7_PCV7_10, VT7_PCV7_13, VT7_PCV10, VT7_PCV13, VT10_PCV7_10, VT10_PCV7_13, VT10_PCV10, VT10_PCV13, VT13_PCV7_10, VT13_PCV7_13, VT13_PCV10, VT13_PCV13,
                 NVT_PCV7_10, NVT_PCV7_13, NVT_PCV10, NVT_PCV13)
  OR_df <- OR_df[complete.cases(OR_df),]
  
  return(OR_df)
}

OR_plot <- function(df_OR, country) { ## Plot OR by age group
  df_OR$Age.group <- factor(df_OR$Age.group, levels = rev(levels(df_OR$Age.group)))
  ggplot(data = df_OR, aes(x = Age.group, y = OR, ymin = CI.low, ymax = CI.high, shape = factor(Vaccine, levels = rev(levels(Vaccine))))) + 
    geom_pointrange(aes(col = factor(Category, levels = rev(levels(Category)))), position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c(15,16,17)) +
    geom_hline(yintercept = 1, color = "black") + scale_y_continuous(trans = 'log10') +
    labs(x = "Age group", y = "Odds Ratio", color = "Category", shape = "Vaccine", title = country)+#paste(country, "age-stratified odds ratio", sep = " ")) +
    coord_flip() + #flip axes/graph
    guides(colour = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
    scale_color_manual(values=c(VT7 = "#35B779FF", VT10= "#31688EFF", VT13 = "#E69F00", NVT = "#440154FF")) +
      #VT7 = "#70AD47", VT10 = "#4472C4", VT13 = "#ED7D31", NVT = "darkgrey")) + # PCV7 green, PCV10 blue, PCV13 orange
    theme_bw() + theme_light()
  
  #ggplot(data = df_OR) +
  #geom_errorbar(aes(x = Age.group, ymin = CI.low, ymax = CI.high, color = Category), width =0.01) + #error bars
  #geom_point(mapping = aes(x = Age.group, y = OR, shape = Vaccine, color = Category), size = 2.5) +
  #scale_shape_manual(values= c(15, 16, 17)) + #change the shapes of Serogroup
  #geom_hline(yintercept = 1, color = "black") + scale_y_continuous(trans = 'log10') + #draw a horizontal line at y=1 for the OR to indicate signif or not
  #labs(x = "Age group", y = "Odds Ratio", title = paste(country, "age-stratified odds ratio", sep = " ")) +
  #coord_flip() + #flip axes/graph
  #scale_color_manual(values=c("#70AD47", "#4472C4", "#ED7D31", "black")) + # PCV7 green, PCV10 blue, PCV13 orange
  #theme_bw() + theme_light()
}

OR_plot <- function(df_OR, country) { ## Plot OR by age group
  df_OR$Age.group <- factor(df_OR$Age.group, levels = rev(levels(df_OR$Age.group)))
  ggplot(data = df_OR, aes(x = OR, y = Age.group, col = factor(Category, levels = rev(levels(Category))), group = factor(Category, levels = rev(levels(Category))))) + 
    geom_errorbarh(aes(xmin = CI.low, xmax = CI.high), height = 0, position = position_dodgev(height = 0.5)) +
    geom_point(aes(shape = factor(Vaccine, levels = rev(levels(Vaccine)))), position = position_dodgev(height = 0.5)) +
    scale_shape_manual(values = c(PCV13 = 15, PCV7 = 1, PCV10 = 17)) +
    geom_vline(xintercept = 1, color = "black") + scale_x_continuous(trans = 'log10') +
    labs(y = "Age group", x = "Odds Ratio", color = "Category", shape = "Vaccine", title = country)+
    guides(colour = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
    scale_color_manual(values= cols)+#c(VT7 = "#70AD47", VT10 = "#4472C4", VT13 = "#ED7D31", NVT = "darkgrey")) + # PCV7 green, PCV10 blue, PCV13 orange
    theme_bw() + theme_light()
}


australia_OR <- eachdf_OR(australia_agestrat, "Australia")
australia_plot <- OR_plot(australia_OR, "Australia") # pdf dim 5 x 7 in

finland_OR <- eachdf_OR(finland_agestrat, "Finland")
finland_plot <- OR_plot(finland_OR, "Finland") # pdf dim 5 x 7 in

fakeDSforleg <- rbind(finland_OR, australia_OR)
fakeplot <- OR_plot(fakeDSforleg, "Australia") 
legendgrob <- g_legend(fakeplot+theme(legend.box = "horizontal"))

france_OR <- eachdf_OR(france_agestrat, "France")
france_plot <- OR_plot(france_OR, "France") # pdf dim 5 x 7 in

#italy_OR <- eachdf_OR(italy_agestrat, "Italy")
#italy_plot <- OR_plot(italy_OR, "Italy") # pdf dim 5 x 7 in

norway_OR <- eachdf_OR(norway_agestrat, "Norway")
norway_plot <- OR_plot(norway_OR, "Norway") # pdf dim 5 x 7 in

usa_OR <- eachdf_OR(usa_agestrat, "USA")
usa_plot <- OR_plot(usa_OR, "USA") # pdf dim 5 x 7 in

grid.arrange(australia_plot + theme(legend.position = "none"), finland_plot+ theme(legend.position = "none"), 
             france_plot+ theme(legend.position = "none"), #italy_plot, 
             norway_plot+ theme(legend.position = "none"), usa_plot+ theme(legend.position = "none"), 
             legend= legendgrob, ncol = 3, nrow = 2) # pdf dim 6 x 9

##########################################################################################################################################################
########### SIMPSON'S DIVERSITY INDEX ####################################################################################################################

setwd("Q:/Technical/R/SeroReplacement/RankFreq - Simpson Diversity")

calc_SimpsDiv <- function(isolates){ # function to calculate simpson's diversity given the list of isolates
  N <- sum(isolates)  #sum of all isolates
  freq <- isolates/N #as.data.frame(lapply(isolates, function(x) {x/N}))
  sq <- freq^2 #as.data.frame(lapply(freq, function(x) {x^2})) # square each frequency
  totsq <- sum(sq) # sum of squares
  cb <- freq^3 #as.data.frame(lapply(freq, function(x) {x^3})) # cube each frequency
  totcb <- sum(cb) # sum of cubes
  num <- isolates*(isolates-1) #as.data.frame(lapply(isolates, function(x) {x*(x-1)}))   # numerator = n*(n-1)
  var <- (4/N) * (totcb - (totsq ^ 2)) # variance^2 as per Grundmann et al 2001 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC88515/
  div_ind <- 1 - (sum(num)/(N*(N-1))) # diversity index
  CI_high <- div_ind + 2*(sqrt(var)) # upper CI
  CI_low <- div_ind - 2*(sqrt(var)) # lower CI
  simpsdivindx <- data.frame(div_ind, CI_low, CI_high, var)
  return(simpsdivindx)
}

plot_SDI_strat <- function(df, country) { # plot Simpson's Diversity Index for age group over time
  max.yr <- df$Year[length(df$Year)]
  ggplot(data = df, aes(x = Year, y = div_ind)) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.01) +
    geom_point() +
    geom_point(data = df[which(df$Year == PCVintro[country,]$PCV7intro),], aes(x = Year, y = div_ind), colour = "#35B779FF") + #"#70AD47") +
    geom_point(data = df[which(df$Year == PCVintro[country,]$PCV10intro),], aes(x = Year, y = div_ind), colour = "#31688EFF") + #"#4472C4") +
    geom_point(data = df[which(df$Year == PCVintro[country,]$PCV13intro),], aes(x = Year, y = div_ind), colour = "#E69F00") + #"#ED7D31") +
    labs(x = "Year", y = "Simpson's Diversity Index", title = paste(country, df$Age.group[1], "Simpson Diversity Index", sep = " ")) +
    theme_light()
}

plot_SDI_strat_combo <- function(df, country) { # plot Simpson's Diversity Index for age group over time
  df$Year <- as.numeric(levels(df$Year))[df$Year]
  startyr <- min(df$Year)
  max.yr <- max(df$Year)
  PCV7intro <- PCVintro[country,]$PCV7intro
  PCV10intro <- PCVintro[country,]$PCV10intro
  PCV13intro <- PCVintro[country,]$PCV13intro
  df$Age.group <- factor(df$Age.group, labels = c('Children', 'Adults')) #comment out for age stratified
  
  rectang <- data.frame(x1 = c(PCV7intro, PCV7intro, PCV10intro, PCV13intro), 
                        x2 = c(PCV10intro, PCV13intro, max.yr, max.yr), 
                        y1 = c(-Inf, -Inf, -Inf, -Inf), y2 = c(Inf, Inf, Inf, Inf), 
                        t = c("PCV7", "PCV7", "PCV10", "PCV13"))
  rectang$t <- factor(rectang$t, levels = c("PCV7", "PCV10", "PCV13"), labels = c('PCV7', 'PCV10/13', 'PCV10/13'))
  rectang <- na.omit(rectang)
  rectang.col <- c("PCV7"="grey65", "PCV7"="grey65", "PCV10"="grey20", "PCV13"="grey20", "PCV10/13" = "grey20")
  
  ggplot() +
    geom_rect(data = rectang, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = t), alpha = 0.3) + 
    scale_fill_manual("PCV era", values = c(rectang.col)) + 
    geom_errorbar(data = df, aes(x = Year, ymin = CI_low, ymax = CI_high, colour = Age.group), width = 0.01) +
    geom_point(data = df, aes(x = Year, y = div_ind, colour = Age.group, shape = Age.group)) +
    #scale_shape_manual(values = c(0, 1, 2, 15, 16, 17)) +
    scale_colour_manual(values = c("#E69F00", "deepskyblue4"))+
    labs(x = "Year", y = "Simpson's Diversity Index", title = paste(country, "Simpson Diversity Index", sep = " "), shape = "Age Group", colour = "Age Group",
         fill = "PCV era") +
    coord_cartesian(xlim = c(startyr, max.yr))+ 
    theme_light()
}

age_strat <- function(df, country) { # calculate Simpson's Diversity index for each age group and year
  
  num.agegrp <- length(unique(df$Age.group))
  
  df_nest <- df %>% group_by(Age.group) %>% nest()
  simpsdiv.strat <- lapply(1:nrow(df_nest), function(r) lapply(2:ncol(df_nest$data[[r]]), function(x) calc_SimpsDiv(df_nest$data[[r]][,x])))
  
  newvarz <- lapply(1:num.agegrp, function(x) bind_rows(simpsdiv.strat[[x]]))
  yrs <- gsub("^X", "", colnames(df)[3:ncol(df)])
  newvarz <- lapply(1:length(newvarz), function(x) cbind(newvarz[[x]], Year = yrs, Age.group = unique(df$Age.group)[x]))
  
  return(newvarz)
  
}

# Australia
oz_agestratSDI <- age_strat(australia_agestrat, "Australia")
oz_agestrat_plots <- lapply(oz_agestratSDI, function(x) plot_SDI_strat(x, "Australia"))

# < 18 & >= 18 yrs age groups
oz_18_SDI.d <- age_strat(oz_18, "Australia")
oz_18_SDI <- bind_rows(oz_18_SDI.d)
oz_18_SDI_plot <- plot_SDI_strat_combo(oz_18_SDI, "Australia")

# all age groups in one graph
oz_agestratSDI_combo <- bind_rows(oz_agestratSDI)
oz_allage <- plot_SDI_strat_combo(oz_agestratSDI_combo, "Australia") # pdf dim 7 x 10 in

# SDI consolidated all age groups
australia_SDI <- age_strat(australia_eachyr, "Australia")
australia_plot <- plot_SDI_strat(australia_SDI[[1]], "Australia")

# Finland
PCVintro$PCV7intro[which(rownames(PCVintro) == "Finland")] <- 2009
finland_agestratSDI <- age_strat(finland_agestrat, "Finland")
finland_agestrat_plots <- lapply(finland_agestratSDI, function(x) plot_SDI_strat(x, "Finland"))

# < 18 & >= 18 yrs age groups
finland_18_SDI.d <- age_strat(finland_18, "finland")
finland_18_SDI <- bind_rows(finland_18_SDI.d)
finland_18_SDI_plot <- plot_SDI_strat_combo(finland_18_SDI, "Finland")

# all age groups in one graph
finland_agestratSDI_combo <- bind_rows(finland_agestratSDI)
finland_allage <- plot_SDI_strat_combo(finland_agestratSDI_combo, "Finland") # pdf dim 7 x 10 in

finland_SDI <- age_strat(finland_eachyr, "Finland")
finland_plot <- plot_SDI_strat(finland_SDI[[1]], "Finland")

# France
france_agestratSDI <- age_strat(france_agestrat, "France")
france_agestrat_plots <- lapply(france_agestratSDI, function(x) plot_SDI_strat(x, "France"))

# all age groups in one graph
france_agestratSDI_combo <- bind_rows(france_agestratSDI)
france_allage <- plot_SDI_strat_combo(france_agestratSDI_combo, "France") # pdf dim 7 x 10 in

france_SDI <- age_strat(france_eachyr, "France")
france_plot <- plot_SDI_strat(france_SDI[[1]], "France")

# # Italy
# italy_agestratSDI <- age_strat(italy_agestrat, "Italy")
# italy_agestrat_plots <- lapply(italy_agestratSDI, function(x) plot_SDI_strat(x, "Italy"))
# 
# # all age groups in one graph
# italy_agestratSDI_combo <- bind_rows(italy_agestratSDI)
# italy_allage <- plot_SDI_strat_combo(italy_agestratSDI_combo, "Italy") # pdf dim 7 x 10 in
# 
# italy_SDI <- age_strat(italy_eachyr, "Italy")
# italy_plot <- plot_SDI_strat(italy_SDI[[1]], "Italy")

# Norway
PCVintro$PCV7intro[which(rownames(PCVintro) == "Norway")] <- 2006
norway_agestratSDI <- age_strat(norway_agestrat, "Norway")
norway_agestrat_plots <- lapply(norway_agestratSDI, function(x) plot_SDI_strat(x, "Norway"))

# < 18 & >= 18 yrs age groups
norway_18_SDI.d <- age_strat(norway_18, "Norway")
norway_18_SDI <- bind_rows(norway_18_SDI.d)
norway_18_SDI_plot <- plot_SDI_strat_combo(norway_18_SDI, "Norway") # pdf dim 7 x 10 in

# all age groups in one graph
norway_agestratSDI_combo <- bind_rows(norway_agestratSDI)
norway_allage <- plot_SDI_strat_combo(norway_agestratSDI_combo, "Norway") # pdf dim 7 x 10 in

norway_SDI <- age_strat(norway_eachyr, "Norway")
norway_plot <- plot_SDI_strat(norway_SDI[[1]], "Norway")

# USA
usa_agestratSDI <- age_strat(usa_agestrat, "USA")
usa_agestrat_plots <- lapply(usa_agestratSDI, function(x) plot_SDI_strat(x, "USA"))

# < 18 & >= 18 yrs age groups
usa_18_SDI.d <- age_strat(usa_18, "USA")
usa_18_SDI <- bind_rows(usa_18_SDI.d)
usa_18_SDI_plot <- plot_SDI_strat_combo(usa_18_SDI, "USA")

# all age groups in one graph
usa_agestratSDI_combo <- bind_rows(usa_agestratSDI)
usa_allage <- plot_SDI_strat_combo(usa_agestratSDI_combo, "USA") # pdf dim 7 x 10 in

usa_SDI <- age_strat(usa_eachyr, "USA")
usa_plot <- plot_SDI_strat(usa_SDI[[1]], "USA")

## all combo plots into one
grid.arrange(oz_allage, finland_allage, france_allage, #italy_allage, 
             norway_allage, usa_allage, ncol = 2, nrow = 3) # pdf dim 9x17 in
grid.arrange(oz_18_SDI_plot+theme(legend.position = 'none'), finland_18_SDI_plot+theme(legend.position = 'none'), 
             france_allage+theme(legend.position = 'none'), norway_18_SDI_plot+theme(legend.position = 'none'), 
             usa_18_SDI_plot+theme(legend.position = 'none'), 
             ncol = 2, nrow = 3) # pdf dim 9x17 in

sharedleg <- get_legend(oz_18_SDI_plot+theme(legend.box = 'horizontal'))
SDI.path <- file.path("Q:","Technical","R","SeroReplacement","RankFreq - Simpson Diversity", paste("Fig2 - SDI children adults- leg2", ".pdf", sep=""))
pdf(SDI.path, width = 7, height = 6)
grid.arrange(oz_18_SDI_plot+theme(legend.position = 'none') + labs(y = ""), 
             finland_18_SDI_plot+theme(legend.position = 'none') + labs(y = ""), 
             france_allage+theme(legend.position = 'none'), 
             norway_18_SDI_plot+theme(legend.position = 'none') + labs(y = ""), 
             usa_18_SDI_plot+theme(legend.position = 'none') + labs(y = ""), sharedleg, ncol = 2, nrow = 3)
dev.off()

# consolidated SDI for Australia + France + Norway post-PCV7 and post-PCV13 comparison
combo_aus <- oz_18
combo_aus$country <- "Australia"
combo_aus$Age.group[combo_aus$Age.group == "< 18 yrs"] <- "children"
combo_aus$Age.group[combo_aus$Age.group == "> 18 yrs"] <- "adults"
combo_aus <- melt(combo_aus)
colnames(combo_aus) <- c("Serotype", "Age.group", "country", "Year", "Value")
combo_aus$Year <- as.numeric(gsub("^X", "", combo_aus$Year))
oz_yrs_rel_PCV7 <- combo_aus$Year - PCVintro["Australia",]$PCV7intro
oz_PCV7_first4 <- combo_aus[oz_yrs_rel_PCV7 >= 0 & oz_yrs_rel_PCV7 < 4,]
oz_PCV7_first4$Year <- oz_PCV7_first4$Year - PCVintro["Australia",]$PCV7intro
oz_yrs_rel_PCV13 <- combo_aus$Year - PCVintro["Australia",]$PCV13intro
oz_PCV13_first4 <- combo_aus[oz_yrs_rel_PCV13 >= 0 & oz_yrs_rel_PCV13 < 4,]
oz_PCV13_first4$Year <- oz_PCV13_first4$Year - PCVintro["Australia",]$PCV13intro

combo_france <- france_agestrat
combo_france$country <- "France"
combo_france$Age.group <- as.character(combo_france$Age.group)
combo_france$Age.group[combo_france$Age.group == "<16 yrs"] <- "children"
combo_france$Age.group[combo_france$Age.group == ">16 yrs"] <- "adults"
combo_france <- melt(combo_france)
colnames(combo_france) <- c("Serotype", "Age.group", "country", "Year", "Value")
combo_france$Year <- as.numeric(gsub("^X", "", combo_france$Year))
fra_yrs_rel_PCV7 <- combo_france$Year - PCVintro["France",]$PCV7intro
fra_PCV7_first4 <- combo_france[fra_yrs_rel_PCV7 >= 0 & fra_yrs_rel_PCV7 < 4,]
fra_PCV7_first4$Year <- fra_PCV7_first4$Year - PCVintro["France",]$PCV7intro
fra_yrs_rel_PCV13 <- combo_france$Year - PCVintro["France",]$PCV13intro
fra_PCV13_first4 <- combo_france[fra_yrs_rel_PCV13 >= 0 & fra_yrs_rel_PCV13 < 4,]
fra_PCV13_first4$Year <- fra_PCV13_first4$Year - PCVintro["France",]$PCV13intro

combo_nor <- norway_18
combo_nor$country <- "Norway"
combo_nor$Age.group[combo_nor$Age.group == "< 18 yrs"] <- "children"
combo_nor$Age.group[combo_nor$Age.group == "> 18 yrs"] <- "adults"
combo_nor <- melt(combo_nor)
colnames(combo_nor) <- c("Serotype", "Age.group", "country", "Year", "Value")
combo_nor$Year <- as.numeric(gsub("^X", "", combo_nor$Year))
nor_yrs_rel_PCV7 <- combo_nor$Year - 2006 # Norway PCV7 intro year
nor_PCV7_first4 <- combo_nor[nor_yrs_rel_PCV7 >= 0 & nor_yrs_rel_PCV7 < 4,]
nor_PCV7_first4$Year <- nor_PCV7_first4$Year - 2006
nor_yrs_rel_PCV13 <- combo_nor$Year - PCVintro["Norway",]$PCV13intro
nor_PCV13_first4 <- combo_nor[nor_yrs_rel_PCV13 >= 0 & nor_yrs_rel_PCV13 < 4,]
nor_PCV13_first4$Year <- nor_PCV13_first4$Year - PCVintro["Norway",]$PCV13intro

# children post-PCV7
pooled_childrenSDI_PCV7 <- rbind(oz_PCV7_first4[oz_PCV7_first4$Age.group == "children",], fra_PCV7_first4[fra_PCV7_first4$Age.group == "children",],
                               nor_PCV7_first4[nor_PCV7_first4$Age.group == "children",])
pooled_children_PCV7_agg <- aggregate(pooled_childrenSDI_PCV7$Value, by = list(Serotype = pooled_childrenSDI_PCV7$Serotype, Year = pooled_childrenSDI_PCV7$Year), FUN = sum)
children_PCV7_nest <- pooled_children_PCV7_agg %>% group_by(Year) %>% nest()
children_PCV7_SDI <- bind_rows(lapply(1:nrow(children_PCV7_nest), function(r) calc_SimpsDiv(children_PCV7_nest$data[[r]][,2])))
children_PCV7_SDI$Age.group <- "children"
children_PCV7_SDI$Vaccination <- "PCV7"
children_PCV7_SDI$Year_post_vaccination <- 0:(nrow(children_PCV7_SDI)-1)

# children post-PCV13
pooled_childrenSDI_PCV13 <- rbind(oz_PCV13_first4[oz_PCV13_first4$Age.group == "children",], fra_PCV13_first4[fra_PCV13_first4$Age.group == "children",],
                                 nor_PCV13_first4[nor_PCV13_first4$Age.group == "children",])
pooled_children_PCV13_agg <- aggregate(pooled_childrenSDI_PCV13$Value, by = list(Serotype = pooled_childrenSDI_PCV13$Serotype, Year = pooled_childrenSDI_PCV13$Year), 
                                       FUN = sum)
children_PCV13_nest <- pooled_children_PCV13_agg %>% group_by(Year) %>% nest()
children_PCV13_SDI <- bind_rows(lapply(1:nrow(children_PCV13_nest), function(r) calc_SimpsDiv(children_PCV13_nest$data[[r]][,2])))
children_PCV13_SDI$Age.group <- "children"
children_PCV13_SDI$Vaccination <- "PCV13"
children_PCV13_SDI$Year_post_vaccination <- 0:(nrow(children_PCV13_SDI)-1)

# adults post-PCV7
pooled_adultsSDI_PCV7 <- rbind(oz_PCV7_first4[oz_PCV7_first4$Age.group == "adults",], fra_PCV7_first4[fra_PCV7_first4$Age.group == "adults",],
                                 nor_PCV7_first4[nor_PCV7_first4$Age.group == "adults",])
pooled_adults_PCV7_agg <- aggregate(pooled_adultsSDI_PCV7$Value, by = list(Serotype = pooled_adultsSDI_PCV7$Serotype, Year = pooled_adultsSDI_PCV7$Year), FUN = sum)
adults_PCV7_nest <- pooled_adults_PCV7_agg %>% group_by(Year) %>% nest()
adults_PCV7_SDI <- bind_rows(lapply(1:nrow(adults_PCV7_nest), function(r) calc_SimpsDiv(adults_PCV7_nest$data[[r]][,2])))
adults_PCV7_SDI$Age.group <- "adults"
adults_PCV7_SDI$Vaccination <- "PCV7"
adults_PCV7_SDI$Year_post_vaccination <- 0:(nrow(adults_PCV7_SDI)-1)

# adults post-PCV13
pooled_adultsSDI_PCV13 <- rbind(oz_PCV13_first4[oz_PCV13_first4$Age.group == "adults",], fra_PCV13_first4[fra_PCV13_first4$Age.group == "adults",],
                                  nor_PCV13_first4[nor_PCV13_first4$Age.group == "adults",])
pooled_adults_PCV13_agg <- aggregate(pooled_adultsSDI_PCV13$Value, by = list(Serotype = pooled_adultsSDI_PCV13$Serotype, Year = pooled_adultsSDI_PCV13$Year), FUN = sum)
adults_PCV13_nest <- pooled_adults_PCV13_agg %>% group_by(Year) %>% nest()
adults_PCV13_SDI <- bind_rows(lapply(1:nrow(adults_PCV13_nest), function(r) calc_SimpsDiv(adults_PCV13_nest$data[[r]][,2])))
adults_PCV13_SDI$Age.group <- "adults"
adults_PCV13_SDI$Vaccination <- "PCV13"
adults_PCV13_SDI$Year_post_vaccination <- 0:(nrow(adults_PCV13_SDI)-1)

SDI_combined_df <- bind_rows(children_PCV7_SDI, children_PCV13_SDI, adults_PCV7_SDI, adults_PCV13_SDI)
SDI_combined_df$Age.group <- factor(SDI_combined_df$Age.group, levels = c("children", "adults"), labels = c("Children", "Adults"))
SDI_combined_df$Vaccination <- factor(SDI_combined_df$Vaccination, levels = c("PCV7", "PCV13"))
pooledcountries <- ggplot(SDI_combined_df) + geom_point(aes(x = Year_post_vaccination, y = div_ind, colour = Age.group, shape = Vaccination), size = 2) + 
  geom_errorbar(aes(x = Year_post_vaccination, ymin = CI_low, ymax = CI_high, colour = Age.group), width = 0.05) + 
  geom_line(aes(x = Year_post_vaccination, y = div_ind, colour = Age.group, linetype = Vaccination)) + 
  scale_colour_manual(values = c("#E69F00", "deepskyblue4")) + labs(colour = "Age group", shape = 'Vaccine', linetype = 'Vaccine') +
  theme_minimal() + theme_bw() + xlab("Years post-vaccination") + ylab("Simpson's Diversity Index") 

pooled.path <- file.path("Q:","Technical","R","SeroReplacement","RankFreq - Simpson Diversity", paste("Fig4 - pooled countries", ".pdf", sep=""))
pdf(pooled.path, width = 5, height = 3)
print(pooledcountries)
dev.off()

##########################################################################################################################################################
########### COMPARISON MAJOR IPD-CAUSING SEROTYPES #######################################################################################################

serorank <- function(df, country){ # function that takes dataframe and age group and shows ranking and change in rank over vaccination periods
  colnames(df) <- c("Serotype", "Age.group", gsub("^X", "", colnames(df)[3:ncol(df)]))
  
  PCV7intro <- PCVintro[country,]$PCV7intro
  PCVhigh <- as.numeric(PCVintro[country, 2:3])
  PCVhigh.intro <- PCVhigh[!is.na(PCVhigh)]
  max.yr <- colnames(df[ncol(df)])
  
  yrs <- c(PCV7intro, PCVhigh.intro, max.yr)
  yrs <- yrs[!is.na(yrs)]
  pcv_ind <- unlist(lapply(yrs, function(x) which(colnames(df) == x)))

  agegrp_vcn <- data.frame(Serotype = df$Serotype, 
                           pre.vaccn = rowSums(df[3:pcv_ind[1]]), 
                           pre.pcvhigh = rowSums(df[(pcv_ind[1]+1):pcv_ind[2]]))
  if(length(pcv_ind) > 2) {agegrp_vcn$post.pcvhigh = rowSums(df[(pcv_ind[2]+1):ncol(df)])}
  agegrp_vcn <<- agegrp_vcn
  
  top10 <- lapply(2:ncol(agegrp_vcn), function(x) agegrp_vcn[order(-agegrp_vcn[,x]),]$Serotype)
  topSero <- lapply(top10, function(x) head(x, 10))
  if((ncol(agegrp_vcn)-1)<3) {topSero[[3]] <- NA}
  
  o <- 0.05

  PCV_eras <- PCVera[,which(colnames(PCVera) == country)]
  PCV_eras <- PCV_eras[!is.na(PCV_eras)]
  if (is.na(PCV7intro) == FALSE & PCV7intro < colnames(df)[3]) {PCV_eras <- PCV_eras[-1]} # remove pre-PCV era if data starts after PCV implemented
  PCV_eras <<- PCV_eras
  
  max.len <- max(length(topSero[[1]]), length(topSero[[2]]),length(topSero[[3]]))
  topSero[[1]] <- c(as.character(PCV_eras[1]), as.character(topSero[[1]]), rep(NA, max.len - length(topSero[[1]])))
  topSero[[2]] <- c(as.character(PCV_eras[2]), as.character(topSero[[2]]), rep(NA, max.len - length(topSero[[2]])))
  topSero[[3]] <- c(as.character(PCV_eras[3]), as.character(topSero[[3]]), rep(NA, max.len - length(topSero[[3]])))
  
  # x-axis alignment of text
  x_ind <- c(1.4, 2, 2.6)
  if (length(PCV_eras) > 2) {x_ind} else {x_ind[2] <- x_ind[3]}
    
  DF <- data.frame(x = c(rep(x_ind[1], length(topSero[[1]])), rep(x_ind[2], length(topSero[[2]])), rep(x_ind[3], length(topSero[[3]]))),
                   x1 = c(rep(x_ind[1]+o, length(topSero[[1]])), rep(x_ind[2], length(topSero[[2]])), rep(x_ind[3]-o, length(topSero[[3]]))),
                   y = c(rev(seq_along(topSero[[1]])), rev(seq_along(topSero[[2]])), rev(seq_along(topSero[[3]]))),
                   g = c(as.character(topSero[[1]]), as.character(topSero[[2]]),as.character(topSero[[3]])))
  DF$x1[is.na(DF$g)] <- NA
  DF$x1[c(1, 12, 23)] <- NA
  DF$cat[DF$g %in% VT7] <- "VT7"
  DF$cat[DF$g %in% VT10] <- "VT10"
  DF$cat[DF$g %in% VT13] <- "VT13"
  DF$cat[!(DF$g %in% c(VT7, VT10, VT13))] <- "NVT"
  DF$cat[DF$g %in% PCV_eras] <- "pcv_eras"
  
  #path.cols <- c("VT7" = "#70AD47", "VT10" = "#4472C4", "VT13" = "#ED7D31", "NVT" = "black", "pcv_eras" = "white")
  path.cols <- c("VT7" = "#35B779FF", "VT10"= "#31688EFF", "VT13" = "#E69F00", "NVT" = "#440154FF", "pcv_eras" = "white")
  
  titleseq <- data.frame(x = unique(c(x_ind[1], x_ind[2], x_ind[3])), y = c(rep(DF$y[1], length(PCV_eras))))
  
  ggplot(DF, aes(x=x, y=y, group=g, label=g)) +
    geom_path(aes(x=x1, colour = cat), size=1.2, linetype = "dotted") +
    scale_colour_manual(values = path.cols) +
    
    geom_label(aes(fill = cat),colour = "white", fontface = "bold", hjust = "inward") +
    scale_fill_manual(values = path.cols) +
    
    theme_minimal() +
    ggtitle(paste(country, unique(df$Age.group), sep = " ")) +
    scale_y_continuous(minor_breaks = seq(0,11, 1.5)) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    annotate("text", x = titleseq$x, y = titleseq$y, label = PCV_eras, hjust = "inward", fontface = "bold")
  }
  
rankfreqplot <- function(df, country, agegroup, PCV_eras, positionheight) { # function that plots rank frequency distribution of IPD cases
  df$cat[df$Serotype %in% VT7] <- "VT7"
  df$cat[df$Serotype %in% VT10] <- "VT10"
  df$cat[df$Serotype %in% VT13] <- "VT13"
  df$cat[!(df$Serotype %in% c(VT7, VT10, VT13))] <- "NVT"
  df$cat <- factor(df$cat, levels = c("VT7", "VT10", "VT13", "NVT"))
  df <- df[,c(1,ncol(df),2:(ncol(df)-1))]
  
  rank <- lapply(3:ncol(df), function(x) round(rank(-df[,x])))
  
  df_rank <- lapply(1:length(rank), function(x) data.frame(serotype = df[,1], cat = df[,2], rank = rank[[x]], IPD.cases = df[,x+2]))
  df_rank <- lapply(1:length(rank), function(x) dplyr::arrange(df_rank[[x]], rank))
  
  cumfreq <- lapply(1:length(rank), function(x) cumsum(df_rank[[x]]$IPD.cases)/sum(df_rank[[x]]$IPD.cases))
  df_rank <- mapply(cbind, df_rank, "cumfreq" = cumfreq, SIMPLIFY = F)
  
  df_rank <- lapply(1:length(rank), function(x) rbind(c(NA, NA, rank = 0, IPD.cases = NA, cumfreq = 0), df_rank[[x]]))
  ymax <- unlist(lapply(1:length(rank), function(x) df_rank[[x]]$IPD.cases[2]))
  
  newrankfreqplots <- lapply(1:length(rank), function(x) {
    repeatedranks <- as.data.frame(table(df_rank[[x]][,3]))[as.data.frame(table(df_rank[[x]][,3]))$Freq > 1,]$Var1
    jitteredpts <- df_rank[[x]][df_rank[[x]][,3] %in% repeatedranks,]
    nonjitteredpts <- df_rank[[x]][!(df_rank[[x]][,3] %in% repeatedranks),]
    ggplot() + #df_rank[[x]])
      geom_point(data = nonjitteredpts, aes(x = rank, y = IPD.cases, colour = cat, group = cat)) + 
      geom_point(data = jitteredpts, aes(x = rank, y = IPD.cases, colour = cat, group = cat), position=position_dodgev(height=positionheight)) +
      geom_line(data = df_rank[[x]], aes(x = rank, y = as.numeric(cumfreq)*max(ymax)), size = 0.8, alpha = 0.7, linetype = "dashed") + 
      scale_colour_manual(breaks = c("VT7","VT10","VT13","NVT"), values = cols) + #c("VT7" = "#70AD47", "VT10" = "#4472C4", "VT13" = "#ED7D31", "NVT" = "black")) +
      scale_y_continuous(sec.axis = sec_axis(~./max(ymax))) + theme_light() + theme_bw() + # name = "Cumulative Frequency"
      labs(x = "Rank", y = "IPD cases", title = paste(country, agegroup, PCV_eras[x], sep = " "), colour = "Category")})
  
  return(newrankfreqplots)
}


# australia
aus_less18 <- oz_lessthan18 %>% filter(Serotype != "NON TYPABLE") %>% filter(Serotype != "UNTYPED")
ausless18.plot <- serorank(aus_less18, country = "Australia") # pdf dim 5x4 in
ausless18.uniquesero <- unique(ausless18.plot$data[,4])
ausless18.rank <- agegrp_vcn
ausless18.era <- PCV_eras
ausless18.rankfreqplot <- rankfreqplot(df = ausless18.rank, country = "Australia", agegroup = "< 18 yrs", PCV_eras = ausless18.era, positionheight = 50)
ausless18.alleras <- ggarrange(ausless18.rankfreqplot[[1]]+ coord_cartesian(xlim = c(0,50)), 
             ausless18.rankfreqplot[[2]]+ labs(title = ausless18.era[2]) +
               coord_cartesian(xlim = c(0,50), ylim = c(0,1000)), 
             ausless18.rankfreqplot[[3]]+ labs(title = ausless18.era[3]) + 
               coord_cartesian(xlim = c(0,50), ylim = c(0,1000)), 
             common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1) # pdf dim 3x12 in

aus_more18 <- oz_morethaneq18 %>% filter(Serotype != "NON TYPABLE") %>% filter(Serotype != "UNTYPED")
ausmore18.plot <- serorank(aus_more18, country = "Australia")
ausmore18.uniquesero <- unique(ausmore18.plot$data[,4])
ausmore18.rank <- agegrp_vcn
ausmore18.era <- PCV_eras
ausmore18.rankfreqplot <- rankfreqplot(df = ausmore18.rank, country = "Australia", agegroup = "> 18 yrs", PCV_eras = ausmore18.era, positionheight = 50)
ausmore18.alleras <- ggarrange(ausmore18.rankfreqplot[[1]]+coord_cartesian(xlim = c(0,50), ylim = c(0,1250)), 
                               ausmore18.rankfreqplot[[2]]+ labs(title = ausmore18.era[2]) +
                                 coord_cartesian(xlim = c(0,50), ylim = c(-25,1250)), 
                               ausmore18.rankfreqplot[[3]]+ labs(title = ausmore18.era[3]) +
                                 coord_cartesian(xlim = c(0,50), ylim = c(-25,1250)), 
                               common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1)# pdf dim 3x12 in

ausranks <- grid.arrange(ausless18.plot, ausmore18.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in
ausrankfreqgrob <- ggarrange(ausless18.rankfreqplot[[1]]+ labs(x = "", y = "") + coord_cartesian(xlim = c(0,50), ylim = c(0,1000)), 
          ausless18.rankfreqplot[[2]]+ labs(title = ausless18.era[2],x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(0,1000)), 
          ausless18.rankfreqplot[[3]]+ labs(title = ausless18.era[3],x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,1000)), 
          ausmore18.rankfreqplot[[1]]+ labs(x = "", y = "")+ coord_cartesian(xlim = c(0,50), ylim = c(0,1250)), 
          ausmore18.rankfreqplot[[2]]+ labs(title = ausmore18.era[2], x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(0,1250)), 
          ausmore18.rankfreqplot[[3]]+ labs(title = ausmore18.era[3], x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(0,1250)), 
          common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 2) # pdf dim 6x12 in
annotate_figure(ausrankfreqgrob, 
                left = text_grob("IPD cases", rot = 90),
                right = text_grob("Cumulative frequency", rot = 270))#,
                #fig.lab = "A", fig.lab.face = "bold") # pdf dim 4.5 x 9 in

# print cumulative freq to CSV - need df_rank to be assigned to global environment within rankfreqplot fn
#lapply(df_rank, function(x) write.table(data.frame(x), 'test.csv', append = T, sep =','))

# finland
fin_less18 <- finland_lessthan18 %>% filter(Serotype != "Unknown")
finless18.plot <- serorank(fin_less18, country = "Finland") # for finland, CHANGE PCV INTRO YEARS TO JUST HAVE PRE/POST PCV10
finless18.uniquesero <- unique(finless18.plot$data[,4])
finless18.rank <- agegrp_vcn
finless18.era <- PCV_eras
finless18.rankfreqplot <- rankfreqplot(df = finless18.rank, country = "Finland", agegroup = "< 18 yrs", PCV_eras = finless18.era, positionheight = 15)
finless18.alleras <- ggarrange(finless18.rankfreqplot[[1]] + 
                                    coord_cartesian(xlim = c(0,50), ylim = c(-25,190)), 
                                  finless18.rankfreqplot[[2]]+ labs(title = finless18.era[2]) + 
                                    coord_cartesian(xlim = c(0,50), ylim = c(-25,190)), 
                                  common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)# pdf dim 3x12 in

fin_more18 <- finland_morethaneq18 %>% filter(Serotype != "Unknown")
finmore18.plot <- serorank(fin_more18, country = "Finland")
finmore18.uniquesero <- unique(finmore18.plot$data[,4])
finmore18.rank <- agegrp_vcn
finmore18.era <- PCV_eras
finmore18.rankfreqplot <- rankfreqplot(df = finmore18.rank, country = "Finland", agegroup = "> 18 yrs", PCV_eras = finmore18.era, positionheight = 50)
finmore18.alleras <- ggarrange(finmore18.rankfreqplot[[1]] +
                                  coord_cartesian(xlim = c(0,50), ylim = c(-25,800)), 
                               finmore18.rankfreqplot[[2]]+ labs(title = finmore18.era[2]) + 
                                    coord_cartesian(xlim = c(0,50), ylim = c(-25,800)), 
                               common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)# pdf dim 3x12 in

finranks <- grid.arrange(finless18.plot, finmore18.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in
finrankfreqgrob <- ggarrange(finless18.rankfreqplot[[1]] + labs(x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(-25,190)), 
          finless18.rankfreqplot[[2]]+ labs(title = finless18.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(-25,190)), 
          finmore18.rankfreqplot[[1]] + labs(x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(-25,700)), 
          finmore18.rankfreqplot[[2]]+ labs(title = finmore18.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(-25,700)), 
          common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 2)
annotate_figure(finrankfreqgrob, 
                left = text_grob("IPD cases", rot = 90),
                right = text_grob("Cumulative frequency", rot = 270))#,
                #fig.lab = "B", fig.lab.face = "bold")

# france
france_less16 <- france_lessthan16 %>% filter(Serotype != "Other")
franceless16.plot <- serorank(france_less16, country = "France")
franceless16.uniquesero <- unique(franceless16.plot$data[,4])
franceless16.rank <- agegrp_vcn
franceless16.era <- PCV_eras
franceless16.rankfreqplot <- rankfreqplot(df = franceless16.rank, country = "France", agegroup = "< 16 yrs", PCV_eras = franceless16.era, positionheight = 30)
franceless16.alleras <- ggarrange(franceless16.rankfreqplot[[1]]+
                                    coord_cartesian(xlim = c(0,50), ylim = c(0,400)), 
                                  franceless16.rankfreqplot[[2]] + 
                                    labs(title = franceless16.era[2]) + 
                                      coord_cartesian(xlim = c(0,50), ylim = c(0,400)), 
                                  franceless16.rankfreqplot[[3]]+ labs(title = franceless16.era[3]) + 
                                      coord_cartesian(xlim = c(0,50), ylim = c(0,400)), 
                                  common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1) # pdf dim 3x12 in

france_more16 <- france_morethan16 %>% filter(Serotype != "Other")
francemore16.plot <- serorank(france_more16, country = "France")
francemore16.uniquesero <- unique(francemore16.plot$data[,4])
francemore16.rank <- agegrp_vcn
francemore16.era <- PCV_eras
francemore16.rankfreqplot <- rankfreqplot(df = francemore16.rank, country = "France", agegroup = "> 16 yrs", PCV_eras = francemore16.era, positionheight = 50)
francemore16.alleras <- ggarrange(francemore16.rankfreqplot[[1]] +
                                    coord_cartesian(xlim = c(0,50), ylim = c(0,650)),
                                  francemore16.rankfreqplot[[2]] +
                                    labs(title = francemore16.era[2]) + 
                                    coord_cartesian(xlim = c(0,50), ylim = c(0,650)), 
                                  francemore16.rankfreqplot[[3]]+ labs(title = francemore16.era[3]) + 
                                    coord_cartesian(xlim = c(0,50), ylim = c(0,650)), 
                                  common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1) # pdf dim 3x12 in

franranks <- grid.arrange(franceless16.plot, francemore16.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in
franrankfreqgrob <- ggarrange(franceless16.rankfreqplot[[1]] + labs(x = "", y = "") +
            coord_cartesian(xlim = c(0,50), ylim = c(0,400)), 
          franceless16.rankfreqplot[[2]] + 
            labs(title = franceless16.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,400)), 
          franceless16.rankfreqplot[[3]]+ labs(title = franceless16.era[3], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,400)),
          francemore16.rankfreqplot[[1]] + labs(x = "", y = "")+
            coord_cartesian(xlim = c(0,50), ylim = c(0,650)),
          francemore16.rankfreqplot[[2]] +
            labs(title = francemore16.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,650)), 
          francemore16.rankfreqplot[[3]]+ labs(title = francemore16.era[3], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,650)), 
          common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 2) # pdf dim 6x12 in
annotate_figure(franrankfreqgrob, 
                left = text_grob("IPD cases", rot = 90),
                right = text_grob("Cumulative frequency", rot = 270))#,
                #fig.lab = "C", fig.lab.face = "bold") # pdf dim 4.5 x 9
# # italy
# italy_infants <- italy_inf %>% filter(Serotype != "Others")
# italy_infants.plot <- serorank(italy_infants, country = "Italy")
# italyinf.uniquesero <- unique(italy_infants.plot$data[,4])
# italy_infants.rank <- agegrp_vcn
# italy_infants.era <- PCV_eras
# italy_infants.rankfreqplot <- rankfreqplot(df = italy_infants.rank, country = "Italy", agegroup = "Infants", PCV_eras = italy_infants.era)
# italy_infants.alleras <- grid.arrange(italy_infants.rankfreqplot[[1]]+theme(legend.position = "none") + 
#                                        coord_cartesian(xlim = c(0,50), ylim = c(0,50)), 
#                                       italy_infants.rankfreqplot[[2]] + theme(legend.position = "none") +
#                                        labs(title = italy_infants.era[2], y = element_blank()) + 
#                                        coord_cartesian(xlim = c(0,50), ylim = c(0,50)), 
#                                       italy_infants.rankfreqplot[[3]]+ labs(title = italy_infants.era[3], y = element_blank()) + 
#                                        coord_cartesian(xlim = c(0,50), ylim = c(0,50)), 
#                                      ncol = 3, nrow = 1) # pdf dim 3x12 in
# 
# # italy_elderly plot doesn't work because there's no pre/post vaccination data for the elderly --> potentially take out?
# #italy_elderly <- italy_eld %>% filter(Serotype != "Others")
# #italy_elderly.plot <- serorank(italy_elderly, country = "Italy")
# #grid.arrange(italy_infants.plot, italy_elderly.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in

# norway
nor_less18 <- norway_lessthan18 %>% filter(Serotype != "Empty") %>% filter(Serotype != "NT") 
norless18.plot <- serorank(nor_less18, country = "Norway")
norless18.uniquesero <- unique(norless18.plot$data[,4])
norless18.rank <- agegrp_vcn
norless18.era <- PCV_eras
norless18.rankfreqplot <- rankfreqplot(df = norless18.rank, country = "Norway", agegroup = "< 18 yrs", PCV_eras = norless18.era, positionheight = 5)
norless18.alleras <- ggarrange(norless18.rankfreqplot[[1]] + coord_cartesian(xlim = c(0,60), ylim = c(0,60)), 
                               norless18.rankfreqplot[[2]] + labs(title = norless18.era[2]) + coord_cartesian(xlim = c(0,50), ylim = c(0,50)),
                               common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1) # pdf dim 3x12 in

nor_more18 <- norway_morethaneq18 %>% filter(Serotype != "Empty") %>% filter(Serotype != "NT")  %>% filter(Serotype != "ROUGH")
normore18.plot <- serorank(nor_more18, country = "Norway")
normore18.uniquesero <- unique(normore18.plot$data[,4])
normore18.rank <- agegrp_vcn
normore18.era <- PCV_eras
normore18.rankfreqplot <- rankfreqplot(df = normore18.rank, country = "Norway", agegroup = "> 18 yrs", PCV_eras = normore18.era, positionheight = 50)
normore18.alleras <- ggarrange(normore18.rankfreqplot[[1]] + coord_cartesian(xlim = c(0,50), ylim = c(0,500)), 
                               normore18.rankfreqplot[[2]] + labs(title = normore18.era[2]) + coord_cartesian(xlim = c(0,50), ylim = c(0,500)),
                               common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1) # pdf dim 3x12 in

norranks <- grid.arrange(norless18.plot, normore18.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in
norrankfreqgrob <- ggarrange(norless18.rankfreqplot[[1]] + labs(x = "", y = "") +
                               coord_cartesian(xlim = c(0,50), ylim = c(0,50)), 
          norless18.rankfreqplot[[2]] + labs(title = norless18.era[2],x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,50)),
          normore18.rankfreqplot[[1]] + labs(x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,500)), 
          normore18.rankfreqplot[[2]] + labs(title = normore18.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,50), ylim = c(0,500)),
          common.legend = TRUE, legend = "bottom", ncol =2, nrow = 2) # pdf dim 6x10 in
annotate_figure(norrankfreqgrob, 
                left = text_grob("IPD cases", rot = 90),
                right = text_grob("Cumulative frequency", rot = 270))#,
                #fig.lab = "D", fig.lab.face = "bold") # pdf dim 4.5 x 9


# usa
usaless18.plot <- serorank(usa_lessthan18, country = "USA")
usaless18.uniquesero <- unique(usaless18.plot$data[,4])
usaless18.rank <- agegrp_vcn
usaless18.era <- PCV_eras
usaless18.rankfreqplot <- rankfreqplot(df = usaless18.rank, country = "USA", agegroup = "< 18 yrs", PCV_eras = usaless18.era, positionheight = 100)
usaless18.alleras <- ggarrange(usaless18.rankfreqplot[[1]]+coord_cartesian(xlim = c(0,25), ylim = c(0,1150)),
                               usaless18.rankfreqplot[[2]]+ labs(title = usaless18.era[2]) + coord_cartesian(xlim = c(0,25), ylim = c(0,1150)), 
                               common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)# pdf dim 3x12 in

usamore18.plot <- serorank(usa_morethaneq18, country = "USA")
usamore18.uniquesero <- unique(usamore18.plot$data[,4])
usamore18.rank <- agegrp_vcn
usamore18.era <- PCV_eras
usamore18.rankfreqplot <- rankfreqplot(df = usamore18.rank, country = "USA", agegroup = "> 18 yrs", PCV_eras = usamore18.era, positionheight = 100)
usamore18.alleras <- ggarrange(usamore18.rankfreqplot[[1]]+coord_cartesian(xlim = c(0,25), ylim = c(0,2700)), 
                               usamore18.rankfreqplot[[2]]+ labs(title = usamore18.era[2]) + coord_cartesian(xlim = c(0,25), ylim = c(0,2700)), 
                               common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)# pdf dim 3x12 in

usaranks <- grid.arrange(usaless18.plot, usamore18.plot, ncol = 2, nrow = 1) # pdf dim 5x7 in
usarankfreqgrob <- ggarrange(usaless18.rankfreqplot[[1]]+ labs(x = "", y = "") + 
                               coord_cartesian(xlim = c(0,25), ylim = c(0,1150)),
          usaless18.rankfreqplot[[2]]+ labs(title = usaless18.era[2],x = "", y = "") + 
            coord_cartesian(xlim = c(0,25), ylim = c(0,1150)), 
          usamore18.rankfreqplot[[1]]+ labs(x = "", y = "") + coord_cartesian(xlim = c(0,25), ylim = c(0,2700)), 
          usamore18.rankfreqplot[[2]]+ labs(title = usamore18.era[2], x = "", y = "") + 
            coord_cartesian(xlim = c(0,25), ylim = c(0,2700)), 
          common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 2) # pdf dim 6x12 in
annotate_figure(usarankfreqgrob, 
                left = text_grob("IPD cases", rot = 90),
                right = text_grob("Cumulative frequency", rot = 270))#,
                #fig.lab = "E", fig.lab.face = "bold") # pdf dim 4.5 x 9

grid.arrange(ausranks, finranks, franranks, norranks, usaranks, #labels = c("A","B", "C", "D", "E"), 
             ncol = 2, nrow = 3) # pdf dim 

ggarrange(ausranks, finranks, franranks, norranks, usaranks, labels = c("A","B", "C", "D", "E"), 
             ncol = 2, nrow = 3)

major.sero <- c(levels(ausless18.uniquesero), levels(ausmore18.uniquesero), levels(finless18.uniquesero), levels(finmore18.uniquesero), 
                levels(franceless16.uniquesero), levels(francemore16.uniquesero), #levels(italyinf.uniquesero), 
                levels(norless18.uniquesero), levels(normore18.uniquesero), levels(usaless18.uniquesero), levels(usamore18.uniquesero))
major.sero <- unique(major.sero)
major.sero <- major.sero[!major.sero %in% c("Pre-PCV", "Pre-PCV13", "Pre-PCV10", "Post-PCV")]

##########################################################################################################################################################
########### FREQUENCY OVER TIME PLOTS ####################################################################################################################
tot.freq$country <- factor(tot.freq$country)

plot.freq <- function(df) {
  ggplot(df) + 
    geom_line(aes(x = Year, y = freq, colour = country, group = country)) + 
    geom_point(aes(x = Year, y = freq, colour = country, shape = pcv_era), size = 2) +
    scale_shape_manual(values = c(1, 15, 2)) +
    scale_colour_discrete(drop = FALSE) + 
    labs(y = "Percent of total IPD", title = paste(df$agecat, ": Serotype ", df$Serotype, sep = ""), colour = "Country", shape = "PCV era") +
    theme_bw()
}

pdf(file = "freqOTplots_children.pdf", width = 6, height = 4)
lapply(major.sero, function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% children)))
dev.off() 

pdf(file = "freqOTplots_adults.pdf", width = 6, height = 4)
lapply(major.sero, function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% adults)))
dev.off() 

VT_freqOT_kidz <- lapply(major.sero[major.sero %in% c(VT7, VT10, VT13)], function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% children)) 
                          + theme(legend.position = "none"))
VT_freqOT_kidz_leg <- g_legend(VT_freqOT_kidz[[1]] + theme(legend.position = "bottom"))
grid.arrange(VT_freqOT_kidz[[8]], VT_freqOT_kidz[[10]], VT_freqOT_kidz[[12]], VT_freqOT_kidz[[2]], VT_freqOT_kidz[[3]], VT_freqOT_kidz[[5]], VT_freqOT_kidz[[6]],
             VT_freqOT_kidz[[1]], VT_freqOT_kidz[[11]], VT_freqOT_kidz[[7]], VT_freqOT_kidz[[9]], VT_freqOT_kidz[[4]],  
             bottom = VT_freqOT_kidz_leg, nrow = 3) # pdf dim 8 x 14

NVT_freqOT_kidz <- lapply(major.sero[!major.sero %in% c(VT7, VT10, VT13)], function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% children)) 
                          + theme(legend.position = "none"))
NVT_freqOT_kidz_leg <- g_legend(NVT_freqOT_kidz[[2]] + theme(legend.position = "bottom"))
grid.arrange(NVT_freqOT_kidz[[1]], NVT_freqOT_kidz[[2]], NVT_freqOT_kidz[[3]], NVT_freqOT_kidz[[4]], NVT_freqOT_kidz[[5]], NVT_freqOT_kidz[[6]], NVT_freqOT_kidz[[7]],
             NVT_freqOT_kidz[[8]], NVT_freqOT_kidz[[9]], NVT_freqOT_kidz[[10]], NVT_freqOT_kidz[[11]], NVT_freqOT_kidz[[12]], NVT_freqOT_kidz[[13]], NVT_freqOT_kidz[[14]],
             NVT_freqOT_kidz[[15]], NVT_freqOT_kidz[[16]], NVT_freqOT_kidz[[17]], bottom = NVT_freqOT_kidz_leg, nrow = 4, ncol = 5) #pdf dim 8 x 14

VT_freqOT_aduz <- lapply(major.sero[major.sero %in% c(VT7, VT10, VT13)], function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% adults))
                           + theme(legend.position = "none"))
VT_freqOT_aduz_leg <- g_legend(VT_freqOT_aduz[[1]] + theme(legend.position = "bottom"))
grid.arrange(VT_freqOT_aduz[[8]], VT_freqOT_aduz[[10]], VT_freqOT_aduz[[12]], VT_freqOT_aduz[[2]], VT_freqOT_aduz[[3]], VT_freqOT_aduz[[5]], VT_freqOT_aduz[[6]],
             VT_freqOT_aduz[[1]], VT_freqOT_aduz[[11]], VT_freqOT_aduz[[7]], VT_freqOT_aduz[[9]], VT_freqOT_aduz[[4]],  
             bottom = VT_freqOT_aduz_leg, nrow = 3) # pdf dim 8 x 14

NVT_freqOT_aduz <- lapply(major.sero[!major.sero %in% c(VT7, VT10, VT13)], function(x) plot.freq(tot.freq %>% filter(Serotype == x) %>% filter(Age.group %in% adults)) 
                          + theme(legend.position = "none"))
NVT_freqOT_aduz_leg <- g_legend(NVT_freqOT_aduz[[6]] + theme(legend.position = "bottom"))
grid.arrange(NVT_freqOT_aduz[[1]], NVT_freqOT_aduz[[2]], NVT_freqOT_aduz[[3]], NVT_freqOT_aduz[[4]], NVT_freqOT_aduz[[5]], NVT_freqOT_aduz[[6]], NVT_freqOT_aduz[[7]],
             NVT_freqOT_aduz[[8]], NVT_freqOT_aduz[[9]], NVT_freqOT_aduz[[10]], NVT_freqOT_aduz[[11]], NVT_freqOT_aduz[[12]], NVT_freqOT_aduz[[13]], NVT_freqOT_aduz[[14]],
             NVT_freqOT_aduz[[15]], NVT_freqOT_aduz[[16]], NVT_freqOT_aduz[[17]], bottom = NVT_freqOT_aduz_leg, nrow = 4, ncol = 5) #pdf dim 8 x 14

##########################################################################################################################################################
########### INCIDENCE OVER TIME PLOTS ####################################################################################################################
melt.tot$country <- factor(melt.tot$country)

plot.inc <- function(df){
  ggplot(df) + 
    geom_line(aes(x = Year, y = (dis.cases/population)*100000, colour = country, group = country)) +
    geom_point(aes(x = Year, y = (dis.cases/population)*100000, colour = country, shape = pcv_era), size = 2) +
    scale_shape_manual(values = c(1, 15, 2)) + 
    scale_colour_discrete(drop = FALSE) + 
    labs(y = "Incidence per 100,000 ppl", colour = "Country", shape = "PCV era", title = paste(df$agecat, ": Serotype ", df$Serotype, sep = "")) +
    theme_bw()
}

countryplot.inc <- function(df){
  ggplot(df) + 
    geom_line(aes(x = Year, y = (dis.cases/population)*100000, colour = Serotype, group = Serotype)) +
    geom_point(aes(x = Year, y = (dis.cases/population)*100000, colour = Serotype, shape = pcv_era), size = 2) +
    scale_shape_manual(values = c(1, 15, 2)) + 
    labs(y = "Incidence per 100,000 ppl", colour = "Serotype", shape = "PCV era", title = paste(df$country, df$Serotype, sep = "")) +
    theme_bw()
}

pdf(file = "incOTplots_children.pdf", width = 6, height = 4)
lapply(major.sero, function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% children)))
dev.off() 

pdf(file = "incOTplots_adults.pdf", width = 6, height = 4)
lapply(major.sero, function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% adults)))
dev.off() 

VT_incOT_kidz <- lapply(major.sero[major.sero %in% c(VT7, VT10, VT13)], function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% children)) 
                         + theme(legend.position = "none"))
VT_incOT_kidz_leg <- g_legend(VT_incOT_kidz[[1]] + theme(legend.position = "bottom"))
grid.arrange(VT_incOT_kidz[[8]], VT_incOT_kidz[[10]], VT_incOT_kidz[[12]], VT_incOT_kidz[[2]], VT_incOT_kidz[[3]], VT_incOT_kidz[[5]], VT_incOT_kidz[[6]],
             VT_incOT_kidz[[1]], VT_incOT_kidz[[11]], VT_incOT_kidz[[7]], VT_incOT_kidz[[9]], VT_incOT_kidz[[4]],  
             bottom = VT_incOT_kidz_leg, nrow = 3) # pdf dim 8 x 14

NVT_incOT_kidz <- lapply(major.sero[!major.sero %in% c(VT7, VT10, VT13)], function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% children)) 
                          + theme(legend.position = "none"))
NVT_incOT_kidz_leg <- g_legend(NVT_incOT_kidz[[1]] + theme(legend.position = "bottom"))
grid.arrange(NVT_incOT_kidz[[1]], NVT_incOT_kidz[[2]], NVT_incOT_kidz[[3]], NVT_incOT_kidz[[4]], NVT_incOT_kidz[[5]], NVT_incOT_kidz[[6]], NVT_incOT_kidz[[7]],
             NVT_incOT_kidz[[8]], NVT_incOT_kidz[[9]], NVT_incOT_kidz[[10]], NVT_incOT_kidz[[11]], NVT_incOT_kidz[[12]], NVT_incOT_kidz[[13]], NVT_incOT_kidz[[14]],
             NVT_incOT_kidz[[15]], NVT_incOT_kidz[[16]], NVT_incOT_kidz[[17]], bottom = NVT_incOT_kidz_leg, nrow = 4, ncol = 5) #pdf dim 8 x 14

VT_incOT_aduz <- lapply(major.sero[major.sero %in% c(VT7, VT10, VT13)], function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% adults))
                         + theme(legend.position = "none"))
VT_incOT_aduz_leg <- g_legend(VT_incOT_aduz[[1]] + theme(legend.position = "bottom"))
grid.arrange(VT_incOT_aduz[[8]], VT_incOT_aduz[[10]], VT_incOT_aduz[[12]], VT_incOT_aduz[[2]], VT_incOT_aduz[[3]], VT_incOT_aduz[[5]], VT_incOT_aduz[[6]],
             VT_incOT_aduz[[1]], VT_incOT_aduz[[11]], VT_incOT_aduz[[7]], VT_incOT_aduz[[9]], VT_incOT_aduz[[4]],  
             bottom = VT_incOT_aduz_leg, nrow = 3) # pdf dim 8 x 14

NVT_incOT_aduz <- lapply(major.sero[!major.sero %in% c(VT7, VT10, VT13)], function(x) plot.inc(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% adults)) 
                          + theme(legend.position = "none"))
NVT_incOT_aduz_leg <- g_legend(NVT_incOT_aduz[[1]] + theme(legend.position = "bottom"))
grid.arrange(NVT_incOT_aduz[[1]], NVT_incOT_aduz[[2]], NVT_incOT_aduz[[3]], NVT_incOT_aduz[[4]], NVT_incOT_aduz[[5]], NVT_incOT_aduz[[6]], NVT_incOT_aduz[[7]],
             NVT_incOT_aduz[[8]], NVT_incOT_aduz[[9]], NVT_incOT_aduz[[10]], NVT_incOT_aduz[[11]], NVT_incOT_aduz[[12]], NVT_incOT_aduz[[13]], NVT_incOT_aduz[[14]],
             NVT_incOT_aduz[[15]], NVT_incOT_aduz[[16]], NVT_incOT_aduz[[17]], bottom = NVT_incOT_aduz_leg, nrow = 4, ncol = 5) #pdf dim 8 x 14

##########################################################################################################################################################
########### SEROTYPE INCIDENCE GROWTH/DECLINE RATES ######################################################################################################
inc.rates <- function(df, agegrp) { # function that estimates crude rate of growth/decline in incidence for each PCV era
  df.cases <- df %>% filter(agecat %in% agegrp)
  df.tibs <- df.cases %>% group_by(Serotype, Age.group, agecat, country, pcv_era) %>% nest()
  #df.tibs.mod <- df.tibs %>% mutate(model = purrr::map(data, ~glm(dis.cases ~ Year + offset(log(population/100000)), family = poisson, data = .)))
  df.tibs.mod <- df.tibs %>% mutate(model = purrr::map(data, ~lm((dis.cases/population)*100000 ~ Year, data = .)))
  #df.tibs.summ <- df.tibs.mod %>% unnest(model %>% map(broom::glance))
  #df.tibs.coeff <- data.frame(df.tibs.summ %>% unnest(model %>% purrr::map(broom::tidy)))
  #bh <- p.adjust(df.tibs.coeff$p.value, method="BH") # Benjamin & Hochberg method - method to control false discovery rate
  
  df.tibs.summ <- bind_rows(df.tibs.mod$model %>% map(broom::glance))
  df.tibs.coeff <- bind_rows(df.tibs.mod$model %>% map(broom::tidy))
    
  rates <- df.tibs.coeff %>% filter(!grepl('(Intercept)', term))
  
  bh <- p.adjust(rates$p.value, method="BH") # Benjamin & Hochberg method - method to control false discovery rate
  rates$p.value <- bh
  #mean(bh, na.rm = TRUE) # 0.2797367
  rates <- bind_cols(df.tibs, rates)
  return(rates)
}

rateplotz <- function(countrydf) { ## function to plot each country's serotype rates before and after vaccination ##
  ggplot(countrydf, aes(x = country, y = estimate, fill = pcv_era)) + 
    geom_bar(stat = "identity", position = position_dodge(0.5), width = 0.5) + theme_minimal() + 
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
    geom_errorbar(aes(ymin = countrydf$estimate - countrydf$std.error, ymax = countrydf$estimate + countrydf$std.error),
                  position = position_dodge(.5), width = 0.2) +
    labs(fill = "PCV era", x = "Country", y = "Incidence growth rate", title = countrydf$Serotype[1]) +
    theme_bw() + theme_light()
}

country.rateplotz <- function(countrydf) { ## function to plot each country's serotype rates before and after vaccination ##
  countrydf$cat[countrydf$Serotype %in% VT7] <- "VT7"
  countrydf$cat[countrydf$Serotype %in% VT10] <- "VT10"
  countrydf$cat[countrydf$Serotype %in% VT13] <- "VT13"
  countrydf$cat[!(countrydf$Serotype %in% c(VT7, VT10, VT13))] <- "NVT"
  countrydf$cat <- factor(countrydf$cat, levels = c("VT7", "VT10", "VT13", "NVT"))
  
  countrydf$cols[countrydf$Serotype %in% VT7] <- "#35B779FF" #"#70AD47"
  countrydf$cols[countrydf$Serotype %in% VT10] <- "#31688EFF" #"#4472C4"
  countrydf$cols[countrydf$Serotype %in% VT13] <- "#E69F00" #"#ED7D31"
  countrydf$cols[!(countrydf$Serotype %in% c(VT7, VT10, VT13))] <- "#440154FF"#"black"
  
  countrydf.nest <- countrydf %>% group_by(Serotype, cols) %>% nest()
  countrydf.nest <- countrydf.nest[,-3]
  countrydf.nest <- countrydf.nest[order(countrydf.nest$Serotype),]
  
  ggplot(countrydf, aes(x = Serotype, y = estimate, colour = pcv_era)) + 
    geom_bar(stat = "identity", fill = "white", position = position_dodge(width = 0.7), width = 0.7) + 
    scale_colour_manual(values = c("Pre-PCV" = "#999999", "Pre-PCV10/13" = "#E69F00", "Post-PCV" = "#56B4E9")) +
    geom_errorbar(aes(ymin = countrydf$estimate - countrydf$std.error, ymax = countrydf$estimate + countrydf$std.error),
                  position = position_dodge(0.7), width = 0.2) +
    labs(colour = "PCV era", x = "Serotype", y = "Incidence growth rate", title = countrydf$country[1]) +
    theme(axis.text.x = element_text(colour = countrydf.nest$cols), 
          panel.background = element_rect(fill = "white", colour = "grey"), panel.grid = element_line(colour = "grey92"))
}

#### rates for children 
sero.rates.children <- lapply(major.sero, function(x) inc.rates(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% children), "Children"))
serorate.plots.children <- lapply(sero.rates.children, function(x) rateplotz(x) + theme(legend.position = "none"))
sero.legend.children <- g_legend(rateplotz(sero.rates.children[[1]])+ theme(legend.position = "bottom"))

VT.rates.children <- grid.arrange(serorate.plots.children[[13]], serorate.plots.children[[15]], serorate.plots.children[[18]], serorate.plots.children[[3]], 
                         serorate.plots.children[[5]], serorate.plots.children[[7]], serorate.plots.children[[10]], serorate.plots.children[[1]], 
                         serorate.plots.children[[17]], serorate.plots.children[[11]], serorate.plots.children[[14]], serorate.plots.children[[6]], 
                         bottom = sero.legend.children, nrow = 3, ncol = 4) #pdf dim 8 x 15

sero.legend.children <- g_legend(rateplotz(sero.rates.children[[1]]))
NVT.rates.children <- grid.arrange(serorate.plots.children[[2]], serorate.plots.children[[4]], serorate.plots.children[[8]], serorate.plots.children[[9]], 
                                   serorate.plots.children[[12]], serorate.plots.children[[16]], serorate.plots.children[[19]], serorate.plots.children[[20]], 
                                   serorate.plots.children[[21]], serorate.plots.children[[22]], serorate.plots.children[[23]], serorate.plots.children[[24]], 
                                   serorate.plots.children[[25]], serorate.plots.children[[26]], serorate.plots.children[[27]], serorate.plots.children[[28]], 
                                   serorate.plots.children[[29]], sero.legend.children, nrow = 3, ncol = 6) #pdf dim 8 x 18

tot.rates.children <- bind_rows(sero.rates.children)

country.rates.children <- lapply(unique(tot.rates.children$country), function(x) country.rateplotz(tot.rates.children %>% filter(country == x)))

pdf(file = "inc_growthrates_children.pdf", width = 15, height = 4)
country.rates.children
dev.off()


#### rates for adults
sero.rates.adults <- lapply(major.sero, function(x) inc.rates(melt.tot %>% filter(Serotype == x) %>% filter(Age.group %in% adults), "Adults"))
serorate.plots.adults <- lapply(sero.rates.adults, function(x) rateplotz(x) + theme(legend.position = "none"))
sero.legend.adults <- g_legend(rateplotz(sero.rates.adults[[1]])+ theme(legend.position = "bottom"))

VT.rates.adults <- grid.arrange(serorate.plots.adults[[13]], serorate.plots.adults[[15]], serorate.plots.adults[[18]], serorate.plots.adults[[3]], 
                                  serorate.plots.adults[[5]], serorate.plots.adults[[7]], serorate.plots.adults[[10]], serorate.plots.adults[[1]], 
                                  serorate.plots.adults[[17]], serorate.plots.adults[[11]], serorate.plots.adults[[14]], serorate.plots.adults[[6]], 
                                  bottom = sero.legend.adults, nrow = 3, ncol = 4) #pdf dim 8 x 15

sero.legend.adults <- g_legend(rateplotz(sero.rates.adults[[1]]))
NVT.rates.adults <- grid.arrange(serorate.plots.adults[[2]], serorate.plots.adults[[4]], serorate.plots.adults[[8]], serorate.plots.adults[[9]], 
                                   serorate.plots.adults[[12]], serorate.plots.adults[[16]], serorate.plots.adults[[19]], serorate.plots.adults[[20]], 
                                   serorate.plots.adults[[21]], serorate.plots.adults[[22]], serorate.plots.adults[[23]], serorate.plots.adults[[24]], 
                                   serorate.plots.adults[[25]], serorate.plots.adults[[26]], serorate.plots.adults[[27]], serorate.plots.adults[[28]], 
                                   serorate.plots.adults[[29]], sero.legend.adults, nrow = 3, ncol = 6) #pdf dim 8 x 18

tot.rates.adults <- bind_rows(sero.rates.adults)

country.rates.adults <- lapply(unique(tot.rates.adults$country), function(x) country.rateplotz(tot.rates.adults %>% filter(country == x)))

pdf(file = "inc_growthrates_adults.pdf", width = 15, height = 4)
country.rates.adults
dev.off()

totrates <- bind_rows(tot.rates.children, tot.rates.adults)

aus.rateplots <- ggarrange(country.rates.children[[1]]+labs(x = NULL, title = "Australia Children"), 
                              country.rates.adults[[1]] + labs(title = "Australia Adults"), nrow = 2,
                              common.legend = TRUE, legend = "right") # pdf dim 6 x 12

fin.rateplots <- ggarrange(country.rates.children[[2]]+labs(x = NULL, title = "Finland Children"), 
                              country.rates.adults[[2]] + labs(title = "Finland Adults"), nrow = 2,
                              common.legend = TRUE, legend = "right")

fra.rateplots <- ggarrange(country.rates.children[[3]]+labs(x = NULL, title = "France Children"),
                           country.rates.adults[[3]] + labs(title = "France Adults"), nrow = 2,
                           common.legend = TRUE, legend = "right")

nor.rateplots <- ggarrange(country.rates.children[[4]]+labs(x = NULL, title = "Norway Children"), 
                              country.rates.adults[[4]] + labs(title = "Norway Adults"), nrow = 2,
                           common.legend = TRUE, legend = "right")

usa.rateplots <- ggarrange(country.rates.children[[5]]+labs(x = NULL, title = "USA Children"), 
                              country.rates.adults[[5]] + labs(title = "USA Adults"), nrow = 2,
                           common.legend = TRUE, legend = "right")

#ita.rateplots <- country.rates.children[[4]] + labs(title = "Italy Children")
