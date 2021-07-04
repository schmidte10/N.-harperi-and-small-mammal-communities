# N.-harperi-and-small-mammal-communities
## set working directory 
```
setwd()
```
## Necessary packages that need to be loaded
```
#list.files(path = ".")
library(readr)
library(readxl) 
library(grid)
library(lubridate)
library(janitor)
library(brms)
library(plyr)
library(dplyr)
library(rstan)
library(reshape2)
library(modelr)
library(emmeans)
library(bayestestR) 
library(tidyverse) 
library(tidybayes) 
library(patchwork) 
library(ggpubr)
library(cowplot)
library(HDInterval)
library(ggridges)
```
## Importing data 
```
orangemitedata_2020 <- read_delim("orangemitedata_2020.txt", 
                                  "\t", escape_double = FALSE, col_types = cols(Date = col_datetime(format = "%d/%m/%Y"), 
                                   Trapline = col_factor(levels = c("1","2", "3", "4", "5", "7", "8", "9")), 
                                   Species = col_factor(levels = c("1", "2", "3")), 
                                   Sex = col_factor(levels = c("1","2")), 
                                   Year = col_factor(levels = c("1","2", "3", "4"))), trim_ws = TRUE)

smallmammaldata3 <- read_delim("smallmammaldata_v3.txt", 
                               "\t", 
                               escape_double = FALSE, 
                               col_types = cols(Date = col_date(format = "%Y-%m-%d"), 
                               Trapline = col_factor(levels = c("1","2", "3", "4", "5", "6", "7","8", "9")), 
                               Species = col_factor(levels = c("1","2", "3")), Sex = col_factor(levels = c("1", "2")),
                               Year = col_factor(levels = c("1","2", "3")), 
                               X14 = col_skip(), X15 = col_skip()), 
                               trim_ws = TRUE)
abund <- read_excel("smallmammaldata2016-2018_all.xlsx") #import dataframes
abund_2020 <- read_delim("smallmammaldata2020_all.txt", 
                         "\t", 
                         escape_double = FALSE, 
                         col_types = cols(Date = col_date(format = "%d/%m/%Y")), 
                         trim_ws = TRUE)                       
colnames(smallmammaldata3)[10] <- 'OrangeMites'
colnames(orangemitedata_2020)[10] <- 'OrangeMites'
``` 
## Data cleaning and organising 
```
orangemitedata_v2 <- dplyr::bind_rows(smallmammaldata3, orangemitedata_2020)           ##joining dataframes
orangemitedata_v2$date.julian <- as.POSIXct(orangemitedata_v2$Date, format="%Y-%m-%d") ##convert to POSIXct date
orangemitedata_v2$date.julian <- format(orangemitedata_v2$date.julian, "%j")           ##convert to julian date
orangemitedata_v2$date.julian <- as.numeric(orangemitedata_v2$date.julian)             ##convert to numeric
orangemitedata_v2$Species  <- relevel(orangemitedata_v2$Species , ref = "2")           ##change reference species to '2' (red backed voles)
orangemitedata_v2$Trapline  <- factor(orangemitedata_v2$Trapline, levels = c('4','1','2','3','5','7','8','9')) 
orangemitedata_v2$Trapline  <- relevel(orangemitedata_v2$Trapline , ref = "4")         ##change reference habitat to '4' (conifer)
orangemitedata_v2$OrangeMites <- as.factor(orangemitedata_v2$OrangeMites)              ##Change orange mite variable to factor
orangemitedata_v2$date.julian <- scale(orangemitedata_v2$date.julian, center = T, scale = T) ##center and scale julian date
orangemitedata_v2$Year <- as.factor(orangemitedata_v2$Year)                            ##change variable 'Year' to factor
orangemitedata_v2$TagLeft <- ifelse(is.na(orangemitedata_v2$TagLeft), orangemitedata_v2$TagRight, orangemitedata_v2$TagLeft) 
                                                                                       ## ^ If there is no right ear tag replace value with same value as left ear tag
orangemitedata_v2 <- orangemitedata_v2[!with(orangemitedata_v2, is.na(TagLeft)&is.na(TagRight)),] ##remove samples missing ear tags
names(orangemitedata_v2[c(3,4,8,13,14)])                                               ##names of important variables to be included in model


abund_1 <- rbind(abund,abund_2020) %>%                                                 ##join data frames 
  clean_names()%>%                                                                     ##make sure that column names are formatted 
  mutate(orange_mites_y_n = toupper(orange_mites_y_n),                                 ##make 'Y/N' in orange_mites_y_n column uppercase
         species = toupper(species),                                                   ##make species codes uppercase
         orange_mites_y_n = replace_na(orange_mites_y_n, "N"),                         ##replace NA's with a 'N'                                      
         year = year(date),                                                            ##make separate column that contains only the year          
         month = month(date),                                                          ##make separate colum that contains only the month 
         tag_left = coalesce(tag_left,tag_right),                                      ##if left ear tage is missing replace with right ear tag
         tag_left = ifelse(tag_left=="RIP", tag_right, tag_left)) %>%                  ## if left ear tag says "RIP" replace with right ear tag 
  group_by(year)%>% 
  mutate(first_om_appearence = min(date[orange_mites_y_n=="Y"]))%>%                    ##determine the day mites appeared first for each year
  ungroup() %>% 
  subset(date > first_om_appearence)%>%                                                ##remove observations that occurred before mites first appeared for each year
  mutate(month_name = month.abb[month],                                                ##column that has abbreviated month names - will be used later 
                             species = case_when((species == "SM") ~ "DM",             ##there is a species listed as 'SM' in the dataframe, since we do not have a species                                                                                            ##designated to the code 'sm' I assume the observer meant 'DM'. D and S are beside each                                                                                          ##other on the keybroad
                             (species == "NFLYER") ~ "FLYER",                          ##change species called 'NFLYER' to "FLYER'
                             TRUE ~ as.character(.$species)),                          ##remaining species retain name in currently in column
         trapline_group = as.numeric(trapline),                                        ##make trapline values numeric
         trapline_group = case_when((trapline >100 & trapline<200) ~ 100,              ##make a separate group that groups traplines based on habitat type
                                    (trapline >200 & trapline<300) ~ 200,
                                    (trapline >300 & trapline<400) ~ 300,
                                    (trapline >400 & trapline<500) ~ 400,
                                    (trapline >500 & trapline<600) ~ 500,
                                    (trapline >700 & trapline<800) ~ 700,
                                    (trapline >800 & trapline<900) ~ 800,
                                    TRUE ~ 900)) %>% 
  relocate(trapline_group, .after=trapline) %>%                                        ##place new column next to the 'trapline' column to make sure it worked
  unite("group_ID", year, month_name, trapline_group, remove = F)                      ##create groups based on 'year','month', and 'trapline_group' of each observation


abund_2 <-  abund_1 %>%                                                                
  group_by(group_ID) %>%                                                               
  distinct(tag_left, .keep_all = T) %>%                                                ##only keep unique left tag observations
  ungroup() %>%
  reshape2::dcast(group_ID~species, length) %>%                                        ##reorganise dataframe
  select(-c(`TRAP MALFUNCTION`, `TRAP DISTURBANCE`)) %>%                               ##remove oberservations recorded as 'TRAP MALFUNCTION' or 'TRAP DISTURBANCE' 
  mutate(total = rowSums(.[2:8]),                                                      ##get the sum of rows for rows 2:8 
         perc_rbv = RBV/total)                                                         ##divide the number of red backed voles by the total number of animals to get the                                                                                                ##percentage of red backed voles (perc_rbv) within the small mammal community
         
              
                                                     


temp1 <- separate(abund_2, group_ID, c("year","month","trapline.no"))%>%               ##separate group_ID into 'year','month', and 'trapline number'
  filter(year=="2016"|year=="2017")                                                    ##subset observations from years 2016 and 2017
summarise_all(temp1[4:14],funs(sum))/sum(temp1[14])                                    ##summary of community composition 

temp2 <- separate(abund_2, group_ID, c("year","month","trapline.no"))%>%               ##see description above ^
  filter(year=="2018"|year=="2020")
summarise_all(temp2[4:14],funs(sum))/sum(temp2[14])


##### need to combine the abund_2 dataframe with previous data frame 
#### so that the model OrangeMites ~ species*rbc_perc can be run 

orangemitedata_v3 <- orangemitedata_v2 %>% 
  mutate(year = year(Date), 
         month = month(Date),
         month_name = month.abb[month], 
         trapline_group= as.numeric(as.character(Trapline))*100) %>% 
  unite("group_ID", year,month_name,trapline_group, remove=F) %>%
  #join 'abund_2' dataframe to 'orangemitedata_v3' dataframe by the group_ID column 
  #so that the percentage of voles within a community (per month) can be recorded for each observation 
  left_join(abund_2[,c(1,13)], by="group_ID")
  ``` 
  # Community composition 1 (High red backed vole year) model 
  ### Priors
  ``` 
p4 <- get_prior(OrangeMites~Species+Trapline+Year+date.julian+Species:Trapline+Species:date.julian+(1|TagLeft), 
                family=bernoulli(), 
                data=orangemitedata_v2[!orangemitedata_v2$Trapline==9,][c(3,4,8,10,13,14)] %>% 
                  subset(Year == 1 | Year == 2));p4
p4 <- c(set_prior("normal(0,100)", class="b", coef="Species1"), 
        set_prior("normal(0,100)", class="b", coef="Species3"), 
        set_prior("normal(0,100)", class="b", coef="Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Trapline8"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline8"),
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline8"), 
        set_prior("normal(0,100)", class="b", coef="date.julian"),
        set_prior("normal(0,100)", class="b", coef="Species1:date.julian"),
        set_prior("normal(0,100)", class="b", coef="Species3:date.julian"), 
        set_prior("normal(0,100)", class="b", coef="Year2"));p4
``` 
### Model 1 (High red backed vole year)
``` 
comm.comp1 <- brm(OrangeMites~Species+Trapline+Year+date.julian+Species:Trapline+Species:date.julian+(1|TagLeft), 
                  data=na.omit(orangemitedata_v2[!orangemitedata_v2$Trapline==9,][c(3,4,8,10,13,14)]) %>% 
                    subset(Year == 1 | Year == 2), 
                  family=bernoulli(),
                  prior = p4,
                  warmup = 4000,
                  iter =8000, 
                  seed=123,
                  thin=2,
                  chains = 4,
                  cores = 2, 
                  save_all_pars=T, 
                  control = list(adapt_delta=0.99, max_treedepth=15))
saveRDS(comm.comp1, "comm_comp1.RDS")                                               ## save model 
comm.comp1<- readRDS("comm_comp1.RDS")                                              ## load model for future use 
``` 
# Community composition 2 (Low red backed vole years) 
### Priors 
``` 
p4 <- c(set_prior("normal(0,100)", class="b", coef="Species1"),
        set_prior("normal(0,100)", class="b", coef="Species3"),
        set_prior("normal(0,100)", class="b", coef="Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Trapline8"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Species1:Trapline8"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline1"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline2"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline3"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline5"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline7"), 
        set_prior("normal(0,100)", class="b", coef="Species3:Trapline8"), 
        set_prior("normal(0,100)", class="b", coef="date.julian"),
        set_prior("normal(0,100)", class="b", coef="Species1:date.julian"),
        set_prior("normal(0,100)", class="b", coef="Species3:date.julian"), 
        set_prior("normal(0,100)", class="b", coef="Year4"));p4
 ``` 
 ### Model 2 (Low red backed vole years) 
 ```
 comm.comp2 <- brm(OrangeMites~Species+Trapline+Year+date.julian+Species:Trapline+Species:date.julian, 
                  data=orangemitedata_v2[!orangemitedata_v2$Trapline==9,][c(3,4,8,10,13,14)] %>% 
                    subset(Year == 3 | Year == 4), 
                  family=bernoulli(),
                  prior = p4,
                  warmup = 4000,
                  iter =8000, 
                  thin=2,
                  seed=123,
                  chains = 4,
                  cores = 2, 
                  save_all_pars=T, 
                  control = list(adapt_delta=0.99, max_treedepth=15))
saveRDS(comm.comp2, "comm_comp2.RDS")                                               ## save model
comm.comp2<- readRDS("comm_comp2.RDS")                                              ## load model (more future use) 
``` 
# Influence of red backed vole abundance on infection probability 
### Priors 
``` 
p4 <- get_prior(OrangeMites~perc_rbv+Species*perc_rbv+(1|TagLeft), 
                family=bernoulli(), 
                data=orangemitedata_v3[!orangemitedata_v3$Trapline==9,]);p4
#%>% 
#subset(Trapline == 2|Trapline == 4|Trapline == 5|Trapline == 6|Trapline == 7)%>%
#group_by(TagLeft) %>% 
#top_n(1, date.julian))
p4 <- c(set_prior("normal(0,100)", class="b", coef="Species1"),
        set_prior("normal(0,100)", class="b", coef="Species3"),
        set_prior("normal(0,100)", class="b", coef="perc_rbv"), 
        set_prior("normal(0,100)", class="b", coef="perc_rbv:Species1"), 
        set_prior("normal(0,100)", class="b", coef="perc_rbv:Species3"));p4
``` 
### Model 
```
#%>% subset(Species == 1|Species==2|Species)
rbv_pop_model <- brm(OrangeMites~perc_rbv+Species*perc_rbv+(1|TagLeft), 
                     data=orangemitedata_v3[!orangemitedata_v3$Trapline==9,], 
                     family=bernoulli(),
                     prior = p4,
                     warmup = 4000,
                     iter =8000, 
                     thin=2,
                     seed=123,
                     chains = 4,
                     cores = 2, 
                     save_all_pars=T, 
                     control = list(adapt_delta=0.99, max_treedepth=15))
saveRDS(rbv_pop_model, "rbv_pop_model.RDS")                                       ## save model
rbv_pop_model <- readRDS("rbv_pop_model.RDS")                                     ## load model (for future use)
``` 
# Figures 
## Figure 1 

 
