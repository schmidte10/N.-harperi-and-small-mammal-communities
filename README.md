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

Before we can start plotting we need to make sure we are plotting with the correct values the code below take you through the steps on how to extract the proper values to create the dataframe necessary to begin plotting. Note: 'comm.comp1' NEEDS TO BE LOADED!!!! (or whatever you called your brms model). 
```
var <- get_variables(comm.comp1)                                                ## character vector of the names of the variables in a variety of fitted Bayesian model types
vari <- get_variables(comm.comp1)[c(1,4:9,12:23)]                               ## select only the variables needed for figure [trapline and trapline:species interactions]

comm.comp1_int_draws2 <- comm.comp1 %>% spread_draws(!!!syms(vari))             ## extract posterior draws for the variables selected above ^

comm.comp1_draws <- comm.comp1 %>% gather_draws(b_Intercept, b_Species1, b_Species3) %>%  ## get the posterior draws for the intercept (reference value [RBV]), and other values                                                                                           ## (in this case other species [1 = deer mouse]; [2 = woodland jumping mouse]
  left_join(comm.comp1_int_draws2, by = c(".chain",".iteration",".draw")) %>%             ## join collected draws from above ^ with posterior draws located in                                                                                                             ##  comm.comp1_int_draws2
  
  mutate(Trapline4_mean = case_when(`.variable` == 'b_Intercept' ~                        ### to get the mean probability of infection for each species on each trapline it is 
                                      `.value`,                                           ### necessary to add values from the interecept (RBV on trapline 2) to those 
                                    `.variable` != 'b_Intercept' ~                        ### pertaining to other factors in order to get the correct values. See below for an
                                      `.value` + b_Intercept),                            ### example
         
         Trapline1_mean = case_when(`.variable` == 'b_Intercept' ~                        ### For example to get the correct values for species 1 on trapleine 1 we would need 
                                      `.value` + b_Trapline1,                             ### to added the intercept (b_intercept) to the value related to species 1 [DM;.value]
                                    `.variable` == 'b_Species1' ~                         ### to the value related to trapline 1 (b_Trapline1), to the interaction between 
                                      `.value` + b_Intercept +                            ### species 1:trapline1 (b_species1:Trapline1). This process needs to be repeated for
                                      b_Trapline1 + `b_Species1:Trapline1`,               ### each speecies in each trapline. 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept +                            ### Remember that species 2 [RBV] is the reference factor. 
                                      b_Trapline1 + `b_Species3:Trapline1`), 
         
         Trapline2_mean = case_when(`.variable` == 'b_Intercept' ~ 
                                      `.value` + b_Trapline2, 
                                    `.variable` == 'b_Species1' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline2 + `b_Species1:Trapline2`, 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline2 + `b_Species3:Trapline2`), 
         
         Trapline3_mean = case_when(`.variable` == 'b_Intercept' ~ 
                                      `.value` + b_Trapline3, 
                                    `.variable` == 'b_Species1' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline3 + `b_Species1:Trapline3`, 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline3 + `b_Species3:Trapline3`), 
         
         Trapline5_mean = case_when(`.variable` == 'b_Intercept' ~ 
                                      `.value` + b_Trapline5, 
                                    `.variable` == 'b_Species1' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline5 + `b_Species1:Trapline5`, 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline5 + `b_Species3:Trapline5`), 
         
         Trapline7_mean = case_when(`.variable` == 'b_Intercept' ~ 
                                      `.value` + b_Trapline7, 
                                    `.variable` == 'b_Species1' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline7 + `b_Species1:Trapline7`, 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline7 + `b_Species3:Trapline7`), 
         
         Trapline8_mean = case_when(`.variable` == 'b_Intercept' ~ 
                                      `.value` + b_Trapline8, 
                                    `.variable` == 'b_Species1' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline8 + `b_Species1:Trapline8`, 
                                    `.variable` == 'b_Species3' ~ 
                                      `.value` + b_Intercept + 
                                      b_Trapline8 + `b_Species3:Trapline8`), 
         
         transformed_trapline4 = exp(Trapline4_mean)/(1+exp(Trapline4_mean)),           ## Bernoulli distribution applied a logit transformation on values, therefore 
         transformed_trapline1 = exp(Trapline1_mean)/(1+exp(Trapline1_mean)),           ## to get the true values we need to reverse transform the values we calculated
         transformed_trapline2 = exp(Trapline2_mean)/(1+exp(Trapline2_mean)), 
         transformed_trapline3 = exp(Trapline3_mean)/(1+exp(Trapline3_mean)), 
         transformed_trapline5 = exp(Trapline5_mean)/(1+exp(Trapline5_mean)), 
         transformed_trapline7 = exp(Trapline7_mean)/(1+exp(Trapline7_mean)), 
         transformed_trapline8 = exp(Trapline8_mean)/(1+exp(Trapline8_mean))) 
```
## Almost there!
```
comm.comp1_draws_plotting <-  comm.comp1_draws %>%                                      
  pivot_longer(cols = starts_with("transformed"),                                       ## pivot transformed values within the data frame
               names_to = "Trapline",                                                   ## change column name to Trapline (values in column should be 'trapline1;trapline2...'
               values_to = "infect_prob") %>%                                           ## change 'value' column name to 'infect_prob'...or whatever you want to call it
  transmute(species = case_when(`.variable` == "b_Intercept" ~ "RBV",                   ## changing variable names 
                                `.variable` == "b_Species1" ~ "DM", 
                                `.variable` == "b_Species3" ~ "WJM"), 
            trapline = case_when(Trapline == 'transformed_trapline4' ~ "trapline4",    ## changing variable names in Trapline column 
                                 Trapline == 'transformed_trapline1' ~ "trapline1", 
                                 Trapline == 'transformed_trapline2' ~ "trapline2",
                                 Trapline == 'transformed_trapline3' ~ "trapline3",
                                 Trapline == 'transformed_trapline5' ~ "trapline5",
                                 Trapline == 'transformed_trapline7' ~ "trapline7",
                                 Trapline == 'transformed_trapline8' ~ "trapline8"), 
            infect_prob = infect_prob,                                                 ## making new columns
            speciesb = species, 
            traplineb = trapline, 
            chain = `.chain`, 
            iteration = `.iteration`, 
            draw_num = `.draw`) %>% 
  unite("id",speciesb,traplineb,sep = "_") %>%                                         ## merging columns 'speciesb', and 'traplineb' to creat new variable names stored in 
                                                                                       ## column called 'id'
 mutate(id = factor(id, levels = c("RBV_trapline4","DM_trapline4","WJM_trapline4",     ## setting factor levels
                                    "RBV_trapline1","DM_trapline1","WJM_trapline1", 
                                    "RBV_trapline2","DM_trapline2","WJM_trapline2", 
                                    "RBV_trapline3","DM_trapline3","WJM_trapline3",
                                    "RBV_trapline5","DM_trapline5","WJM_trapline5",
                                    "RBV_trapline7","DM_trapline7","WJM_trapline7",
                                    "RBV_trapline8","DM_trapline8","WJM_trapline8")))
 ``` 
You have now created your dataframe that wil be used for plotting!! Congratulations this has been a lot of work. Perhaps time to take a bit of a break, pat your self on the back, and grab a drink. Next step is plotting, ready whenever you are. Quick note the dataframe 'comm.comp1_draws_plotting' will also be used to creat figures 1,3, & 4; Tables 1 & 2. 

## Figure 1 [place holder for actually figure number in the manuscript] 
Plotting will be down in ggplot2. For the first plot we will make use of the package 'ggridges' to get a nice look at our posterior distributions. The package 'ggridges' should be loaded from the beginning. 'ggridges' is a really useful plot for examining posterior draws. 
### Lets go!
 ```  
 ridge_plot2 <- comm.comp1_draws_plotting %>%                       
  ggplot(aes(infect_prob, trapline, fill = species)) +                                              ## x-axis = infect_prob; y-axis= trapline
  scale_fill_manual(values = alpha(c("steelblue1","orange1",  "seagreen3"),0.7),                    ## sitting colors for different species, probably want some transparency for                                                                                                     ## overlapping ridges
                    labels=c("Deer mouse", "Red-backed voles","Woodland jumping mouse")) +          ## labels for species names
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +                                  ## scale determines how much overlap will occur between ridges;                                                                                                                 ## rel_min_height refers to the tail cut-off, adjust as necessary
  theme_bw()+                                                                                       ## choosen theme
  scale_x_continuous(breaks = seq(0, 1, by = .2))+                                                  ## x-axis numerical limit and break
  xlab("P(Infection)")+ylab("Habitat type")+                                                        ## axis labels
  scale_y_discrete(labels = c("Sugar maple hardwood", " Cut-over mixed-wood", "Dense mixed-wood",   ## y-axis tick labels 
                              "Conifer", "White pine/white spruce", "Black spruce/aspen", 
                              "White/red pine"))+
  theme(axis.text=element_text(size=18),                                                            ## size of habitat names (y-axis); numbers (x-axis)
        axis.title=element_text(size=25,face="bold"),                                               ## size of axis label names
        legend.position = "bottom",                                                                 ## position of legend
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),                     ## gridlines (or lack thereof) 
        legend.text=element_text(size=20))+                                                         ## legend text size 
  guides(fill=guide_legend(title=""));ridge_plot2                                                   ## legend title (or lack thereof) 

pdf(file = "Figure2_ridge_plot2.pdf", height = 14, width=14)                                        ## save file as .pdf
ridge_plot2                                                                                         ## plot
dev.off()
``` 
### Congrats! 

## Figure 2 [place holder for actually figure number in the manuscript] 
Once again before starting make sure that your model is loaded
``` 
prp <- orangemitedata_v3 %>%                                                                ## dataframe that the model was run on 
  group_by(Species) %>%                                                                     ## group by 'species' variable
  data_grid(perc_rbv=seq(from = min(orangemitedata_v3$perc_rbv),                            ## make a dataframe that contains the lowest and highester values for x-axis value
                         to=max(orangemitedata_v3$perc_rbv),                                ## dataframe should contain 101 values (per group - species were grouped in the
                         length.out = 101))%>%                                              ## previous section of code
  add_fitted_draws(rbv_pop_model, n=100,                                                    ## Add 100 observations from fitted posterior draws (obtained from the brms model)                                                                                               ## for each perc_rbv value created in the 
                   re_formula = NA)%>%                                                      ## dateframe above
  ggplot(aes(x=perc_rbv, color=Species))+                                                   ## begin plotting
  geom_line(aes(y=.value, group = paste(Species, .draw)), alpha=0.2)+                       ## Draw lines based on the posterior draws for each species 
  geom_smooth(aes(y=.value, color=Species), se=F, size=1.5)+theme_classic()+                ## Draw a thicker line that shows the means (of draw) value for each species
  xlab("% red backed vole in community composition")+ylab("P(Infection)")+                  ## axis labels
  scale_color_manual(labels=c("Red backed vole", "Deer mouse", "Woodland jumping mouse"),   ## manuallt changing color labels
                     values=c("orange1", "steelblue1", "seagreen3"))+
  theme(legend.position = "bottom",                                                         ## legend position
        legend.title = element_blank(),                                                     ## legend title (none)
        axis.title.x = element_text(size=25),                                               ## x-axis label size
        axis.text.x = element_text(size=18),                                                ## x-axis text size                                                
        axis.title.y = element_text(size=25),                                               ## y-axis label size 
        axis.text.y = element_text(size=18),                                                ## y-axis text size 
        legend.text = element_text(size=15),                                                ## legend text size 
        axis.line = element_line(size=1.2))+                                                ## axis line size
  geom_vline(xintercept = 0.2, color="black", alpha=0.4, size=1); prp                       ## vertical line place on x-axis 

pdf(file = "perc_rbv_pop_model.pdf", height = 10, width=14)                                 ## save file as .pdf
prp                                                                                         ## plot 
dev.off()
``` 
## Figure 3 and Table 1 [place holder for actually figure and table number in the manuscript]
``` 
p_rbv <- comm.comp1_draws_plotting[comm.comp1_draws_plotting$species=="RBV",]     ## subset posterior draws that are red backed voles [RBV]
p_dm <- comm.comp1_draws_plotting[comm.comp1_draws_plotting$species=="DM",]       ## subset posterior draws that are deer mice [DM] 
p_wjm <- comm.comp1_draws_plotting[comm.comp1_draws_plotting$species=="WJM",]     ## subset posterior draws that are woodland jumpin mice [wjm]
p <- left_join(p_rbv, p_dm, by=c("trapline","draw_num")) %>%                      ## left join p_rbv, p_dm, p_wjm; simply reorganizing the dataframe
  left_join(p_wjm, by=c("trapline", "draw_num")) %>%                              
  mutate( 
    rbv.dm = infect_prob.x - infect_prob.y,              ## substract infection probability of red backed voles (infect_prob.x) from that of deer mice (infect_prob.y) 
    rbv.wjm = infect_prob.x - infect_prob,               ## substract infection probability of red backed voles (infect_prob.x) from that of woodland jumping mice (infect_prob)
    dm.wjm = infect_prob.y - infect_prob)                ## substract infection probability of deer mice (infect_prob.y) from that of woodland jumping mice (infect_prob)
p_rbvdm <- p[,c(3,21)] %>%                               ## subset values from the rbv.dm comparison and column that lists the traplines
  mutate(comparison = "rbv.dm") %>%                      ## create a new column that lists the species comparison that is beinmade (this will be useful later on) 
  rename(infect_prob = rbv.dm)                           ## rename the rbv.dm column 'infect_prob' - now repear the process for the different comparisons made. 
p_rbvwjm <- p[,c(3,22)] %>% 
  mutate(comparison = "rbv.wjm") %>% 
  rename(infect_prob = rbv.wjm) 
p_dmwjm <- p[,c(3,23)] %>% 
  mutate(comparison = "dm.wjm") %>% 
  rename(infect_prob = dm.wjm) 

## Code below is used to obtain the median as well as lower and upper HPD interval for differences in infection probability between different species

db_rbvdm <- as.data.frame(t(aggregate(p_rbvdm$infect_prob ~ p_rbvdm$trapline, FUN = function(i)c(median = median(i), HPDinterval = hdi(i)))))%>% 
  row_to_names(row_number = 1)
db_rbvwjm <- as.data.frame(t(aggregate(p_rbvwjm$infect_prob ~ p_rbvwjm$trapline, FUN = function(i)c(median = median(i), HPDinterval = hdi(i)))))%>% 
  row_to_names(row_number = 1)
db_dmwjm <- as.data.frame(t(aggregate(p_dmwjm$infect_prob ~ p_dmwjm$trapline, FUN = function(i)c(median = median(i), HPDinterval = hdi(i)))))%>% 
  row_to_names(row_number = 1)
table1 <- rbind(db_rbvdm, db_rbvwjm, db_dmwjm)                ## combined all the values into one table
write.csv(table1, "table1_commcomp1.csv", row.names = T)      ## convert the table to a .csv file and save it to the working directory


p2 <- bind_rows(p_rbvdm,p_rbvwjm,p_dmwjm) %>%                                                   ## combine p_rbvdm, p_rbvwjm, and p_dmwjm 
  mutate(trapline_proper = case_when((trapline == "trapline1") ~ "Sugar maple hardwood",        ## Make a new column that will change trapline names so there correspond the 
                                     (trapline == "trapline2") ~ "Cut-over mixed-wood",         ## appropriate habitat type       
                                     (trapline == "trapline3") ~ "Dense mixed-wood",
                                     (trapline == "trapline4") ~ "Conifer",
                                     (trapline == "trapline5") ~ "White pine/white spruce",
                                     (trapline == "trapline7") ~ "Black spruce/aspen",
                                     TRUE ~ "White/red pine")) %>% 
  mutate(trapline_proper = factor(trapline_proper, levels = c("Sugar maple hardwood", "Cut-over mixed-wood", "Dense mixed-wood", #convert column to factor
                                                              "Conifer", "White pine/white spruce", "Black spruce/aspen", 
                                                              "White/red pine"))) %>% 
  mutate(comparison = factor(comparison, levels = c("dm.wjm","rbv.wjm","rbv.dm")))                                               # convert column to factor


sp_comp <- ggplot(p2, aes(infect_prob, comparison, fill=comparison))+                   ## begin plotting figure
  geom_vline(xintercept=0, lty=2, color="grey28")+                                      ## create verticle line on the x-axis at 0
  stat_halfeye(interval_colour="red", point_colour="darkred", point_fill="red")+        ## show posterior distribution using stat_halfeye
  facet_wrap(~trapline_proper, scales = "free_y")+theme_classic()+                      ## facet wrap figure based on habitat type (`trapline proper`)
  scale_fill_manual(values = c('grey10', 'grey10', 'grey10'),guide=F)+                  ## Manually fill in colors; no legend
  xlab("P(Infection)")+ylab("Comparisons")+                                             ## axis labels
  theme(axis.title.x = element_text(size=25),                                           ## x-axis label title size
        axis.text.x = element_text(size=18),                                            ## x-axis text size
        axis.title.y = element_text(size=25),                                           ## y-axis label title size
        axis.text.y = element_text(size=18),                                            ## y-axis text size
        strip.text = element_text(size=20));sp_comp                                     ## facet plot title size
  
pdf(file = "species_comparisons.pdf", height = 10, width=14)                            ## save as .pdf
sp_comp                                                                                 ## plot
dev.off()
```
## Figure 4 and Table 2 [place holder for actually figure and table number in the manuscript]
Remember you will need to have created the dataframe comm.comp1_draws_plotting (see above) to create this figure
```
p <- comm.comp1_draws_plotting  %>% 
  group_by(trapline) %>%                                                  ## group by trapline
  group_split()%>%                                                        ## split into different dataframes based in trapline (habitat)
  reduce(left_join, by=c("species","chain","iteration","draw_num")) %>%   ## combines dataframes (simply to left_join); reorganization of dataframes
  ungroup() %>%                                                           ## ungroup dataframe
  group_by(species) %>%                                                   ## group by species
  group_split()                                                           ## split dataframe apart based on species

hd <- list()                                                              ## create an empty list called 'hd'
habitat <- list()                                                         ## create an empty list called 'habitat'
for (i in 1:3) {                                                          ## run a for loop that compares infection probability between habitat types for each species
  habitat[[i]] <- p[[i]] %>% 
    mutate(`4-1` = infect_prob.y.y - infect_prob.x, 
           `4-2` = infect_prob.y.y - infect_prob.y, 
           `4-3` = infect_prob.y.y - infect_prob.x.x, 
           `4-5` = infect_prob.y.y - infect_prob.x.x.x, 
           `4-7` = infect_prob.y.y - infect_prob.y.y.y, 
           `4-8` = infect_prob.y.y - infect_prob, 
           `1-2` = infect_prob.x - infect_prob.y, 
           `1-3` = infect_prob.x - infect_prob.x.x,
           `1-5` = infect_prob.x - infect_prob.x.x.x, 
           `1-7` = infect_prob.x - infect_prob.y.y.y, 
           `1-8` = infect_prob.x - infect_prob, 
           `2-3` = infect_prob.y - infect_prob.x.x, 
           `2-5` = infect_prob.y - infect_prob.x.x.x, 
           `2-7` = infect_prob.y - infect_prob.y.y.y, 
           `2-8` = infect_prob.y - infect_prob, 
           `3-5` = infect_prob.x.x - infect_prob.x.x.x,
           `3-7` = infect_prob.x.x - infect_prob.y.y.y, 
           `3-8` = infect_prob.x.x - infect_prob, 
           `5-7` = infect_prob.x.x.x - infect_prob.y.y.y, 
           `5-8` = infect_prob.x.x.x - infect_prob, 
           `7-8` = infect_prob.y.y.y - infect_prob) %>% 
    select(c(2,6:8,33:53))                                                ## select columns for 'species', 'chain', 'iteration', 'draw_num', and all habitat comparisons
  
  hd[[i]] <- habitat[[i]][5:25] %>%                                       ## place habitat comparisons within a new list 'hd', only include habitat comparisons
    summarise_each(list(mean = mean, hdi = hdi)) %>%                      ## summarise comparisons by getting the mean and HDI (lower and upper)
    as.data.frame()%>% t() %>%                                            ## create a dataframe for habitat comaprisons for each species; transpose
    as.data.frame()                                                       ## make dataframe
}
hd[[1]] <- rename(hd[[1]], DM.lower = "V1", DM.upper = "V2")              ## rename dataframe columns as DM.lower and DM.upper (repeat for other dataframes/species
hd[[2]] <- rename(hd[[2]], RBV.lower = "V1", RBV.upper = "V2")
hd[[3]] <- rename(hd[[3]], WJM.lower = "V1", WJM.upper = "V2")
hrdw <- cbind(hd[[1]],hd[[2]],hd[[3]])                                    ## combined columns 
                                                                          ## note that in the table create the rows that contain the mean will be the same for .lower and .upper                                                                           ## columns, however the HDI lower and upper values shoud be different 
write.csv(hrdw, "table2_commcomp1.csv", row.names = T)                    ## write table as a .csv value; will be placed in your working directory 

p2 <- list()                                                              ## create list called p2
for (i in 1:3) {                                                          ## create 'for' loop
  p2[[i]] <-  habitat[[i]][5:25] %>% gather()%>%                          ## gether posterior draws from 'habitat' list for all habitat (trapline) comparisons
    mutate(species = paste(unique(habitat[[i]][1]))) %>%                  ## create a new column that contains the species name
    rename(comparison = "key",                                            ## rename columns
           infect_prob = "value", 
           species = "species")
}



p2_plot <- rbind(p2[[1]],p2[[2]],p2[[3]]) %>%                             ## rejoin species data
  as.data.frame()                                                         ## make a dataframe


habitat_comp <- ggplot(p2_plot, aes(infect_prob, comparison))+                              ## begin plotting
  geom_vline(xintercept=0, lty=2, color="grey28")+                                          ## create a vertical line that will at value 0 on the x-axis
  stat_pointinterval(interval_colour="red", point_colour="darkred", point_fill="red")+      ## use stat_pointinterval for plotting (too many comparisons for stat_halfeye)
  facet_wrap(~species)+theme_classic()+                                                     ## facet wrap plots
  scale_fill_manual(values = c('grey10', 'grey10', 'grey10'),guide=F)+                      ## fill in colors (choosing them to be all the same) 
  xlab("P(Infection)")+ylab("Comparisons")+                                                 ## axis labels
  theme(axis.title.x = element_text(size=25),                                               ## x-axis label title size
        axis.text.x = element_text(size=18),                                                ## x-axis text size
        axis.title.y = element_text(size=25),                                               ## y-axis label title size
        axis.text.y = element_text(size=18),                                                ## y-axis text size
        strip.text = element_text(size=20),                                                 ## facet plot titles size
        panel.spacing = unit(2, "lines")); habitat_comp                                     ## spacing between panels

pdf(file = "habitat_comp.pdf", height = 10, width=14)                                       ## save as .pdf
habitat_comp                                                                                ## plot
dev.off()
``` 
All of these figures and tables were created for community compositons 1 & 2 (high and low abundance red backed vole years). Different dataframes were used for each year however the process for making the figures remained the same. Please feel free to contact me using the contact information on my GitHub profile is you have any problems, questions, or concerns regardin the data manipulation/organisation or creation of figures. Hope you were able to find the code above helpful. 
