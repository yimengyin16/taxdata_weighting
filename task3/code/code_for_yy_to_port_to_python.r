
#****************************************************************************************************
#                globals - change as needed ####
#****************************************************************************************************
# target_state <- "NY"
#dbox <- "C:/Users/donbo/Dropbox/"
yy_dir <- "task3/"  # paste0(dbox, "state_puf_shared/yy/")

includes_dir <- paste0(yy_dir, "code/")
input_data_dir <- paste0(yy_dir, "input_files/")
log_dir <- paste0(yy_dir, "logfiles/") # for ipopt logfiles
interim_results_dir <- paste0(yy_dir, "interim_results/")
outfiles_dir <- paste0(yy_dir, "output_files/")


#****************************************************************************************************
#                libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library("knitr")

library("ipoptr")


#****************************************************************************************************
#                includes ####
#****************************************************************************************************
source(paste0(includes_dir, "functions_for_yy_porting.r"))


#****************************************************************************************************
#                get data ####
#****************************************************************************************************
targets_state <- read_csv(paste0(input_data_dir, "targets_state.csv"))
targets_state # note that target_num uniquely identifies each target - unique combination of stabbr, constraint_name, constraint_type

# get national puf using Boyd 2017 grow factors, with selected 2017 Tax-Calculator output variables, and wtus_2017 weight
pufbase_all <- read_csv(paste0(input_data_dir, "puf2017_weighted.csv"))

ratio_adjusted_weights_state <- read_csv(paste0(input_data_dir, "ratio_adjusted_weights_state.csv"))


#****************************************************************************************************
#                prepare PUF base file - with initial state weight and essential variables ####
#****************************************************************************************************
pufbase_state <- pufbase_all %>%
  select(RECID, MARS, one_of(unique(targets_state$puf_link))) %>% # harmless warning given for MARS1..4, which we will create in a minute
  left_join(ratio_adjusted_weights_state %>% 
              select(RECID, AGI_STUB, stabbr, weight_state),
            by="RECID") %>%
  mutate(weight_initial=weight_state) %>% # we need to keep a copy of the state weight
  # create a column for each MARS category that has the initial state weight if the record is of that category
  # so that we have columns that match up against targets for number of returns by MARS category
  pivot_wider(names_from = MARS, names_prefix = "MARS", values_from = weight_state, values_fill = list(weight_state=0))
glimpse(pufbase_state)
ns(pufbase_state)
saveRDS(pufbase_state, paste0(interim_results_dir, "pufbase_state.rds"))


#****************************************************************************************************
#                prepare nonzero constraint coefficients ####
#****************************************************************************************************
# make a long sparse file that has the initial weight and the variable value for each numeric variable
# we will make constraint coefficients from this
long <- pufbase_state %>%
  pivot_longer(cols=-c(RECID, stabbr, AGI_STUB, weight_initial), names_to="puf_link", values_to = "value") %>%
  filter(value!=0) # if the variable value is 0 then the constraint coefficient will always be 0 and we don't need it
ht(long)
count(long, AGI_STUB)
count(long, puf_link)

# full join against the targets file allowing multiple targets (typically 2) per variable -- # of nonzeros, and amount
long_full <- long %>% 
  full_join(targets_state, by=c("stabbr", "AGI_STUB", "puf_link")) %>%  # we don't need to keep target but it is convenient to have
  filter(!is.na(RECID)) # RECID will be NA in a full join if we have a constraint for a category that has no records in the long file
# example:
#   we won't have any returns where variable c05800 tax liability is nonzero and the AGI_STUB==1 (AGI < $1)
#   even though the historical table 2 data have about 910 such returns, with total taxbc of ~$31m
#   because all returns in this group in the puf have liability 0 and we dropped them from the long file

# keep track of the good constraints - those that are found in long_full
good_con <- unique(long_full$target_num)
targets_state %>% filter(target_num %in% good_con) # view full details on good constraints
targets_state %>% filter(!target_num %in% good_con) # view full details on bad constraints
# not surprisingly, we lost a few constraints in the lowest AGI ranges, where no PUF records have nonzero values for these

targets_state_good <- targets_state %>%
  filter(target_num %in% good_con)
# this is what we will use for targeting and for constraint coefficients

count(targets_state_good, constraint_type) # the types of constraint coefficients we will need

# now we can create a data frame of nonzero constraint coefficients
nzcc <- targets_state_good %>%
  left_join(long_full %>% select(-c(puf_link, year, stabbr, AGI_STUB, table_desc, target,
                                    constraint_name, constraint_type, constraint_base_name)), 
            by="target_num") %>%
  mutate(nzcc=case_when(constraint_type=="amount" ~ value * weight_initial,
                        constraint_type=="n_nonzero" ~ weight_initial,
                        constraint_type=="n_exempt" ~ value * weight_initial, # number of exemptions times record weight
                        constraint_type=="n_returns" ~ weight_initial,
                        TRUE ~ NA_real_))
glimpse(nzcc) # we've kept a lot of unnecessary variables but some are nice to have
ht(nzcc)


#****************************************************************************************************
#                set tolerances so that we can calculate constraint bounds ####
#****************************************************************************************************
# compute the starting point (the value on the file using weight_initial) for each target and use it to set constraint bounds
starting_point <- nzcc %>%
  group_by(AGI_STUB, target_num, constraint_name, constraint_type, table_desc) %>%
  summarise(target=first(target), # all records within a target group will have the same value
            file=sum(nzcc),
            diff=file - target,
            pdiff=diff / target * 100) %>% # pdiff will help us decide on tolerances
  ungroup

starting_point %>%
  filter(AGI_STUB==2) %>%
  # filter(str_detect(table_desc, coll("unemp", ignore_case = TRUE))) %>%
  mutate(table_desc=str_remove(table_desc, "Number of") %>% str_sub(., 1, 35)) %>% # shorten the name for better readability
  # note that we do not have target A05800 as a target in AGI_STUB 1 but we do in other stubs
  kable(digits=c(rep(0, 8), 1), format="rst", format.args = list(big.mark=","))


# create a few priority levels
priority1 <- c("A00100", "A00200", "A05800", "A09600", "A18500", "N2", "MARS1", "MARS2", "MARS4")
priority2 <- c("N00100", "N00200", "N05800", "N09600", "N18500")
tolerances <- starting_point %>%
  mutate(tol_default=case_when(constraint_name %in% c(priority1, priority2) ~ .005,
                               TRUE ~ abs(pdiff/100) * .10))

tolerances %>%
  filter(AGI_STUB==2) %>%
  select(-c(target_num, file)) %>%
  mutate(tol_default=tol_default * 100) %>%
  mutate(table_desc=str_remove(table_desc, "Number of") %>% str_sub(., 1, 35)) %>%
  kable(digits=c(rep(0, 6), 1, 1), format = "rst", format.args = list(big.mark = ","))


#****************************************************************************************************
#                run the optimization on one or more agi stubs ####
#****************************************************************************************************
stubs <- 1:10 # to run all stubs
# stubs <- 2 # for test run, to compare to boyd_stub02.out
xdf <- ldply(stubs, runstub, tolerances, nzcc, pufbase_state, log_dir, interim_results_dir, .progress="text")
glimpse(xdf)
count(xdf, AGI_STUB)
ht(xdf)
quantile(xdf$x, probs=c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))
quantile(xdf$x * xdf$wt_init, probs=c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))

saveRDS(xdf, paste0(outfiles_dir, "x_state.rds"))


#****************************************************************************************************
#                create the final file ####
#****************************************************************************************************
pufbase_state <- readRDS(paste0(interim_results_dir, "pufbase_state.rds"))
xdf <- readRDS(xdf, file=paste0(outfiles_dir, "x_state.rds"))

puf_state <- pufbase_state %>%
  left_join(xdf %>% select(RECID, x), by="RECID") %>%
  mutate(weight_state=weight_initial * x)

sum(puf_state$weight_initial)
sum(puf_state$weight_state)
sum(puf_state$weight_state * puf_state$E00200) / 1e9 # wages in billions
# $550.8716 billion is the number that I (Boyd) get, which is within 0.2% of the Historical Table 2 target

# save file if desired



