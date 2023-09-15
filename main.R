### Title: Functional Limitations in Stroke Survivors: Pre-Stroke Physical Activity Matters
### Authors: Dan Orsholits, Matthieu Boisgontier

### Packages
## To import the data
library("haven")
## Use tidyverse stuff
library("dplyr")
library("tidyr")
library("magrittr")
## For matching
library("MatchIt")
## To estimate models
library("lme4")
library("lmerTest")
## For figures
library("ggplot2")
library("ggeffects")

### Import data

## We are not interested in the retrospective waves
waves <- c(1:2, 4:8)

## The modules we want
modules <- c(
  physical_health = "ph", cover_screen = "cv_r", behavioural_risk = "br",
  demographics = "dn", education = "gv_isced", health_generated = "gv_health",
  household_income = "hh"
)

data_list <- vector(mode = "list", length = length(waves))
names(data_list) <- paste0("wave_", waves)

## Set the directory where the data is located
data_dir <- "."
## It is assumed that the directory structure from SHARE's original zip files
## is preserved

## Import SHARE master file
master_file <- read_stata(file = paste0(data_dir, "/sharewX_rel8-0-0_gv_allwaves_cv_r.dta"))

for (wave in waves) {
    files <- list.files(
        path = data_dir, pattern = paste0("^sharew", wave), full.names = TRUE,
        recursive = TRUE
    )
    cat("Currently on wave ", wave, "\n", sep = "")

    ## Convert directly from tibble to data.table
    for (module in modules) {
        data_list[[paste0("wave_", wave)]][[module]] <-
            read_stata(file = files[grepl(paste0(
                                      ".*_", module, "\\.dta"
                                  ), files)])
        cat("Finished reading module: \"", module, "\"\n", sep = "")
    }
}

### Convert wide format master to long
ivw_ages <- master_file %>%
    select(matches("^mergeid$|^country$|^gender$|age_int_w[[:digit:]]$|^interview_w[[:digit:]]$"))
ivw_ages_long <- ivw_ages %>%
    pivot_longer(cols = !c("mergeid", "country", "gender"),
                 names_sep = "_w",
                 names_to = c(".value", "wave")) %>%
    mutate(wave = as.numeric(wave)) %>%
    arrange(mergeid, wave) %>%
    group_by(mergeid) %>%
    mutate(firstivw =
               ifelse(is.na(match(1, interview)),
                      NA,
                      nth(wave, match(1, interview)))) %>%
    ungroup()

### A simple function to be used with mapply to add a wave variable to
### tibbles
add_wave_mapply <- function(tib, wave) tib %>% mutate(wave = wave)

### Get the variables we need

## Variables that should always be selected
always_keep <- c("mergeid", "wave")

## These come from the ph modules
stroke_vars <- c(
    "ph067_2", "ph071_2", "ph072_2",
    "ph006d4", "ph073_2", "ph074_2",
    "ph075_2"
)

## Other raw health variables
raw_health_vars <- c(
   "ph049",
   "ph004",
   "ph006"
)

## These come from the br modules
phys_act_vars <- c(vig_act = "br015_", mod_act = "br016_")

## These come from the gv_health modules
adl_vars <- c("adl", "iadl", "bmi", "bmi2")
chronic_cond <- c("chronic", "chronic2")

## These come from the gv_isced modules
educ_vars <- c("isced1997_r")

## This function is a wrapper to simplify extracting columns/variables from a
## list of data.frames or objects that inherit from data.frame (data.table,
## tibble, etc.). It uses "grepl" meaning that if there are no matches, an empty
## object will be returned. Also this function will do partial matching and
## *not* exact matching.
select_vars_from_data <- function(data_list, module_name, variables) {
  ret <- lapply(data_list, "[[", module_name)
  grep_vars <- paste0("^", variables, collapse = "|")
  ret <- lapply(ret, function(x) {
    sub <- x[, grepl(grep_vars, names(x)), with = FALSE]

    return(sub)
  })
  return(ret)
}

### Select the relevant variables from the modules

## Physical health
ph_select <-
    select_vars_from_data(data_list, "ph",
                          c(always_keep, stroke_vars, raw_health_vars))
ph_select <-
    mapply(FUN = add_wave_mapply, ph_select,
           c(1, 2, 4:8), SIMPLIFY = FALSE)

ph_select <- lapply(ph_select, function(x) {
  ret <- x
  names(ret) <- gsub("ph071_2", "n_strokes_last_int", names(ret), fixed = TRUE)
  names(ret) <- gsub("^ph067_2$|ph072_2", "stroke_last_int", names(ret))
  return(ret)
})

ph_select <- bind_rows(ph_select)

## Generated health variables
gv_health_select <-
    select_vars_from_data(data_list,
                          "gv_health", c(always_keep, adl_vars, chronic_cond))

gv_health_select <-
  lapply(gv_health_select,
         function(x) {
           names(x) <- gsub("chronic2.*", "chronic2", names(x))
           names(x) <- gsub("chronicw.*", "chronic", names(x))
          return(x)
         })

gv_health_select <- mapply(
  FUN = add_wave_mapply, gv_health_select, c(1, 2, 4:8),
  SIMPLIFY = FALSE
)
gv_health_select <- bind_rows(gv_health_select)

## Physical activity
phys_select <-
    select_vars_from_data(data_list, "br", c(always_keep, phys_act_vars))
phys_select <-
    mapply(FUN = add_wave_mapply, phys_select, c(1, 2, 4:8), SIMPLIFY = FALSE)
phys_select <- bind_rows(phys_select)

## Education
educ_select <-
    select_vars_from_data(data_list, "gv_isced", c(always_keep, educ_vars))
educ_select <-
    mapply(FUN = add_wave_mapply, educ_select, c(1, 2, 4:8), SIMPLIFY = FALSE)
educ_select <- bind_rows(educ_select)

## Merge together all the modules using Reduce. We use the full_join
## function from dplyr.
reduc_merge <- dplyr::full_join
model_data <- Reduce(reduc_merge, list(
  ivw_ages_long, ph_select, gv_health_select,
  phys_select, educ_select
))

## Clean up education (set missing values) and find highest attained level for
## each respondent
model_data <- model_data %>%
    mutate(isced1997_r = if_else(isced1997_r < 0 & !is.na(isced1997_r),
                                 haven::labelled(NA_real_),
                                 isced1997_r)) %>%
    mutate(isced1997_r = if_else(isced1997_r %in% c(95, 97),
                                 haven::labelled(NA_real_),
                                 isced1997_r)) %>%
    group_by(mergeid) %>%
    mutate(max_edu = if_else(all(is.na(isced1997_r)),
                             haven::labelled(NA_real_),
                             max(isced1997_r, na.rm = TRUE))) %>%
    ungroup()

## Recode education into 3 groups
# Primary or less
# Secondary
# Tertiary
model_data <- model_data %>%
  mutate(max_edu = factor(max_edu, levels = 0:6,
         labels = c(rep("Primary or less", 2),
                    rep("Secondary", 3),
                    rep("Tertiary", 2)))) %>%
  mutate(max_edu = relevel(max_edu, "Secondary"))

## Clean up adl & iadl below 0 recoded as missing
model_data <- model_data %>%
    mutate(adl = ifelse(adl < 0, haven::labelled(NA_real_), adl)) %>%
    mutate(iadl = ifelse(iadl < 0, haven::labelled(NA_real_), iadl))

## Use raw ADL/IADL info
adl_cols <- paste0("ph049d", 1:6)
## Dummies 14 & 15 introduced from wave 6
iadl_cols <- paste0("ph049d", 7:13)
model_data <- model_data %>%
  mutate(across(starts_with("ph049"),
                ~ if_else(.x < 0, haven::labelled(NA_real_), .x))) %>%
  mutate(adl_raw = rowSums(.[adl_cols])) %>%
  mutate(iadl_raw = rowSums(.[iadl_cols]))

## See if they had a stroke before the first interview
model_data <- model_data %>%
    group_by(mergeid) %>%
    mutate(stroke_before_study = if_else(wave == firstivw, ph006d4,
                                         haven::labelled(NA_real_))) %>%
    mutate(stroke_before_study = if_else(stroke_before_study < 0,
                                         haven::labelled(NA_real_),
                                         stroke_before_study)) %>%
    mutate(stroke_before_study = if_else(any(stroke_before_study %in% 1),
                                         1L, 0L)) %>%
    ungroup()

## For people who have wave 3 as first interview, find out if they were ever
## diagnosed with a stroke
wave3_stroke <- read_stata(file = "./sharew3/sharew3_rel8-0-0_hs.dta")
wave3_stroke <- left_join(wave3_stroke,
                          master_file[, c("mergeid", "int_year_w2")],
                          by = "mergeid")
wave3_stroke <- wave3_stroke %>%
    mutate(had_stroke_before =
               (sl_hs055d7_1 %in% 1) |
               (sl_hs055d7_2 %in% 1) |
               (sl_hs055d7_3 %in% 1))
wave3_stroke_ids <- wave3_stroke %>%
    filter(had_stroke_before == TRUE) %>%
    pull(mergeid)
model_data <- model_data %>%
    ungroup() %>%
    mutate(stroke_before_study =
               if_else(mergeid %in% wave3_stroke_ids & firstivw == 3,
                       haven::labelled(1), stroke_before_study))

## Use wave 3 information to figure out if stroke happened between wave 2 and
## wave 3 interview
wave3_stroke <- wave3_stroke %>%
    mutate(has_stroke_w3 =
               case_when(sl_hs055d7_1 %in% 1 & int_year_w2 > 0 & sl_hs059_1 >= int_year_w2 ~ 1L,
                         sl_hs055d7_1 %in% 1 & int_year_w2 > 0 & sl_hs059_1 < int_year_w2 ~ 0L,
                         sl_hs055d7_2 %in% 1 & int_year_w2 > 0 & sl_hs059_2 >= int_year_w2 ~ 1L,
                         sl_hs055d7_2 %in% 1 & int_year_w2 > 0 & sl_hs059_2 < int_year_w2 ~ 0L,
                         sl_hs055d7_3 %in% 1 & int_year_w2 > 0 & sl_hs059_3 >= int_year_w2 ~ 1L,
                         sl_hs055d7_3 %in% 1 & int_year_w2 > 0 & sl_hs059_3 < int_year_w2 ~ 0L))

## Add this information back into the model data
model_data[model_data$wave == 3 & model_data$mergeid %in% wave3_stroke$mergeid, "stroke_last_int"] <-
    wave3_stroke[, "has_stroke_w3"]

## Make lagged versions of adl, iadl, moderate physical activity (br016_) and
## vigorous physical activity (br015_)
model_data <- model_data %>%
    group_by(mergeid) %>%
    mutate(pa_mod = lag(br016_),
           pa_vig = lag(br015_),
           pa_mod_1 = br016_,
           pa_vig_1 = br015_,
           adl_1 = lag(adl),
           iadl_1 = lag(iadl)) %>%
    ungroup()

## Set up missing values for age
model_data <- model_data %>%
    mutate(age = ifelse(age_int < 0,
                        haven::labelled(NA_real_),
                        age_int))

## Recode binary stroke variable (i.e. stroke since last interview) so that no
## stroke since last interview is 0 and not 5
model_data <- model_data %>%
    mutate(stroke_last_int = if_else(stroke_last_int %in% 5,
                                     haven::labelled(0),
                                     stroke_last_int)) %>%
    mutate(stroke_last_int = if_else(!is.na(stroke_last_int) & stroke_last_int < 0,
                                     haven::labelled(NA_real_),
                                     stroke_last_int))

## Make country variable a factor
model_data <- model_data %>%
    mutate(country = factor(country,
                            levels = attr(country, "labels", exact = TRUE),
                            labels = names(attr(country, "labels", exact = TRUE))))

## Set missing values for lagged physical activity
model_data <- model_data %>%
    mutate(pa_mod = if_else(!is.na(pa_mod) & pa_mod < 0,
                            haven::labelled(NA_real_),
                            pa_mod),
           pa_vig = if_else(!is.na(pa_vig) & pa_vig < 0,
                            haven::labelled(NA_real_),
                            pa_vig),
           pa_mod_1 = if_else(!is.na(pa_mod_1) & pa_mod_1 < 0,
                            haven::labelled(NA_real_),
                            pa_mod_1),
           pa_vig_1 = if_else(!is.na(pa_vig_1) & pa_vig_1 < 0,
                            haven::labelled(NA_real_),
                            pa_vig_1))

## Reverse the physical activity variable
model_data <- model_data %>%
    mutate(pa_mod_rev = -(pa_mod - max(pa_mod, na.rm = TRUE)),
           pa_vig_rev = -(pa_vig - max(pa_vig, na.rm = TRUE)),
           pa_mod_rev_1 = -(pa_mod_1 - max(pa_mod_1, na.rm = TRUE)),
           pa_vig_rev_1 = -(pa_vig_1 - max(pa_vig_1, na.rm = TRUE)))

## Binary version of adl and iadl
model_data <- model_data %>%
    mutate(adl_bin = if_else(adl > 0, 1, 0),
           iadl_bin = if_else(iadl > 0, 1, 0),
           adl_1_bin = if_else(adl_1 > 0, 1, 0),
           iadl1_bin = if_else(iadl_1 > 0, 1, 0))

## Clean up body mass index (BMI)
model_data <- model_data %>%
    mutate(bmi_clean = if_else(!is.na(bmi) & bmi < 0,
                         haven::labelled(NA_real_),
                         bmi),
           bmi2_clean = if_else(!is.na(bmi2) & bmi2 < 0,
                         haven::labelled(NA_real_),
                         bmi2))

model_data <- model_data %>%
    mutate(bmi2_clean = factor(bmi2_clean,
                         levels = as.character(attr(x = bmi2, which = "labels", exact = TRUE)),
                         labels = names(attr(x = bmi2, which = "labels", exact = TRUE))))

model_data <- model_data %>%
    mutate(bmi2_clean = relevel(bmi2_clean, ref = "18.5-24.9 - normal"))

## Remove unused factor levels
model_data <- model_data %>%
    mutate(bmi2_clean = factor(bmi2_clean))

## Recode physical activity into binary variables
## Max a low and high cut-off version
model_data <- model_data %>%
    mutate(pa_bin_low_cut = if_else(pa_mod_rev == 0 | pa_vig_rev == 0, 1, 0),
           pa_bin_high_cut = if_else(pa_mod_rev == 3 | pa_vig_rev == 3, 0, 1),
           pa_bin_low_cut_1 = if_else(pa_mod_rev_1 == 0 | pa_vig_rev_1 == 0, 1, 0),
           pa_bin_high_cut_1 = if_else(pa_mod_rev_1 == 3 | pa_vig_rev_1 == 3, 0, 1))

## Clean up chronic conditions
model_data <- model_data %>%
  mutate(chronic = if_else(!is.na(chronic) & chronic < 0, NA_real_,
                           chronic),
         chronic2 = if_else(!is.na(chronic2) & chronic2 < 0, NA_real_,
                            chronic2)) %>%
  mutate(chronic1 = case_when(chronic %in% 0 ~ 0,
                              !(chronic %in% 0) ~ 1,
                              is.na(chronic) ~ NA))

## Binary variable if there was ever any stroke during follow-up
model_data <- model_data %>%
    group_by(mergeid) %>%
    mutate(ever_stroke = any(stroke_last_int %in% 1)) %>%
    ungroup()

### Matching of stroke survivors with stroke-free controls

## Function to get the value of a variable at the 1st interview
first_ivw_value <- function(first_ivw, wave_var, var) {
  stopifnot(all.equal(length(wave_var), length(var)),
            isTRUE(length(first_ivw) == 1))

  pos <- match(first_ivw, wave_var)

  if (is.na(pos)) {
    ret <- NA
  } else {
    ret <- var[[pos]]
  }
  return(ret)
}

## Function to simplify making baseline variables
make_baseline <- function(first_ivw, wave_var, var) {
  ifelse(all(is.na(first_ivw)), NA,
         first_ivw_value(unique(first_ivw), wave_var, var))
}

match_data <- model_data
master_file_m <- master_file %>%
  mutate(country = factor(country,
                            levels = attr(country, "labels", exact = TRUE),
                            labels = names(attr(country, "labels", exact = TRUE)))) %>%
  select(mergeid, matches("^interview_w[[:digit:]]$"))

match_data <- left_join(match_data, master_file_m)
match_data <- match_data %>%
  group_by(mergeid) %>%
  mutate(baseline_age = make_baseline(firstivw, wave, age)) %>%
  mutate(baseline_bmi = make_baseline(firstivw, wave, bmi2_clean)) %>%
  mutate(baseline_chronic1 = make_baseline(firstivw, wave, chronic1)) %>%
  mutate(baseline_chronic2 = make_baseline(firstivw, wave, chronic2)) %>%
  mutate(baseline_pa_bin_low_cut_1 = make_baseline(firstivw, wave, pa_bin_low_cut_1)) %>%
  mutate(baseline_pa_bin_high_cut_1 = make_baseline(firstivw, wave, pa_bin_high_cut_1)) %>%
  mutate(baseline_adl = make_baseline(firstivw, wave, adl_raw)) %>%
  mutate(baseline_iadl = make_baseline(firstivw, wave, iadl_raw)) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Drop respondents who had a stroke before the study
match_data <- match_data %>%
  filter(stroke_before_study == 0)

## Make factors
match_data <- match_data %>%
  mutate(gender = factor(gender,
                           levels = attr(gender, "labels", exact = TRUE),
                           labels = names(attr(gender, "labels", exact = TRUE))))

match_data <- match_data %>%
    mutate(across(c(
        "interview_w1", "interview_w2", "interview_w3", "interview_w4",
        "interview_w5", "interview_w6", "interview_w7", "interview_w8"
    ),
    ~ if_else(.x < 0, labelled(0), .x))) %>%
    mutate(across(c(
        "interview_w1", "interview_w2", "interview_w3", "interview_w4",
        "interview_w5", "interview_w6", "interview_w7", "interview_w8"
    ),
    ~ relevel(factor(.x), "1"))) %>%
    group_by(mergeid) %>%
    mutate(n_invw = sum(across(.cols = c("interview_w1", "interview_w2",
                                 "interview_w3", "interview_w4",
                                 "interview_w5", "interview_w6",
                                 "interview_w7", "interview_w8"), .fns = ~ .x %in% "1"))) %>%
    ungroup() %>%
    filter(n_invw >= 4)

match_data <- match_data %>%
    filter(complete.cases(match_data %>% select(ever_stroke, country, baseline_adl, baseline_iadl,
                                                gender, baseline_age,
                                                firstivw, n_invw,
                                                baseline_bmi, baseline_chronic2)))

match_res <- matchit(ever_stroke ~ country + baseline_adl + baseline_iadl + gender + baseline_age + firstivw +
                       n_invw + baseline_bmi + baseline_chronic2,
                      data = match_data, ratio = 5)
matched_sample <- match.data(match_res)

## Find the occurrence of stroke for matched
model_data <-
    model_data %>%
    group_by(mergeid) %>%
    mutate(when_stroke = ifelse(!is.na(match(1, stroke_last_int)),
                                wave[[match(1, stroke_last_int)]],
                                NA)) %>%
    ungroup()

model_data_matched <-
    model_data %>%
    filter(mergeid %in% matched_sample$mergeid)

model_data_matched <-
    left_join(model_data_matched,
              matched_sample[, c("mergeid", "weights", "subclass", "baseline_adl", "baseline_iadl", "baseline_age", "baseline_bmi", "baseline_chronic1", "baseline_chronic2",
                                 "baseline_pa_bin_low_cut_1", "baseline_pa_bin_high_cut_1")])
model_data_matched <-
    model_data_matched %>%
    group_by(subclass) %>%
    mutate(when_stroke_match = unique(when_stroke[!is.na(when_stroke)])) %>%
    ungroup()

model_data_matched <-
    model_data_matched %>%
    mutate(adl = if_else(adl < 0, NA_real_, adl)) %>%
    mutate(iadl = if_else(iadl < 0, NA_real_, iadl)) %>%
    mutate(stroke_1_0 = factor(!ever_stroke)) %>%
    mutate(ever_stroke = factor(ever_stroke)) %>%
    mutate(time_stroke = wave - when_stroke_match)

model_data_matched <-
    model_data_matched %>%
    mutate(wave = wave - 1)

#### Estimate main models
# Main models testing the effect of physical activity (baseline_pa_bin_high) on ADLs
res_dir <- as.character(Sys.Date())
dir.create(res_dir)

mod_2_adl_pa_high <- lmer(adl_raw ~ ever_stroke * baseline_pa_bin_high_cut_1 + wave +
                            I(wave^2) + baseline_age + gender +
                            max_edu + chronic2 +
                            (wave + I(wave^2) | mergeid) +
                            (1 | subclass),
                          data = model_data_matched,
                          weights = model_data_matched$weights
)

data_fit_adl_pa_high <- mod_2_adl_pa_high@frame

mod_1_adl_pa_high <- lmer(adl_raw ~ever_stroke * baseline_pa_bin_high_cut_1 + wave +
                            I(wave^2) +
                            (wave + I(wave^2) | mergeid) +
                            (1 | subclass),
                          data = data_fit_adl_pa_high,
                          weights = data_fit_adl_pa_high$weights
)

write.csv(summary(mod_1_adl_pa_high)$coefficients,
          file = file.path(res_dir, "mod_1_adl_pa_high_coefs.csv"))
write.csv(summary(mod_2_adl_pa_high)$coefficients,
          file = file.path(res_dir, "mod_2_adl_pa_high_coefs.csv"))

# Main models testing the effect of physical activity (baseline_pa_bin_high) on IADLs
mod_2_iadl_pa_high <- lmer(iadl_raw ~ever_stroke * baseline_pa_bin_high_cut_1 + wave +
                             I(wave^2) + baseline_age +
                             gender + max_edu + chronic2 +
                             (wave + I(wave^2) | mergeid) +
                             (1 | subclass),
                           data = model_data_matched,
                           weights = model_data_matched$weights
)

data_fit_iadl_pa_high <- mod_2_iadl_pa_high@frame

mod_1_iadl_pa_high <- lmer(iadl_raw ~ever_stroke * baseline_pa_bin_high_cut_1 + wave +
                             I(wave^2) +
                             (wave + I(wave^2) | mergeid) +
                             (1 | subclass),
                           data = data_fit_iadl_pa_high,
                           weights = data_fit_iadl_pa_high$weights
)

write.csv(summary(mod_1_iadl_pa_high)$coefficients,
          file = file.path(res_dir, "mod_1_iadl_pa_high_coefs.csv"))
write.csv(summary(mod_2_iadl_pa_high)$coefficients,
          file = file.path(res_dir, "mod_2_iadl_pa_high_coefs.csv"))

# Plot main analyses (pa_high)
plot_pa_high_data <- lapply(mget(ls(pattern = "mod.*high$")), ggeffects::ggpredict,
                            terms = c("wave[all]", "ever_stroke", "baseline_pa_bin_high_cut_1 [0, 1]"),
                            type = "fe"
)

plots_pa_high <- lapply(plot_pa_high_data, function(data) {
  ggplot(data = data, mapping = aes(x = x, y = predicted, colour = facet, linetype = group)) +
    geom_line()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet), alpha = 0.2) +
    scale_colour_hue (aesthetics = c("colour", "fill"), direction = -1, labels = c("Stroke-Free Controls", "Stroke Survivors")) +
    labs(x = "Survey Wave", y= "Functional Limitation", colour = "", linetype = "Low Physical Activity")
})

ggsave(plot = plots_pa_high[[1]] + ylab("ADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_1_adl_pa_high.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_high[[3]] + ylab("ADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_2_adl_pa_high.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_high[[2]] + ylab("IADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_1_iadl_pa_high.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_high[[4]] + ylab("IADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_2_iadl_pa_high.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

plots_pa_high

#### Estimate sensitivity models
# Using low cut-off  of physical activity instead of a high cut-off
# i.e., replacing baseline_pa_bin_high_cut_1 with baseline_pa_bin_low_cut_1

# Sensitivity models testing the effect of physical activity (baseline_pa_bin_low) on ADLs
mod_2_adl_pa_low <- lmer(adl_raw ~ever_stroke * baseline_pa_bin_low_cut_1 + wave +
                  I(wave^2) + baseline_age + gender + chronic2 + max_edu + 
                  (wave + I(wave^2) | mergeid) + (1 | subclass),
                  data = model_data_matched,
                  weights = model_data_matched$weights
)

data_fit_adl_pa_low <- mod_2_adl_pa_low@frame

mod_1_adl_pa_low <-
  lmer(adl_raw ~ever_stroke * baseline_pa_bin_low_cut_1 + wave +
         I(wave^2) +
         (wave + I(wave^2) | mergeid) +
         (1 | subclass),
       data = data_fit_adl_pa_low,
       weights = data_fit_adl_pa_low$weights
  )

write.csv(summary(mod_1_adl_pa_low)$coefficients,
          file = file.path(res_dir, "mod_1_adl_pa_low_coefs.csv"))
write.csv(summary(mod_2_adl_pa_low)$coefficients,
          file = file.path(res_dir, "mod_2_adl_pa_low_coefs.csv"))

# Sensitivity models testing the effect of physical activity (baseline_pa_bin_low) on IADLs
mod_2_iadl_low <- lmer(iadl_raw ~ever_stroke * baseline_pa_bin_low_cut_1 + wave +
                    I(wave^2) + baseline_age + gender + chronic2 + max_edu + 
                    (wave + I(wave^2) | mergeid) + (1 | subclass),
                    data = model_data_matched,
                    weights = model_data_matched$weights
)

data_fit_iadl_low <- mod_2_iadl_low@frame

mod_1_iadl_low <- lmer(iadl_raw ~ever_stroke * baseline_pa_bin_low_cut_1 + wave +
                     I(wave^2) +
                     (wave + I(wave^2) | mergeid) +
                     (1 | subclass),
                   data = data_fit_iadl_low,
                   weights = data_fit_iadl_low$weights
)

write.csv(summary(mod_1_iadl_low)$coefficients,
          file = file.path(res_dir, "mod_1_iadl_pa_low_coefs.csv"))
write.csv(summary(mod_2_iadl_low)$coefficients,
          file = file.path(res_dir, "mod_2_iadl_pa_low_coefs.csv"))

### Plot sensitivity models (pa_low)
plot_data <- lapply(mget(ls(pattern = "mod.*low$")), ggeffects::ggpredict,
                    terms = c("wave[all]", "ever_stroke", "baseline_pa_bin_low_cut_1 [0, 1]"),
                    type = "fe"
)

plots_pa_low <- lapply(plot_data, function(data) {
  ggplot(data = data, mapping = aes(x = x, y = predicted, colour = facet, linetype = group)) +
    geom_line()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet), alpha = 0.2) + # default conf.level = .95
    scale_colour_hue (aesthetics = c("colour", "fill"), direction = -1, labels = c("Stroke-Free Controls", "Stroke Survivors")) +
    labs(x = "Survey Wave", y= "Functional Limitation", colour = "", linetype = "Low Physical Activity")
})

ggsave(plot = plots_pa_low[[1]] + ylab("ADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_1_adl_pa_low.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_low[[3]] + ylab("ADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_2_adl_pa_low.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_low[[2]] + ylab("IADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_1_iadl_pa_low.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(plot = plots_pa_low[[4]] + ylab("IADL"),
       filename = file.path(res_dir, "stroke_test_plot_mod_2_iadl_pa_low.png"),
       device = "png",
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

plots_pa_low


##### Descriptive table ######
library("openxlsx")
ids <- Reduce(`union`, lapply(mget(ls(pattern = "data_fit")), `[[`, "mergeid"))

desc_data <- model_data_matched %>%
  filter(mergeid %in% !!ids)

## Missing control vs stroke
missing_adl_iadl <- desc_data %>%
  filter(!(wave %in% 2)) %>%
  group_by(mergeid) %>%
  summarise(n_obs_adl = sum(!is.na(adl)),
            n_obs_iadl = sum(!is.na(iadl)),
            n_obs_adl_raw = sum(!is.na(adl_raw)),
            n_obs_iadl_raw = sum(!is.na(iadl_raw)),
            had_stroke = unique(ever_stroke)) %>%
  select(-mergeid) %>%
  group_by(had_stroke) %>%
  summarise(avg_obs_adl = mean(n_obs_adl),
            sd_obs_adl = sd(n_obs_adl),
            avg_obs_iadl = mean(n_obs_iadl),
            sd_obs_iadl = sd(n_obs_iadl),
            avg_obs_adl_raw = mean(n_obs_adl_raw),
            sd_obs_adl_raw = sd(n_obs_adl_raw),
            avg_obs_iadl_raw = mean(n_obs_iadl_raw),
            sd_obs_iadl_raw = sd(n_obs_iadl_raw))

export_desc_vars <- createWorkbook()
addWorksheet(export_desc_vars, "Avg. Obs.")
writeData(export_desc_vars, sheet = "Avg. Obs.", x = missing_adl_iadl)

## Time constant variables
time_constant_vars <-
  c("ever_stroke",
    "baseline_adl",
    "baseline_iadl",
    "baseline_pa_bin_low_cut_1",
    "baseline_pa_bin_high_cut_1",
    "baseline_age",
    "gender",
    "country",
    "max_edu")

make_desc <- function(x) {
  if (isTRUE(length(unique(x)) <= 5) | is.factor(x)) {
    tmp <- table(x)
    empty <- tmp %in% 0
    tmp <- tmp[!empty]
    tmp1 <- sprintf("%.1f%%", prop.table(tmp) * 100)

    ret <- paste0(tmp, " (", tmp1, ")")

    dim(ret) <- c(length(ret), 1)
    ret <- cbind(names(tmp), ret)

    return(ret)
  } else {
    tmp <- sprintf("%.2f", mean(x, na.rm = TRUE))
    tmp1 <- sprintf("%.2f", sd(x, na.rm = TRUE))

    ret <- paste0(tmp, " (", tmp1, ")")
    ret <- cbind("", ret)

    return(ret)
  }
}

desc_data_constant <-
  lapply(time_constant_vars,
         function(x, data, func) {
           ret <- data %>% reframe(func(.data[[x]]))
           names(ret)[-1] <- "vals"

           return(ret)
         },
         data = desc_data %>%
           group_by(mergeid) %>%
           filter(row_number() == 1) %>%
           ungroup() %>%
           select(!!time_constant_vars) %>%
           group_by(ever_stroke),
         func = make_desc
  )

names(desc_data_constant) <- time_constant_vars

desc_data_constant <-
  cbind(rep(
    names(desc_data_constant),
    times = vapply(desc_data_constant, nrow, FUN.VALUE = integer(1L))
  ),
  Reduce(`rbind`, desc_data_constant))

desc_data_constant <-
  lapply(desc_data_constant,
         function(x) if (is.matrix(x)) as.data.frame(x) else x)
desc_data_constant <- data.frame(desc_data_constant)

prep_excel_export_constant <- desc_data_constant
prep_excel_export_constant[-match(unique(prep_excel_export_constant[, 1]),
                                  prep_excel_export_constant[, 1]), 1] <- ""
cell_merge_list_constant <-
  list(
    start_row = match(unique(desc_data_constant[, 1]), desc_data_constant[, 1]),
    rows_to_span = table(match(
      desc_data_constant[, 1], unique(desc_data_constant[, 1])
    ))
  )

addWorksheet(export_desc_vars, "Time Constant")
writeData(export_desc_vars, sheet = "Time Constant", x = prep_excel_export_constant, colNames = FALSE)
invisible(mapply(
  function(start_row, rows_to_span, wb, sheet) {
    row_calc <- seq_len(rows_to_span) - 1 + start_row
    mergeCells(
      wb = wb,
      sheet = sheet,
      cols = 1,
      rows = row_calc
    )
  },
  cell_merge_list_constant$start_row,
  cell_merge_list_constant$rows_to_span,
  MoreArgs = list(wb = export_desc_vars, sheet = "Time Constant")
))

#### Time-varying vars ####
time_varying_vars <- c("bmi2_clean", "chronic2")
desc_data_varying <-
  lapply(time_varying_vars,
         function(x, data, func) {
           ret <- data %>% reframe(func(.data[[x]]))
           names(ret)[-1] <- "vals"

           return(ret)
         },
         data = desc_data %>%
           group_by(ever_stroke) %>%
           select(!!time_varying_vars),
         func = make_desc
  )

names(desc_data_varying) <- time_varying_vars

desc_data_varying <-
  cbind(rep(
    names(desc_data_varying),
    times = vapply(desc_data_varying, nrow, FUN.VALUE = integer(1L))
  ),
  Reduce(`rbind`, desc_data_varying))

desc_data_varying <-
  lapply(desc_data_varying,
         function(x) if (is.matrix(x)) as.data.frame(x) else x)
desc_data_varying <- data.frame(desc_data_varying)

prep_excel_export_varying <- desc_data_varying
prep_excel_export_varying[-match(unique(prep_excel_export_varying[, 1]),
                                  prep_excel_export_varying[, 1]), 1] <- ""
cell_merge_list_varying <-
  list(
    start_row = match(unique(desc_data_varying[, 1]), desc_data_varying[, 1]),
    rows_to_span = table(match(
      desc_data_varying[, 1], unique(desc_data_varying[, 1])
    ))
  )

addWorksheet(export_desc_vars, "Time Varying")
writeData(export_desc_vars, sheet = "Time Varying", x = prep_excel_export_varying, colNames = FALSE)
invisible(mapply(
  function(start_row, rows_to_span, wb, sheet) {
    row_calc <- seq_len(rows_to_span) - 1 + start_row
    mergeCells(
      wb = wb,
      sheet = sheet,
      cols = 1,
      rows = row_calc
    )
  },
  cell_merge_list_varying$start_row,
  cell_merge_list_varying$rows_to_span,
  MoreArgs = list(wb = export_desc_vars, sheet = "Time Varying")
))

saveWorkbook(export_desc_vars, file = "test.xlsx", overwrite = TRUE)

## Before and after PA
## PA variables available
# - pa_mod --- lagged version
# - pa_vig --- lagged version
# - pa_mod_1
# - pa_vig_1
# - pa_mod_rev --- scale reversed, higher values = more PA; lagged version
# - pa_vig_rev --- scale reversed, higher values = more PA; lagged version
# - pa_mod_rev_1 --- scale reversed, higher values = more PA
# - pa_vig_rev_1 --- scale reversed, higher values = more PA
# - pa_bin_low_cut -- lagged version
# - pa_bin_high_cut -- lagged version
# - pa_bin_low_cut_1
# - pa_bin_high_cut_1
# - baseline_pa_bin_low_cut_1
# - baseline_pa_bin_high_cut_1

ids <- Reduce(`union`, lapply(mget(ls(pattern = "data_fit")), `[[`, "mergeid"))
before_after <- model_data_matched %>%
  filter(mergeid %in% !!ids) %>%
  group_by(mergeid)

before_after <-
  before_after %>%
  mutate(baseline_pa_mod_rev_1 = nth(pa_mod_rev_1, match(unique(firstivw) - 1, wave))) %>%
  mutate(baseline_pa_vig_rev_1 = nth(pa_vig_rev_1, match(unique(firstivw) - 1, wave)))

before_after <-
  before_after %>%
  filter(ever_stroke %in% TRUE) %>%
  mutate(stroke_wave = nth(wave, unique(when_stroke))) %>%
  mutate(pa_before_low_cut = if_else(
    stroke_wave != 3,
    nth(pa_bin_low_cut_1, unique(when_stroke) - 1),
    nth(pa_bin_low_cut_1, unique(when_stroke) - 2)
  )) %>%
  mutate(pa_after_low_cut = nth(pa_bin_low_cut_1, unique(when_stroke))) %>%
  mutate(pa_before_high_cut = if_else(
    stroke_wave != 3,
    nth(pa_bin_high_cut_1, unique(when_stroke) - 1),
    nth(pa_bin_high_cut_1, unique(when_stroke) - 2)
  )) %>%
  mutate(pa_after_high_cut = nth(pa_bin_high_cut_1, unique(when_stroke))) %>%
  filter(row_number() == 1) %>%
  select(
    mergeid,
    pa_after_high_cut,
    pa_before_high_cut,
    pa_after_low_cut,
    pa_before_low_cut,
    when_stroke,
    stroke_wave,
    baseline_pa_bin_low_cut_1,
    baseline_pa_bin_high_cut_1,
    baseline_pa_mod_rev_1,
    baseline_pa_vig_rev_1
  )
