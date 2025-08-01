setwd('') # Set working dir
library(tidyverse)
# Related constants functions -------------------------------------------------------
N_G_POWER = 921 # Required sample size estimated through Gpower analysis
DESCRIPTION_ALL_PHE <<- read_csv('description/descriptions_all_phenotypes.csv')
PATH.GGSEG.HO <- 'description/ggsesg/ggseg_HO.csv'
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}

index_measure_2_description <- function(index, measure) {
  DESCRIPTION_ALL_PHE %>% filter(index == index & measure == measure) %>% pull(description) %>% unique()
}

index2label <- function(index) {
  if(index > 0) {
    return(acronym.schafer216[index])
  }
  if(index == 0) {
    return('whole-brain')
  }
  if(index <= -100) {
    return(label.intermodules[-index-99])
  }
  if(index <= -10) {
    return(c("DAN", "DMN", "ECN", "LIM", "SMT", "VAN", "VIS", "SUB")[-(index/10)])
  }
  return(NA_character_)
}

remove_leading_zero <- function(x) {
  # Use a regular expression to match the leading zero(s) in the number
  # and replace it with an empty string if it's directly before a decimal point,
  # taking into account possible negative sign
  sub("^(-?)0+(\\.)", "\\1\\2", x)
}

# Function for conditional number formatting
format_number <- function(x, decimal = 2, leadingZero = T) {
  if(!leadingZero) {
    case_when(
      abs(x) < 0.0099 ~ formatC(x, format = "e", digits = decimal - 1), # scientific notation for very small numbers. e.g. -1.9e-04
      abs(x) < 0.099 ~ round(x, digits = decimal + 1) %>% as.character() %>% remove_leading_zero(), # e.g. 0.012
      abs(x) <= 1 ~ round(x, digits = decimal) %>% as.character() %>% remove_leading_zero(),
      abs(x) < 9999 ~ round(x, digits = decimal) %>% as.character(),
      TRUE ~ formatC(x, format = "e", digits = decimal - 1)
    )
  }
  else {
    case_when(
      abs(x) < 0.0099 ~ formatC(x, format = "e", digits = decimal), # scientific notation for very small numbers. e.g. -1.9e-04
      abs(x) < 0.099 ~ round(x, digits = decimal + 1) %>% as.character(), # e.g. 0.012
      abs(x) <= 1 ~ round(x, digits = decimal) %>% as.character(),
      abs(x) < 9999 ~ round(x, digits = decimal) %>% as.character(),
      TRUE ~ formatC(x, format = "e", digits = decimal - 1)
    )
  }
}

# Tidy data ---------------------------------------------------------------

# 1. Prepare UK Biobank phenotype data (cognition and IDPs) as a long table with the following columns:
# - eid: Participant identifier (UK Biobank ID)
# - atlas: Variable domain or anatomical atlas, e.g., "cognition", "HO", "diffusion"
# - weight: Measurement type within atlas, e.g., "CFA", "volume", "diffusion"
# - index: Index of brain regions or cognitive functions
# - measure: Specific measure type, e.g., "cognition", "volume_first_normalised",
#            "volume_whole_brain_normalised", "weighted_fa", "mean_icvf", etc.
# - description: Phenotype description (from 'description/descriptions_all_phenotypes.csv'),
#                e.g., "Mean ISOVF in cingulum cingulate gyrus on FA skeleton (right)"
# - value: Phenotype value
#
# Note: UK Biobank data is available upon approved application:
#       https://www.ukbiobank.ac.uk
#
# This processed dataset is saved as 'data/tidy_imaging_cognition.rds' and used in subsequent scripts.


# 2. Prepare demographic data in following cols:
# Format demographic and genetic covariate data into the following columns:
# - eid: Participant identifier
# - ethnicity: Ethnic background, e.g., "UKB-UC"  "UKB-MIX" "UKB-CC" 
# - apoe: APOE genotype, e.g., "e2e3" "e3e4" "e3e3" "e2e2" "e4e4"
# - PRS-UKB, PRS-APin, PRS-APex: Polygenic risk scores
# - age: Precise age at scan, e.g. 72.23454
# - sex: Biological sex
# - education: Education score
# - BMI: Body mass index
# - townsend: Townsend deprivation index
# - date.scanning: MRI/PET scan date
# - center.newcastle / center.cheadle / center.reading: Scanning site dummy variables
# - hypertensive: Hypertension status
# - diabetic: Diabetes status
# - genePC1–genePC40: Genetic principal components for population stratification
#
# Note: UK Biobank data is available upon approved application:
#       https://www.ukbiobank.ac.uk
#
# This processed dataset is saved as 'data/demo.csv' and used in subsequent scripts.


# 3. Inspect sample sizes by age window length for subsequent SW analysis
# VERY IMPORTANT
df.age.uc <- read_csv('data/demo.csv') %>% 
  filter(ethnicity == 'UKB-UC') %>% 
  select(eid, age, PRS = `PRS-UKB`) %>% 
  mutate(tertile = ntile(PRS, 3)) %>%
  mutate(tertile = if_else(tertile == 1, 'Low', if_else(tertile == 2, 'Medium', 'High'))) %>% 
  mutate(tertile = factor(tertile, levels = c('Low', 'Medium', 'High')))

sample_sizes_uc <- map_dfr(
  sort(c(seq(0.5, 5, by=0.1))),
  function(WIDTH) {
    age.window <- tibble(
      age.min = seq(45, 85, by = 0.1),
      age.max = age.min + WIDTH
    ) %>% filter(age.max < 85)
    map2_dfr(age.window$age.min, age.window$age.max, function(age.min, age.max){
      DF <- df.age.uc %>% filter(between_interval(age, age.min, age.max))
      tibble(
        width = WIDTH,
        age = median(DF$age),
        age.min = age.min,
        age.max = age.max,
        `Total group` = nrow(DF),
        `Low-risk group` = nrow(DF %>% filter(tertile == 'Low')),
        `Med-risk group` = nrow(DF %>% filter(tertile == 'Medium')),
        `High-risk group` = nrow(DF %>% filter(tertile == 'High'))
      )
    })
  }
)

age_span <- sample_sizes_uc %>% 
  filter(`Total group` >= N_G_POWER) %>% 
  group_by(width) %>% 
  summarise(
    `Min(age)` = min(age),
    `Max(age)` = max(age)
  ) %>% 
  ungroup()

sample_sizes_uc <- na.omit(sample_sizes_uc)
write_csv(sample_sizes_uc, 'data/sample_size_by_window.csv')


# Effects after 75 --------------------------------------------------------
create_when_absent('results/201-IDP-PRS-after-75')
if(!file.exists('results/201-IDP-PRS-after-75/result_report_all.rds')) {
  data.demo <- read_csv('data/demo.csv') %>% 
    filter(age>=75) %>% 
    filter(ethnicity=='UKB-UC') %>% na.omit() %>% 
    mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
    pivot_longer(c('PRS-UKB', 'PRS-APin', 'PRS-APex'), names_to = 'PRS.type', values_to = 'PRS') 
  data.idp <- read_rds('data/tidy_imaging_cognition.rds') %>% 
    filter(eid %in% data.demo$eid) %>% 
    left_join(data.demo) %>% 
    group_by(atlas, weight, index, measure, PRS.type) %>% 
    nest()
  
  future::plan(future::multisession, workers = 8)
  lm.result <- furrr::future_pmap_dfr(
    list(
      data.idp$atlas, data.idp$weight, data.idp$index, data.idp$measure, data.idp$PRS.type, data.idp$data
    ), 
    function(iAtlas, iWeight, iIndex, iMeasure, iPRS, iData){
      df.std <- iData %>% 
        na.omit() %>% 
        mutate(across(c(value, age, education, BMI, townsend, date.scanning, genePC1:genePC40, PRS), ~scale(.)[,]))
      # Try fitting the linear model and handling errors
      fit.model <- tryCatch(
        lm(df.std, formula = value ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
        error = function(e) NULL
      )
      # If lm() fails, return missing values for linear model fitting
      if (is.null(fit.model)) {
        result.lm <- tibble(
          lm_p = NA_real_,
          lm_beta = NA_real_,
          lm_se = NA_real_,
          lm_statistic = NA_real_
        )
      } else {
        result.lm <- broom::tidy(fit.model) %>% 
          filter(term == 'PRS') %>% 
          select(lm_p = p.value, lm_beta = estimate, lm_se = std.error, lm_statistic = statistic) 
      }
      cbind(
        atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, PRS = iPRS,
        result.lm
      )
    })
  
  result.report <- lm.result %>% 
    mutate(report = glue::glue('β = {format_number(lm_beta)}, p = {format_number(lm_p)}'))
    
  result.report <- lm.result %>% mutate(
    report = glue::glue('β = {format_number(lm_beta)}, p = {format_number(lm_p)}')
  ) %>% 
    mutate(description = map2_chr(index, measure, function(x, y){DESCRIPTION_ALL_PHE %>% filter(index == x, measure == y) %>% pull(description)}))
  write_rds(result.report, 'results/201-IDP-PRS-after-75/result_report_all.rds')
} else {
  result.report <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds')
}


# ggsesg-visualisation based on normalized FAST-------
temp <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds') %>% filter(measure == 'volume_fast_normalised') %>% filter(lm_p<0.05) %>% filter(PRS == 'PRS-UKB')
temp <- result.report %>% filter(measure == 'volume_first_normalised') %>% filter(lm_p<0.05) %>% filter(PRS == 'PRS-UKB')
temp <- result.report %>% filter(measure == 'cognition') %>% filter(lm_p<0.05) %>% filter(PRS == 'PRS-UKB')

# Cortical with ggseg
df.vis.ggseg <- temp %>% mutate(
  beta = lm_beta,
  logp = -log10(lm_p),
  index = index,
  .keep = 'none'
) %>% right_join(read_csv(PATH.GGSEG.HO))
library(ggseg)
library(ggsegHO)
ggseg(
  .data = df.vis.ggseg,
  atlas = 'hoCort',
  mapping = aes(fill=logp),
  size = .1,
  color = 'black',
  position = 'stacked'
) +
  scale_fill_gradient2(high = "#093F9E", low = '#C30600', na.value = 'white')


# Subcortical with points ---------
df.vis.sub <- temp %>% filter(index>96)

index.sub <- df.vis.sub %>% pull(index)
df.vis.sub <- result.report %>% filter(measure == 'volume_fast_normalised') %>% 
  #filter(lm_p<0.05) %>% 
  filter(index %in% index.sub) %>% 
  mutate(region = description %>%
           str_remove("^Volume of grey matter in ") %>%  # Remove prefix
           str_remove(" \\(normalised\\)$"))             # Remove suffix
  

ggplot(df.vis.sub, aes(x = region, y = lm_beta, color = PRS)) +
  geom_point() +
  geom_errorbar(aes(ymin = lm_beta - lm_se, ymax = lm_beta + lm_se), width = 0.2) +
  theme_classic() +
  facet_wrap(~PRS) +
  labs(x = 'Subcortical region', y = 'β after 75 years', fill = 'Weight') +
  theme(legend.position = 'right') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  #scale_fill_gradient2(low = "#C30600", mid = 'white', high = '#093F9E', na.value = 'white') +
  theme(
    axis.line = element_line(colour = 'black', linewidth = 0.3),
    axis.ticks = element_line(colour = 'black', linewidth = 0.3),
    axis.text.x = element_text(angle = 45, color = 'black', size = 8, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8),
    title = element_text(color = 'black', size = 6),
    axis.title.x = element_text(color = 'black', size = 8),
    axis.title.y = element_text(color = 'black', size = 8)
  ) +
  ggsci::scale_color_npg()




# Sensitivity-analysis-FIRST ---------
temp <- result.report %>% filter(measure == 'volume_first_normalised') %>% filter(lm_p<0.05)

df.vis.sub <- result.report %>% filter(measure == 'volume_first_normalised') %>% 
  #filter(lm_p<0.05) %>% 
  filter(index %in% temp$index) %>% 
  mutate(region = description %>%
           str_remove("^Volume of ") %>%  # Remove prefix
           str_remove(" \\(normalised\\)$"))             # Remove suffix
ggplot(df.vis.sub, aes(x = region, y = lm_beta, color = PRS)) +
  geom_point() +
  geom_errorbar(aes(ymin = lm_beta - lm_se, ymax = lm_beta + lm_se), width = 0.2) +
  theme_classic() +
  facet_wrap(~PRS) +
  labs(x = 'Subcortical region', y = 'β after 75 years', fill = 'Weight') +
  theme(legend.position = 'right') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme(
    axis.line = element_line(colour = 'black', linewidth = 0.3),
    axis.ticks = element_line(colour = 'black', linewidth = 0.3),
    axis.text.x = element_text(angle = 45, color = 'black', size = 8, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8),
    title = element_text(color = 'black', size = 6),
    axis.title.x = element_text(color = 'black', size = 8),
    axis.title.y = element_text(color = 'black', size = 8)
  ) +
  ggsci::scale_color_npg()


# Heatmap for mean diffusion measures ------------------
df.vis.diffusion <- result.report %>% 
  filter(measure %in% c("mean_fa", "mean_md", "mean_mo", "mean_icvf", "mean_od", "mean_isovf")) %>% 
  mutate(tract = map2_chr(
    index, measure, 
    function(x, y){
      DESCRIPTION_ALL_PHE %>% 
        filter(index == x & measure == y) %>% 
        mutate(hemi = ifelse(
          hemi %in% c('L', 'R'),
          str_c(' (', hemi, ')'),
          ''
        )) %>% 
        mutate(tract = str_c(acronym, hemi)) %>% 
        pull(tract)
    }
  )) %>% 
  mutate(measure = toupper(str_replace(measure, "\\mean_", ""))) %>% 
  mutate(measure = factor(measure, levels=c("FA", "MD", "ICVF", "ISOVF", "MO", "OD")))

df.vis.diffusion.star <- df.vis.diffusion %>% mutate(
  star = case_when(
    lm_p < 0.0001 ~ '****',
    lm_p < 0.001 ~ '***',
    lm_p < 0.01 ~ '**',
    lm_p < 0.05 ~ '*',
    TRUE ~ ''
  )
)


ggplot(df.vis.diffusion, aes(x = tract, y = measure, fill = lm_beta)) +
  geom_tile() +
  geom_text(data = df.vis.diffusion.star, mapping = aes(label = star), size = 3.5, nudge_y = -0.15, color = 'white') +
  scale_fill_gradientn(colors = c('darkred', 'lightpink', "white", "lightblue", "#232E73"), 
                       values = scales::rescale(c(-1, -0.75, -0.2, 0, 0.2, 0.75, 1))) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, color = "black",vjust = 0.8)
  ) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20)) +
  facet_wrap(~PRS, ncol = 1)


# Heatmap for weighted diffusion measures ------------------
df.vis.diffusion <- result.report %>% 
  filter(measure %in% c("weighted_fa", "weighted_md", "weighted_mo", "weighted_icvf", "weighted_od", "weighted_isovf")) %>% 
  mutate(tract = map2_chr(
    index, measure, 
    function(x, y){
      DESCRIPTION_ALL_PHE %>% 
        filter(index == x & measure == y) %>% 
        mutate(hemi = ifelse(
          hemi %in% c('L', 'R'),
          str_c(' (', hemi, ')'),
          ''
        )) %>% 
        mutate(tract = str_c(acronym, hemi)) %>% 
        pull(tract)
    }
  )) %>% 
  mutate(measure = toupper(str_replace(measure, "\\weighted_", ""))) %>% 
  mutate(measure = factor(measure, levels=c("FA", "MD", "ICVF", "ISOVF", "MO", "OD")))

df.vis.diffusion.star <- df.vis.diffusion %>% mutate(
  star = case_when(
    lm_p < 0.0001 ~ '****',
    lm_p < 0.001 ~ '***',
    lm_p < 0.01 ~ '**',
    lm_p < 0.05 ~ '*',
    TRUE ~ ''
  )
)


ggplot(df.vis.diffusion, aes(x = tract, y = measure, fill = lm_beta)) +
  geom_tile() +
  geom_text(data = df.vis.diffusion.star, mapping = aes(label = star), size = 3.5, nudge_y = -0.15, color = 'white') +
  scale_fill_gradientn(colors = c('darkred', 'lightpink', "white", "lightblue", "#232E73"), 
                       values = scales::rescale(c(-1, -0.75, -0.2, 0, 0.2, 0.75, 1))) +
  #geom_text(data = data.highlight, aes(label = round(value, 2)), color = 'white', size = 2.5) +
  #geom_rect(data= data.highlight, size=1, fill=NA, colour="black") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, color = "black",vjust = 0.8)
  ) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20)) +
  facet_wrap(~PRS, ncol = 1)
# Generates aging trajectories --------------------------------------------

path.plot <- '' %>% create_when_absent() # Customized output dir

# Gets the age windows along the age distribution
seq_along_age <- function(vector, interval, start = 1) {
  # Initialize the subvector with the first element
  subvector <- c(vector[start])
  prev_element <- vector[start]
  
  # Iterate through the vector to find elements nearest the specified interval
  for (i in (start + 1):length(vector)) {
    if ( abs(vector[i] - prev_element) >= interval ) {
      subvector <- c(subvector, vector[i])
      prev_element <- vector[i]  # Update prev_element
    }
  }
  
  return(subvector)
}

between_interval <- function(x, a, b) {
  dplyr::between(x, a, b) | dplyr::between(x, b, a)
}

generate_age_bins <- function(age, start = 1, bin_step, size.window) {
  tibble(age.left = age) %>% 
    arrange(desc(age.left)) %>% # Change it to arrange(age) if your study focus on younger participants 
    mutate(age.right = lead(age.left,n=size.window)) %>% 
    filter(age.left %in% seq_along_age(.$age.left, interval = bin_step)) %>% # Gets the age windows along the age distribution
    mutate(rn = row_number()) %>% # row_number
    filter(rn <= ( nrow(drop_na(.)) + 1) ) %>% 
    mutate(age.right = ifelse(is.na(age.right), min(age), age.right)) %>% # Change it to max(age) if your study focus on younger participants 
    select(age.left, age.right) %>% 
    # only the first observation for each unique age.left value will be kept, and among observations with the same age.left, the one with the smallest age.right will be retained
    arrange(age.left, age.right) %>%
    distinct(age.left, .keep_all = T) 
}

gaussian_smooth <- function(vec, n) {
  width.kernel = 2*n+1
  # Normalized Gaussian kernel of 2n+1
  # Scaling factors inside the kernel are set as suggested by references
  kernel.gaussian <- exp(-(3 * seq(-n, n) / n)^2) / sum(exp(-(3 * seq(-n, n) / n)^2)) 
  vec.pad <- c(rep(vec[1], n), vec, rep(vec[length(vec)], n))
  slider::slide(vec.pad, ~.x, .before = width.kernel - 1, .complete = TRUE) %>% 
    keep(~!is.null(.)) %>% 
    map_dbl(function(x) sum(x*kernel.gaussian))
}

sliding_window <- function(df, size.window, move.age.bin, n.kernel) {
  # Generates age bins  
  bins <- generate_age_bins(df$age, bin_step = move.age.bin, size.window = size.window)
  result <- map2_dfr(bins$age.left, bins$age.right, function(x, y) {
    temp <- df %>% filter(between_interval(age, x, y))
    return(
      tibble(
        age = median(temp$age),
        value.raw = mean(temp$value, na.rm = T),
        se.raw = sd(temp$value, na.rm = T) / sqrt(length(temp$value))
      ) %>% na.omit() 
    )
  }) %>% 
    # Gaussian smoothing
    mutate(
      value = gaussian_smooth(value.raw, n = n.kernel),
      se = gaussian_smooth(se.raw, n = n.kernel),
      .after = age
    )
  return(result)
}
demo.uc <- read_csv('data/demo.csv') %>% 
  filter(ethnicity=='UKB-UC') %>% na.omit() %>% 
  mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
  pivot_longer(c('PRS-UKB', 'PRS-APin', 'PRS-APex'), names_to = 'PRS.type', values_to = 'PRS') %>% 
  group_by(PRS.type) %>% 
  mutate(tertile = ntile(PRS, 3)) %>%
  mutate(tertile = if_else(tertile == 1, 'Low', if_else(tertile == 2, 'Medium', 'High'))) %>% 
  mutate(tertile = factor(tertile, levels = c('Low', 'Medium', 'High'))) %>% 
  ungroup()

data.idp.nest <- read_rds('data/tidy_imaging_cognition.rds') %>% 
  filter(eid %in% demo.uc$eid) %>% 
  select(-description) %>% 
  na.omit() %>%
  left_join(demo.uc %>% select(eid, age, PRS.type, tertile)) %>%
  group_by(atlas, weight, measure, index, PRS.type, tertile) %>% nest()
span = 0.1 # Span of sliding window
step = 0.5 # Move of age bins
kernel = 20 # Gaussian kernel of smoothing

future::plan(future::multisession, workers = 10)
data.trajectories <- furrr::future_pmap_dfr(
  list(data.idp.nest$atlas, data.idp.nest$weight, data.idp.nest$measure, data.idp.nest$index, data.idp.nest$PRS.type, data.idp.nest$tertile, data.idp.nest$data),
  function(iAtlas, iWeight, iMeasure, iIndex, iPRS, iTertile, iDF) {
    size_window <- ceiling( nrow(iDF) / 3 * span)
    result <- iDF %>% select(age, value) %>% sliding_window(move.age.bin = step, n.kernel = kernel, size.window = size_window)
    return(
      cbind(
        atlas = iAtlas, weight = iWeight, measure = iMeasure, index = iIndex, PRS = iPRS, tertile = iTertile,
        result
      )
    )      
  }
)

write_rds(data.trajectories, '') # Customized output dir
data.trajectories.nest <- data.trajectories %>% 
  group_by(atlas, weight, measure, index) %>% nest()

future::plan(future::multisession, workers = 10)
furrr::future_pwalk(
  list(data.trajectories.nest$atlas, data.trajectories.nest$weight, data.trajectories.nest$measure, data.trajectories.nest$index, data.trajectories.nest$data),
  function(iAtlas, iWeight, iMeasure, iIndex, iDF) {
    title <- DESCRIPTION_ALL_PHE %>% filter(index == iIndex & measure == iMeasure) %>% pull(description)
    p <- ggplot(data = iDF, mapping = aes(x = age, y = value, color = tertile, fill = tertile)) +
      geom_line(size = 0.5) +
      geom_ribbon(mapping = aes(ymin = value-se, ymax = value+se), alpha = 0.05, color = NA, position = "identity") +
      geom_point(mapping = aes(y = value), size = 0.3) +
      theme_classic() +
      scale_fill_manual(values = c('High' = '#CB3746', 'Medium' = '#F3D01F', 'Low' = '#18489D')) +
      scale_color_manual(values = c('High' = '#CB3746', 'Medium' = '#F3D01F', 'Low' = '#18489D')) +
      geom_vline(xintercept = 75, linetype = "dashed", color = "red") +
      facet_wrap(~PRS, nrow = 1, scales = 'free') +
      ylab(title) +
      xlab('Age (years)') +
      theme(
        axis.line = element_line(colour = 'black', linewidth = 0.3),
        axis.ticks = element_line(colour = 'black', linewidth = 0.3),
        axis.text.x = element_text(color = 'black', size = 8),
        axis.text.y = element_text(color = 'black', size = 8),
        title = element_text(color = 'black', size = 6),
        axis.title.x = element_text(color = 'black', size = 8),
        axis.title.y = element_text(color = 'black', size =4)
      ) +
      scale_x_continuous(breaks = seq(50, 80, by = 5)) 
    ggsave(plot = p, filename = str_c(path.plot, '/', iAtlas, '_', iWeight, '_', iMeasure, '_', iIndex, '.pdf'), width = 6.7, height = 2.1)
  }
)


# Generates SW statistics ---------------------------------------
path.result <- create_when_absent('') # Customized output dir
if(file.exists(str_c(path.result, '/intra_age_statistics.rds'))) {
  intra_age_statistics <- read_rds(str_c(path.result, '/intra_age_statistics.rds'))
} else {
  demo.uc <- read_csv('data/demo.csv') %>% 
    filter(ethnicity=='UKB-UC') %>% na.omit() %>% 
    mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
    mutate(Age = age) %>% # for sliding window
    mutate(across(c(age, education, BMI, townsend, date.scanning, genePC1:genePC40, 'PRS-UKB', 'PRS-APin', 'PRS-APex'), ~scale(.)[,])) %>% 
    pivot_longer(c('PRS-UKB', 'PRS-APin', 'PRS-APex'), names_to = 'PRS.type', values_to = 'PRS') %>% 
    group_by(PRS.type) %>% 
    mutate(tertile = ntile(PRS, 3)) %>%
    mutate(tertile = if_else(tertile == 1, 'Low', if_else(tertile == 2, 'Medium', 'High'))) %>% 
    mutate(tertile = factor(tertile, levels = c('Low', 'Medium', 'High'))) %>% 
    ungroup()
  
  data.idp <- read_rds('data/tidy_imaging_cognition.rds') %>% 
    filter(eid %in% demo.uc$eid)
  
  data.idp.nest <- data.idp %>% 
    na.omit() %>%
    group_by(atlas, weight, measure, index) %>% nest()
  
  # Computes intra-age statistics
  window_age_SW <- read_csv('data/sample_size_by_window.csv') %>% 
    filter(width == 3) %>% 
    filter(`Total group` > N_G_POWER) %>% 
    select(width, age, age.min, age.max)
  
  future::plan(future::multisession, workers = 10)
  intra_age_statistics <- furrr::future_pmap_dfr(
    list(data.idp.nest$atlas, data.idp.nest$weight, data.idp.nest$index, data.idp.nest$measure, data.idp.nest$data),
    function(iAtlas, iWeight, iIndex, iMeasure, iData){
      iDF <- iData %>% left_join(demo.uc) %>% group_by(PRS.type) %>% nest()
      map2_dfr(iDF$PRS.type, iDF$data, function(type.PRS, df){
        df <- df %>% mutate(value.std = scale(value)[,])
        pmap_dfr(
          list(window_age_SW$width, window_age_SW$age, window_age_SW$age.min, window_age_SW$age.max),
          function(iWidth, iAge, iAge.min, iAge.max){
            df.windowed <- df %>% filter(Age >= iAge.min & Age <= iAge.max)
            # First run linear models
            # Try fitting the linear model and handling errors
            fit.model <- tryCatch(
              lm(df.windowed, formula = value.std ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
              error = function(e) NULL
            )
            # If lm() fails, return missing values for linear model fitting
            if (is.null(fit.model)) {
              result.lm <- tibble(
                lm_p = NA_real_,
                lm_beta = NA_real_,
                lm_se = NA_real_,
                lm_statistic = NA_real_
              )
            } else {
              result.lm <- broom::tidy(fit.model) %>% 
                filter(term == 'PRS') %>% 
                select(lm_p = p.value, lm_beta = estimate, lm_se = std.error, lm_statistic = statistic) 
            }
            
            # Next compare between-group differences
            df.diff <- df.windowed %>% select(value, tertile) %>% filter(value > 0) %>%  na.omit()
            value.low <- df.diff %>% filter(tertile == 'Low') %>% pull(value)
            value.mid <- df.diff %>% filter(tertile == 'Medium') %>% pull(value)
            value.high <- df.diff %>% filter(tertile == 'High') %>% pull(value)
            result.diff <- map2_dfr(
              list(
                list(value.high, value.low),
                list(value.high, value.mid),
                list(value.mid, value.low)
              ),
              list('HL', 'HM', 'ML'),
              function(pair, label_pair) {
                a <- pair[[1]]
                b <- pair[[2]]
                a = keep(a, ~!is.na(.))
                b = keep(b, ~!is.na(.))
                t.test.coef <- t.test(a, b) %>% broom::tidy()
                cohen.d <- effectsize::cohens_d(a, b) # Cohen's d with pooled sd, 
                response.ratio <- effectsize::means_ratio(a, b) # Response ratio
                tibble(
                  label = label_pair,
                  mean_a = mean(a),
                  mean_b = mean(b),
                  cohen_d = cohen.d %>% pull(Cohens_d),
                  cohen_d_low = cohen.d %>% pull(CI_low),
                  cohen_d_high = cohen.d %>% pull(CI_high),
                  response_ratio = response.ratio %>% pull(Means_ratio_adjusted),
                  response_ratio_low = response.ratio %>% pull(CI_low),
                  response_ratio_high = response.ratio %>% pull(CI_high),
                  relative_diff = 1 - response_ratio,
                  t_test_p_value = t.test.coef %>% pull(p.value),
                  t_test_df = t.test.coef %>% pull(parameter),
                  t_test_t = t.test.coef %>% pull(statistic)
                )
                
              }
            ) %>% 
              pivot_longer(cols = -label, names_to = "variable", values_to = "value") %>%
              pivot_wider(names_from = c("label", "variable"), values_from = "value")
            
            # return the result
            cbind(
              atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure,
              age_window = iWidth, age = iAge,
              PRS = type.PRS,
              result.lm,
              result.diff
            ) %>% return()
          }
        )
      })
    }
  )
  
  write_rds(intra_age_statistics, '') # Customized output dir
}
# Fitting age related statistics with GAM ---------------------------------
# Model edge-wise age-related changes with GAM
path.result <- create_when_absent('') # Customized output dir

if(file.exists(str_c(path.result, '/GAM_fitting_raw.rds'))) {
  result.GAM <- read_rds(str_c(path.result, '/GAM_fitting_raw.rds'))
} else {
  result.SW <- read_rds('results/intra_age_statistics.rds')
  sw.nest <- result.SW %>%
    select(atlas, weight, index, measure, age_window, age, PRS, lm_beta, HL_cohen_d, HL_relative_diff, ML_cohen_d, ML_relative_diff, HM_cohen_d, HM_relative_diff) %>% 
    na.omit() %>% 
    pivot_longer(lm_beta:HM_relative_diff, names_to = 'statistics_type', values_to = 'statistics') %>% 
    group_by(atlas, weight, index, measure, age_window, PRS, statistics_type) %>% 
    nest()
  future::plan(future::multisession, workers = 10)
  result.GAM <- furrr::future_pmap_dfr(
    list(sw.nest$atlas, sw.nest$weight, sw.nest$index, sw.nest$measure, sw.nest$age_window, sw.nest$PRS, sw.nest$statistics_type, sw.nest$data),
    function(iAtlas, iWeight, iIndex, iMeasure, iAge_Window, iPRS_Type, iStatistics, iDF) {
      iDF <- na.omit(iDF) %>% distinct()
      # 1. Perform Mann-Kendall trend test
      mk <- trend::mk.test(iDF %>% arrange(age) %>% pull(statistics))
      mk.z <- mk$statistic %>% unname()
      mk.p <- mk$p.value
      mk.tau <- (mk$estimates)['tau'] %>% unname()
      
      # 2. Fit the GAMs: statistics ~ age
      model.gam <- mgcv::gam(statistics ~ s(age, k = 4, bs = 'ps'), data = iDF)
      model.sum <- summary(model.gam)
      gam.rsq <- model.sum$r.sq
      gam.age.p <- model.sum$s.table %>% as_tibble() %>% pull('p-value') # F and p of s(age)
      gam.age.F <- model.sum$s.table %>% as_tibble() %>% pull('F')
      gam.age.edf <- model.sum$s.table %>% as_tibble() %>% pull('edf') %>% round(digits = 0)
      
      # 3. Find out EOA (ESTIMATED Onset Age), i.e. since when RD is greater than 0 and increasing and inflection point
      model.gam.derv <- gratia::derivatives(model.gam, interval = "confidence", n = 1000) # from gratia. "confidence" for point-wise intervals, 95% CI
      age_span <- model.gam.derv %>% select(age)
      
      coef_gam <- coef(model.gam)
      second_derivatives <- predict(model.gam, age_span, type = "lpmatrix", deriv = 2) %>% 
        as_tibble() %>% 
        mutate(age = age_span$age) %>% 
        mutate(second_deriv = `s(age).1` * coef_gam["s(age).1"] + `s(age).2` * coef_gam["s(age).2"] + `s(age).3` * coef_gam["s(age).3"] ) # Calculate the Total Second Derivative
      
      pred_statistics <- predict(model.gam, age_span, type = 'response', se.fit = T) %>% 
        as_tibble() %>% 
        mutate(age = age_span$age) %>% 
        left_join(model.gam.derv %>% select(age, derivative = .derivative, derivative.lower = .lower_ci, derivative.upper = .upper_ci)) %>% 
        #left_join(model.gam.derv %>% select(age = data, derivative = derivative, derivative.lower = lower, derivative.upper = upper)) %>%
        select(age, derivative, derivative.lower, derivative.upper, statistics.pred = fit, statistics.pred.se = se.fit) %>% 
        mutate(statistics.lower = statistics.pred - 1.96 * statistics.pred.se) %>% # 95% CI
        mutate(statistics.upper = statistics.pred + 1.96 * statistics.pred.se)
      
      if(mk.tau > 0) {
        temp.EOA <- pred_statistics %>% filter(statistics.lower > 0 & derivative.lower > 0)
      } else {
        temp.EOA <- pred_statistics %>% filter(statistics.upper < 0 & derivative.upper < 0)
      }
      
      EOA <- ifelse(
        nrow(temp.EOA) == 0, NA_real_,
        min(temp.EOA$age)
      )
      
      # 4. Predicts statistics at multiple age points 
      statistics.pred.age <- predict(model.gam, tibble(age=seq(55, 85, by = 5)), type = 'response', se.fit = T) %>%
        as_tibble() %>%
        mutate(fit = as.numeric(fit), se.fit = as.numeric(se.fit)) %>% 
        mutate(age=seq(55, 85, by = 5)) 
      statistics.age <- statistics.pred.age %>% 
        select(age, value = fit) %>% 
        pivot_wider(names_prefix = 'statistics_at_age_', names_from = 'age', values_from = 'value')
      statistics.se.age <- statistics.pred.age %>% 
        select(age, value = se.fit) %>% 
        pivot_wider(names_prefix = 'statistics_se_at_age_', names_from = 'age', values_from = 'value')
      
      # 5. Run statistics ~ age by age quantiles
      B.quartile <- iDF %>% 
        mutate(quartile = ntile(age, 4)) %>%
        mutate(quartile = case_when(
          quartile == 1 ~ 'Q1',
          quartile == 2 ~ 'Q2',
          quartile == 3 ~ 'Q3',
          quartile == 4 ~ 'Q4',
          TRUE ~ NA_character_
        )) %>% 
        mutate(quartile = factor(quartile, levels = c('Q1', 'Q2', 'Q3', 'Q4'))) %>% 
        group_by(quartile) %>% 
        nest() %>% 
        mutate(map_dfr(data, function(df.quartile){
          lm(statistics~age, data=df.quartile) %>% 
            broom::tidy() %>% 
            filter(term == 'age') %>% 
            select(B_quartile = estimate, se_quartile = std.error, p_quartile = p.value)
        })) %>% select(-data)
      
      diff.B.quartile <- map2_dfr(
        list(list('Q4', 'Q1'), list('Q4', 'Q2'), list('Q4', 'Q3'), list('Q3', 'Q1'), list('Q3', 'Q2'), list('Q2', 'Q1')),
        list('Q4-Q1', 'Q4-Q2', 'Q4-Q3', 'Q3-Q1', 'Q3-Q2', 'Q2-Q1'),
        function(pair, label) {
          A <- pair[[1]]
          B <- pair[[2]]
          est_A <- B.quartile %>% filter(quartile == A) %>% pull(B_quartile)
          est_B <- B.quartile %>% filter(quartile == B) %>% pull(B_quartile)
          se_A <- B.quartile %>% filter(quartile == A) %>% pull(se_quartile)
          se_B <- B.quartile %>% filter(quartile == B) %>% pull(se_quartile)
          diff_B <- est_A - est_B
          z <- diff_B / sqrt((se_A)^2 + (se_B)^2)
          p <- 2 * (1 - pnorm(abs(z)))
          tibble(pair = label, diff_B, z, p)
        }
      )
      
      statistics.quartile.diff <- cbind(
        B.quartile %>% 
          pivot_longer(cols = -quartile, names_to = "variable", values_to = "value") %>%
          pivot_wider(names_from = c("quartile", "variable"), values_from = "value"),
        diff.B.quartile %>% 
          pivot_longer(cols = -pair, names_to = "variable", values_to = "value") %>%
          pivot_wider(names_from = c("pair", "variable"), values_from = "value")
      )
      
      
      # Combine the statistics
      tibble(
        atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, age_window = iAge_Window, PRS = iPRS_Type, statistics = iStatistics,
        mk.tau, mk.p, mk.z,
        gam.age.p, gam.age.F, gam.age.edf, gam.rsq,
        EOA
      ) %>% cbind(statistics.age, statistics.se.age, statistics.quartile.diff) 
    }
  )
  write_rds(result.GAM, str_c(path.result, '/GAM_fitting_raw.rds'))
}

# PCA analysis on volumetric IDPs------------------------------------------------------------
demo.uc <- read_csv('data/demo.csv') %>% 
  filter(ethnicity == 'UKB-UC')
path.pca <- create_when_absent('results/600-PCA')

idp.input <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds') %>% 
  filter(PRS == 'PRS-UKB') %>%
  filter(lm_p < 0.05) %>% 
  filter(measure %in% c("volume_fast_normalised", "volume_first_normalised")) %>% 
  filter(index <= 96) %>% 
  select(atlas, weight, index, measure) %>%
  left_join(read_rds('data/tidy_imaging_cognition.rds'), by = c('atlas', 'weight', 'index', 'measure')) %>% 
  mutate(description = description %>%
           str_remove("^Volume of grey matter in ") %>%  # Remove prefix
           str_remove("^Volume of ") %>%  # Remove prefix
           str_remove(" \\(normalised\\)$")) %>%            # Remove suffix
  select(eid, description, value) %>% 
  pivot_wider(names_from = description, values_from = value) %>% 
  filter(eid %in% demo.uc$eid)

# Parallel analysis to determine PC nums
pc.parallel <- psych::fa.parallel(idp.input %>% select(-eid), fa='pc', n.iter=10000, show.legend=FALSE)
df.parallel <- 
  tibble(
    x = 1:length(pc.parallel$fa.values),
    `Principle values` = pc.parallel$pc.values,
    `Simulated data` = pc.parallel$pc.sim,  
    `Resampled data` = pc.parallel$pc.simr) %>% 
  pivot_longer(cols = `Principle values`:`Resampled data`,
               values_to = 'pc',
               names_to = 'type')
plot.parallel <- ggplot(df.parallel, aes(x, pc)) + 
  geom_point(aes(colour = type, shape = type), size = 1) + 
  geom_line(aes(colour = type)) + 
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = as.numeric(pc.parallel$ncomp), linetype = 'dashed', color = 'blue') + 
  labs(
    title = str_c('Parallel analysis suggests ', pc.parallel$ncomp, ' components'),
    x = 'Principle component numbers',
    y = 'Eigen values'
  )  + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks = element_line(colour = 'black', size = 0.3),
    axis.title = element_text(size = 8, color = 'black'),
    axis.text = element_text(size = 8, color = 'black')
  ) + # Mediate the title
  ggsci::scale_color_bmj()
plot.parallel # 5 5

# Do pca and saves the scores and weights
NCOMP = pc.parallel$ncomp
pca.result <- psych::principal(idp.input %>% select(-eid), nfactors=NCOMP, rotate='promax', scores=TRUE)
saveRDS(pca.result, glue::glue('{path.pca}/pca_GM_norm.rds'))
# promax > geominT > varimax
# Visualize the weights
pca.weights <- 
  as.matrix(pca.result$weights) %>% 
  reshape::melt.matrix() %>% 
  set_names(c("x", "y", "value")) %>% 
  mutate(y = factor(y, levels = c(str_c('RC', 1:NCOMP)))) %>% 
  rename(region = x) %>% 
  mutate(
    y = case_when(
      y == 'RC1' ~ 'subcortical',
      y == 'RC2' ~ 'vmPFC', # A region combining the frontal orbital and subcallosal areas can be referred to as the ventromedial prefrontal cortex 
      y == 'RC3' ~ 'paraHipp-A',
      y == 'RC4' ~ 'frontal',
      y == 'RC5' ~ 'paraHipp-P',
      y == 'RC6' ~ 'IPL', # Both AG and SMG together form the IPL, which plays a key role in language, spatial awareness, and numerical cognition
      y == 'RC7' ~ 'insular'
    )
  )

ggplot(pca.weights, aes(x = y, y = region, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = 'white') +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(angle = 45, color = "black", vjust = 0, size = 8),
    axis.text.x = element_text(angle = 45, color = "black", size = 8),
    axis.line = element_line(colour = "black", size = 0.3),
    axis.ticks = element_line(colour = "black", size = 0.3)
  ) 
plot.pca.loading
ggsave(str_c(path.result, '/loading.pdf'), plot.pca.loading, width=4.7, height=8.9)

data.pca <- 
  cbind(eid=idp.input$eid, pca.result$scores) %>% 
  as_tibble() %>% 
  pivot_longer(-eid, names_to = 'label_PCA', values_to = 'value_PCA') %>% 
  mutate(
    label_PCA = case_when(
      label_PCA == 'RC1' ~ 'subcortical',
      label_PCA == 'RC2' ~ 'vmPFC',
      label_PCA == 'RC3' ~ 'paraHipp-A',
      label_PCA == 'RC4' ~ 'frontal',
      label_PCA == 'RC5' ~ 'paraHipp-P',
      label_PCA == 'RC6' ~ 'IPL',
      label_PCA == 'RC7' ~ 'insular'
    )
  ) %>% 
  mutate(component_type = 'GM_normalised')

# Make sure they all decline with age
pca.qc <- data.pca %>% left_join(demo.uc %>% select(eid, age))

ggplot(pca.qc, aes(x = age, y = value_PCA, color = label_PCA)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~label_PCA, scales = 'free') 
write_csv(data.pca, glue::glue('{path.pca}/pca_score_GM.csv'))

# PCA analysis on diffusion IDPs------------------------------------------------------------
demo.uc <- read_csv('data/demo.csv') %>% 
  filter(ethnicity == 'UKB-UC')
path.pca <- create_when_absent('results/600-PCA')
DESCIRPTION_MEAN_TRACT <- DESCRIPTION_ALL_PHE %>% 
  filter(measure == 'mean_fa') %>% select(index, acronym, hemi) %>% 
  mutate(hemi = ifelse(
    hemi %in% c('L', 'R'),
    str_c('-', hemi),
    ''
  )) %>% 
  mutate(tract = str_c(acronym, hemi)) %>% 
  arrange(index) %>% 
  pull(tract)

idp.input <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds') %>% 
  filter(PRS == 'PRS-UKB') %>%
  filter(lm_p < 0.05) %>% 
  filter(measure %in% c("mean_fa", "mean_md", "mean_mo", "mean_icvf", "mean_od", "mean_isovf")) %>%
  select(atlas, weight, index, measure) %>%
  left_join(read_rds('data/tidy_imaging_cognition.rds') %>% filter(eid %in% demo.uc$eid), by = c('atlas', 'weight', 'index', 'measure')) %>% 
  mutate(tract = map_chr(index, ~DESCIRPTION_MEAN_TRACT[[.]])) %>% 
  mutate(measure = toupper(str_replace(measure, "\\mean_", ""))) %>% 
  mutate(measure = factor(measure, levels=c("FA", "MD", "MO", "ICVF", "ISOVF", "OD"))) %>% 
  mutate(description = str_c(tract, '-', measure)) %>%
  select(eid, description, value) %>% 
  pivot_wider(names_from = description, values_from = value)

# Parallel analysis to determine PC nums
pc.parallel <- psych::fa.parallel(idp.input %>% select(-eid), fa='pc', n.iter=10000, show.legend=FALSE)
df.parallel <- 
  tibble(
    x = 1:length(pc.parallel$fa.values),
    `Principle values` = pc.parallel$pc.values,
    `Simulated data` = pc.parallel$pc.sim,  
    `Resampled data` = pc.parallel$pc.simr) %>% 
  pivot_longer(cols = `Principle values`:`Resampled data`,
               values_to = 'pc',
               names_to = 'type')
ggplot(df.parallel, aes(x, pc)) + 
  geom_point(aes(colour = type, shape = type), size = 1) + 
  geom_line(aes(colour = type)) + 
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = as.numeric(pc.parallel$ncomp), linetype = 'dashed', color = 'blue') + 
  labs(
    title = str_c('Parallel analysis suggests ', pc.parallel$ncomp, ' components'),
    x = 'Principle component numbers',
    y = 'Eigen values'
  )  + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks = element_line(colour = 'black', size = 0.3),
    axis.title = element_text(size = 8, color = 'black'),
    axis.text = element_text(size = 8, color = 'black')
  ) + # Mediate the title
  ggsci::scale_color_bmj()


# Do pca and saves the scores and weights
NCOMP = pc.parallel$ncomp
pca.result <- psych::principal(idp.input %>% select(-eid), nfactors=NCOMP, rotate='promax', scores=TRUE)
saveRDS(pca.result, glue::glue('{path.pca}/pca_GM_norm.rds'))
pca.weights <- 
  as.matrix(pca.result$weights) %>% 
  reshape::melt.matrix() %>% 
  set_names(c("x", "y", "value")) %>% 
  mutate(y = factor(y, levels = c(str_c('RC', 1:NCOMP)))) %>% 
  rename(region = x) %>% 
  mutate(
    y = case_when(
      y == 'RC1' ~ 'subcortical',
      y == 'RC2' ~ 'vmPFC', # A region combining the frontal orbital and subcallosal areas can be referred to as the ventromedial prefrontal cortex 
      y == 'RC3' ~ 'paraHipp-A',
      y == 'RC4' ~ 'frontal',
      y == 'RC5' ~ 'paraHipp-P',
      y == 'RC6' ~ 'IPL', # Both AG and SMG together form the IPL, which plays a key role in language, spatial awareness, and numerical cognition
      y == 'RC7' ~ 'insular'
    )
  )

ggplot(pca.weights, aes(x = y, y = region, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = 'white') +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(angle = 45, color = "black", vjust = 0, size = 4),
    axis.text.x = element_text(angle = 45, color = "black", size = 4),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1)
  ) 
plot.pca.loading
ggsave(str_c(path.result, '/loading.pdf'), plot.pca.loading, width=4.7, height=8.9)

data.pca <- 
  cbind(eid=idp.input$eid, pca.result$scores) %>% 
  as_tibble() %>% 
  pivot_longer(-eid, names_to = 'label_PCA', values_to = 'value_PCA') %>% 
  mutate(component_type = 'diffusion_measures')

# Make sure they all decline with age
pca.qc <- data.pca %>% left_join(demo.uc %>% select(eid, age))

ggplot(pca.qc, aes(x = age, y = value_PCA, color = label_PCA)) +
  #geom_point() +
  geom_smooth() +
  facet_wrap(~label_PCA, scales = 'free') 
data.pca <- data.pca %>% 
  mutate(value_PCA = case_when(
    label_PCA %in% c('RC2', 'RC3', 'RC5', 'RC10', 'RC11', 'RC12', 'RC14', 'RC15') ~ -value_PCA,
    TRUE ~ value_PCA
  ))
 
write_csv(data.pca, glue::glue('{path.pca}/pca_score_diffusion.csv'))

# Mediation analysis --------------

# Combine PCA scores
PCA.score <- map_dfr(
  c(
    # List of PCA scores generated earlier
  ),
  ~read_csv(.)
)

write_rds(PCA.score, '') # Customised dir

# Cognition data
cognition <- read_rds('data/tidy_imaging_cognition.rds') %>% 
  filter(measure == 'cognition') %>% 
  mutate(eid, cognition_label = description, cognition = value, .keep = 'none')

demo <- read_csv('demo_combined.csv') %>% 
  mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
  mutate(age_group = case_when(
    age < 55 ~ "<55",
    age >= 55 & age < 60 ~ "55-60",
    age >= 60 & age < 65 ~ "60-65",
    age >= 65 & age < 70 ~ "65-70",
    age >= 70 & age < 75 ~ "70-75",
    age >= 75 ~ ">75",
    TRUE ~ NA_character_
  ), .after = age) %>% 
  pivot_longer(c('PRS-UKB', 'PRS-APin', 'PRS-APex'), names_to = 'PRS.type', values_to = 'PRS') %>% 
  mutate(age_group = factor(age_group, levels = c("<55", "55-60", "60-65", "65-70", "70-75", ">75"), ordered = TRUE)) 
PCA.score <- read_rds('PCA_score.rds')

data.nest <- demo %>% left_join(PCA.score, by = 'eid') %>% 
  left_join(read_csv('cognition.csv')) %>% 
  na.omit() %>% 
  group_by(PRS.type, age_group, component_type, label_PCA, cognition_label) %>% 
  nest()

NUM_BOOT_MEDIATION = 10000 # Customized boot num
future::plan(future::multisession, workers = 90)
mediation.statistics <- furrr::future_pmap_dfr(
  list(data.nest$component_type, data.nest$label_PCA, data.nest$PRS.type, data.nest$age_group, data.nest$cognition_label, data.nest$data),
  function(iComponent_type, iLabel_PCA, iPRS, iAge_GROUP, iLabel_cognition, iDF) {
    # Scale the variables to get standardised coefficients
    df_scaled <- iDF %>% 
      mutate(across(c(PRS, cognition, value_PCA, age, education, BMI, townsend, date.scanning, genePC1:genePC40), ~scale(.)[,]))
    # Fit the model for mediator and X
    model.mediator <- lm(value_PCA ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40, data = df_scaled)
    coef.mediator <- broom::tidy(model.mediator)
    a <- coef.mediator %>% filter(term == 'PRS') %>% pull('estimate')
    a.se <- coef.mediator %>% filter(term == 'PRS') %>% pull('std.error')
    a.pval <- coef.mediator %>% filter(term == 'PRS') %>% pull('p.value')
    
    # Fit regression model for the total effect (Y, X and M)
    model.total <- lm(cognition ~ PRS + value_PCA + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40, data = df_scaled)
    coef.total <- broom::tidy(model.total)
    b <- coef.total %>% filter(term == 'value_PCA') %>% pull('estimate')
    b.se <-  coef.total %>% filter(term == 'value_PCA') %>% pull('std.error')
    b.pval <-  coef.total %>% filter(term == 'value_PCA') %>% pull('p.value')
    
    # Conduct mediation analysis
    mediate.result <- mediation::mediate(model.mediator, model.total,  treat = 'PRS', mediator = 'value_PCA', sims = NUM_BOOT_MEDIATION, boot = T)
    mediate.sum <- summary(mediate.result)
    
    # ACME stands for average causal mediation effects. This is the indirect effect of the X on the Y that goes through the M.
    acme <- mediate.sum$d0
    acme.pval <- mediate.sum$d0.p
    acme.low <- mediate.sum$d0.ci[[1]]
    acme.up <- mediate.sum$d0.ci[[2]]
    
    # ADE stands for average direct effects. It describes the direct effect of the IV on the DV.
    ade <- mediate.sum$z0
    ade.pval <- mediate.sum$z0.p
    ade.low <- mediate.sum$z0.ci[[1]]
    ade.up <- mediate.sum$z0.ci[[2]]
    
    # Total Effect stands for the total effect (direct + indirect) of the IV onto the DV.
    total <- mediate.sum$tau.coef
    total.pval <- mediate.sum$tau.p
    total.low <- mediate.sum$tau.ci[[1]]
    total.up <- mediate.sum$tau.ci[[2]]
    
    # Prop. Mediated describes the proportion of the effect of the IV on the DV that goes through the mediator.
    
    # Save the result of mediation
    return(tibble(
      age_group = iAge_GROUP, PRS = iPRS, PCA = iLabel_PCA, PCA_component = iComponent_type, cognition = iLabel_cognition,
      a, a.se, a.pval,
      b, b.se, b.pval,
      acme, acme.pval, acme.low, acme.up,
      ade, ade.pval, ade.low, ade.up,
      total, total.pval, total.low, total.up
    ))
  }
)


# Sensitivity - IA-beta and IARD By PRS  ------------------------------------------

df.vis <-  read_rds('results/402-SW-analysis/intra_age_statistics.rds') %>%
  filter(
    (measure == 'volume_first_normalised' & index %in% c(9) ) |
      (measure == 'mean_fa' & index %in% c(6) ) |
      (measure == 'volume_whole_brain_normalised' & index %in% c(-102) )
  )  %>% 
  select(index, measure, age, beta = lm_beta, se = lm_se, PRS) %>% 
  na.omit() %>% 
  mutate(region = map2_chr(index, measure, ~DESCRIPTION_ALL_PHE %>% filter(index == .x & measure == .y) %>% pull(description)) ) %>% 
  mutate(region = region %>%
           str_remove("^Volume of grey matter in ") %>%  # Remove prefix
           str_remove("^Volume of ") %>%  # Remove prefix
           str_remove(" \\(normalised\\)$")) %>%            # Remove suffix
  mutate(region = ifelse(region == 'Mean FA in fornix on FA skeleton', 'Fornix-FA', region))

p1 <- ggplot(
  data = df.vis %>% 
    mutate(region = factor(region, levels = c('hippocampus (left)', 'ventricular CSFV', 'Fornix-FA'))), 
  mapping = aes(x = age, y = beta, color = region, fill = region)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_ribbon(aes(ymin = beta-se, ymax=beta+se), alpha = 0.1, color = NA, position = "identity") +
  geom_smooth(method = 'gam', formula = y~s(x, k=4, bs = 'ps'), linewidth = 0.5, alpha = 0.2) +
  theme_classic(base_size = 10) + 
  scale_x_continuous(breaks = seq(50, 80, by = 5)) +
  labs(
    x = 'Age (years)',
    y = 'Beta'
  ) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  ggsci::scale_color_aaas() +
  ggsci::scale_fill_aaas() +
  theme(
    strip.placement = "outside",
    strip.clip = "off",
    strip.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    title = element_text(size = 8, color = 'black'),
    axis.text = element_text(size = 8, color = 'black')) +
  facet_wrap(PRS~region, scales = 'free') +
  ggeasy::easy_remove_legend()


df.vis <-  read_rds('results/402-SW-analysis/intra_age_statistics.rds') %>%
  filter(
    (measure == 'volume_first_normalised' & index %in% c(9) ) |
      #(measure == 'volume_fast_normalised' & index %in% c(67, 68, 69, 70) ) |
      (measure == 'mean_fa' & index %in% c(6) ) |
      (measure == 'volume_whole_brain_normalised' & index %in% c(-102) )
  )  %>% 
  select(index, measure, age, HL = HL_relative_diff, ML = ML_relative_diff, HM = HM_relative_diff, PRS) %>% 
  na.omit() %>% 
  mutate(region = map2_chr(index, measure, ~DESCRIPTION_ALL_PHE %>% filter(index == .x & measure == .y) %>% pull(description)) ) %>% 
  mutate(region = region %>%
           str_remove("^Volume of grey matter in ") %>%  # Remove prefix
           str_remove("^Volume of ") %>%  # Remove prefix
           str_remove(" \\(normalised\\)$")) %>%            # Remove suffix
  mutate(region = ifelse(region == 'Mean FA in fornix on FA skeleton', 'Fornix-FA', region)) %>% 
  pivot_longer(HL:HM, names_to = 'Risk pair', values_to = 'IARD') %>% 
  mutate(`Risk pair` = case_when(
    `Risk pair` == 'HL' ~ 'high to low',
    `Risk pair` == 'ML' ~ 'medium to low',
    `Risk pair` == 'HM' ~ 'high to medium',
    TRUE ~ NA_character_
  ))


p2 <- ggplot(
  data = df.vis %>%
    filter(PRS == 'PRS-UKB') %>% 
    mutate(region = factor(region, levels = c('hippocampus (left)', 'ventricular CSFV', 'Fornix-FA'))), 
  mapping = aes(x = age, y = IARD, color = `Risk pair`, fill = `Risk pair`)) +
  #geom_point(size = 0.1, alpha = 0.5) +
  #geom_ribbon(aes(ymin = beta-se, ymax=beta+se), alpha = 0.1, color = NA, position = "identity") +
  geom_smooth(method = 'gam', formula = y~s(x, k=4, bs = 'ps'), linewidth = 0.5, alpha = 0.2) +
  theme_classic(base_size = 10) + 
  scale_x_continuous(breaks = seq(50, 80, by = 5)) +
  labs(
    x = 'Age (years)',
    y = 'IARD'
  ) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg() +
  theme(
    strip.placement = "outside",
    strip.clip = "off",
    strip.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    title = element_text(size = 8, color = 'black'),
    axis.text = element_text(size = 8, color = 'black')) +
  facet_wrap(~region, scales = 'free') +
  ggeasy::easy_remove_legend()

ggpubr::ggarrange(p1, p2, ncol = 1, heights = c(3,1))
# Sensitivity analysis - APex * APOE4 -------------------------------------
data.demo.by.apoe <- read_csv('data/demo.csv') %>% 
  filter(ethnicity=='UKB-UC') %>% na.omit() %>% 
  mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
  select(-`PRS-UKB`, -`PRS-APin`) %>% 
  rename(PRS = `PRS-APex`) %>% 
  mutate(carrier = case_when(
    apoe %in% c('e3e4', 'e4e4', 'e2e4') ~ 'Carrier',
    TRUE ~ 'Non-Carrier'
  )) %>% 
  mutate(apoe = case_when(
    apoe %in% c('e3e4', 'e4e4', 'e2e4') ~ 'Carrier',
    TRUE ~ 'Non-Carrier'
  )) %>% 
  mutate(age_group = case_when(
    age < 55 ~ "<55",
    age >= 55 & age < 60 ~ "55-60",
    age >= 60 & age < 65 ~ "60-65",
    age >= 65 & age < 70 ~ "65-70",
    age >= 70 & age < 75 ~ "70-75",
    age >= 75 ~ ">75",
    TRUE ~ NA_character_
  ), .after = age) %>% 
  mutate(age_group = factor(age_group, levels = c("<55", "55-60", "60-65", "65-70", "70-75", ">75"), ordered = TRUE)) 

data.demo.combined <- data.demo.by.apoe %>% mutate(carrier = 'Combined')
data.demo <- rbind(data.demo.by.apoe, data.demo.combined) %>% 
  mutate(apoe = factor(apoe, levels = c('Non-Carrier', 'Carrier')))

data.idp <- read_rds('data/tidy_imaging_cognition.rds') %>% 
  filter(eid %in% data.demo$eid) %>% 
  group_by(atlas, weight, index, measure) %>% 
  nest()

future::plan(future::multisession, workers = 8)
lm.result <- furrr::future_pmap_dfr(
  list(
    data.idp$atlas, data.idp$weight, data.idp$index, data.idp$measure, data.idp$data
  ), 
  function(iAtlas, iWeight, iIndex, iMeasure, iData){
    iData <- iData %>% left_join(data.demo) %>% group_by(carrier, age_group) %>% nest()
    result <- pmap_dfr(
      list(iData$carrier, iData$age_group, iData$data),
      function(iCarrier, iAge_group, iDf) {
        df.std <- iDf %>% 
          na.omit() %>% 
          mutate(across(c(value, age, education, BMI, townsend, date.scanning, genePC1:genePC40, PRS), ~scale(.)[,]))
        # Fit a linear model with interaction for combined (carrier + non-carrier)
        if(iCarrier == 'Combined') {
          # Try fitting the linear model and handling errors
          fit.model <- tryCatch(
            lm(df.std, formula = value ~ PRS*apoe + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
            error = function(e) NULL
          )
          # If lm() fails, return missing values for linear model fitting
          if (is.null(fit.model)) {
            #print('A')
            return(
              tibble(
                atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
                beta_PRS = NA_real_,
                p_PRS = NA_real_,
                se_PRS = NA_real_,
                beta_interaction = NA_real_,
                p_interaction = NA_real_,
                se_interaction = NA_real_
              ))
          } else {
            result.lm <- broom::tidy(fit.model) %>% 
              filter(term %in% c('PRS', 'PRS:apoeCarrier')) 
            return(
              tibble(
                atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
                beta_PRS = result.lm %>% filter(term == 'PRS') %>% pull('estimate'),
                p_PRS = result.lm %>% filter(term == 'PRS') %>% pull('p.value'),
                se_PRS = result.lm %>% filter(term == 'PRS') %>% pull('std.error'),
                beta_interaction = result.lm %>% filter(term == 'PRS:apoeCarrier') %>% pull('estimate'),
                p_interaction = result.lm %>% filter(term == 'PRS:apoeCarrier') %>% pull('p.value'),
                se_interaction = result.lm %>% filter(term == 'PRS:apoeCarrier') %>% pull('std.error')
              ))
          }
        } else {
          # Try fitting the linear model and handling errors
          fit.model <- tryCatch(
            lm(df.std, formula = value ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
            error = function(e) NULL
          )
          # If lm() fails, return missing values for linear model fitting
          if (is.null(fit.model)) {
            return(
              tibble(
                atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
                beta_PRS = NA_real_,
                p_PRS = NA_real_,
                se_PRS = NA_real_,
                beta_interaction = NA_real_,
                p_interaction = NA_real_,
                se_interaction = NA_real_
              ))
          } else {
            result.lm <- broom::tidy(fit.model) %>% 
              filter(term %in% c('PRS'))
            return(
              tibble(
                atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
                beta_PRS = result.lm %>% filter(term == 'PRS') %>% pull('estimate'),
                p_PRS = result.lm %>% filter(term == 'PRS') %>% pull('p.value'),
                se_PRS = result.lm %>% filter(term == 'PRS') %>% pull('std.error'),
                beta_interaction = NA_real_,
                p_interaction = NA_real_,
                se_interaction = NA_real_
              ))
          }
        }
      }
    )
  })

result.sensitivity <- lm.result %>% mutate(region = map2_chr(index, measure, ~DESCRIPTION_ALL_PHE %>% filter(index == .x & measure == .y) %>% pull(description)) )
write_rds(result.sensitivity, '') # Customised output dir

# Sensitivity analysis - Race ---------------------------------------------
data.demo <- read_csv('data/demo.csv') %>% 
  mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
  mutate(age_group = case_when(
    age < 55 ~ "<55",
    age >= 55 & age < 60 ~ "55-60",
    age >= 60 & age < 65 ~ "60-65",
    age >= 65 & age < 70 ~ "65-70",
    age >= 70 & age < 75 ~ "70-75",
    age >= 75 ~ ">75",
    TRUE ~ NA_character_
  ), .after = age) %>% 
  mutate(age_group = factor(age_group, levels = c("<55", "55-60", "60-65", "65-70", "70-75", ">75"), ordered = TRUE)) %>% 
  pivot_longer(c('PRS-UKB', 'PRS-APin', 'PRS-APex'), names_to = 'PRS.type', values_to = 'PRS') 

# 109 IDPs
data.idp <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds') %>% 
  filter(lm_p < 0.05) %>% 
  filter(PRS == 'PRS-UKB') %>% 
  filter(measure %in% c(
    'volume_fast_normalised', 'volume_first_normalised', 'volume_whole_brain_normalised',
    'mean_fa', 'mean_md', 'mean_mo', 'mean_icvf', 'mean_od', 'mean_isovf'
  )) %>% 
  select(atlas, weight, index, measure, description) %>% 
  distinct() %>%
  filter(index <= 96) %>% 
  left_join(read_rds('data/tidy_imaging_cognition.rds')) %>% 
  left_join(data.demo) %>% 
  group_by(atlas, weight, index, measure, ethnicity, age_group, PRS.type) %>% 
  nest()

future::plan(future::multisession, workers = 8)
lm.result <- furrr::future_pmap_dfr(
  list(
    data.idp$atlas, data.idp$weight, data.idp$index, data.idp$measure, data.idp$ethnicity, data.idp$age_group, data.idp$PRS.type, data.idp$data
  ), 
  function(iAtlas, iWeight, iIndex, iMeasure, iEthnicity, iAge_group, iPRS, iData){
    df.std <- iData %>% 
      na.omit() %>% 
      mutate(across(c(value, age, education, BMI, townsend, date.scanning, genePC1:genePC40, PRS), ~scale(.)[,]))
    # Try fitting the linear model and handling errors
    fit.model <- tryCatch(
      lm(df.std, formula = value ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
      error = function(e) NULL
    )
    # If lm() fails, return missing values for linear model fitting
    if (is.null(fit.model)) {
      return(
        tibble(
          atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, age_group = iAge_group, ethnicity = iEthnicity, PRS = iPRS,
          beta_PRS = NA_real_,
          p_PRS = NA_real_,
          se_PRS = NA_real_
        ))
    } else {
      result.lm <- broom::tidy(fit.model) %>% 
        filter(term %in% c('PRS'))
      return(
        tibble(
          atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, age_group = iAge_group, ethnicity = iEthnicity, PRS = iPRS,
          beta_PRS = result.lm %>% filter(term == 'PRS') %>% pull('estimate'),
          p_PRS = result.lm %>% filter(term == 'PRS') %>% pull('p.value'),
          se_PRS = result.lm %>% filter(term == 'PRS') %>% pull('std.error')
        ))
      }
  })

result.lm <- lm.result %>% mutate(region = map2_chr(index, measure, ~DESCRIPTION_ALL_PHE %>% filter(index == .x & measure == .y) %>% pull(description)) )
write_rds(result.lm, '') # Customised output dir

# Sensitivity - Other PRS --------------------------------------------------------------
path.result <- create_when_absent('results/502-PRS-interaction')
demo.uc <- read_csv('data/demo.csv') %>% 
  filter(ethnicity=='UKB-UC') %>% na.omit() %>% 
  mutate(sex = ifelse(sex == 'M', 1, 0)) %>% 
  select(eid, PRS.AD = `PRS-UKB`, age:genePC40) %>% 
  mutate(age_group = case_when(
    age < 55 ~ "<55",
    age >= 55 & age < 60 ~ "55-60",
    age >= 60 & age < 65 ~ "60-65",
    age >= 65 & age < 70 ~ "65-70",
    age >= 70 & age < 75 ~ "70-75",
    age >= 75 ~ ">75",
    TRUE ~ NA_character_
  ), .after = age) %>% 
  mutate(age_group = factor(age_group, levels = c("<55", "55-60", "60-65", "65-70", "70-75", ">75"), ordered = TRUE))

data.base <- read_csv('data/PRS_all_disease.csv') %>% 
  filter(eid %in% demo.uc$eid) %>% 
  select(eid, asthma:venous_thromboembolic_disease) %>% 
  pivot_longer(asthma:venous_thromboembolic_disease, names_to = 'disease', values_to = 'PRS.disease') %>% 
  left_join(demo.uc)

data.idp.ROI.nest <- read_rds('results/201-IDP-PRS-after-75/result_report_all.rds') %>% 
  filter(PRS == 'PRS-UKB') %>% 
  filter(lm_p < 0.05) %>% 
  select(atlas, weight, index, measure) %>% 
  left_join((
    read_rds('data/tidy_imaging_cognition.rds') %>% 
      filter(eid %in% demo.uc$eid)
  ), by = c('atlas', 'weight', 'index', 'measure')) %>% 
  na.omit() %>% 
  group_by(atlas, weight, index, measure) %>% nest()

future::plan(future::multisession, workers = 10)
result.between.PRS <- furrr::future_pmap_dfr(
  list(data.idp.ROI.nest$atlas, data.idp.ROI.nest$weight, data.idp.ROI.nest$index, data.idp.ROI.nest$measure, data.idp.ROI.nest$data),
  function(iAtlas, iWeight, iIndex, iMeasure, iDF) {
    iDF <- left_join(iDF, data.base, by = 'eid') %>% group_by(age_group, disease) %>% nest()
    pmap_dfr(
      list(iDF$disease, iDF$age_group, iDF$data),
      function(iDisease, iAgeGroup, df) {
        #df <- iDF$data[[1]]
        df.std <- df %>% 
          na.omit() %>% 
          mutate(across(c(value, age, education, BMI, townsend, date.scanning, genePC1:genePC40, PRS.AD, PRS.disease), ~scale(.)[,]))
        fit.model <- lm(df.std, formula = value ~ PRS.AD * PRS.disease + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40)
        result.lm <- broom::tidy(fit.model) %>% 
          filter(term %in% c('PRS.AD', 'PRS.disease', 'PRS.AD:PRS.disease'))
        cbind(
          atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, disease = iDisease, age_group = iAgeGroup,
          result.lm
        )
      }
    )
  }
)

write_rds(result.between.PRS, str_c('')) # Customised output dir

# Sensitivity - AD-PRS with other thresh ----------------------------------------
path.result <- create_when_absent('results/800-sensitive/')
demo <- read_csv('data/demo.csv') %>% 
  select(-(`PRS-UKB`:`PRS-APex`)) %>% 
  filter(ethnicity == 'UKB-UC')

PRS.CT <- read_csv('data/raw/PRS_DXY_54386_PLINK.csv') %>% 
  filter(eid %in% demo$eid) %>% 
  pivot_longer(-eid, names_to = 'PRS.type', values_to = 'PRS')

demo.PRS <- demo %>% 
  left_join(PRS.CT, by = 'eid') %>% 
  mutate(
    age_group = case_when(
      age < 55 ~ "<55",
      age >= 55 & age < 60 ~ "55-60",
      age >= 60 & age < 65 ~ "60-65",
      age >= 65 & age < 70 ~ "65-70",
      age >= 70 & age < 75 ~ "70-75",
      age >= 75 ~ ">75",
      TRUE ~ NA_character_
    ), .after = age
  ) %>% 
  mutate(age_group = factor(age_group, levels = c("<55", "55-60", "60-65", "65-70", "70-75", ">75"), ordered = TRUE))


data.idp <- read_rds('data/tidy_imaging_cognition.rds') %>% 
  filter(eid %in% demo$eid) %>% 
  group_by(atlas, weight, index, measure) %>% 
  nest()

# Run linear regression models within age groups and PRSs
future::plan(future::multisession, workers = 8)
lm.result <- furrr::future_pmap_dfr(
  list(
    data.idp$atlas, data.idp$weight, data.idp$index, data.idp$measure, data.idp$data
  ), 
  function(iAtlas, iWeight, iIndex, iMeasure, iData){
    #iData <- data.idp$data[[1]]
    iData <- iData %>% left_join(demo.PRS) %>% group_by(PRS.type, age_group) %>% nest()
    result <- pmap_dfr(
      list(iData$PRS.type, iData$age_group, iData$data),
      function(iPRS_type, iAge_group, iDf) {
        df.std <- iDf %>% 
          na.omit() %>% 
          mutate(across(c(value, age, education, BMI, townsend, date.scanning, genePC1:genePC40, PRS), ~scale(.)[,]))
        # Try fitting the linear model and handling errors
        fit.model <- tryCatch(
          lm(df.std, formula = value ~ PRS + age + sex + education + townsend + BMI + hypertensive + diabetic + date.scanning + center.newcastle + center.cheadle + center.reading + genePC1 + genePC2 + genePC3 + genePC4 + genePC5 + genePC6 + genePC7 + genePC8 + genePC9 + genePC10 + genePC11 + genePC12 + genePC13 + genePC14 + genePC15 + genePC16 + genePC17 + genePC18 + genePC19 + genePC20 + genePC21 + genePC22 + genePC23 + genePC24 + genePC25 + genePC26 + genePC27 + genePC28 + genePC29 + genePC30 + genePC31 + genePC32 + genePC33 + genePC34 + genePC35 + genePC36 + genePC37 + genePC38 + genePC39 + genePC40),
          error = function(e) NULL
        )
        # If lm() fails, return missing values for linear model fitting
        if (is.null(fit.model)) {
          return(
            tibble(
              atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
              beta_PRS = NA_real_,
              p_PRS = NA_real_,
              se_PRS = NA_real_,
              beta_interaction = NA_real_,
              p_interaction = NA_real_,
              se_interaction = NA_real_
            ))
        } else {
          result.lm <- broom::tidy(fit.model) %>% 
            filter(term %in% c('PRS'))
          return(
            tibble(
              atlas = iAtlas, weight = iWeight, index = iIndex, measure = iMeasure, carrier = iCarrier, age_group = iAge_group,
              beta_PRS = result.lm %>% filter(term == 'PRS') %>% pull('estimate'),
              p_PRS = result.lm %>% filter(term == 'PRS') %>% pull('p.value'),
              se_PRS = result.lm %>% filter(term == 'PRS') %>% pull('std.error'),
              beta_interaction = NA_real_,
              p_interaction = NA_real_,
              se_interaction = NA_real_
            ))
        }
      }
    )
  })

write_rds(lm.result, str_c('')) # Customized output dir


