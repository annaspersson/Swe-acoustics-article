# Constants ------------------------------------------------------------------
theme_set(
  theme_bw(base_size = 11) + 
    theme(
      panel.grid = element_blank()))

levels.vowel.IPA = c("[i]", "[ɪ]", "[ɛ]", "[æ]", "[ʌ]", "[ʊ]", "[u]", "[ɑ]")
levels.vowel.Arpabet = c("iy1", "ih1", "eh1", "ae1", "ah1", "uh1", "uw1", "aa1")
levels.vowel.SweFA = c("ii1", "ih1", "yy1", "yh1", "uu1", "uh1", "ee1", "eh1", "ae1", "{:", "{", "oe1", "oeh1", "9:", "9", "aa1", "ah1", "oa1", "oah1", "oo1", "oh1")
levels.vowel.IPA.swe = c("[iː]", "[ɪ]", "[yː]", "[ʏ]", "[ʉː]", "[ɵ]", "[eː]", "[ɛ]", "[ɛː]", "[æː]", "[æ]", "[øː]", "[ø]", "[œː]", "[œ]", "[ɑː]", "[a]", "[oː]", "[ɔ]", "[uː]", "[ʊ]")
levels.vowel.IPA.swe.long = c("[iː]", "[yː]", "[ʉː]", "[eː]", "[ɛː]", "[æː]", "[øː]", "[œː]", "[ɑː]", "[oː]", "[uː]")
levels.vowel.IPA.swe.long.noallo = c("[iː]", "[yː]", "[ʉː]", "[eː]", "[ɛː]", "[øː]", "[ɑː]", "[oː]", "[uː]")
levels.vowel.IPA.swe.short = c("[ɪ]", "[ʏ]", "[ɵ]", "[ɛ]", "[æ]", "[ø]", "[œ]", "[a]", "[ɔ]", "[ʊ]")
levels.vowel.IPA.swe.short.noallo = c("[ɪ]", "[ʏ]", "[ɵ]", "[ɛ]", "[ø]", "[a]", "[ɔ]", "[ʊ]")
levels.vowel.SweFA.wrongEncod = c("ii1", "ih1", "yy1", "yh1", "uu1", "uh1", "ee1", "eh1", "ae1", "{:", "{", "oe1", "oeh1", "09:00", "9", "aa1", "ah1", "oa1", "oah1", "oo1", "oh1")
levels.vowel.Swe.word = c("hid", "hidd", "hyd", "hydd", "hud", "hudd", "hed", "hedd", "häd", "härd", "härr", "höd", "hödd", "hörd", "hörr", "had", "hadd", "håd", "hådd", "hod", "hodd")

# While R can handle unicode, stan cannot. Thus using IPA as values does not work for models
levels.response.natural <- c("heed", "hid", "head", "had", "hut", "hood", "who'd", "odd")
levels.response <- c("heed", "hid", "head", "had", "hud", "hood", "whod", "hod")
levels.response.natural.swe.long <- c("hid", "hyd", "hud", "hed", "häd", "härd", "höd", "hörd", "had", "håd", "hod")
levels.response.natural.swe.long.noallo <- c("hid", "hyd", "hud", "hed", "häd", "höd", "had", "håd", "hod")

levels.cue.names.transform = c("F1_Hz", "F2_Hz", "F3_Hz", "F1_Mel", "F2_Mel", "F3_Mel", "F1_Bark", "F2_Bark", "F3_Bark", "F1_ERB", "F2_ERB", "F3_ERB", "F1_semitones", "F2_semitones", "F3_semitones")
levels.normalization = c("r_Hz", "r_log", "r_Mel", "r_Bark", "r_ERB", "r_semitones", "SyrdalGopal_Bark", "SyrdalGopal2_Bark", "Miller_log", "Nearey2_log", "Nearey1_log", "CCuRE_Hz", "CCuRE_Mel", "CCuRE_Bark", "CCuRE_ERB", "CCuRE_semitones", "Gerstman_Hz", "Lobanov_Hz")
labels.normalization = c("no normalization (Hz)", "transformed (log)", "transformed (Mel)", "transformed (Bark)", "transformed (ERB)","transformed (semitones)", "SyrdalGopal (Bark)", "SyrdalGopal2 (Bark)", "Miller (log)", "Nearey's uniform scaling (log)", "Nearey's formantwise mean (log)", "C-CuRE (Hz)", "C-CuRE (Mel)", "C-CuRE (Bark)", "C-CuRE (ERB)", "C-CuRE (semitones)", "Gerstman (Hz)", "Lobanov (Hz)")
levels.norm = c("r", "SyrdalGopal", "SyrdalGopal2", "Miller", "CCuRE", "Nearey1", "Nearey2", "Gerstman", "Lobanov")

# Color codes for plotting
colors.vowel <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
names(colors.vowel) <- levels.vowel.IPA
colors.vowel.word <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
names(colors.vowel.word) <- levels.response.natural
# repeat the same code to assign the same colours to IO.Vowel in order to plot posteriors in the same colours
colors.vowel.IO <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
names(colors.vowel.IO) <- levels.response
# Colors for plotting IOs categorization accuracy
colScale <-   scale_colour_manual(name = "Vowel",values = colors.vowel)
colors.vowel.swe <- c("#66C2A5", "#1B9E77", "#FC8D62", "#D95F02", "#8DA0CB", "#7570B3", "#FB9A99", "#E78AC3", "#E7298A", "#A6D854", "#66A61E", "#FFD92F", "#E6AB02", "#E5C494", "#A6761D", "#B3B3B3", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
names(colors.vowel.swe) <- levels.vowel.IPA.swe

vowelPlot_components <- list(
  axis = c(scale_x_reverse(), scale_y_reverse()),
  scale = colScale, 
  theme(legend.position="right"), 
  guides(color = guide_legend(order=1), shape = guide_legend(order=2)))

# Using grob.element
element_textbox_highlight <- function(
    ...,
    hi.labels = NULL, hi.fill = NULL,
    hi.col = NULL, hi.box.col = NULL,
    hi.labels2 = NULL, hi.fill2 = NULL,
    hi.col2 = NULL, hi.box.col2 = NULL,
    hi.labels3 = NULL, hi.fill3 = NULL,
    hi.col3 = NULL, hi.box.col3 = NULL
) {
  structure(
    c(ggtext::element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col,
           hi.labels2 = hi.labels2, hi.fill2 = hi.fill2, hi.col2 = hi.col2, hi.box.col2 = hi.box.col2,
           hi.labels3 = hi.labels3, hi.fill3 = hi.fill3, hi.col3 = hi.col3, hi.box.col3 = hi.box.col3)),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element",
              "element_textbox_highlight", "element_textbox", "element_text", "element",
              "element_textbox_highlight", "element_textbox", "element_text", "element"))
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
  }
  if (label %in% element$hi.labels2) {
    element$fill <- element$hi.fill2 %||% element$fill
    element$colour <- element$hi.col2 %||% element$colour
    element$box.colour <- element$hi.box.col2 %||% element$box.colour
  }
  if (label %in% element$hi.labels3) {
    element$fill <- element$hi.fill3 %||% element$fill
    element$colour <- element$hi.col3 %||% element$colour
    element$box.colour <- element$hi.box.col3 %||% element$box.colour
  }
  NextMethod()
}


## STAN / BRMS constants ---------------------------------------------------------------
chains <- 4
options(
  width = 110,
  mc.cores = min(chains, parallel::detectCores()))

require(brms)
my_priors <- c(
  prior(student_t(3, 0, 2.5), class = "b"),
  prior(cauchy(0, 2.5), class = "sd"),
  prior(lkj(1), class = "cor")
)

# Functions ------------------------------------------------------------------
myGplot.defaults = function(
    type = c("paper","poster","slides")[1],
    base_size = if (type == "paper") { 10 } else if (type == "slides") { 32 } else if (type == "poster") { 36 } else { 10 },
    margin=c("t" = 0.6, "r" = 0.5, "b" = 0.5, "l" = 0.3),
    set_theme = T
)
{
  require(ggplot2)
  t <- theme(axis.text.x = element_text(size=base_size, vjust=1),
             axis.text.y = element_text(size=base_size, hjust=1, vjust=.5),
             axis.title.x = element_text(size=base_size , vjust=0, hjust=0.5, face = "bold"),
             axis.title.y = element_text(size=base_size, hjust= 0.5, vjust=0.5, face = "bold"),
             strip.text = element_text(size=base_size, color = "white"),
             strip.background = element_rect(fill = "black", color = "black"),
             legend.title = element_text(size=base_size, face = "bold", hjust= 0),
             legend.text = element_text(size=base_size),
             plot.margin = unit(margin, "lines"),
             aspect.ratio = 1)
  
  if (set_theme) theme_set(theme_bw(base_size=base_size) + t) else return(t)
}

## General functions ------------------------------------------------------------------
se <- function(x) sqrt(var(x) / (length(x) - 1))

scale_Gelman <- function(x) {
  (x - mean(x)) / (2 * sd(x))
}

my_unscale_Gelman <- function(d, scale_summary) {
  d %>%
    mutate(
      Item.Height = sItem.Height * 2 * scale_summary$Item.Height_sd + scale_summary$Item.Height_mean,
      Item.Backness = sItem.Backness * 2 * scale_summary$Item.Backness_sd + scale_summary$Item.Backness_mean)
}

geometric.mean = psych::geometric.mean

get_vowel_stats <-
  . %>%
  select(-F0_gm) %>%
  summarise(
    across(ends_with("gm"),
           .fns = list("mean" = mean, "var" = var)),
    F1F2_cov = cov(F1_gm, F2_gm), 
    F1F3_cov = cov(F1_gm, F3_gm),
    F2F3_cov = cov(F2_gm, F3_gm),
    heightbackness_cov = cov(height_gm, backness_gm),
    synthtalkerF0_gm = synthtalker_F0)

split_data <- function(
  data,
  proportion_data = .2
) {
  n_fold <- round(1 / proportion_data)
  
  data %>%
    group_by(Talker, category) %>%
    mutate(
      fold = base::sample(
        x = rep(1:n_fold, first(length(Talker)) %/% n_fold + 1), 
        size = first(length(Talker)), 
        replace = F))
}

## Normalization functions ------------------------------------------------------------
F0_to_SR = function(F0) {
  168 * (F0 / 168)^(1/3)
}

F1_to_height = function(F1, SR) {
  log(F1 / SR)
}

F2_to_backness = function(F1, F2) {
  log(F2 / F1)
}

height_to_F1 = function(height, SR) {
  exp(height) * SR
}

backness_to_F2 = function(backness, height, SR) {
  exp(backness) * height_to_F1(height, SR)
}

# Aliases
get_F1_Hz_from_height <- height_to_F1

get_F2_Hz_from_backness <- backness_to_F2

add_Hz_from_height_backness <- function(
  data,
  var.height = "Item.Height",
  var.backness = "Item.Backness",
  var.SR = "SR"
) {
  require(rlang)
  
  data %>%
    mutate(
      Item.Cue_Hz_F1 = height_to_F1(!! sym(var.height), !! sym(var.SR)),
      Item.Cue_Hz_F2 = backness_to_F2(!! sym(var.backness), !! sym(var.height), !! sym(var.SR)))
}

add_Mel_from_Hz <- function(
  data,
  var.F1 = "Item.Cue_Hz_F1",
  var.F2 = "Item.Cue_Hz_F2",
  var.F3 = "Item.Cue_Hz_F3"
) {
  require(rlang)
  
  data %>%
    mutate(
      Item.Cue_Mel_F1 = phonR::normMel(!! sym(var.F1)),
      Item.Cue_Mel_F2 = phonR::normMel(!! sym(var.F2)),
      Item.Cue_Mel_F3 = phonR::normMel(!! sym(var.F3)))
}

add_normalized <- function(
  data,
  var.cue.names = levels.cue.names
) {
  data %>%
    mutate_at(
      var.cue.names, 
      .funs = 
        c("r" = function(x) x,
          "c" = function(x) x - mean(x, na.rm = T),
          "s" = function(x) as.numeric(scale(x)))) %>%
    select(-(all_of(var.cue.names)))
}

get_transformation <- function(
    data,
    cues = c("F0", "F1", "F2", "F3")
) {
  data %>%
    mutate(
      across(
        c(!!! syms(cues)),
        c("log" = log,
          "Mel" = function(x) phonR::normMel(x),
          "Bark" = function(x) phonR::normBark(x),
          "ERB" = function(x) phonR::normErb(x),
          "semitones" = function(x) 12 * log2(x/100))))
}

# Takes a normalization name as argument, along with a named list of parameter vectors.
# The name of each list element is the name of the statistics that the parameters sum-
# marizes (e.g., "mean", "min", etc.). The elements of the vector are the values for 
# that statistic for all formants that are meant to be normalized (e.g., a 3-element 
# vector if F1-F3 are meant to be normalized by the returned function). 
# 
# The output of this function is a function that can take vectors of cues values as 
# input and returns vectors of normalized cue values.
# (Could potentially be extended to also take *lists* of cue values (i.e., multiple 
# observations) and then return lists of normalized cue values.)
get_normalization_function_from_parameters <- function(normalization, ..., verbose = F) {
  normalization_params <- list(...)[[1]]
  
  # if (verbose) {
  #   message("Returning normalization function for ", normalization, " with ", 
  #           length(normalization_params), " parameters ",
  #           paste(normalization_params, collapse = " and "))
  # }
  
  if (normalization %in% c("CCuRE", "Nearey1", "Nearey2")) {
    function(x) x - unlist(normalization_params[[1]])
  } else if (normalization == "Gerstman") { 
    function(x) 999 * (x - unlist(normalization_params[[2]])) / 
      (unlist(normalization_params[[1]]) - unlist(normalization_params[[2]]))
  } else if (normalization == "Lobanov") {
    function(x) (x - unlist(normalization_params[[1]])) / unlist(normalization_params[[2]])
  } else message("Unrecognized normalization procedure: ", normalization)
}

# Assumes that F0, F1, F2, and F3 exist in the data.
# (if Mel, Bark, or alike are meant to be normalized, then the F0-F3_* columns need to be renamed to F0-F3)
get_normalization_functions_from_data <- function(
    data, 
    normalize_based_on_fold_types = c("training")
) {
  message(paste("Making classic formant normalization functions based on data in", 
                paste(normalize_based_on_fold_types, collapse = ", "), "folds."))
  talker_statistics <-
    data %>%
    filter(fold_type %in% normalize_based_on_fold_types) %>% 
    group_by(Talker) %>%
    summarise(
      formants_mean = list(c(mean(F0, na.rm = T), mean(F1, na.rm = T), mean(F2, na.rm = T), mean(F3, na.rm = T))),
      formants_min = list(c(min(F0, na.rm = T), min(F1, na.rm = T), min(F2, na.rm = T), min(F3, na.rm = T))),
      formants_max = list(c(max(F0, na.rm = T), max(F1, na.rm = T), max(F2, na.rm = T), max(F3, na.rm = T))),
      formants_sd = list(c(sd(F0, na.rm = T), sd(F1, na.rm = T), sd(F2, na.rm = T), sd(F3, na.rm = T))),
      formants_mean_log = list(c(mean(log(F0), na.rm = T), mean(log(F1), na.rm = T), mean(log(F2), na.rm = T), mean(log(F3), na.rm = T))),
      # For Nearey overall logmean *both* formants are normalized by the same quantity: 
      # the mean of the three means of log-transformed F1-F3. We thus make a 4-element 
      # vector with the same mean of log formant means.
      formants_overall_mean_log = list(rep(mean(c(log(F1), log(F2), log(F3)), na.rm = T), 4)))
  
  # use constants for normalization functions, when needed
  f <- function(newdata) {
    message("Applying formant normalization functions to data.")
    newdata %<>%
      left_join(talker_statistics, by = "Talker") %>%
      mutate(
        formants = pmap(
          list(F0, F1, F2, F3), cbind),
        log_Nearey1 = map2(
          formants, formants_mean_log,
          ~ log(.x) - .y),
        log_Nearey2 = map2(
          formants, formants_overall_mean_log,
          ~ log(.x) - .y),
        Hz_Gerstman = pmap(
          list(formants, formants_min, formants_max),
          function(x, y, z) 999 * ((x - y) / (z - y))),
        Hz_Lobanov = pmap(
          list(formants, formants_mean, formants_sd), 
          function(x, y, z) (x - y) / z), 
        across(
          c(log_Nearey1, log_Nearey2, Hz_Gerstman, Hz_Lobanov),
          list("F0" = ~ unlist(map(.x, function(x) x[1])),
               "F1" = ~ unlist(map(.x, function(x) x[2])),
               "F2" = ~ unlist(map(.x, function(x) x[3])),
               "F3" = ~ unlist(map(.x, function(x) x[4]))), 
          .names = "{.fn}_{.col}"),
        # Add Miller 
        # (since this is intrinsic, normalizing params are taken from newdata rather than training_talker_statistics)
        F0_gm = geometric.mean(F0),
        F0_log_Miller = log(F0),
        F1_log_Miller = F1_to_height(F1 = .data$F1, SR = 168 * (.data$F0_gm / 168) ^ (1/3)),
        F2_log_Miller = F2_to_backness(.data$F1, .data$F2),
        F3_log_Miller = log(.data$F3 / .data$F2),
        # Add Syrdal-Gopal's Bark-distance model 
        # (since this is intrinsic, normalizing params are taken from newdata rather than training_talker_statistics)
        F0_Bark_SyrdalGopal = F0_Bark,
        F0_Bark_SyrdalGopal2 = F0_Bark,
        F1_Bark_SyrdalGopal = .data$F1_Bark - .data$F0_Bark,
        F1_Bark_SyrdalGopal2 = .data$F1_Bark - .data$F0_Bark, 
        F2_Bark_SyrdalGopal = .data$F2_Bark - .data$F1_Bark,
        F2_Bark_SyrdalGopal2 = .data$F3_Bark - .data$F2_Bark) %>%
      select(-c(F0_gm, Hz_Lobanov, log_Nearey1, log_Nearey2, Hz_Gerstman))
    
    return(newdata)
  }
  
  return(f)
}

# For legacy use only
get_normalization <- get_normalization_functions_from_data

# Assumes that the two formants are called F1 and F2 and that F0 and F3 exists, too.
# (if Mel, Bark, or alike are meant to be normalized, then the F0-F3_* columns need to be renamed to F0-F3)
get_C_CuRE_function <- function(
    data, 
    cues = levels.cue.names.transform
) {
  data %<>% 
    group_by(Talker) %>%
    summarise(
      across(
        .cols = .env$cues,
        .fns = list("overall_mean_for_CCuRE" = ~ mean(.x, na.rm = T)),
        .names = "{.fn}_{.col}"))
  
  f <- function(newdata) {
    require(glue)
    message(paste("Applying C-CuRE normalization to cues:", paste(cues, collapse = ",")))
    
    newdata %<>%
      left_join(data, by = "Talker") %>%
      mutate(
        across(
          .cols = .env$cues,
          .fns = list("CCuRE" = ~ .x - get(glue("overall_mean_for_CCuRE_{cur_column()}"))),
          .names = "{.col}_{.fn}")) %>%
      select(-starts_with("overall_mean_for_CCuRE_"))
    
    return(newdata)
  }
  
  return(f)
}

# Apply transformations and normalizations to the test data 
apply_all_transformations_and_normalization <- function(
    data, 
    normalize_based_on_fold_types
) {
  data %>%
    get_transformation() %>%
    get_normalization_functions_from_data(
      data = ., 
      normalize_based_on_fold_types = normalize_based_on_fold_types)() %>%
    rename(F0_Hz = F0, F1_Hz = F1, F2_Hz = F2, F3_Hz = F3) %>%
    add_C_CuRE(
      data = .,
      cues = c(
        "F0_Hz", "F1_Hz", "F2_Hz", "F3_Hz", 
        "F0_Mel", "F1_Mel", "F2_Mel", "F3_Mel", 
        "F0_ERB", "F1_ERB", "F2_ERB", "F3_ERB", 
        "F0_Bark", "F1_Bark", "F2_Bark", "F3_Bark", 
        "F0_semitones", "F1_semitones", "F2_semitones", "F3_semitones", 
        "Duration"),
      normalize_based_on_fold_types = normalize_based_on_fold_types)() %>%
    # Add '_r' for 'raw' to columns with scale-transformed data in order for pivoting to work in next chunk
    rename_with(
      .fn = ~ paste(.x, "r", sep = "_"), 
      .cols = matches(c("F[0-3].*Hz$", "F[0-3].*log$", "F[0-3].*Mel$", "F[0-3].*Bark$", "F[0-3].*ERB$", "F[0-3].*semitones$")))
}


# Keep around for legacy use but get rid of it once it's not used anymore (to avoid confusion)
add_C_CuRE <- function(
    data, 
    normalize_based_on_fold_types = c("training"), 
    cues = levels.cue.names.transform
) {
  message(paste("Making C-CuRE normalization functions based on data in", normalize_based_on_fold_types, "folds."))
  return(get_C_CuRE_function(data %>% filter(fold_type %in% normalize_based_on_fold_types), cues = cues))
}

## Outlier detection and correction --------------------------------------------------
get_cumulative_probability = function(x1, x2, mean, sigma) {
  pmvnorm(lower = -Inf, upper = c(x1, x2), mean = mean, sigma = sigma)
}

get_cumulative_probability_allCues = function(x1, x2, x3, x4, x5, mean, sigma) {
  pmvnorm(lower = -Inf, upper = c(x1, x2, x3, x4, x5), mean = mean, sigma = sigma)
}

get_cumulative_probability_univariate = function(x, mean, sd) {
  pnorm(x, mean = mean, sd = sd)
}

is_outlier = function(x, cutoff = outlier_probability_cutoff) {
  !between(
    x, 
    cutoff / 2, 
    1 - cutoff / 2)
}

obtain_densities <- . %>%
  group_by(Talker, category) %>%
  nest() %>%
  mutate(
    # density based on F1 and F2
    x_F1F2 = map(data, ~ cbind(.$F1, .$F2)),
    mean_F1F2 = map(x_F1F2, ~ colMeans(.x)),
    cov_F1F2 = map(x_F1F2, ~ cov(.x)),
    x_F1F2 = NULL) %>%
  unnest(data) %>%
  # normalize densities within each talker and vowel
  mutate(
    cumulative_probability_F1F2 = pmap(
      .l = list(F1, F2, mean_F1F2, cov_F1F2), 
      .f = get_cumulative_probability)) %>%
  mutate_at(
    vars(cumulative_probability_F1F2),
    unlist) %>%
  ungroup()

# obtain_densities_univariates <- function(
#   data,
#   cue = "cue"
# ) {
#   data %>%
#     group_by(Talker, Vowel) %>%
#     nest() %>%
#     mutate(
#       x_cue = map(data, ~cbind(!!! syms(cue))),
#       mean = map(x_cue, ~mean(.x, na.rm = T)),
#       sd = map(x_cue, ~sd(.x))) %>%
#     unnest(data) %>%
#     mutate(
#       "cumulative_probability_" = pmap(
#         .l = list(x, mean, sd),
#         .f = get_cumulative_probability_univariate)) %>%
#     mutate_at(
#       vars(starts_with("cumulative_probability")),
#       unlist) %>%
#     ungroup()
# }

#Obtain univariate densities for outlier correction of all cues in raw Hz
obtain_densities_univariates <- . %>%
  group_by(Talker, Vowel) %>%
  nest() %>%
  mutate(
    # density based on F0
    x_F0 = map(data, ~ cbind(.$F0)),
    mean_F0 = map(x_F0, ~ mean(.x, na.rm = TRUE)),
    sd_F0 = map(x_F0, ~ sd(.x)),
    # density based on F1
    x_F1 = map(data, ~ cbind(.$F1)),
    mean_F1 = map(x_F1, ~ mean(.x)),
    sd_F1 = map(x_F1, ~ sd(.x)),
    # density based on F2
    x_F2 = map(data, ~ cbind(.$F2)),
    mean_F2 = map(x_F2, ~ mean(.x)),
    sd_F2 = map(x_F2, ~ sd(.x)),
    # density based on F3
    x_F3 = map(data, ~ cbind(.$F3)),
    mean_F3 = map(x_F3, ~ mean(.x)),
    sd_F3 = map(x_F3, ~ sd(.x)),
    # density based on Duration
    x_Duration = map(data, ~ cbind(.$Duration)),
    mean_Duration = map(x_Duration, ~ mean(.x)),
    sd_Duration = map(x_Duration, ~ sd(.x)),
    x_F0 = NULL,
    x_F1 = NULL,
    x_F2 = NULL,
    x_F3 = NULL,
    x_Duration = NULL) %>%
  unnest(data) %>%
  # normalize densities within each talker and vowel
  mutate(
    cumulative_probability_F0 = pmap(
      .l = list(F0, mean_F0, sd_F0), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F1 = pmap(
      .l = list(F1, mean_F1, sd_F1), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F2 = pmap(
      .l = list(F2, mean_F2, sd_F2),
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F3 = pmap(
      .l = list(F3, mean_F3, sd_F3), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F3 = pmap(
      .l = list(F3, mean_F3, sd_F3), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_Duration = pmap(
      .l = list(Duration, mean_Duration, sd_Duration), 
      .f = get_cumulative_probability_univariate)) %>%
  mutate_at(
    vars(starts_with("cumulative_probability")),
    unlist) %>%
  ungroup()

#Obtain multivariate densities for outlier correction of all cues in raw Hz
obtain_densities_allCues <- function(data) {
  data %>%
  group_by(Talker, Vowel) %>%
  nest() %>%
  mutate(
    # density based on all cues
    x_cues = map(data, ~ cbind(.$F0, .$F1, .$F2, .$F3, .$Duration)),
    mean_cues = map(x_cues, ~ colMeans(.x)),
    cov_cues = map(x_cues, ~ cov(.x))) %>%
  unnest(data) %>%
  # normalize densities within each talker and vowel
  mutate(
    cumulative_probability_allCues = pmap(
      .l = list(F0, F1, F2, F3, Duration, mean_cues, cov_cues), 
      .f = get_cumulative_probability_allCues)) %>%
  mutate_at(
    vars(cumulative_probability_allCues),
    unlist) %>%
  ungroup()
}

# Vowel plot function
plot_vowels <- function(data, x, y) {
  ggplot(data,
        aes(
          x = F2_gm, 
          y = F1_gm,
          color = Vowel)) +
    geom_point(alpha = .5) +
    scale_x_reverse() +
    scale_y_reverse() +
    vowelPlot_components
}

