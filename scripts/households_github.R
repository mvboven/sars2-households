######################################################################################
# SARS-CoV-2 household analyses for the Dutch CoKids study. Code is copyrighted by   #
# Michiel van Boven and licensed with BSD3 clause (use but mention). Created: 7/2022 #
######################################################################################


# load packages
library(tidyverse)
library(cmdstanr)
set_cmdstan_path("~/.cmdstanr/cmdstan-2.29.1/")
library(loo)

# colors and folders
cbPalette <-
  c("#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7")
dir_figures = "figures/"
dir_data = "data/"
dir_results = "output/"
dir_cache = "cmdstan-cache/"

# some preparations
num_households <- 307
num_households_infected <- 59
num_persons <- 1209
date.start <- as.Date("2020-08-24") # Stan needs 1-indexing
date.end <- as.Date("2021-07-29")

# read stan data files
data.finalsize.stan <- read_csv("data/data_finalsize_stan_05072022")
data.escape.stan <- read_csv("data/data_escape_stan_05072022")

# import hospitalisations
library(jsonlite)
hospitalisations <-
  fromJSON(
    "https://data.rivm.nl/covid-19/COVID-19_ziekenhuisopnames_tm_03102021.json",
    flatten = TRUE
  )
hospitalisations <- hospitalisations %>%
  select(Date_of_statistics, Municipality_code, Hospital_admission) %>%
  group_by(Date_of_statistics) %>%
  summarise(total = sum(Hospital_admission)) %>%
  filter(Date_of_statistics >= date.start,
         Date_of_statistics <= date.end)


# Stan preparations
num_days <- as.numeric(date.end - date.start + 1)
num_knots <- 50
spline_degree <- 3
knots <-
  seq(from = 1,
      to = num_days,
      by = (num_days - 1) / (num_knots - 1))
num_types <- 3
household.data.list <- list(
  num_types = num_types,
  num_days = num_days,
  num_households = num_households,
  num_households_infected = num_households_infected,
  num_persons = num_persons,
  ts = 1:num_days,
  spline_degree = spline_degree,
  num_knots = num_knots,
  knots = knots,
  J = as.matrix(data.finalsize.stan[, c("j1", "j2", "j3")]),
  A = as.matrix(data.finalsize.stan[, c("a1", "a2", "a3")]),
  N = as.matrix(data.finalsize.stan[, c("n1", "n2", "n3")]),
  C = as.matrix(data.finalsize.stan[, c("c1", "c2", "c3")]),
  D = as.matrix(data.finalsize.stan[, c("outbreak_start", "outbreak_end")]),
  conditioning = as.vector(pull(data.finalsize.stan[, c("conditioning")])),
  hazard_times = as.matrix(data.escape.stan[, c("time_start", "time_end", "Type", "ext_infected")]),
  id_household = as.vector(pull(data.escape.stan[, c("household")])),
  id_infected_household = as.vector(pull(data.finalsize.stan[, c("household")])),
  number_subjects = nrow(data.escape.stan),
  mode = 0 # 0=normal sampling; 1=WBIC sampling
)

# initial values
initials = function() {
  return(
    list(
      rel_susceptibility = c(1.5, 0.7),
      infectivity = c(0.3, 0.1, 0.1),
      extra_trans = 1.0,
      RWvar = 0.5,
      ext_hazard_children = 0.6,
      ext_hazard_adolescents = 1.1,
      ext_hazard_weights = matrix(-1.0, nrow = 1, ncol = num_knots + spline_degree - 1)
    )
  )
}


# run Stan model
set_cmdstan_path("~/.cmdstanr/cmdstan-2.29.1/")

stan_file <- "scripts/households_github.stan"

stan_model <- cmdstan_model(stan_file,
                            dir  = "~/.cmdstanr/cmdstan-2.29.1/",
                            force_recompile = TRUE)

chains = 10
iter_sampling = 1000
thin = 10
stan_fit <- stan_model$sample(
  data = household.data.list,
  init = initials,
  seed = 123,
  chains = chains,
  parallel_chains = 10,
  iter_warmup = 1000,
  iter_sampling = iter_sampling,
  thin = thin,
  save_warmup = FALSE,
  refresh = 100,
  output_dir = "cmdstan-cache",
  adapt_delta = 0.99,
  max_treedepth = 15
)

# transform to rstan object and extract parameters
stan_fit <- rstan::read_stan_csv(stan_fit$output_files())
params = rstan::extract(stan_fit)
n_samples = chains * iter_sampling / thin

# model selection using WBIC
print(stan_fit, pars = "WBIC", digits = 4) # nb set mode=1

# model selection using LOO_IC
loo_unstr <- loo(stan_fit) # mode == 0

# some checks
pp <-
  c(
    "beta",
    "rel_infectivity[1]",
    "rel_infectivity[2]",
    "rel_susceptibility[1]",
    "rel_susceptibility[2]",
    "extra_trans",
    "ext_hazard_children",
    "ext_hazard_adolescents",
    "RWvar"
  )
rstan::traceplot(stan_fit, pars = pp)
pairs(stan_fit, pars = pp)
print(
  stan_fit,
  pars = pp,
  digits = 4,
  probs = c(0.025, 0.5, 0.975)
)

# figure 1
cokids.lexis <- data.escape %>%
  group_by(household_id) %>%
  mutate(
    ext_infected = case_when(
      ext_infected == 'primary' ~ as.numeric(1),
      ext_infected == 'coprimary' ~ as.numeric(1),
      ext_infected == 'at risk' ~ as.numeric(0),
      ext_infected == 'secondary case' ~ as.numeric(0)
    )
  ) %>%
  summarise(
    housedhold_id = household_id,
    inf = ifelse(sum(ext_infected) > 0, 1, 0),
    start = esc_time_start,
    end = min(esc_time_end),
    duration = as.numeric(end - start)
  ) %>%
  unique()

figure_1 <- ggplot(cokids.lexis) +
  geom_segment(aes(
    x = start,
    y = 0,
    xend = end,
    yend = duration,
    color = as.factor(inf)
  ),
  alpha = 0.1) +
  geom_point(aes(x = end, y = duration, color = as.factor(inf)),
             alpha = 0.5,
             size = 2) +
  scale_color_manual(
    name = "",
    labels = c("Uninfected", "Infected"),
    breaks = c(0, 1),
    values = c(cbPalette[5], cbPalette[6])
  ) +
  scale_y_continuous(limits = c(0, 165), expand = c(0, 5)) +
  xlab("Date") +
  ylab("Time in study (days)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

figure_1

# figure 2
h3_med <- apply(params$ext_hazards[, , 3], 2, median)
h3_lb <- apply(params$ext_hazards[, , 3], 2, quantile, probs = 0.05)
h3_ub <- apply(params$ext_hazards[, , 3], 2, quantile, probs = 0.95)

df_hazards3 <- data.frame(
  date = seq(date.start, by = "day", length.out = num_days),
  med = h3_med,
  low = h3_lb,
  upp = h3_ub,
  hospitalisations = as.numeric(hospitalisations$total)
) %>%
  mutate(hospsrolling = zoo::rollapply(hospitalisations, 7, mean, fill = "extend")) %>%
  select(-hospitalisations)

color <- cbPalette[2]
secaxisfac <- 200000
figure_2 <-
  ggplot(df_hazards3, aes(
    x = date,
    y = med,
    ymin = low,
    ymax = upp
  )) +
  geom_ribbon(alpha = 0.25) +
  geom_line(color = color) +
  geom_point(
    aes(y = hospsrolling / secaxisfac),
    size = 1,
    color = cbPalette[1],
    alpha = 0.4
  ) +
  labs(x = "Date", y = expression(paste("Introduction hazard (", day ^ -1, ")", sep = ""))) +
  coord_cartesian(ylim = c(0, 0.002)) +
  scale_y_continuous(
    breaks = c(0, 0.0005, 0.0010, 0.0015, 0.002),
    sec.axis = sec_axis( ~ . * secaxisfac, name = expression(
      paste("Hospital admissions (", day ^ -1, ")", sep = "")
    ))
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title.y = element_text(color = cbPalette[2]),
    axis.title.y.right = element_text(color = cbPalette[1]),
    plot.title = element_text(color = color),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

figure_2

# figure 4
cbPalette <- c("#0072B2", "#E69F00", "#009E73")

simulation.data <- function(probs, letter) {
  probs %>%
    as.data.frame() %>%
    setNames(lev) %>%
    pivot_longer(
      col = lev,
      names_to = c("Primary", "Household"),
      names_sep = "_",
      values_to = "Probability"
    ) %>%
    mutate(Household = factor(Household, levels = c("two", "four", "six"))) %>%
    ggplot(aes(
      x = factor(Primary, levels = c("child", "adolescent", "adult")),
      y = Probability,
      fill = Household
    )) +
    geom_violin(aes(fill = Household), color = NA, alpha = 0.5) +
    stat_summary(
      fun = median,
      geom = "point",
      color = "black",
      alpha = 0.3,
      position = position_dodge(0.9)
    ) +
    geom_text(
      x = Inf,
      y = Inf,
      hjust = 1.5,
      vjust = 1.5,
      label = letter,
      size = 12
    ) +
    scale_y_continuous(limits = c(0, 0.80)) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x = element_text(
        angle = 0,
        vjust = 0.5,
        hjust = 0.5,
        size = 15
      ),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    labs(x = ifelse(letter == "C", "Primary case", " "), y = "Prob(infection)") +
    scale_fill_manual(
      values = cbPalette,
      name = "Household",
      labels = c(
        "1 adult, 1 other",
        "2 adults, 1 adolescent, 1 child",
        "2 adults, 2 adolescents, 2 children"
      )
    )
}

figure_4a <- simulation.data(params$outbreak_probs, "A")
figure_4b <- simulation.data(params$outbreak_probs_vacadults, "B")
figure_4c <- simulation.data(params$outbreak_probs_vacaduladols, "C")

gridExtra::grid.arrange(figure_4a,
                        figure_4b,
                        figure_4c,
                        ncol = 3,
                        widths = c(1, 1, 1))
