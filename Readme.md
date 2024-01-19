Estimation of introduction and transmission rates of SARS-CoV-2 in a
prospective household study
================

## Introduction

This page contains an overview of the data, model, and main analyses presented in our 
[medRxiv preprint](https://doi.org/10.1101/2023.06.02.23290879).
Full information is given in the R and Stan scripts, which can be found
in the ‘scripts’ directory. Data are available in the ‘data’ directory,
and are sufficient to reproduce all results in the manuscript. Here we
present the analyses for the model with proportional mixing in the
household and a separate transmission rate for child-to-child
transmission. Small modifications are needed to run all model scenarios
and perform model selection using WBIC and LOO_IC.

## Data

The main data are available in two files called
‘data_finalsize_stan_05072022’ and ‘data_escape_stan_05072022’. The
first contains information on the final size of household outbreaks, and
the second contains information on the period that each person has been
at risk for introduction SARS-CoV-2 into an as yet uninfected household.

``` r
# read stan data files
data.finalsize.stan <- read_csv("data/data_finalsize_stan_05072022")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   household = col_double(),
    ##   j1 = col_double(),
    ##   j2 = col_double(),
    ##   j3 = col_double(),
    ##   a1 = col_double(),
    ##   a2 = col_double(),
    ##   a3 = col_double(),
    ##   n1 = col_double(),
    ##   n2 = col_double(),
    ##   n3 = col_double(),
    ##   c1 = col_double(),
    ##   c2 = col_double(),
    ##   c3 = col_double(),
    ##   outbreak_start = col_double(),
    ##   outbreak_end = col_double(),
    ##   conditioning = col_double()
    ## )

``` r
head(data.finalsize.stan)
```

    ## # A tibble: 6 x 16
    ##   household    j1    j2    j3    a1    a2    a3    n1    n2    n3    c1    c2
    ##       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1         1     0     0     0     0     0     1     2     2     1     0     0
    ## 2         3     0     0     0     0     0     1     3     1     1     0     0
    ## 3         5     0     0     0     0     1     1     0     0     1     0     0
    ## 4         7     0     0     0     0     1     0     0     1     2     0     0
    ## 5        13     0     0     0     0     0     1     0     1     1     0     0
    ## 6        24     0     0     0     0     1     0     1     0     2     0     0
    ## # … with 4 more variables: c3 <dbl>, outbreak_start <dbl>, outbreak_end <dbl>,
    ## #   conditioning <dbl>

This file contains the so-called final size of the outbreak in the
**a**, **n**, **j** format, where **a** is the vector of type-specific
primary cases, **n** contains the non-primary persons in the household,
and **j** contains the number of household infections. See
<https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/distribution-of-general-final-state-random-variables-for-stochastic-epidemic-models/D685556D125E03466B82E73382380C51>
for details. Throughout, we distinguish between children (0-12 years),
adolescents (12-18 years), and adults (over 18 years).

``` r
# read stan data files
data.escape.stan <- read_csv("data/data_escape_stan_05072022")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   household = col_double(),
    ##   time_start = col_double(),
    ##   time_end = col_double(),
    ##   Type = col_double(),
    ##   ext_infected = col_double()
    ## )

``` r
head(data.escape.stan)
```

    ## # A tibble: 6 x 5
    ##   household time_start time_end  Type ext_infected
    ##       <dbl>      <dbl>    <dbl> <dbl>        <dbl>
    ## 1         1         39      102     3            0
    ## 2         1         39      102     3            1
    ## 3         1         39      102     2            0
    ## 4         1         39      102     2            0
    ## 5         1         39      102     2            0
    ## 6         1         39      102     1            0

This file contains, for each person in the study, information of the
start and end of the at-risk period, the household to which he/she
belongs, and the outcome.

In addition, Dutch hospital admission data are downloaded from
<https://data.rivm.nl/covid-19/> for comparison with the estimated
hazards of introduction.

``` r
# import hospitalisations
library(jsonlite)
```

    ## 
    ## Attaching package: 'jsonlite'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

``` r
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
head(hospitalisations)
```

    ## # A tibble: 6 x 2
    ##   Date_of_statistics total
    ##   <chr>              <int>
    ## 1 2020-08-24            11
    ## 2 2020-08-25            19
    ## 3 2020-08-26            12
    ## 4 2020-08-27             9
    ## 5 2020-08-28            13
    ## 6 2020-08-29             4

## Stan model

Inference is based on the final size of the household epidemics. These
are calculated iteratively in the function ‘prob_infect_pattern’:

      /* iterative calculation of the household infection probabilities */
      real prob_infect_pattern(array[] int j, array[] int a, array[] int n, matrix beta, vector b) {
        int num_types = num_elements(n);
        array[num_types] int jp1 = add_one(j);
        int dimP = prod(jp1); // dimension of cache array
        // cache array: we need space for 0,...,j elements (i.e. j+1)
        array[dimP] real P = rep_array(0.0, dimP);
        // now compute all elements of P iteratively
        for ( idxP1 in 1:dimP ) {
          // convert index of P to a multi-index (or household state) omega
          array[num_types] int omega = unravel_idx(idxP1-1, jp1);
          // pre-compute phi(beta (n-omega))
          vector[num_types] phi_omega = phi_vec(omega, n, beta);
          real S = vec_power(phi_omega, a) * vec_power(b, n) / vec_power(b, omega);
          // we now have to compute prod(omega+1) terms using cached values in P
          for ( idx in 1:prod(add_one(omega))-1 ) {
            array[num_types] int wau = unravel_idx(idx-1, add_one(omega));
            // we have to find the correct flat index in P
            int idxP2 = ravel_idx(wau, add_one(j)) + 1;
            S -= P[idxP2] * binom_prod(omega, wau) / vec_power(phi_omega, wau);
          }
          P[idxP1] = S * vec_power(phi_omega, omega);
        }
        return P[dimP] * binom_prod(n, j);
      }

For estimation of the hazards of introduction of SARS-CoV-2 into the
households we employ penalised splines with first-order penalisation.
See <https://www.jstor.org/stable/1391151> for details. Our
implementation is based on code by Milad Kharratzadeh
(<https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html>).

Estimated parameters are defined in the ‘parameters’ block of the Stan
program, which in our case is given by

    parameters {
      vector<lower=0>[num_types-1] rel_susceptibility;          // relative susceptibilies (compared to reference class) 
      vector<lower=0>[num_types] infectivity;                   // infectivity parameters
      real<lower=0> extra_trans;                                // to accomodate potential additional specific transmsision 
      matrix[1, num_basis] ext_hazard_weights;                  // spline weights, 1 spline, other groups with multiplier
      real<lower = 0> RWvar;                                    // variance of RW1/2 parameter, see Lang & Brezger (2004)
      real<lower = 0> ext_hazard_children;                      // 1 spline for adults, other types multiplicative to adults
      real<lower = 0> ext_hazard_adolescents;
    }

Here we estimate susceptibility of children and adolescents relative to
adults (‘rel_susceptibilty’), and absolute infectivity of all
person-types (‘infectivity’). A separate parameter is estimated for
child-to-child transmission (‘extra_trans’). The weights of the spline
defining the introduction hazard for adults are given by
‘ext_hazard_weights’, and the relative hazards for children and
adolescents relative to adults are given by ‘ext_hazard_children’ and
‘ext_hazard_adolescents’.

The log-likelihood contributions of all households are specified in the
‘transformed parameters’ block, and include the household outbreaks (for
those households that are infected), and survival hazards of all persons
in the households. Notice that external infections are also allowed
during the household outbreaks. The hazards of external infection are
estimated to be very small relative to the infection rates within the
households, and therefore affect the results only marginally.

      /* log-likelihood contributions per household (survival and outbreak) */
      for (i in 1 : num_households) {
         log_lik[i] = 0.0;
        }
      
      for (i in 1: num_households_infected) { // household outbreaks
         int hh_size = sum(A[i,:]) + sum(N[i,:]);
         vector[num_types] ext_esc;
        for ( k in 1 : num_types) {
            ext_esc[k] = exp(-sum(ext_hazards[D[i,1] : D[i,2], k])); 
           };   
         /* notice that it is possible to switch between density- and frequency-dependent models */
         /* by dividing transmission matrix function call by hh_size (frequency-dependence)      */
         /* or by 1 (density-dependence). In the former, transmision rates are per infectious    */
         /* period, and in the latter rates are per infectious period per person.                */
         /* also notice that continous time hazards are replaced by fixed hazards per day.       */
         /* this leads to faster code but with small differences with fully continuous model     */
         log_lik[id_infected_household[i]] = log(prob_infect_pattern(J[i,:], A[i,:], N[i,:], transmission_rate, ext_esc)); 
        }
        
      for ( i in 1 : num_persons ) { // escapes
         int j = id_household[i]; 
         real total_hazard = sum(ext_hazards[hazard_times[i,1] : hazard_times[i,2]-1, hazard_times[i,3]]);
         if ( hazard_times[i,4] == 0 ) { // no infection on last day
            log_lik[j] += -ext_hazards[hazard_times[i,2], hazard_times[i,3]] - total_hazard;  
           } else { // primary infection on last day
            log_lik[j] += log(ext_hazards[hazard_times[i,2], hazard_times[i,3]]) - total_hazard; 
           };           
        }

Finally, the ‘generated quantities’ block contains code to calculate
secondary attack rates for households of different compositions and
different numbers and types of primary cases. The R script also contains
code to reproduce the figures of the manuscript.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
