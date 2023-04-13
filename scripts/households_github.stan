/* Household final size model. Notation is mostly taken from Ball (1986).    */
/*                                                                           */
/* Let a, n and j be the number of initial infections, the number of         */
/* uninfected individuals in the household, and the number of new            */
/* infections after the outbreak. The model allows for an arbitrary number   */
/* of types / classes of individuals (age classes, sex, etc).                */
/* further, let q(j|a,n,beta) denote the probability of final size j+a,      */
/* given a and n and parameters beta etcetera, then                          */
/* q(j|a,n,beta) = \binom{n}{j} p(j|a,n,beta) where p is defined by          */
/*                                                                           */
/* 1 = \sum_{\omega \leq j} \binom{j}{\omega} p(\omega | a, n) \times        */
/* \prod_{k=1}^K \phi(\sum_{l=1}^K \beta_{lk}(n_l-j_l)/|a+n|)^{-omega_i-a_i} */
/*                                                                           */
/* Here, K denotes the number of classes, and phi is the Laplace transform   */
/* of the pdf of the infectious period.                                      */
/* For instance, in case of a fixed infections period of 1 time unit, we get */
/* \phi(x) = \int_0^{\infty} e^{-tx} \delta_{t-1} dt = e^{-x}                */
/*                                                                           */
/* Also included is a probability of external infection (or rather escape).  */
/* The probabilities p can be computed iteratively as follows:               */
/*                                                                           */
/* p(\omega | a, n) = b^{n-\omega) \phi(beta'(n-\omega))^{a+\omega} -        */
/*      \sum_{\wau < \omega} p(\wau | a, n) \binom{\omega}{\wau} \times      */    
/*              \phi(\beta'(n-\omega))^{\omega - \wau}                       */
/*                                                                           */
/*                                                                           */
/* Code is drafted by Chris van Dorp and Michiel van Boven (october 2020)    */
/* and licensed with BSD3 clause (use but mention).                          */
/*                                                                           */
/* May 2022: This update uses prospective cohort data. Here final size       */
/* analyses apply to the households whenever an introduction has occurred.   */
/* At the same time, a p-spline is fitted to estimate the hazard of          */
/* introduction into households Hazards are related to infection             */
/* probabilities through p(infection(t)) = 1 - exp(-lambda(t)).              */
/* Notice that hazards are type-specific.                                    */

functions {
  /* Laplace transform of the scaled infectious period (E(T_I)=1) with       */
  /* realistic gamma distribution (alpha=beta=50; alt par k=50, theta=1/50). */
  /* With ~95% coverage of 0.75-1.25, eg 6-10 days when mean is 8 days       */
  real phi(real x) {
    return (1+x/50)^(-50);                                  // gamma distributed infectious period with mean 1
    //return (1+x/20)^(-20);                                // variant gamma distributed infectious period	
    //return exp(-x);	                                    // fixed infectious period with mean 1 
  }
  
  /* probability of n-j individuals not being infected by an infected household member */
  vector phi_vec(array[] int j, array[] int n, matrix beta) {
    int num_types = num_elements(n);
    vector[num_types] res;
    for ( k in 1:num_types ) {
      real Bx = 0;
      for ( l in 1:num_types ) {
        Bx += beta[l,k]*(n[l] - j[l]);
      }
      res[k] = phi(Bx);
    }
    return res;
  }
  
  /* compute product of power functions: x^m = x_1^{m_1} * x_2^{m_2} * ... */
  real vec_power(vector x, array[] int m) {
    return prod(pow(x, to_vector(m))); //stan 2.29 with vectorised pow
  }
      
  /* take binomial coefficient for integer vectors as product of components */
  int binom_prod(array[] int n, array[] int k) {
    int B = 1;
    int num_types = num_elements(n);
    for ( i in 1:num_types ) {
      B *= choose(n[i], k[i]);
    }
    return B;
  }
  
  /* adds ones to an array of integers */
  array[] int add_one(array[] int x) {
    int n = num_elements(x);
    array[n] int xp1;
    for ( i in 1:n ) {
        xp1[i] = x[i] + 1;
    }
    return xp1;
  }

  /* convert index on 1d array to multi-index of nd-array
   * NB: both indices are 0-based!
   */
  array[] int unravel_idx(int idx, array[] int shape) {
      int num_types = num_elements(shape);
      array[num_types] int idxs;
      for ( tp in 1:num_types ) {
          idxs[tp] = (idx %/% prod(shape[1:tp-1])) % shape[tp]; //stan 2.29
          //idxs[tp] = (idx / prod(shape[1:tp-1])) % shape[tp]; //stan 2.21
      }
      return idxs;
  }

  /* convert multi-index of nd array to index of 1d-array
   * NB: both indices are 0-based
   */
  int ravel_idx(array[] int idxs, array[] int shape) {
      int idx = 0;
      int num_types = num_elements(shape);
      for ( tp in 1:num_types ) {
          idx += prod(shape[1:tp-1]) * idxs[tp];
      }
      return idx;
  }
  
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
  
  /* construct b-spline basis functions; taken from Kharratzadeh's example         */
  /* see https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html */
  /* this is used to define the hazard of infection from outside the households    */
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  } 
  
  /* for calculation of infection of specific individuals in simulations with 2 adults */
  real pdouble(int j) {  // MEH
    if ( j == 1 ) { 
      return 0.5;
    } else { 
      return 1.0;
    }
  }  
}

data {
  int<lower=1> num_types;                                   // number of types in the household
  int<lower = 1> num_days;                                  // number of days 
  int<lower = 1> num_households;                            // number of observations (all households) 
  int<lower = 1> num_households_infected;                   // number of infected households 
  array[num_households_infected, num_types] int<lower=0> J; // number of household infections
  array[num_households_infected, num_types] int<lower=0> A; // number of initial cases
  array[num_households_infected, num_types] int<lower=0> N; // initially uninfected individuals
  array[num_households_infected, 2] int<lower=0> D;         // start and end of household outbreaks 
  int<lower=1> num_persons;                                 // total number of persons in the study
  array[num_persons, 4] int<lower=0> hazard_times;          // indvidual hazards data (1=start, 2=end, 3=type, 4=infected)
  array[num_persons] int id_household;                      // household id for each of the participants
  array[num_households_infected] int id_infected_household; // infected household ids
  int num_knots;                                            // number of knots
  int spline_degree;                                        // spline degree (order - 1)
  array[num_days] real ts;                                  // points at which hazards are calculated
  vector[num_knots] knots;                                  // sequence of knots 
  int<lower = 0, upper = 1> mode;                           // sampling: 0 = estimation; 1 = WBIC 
}

transformed data {
  int num_basis = num_knots + spline_degree - 1;            // number of B-splines
  matrix[num_basis, num_days] B;                            // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2 * spline_degree + num_knots] ext_knots;          // extended knots
  row_vector[num_days] u_row = rep_vector(1.0, num_days)';           
  real<lower = 0, upper = 1> watanabe_beta;                 // sampling: 0 = normal; 1 = WBIC
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1 : num_basis){
    B[ind,:] = to_row_vector(build_b_spline(ts, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, num_days] = 1;                               // not needed if ts is sorted 
  if ( mode == 0 ) {
    watanabe_beta = 1.0;
  }
  else { // mode == 1
    watanabe_beta = 1.0/log(num_households);     
  }
}

parameters {
  vector<lower=0>[num_types-1] rel_susceptibility;          // relative susceptibilies (compared to reference class) 
  vector<lower=0>[num_types] infectivity;                   // infectivity parameters
  real<lower=0> extra_trans;                                // to accomodate potential additional specific transmsision 
  matrix[1, num_basis] ext_hazard_weights;                  // spline weights, 1 spline, other groups with multiplier
  real<lower = 0> RWvar;                                    // variance of RW1/2 parameter, see Lang & Brezger (2004)
  real<lower = 0> ext_hazard_children;                      // 1 spline for adults, other types multiplicative to adults
  real<lower = 0> ext_hazard_adolescents;
}

transformed parameters {
  matrix<lower=0>[num_types, num_types] transmission_rate;  // transmission rates (unit: per infectious period)
  vector<lower=0>[num_types] susceptibility;                // susceptibility with 1 padded for reference class
  real<lower = 0> beta;                                     // transmission rate in reference class
  vector<lower=0>[num_types-1] rel_infectivity;             // infectivity relative to reference type (last class)
  matrix<lower=0>[num_days, num_types] ext_hazards;         // estimated external infection hazard(s) - for flexible extension
  vector[num_households] log_lik;                           // log-likelihood contributions (household and external)
  matrix[num_types, num_basis] weights;                     // regression coefficients for the p-spline 

  /* calculation of spline weights */  
  weights[3,1] = ext_hazard_weights[1,1];                   // now assuming 1 p-spline for adults (group 3)
  for (i in 2 : num_basis) {				                // RW1 smoothing prior (Lang & Brezger, 2004) 
      weights[3,i] = weights[3,i-1] + ext_hazard_weights[1,i] * sqrt(RWvar); // RW1 prior
    }
  for (i in 1 : num_basis) {				                // exponentiated weights
     weights[3,i] = exp(weights[3,i]);  
    }
  weights[1,] = ext_hazard_children * weights[3,];          // other groups have type-specific multiplier relative to adults
  weights[2,] = ext_hazard_adolescents * weights[3,];

  //ext_hazards = (hazards * u_row)';                       // fixed type-specific hazards
  ext_hazards = (weights * B)';                             // p-spline type-specific hazards      

  /* calculate transmission rates from underlying assumptions   */
  /* notice that the last class (adults) is the reference group */
  susceptibility = append_row(rel_susceptibility, 1.0);
  transmission_rate = susceptibility * infectivity';        // type-to-type transmission rates per infectious period
  transmission_rate[1, 1] = extra_trans;                    // 1=children; 2=adolescents; 3-adults
  
  print(median(D[,2] - D[,1]));
  
  /* fitting the model in terms of infectivity and rel_susceptibility is faster and       */
  /* easier than fitting in terms of beta and relative infectivities and relative         */
  /* susceptibilities because of strong parameter correlations in the latter case. beta   */
  /* and rel_infectivity are still interesting outputs in their own right. notice that    */
  /* proportional mixing models with n types contain 2n-1 independent parameters (not 2n!)*/
  /* the code below assumes that the last type (e.g., adults) is the reference class.     */
  beta = infectivity[num_types];
  rel_infectivity = infectivity[1 : num_types-1]/beta;

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
}

model {
  /* prior distributions */
  ext_hazard_weights[1,1] ~ normal(-7.5, 2.5);              // first weight; note: log-scale
  ext_hazard_weights[1,2:num_basis] ~ normal(0, 10);        // note: log-scale
  RWvar ~ inv_gamma(1, 0.0005);                             // Lang & Brezger (2004) 
 
  /* log-likelihood */
  target += sum(log_lik);                                   // log-likelihood
}

generated quantities {
  real WBIC = -2*sum(log_lik);                              // NB set mode=1
  array[8] real outbreak_probs;
  array[8] real outbreak_probs_vacadults;
  array[8] real outbreak_probs_vacaduladols;
  array[11] real outbreak_size;
  array[11] real outbreak_size_vacadults;
  array[11] real outbreak_size_vacaduladols;
  vector[num_types] ext_esc = rep_vector(1.0, num_types);
  matrix[num_types, num_types] VE_S_adults = [[0, 0, 0], [0, 0, 0], [0.9, 0.9, 0.9]]; //move to data
  matrix[num_types, num_types] VE_I_adults = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];  
  matrix[num_types, num_types] VE_S_aduladols = [[0, 0, 0], [0.9, 0.9, 0.9], [0.9, 0.9, 0.9]]; 
  matrix[num_types, num_types] VE_I_aduladols = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; 
  array[num_types] int nn;
  array[num_types] int aa;
  int ll;
  real hh;

  /* probabilitiers of infection of named adults - Figure 4 in manuscript */
  for (i in 1 : size(outbreak_probs)) {
	outbreak_probs[i] = 0;
	outbreak_probs_vacadults[i] = 0;
	outbreak_probs_vacaduladols[i] = 0;
  };
  
  /* one child and one adult, primary case child */
  hh = 1;
  outbreak_probs[1] += prob_infect_pattern({0,0,1}, {1,0,0}, {0,0,1}, transmission_rate/hh, ext_esc);	  
  outbreak_probs_vacadults[1] += prob_infect_pattern({0,0,1}, {1,0,0}, {0,0,1}, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
  outbreak_probs_vacaduladols[1] += prob_infect_pattern({0,0,1}, {1,0,0}, {0,0,1}, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);	
  
  /* one adolescent and one adult, primary case adolescent */
  outbreak_probs[2] += prob_infect_pattern({0,0,1}, {0,1,0}, {0,0,1}, transmission_rate/hh, ext_esc);	  
  outbreak_probs_vacadults[2] += prob_infect_pattern({0,0,1}, {0,1,0}, {0,0,1}, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
  outbreak_probs_vacaduladols[2] += prob_infect_pattern({0,0,1}, {0,1,0}, {0,0,1}, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);	
  
  /* one child, one adolescent, two adults, primary case child */
  hh = 1;
  for (i in 0 : 1) {
    for (j in 1 : 2) { // 0-2 as check that probs sum to 1; here focus on adult infections
	  outbreak_probs[3] += pdouble(j) * prob_infect_pattern({0,i,j}, {1,0,0}, {0,1,2}, transmission_rate/hh, ext_esc);	  
	  outbreak_probs_vacadults[3] += pdouble(j) * prob_infect_pattern({0,i,j}, {1,0,0}, {0,1,2}, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	  outbreak_probs_vacaduladols[3] += pdouble(j) * prob_infect_pattern({0,i,j}, {1,0,0}, {0,1,2}, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);	  
	}
  }
  
  /* one child, one adolescent, two adults, primary case adolescent */
  for (i in 0 : 1) {
    for (j in 1 : 2) { 
	  outbreak_probs[4] += pdouble(j) * prob_infect_pattern({i,0,j}, {0,1,0}, {1,0,2}, transmission_rate/hh, ext_esc);	  
	  outbreak_probs_vacadults[4] += pdouble(j) * prob_infect_pattern({i,0,j}, {0,1,0}, {1,0,2}, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	  outbreak_probs_vacaduladols[4] += pdouble(j) * prob_infect_pattern({i,0,j}, {0,1,0}, {1,0,2}, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);	  
	}
  }
  
  /* one child, one adolescent, two adults, primary case adult */
  for (i in 0 : 1) {
    for (j in 0 : 1) { 
	  outbreak_probs[5] += prob_infect_pattern({i,j,1}, {0,0,1}, {1,1,1}, transmission_rate/hh, ext_esc);	  
	  outbreak_probs_vacadults[5] += prob_infect_pattern({i,j,1}, {0,0,1}, {1,1,1}, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	  outbreak_probs_vacaduladols[5] += prob_infect_pattern({i,j,1}, {0,0,1}, {1,1,1}, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);	  
	}
  }  
    
  /* two children, two adolescents, two adults, primary case child */
  hh = 1;
  nn = {1,2,2};
  aa = {1,0,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 1 : nn[3] ) {
	    outbreak_probs[6] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_probs_vacadults[6] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_probs_vacaduladols[6] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }
    
  /* two children, two adolescents, two adults, primary case adolescent */
  nn = {2,1,2};
  aa = {0,1,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 1 : nn[3] ) {
	    outbreak_probs[7] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_probs_vacadults[7] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_probs_vacaduladols[7] += pdouble(k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }
    
  /* two children, two adolescents, two adults, primary case adult */
  nn = {2,2,1};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 1 : nn[3] ) { 
	    outbreak_probs[8] += prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_probs_vacadults[8] += prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_probs_vacaduladols[8] += prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  
  /* outbreak sizes - Table 3 in manuscript */
  for (i in 1 : size(outbreak_size) ) {
	outbreak_size[i] = 0;
	outbreak_size_vacadults[i] = 0;
	outbreak_size_vacaduladols[i] = 0;
  };
   
  /* fill in */
  hh = 1;

  ll = 1;
  nn = {0,0,1};
  aa = {1,0,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 2;
  nn = {1,0,0};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 3;
  nn = {0,0,1};
  aa = {0,1,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 4;
  nn = {0,1,0};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 5;
  nn = {1,0,2};
  aa = {1,0,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 6;
  nn = {2,0,1};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 7;
  nn = {0,1,2};
  aa = {0,1,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 8;
  nn = {0,2,1};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 9;
  nn = {1,2,2};
  aa = {1,0,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 10;
  nn = {2,1,2};
  aa = {0,1,0};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);

  ll = 11;
  nn = {2,2,1};
  aa = {0,0,1};
  for (i in 0 : nn[1]) { 
    for (j in 0 : nn[2]) { 
	  for (k in 0 : nn[3] ) { 
	    outbreak_size[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, transmission_rate/hh, ext_esc);	  
	    outbreak_size_vacadults[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_adults) .* transmission_rate/hh .* (1-VE_I_adults), ext_esc);	  
	    outbreak_size_vacaduladols[ll] += (i+j+k) * prob_infect_pattern({i,j,k}, aa, nn, (1-VE_S_aduladols) .* transmission_rate/hh .* (1-VE_I_aduladols), ext_esc);
      }		
	}
  }  
  outbreak_size[ll] /= sum(nn);
  outbreak_size_vacadults[ll] /= sum(nn);
  outbreak_size_vacaduladols[ll] /= sum(nn);
}

