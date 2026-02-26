# working environment ####

# cleaning R start
rm(list=ls()); gc()

# loading libraries
libraries = c( 
  # general purpose
  'here','tidyverse','rtables','gt','cowplot',
  
  # frequentist implementations
  'survival','cmprsk','tidycmprsk','ggsurvfit',
  
  # Bayesian implementations
  'cmdstanr', 'coda','posterior','bayesplot'
) 
# sapply(libraries, install.packages, character.only=T, dependencies=T)
sapply(libraries, require, character.only=T)




# Source data ####

# load data
rd  = file.path( here(), 'data', 'liver transplant example.dat')
ds = rd %>%
  read.table( header=T )

              
# data structure
str(ds)
# 
# Note:
#   patient: patient number
#   age: age at baseline
#   treatment: 1=Solution A, 2=Solution B
#   disease: primary liver disease, 1=PBC, 2=PSC, 3=ALD
#   time: time to cause of graft failure
#   status: 1=event, 0=censored
#   cof: cause of graft failure, 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other


# convert to factor
ds = ds %>%
  mutate(
    treatment_label = treatment,
    disease_label = disease,
    cof_label = cof ) %>%
  mutate_at( 
    .vars = vars(treatment_label), 
    .funs = factor,
    levels = 1:2, 
    labels = c('Solution A','Solution B')
  ) %>%
  mutate(
    treatment_label = 
      fct_relevel( treatment_label, 
                   c('Solution B','Solution A'))) %>%
  mutate_at( 
    .vars = vars(disease_label), 
    .funs = factor,
    levels = 1:3, 
    labels = c('PBC','PSC','ALD') ) %>%
  mutate(
    disease_label = 
      fct_relevel( disease_label, 
                   c('ALD','PBC','PSC'))) %>%
  mutate_at( 
    .vars = vars(cof_label), 
    .funs = factor,
    levels = 0:4, 
    labels = c('Functioning graft', 
               'Organ Rejection', 
               'Hepatic Artery Thrombosis', 
               'Recurrent Disease', 
               'Other') )




# ADTTE? ####

## Frequentist ####

# covariates
X = with(ds, model.matrix( ~ age + treatment_label + disease_label ))[,-1] 

## Bayesian ####

# modifying 'immediate failures'
ds$time[ds$time==0] = 1e-10

# defining status per cause: 
#   status=1 if cof==j; status=0 otherwise
status = c()
for( j in 1:max(ds$cof) ){ # 4 causes
  status = cbind(status, statusX=as.integer(ds$cof==j) )  
}
attr(status, "dimnames") = NULL

# static weight calculation
fit_cens = survfit( Surv(time, status==0) ~ 1, data=ds )
summary_fit = summary(fit_cens, times=ds$time, extend=T)
w = summary_fit$surv

# additional data for simulation
#   XS defines a patient with 55 years, male, and PSC
t = seq(1, max(ds$time), by=3)
XS = matrix( c(55,1,0,1), nrow=1) 


# data list
dlist = list( 
  # fitting data
  m = max(ds$cof),
  n = nrow(ds),
  q = ncol(X),
  t = ds$time,
  s = status,
  X = X,
  w = w,
  # for MC integration
  ms = 500,
  ns = length(t),
  ts = t,
  XS = XS )





# ARD? ####

## Frequentist ####

### Fine-Gray #####
# require(cmprsk)

for( j in 1:max(ds$cof)){
  
  # model fit
  fit_cof = with( 
    ds,
    cmprsk::crr( ftime=time, fstatus=cof, cov1=X, cencode=0, failcode=j) )
  saveRDS( fit_cof, file.path( here(),'results', paste0('fit_F_FG_',j,'.RDS') ) )
  
  # model summary
  summary_cof = fit_cof %>%
    summary()
  
  # Hazard ratios
  ARD_mom = with( summary_cof,
    c( exp(coef[,1]),
       exp(coef[,1]+qnorm(0.05)*coef[,3]),
       exp(coef[,1]+qnorm(0.95)*coef[,3]),
       c( NA, pchisq( q=logtest[1], df=logtest[2] ), rep(NA,2) ) )
    ) %>%
    matrix( 
      ncol=4, 
      byrow=F,
      dimnames = list(
        NULL,
        c('Hazard ratio', 'lower', 'upper', 'p-value'))
      ) %>%
    as_tibble() %>%
    mutate(
      'Effect/Covariate Included in the model' = 
        c('Age', 
          'Solution A vs control (Solution B)', 
          'Disease PBC vs ALD', 
          'Disease PSC vs ALD')) %>%
    arrange( match(`Effect/Covariate Included in the model`, 
                   c( 'Solution A vs control (Solution B)',
                      'Age', 
                      'Disease PBC vs ALD', 
                      'Disease PSC vs ALD'))) %>%
    mutate( cof = j ) %>%
    mutate( cof_label = cof ) %>% 
    mutate_at( 
      .vars = vars(cof_label), 
      .funs = factor,
      levels = 0:4, 
      labels = c('Functioning graft', 
                 'Organ Rejection', 
                 'Hepatic Artery Thrombosis', 
                 'Recurrent Disease', 
                 'Other') )
  
  # save all causes
  if( j==1 ){
    ARD_cof = ARD_mom
  } else{
    ARD_cof = ARD_cof %>%
      add_row(ARD_mom)
  }
}

# save ARD object
saveRDS( ARD_cof, file.path( here(),'results', paste0('ARD_F_FG.RDS') ) )




## Bayesian ####

### Fine-Gray #####

# stan model
StanModel = "
data {
  int<lower=0> m; // number of causes
  int<lower=0> n; // num obs
  int<lower=0> q; // num covariates
  vector[n] t;    // time
  matrix[n,m] s;  // cause of failure
  matrix[n,q] X;  // covariates
  vector[n] w;    // weights
}
parameters {
  matrix[q,m] beta; // beta[covariate,cause]
}
model {
  // priors
  to_vector(beta) ~ normal(0,2);
  
  // risk score
  matrix[n,m] rs = X*beta;
  
  // partial log-likelihood
  for( j in 1:m ){
    for( i in 1:n ) {
      if( s[i,j]==1 ){
        real log_den = log_sum_exp( log(w[1:i]) + rs[1:i, j] );
        target += rs[i,j] - log_den;
      }
    }
  }
}
"

# saving model
model_nam = 'Bayesian_Fine-Gray.stan'
writeLines( StanModel, con=file.path( here(),'code', model_nam) )

# loading model
set_cmdstan_path( '/home/josema/.cmdstan/cmdstan-2.37.0')
mod = cmdstan_model( file.path( here(),'code', model_nam) )

# sampling parameters
fit = mod$sample( 
  data = dlist, 
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  chains = 4, 
  parallel_chains = 4,
  seed = 12345 ) 

# saving parameter samples
fit$save_object( 
  file = file.path( 
    here(),
    'results', 
    str_replace(model_nam,'.stan','.RDS') ) )


# load parameter samples
fit = readRDS( 
  file = file.path( here(),'results', str_replace(model_nam,'.stan','.RDS')) )

# parameter posterior distributions
post = as_draws_df( fit ) 

# hazard ratio posterior distribution
HR = post %>%
  subset_draws( variable='beta') %>%
  exp()

# hazard ratio estimate summaries
HRsumm = HR %>%
  summarize_draws( ) 

# hazard ratio "p-value" estimate summary
thr = c(0.95, 1.1) # rope
p_value = apply( HR>thr[1] & HR<thr[2], 2, as.integer) %>%
  as_tibble() %>%
  summarize_draws( ) %>%
  select( variable, mean ) %>%
  rename( 'P(0.95< HR <1.1)' = mean )

# join hazard ratio
by = join_by( variable )
HR_join = left_join( HRsumm, p_value, by )


# create ARD object
ARD_cof = HR_join %>%
  mutate(
    cof = str_extract( variable, "[:digit:][:punct:]$" )
  ) %>%
  mutate(
    cof = str_extract( cof, "[:digit:]" )
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]1[:punct:][:digit:][:punct:]$'),
      'Age' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]2[:punct:][:digit:][:punct:]$'),
      'Solution A vs control (Solution B)' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]3[:punct:][:digit:][:punct:]$'),
      'Disease PBC vs ALD' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]4[:punct:][:digit:][:punct:]$'),
      'Disease PSC vs ALD' ) 
  ) %>%
  arrange( 
    cof,
    match( variable, 
           c( 'Solution A vs control (Solution B)',
              'Age', 
              'Disease PBC vs ALD', 
              'Disease PSC vs ALD'))) %>%
  rename( 
    'Effect/Covariate Included in the model' = variable,
    'Hazard ratio' = mean,
    lower = q5,
    upper = q95
  ) %>%
  mutate( cof_label = cof ) %>% 
  mutate_at( 
    .vars = vars(cof_label), 
    .funs = factor,
    levels = 0:4, 
    labels = c('Functioning graft', 
               'Organ Rejection', 
               'Hepatic Artery Thrombosis', 
               'Recurrent Disease', 
               'Other') )

# save ARD object
saveRDS( 
  ARD_cof, 
  file.path( 
    here(),
    'results', 
    paste0('ARD_B_FG.RDS') ) )




### Jeong-Fine #####

# stan model
StanModel = "
functions{
  // cause-specific log hazard function
  real log_hij( real t, real rs, real lambda, real gamma ){
    return( rs + log(lambda) + log(gamma) + (gamma-1)*log(t) );
  }
  // cause-specific cummulative hazard function
  real Hij( real t, real rs, real lambda, real gamma ){
    return( exp( rs ) * lambda * t^gamma ); 
  }
  // overall log-survival function
  real log_S( real t, vector rs, vector lambda, vector gamma, int m ) {
    real oSt = 0;
    for (j in 1:m) {
      oSt += Hij(t, rs[j], lambda[j], gamma[j]);
    }
    return(-oSt);
  }
  // cause specific log-density
  real log_fij( real t, int j, vector rs, vector lambda, vector gamma, int m ){
    return ( log_hij( t, rs[j], lambda[j], gamma[j] ) + 
             log_S( t, rs, lambda, gamma, m ) );
  }
  // Integrand for the density function
  real fij(real t, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    // parameters and data
    int j = x_i[1];
    int m = x_i[2];
    vector[m] rs = to_vector( theta[1:m] );
    vector[m] lambda = to_vector( theta[(m+1):(2*m)] );
    vector[m] gamma =  to_vector( theta[(2*m+1):(3*m)] );
    return( exp( log_fij( t, j, rs, lambda, gamma, m ) ) );
  }
}
data {
  // fitting data
  int<lower=0> m;       // number of causes
  int<lower=0> n;       // num obs
  int<lower=0> q;       // num covariates
  vector[n] t;          // time
  matrix[n,m] s;        // cause of failure
  matrix[n,q] X;        // covariates
  
  // For Marginal Approximation
  int<lower=1> ms;      // Number of Monte Carlo draws (e.g., 100)
  int<lower=0> ns;      // number of time samples
  vector[ns] ts;         // time samples
  matrix[1,q] XS;       // covariates for simulation
}
parameters {
  matrix[q,m] beta;           // beta[covariate,cause]
  vector<lower=0>[m] lambda;  // lambda[cause]
  vector<lower=0>[m] gamma;   // gamma[cause]
}
model {
  // priors
  to_vector(beta) ~ normal(0,2);
  lambda ~ exponential(2);
  gamma ~ gamma(8,5); 
  
  // risk score
  matrix[n,m] rs = X*beta;
  
  // log-likelihood
  for( j in 1:m ){
    for( i in 1:n ){
      target += s[i,j] * log_hij( t[i], rs[i,j], lambda[j], gamma[j] ) -
                Hij( t[i], rs[i,j], lambda[j], gamma[j] );
    }
  }
}
generated quantities {
  matrix[ns,m] Fij_MC;
  matrix[ns,m] Fij_I;
  
  // risk score
  vector[m] rs = to_vector(XS*beta);
  
  // flatten parameters 
  array[3*m] real flat_theta;
  for (j in 1:m) {
    flat_theta[j] = rs[j];
    flat_theta[m+j] = lambda[j];
    flat_theta[2*m+j] = gamma[j];
  }
  
  for (i in 1:ns) {
    for (j in 1:m) {
    
      // Monte Carlo integration
      real fij_sum = 0;
      for (k in 1:ms) {
        real u = uniform_rng(0, ts[i]);
        fij_sum += exp(log_fij(u, j, rs, lambda, gamma, m));
      }
      Fij_MC[i,j] = ts[i] * (fij_sum / ms); 
      // Fij = t[i] * p(t[i])
      
      // Double Exponential Quadrature integration
      Fij_I[i,j] = integrate_1d( fij, 0.0, ts[i], flat_theta, {0.0}, {j,m} );
      
    }
  }
}
"

# saving
model_nam = 'Bayesian_Jeong-Fine.stan'
writeLines( StanModel, con=file.path( here(),'code', model_nam) )

# loading model
set_cmdstan_path( '/home/josema/.cmdstan/cmdstan-2.37.0')
mod = cmdstan_model( file.path( here(),'code', model_nam) )

# sampling parameters
fit = mod$sample( 
  data = dlist, 
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  chains = 4, 
  parallel_chains = 4, 
  seed = 28947 )

# saving parameter samples
fit$save_object( 
  file = file.path( 
    here(),
    'results', 
    str_replace(model_nam,'.stan','.RDS') ) )


# load parameter samples
fit = readRDS( 
  file = file.path( here(),'results', str_replace(model_nam,'.stan','.RDS')) )

# parameter posterior distributions
post = as_draws_df( fit ) 

# hazard ratio posterior distribution
HR = post %>%
  subset_draws( variable='beta') %>%
  exp()

# hazard ratio estimate summaries
HRsumm = HR %>%
  summarize_draws( ) 

# hazard ratio "p-value" estimate summary
thr = c(0.95, 1.1) # rope
p_value = apply( HR>thr[1] & HR<thr[2], 2, as.integer) %>%
  as_tibble() %>%
  summarize_draws( ) %>%
  select( variable, mean ) %>%
  rename( 'P(0.95< HR <1.1)' = mean )

# join hazard ratio
by = join_by( variable )
HR_join = left_join( HRsumm, p_value, by )


# theta posterior distribution
theta = post %>%
  subset_draws( variable=c('lambda','gamma')) %>% 
  summarize_draws( )


# create ARD object
ARD_cof = HR_join %>%
  add_row(theta) %>%
  mutate(
    cof = str_extract( variable, "[:digit:][:punct:]$" )
  ) %>%
  mutate(
    cof = str_extract( cof, "[:digit:]" )
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]1[:punct:][:digit:][:punct:]$'),
      'Age' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]2[:punct:][:digit:][:punct:]$'),
      'Solution A vs control (Solution B)' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]3[:punct:][:digit:][:punct:]$'),
      'Disease PBC vs ALD' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, '[:punct:]4[:punct:][:digit:][:punct:]$'),
      'Disease PSC vs ALD' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, 'lambda'),
      'Scale' ) 
  ) %>%
  mutate(
    variable = replace( 
      variable, 
      str_detect(variable, 'gamma'),
      'Shape' ) 
  ) %>%
  arrange( 
    cof,
    match( variable, 
           c( 'Solution A vs control (Solution B)',
              'Age', 
              'Disease PBC vs ALD', 
              'Disease PSC vs ALD',
              'Scale',
              'Shape'))) %>%
  rename( 
    'Effect/Covariate Included in the model' = variable,
    'Hazard ratio' = mean,
    lower = q5,
    upper = q95
  ) %>%
  mutate( cof_label = cof ) %>% 
  mutate_at( 
    .vars = vars(cof_label), 
    .funs = factor,
    levels = 0:4, 
    labels = c('Functioning graft', 
               'Organ Rejection', 
               'Hepatic Artery Thrombosis', 
               'Recurrent Disease', 
               'Other') )

# save ARD object
saveRDS( 
  ARD_cof, 
  file.path( 
    here(),
    'results', 
    paste0('ARD_B_JF.RDS') ) )



# summaries
Fij_MC = summarize_draws( subset_draws( post, variable='Fij_MC' ) ) 
saveRDS( 
  Fij_MC, 
  file.path( 
    here(),
    'results', 
    paste0('ARD_B_JF_cumincMC.RDS') ) )

Fij_I = summarize_draws( subset_draws( post, variable='Fij_I' ) ) 
saveRDS( 
  Fij_I, 
  file.path( 
    here(),
    'results', 
    paste0('ARD_B_JF_cumincI.RDS') ) )





# TLGs/TLFs ####

## Demographics (DMT01) #####

# calculation and formatting function
f_summary = function(x) {
  if ( is.numeric(x) ) {
    in_rows(
      "n" = rcell( sum(!is.na(x)), format = "xx"),
      "Mean (SD)" = rcell( c(mean(x, na.rm=T), sd(x, na.rm=T)), format = "xx.x (xx.x)"),
      "Median" = rcell( median(x, na.rm=T), format = "xx.x"),
      "Min - Max" = rcell( range(x, na.rm=T), format = "xx.x - xx.x")
    )
  } else if ( is.factor(x) ) {
    tab = c( table(x) )
    do.call( 
      in_rows, 
      c( "n" = rcell( sum(tab), format = "xx" ),
         lapply(tab, rcell, format = "xx" ) ) )
  } else {
    stop("type not supported")
  }
}

# table
basic_table( show_colcounts=T ) %>%
  split_cols_by( 
    var = "treatment_label",
    ref_group = "Solution B") %>%
  add_overall_col("All Patients") %>%
  analyze( 
    vars = c('age','disease_label'), 
    afun = f_summary,
    var_labels = c('Age (yr)','Liver Disease') 
  ) %>%
  build_table( df=ds )



## Deaths (DTHT01) #####

# calculation and formatting functions
event_count = function(x, .N_col) {
  tab = c( table(x) )
  in_rows(
    "Total number of events" = rcell( tab[2] * c(1, 1 / .N_col), format = "xx (xx.x%)" )
  )
}

cof_count = function(x, .N_col) {
  tab = c( table(x) )
  tab = tab[ names(tab) != "Functioning graft" ]
  do.call( 
    in_rows, 
    c( "n" = rcell( sum(tab), format = "xx" ),
       lapply( 
         tab,
         function(xi){ rcell(xi * c(1, 1 / .N_col), format = "xx (xx.x%)") } 
       ) ) )
}


# table
basic_table( show_colcounts=T ) %>%
  split_cols_by( 
    var = "treatment_label", 
    ref_group = "Solution B") %>%
  add_overall_col("All Patients") %>%
  analyze( 
    vars = "status", 
    afun = event_count, 
    show_labels = "hidden" ) %>%
  analyze( 
    vars = "cof_label", 
    afun = cof_count, 
    var_labels = "Cause of graft failure", 
    show_labels = "visible" ) %>%
  build_table( df=ds )




## Frequentist ####

### Fine-Gray regression ####
# (Equivalence: COXT01) 

# load ADR?
hr_cof = readRDS( 
  file.path( 
    here(),
    'results',
    paste0('ARD_F_FG.RDS') ) )


# Hazard tables
hr_cof[,c(5,1:4,7)] %>%
  filter( cof_label == 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','p-value')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )


# Hazard tables
hr_cof[,c(5,1:4,7)] %>%
  filter( cof_label != 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','p-value')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )



### Fine-Gray failures plot ####
# (Equivalence: KMG01)

# require(ggsurvfit)

cof = c('Organ Rejection', 
        'Hepatic Artery Thrombosis', 
        'Recurrent Disease', 
        'Other')

for( j in 1:max(ds$cof) ){
  assign( 
    paste0('F', j),
    cuminc( 
      Surv(time, cof_label) ~ treatment_label, 
      data = ds ) %>%
      ggcuminc( outcome = cof[j] ) +
      add_confidence_interval() +
      add_risktable() +
      scale_ggsurvfit( 
        x_scales = list(breaks=seq(0,max(ds$time), by=500)),
        y_scales = list(limits=c(0,0.07)) )
  )
}

F1

plot_grid( 
  F1, F2, F3, F4,
  ncol=2, nrow=2,
  labels = cof, label_size=10, label_x=0.15, label_y=0.95)




## Bayesian ####

### Fine-Gray regression ####
# (Equivalence: COXT01) 

# load ADR?
hr_cof = readRDS( 
  file.path( 
    here(),
    'results',
    paste0('ARD_B_FG.RDS') ) )


hr_cof %>%
  select( 
    'Effect/Covariate Included in the model', 
    'Hazard ratio', 
    lower, 
    upper, 
    'P(0.95< HR <1.1)', 
    cof_label ) %>%
  filter( cof_label == 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','P(0.95< HR <1.1)')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )


hr_cof %>%
  select( 
    'Effect/Covariate Included in the model', 
    'Hazard ratio', 
    lower, 
    upper, 
    'P(0.95< HR <1.1)', 
    cof_label ) %>%
  filter( cof_label != 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','P(0.95< HR <1.1)')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )


### Jeong-Fine regression ####
# (Equivalence: COXT01) 

# load ADR?
hr_cof = readRDS( 
  file.path( 
    here(),
    'results',
    paste0('ARD_B_JF.RDS') ) )


hr_cof %>%
  select( 
    'Effect/Covariate Included in the model', 
    'Hazard ratio', 
    lower, 
    upper, 
    'P(0.95< HR <1.1)', 
    cof_label ) %>%
  filter( cof_label == 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','P(0.95< HR <1.1)')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )


hr_cof %>%
  select( 
    'Effect/Covariate Included in the model', 
    'Hazard ratio', 
    lower, 
    upper, 
    'P(0.95< HR <1.1)', 
    cof_label ) %>%
  filter( cof_label != 'Organ Rejection' ) %>%
  group_by( cof_label ) %>%
  gt() %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  tab_spanner(
    label = "90% CI",
    columns = c('lower','upper')
  ) %>%
  tab_spanner(
    label = "Treatment effect adjusted for covariate",
    columns = c('Hazard ratio','lower','upper','P(0.95< HR <1.1)')
  ) %>%
  fmt_number( 
    columns = 2:4,
    decimal = 2
  ) %>%
  fmt_number( 
    columns = 5,
    decimal = 4
  ) %>%
  sub_missing( missing_text = " " )



### Jeong-Fine failures plot ####
# (Equivalence: KMG01)

# load ADR?
Fij_MC = readRDS( 
  file.path( 
    here(),
    'results',
    paste0('ARD_B_JF_cumincMC.RDS') ) )

Fij_I = readRDS( 
  file.path( 
    here(),
    'results',
    paste0('ARD_B_JF_cumincI.RDS') ) )

# save additional ADRs
t = seq(1, max(ds$time), by=3)
idx = str_detect( Fij_I$variable, '1[:punct:]$' )
Fij_I_mom = Fij_I[idx,]



# plot 1
plot( t, Fij_I_mom$mean, type='l', lwd=4, lty=1, ylim=c(0,0.07),
      ylab = "Cumulative incidence" ,
      xlab = "Time" )
for(i in 1:length(t)){
  lines( 
    x = rep(t[i],2), 
    y = with(Fij_I_mom[i,], c(q5, q95) ), 
    col=rgb(0,0,0,0.1) )  
}
legend( 'topleft', bty='n', lty=c(1,1,1),  
        col=c('white','black', rgb(0,0,0,0.3)),
        legend=c( paste0('P(C=',1,') = ', round( tail(Fij_I$mean[idx], 1), 2) ),
                  'Mean incidence', 
                  '90% CI') )


# plot 2
par(mfrow=c(2,2))
for( j in 1:4 ){
  idx = str_detect( Fij_MC$variable, paste0(j,'[:punct:]$') )
  plot( t, Fij_I$mean[idx], type='l', lwd=4, lty=1, ylim=c(0,0.15))
  lines( t, Fij_MC$mean[idx], col='red', lty=2)
  legend( 'topleft', bty='n', col=c('white','black','red'), lty=c(1,1,2),
          legend=c( paste0('P(C=',j,') = ', round( tail(Fij_I$mean[idx], 1), 2) ),
                    'Double Exponential Quadrature', 
                    'Monte Carlo approx') )
}
par(mfrow=c(1,1))
