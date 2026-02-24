## example 12.05, Cox cause-specific hazards ####
dLTrans = read.table( 
  file.path(dir, 'Collett (2023) data', "Survival of liver transplant recipients.dat"), 
  header=T)
# str(dLTrans)
# with(dLTrans, table(cof,status,disease))
# 
# Note:
#   gender: 1=male, 2=female
#   disease: 1=PBC, 2=PSC, 3=ALD
#   status : 1=event, 0=censored
#   cof: 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other

# data transform
dLTrans$time = as.numeric(dLTrans$time)
dLTrans$gender = factor(dLTrans$gender)
dLTrans$gender = relevel(dLTrans$gender, ref=2)
dLTrans$disease = factor(dLTrans$disease)
dLTrans$disease = relevel(dLTrans$disease, ref=3)

# Cox model
for( j in 1:max(dLTrans$cof) ){ # 4 causes
  # define status
  statusC = as.integer(dLTrans$cof==j)  
  
  # Cox cause-specific hazards model
  assign( 
    paste0('summary_cox_cof', j),
    summary( coxph( Surv(time, statusC) ~ age+gender+disease, data=dLTrans) ) )
  
  summ_cof = get(paste0('summary_cox_cof', j))$conf.int[,c(1,3:4)]
  summ_cof = round(summ_cof, 2)
  
  print( paste0('cof == ', j) )
  print( summ_cof )
}
# pred_cof1 = t( predict.crr(fit_cof1, cov1=xvars) )
# str(pred_cof1[,1])
# 
# Note:
#   status=1 if cof==j
#   status=0 otherwise





## example 12.XX, Weibull cause-specific hazards ####
dLTrans = read.table( 
  file.path(dir, 'Collett (2023) data', "Survival of liver transplant recipients.dat"), 
  header=T)
# str(dLTrans)
# with(dLTrans, table(cof,status,disease))
# 
# Note:
#   gender: 1=male, 2=female
#   disease: 1=PBC, 2=PSC, 3=ALD
#   status : 1=event, 0=censored
#   cof: 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other

# data transform
dLTrans$time = as.numeric(dLTrans$time)
dLTrans$gender = factor(dLTrans$gender)
dLTrans$gender = relevel(dLTrans$gender, ref=2)
dLTrans$disease = factor(dLTrans$disease)
dLTrans$disease = relevel(dLTrans$disease, ref=3)

# changing time=0
dLTrans$time_mod = dLTrans$time
dLTrans$time_mod[dLTrans$time_mod==0] = 1e-10


# Weibull model
for( j in 1:max(dLTrans$cof) ){ # 4 causes
  # define status
  statusC = as.integer(dLTrans$cof==j)  
  
  # Cox cause-specific hazards model
  assign( 
    paste0('fit_wei_cof', j),
    survreg( Surv(time_mod, statusC) ~ age+gender+disease, 
             dist='weibull', data=dLTrans) )
  
  assign( paste0('summary_wei_cof', j), 
          summary( get(paste0('fit_wei_cof', j)) ) )
  
  mean_g = with( get(paste0('summary_wei_cof', j)), 1/exp(table[6,1]) ) # equivalently: 1/scale
  lower_g = with( get(paste0('summary_wei_cof', j)), 1/exp(table[6,1] + qnorm(0.975)*table[6,2]))
  upper_g = with( get(paste0('summary_wei_cof', j)), 1/exp(table[6,1] + qnorm(0.025)*table[6,2]))
  summ_g = data.frame( mean_g, lower_g, upper_g)
  
  mean_cof = exp(get(paste0('summary_wei_cof', j))$table[-c(1,6),1])^(-mean_g)
  lower_cof = with( get(paste0('summary_wei_cof', j)), exp(table[-c(1,6),1] + qnorm(0.025)*table[-c(1,6),2])^(-mean_g) )
  upper_cof = with( get(paste0('summary_wei_cof', j)), exp(table[-c(1,6),1] + qnorm(0.975)*table[-c(1,6),2])^(-mean_g) )
  summ_cof = data.frame(mean_cof, lower_cof, upper_cof)
  
  mean_l = with( get(paste0('summary_wei_cof', j)), exp(table[1,1])^(-mean_g) )
  lower_l = with( get(paste0('summary_wei_cof', j)), exp(table[1,1] + qnorm(0.975)*table[1,2])^(-mean_g) )
  upper_l = with( get(paste0('summary_wei_cof', j)), exp(table[1,1] + qnorm(0.025)*table[1,2])^(-mean_g) )
  summ_l = data.frame( mean_l, lower_l, upper_l)
  
  print( paste0('cof == ', j) )
  print( round(summ_l, 7) )
  print( round(summ_g, 2) )
  print( round(summ_cof, 2) )
}
# 
# Note:
#   intercepts are very high because they represent very low lambda0's





## example 12.06, Fine-Gray cause-specific incidence ####
dLTrans = read.table( 
  file.path(dir, 'Collett (2023) data', "Survival of liver transplant recipients.dat"), 
  header=T)
# str(dLTrans)
# with(dLTrans, table(cof,status,disease))
# 
# Note:
#   gender: 1=male, 2=female
#   disease: 1=PBC, 2=PSC, 3=ALD
#   status : 1=event, 0=censored
#   cof: 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other

# data transform
dLTrans$gender = factor(dLTrans$gender)
dLTrans$gender = relevel(dLTrans$gender, ref=2)
dLTrans$disease = factor(dLTrans$disease)
dLTrans$disease = relevel(dLTrans$disease, ref=3)

# covariates
X = with(dLTrans, model.matrix(~age+gender+disease))[,-1] # covariates

# model
# require(cmprsk)
for( i in 1:max(dLTrans$cof)){
  assign( 
    paste0('summary_cof', i),
    with( dLTrans, 
          summary(crr( ftime=time, fstatus=cof, cov1=X, cencode=0, failcode=i) ) )
  )
  
  mean_cof = with(get(paste0('summary_cof', i)), exp(coef[,1]) )
  lower_cof = with(get(paste0('summary_cof', i)), exp(coef[,1]+qnorm(0.025)*coef[,3]) )
  upper_cof = with(get(paste0('summary_cof', i)), exp(coef[,1]+qnorm(0.975)*coef[,3]) )
  summ_cof = data.frame( mean_cof, lower_cof, upper_cof) 
  
  print( paste0('cof == ', i) )
  print( round(summ_cof,2) )
}

# pred_cof1 = t( predict.crr(fit_cof1, cov1=X) )
# str(pred_cof1[,1])




### Cause-specific Incidence (Fine-Gray models) ####

#### example 12.06, Cox (static weights) ####

# data
dLTrans = read.table( 
  file.path(dir, 'Collett (2023) data', "Survival of liver transplant recipients.dat"), 
  header=T)
dLTrans = dLTrans[ order(-dLTrans$time), ]
# 
# Note:
#   gender: 1=male, 2=female
#   disease: 1=PBC, 2=PSC, 3=ALD
#   status : 1=event, 0=censored
#   cof: 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other

# data transform
dLTrans$gender = factor(dLTrans$gender)
dLTrans$gender = relevel(dLTrans$gender, ref=2)
dLTrans$disease = factor(dLTrans$disease)
dLTrans$disease = relevel(dLTrans$disease, ref=3)

# covariates
X = with(dLTrans, model.matrix(~age+gender+disease))[,-1] # covariates
attr(X, "dimnames") = NULL

# status
status = c()
for( j in 1:max(dLTrans$cof) ){ # 4 causes
  status = cbind(status, statusX=as.integer(dLTrans$cof==j) )  
}
attr(status, "dimnames") = NULL
# 
# Note:
#   status=1 if cof==j
#   status=0 otherwise

# static weights
fit_cens = survfit( Surv(time, status==0) ~ 1, data=dLTrans )
# with(dLTrans, table(status==0, cof))
# 
# Note:
#   I am estimating the survival for all censored observations, because I use
#   as an event indicator (status==0), which means now that:
#   status=1 if censored and status=0 if event from any cof.
#   Notably, the variable status used here is the one that comes with the 
#   original data.

# Find the probability of NOT being censored at each time point
summary_fit = summary(fit_cens, times=dLTrans$time, extend=T)
w = summary_fit$surv
# all(dLTrans$time == summary_fit$time)
# 
# Note:
#   For Fine-Gray, the weight for subject i at time t is:
#   w_i = S(t) / S( min(t_i, t) )
#   For a parametric model, we often use the S(t) at the observed time


# data list
dlist = list( 
  m = max(dLTrans$cof),
  n = nrow(dLTrans),
  q = ncol(X),
  t = dLTrans$time,
  s = status,
  X = X,
  w = w)
# str(dlist)


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
  to_vector(beta) ~ normal(0, 2);
  
  // risk score
  matrix[n,m] rs = X*beta;
  
  for( j in 1:m ){
    for( i in 1:n ) {
      if( s[i,j]==1 ){
        // risk_set_rs = rs[1:i, j]
        // risk_set_w = w[1:i]
        // log_den = log_sum_exp( log(risk_set_w) + risk_set_rs )
        //         = log_sum_exp( log(w[1:i]) + rs[1:i, j] )
        real log_den = log_sum_exp( log(w[1:i]) + rs[1:i, j] );
        target += rs[i,j] - log_den;
      }
    }
  }
}
generated quantities{
  // Here, following Equation in Collett (2023) pp. 366, I could use the
  // number of deaths at each time to approximate the cause-specific 
  // cumulative hazard subdistribution (Hij), and then convert it to the
  // cause specific cumulative incidence function (Fij). However, there is 
  // NO new learning outcome from doing this because it does not involve the 
  // use of any probability distribution. Moreover, this is easier to do in 
  // the parametric model (see below) 
}
"

# saving
model_nam = '12.06_Cox_Fine_and_Gray_CSFj_CR.stan'
writeLines(StanModel, con=file.path(dir,'Collett (2023) code','models', model_nam))

# loading
set_cmdstan_path( '/home/josema/.cmdstan/cmdstan-2.37.0')
mod = cmdstan_model( file.path( dir,'Collett (2023) code','models', model_nam), ) # to create the C++ model

fit = mod$sample( 
  data = dlist, 
  iter_warmup = 1000, iter_sampling = 1000, 
  chains = 4, parallel_chains = 4) 
# adapt_delta=0.95, max_treedepth=20 ) #, init=0
fit$save_object( 
  file = file.path( dir,'Collett (2023) code','models', str_replace(model_nam,'.stan','.RDS')) )


# load model
fit = readRDS( 
  file = file.path( dir,'Collett (2023) code','models', str_replace(model_nam,'.stan','.RDS')) )
# fit$summary()

# parameter recovery
post = as_draws_df( fit )
# post = subset_draws( post, variable='lambda' )

summarize_draws( exp(subset_draws( post, variable='beta' )) )
# 
# Frequentist estimates
# [1] "cof == 1"
#              mean_cof lower_cof upper_cof
# age          0.97      0.93      1.00
# gender1      0.60      0.27      1.31
# disease1     0.69      0.25      1.91
# disease2     1.10      0.44      2.76
# [1] "cof == 2"
#              mean_cof lower_cof upper_cof
# age          0.98      0.96      1.00
# gender1      1.08      0.67      1.75
# disease1     1.83      1.06      3.17
# disease2     1.59      0.87      2.91
# [1] "cof == 3"
#              mean_cof lower_cof upper_cof
# age          0.98      0.95      1.01
# gender1      0.89      0.45      1.74
# disease1     0.35      0.13      0.98
# disease2     1.68      0.84      3.39
# [1] "cof == 4"
#              mean_cof lower_cof upper_cof
# age          1.03      1.01      1.05
# gender1      0.91      0.59      1.42
# disease1     0.65      0.38      1.12
# disease2     1.30      0.83      2.04
# 
# Note:
#   All cause-specific HR are quite similar to the frequentist HRj estimates 
#   and to the ones reported in the book.




#### example 12.06, Weibull cause-specific hazard approach with CIF approx ####

# data
dLTrans = read.table( 
  file.path(dir, 'Collett (2023) data', "Survival of liver transplant recipients.dat"), 
  header=T)
dLTrans = dLTrans[ order(-dLTrans$time), ]
# 
# Note:
#   gender: 1=male, 2=female
#   disease: 1=PBC, 2=PSC, 3=ALD
#   status : 1=event, 0=censored
#   cof: 0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other

# data transform
dLTrans$gender = factor(dLTrans$gender)
dLTrans$gender = relevel(dLTrans$gender, ref=2)
dLTrans$disease = factor(dLTrans$disease)
dLTrans$disease = relevel(dLTrans$disease, ref=3)

X = with(dLTrans, model.matrix(~age+gender+disease))[,-1] # covariates
attr(X, "dimnames") = NULL

status = c()
for( j in 1:max(dLTrans$cof) ){ # 4 causes
  status = cbind(status, statusX=as.integer(dLTrans$cof==j) )  
}
attr(status, "dimnames") = NULL

# IMPORTANT: 
#   The data set contains 3 zeros to show that some individuals had an 
#   'immediate failure'. However, the Weibull log_hij need to calculate the
#   log(t), which would evaluate to -Inf when t=0. Thus, for those cases we 
#   add a small constant to the data.
dLTrans$time[dLTrans$time==0] = 1e-10
# with( dLTrans, table(time==0, cof) )


# for simulation
t = seq(1, max(dLTrans$time),by=3)
XS = matrix( c(55,1,0,1), nrow=1) # patient with 55 years, male, and PSC

# data list
dlist = list( 
  # fitting data
  m = max(dLTrans$cof),
  n = nrow(dLTrans),
  q = ncol(X),
  t = dLTrans$time,
  s = status,
  X = X,
  # for MC approx. of integral
  ms = 500,
  ns = length(t),
  ts = t,
  xs = XS )
# str(dlist)


# stan model
StanModel = "
functions{
  // Based on Jeong et al (2006) Equations (9), (10), and (13)
  // see also example 12.05, Weibull 
  // cause-specific log hazard function
  //  hij(t) = exp(X*b) * hOj(t)
  //  log[hij(t)] = X*b + log[hOj(t)]
  //              = rs + log[hOj(t)]
  //    where X*b = rs = risk_score
  //  hOj(t) = lambda_j * gamma_j * t^{gamma_j-1}
  //  log[hOj(t)] = log(lambda_j) + log(gamma_j) + (gamma_j-1)*log(t)
  real log_hij( real t, real rs, real lambda, real gamma ){
    return( rs + log(lambda) + log(gamma) + (gamma-1)*log(t) );
  }
  // cause-specific cummulative hazard function
  //  Hij(t) = int_{0}^{t} hij(t)
  //         = int_{0}^{t} exp(X*b) * lambda_j * gamma_j * u^{gamma_j-1} du
  //         = exp(X*b) * lambda_j * [ int_{0}^{t} gamma_j * u^{gamma_j-1} du ]
  //         = exp(X*b) * lambda_j * u^gamma_j |_{0}^{t}
  //         = exp(X*b) * lambda_j * t^gamma_j
  real Hij( real t, real rs, real lambda, real gamma ){
    return( exp( rs ) * lambda * t^gamma ); 
  }
  // overall log-survival function
  //  S(t) = exp[ -sum_{j=1}^{m} Hij(t) ]
  //  log[S(t)] = -sum_{j=1}^{m} Hij(t)
  real log_S( real t, vector rs, vector lambda, vector gamma, int m ) {
    real oSt = 0;
    for (j in 1:m) {
      oSt += Hij(t, rs[j], lambda[j], gamma[j]);
    }
    return(-oSt);
  }
  //
  // Note:
  // Jeong et al (2006) uses the cause-specific pseudo-survival functions Sk(t)
  // but because: Sk(t) = exp[-Hk(t)], the log[Sk(t)] = -Hk(t) in Equation (13)
  //
  // cause specific log-density
  // fij(t) = hij(t) * S(t)
  // log[fij(t)] = log[hij(t)] + log[S(t)]
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
  matrix[1,q] xs;       // covariates for simulation
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
  
  // Pre-calculate for efficiency
  matrix[n,m] rs = X*beta;
  
  // model
  // Following Equation (12.9) from Collett (2023)
  for( j in 1:m ){
    for( i in 1:n ){
      target += s[i,j] * log_hij( t[i], rs[i,j], lambda[j], gamma[j] ) -
                Hij( t[i], rs[i,j], lambda[j], gamma[j] );
    }
  }
}
generated quantities {
  matrix[ns,m] Fij_MC; // Cumulative Incidence Function for 1st cause
  matrix[ns,m] Fij_I;
  
  // risk score
  vector[m] rs = to_vector(xs*beta);
  
  // flatten all parameters into a single array. Size: 3*m = m(rs) + m(lambda) + m(gamma) 
  array[3*m] real flat_theta;
  for (j in 1:m) {
    flat_theta[j] = rs[j];
    flat_theta[m+j] = lambda[j];
    flat_theta[2*m+j] = gamma[j];
  }
  
  for (i in 1:ns) {
    for (j in 1:m) { // activate for all causes
    
      // Integrate using Monte Carlo approximation
      real fij_sum = 0;
      for (k in 1:ms) {
        real u = uniform_rng(0, ts[i]);
        fij_sum += exp(log_fij(u, j, rs, lambda, gamma, m));
      }
      Fij_MC[i,j] = ts[i] * (fij_sum / ms); 
      // Fij = t[i] * p(t[i])
      
      // Integrate using integrate_1d function
      Fij_I[i,j] = integrate_1d( fij, 0.0, ts[i], flat_theta, {0.0}, {j,m} );
      
    } // activate for all causes
  }
}
"

# saving
model_nam = '12.06_Weibull_Jeong_and_Fine_CShk_CR.stan'
writeLines(StanModel, con=file.path(dir,'Collett (2023) code','models', model_nam))

# loading
set_cmdstan_path( '/home/josema/.cmdstan/cmdstan-2.37.0')
mod = cmdstan_model( file.path( dir,'Collett (2023) code','models', model_nam), ) # to create the C++ model

fit = mod$sample( 
  data = dlist, 
  iter_warmup = 1000, iter_sampling = 1000, 
  chains = 4, parallel_chains = 4) 
# adapt_delta=0.95, max_treedepth=20 ) #, init=0
fit$save_object( 
  file = file.path( dir,'Collett (2023) code','models', str_replace(model_nam,'.stan','.RDS')) )
# 
# Note:
#   In this model the priors make more difference than the adapt_delta, 
#   max_treedepth or a seed.


# load model
fit = readRDS( 
  file = file.path( dir,'Collett (2023) code','models', str_replace(model_nam,'.stan','.RDS')) )
# fit$summary()

# parameter recovery
post = as_draws_df( fit )
# post = subset_draws( post, variable='lambda' )

summarize_draws( subset_draws( post, variable=c('lambda','gamma') ) )
summarize_draws( exp(subset_draws( post, variable='beta' )) )
#
# Frequentist estimates
# [1] "cof == 1"
# mean_l    lower_l upper_l
# 0.0041276 0.0005066 0.0336326
# mean_g lower_g upper_g
# 0.47   0.34    0.67
#              mean_cof lower_cof upper_cof
# age          0.97      1.00      0.93
# gender1      0.60      1.40      0.26
# disease1     0.68      1.91      0.24
# disease2     1.15      2.86      0.46
# 
# [1] "cof == 2"
# mean_l   lower_l upper_l
# 0.011582 0.0027942 0.0480084
# mean_g lower_g upper_g
# 0.29   0.23    0.36
#              mean_cof lower_cof upper_cof
# age          0.98      1.00      0.95
# gender1      1.08      1.91      0.61
# disease1     1.79      3.49      0.92
# disease2     1.62      2.96      0.88
# 
# [1] "cof == 3"
# mean_l  lower_l upper_l
# 1.2e-06 2e-07   7.3e-06       0
# mean_g lower_g upper_g
# 1.45   1.13    1.87
#              mean_cof lower_cof upper_cof
# age          0.98      1.01      0.95
# gender1      0.89      1.81      0.44
# disease1     0.34      1.03      0.11
# disease2     1.73      3.45      0.86
# 
# [1] "cof == 4"
# mean_l    lower_l   upper_l
# 0.0008196 0.0002055 0.0032689
# mean_g lower_g upper_g
# 0.42   0.35    0.49
#              mean_cof lower_cof upper_cof
# age          1.03      1.05      1.00
# gender1      0.92      1.44      0.59
# disease1     0.70      1.21      0.41
# disease2     1.38      2.16      0.88
# 
# Note:
#   These are the "same" parameter estimates ("minus" the randomness)


# summaries
Fij_MC = summarize_draws( subset_draws( post, variable='Fij_MC' ) ) 
Fij_I = summarize_draws( subset_draws( post, variable='Fij_I' ) ) 

par(mfrow=c(2,2))
for( j in 1:4 ){
  idx = str_detect( Fij_MC$variable, paste0(j,'[:punct:]$') )
  plot(t, Fij_I$mean[idx], type='l', lwd=4, lty=1, ylim=c(0,0.15))
  lines(t, Fij_MC$mean[idx], col='red', lty=2)
  legend( 'topleft', bty='n', col=c('white','black','red'), lty=c(1,1,2),
          legend=c( paste0('P(C=',j,') = ', round( tail(Fij_I$mean[idx], 1), 2) ),
                    'Double Exponential Quadrature', 
                    'Monte Carlo approx') )
}
par(mfrow=c(1,1))

