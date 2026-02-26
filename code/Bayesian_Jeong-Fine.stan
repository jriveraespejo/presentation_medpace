
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

