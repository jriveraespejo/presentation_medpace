
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

