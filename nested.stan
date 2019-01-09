data {
  int<lower=1>  n_obs;
  int<lower=1>  n_mouse;
  int<lower=1>  n_line;
  int<lower=1, upper=n_line> obs2line[n_obs];
  int<lower=1, upper=n_mouse> obs2mouse[n_obs];
  int<lower=1, upper=n_line> mouse2line[n_mouse];
  real  y[n_obs];
}

parameters {
  real<lower=0> nu_minus_one; 
  real<lower=0> sigma_mouse; 
  vector[n_mouse] mu_mouse;
  real<lower=0> sigma_line; 
  vector[n_line] mu_line;
  real mu0; 
  real sigma0; 
}

transformed parameters {
real<lower=0> nu;
  nu = nu_minus_one+1;
}
model {

  for (i in 1:n_obs){
    y[i] ~ student_t(nu, mu_mouse[obs2mouse[i]], sigma_mouse);  
  }
  for (m in 1:n_mouse){
    mu_mouse[m] ~ normal(mu_line[mouse2line[m]], sigma_line);
  }
  mu_line ~ normal(mu0, sigma0);
  sigma0 ~ cauchy(0,1);
  sigma_line ~ cauchy(0,1);
  sigma_mouse ~ cauchy(0,1);
  nu_minus_one ~ exponential(1/29.0);
}

generated quantities {
  vector[n_line] diff;
  diff[1] = mu_line[1] - mu_line[2];
  diff[2] = mu_line[1] - mu_line[3];
  diff[3] = mu_line[2] - mu_line[3];
}
