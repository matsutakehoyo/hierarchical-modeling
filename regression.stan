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
  real<lower=0> sigma; 
  real a0; 
	real<lower=0> a_line_s;
	vector[n_line] a_line; 
	real<lower=0> a_mouse_s;
	vector[n_mouse] a_mouse;
}

transformed parameters {
  real<lower=0> nu;
  vector[n_obs] mu; 
  for (i in 1:n_obs){
    mu[i] = a0 + a_line[obs2line[i]] + a_mouse[obs2mouse[i]];
  }
  nu = nu_minus_one+1;
}
model {
	sigma ~ cauchy(0,1);
  a0 ~ cauchy(0,1);
	a_line ~ normal(0, a_line_s);
	a_line_s ~ cauchy(0,1);
  a_mouse ~ normal(0, a_mouse_s);
  a_mouse_s ~ cauchy(0,1);
  nu_minus_one ~ exponential(1/29.0);
  y ~ student_t(nu, mu, sigma);

}

generated quantities {
  vector[n_line] diff;
  vector[n_line] mu_line;
  vector[n_mouse] mu_mouse;
  diff[1] = a_line[1] - a_line[2];
  diff[2] = a_line[1] - a_line[3];
  diff[3] = a_line[2] - a_line[3];
  
  mu_line = a0 + a_line; 
  for (m in 1:n_mouse){
    mu_mouse[m] = a0 + a_line[mouse2line[m]] + a_mouse[m];
  }
}
