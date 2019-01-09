library(tidyverse)
library(rstan)
library(tictoc)
source('take.R')
source('matsu_stan.R')
line <- c(-1, 0, 1)

df <- tibble(line_id=rep(1:3, each=10), # 3 lines with 10 mice each
             mouse_id=1:30) %>% 
	mutate(mouse=map_dbl(line[line_id], ~rnorm(1, ., 2))) %>% 
	mutate(obs = map(mouse, ~rnorm(30, ., 1))) %>% # 30 observations per mouse
	unnest()

glimpse(df)

# regression style modeling
model_data <- list(
   n_obs = nrow(df), 
   n_mouse = df %>% distinct(mouse_id) %>% nrow(), 
   n_line = df %>% distinct(line_id) %>% nrow(), 
   obs2line = df$line_id, 
   obs2mouse = df$mouse_id, 
   mouse2line = df %>% distinct(mouse_id, line_id) %>% pull(line_id),
   y = df$obs
)

model_name <- "regression.stan"
tic("sampling")
model_r <- rstan::stan(file=model_name,
              data = model_data,
              control = list(adapt_delta=0.9, max_treedepth=10),
              cores = getOption("mc.cores", 1L))
toc()
model_name <- gsub(".stan", "", model_name)

posterior_r <- tidy_posterior(model_r)
posterior_r %>% distinct(Parameter)
posterior_summary_r <- mcmc_summary(model)
model_diagnose(posterior_r, posterior_summary_r, name = model_name)
model_posterior(posterior_r, name = model_name)
pairs(model_r, pars=c("nu_minus_one", "sigma", "a0", "a_line_s", "a_line", "a_mouse_s", "a_mouse[1]"))

# nested model
model_name <- "nested.stan"
tic("sampling")
model_n <- rstan::stan(file=model_name,
              data = model_data,
              control = list(adapt_delta=0.9, max_treedepth=10),
              cores = getOption("mc.cores", 1L))
toc()
model_name <- gsub(".stan", "", model_name)

posterior_n <- tidy_posterior(model_n)
posterior_n %>% distinct(Parameter)
posterior_summary_n <- mcmc_summary(model_n)
model_diagnose(posterior_n, posterior_summary_n, name = model_name)
model_posterior(posterior_n, name = model_name)

# line estimates
line_posteriors <- ggplot() +
	geom_histogram(aes(x=value), fill="deeppink",
	               data=posterior_r %>% filter(Parameter=="mu_line") %>% mutate(model="regression")) +
	geom_histogram(aes(x=value), fill="deepskyblue", 
               data=posterior_n%>%filter(Parameter=="mu_line")%>%mutate(model="nested")) +
	geom_vline(aes(xintercept=value), data=tibble(value=line, ind=1:3)) +
	facet_grid(model~ind)

# mouse estimates 
mouse_posteriors <- ggplot() +
	geom_violin(aes(x=as.factor(ind), y=value), fill="deeppink", 
	            data = posterior_r %>% filter(Parameter=="mu_mouse") %>% mutate(model="regression")) +
	geom_violin(aes(x=as.factor(ind), y=value), fill="deepskyblue", 
            data = posterior_r %>% filter(Parameter=="mu_mouse") %>% mutate(model="nested")) +
	coord_flip() +
	geom_point(aes(x=as.factor(mouse_id), y=mouse), data=df %>% distinct(mouse_id, mouse)) +
	facet_grid(.~model)

# sigma
simga_mouse_posteriors <- ggplot() +
	geom_histogram(aes(x=value), fill="deeppink",
	               data=posterior_r %>% filter(Parameter=="sigma") %>% mutate(model="regression")) +
	geom_histogram(aes(x=value), fill="deepskyblue", 
               data=posterior_n%>%filter(Parameter=="sigma_mouse")%>%mutate(model="nested")) 
	facet_grid(model~ind)
