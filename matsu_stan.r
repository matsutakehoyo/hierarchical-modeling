library(ggmcmc)
library(HDInterval)


tidy_posterior <- function(fit, parameter=NA){
	if(is.na(parameter)){
		print("tidy everything")
		posterior <- ggs(fit) 

		nChains <- attributes(posterior)$nChains
		nParameters <- attributes(posterior)$nParameters
		nIterations <- attributes(posterior)$nIterations
		nBurnin <- attributes(posterior)$nBurnin
		nThin <- attributes(posterior)$nThin

		posterior <- posterior %>%
			separate(Parameter, into=c("Parameter", "ind"), sep="\\[", fill="right") %>%
			replace_na(list(ind=0)) # %>% select(Parameter, ind) %>% distinct %>% print(n=30)

		posterior <- posterior %>%
			separate(ind, into=c("ind", "ind2"), sep=c(","), fill="right") %>% 
			replace_na(list(ind2=0)) %>% 
			mutate(ind = gsub("([0-9]+).*$", "\\1",  ind) %>% as.integer(),
				   ind2 = gsub("([0-9]+).*$", "\\1", ind2) %>% as.integer()) %>%
			mutate(Parameter = as.factor(Parameter)) %>%
			mutate(ind1 = ind) %>% 
			mutate(ind = group_indices(., Parameter, ind1, ind2)) %>%
			group_by(Parameter) %>%
			mutate(ind = ind - min(ind) + 1) %>%
			ungroup()

	} else {
		#use loop because gss family arguments only takes one element 
		posterior <- list()
		for (p in seq_along(parameter)){
			print(paste("tidy", parameter[p]))
			posterior[[p]] <- ggs(fit, family = parameter[p])	
			#posterior[[p]] <- ggs(fit, family = paste0(parameter[p], "\\[(\\d+)\\]"))				
			#preserve this attributes for later calculations
			
			nChains <- attributes(posterior[[p]])$nChains
			nParameters <- attributes(posterior[[p]])$nParameters
			nIterations <- attributes(posterior[[p]])$nIterations
			nBurnin <- attributes(posterior[[p]])$nBurnin
			nThin <- attributes(posterior[[p]])$nThin
			
			posterior[[p]] <- posterior[[p]] %>%
				separate(Parameter, into=c("Parameter", "ind"), sep="\\[", fill="right") %>%
				replace_na(list(ind=0))


			posterior[[p]] <- posterior[[p]] %>%
				separate(ind, into=c("ind", "ind2"), sep=c(","), fill="right") %>% 
				replace_na(list(ind2=0)) %>% 
				mutate(ind = gsub("([0-9]+).*$", "\\1",  ind) %>% as.integer(),
					   ind2 = gsub("([0-9]+).*$", "\\1", ind2) %>% as.integer()) %>%
				mutate(Parameter = as.factor(Parameter)) %>%
				mutate(ind1 = ind) %>% 
				mutate(ind = group_indices(., Parameter, ind1, ind2)) %>%
				group_by(Parameter) %>%
				mutate(ind = ind - min(ind) + 1) %>%
				ungroup()


		}
		posterior <- bind_rows(posterior) %>%
			mutate(Parameter = fct_relevel(Parameter, parameter)) #reorder Parameters
	
	}
	attr(posterior, "nChains") <- nChains
	attr(posterior, "nParameters") <- nParameters
	attr(posterior, "nIterations") <- nIterations
	attr(posterior, "nBurnin") <- nParameters
	attr(posterior, "nBurnin") <- nBurnin
	attr(posterior, "nThin") <- nThin
	return(posterior)
}

sum_posterior <- function(posterior, hdi1=0.89, hdi2=0.50){
	posterior_sum <- posterior %>%
		group_by(Parameter, ind) %>%			
		summarise(data = list(data.frame(hist(value, breaks = 'Scott', plot = FALSE)[c('mids', 'counts')])),
				  hdi_lower = hdi(value, credMass=hdi1)[1], 
				  hdi_upper = hdi(value, credMass=hdi1)[2], 
				  hdi_lower2 = hdi(value, credMass = hdi2)[1], 
				  hdi_upper2 = hdi(value, credMass = hdi2)[2], 
				  mode = dens_mode(value)) %>% 
		ungroup() %>% 
		mutate(data = map(data, mutate, width = mids[2] - mids[1])) %>%
		mutate(data = map(data, mutate, max_count = max(counts))) %>%
		unnest() %>%
		arrange(Parameter, ind)	
	attr(posterior_sum, "hdi1") <- hdi1
	attr(posterior_sum, "hdi2") <- hdi2
	return(posterior_sum)

}

#calculate a reasonable range for posterior distributions excluding ontliers
posterior_range <- function(posterior, hdi=0.99){
	posterior_range <- posterior %>%
		group_by(Parameter, ind) %>%
		summarise(hdi_lower = hdi(value, credMass=hdi)[1], 
	          hdi_upper = hdi(value, credMass=hdi)[2]) %>%
		ungroup() %>%
		summarise(min = min(hdi_lower),
	    	      max = max(hdi_upper)) %>%
		unlist()
	return(posterior_range)
}

posterior_hist <-function(posterior, ...){
	summary <- sum_posterior(posterior, ...)
	range <- posterior_range(posterior)
	ggplot(summary, aes(fill = factor(ind))) + 
		geom_blank(aes(y = max_count*1.3)) +
		geom_col(aes(x = mids, y = counts, width = width), position="identity") +
		geom_errorbarh(aes(x=mode, y = max_count * 1.15, 
						   xmin=hdi_lower, xmax=hdi_upper, 
						   height = max_count /30 ),
					   stat="unique" , alpha=1/2) +
		geom_point(aes(x=mode, y=max_count * 1.15), stat = "unique", alpha=1/2) +
		geom_text(aes(x = hdi_lower, y = max_count * 1.1, label = round(hdi_lower,2)), 
				  stat="unique", alpha=1/2, vjust="top", hjust=0.7) +
		geom_text(aes(x = hdi_upper, y = max_count * 1.1, label = round(hdi_upper,2)), 
				  stat="unique", alpha=1/2, vjust="top", hjust=0.3) +
		geom_text(aes(x = mode, y = max_count *1.1, label = round(mode,2)), 
				  stat="unique", vjust="top") +
		geom_text(aes(x = mode, y = max_count *1.2), 
					label = paste0(attributes(summary)$hdi1 * 100, "% HDI"), 
					stat="unique", vjust="bottom") +
		facet_grid(ind~Parameter, scales = "free") +
		labs(x="value", y="counts") +
		xlim(range[1], range[2]) +
		theme(legend.position = "none") +
		basic_style
}

#additional hdi1 and hdi2 specify credibility mass [0, 1]
posterior_cata <-function(posterior, ...){
	summary <- sum_posterior(posterior, ...)
	range <- posterior_range(posterior)
	ggplot(summary, aes(x=as.factor(ind))) +
		geom_linerange(aes(ymin=hdi_lower, ymax=hdi_upper), color="dodgerblue", size=1) +
  		geom_linerange(aes(ymin=hdi_lower2, ymax=hdi_upper2), color="dodgerblue3", size=2) +
  		geom_point(aes(y=mode), color="dodgerblue4", size=3) +
  		coord_flip() +
  		facet_grid(Parameter~., space="free", scale="free") +
  		labs(x='Parameter', y='Value') +
		ylim(range[1], range[2]) +
  		basic_style
}

# the order of x labels is upside down…
posterior_violin <- function(posterior, ...){
	summary <- sum_posterior(posterior, ...)
	range <- posterior_range(posterior)
	ggplot(posterior) +
		geom_violin(aes(x=as.factor(ind), y=value), fill="gray", color=NA) +
		geom_linerange(data=summary, aes(x=as.factor(ind), ymin=hdi_lower, ymax=hdi_upper), color="gray60", size=1) +
  		geom_linerange(data=summary, aes(x=as.factor(ind), ymin=hdi_lower2, ymax=hdi_upper2), color="gray30", size=2) +
  		geom_point(data=summary, aes(x=as.factor(ind), y=mode), color="black", size=3) +
  		coord_flip() +
  		facet_grid(Parameter~., space="free", scale="free") +
  		labs(x='Parameter', y='Value') +
		ylim(range[1], range[2]) +
  		basic_style
}

#summary info n_eff, rhat 
mcmc_summary <- function(stan, family=NA){
	if (is.na(family) %>% isTRUE()) {
		summary <- summary(stan)$summary
		summary_df <- tibble(Parameter = summary %>% rownames) %>% bind_cols(as.tibble(summary))
		summary_df <- summary_df %>%
			dplyr::filter(Parameter!="lp__") %>%
			separate(Parameter, into=c("Parameter", "ind"), sep="\\[", fill="right") %>%
			replace_na(list(ind=0))
			
		summary_df <- summary_df %>%
			separate(ind, into=c("ind", "ind2"), sep=c(","), fill="right") %>% 
			replace_na(list(ind2=0)) %>% 
			mutate(ind = gsub("([0-9]+).*$", "\\1",  ind) %>% as.integer(),
				   ind2 = gsub("([0-9]+).*$", "\\1", ind2) %>% as.integer())
		
	} else {
		summary_df <- list()
		for (f in seq_along(family)){
			summary <- summary(stan, pars=family[f])$summary#[, c("n_eff", "Rhat")]	
			summary_df[[f]] <- tibble(Parameter = summary %>% rownames) %>% bind_cols(as.tibble(summary))
			summary_df[[f]] <- summary_df[[f]] %>%
				dplyr::filter(Parameter!="lp__") %>%
				separate(Parameter, into=c("Parameter", "ind"), sep="\\[", fill="right") %>%
				replace_na(list(ind=0)) 
				
			summary_df[[f]] <- summary_df[[f]] %>%
				separate(ind, into=c("ind", "ind2"), sep=c(","), fill="right") %>% 
				replace_na(list(ind2=0)) %>% 
				mutate(ind = gsub("([0-9]+).*$", "\\1",  ind) %>% as.integer(),
					   ind2 = gsub("([0-9]+).*$", "\\1", ind2) %>% as.integer())
		}
		summary_df <- bind_rows(summary_df)	
		summary_df <- summary_df %>% mutate(Parameter = fct_relevel(Parameter, family)) #reorder Parameters
	}
	summary_df <- summary_df %>% 
		mutate(ind1 = ind, 
			   ind = group_indices(., Parameter, ind, ind2)) %>%
		group_by(Parameter) %>%
		mutate(ind = (ind - min(ind) + 1) %>% as.integer()) %>% ungroup()

	lp <- suppressWarnings(get_logposterior(stan, inc_warmup = FALSE))
    SS <- prod(length(lp), length(lp[[1L]]))
	attr(summary_df, "SS") <- SS
	summary_df
}


mcmc_ess <- function(posterior_summary){		
	ggplot(posterior_summary, aes(x=as.factor(ind))) +
	    geom_linerange(aes(ymin = 0, ymax = n_eff/attributes(posterior_summary)$SS)) +
    	geom_point(aes(y = n_eff/attributes(posterior_summary)$SS)) +
	    facet_grid(Parameter~., scale="free", space="free") +
	    ylim(0, 1) +
	    coord_flip() +
	    labs(y="Effective sample size / Sample size", x='') +
	    basic_style
	
}



#modified ggs_rhat
# mcmc_rhat <- function(D, scaling=1.5, greek=FALSE){
# 	# The computations follow BDA, pg 296-297, and the notation tries to be
# 	# consistent with it
# 	# Compute between-sequence variance using psi.. and psi.j	
# 	psi.dot <- D %>%
# 		dplyr::group_by(Parameter, ind, Chain) %>%
# 		dplyr::summarize(psi.dot=mean(value))
# 	psi.j <- D %>%
# 		dplyr::group_by(Parameter, ind) %>%
# 		dplyr::summarize(psi.j=mean(value))
# 	b.df <- dplyr::inner_join(psi.dot, psi.j, by=c("Parameter", "ind"))
# 	glimpse(b.df)
# 	B <- b.df %>%
# 		dplyr::group_by(Parameter, ind) %>%
# 		dplyr::summarize(B=var(psi.j-psi.dot)*attributes(D)$nIterations)
# 	B <- unique(B)
# 	# Compute within-sequence variance using s2j
# 	s2j <- D %>%
# 		dplyr::group_by(Parameter, ind, Chain) %>%
# 		dplyr::summarize(s2j=var(value))
# 	W <- s2j %>%
# 		dplyr::group_by(Parameter, ind) %>%
# 		dplyr::summarize(W=mean(s2j))
# 	# Merge BW and compute the weighted average (wa, var.hat+) and the Rhat
# 	BW <- dplyr::inner_join(B, W, by=c("Parameter", "ind")) %>%
# 		dplyr::mutate(
# 			wa= (((attributes(D)$nIterations-1)/attributes(D)$nIterations )* W) +
# 	    		((1/ attributes(D)$nIterations)*B),
# 			Rhat=sqrt(wa/W))
# 	# For parameters that do not vary, Rhat is Nan. Move it to NA
# 	BW$Rhat[is.nan(BW$Rhat)] <- NA
# 	# Plot
# 	f <- ggplot(BW, aes(x=Rhat, y=as.factor(ind))) + 
# 		geom_point() +
# 		facet_grid(Parameter~., space="free", scales="free") +
# 		xlab(expression(hat("R"))) + 
# 		ggtitle("Potential Scale Reduction Factors") +
# 		basic_style
# 	if (greek) {
# 		f <- f + scale_y_discrete(labels = parse(text = as.character(BW$Parameter)))
# 	}
# 	# If scaling, add the scale
# 	if (!is.na(scaling)) {
# 		# Use the maximum of Rhat if it is larger than the prespecified value
# 		scaling <- ifelse(scaling > max(BW$Rhat, na.rm=TRUE), scaling, max(BW$Rhat, na.rm=TRUE))
# 		f <- f + xlim(min(BW$Rhat, na.rm=TRUE), scaling)
# 	}
#  f
# }

mcmc_rhat <- function(posterior_summary){
	ggplot(posterior_summary) +
	    geom_point(aes(x = as.factor(ind), y = Rhat)) +
	    facet_grid(Parameter~., scale="free", space="free") +
	    ylab(expression(hat("R"))) + xlab("") +
		#ggtitle("Potential Scale Reduction Factors") +
	    coord_flip() +
	    basic_style
}


mcmc_chain <-function(posterior){
	ggplot(posterior, aes(x=Iteration, y=value, color=as.factor(Chain))) +
  		geom_line() +
  		facet_grid(Parameter+ind~Chain, scale="free", labeller = labeller(Chain = label_both)) +
  		labs(x='Iteration', y='Value', color="chain") +
  		scale_fill_brewer(palette="Set3") +
  		basic_style +
  		theme(legend.position = "none")
}

mcmc_density <- function(posterior){
	ggplot(posterior, aes(x=value, fill=as.factor(Chain), color=as.factor(Chain))) +
		geom_density(alpha=1/6)+
		facet_grid(Parameter+ind~., scale="free", labeller = labeller(Chain = label_both)) +
		labs(x='Value', y='Density') +
		scale_fill_brewer(palette="Set3") +
		basic_style +
		theme(legend.position = "none")
}


mcmc_running <- function(posterior){
	chain_mean <- posterior %>% 
		group_by(Parameter, ind, Chain) %>%
		summarise(mean=mean(value))

	posterior <- posterior %>% dplyr::arrange(Parameter, ind, Iteration) %>%
	    dplyr::group_by(Parameter, ind, Chain) %>%
	    dplyr::mutate(running_mean=cumsum(value) / Iteration)

	ggplot() +
		geom_hline(aes(yintercept=mean), chain_mean, colour="gray50") +
		geom_line(aes(x=Iteration, y=running_mean, color=as.factor(Chain)), posterior) +
		facet_grid(Parameter+ind~., scale="free") +
		labs(x='Iteration', y='Running Mean', color='Chain') +
		basic_style +
		theme(legend.position = "none")
}

mcmc_autocorrelation <- function(D, nLags=50, greek=FALSE) {
  nIter <- attr(D, 'nIteration')
  if (nIter < nLags) {
    warning(sprintf('nLags=%d is larger than number of iterations, computing until max possible lag %d', nLags, nIter))
    nLags <- nIter
  }

  # No way to make the following use summarize(), as of dplyr 0.3
  # https://github.com/hadley/dplyr/issues/154
  # Temporary workaround using dplyr 0.2 and do()
  wc.ac <- D %>%
    dplyr::group_by(Parameter, ind, Chain) %>%
    dplyr::do(ac(.$value, nLags))

  # Manage multiple chains
  if (attributes(D)$nChains <= 1) {
    f <- ggplot(wc.ac, aes(x=Lag, y=Autocorrelation)) + 
      geom_bar(stat="identity", position="identity") #+ ylim(-1, 1)
    if (!greek) {
      f <- f + facet_wrap(~ Parameter+ind)
    } else {
      f <- f + facet_wrap(~ Parameter+ind, labeller = label_parsed)
    }
  } else {
    f <- ggplot(wc.ac, aes(x=Lag, y=Autocorrelation, colour=as.factor(Chain), fill=as.factor(Chain))) + 
      geom_bar(stat="identity", position="identity") + #ylim(-1, 1) +
      scale_fill_discrete(name="Chain") + scale_colour_discrete(name="Chain")
    if (!greek ) {
      f <- f + facet_grid(Parameter+ind ~ Chain, labeller = labeller(Chain = label_both))
    } else {
      f <- f + facet_grid(Parameter+ind ~ Chain, labeller = labeller(Chain = label_both, Parameter = label_parsed))
    }
  }

 f + basic_style + theme(legend.position = "none")
}


model_diagnose <- function(posterior, posterior_summary, name=NA){
	# parameters <- posterior %>% select(Parameter) %>% pull() %>% unique()
	# only process Parameters with fewer than 50 indexes
	parameters <- posterior %>% 
		group_by(Parameter) %>% 
		summarise(ind1= max(ind1), 
		          ind2= max(ind2)) %>%
		filter(ind1 < 50 & ind2 < 50) %>%
		select(Parameter) %>% pull()

	print(paste("Parameters:", parameters))
	
	for (p in parameters){
		print(p)
		n_ind <- posterior %>% filter(Parameter==p) %>% select(ind) %>% distinct %>% nrow()	
		n_plots <- 6			
		lay <- rbind(c(1,3,3,3,4),
		             c(1,3,3,3,4),
		             c(2,5,5,5,6),
		             c(2,5,5,5,6))
		plots <- list()
		if (n_ind==1){

			d_rhat <- mcmc_rhat(posterior_summary %>% 
								filter(Parameter==p))
			d_ess <- mcmc_ess(posterior_summary %>% 
								filter(Parameter==p))
			d_chain <- mcmc_chain(posterior %>% 
								filter(Parameter==p))
			d_dens <- mcmc_density(posterior %>% 
								filter(Parameter==p))
			d_ac <- mcmc_autocorrelation(posterior %>% 
								filter(Parameter==p))
			d_rmean <- mcmc_running(posterior %>% 
								filter(Parameter==p))
			plots[[1]] <- list(d_rhat, d_ess, d_chain, d_dens, d_ac, d_rmean)
		} else {
			
			for (i in 1:ceiling(n_ind/n_plots)){
				
				d_rhat <- mcmc_rhat(posterior_summary %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))	
				d_ess <- mcmc_ess(posterior_summary %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))	
				d_chain <- mcmc_chain(posterior %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))
				d_dens <- mcmc_density(posterior %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))
				d_ac <- mcmc_autocorrelation(posterior %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))
				d_rmean <- mcmc_running(posterior %>% 
									filter(Parameter==p) %>%
									filter(ind >= (i-1)*n_plots+1) %>%
									filter(ind <= i*n_plots))
				plots[[i]] <- list(d_rhat, d_ess, d_chain, d_dens, d_ac, d_rmean)
			}
		}

		mplots <- lapply(plots, function(g) gridExtra::marrangeGrob(grobs = g, layout_matrix = lay))
		if(is.na(name)){
			path = paste0("diagnostics")
		} else {
			path = paste0(gsub(".stan", "", name), " diagnostics")
		}
		dir.create(path, recursive=TRUE, showWarnings=FALSE) #save marker data here	
		pdf(paste0(path, "/", p, ".pdf"), width = 20, height = 10)
			print(mplots)
		dev.off()
	}	
}

model_posterior <- function(posterior, name=NA, ...){
	if(is.na(name)){
		path = paste0("posterior")
	} else {
		path = paste0(gsub(".stan", "", name), " posterior")
	}
	
	dir.create(path, recursive=TRUE, showWarnings=FALSE) 
	# parameters <- posterior$Parameter %>% unique
	# only process Parameters with fewer than 50 indexes
	parameters <- posterior %>% 
		group_by(Parameter) %>% 
		summarise(ind1= max(ind1), 
		          ind2= max(ind2)) %>%
		filter(ind1 < 50 & ind2 < 50) %>%
		select(Parameter) %>% pull()
	for (p in parameters){
		print(p)
		n_ind <- posterior %>% filter(Parameter==p) %>% select(ind) %>% distinct %>% nrow()
		if (n_ind<5){
			p_posterior <- posterior_hist(posterior%>%filter(Parameter==p), ...)
		} else if(n_ind < 10 )  {
			p_posterior <- posterior_violin(posterior%>%filter(Parameter==p), ...)
		} else {
			p_posterior <- posterior_cata(posterior%>%filter(Parameter==p), ...)
		}
		ggsave(paste0(path, "/", p, ".pdf"), width = 10, height = 10)
	}
}


#pairs plot to check parameter correlation
mcmc_pairs <- function(model, family=NA, name=NA){
	require(GGally)

	if (is.na(family)){
		df <- ggs(model) %>% 
		spread(Parameter, value) %>%
		select(-Iteration, -Chain)
	} else {
		par_list <- ggs(model) %>% distinct(Parameter) %>% pull()
		df <- ggs(model) %>% 
		filter(Parameter %in% grep(par_list, pattern=family, value=TRUE, perl=TRUE)) %>%
		spread(Parameter, value) %>%
		select(-Iteration, -Chain)
	}
	
	
	ggp <- df %>%
		ggpairs(., diag = 'blank', lower = 'blank') 
		
	print("customizing pairwise comparison plots: diagonal")
	#ggpair diagonal 
	for(i in seq_along(df)) {
		x <- df[,i] %>% unlist() #unlist to untibble it
		p <- ggplot(data.frame(x), aes(x)) + 
			geom_histogram(binwidth=hist_bw(x)[1], fill='grey20') +
			theme(text=element_text(size=14), 
			      axis.text.x=element_text(angle=40, vjust=1, hjust=1)) 
		ggp <- putPlot(ggp, p, i, i) 
	}
	
	#upper panels change the cor display
	zcolat <- seq(-1, 1, length=81)
	zcolre <- c(zcolat[1:40]+1, rev(zcolat[41:81]))
	print("customizing pairwise comparison plots: upper")
	for(i in 1:(ncol(df)-1)) {
		for(j in (i+1):ncol(df)) {
			x <- df[,i] %>% pull() %>% as.numeric()  #because cor requires numeric
			y <- df[,j] %>% pull() #%>% as.numeric()
			if (class(y) != 'factor' ){
				y <- as.numeric(y)
				r <- cor(x, y, method='spearman', use='pairwise.complete.obs')
				zcol <- lattice::level.colors(r, at=zcolat,
											  col.regions=colorRampPalette(c(scales::muted('red'), 'white', scales::muted('blue')), 
											  							 space='rgb')(81))
				textcol <- ifelse(abs(r) < 0.4, 'grey20', 'gray80')
				ell <- ellipse::ellipse(r, level=0.95, type='l', npoints=50, scale=c(.2, .2), centre=c(.5, .5))
				p <- ggplot(data.frame(ell), aes(x=x, y=y))
				p <- p + theme_bw() + 
					theme(
						plot.background=element_blank(),
						panel.grid.major=element_blank(), 
						panel.grid.minor=element_blank(),
						panel.border=element_blank(), axis.ticks=element_blank()
					)
				p <- p + geom_polygon(fill=zcol, color=zcol)
				p <- p + geom_text(data=NULL, x=.5, y=.5, label=100*round(r, 2), size=6, col=textcol)
				ggp <- putPlot(ggp, p, i, j)	
			} else {
				ggp[i, j] <- ggp[i, j] + scale_fill_manual(values = cbPalette)				
			}
		}
	}
	
	#lower panels
	print("customizing pairwise comparison plots: lower")
	for(j in 1:(ncol(df)-1)) {
		for(i in (j+1):ncol(df)) {
			x <- df[,j] %>% pull() 
			y <- df[,i] %>% pull() 
			p <- ggplot(data.frame(x, y), aes(x=x, y=y)) +
				geom_hex() +
				# scale_fill_gradient(low = "gray70", high = "black") +
				# geom_point(alpha=1/100) +
				theme(#text=element_text(size=14), 
				      axis.text.x=element_text(angle=40, vjust=1, hjust=1))
			ggp <- putPlot(ggp, p, i, j)
		}
	}
	ggp

	#save plot
	if(is.na(name)){
		path = paste0("diagnostics")
	} else {
		path = paste0(gsub(".stan", "", name), " diagnostics")
	}
	dir.create(path, recursive=TRUE, showWarnings=FALSE) 
	if (is.na(family)){
		ggsave(paste0(path, "/pairs.pdf"), ggp, width = 10, height = 10)	
	} else {
		ggsave(paste0(path, "/pairs_", family, ".pdf"), ggp, width = 10, height = 10)	
	}
	
}

diff_zero <- function(sample){
	pos <- sum(sample>0)/length(sample)
	if(pos>0.5)
		return(pos)
	else 
		return(1-pos)
}


# calculate the maximum radius of ROPE that would cover all the posteriors (grouped by ind)
max_rope_radius <- function(posterior, comp_val=0.0){
	df <- posterior %>% 
		group_by(Parameter, ind) %>%			
		summarise(max_radius = hdi(value, credMass=0.97) %>% as.numeric() %>% abs() %>% max() * 1.1 - comp_val)
	max(df$max_radius)
}

# function for computing area in ROPEs ()
rope_area <- function(values, max_radius, comp_val=0.0) {
  rope_radius = seq(0, max_radius, length=201 ) 
  area = rep( NA , length(rope_radius) )
  for (r in 1:length(rope_radius)) {
    area[r] = ( sum( values > (comp_val-rope_radius[r]) 
                            & values < (comp_val+rope_radius[r]) ) 
                         / length(values) )
  }
  tibble(radius=rope_radius, area=area)
}

# plot ROPE that would incclude fraction of the posterior distribution
plot_rope <- function(posterior, comp_val=0.0, fraction=0.89){
	max_radius <- max_rope_radius(posterior, comp_val)
	df <- posterior %>% 
		group_by(Parameter, ind) %>%			
		summarise(rope = list(rope_area(value, max_radius, comp_val))) %>% 
		ungroup() %>%
		# radius at which fraction=0.89 of the posterior is included inside
		# mutate(r = map(rope, ~filter(., abs(area - fraction) == min(abs(area - fraction))) %>% .$radius[[1]])) %>%
		mutate(r = map(rope, ~filter(., abs(area - fraction) == min(abs(area - fraction))) %>% pull(radius))) %>%
		unnest(r, .drop=FALSE) %>%
		unnest(rope, .drop=FALSE) 

	ggplot(df) +
		facet_grid(ind~Parameter) +
		geom_line(aes(x=radius, y=area, color=factor(ind)), size=1) +
		geom_segment(aes(x=0, xend=r, y=fraction, yend=fraction), 
		             linetype="dotted", stat="unique", size=0.5) +
		geom_segment(aes(x=r, xend=r, y=0.05, yend=fraction), 
		             linetype="dotted", stat="unique", size=0.5) +
		geom_text(aes(x=r/2, y=fraction, label=paste0(fraction*100, "%")), 
		          stat="unique", vjust=-0.5, hjust=0.5) +
		geom_text(aes(label=paste(format(r, digits=3)), x=r, y=0), 
		          stat="unique", vjust=0.5) +
		geom_point(aes(x=r, y=fraction), stat="unique", color="red") +
		labs(y="Posterior in ROPE", x=paste("ROPE radius around ", comp_val)) +
		scale_x_continuous(expand = c(0, 0)) +
		theme(legend.position = "none") +
		basic_style
}
