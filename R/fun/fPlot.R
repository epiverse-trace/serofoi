
fCombinedPlots <- function(res, dat, lambda_sim = NA, max_lambda = NA, size_text = 25) {
  
  if (is.character (res$fit) == FALSE)  {
    if  (class(res$fit@sim$samples)  != "NULL" ) {
      
      summary_mod <- extract_summary_mod (res, dat)
      foi <- rstan::extract(res$fit, 'foi', inc_warmup = FALSE)[[1]]
      
      
      prev_expanded <- get_prev_expanded(foi, dat)
      plot_prev <- 
        ggplot(prev_expanded) +
        geom_ribbon(aes(x= age, ymin = predicted_prev_lower, ymax = predicted_prev_upper), fill = '#c994c7') +
        geom_line(aes(x= age, y = predicted_prev), colour = '#7a0177') +
        geom_errorbar(aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u), width = 0.1) +
        geom_point(aes(age, p_obs_bin, size = bin_size), fill = '#7a0177', colour = 'black', shape = 21) +
        theme_bw(size_text) +
        coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
        theme(legend.position = 'none') +
        ylab ('Sero-positivity') + xlab("Age") 
      
      
      
      #-------- This bit is to get the actual length of the foi data  
      foi_dat <- res$foi_cent_est
      if (!is.na(lambda_sim)) {
        lambda_mod_length <- NROW(foi_dat) 
        lambda_sim_length <- length(lambda_sim)
        
        if (lambda_mod_length < lambda_sim_length) {
          remove_x_values <- lambda_sim_length - lambda_mod_length
          lambda_sim <- lambda_sim[-c(1:remove_x_values)]
        }
        
        foi_dat$simulated <- lambda_sim
        
      }
      
      #--------  
      foi_dat$medianv[1] <- NA
      foi_dat$lower[1] <- NA
      foi_dat$upper[1] <- NA
      
      
      plot_foi <- 
        ggplot(foi_dat) +
        geom_ribbon(aes( x= year, ymin = lower, ymax = upper), fill = '#41b6c4', alpha = 0.5) +
        geom_line(aes(x = year, y = medianv), colour = '#253494', size = size_text/8) +
        theme_bw(size_text) +
        coord_cartesian(ylim = c(0, max_lambda)) +
        ylab ('Force-of-Infection') + xlab ("Year")
      
      # browser()
      
      if (!is.na(lambda_sim)) {
        lambda_plot <- plot_foi + geom_line(aes(x = year, y = simulated), colour = 'red', size = 1.5) 
      }
      
      rhats <- get_table_rhats (res)
      
      plot_rhats <- ggplot(rhats, aes(year, rhat)) +
        geom_line(colour = 'purple') +
        geom_point() +
        coord_cartesian(ylim = c(0.7, 2)) +
        geom_hline(yintercept = 1.1, colour = 'blue', size = size_text/12) +
        theme_bw(size_text) +
        ylab ('Convergence (R^)') 
      
      
   
      
      # browser() # hasta aquí summary_mod es dataframe---> OK
      
      plot_data <- plot_info_table(t(summary_mod), size_text = size_text) 
      
      
      plot_arrange <- grid.arrange(plot_data, 
                                   plot_prev, 
                                   plot_foi, 
                                   plot_rhats, nrow = 4,
                                   heights = c(1.5, 1, 1, 1))
      # dev.off() #Tuve que colocar esto aquí porque por defecto me imprime el plot en la pantalla y no quiero eso
      res_plots <- list (plots = list(plot_data = plot_data,
                                      plot_prev =plot_prev,
                                      plot_foi  = plot_foi,
                                      plot_rhats = plot_rhats),
                         summary_mod     = summary_mod,
                         rhats           = rhats,
                         prev_expanded = prev_expanded)
      
      
      
    } } else
      
    {
      print ('model did not run')
      print_warning <- 'errors'
      df <- data.frame() 
      
      g0 <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x = 4, y = 5, label = print_warning) +
        theme_bw(25) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
        ylab('') + xlab('') 
      g1 <- g0
      g0 <- g0 +   labs(subtitle = res$model)+
        theme(plot.title = element_text(size=10))
      
      
      plots <- grid.arrange(g0, g1, g1, g1, g1, nrow = 5)
      res_p <- list (plots = plots,
                     loo_fit = c('Not available'),
                     summary_mod     = c('Not available'),
                     prev_expanded = c('Not available')
      )
      
      res_plots <-  list (plots = plots,
                          summary_mod  =   "model did not run")
      
      
    }
  
  
  
  return(res_plots)
  
  
}



#========================= PLOT comparison
get_model_comparison_plot <- function(res_comp) {
  model_comp  <- res_comp$model_comp
  best_model  <- as.character(res_comp$best_model_data1$model)
  best_modelP <- as.numeric(res_comp$best_model_data1$pvalue)
  
  emptyp <- ggplot(data = data.frame()) +
    geom_point() +
    xlim(0,1) + ylim (0,1) + theme_void() +
    annotate('text',  x = .5, y = .6, label = best_model, size = 14) +
    annotate('text',  x = .5, y = .55, label = '(best model)', size = 13) 
  
  
  infot <- filter(model_comp, converged == 'Yes') %>% dplyr::select(model, difference, diff_se, pvalue, best) %>%
    mutate(diff = round(difference, 2), 
           diff_se = round(diff_se, 2),
           pvalue = round(pvalue, 4))
  
  blank <- data.frame(x= 1:10, y = 1:100)
  table_pars <- 
    ggplot(blank, aes(x, y)) +
    geom_blank() + ylab('') + xlab ('') + 
    annotation_custom(tableGrob(d= infot,
                                theme = ttheme_default(base_size = 20))) +
    theme_void(30)
  
  pf <- plot_grid(emptyp, table_pars, nrow = 2 )
  
  return(pf)
  
}


vertical_plot_arrange_per_model <- function(PPC){
  
  pp <- grid.arrange(PPC$plots$plot_data, 
                     PPC$plots$plot_prev, 
                     PPC$plots$plot_foi, 
                     PPC$plots$plot_rhats, 
                     nrow = 4,
                     heights = c(1.5, 1, 1, 1))
  return(pp)
  
}


plot_info_table <- function(info, size_text){
  
  dato <- data.frame(
    y = NROW(info):1,
    text = paste0(rownames(info), ': ',info[,1])
  )
  p <- ggplot(dato, aes(x=1, y=y)) + 
    scale_y_continuous(limits=c(0, NROW(info) +1), breaks=NULL) +
    # scale_x_continuous(breaks=NULL) + 
    theme_void() +
    geom_text(aes(label=text), size = size_text/2.5, fontface = 'bold') 
  
  return(p)
}


