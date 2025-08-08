#' Coverage Check Function
#' 
#' Performs a coverage check analysis by generating data from a model, fitting it multiple times
#' with different fit arguments, and plotting the results.
#' 
#' @param initialData An initial dataset to fit to determine 'true' parameters for further generation. 
#' @param fitting_model A ctModel object used for fitting the data
#' @param niter Number of iterations to run
#' @param fit_args Named list of fit argument sets to test 
#' (e.g., list(boot = list(optimcontrol = list(bootstrapUncertainty = TRUE)), 
#' hess = list(optimcontrol = list(bootstrapUncertainty = FALSE))))
#' @param cores Number of cores to use for parallel processing
#' @param plot_every Print plots every n iterations (default = 10)
#' 
#' @return A list containing the results data.table and final plots

ctModelCoverage_check <- function(initialData, fitting_model, niter, fit_args, 
  cores = 10, plot_every = max(c(10,cores))) {

  
  neededPacks = c('future','future.apply')
  sapply(neededPacks,function(x) {
    if(!requireNamespace(x, quietly = TRUE)) stop(paste("Package", x, "is required but not installed."))
  })
  
  # Set up parallel processing
  future::plan(future::multisession, workers = cores)
  
  # Check if parallel setup worked
  message(paste("Using", future::nbrOfWorkers(), "workers for parallel processing"))
  
  
  # Default fit arguments
  default_fit_args <- list(
    cores = 1
  )
  
  # Ensure fit_args is a named list
  if(is.null(names(fit_args)) || any(names(fit_args) == "")) {
    stop("fit_args must be a named list (e.g., list(boot = list(...), hess = list(...)))")
  }
  
  
  # Fit the model to get true parameters (use first fit_args configuration)
  initial_fit_args <- default_fit_args
  initial_fit <- do.call(ctStanFit, c(list(datalong = initialData, ctstanmodel = fitting_model), initial_fit_args))
  truepars <- initial_fit$stanfit$rawest
  
  # CRITICAL STEP: Generate new data samples from the fitted model
  # This ensures proper parameter specification for all iterations
  message("Generating samples for all iterations using ctStanGenerateFromFit...")
  generated_samples <- ctStanGenerateFromFit(fit = initial_fit, nsamples = niter, cores = cores)

  # Define the function to run for each iteration
  run_iteration <- function(iter_idx) {
    tryCatch({
      # Extract generated data for this iteration from ctStanGenerateFromFit
      y <- array(generated_samples$generated$Y[iter_idx,,], 
        dim = dim(generated_samples$generated$Y[1,,,drop=FALSE])[-1])
      colnames(y) <- fitting_model$manifestNames
      
      # Create data for this iteration using the original data structure
      # but with the newly generated Y values
      dat <- initialData
      dat[, fitting_model$manifestNames] <- y
      
      # Fit with each set of fit arguments
      iter_results <- list()
      
      for(fit_type in names(fit_args)) {
        # Merge default args with specific fit args
        current_fit_args <- modifyList(default_fit_args, fit_args[[fit_type]])
        
        # Fit the model with current arguments
        current_fit <- do.call(ctStanFit, c(list(datalong = dat, ctstanmodel = fitting_model), current_fit_args))
        estimates <- t(apply(current_fit$stanfit$rawposterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))
        
        # Compile results for this fit type
        current_result <- data.frame(
          iteration = iter_idx,
          type = fit_type,
          truepars = truepars,
          estimates,
          par=unlist(ctFitgetparnamesfromraw(current_fit)))
        
        # current_result$par <- rownames(current_result)
        
        iter_results[[fit_type]] <- current_result
      }
      
      # Combine all results for this iteration
      return(do.call(rbind, iter_results))
      
    }, error = function(e) {
      message(paste("Error in iteration", iter_idx, ":", e$message))
      return(NULL)
    })
  }
  
  # Initialize results storage
  all_results <- list()
  
  # Run iterations in batches and plot every plot_every iterations
  for(batch_start in seq(1, niter, by = plot_every)) {
    batch_end <- min(batch_start + plot_every - 1, niter)
    batch_indices <- batch_start:batch_end
    
    message(paste("Running iterations", batch_start, "to", batch_end))
    
    # Run batch in parallel
    batch_results <- future.apply::future_lapply(batch_indices, run_iteration, future.seed = TRUE)
    
    # Remove NULL results (failed iterations)
    batch_results <- batch_results[!sapply(batch_results, is.null)]
    
    # Add to all results
    all_results <- c(all_results, batch_results)
    
    # Create plots if we have any results
    if (length(all_results) > 0) {
      # Combine all results so far
      res <- do.call(rbind, all_results)
      
      # Convert to data.table and calculate coverage statistics
      res2 <- data.table(res)
      res2[, mean50 := mean(X50.), by = 'par']
      res2[, coverage := mean(X2.5. <= truepars & X97.5. >= truepars), by = c('type', 'par')]
      
      # Create plots
      g1 <- ggplot(res2, aes(x = factor(par), y = mean50 - truepars, colour = type)) +
        geom_boxplot(outliers = FALSE) + 
        geom_boxplot(aes(y = X2.5. - truepars), outliers = FALSE) +
        geom_boxplot(aes(y = X97.5. - truepars), outliers = FALSE) +
        theme_bw() +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom') +
        labs(title = paste("Parameter Estimation Bias (Iterations 1 -", batch_end, ")"), 
             y = "Estimate - True Value")
      
      g2 <- ggplot(res2[iteration == 1,], aes(x = par, y = coverage, colour = type)) +
        geom_point(alpha = 0.7) +
        geom_segment(aes(xend = par, yend = 0.95), linetype = 'dotted') +
        theme_bw() +
        geom_hline(yintercept = 0.95, linetype = 'dashed') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom') +
        labs(title = paste("Coverage Probability (Iterations 1 -", batch_end, ")"), 
             y = "Coverage")
      
      # Print plots
      print(gridExtra::grid.arrange(g1, g2, nrow = 2))
    }
  }
  
  # Final results compilation
  if (length(all_results) > 0) {
    res <- do.call(rbind, all_results)
  } else {
    stop("All iterations failed")
  }
  
  # Convert to data.table and calculate final coverage statistics
  res2 <- data.table(res)
  res2[, mean50 := mean(X50.), by = 'par']
  res2[, coverage := mean(X2.5. <= truepars & X97.5. >= truepars), by = c('type', 'par')]
  
  # Create final plots
  g1 <- ggplot(res2, aes(x = factor(par), y = mean50 - truepars, colour = type)) +
    geom_boxplot(outliers = FALSE) + 
    geom_boxplot(aes(y = X2.5. - truepars), outliers = FALSE) +
    geom_boxplot(aes(y = X97.5. - truepars), outliers = FALSE) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom') +
    labs(title = "Parameter Estimation Bias", y = "Estimate - True Value")
  
  g2 <- ggplot(res2[iteration == 1,], aes(x = par, y = coverage, colour = type)) +
    geom_point(alpha = 0.7) +
    geom_segment(aes(xend = par, yend = 0.95),linetype = 'dotted') +
    theme_bw() +
    geom_hline(yintercept = 0.95, linetype = 'dashed') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom') +
    labs(title = "Coverage Probability", y = "Coverage")
  
  # Clean up parallel processing
  future::plan(future::sequential)
  
  # Return results
  return(list(
    results = res2,
    bias_plot = g1,
    coverage_plot = g2,
    successful_iterations = length(all_results),
    total_iterations = niter
  ))
}
