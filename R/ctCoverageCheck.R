#' Coverage Check Function
#' 
#' Performs a coverage check analysis by generating data from a model, fitting it multiple times
#' with different fit arguments, and plotting the results. Works with \code{fitArgs}/
#' \code{fittingModel} combinations that fit with either \code{backend='stan'} (the
#' \code{\link{ctFit}} default) or \code{backend='julia'}.
#'
#' @param initialData An initial dataset to fit to determine 'true' parameters for further generation. 
#' @param fittingModel A ctModel object used for fitting the data
#' @param niter Number of iterations to run
#' @param fitArgs Named list of fit argument sets to test 
#' (e.g., list(boot = list(optimcontrol = list(uncertainty = 'bootstrap')), 
#' hess = list(optimcontrol = list(uncertainty = 'hessian'))))
#' @param cores Number of outer simulation iterations to run in parallel.
#' Inner model fits use \code{fitCores} by default to avoid nested
#' parallelism.
#' @param fitCores Number of cores for each inner \code{\link{ctFit}} call.
#' This is enforced after merging each \code{fitArgs} entry so that outer
#' simulation \code{cores} cannot accidentally be reused by worker fits.
#' @param generateCores Number of cores for \code{\link{ctGenerateFromFit}}
#' when generating replicated datasets. Defaults to \code{fitCores}; set it
#' explicitly to use more cores for generation.
#' @param plotEvery Print plots every n iterations (default = 10)
#' 
#' @return A list containing the results data.table and final plots
#' @export

ctCoverageCheck <- function(initialData, fittingModel, niter, fitArgs, 
  cores = 10, fitCores = 1, generateCores = fitCores,
  plotEvery = max(c(10,cores))) {
  
  ctCoverageFitArgs <- function(default_fit_args, fitArgs, fitCores){
    out <- utils::modifyList(default_fit_args, fitArgs)
    out$cores <- fitCores
    out
  }

  mean50 <- X50. <- coverage <- X2.5. <- X97.5. <- type <- iteration <- NULL   # To avoid R CMD check note about undefined global variable
  
  neededPacks = c('future','future.apply')
  sapply(neededPacks,function(x) {
    if(!requireNamespace(x, quietly = TRUE)) stop(paste("Package", x, "is required but not installed."))
  })
  
  cores <- suppressWarnings(as.integer(cores[1]))
  fitCores <- suppressWarnings(as.integer(fitCores[1]))
  generateCores <- suppressWarnings(as.integer(generateCores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  if(!is.finite(fitCores) || is.na(fitCores) || fitCores < 1) fitCores <- 1L
  if(!is.finite(generateCores) || is.na(generateCores) || generateCores < 1) generateCores <- 1L
  cores <- min(cores, niter)
  
  
  # Default fit arguments. Inner fits are single-core by default to avoid
  # nested parallelism; fitCores is enforced after merging fitArgs below.
  default_fit_args <- list()
  
  # Ensure fitArgs is a named list
  if(is.null(names(fitArgs)) || any(names(fitArgs) == "")) {
    stop("fitArgs must be a named list (e.g., list(boot = list(...), hess = list(...)))")
  }
  fit_args_cores <- vapply(fitArgs, function(x) !is.null(x$cores),
    logical(1))
  if(any(fit_args_cores)) {
    warning('Ignoring top-level cores entries in fitArgs; use fitCores to control inner ctFit cores.',
      call.=FALSE)
  }
  if(cores > 1 && fitCores > 1) {
    warning('ctCoverageCheck is using outer parallelism and inner ctFit cores > 1. ',
      'This can oversubscribe CPUs; prefer cores > 1 with fitCores = 1 unless you have budgeted for nested workers.',
      call.=FALSE)
  }
  
  
  # Fit the model to get true parameters (use first fitArgs configuration)
  initial_fit_args <- ctCoverageFitArgs(default_fit_args, list(), cores)
  initial_fit <- do.call(ctFit, c(list(datalong = initialData, model= fittingModel), initial_fit_args))
  # .ctFitRawEstimate()/.ctFitRawPosterior()/.ctFitRawParNames() (R/ctBackendSummary.R)
  # read `$stanfit$...` or `$estimate$...` depending on which backend produced
  # the fit, so this works whether fittingModel/fitArgs select backend='stan'
  # (the default) or backend='julia'.
  truepars <- .ctFitRawEstimate(initial_fit)
  
  # CRITICAL STEP: Generate new data samples from the fitted model
  # This ensures proper parameter specification for all iterations
  message("Generating samples for all iterations using ctGenerateFromFit...")
  generated_samples <- ctGenerateFromFit(fit = initial_fit, nsamples = niter,
    cores = generateCores)
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add=TRUE)
  
  # Set up parallel processing over generated datasets.
  future::plan(future::multisession, workers = cores)
  
  # Check if parallel setup worked
  message(paste("Using", future::nbrOfWorkers(), "workers for coverage iterations"))

  # Define the function to run for each iteration
  run_iteration <- function(iter_idx) {
    tryCatch({
      # Extract generated data for this iteration from ctGenerateFromFit
      y <- array(generated_samples$generated$Y[iter_idx,,], 
        dim = dim(generated_samples$generated$Y[1,,,drop=FALSE])[-1])
      colnames(y) <- fittingModel$manifestNames
      
      # Create data for this iteration using the original data structure
      # but with the newly generated Y values
      dat <- initialData
      dat[, fittingModel$manifestNames] <- y
      
      # Fit with each set of fit arguments
      iter_results <- list()
      
      for(fit_type in names(fitArgs)) {
        # Merge default args with specific fit args
        current_fit_args <- ctCoverageFitArgs(default_fit_args,
          fitArgs[[fit_type]], fitCores)
        
        # Fit the model with current arguments
        current_fit <- do.call(ctFit, c(list(datalong = dat, model= fittingModel), current_fit_args))
        estimates <- t(apply(.ctFitRawPosterior(current_fit), 2, quantile, probs = c(0.025, 0.5, 0.975)))

        # Compile results for this fit type
        current_result <- data.frame(
          iteration = iter_idx,
          type = fit_type,
          truepars = truepars,
          estimates,
          par=.ctFitRawParNames(current_fit))
        
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
  
  # Run iterations in batches and plot every plotEvery iterations
  for(batch_start in seq(1, niter, by = plotEvery)) {
    batch_end <- min(batch_start + plotEvery - 1, niter)
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
  
  # Return results
  return(list(
    results = res2,
    bias_plot = g1,
    coverage_plot = g2,
    successful_iterations = length(all_results),
    total_iterations = niter
  ))
}
