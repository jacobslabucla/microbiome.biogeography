#' Generate Upset Plots for GBM and GMM 
#'
#' 
#' 
#' 
#' 
#' @author Julianne C. Yang
#' @author Jonathan P. Jacobs
#' @param file_paths filepaths to dataframes containing statistical testing results
#' @param cohort_prefixes strings describing what's in each dataframe
#' @return a dataframe aggregating all results dataframes with cohort prefixes appended
#' @export 
#' 
#' 
process_results_for_upset_plot <- function(file_paths, cohort_prefixes, filter_by="Site") {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here::here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- dplyr::filter(results, metadata == {{filter_by}} & qval<0.05)
  
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    data <- data %>% mutate(Cohort = cohort)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  return(data_all)
}

