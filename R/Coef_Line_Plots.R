#' Data wrangling for the coefficient Line Plots
#' @author Julianne C. Yang
#' @author Jonathan P. Jacobs
#' @param file_paths filepaths to dataframes containing statistical testing results
#' @param cohort_prefixes strings describing what's in each dataframe
#' @param feature_value string descriptor of feature name
#' @param new_coef integer value for the effect size of the reference variable
#' @param new_value string descriptor of the reference variable
#' @return a dataframe aggregating all results dataframes with cohort prefixes appended
#' @export 

process_results_files <- function(file_paths, feature_value, 
                                  new_value="Distal_Colon", new_coef=0, cohort_prefixes,filter_by="Site") {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == {{filter_by}} & feature == feature_value)
    
    # Add a new row to the data frame
    new_row <- list(metadata = "Site",
                    feature = feature_value,
                    value = new_value,
                    coef = new_coef,
                    stderr = 0,
                    N = NA,
                    N.not.0 = NA,
                    pval = NA,
                    qval = 100)
    data <- rbind(data, new_row)
    
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