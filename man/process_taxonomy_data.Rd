% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Taxa_Barplots.R
\name{process_taxonomy_data}
\alias{process_taxonomy_data}
\title{Generate taxa summary plots by phyla
Plot does not have legend - legends were graphed separately}
\usage{
process_taxonomy_data(file_path)
}
\arguments{
\item{file_path}{filepath to the L6.csv file where you can wrangle names}
}
\value{
dataframe with wrangled names
}
\description{
Generate taxa summary plots by phyla
Plot does not have legend - legends were graphed separately
}
\examples{

generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS", 
"Luminal ( > 0.1\% Relative Abundance)",
 ".*g__",
  assign_cols,
   "Site")


}
\author{
Julianne Yang
}
