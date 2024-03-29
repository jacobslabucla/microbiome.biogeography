% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Alpha_Diversity.R
\name{generate_adiv_plots}
\alias{generate_adiv_plots}
\title{Generate violin plots}
\usage{
generate_adiv_plots(input_data, input_metadata, X, Y, fillvariable, min, max)
}
\arguments{
\item{input_data}{dataframe containing "SampleID", and alpha diversity metrics as columns. "SampleID" must be typed as such.}

\item{input_metadata}{dataframe containing "SampleID" and all relevant metadata. Can contain more SampleIDs than are in input_data.}

\item{X}{what you want to plot on X axis. Available options are Site, Type, or Site_General.}

\item{Y}{name of any alpha diversity metric that you want to plot on y-axis.}

\item{fillvariable}{name of any metadata column that you want to use in the fill argument of aes()}

\item{min}{y-axis minimum}

\item{max}{y-axis maximum}
}
\value{
a ggplot2 encoding a violin plot
}
\description{
Can be used for any continuous variable (Y) against X (Site_General, Site, Type)
}
\examples{

generate_adiv_plots(input_dataframe, 
metadata_frame,
Site, observed_otus, Site,0,500) +
ggtitle ("Mucosal") +
theme(plot.title = element_text(hjust = 0.5)) +
stat_compare_means(comparisons = list(c("DC", "PC"),
                                     c("DC", "Cec"),
                                      c("DC", "Ile"),
                                       c("DC", "Jej"),
                                       c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

}
\author{
Julianne C. Yang

Jonathan P. Jacobs
}
