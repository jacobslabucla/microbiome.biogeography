#' Generating heatmaps to match the metabolic maps containing "GMM" - gut metabolic modules
#'
#'
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param path_to_all_results_tsv filepath to Maaslin2/all_results.tsv
#' @param path_to_module_key filepath to the GMM key with annotations for each GMM
#' @param ystring choose between "hierachy_annotations", "feature_annotations", "metabolic_map", or "Map"
#' @param titlestring title of plot, e.g. "Mucosal"
#' @param colorvector character vector of two colors
#' @return a ggplot2 object encoding a heatmap
#' @export


generate_interregional_GBM_barplot <- function(path_to_all_results_tsv, path_to_Module_Key, titlestring, colorvector){
  ##### test function inputs
  #annotation <- read.csv("GBM_Module_Key.csv", header=TRUE)
  #luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")

  #cols <- viridis::viridis(2)
  #####

  luminal<-read.table(path_to_all_results_tsv, header=TRUE)
  luminal <- luminal %>% filter(metadata=="Site_General" & qval<0.05)

  annotation <- read.csv(path_to_Module_Key, header=TRUE)

  data<- (merge(luminal, annotation, by = 'feature'))
  cols=c(colorvector)


  res_plot <- data %>% select(c("coef", "qval","annotation"))
  res_plot <- unique(res_plot)
  res_plot <- res_plot %>%
      mutate(site = ifelse(coef< 0, "Colon", "SI"))

  y = tapply(res_plot$coef, res_plot$annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  res_plot$annotation= factor(as.character(res_plot$annotation), levels = names(y))

  g1<- res_plot %>%
      arrange(coef) %>%
      filter(qval < 0.05, abs(coef) > 0) %>%
      ggplot2::ggplot(aes(coef, annotation, fill = site)) +
      geom_bar(stat = "identity") +
      cowplot::theme_cowplot(12) +
      theme(axis.text.y = element_text(face = "italic")) +
      scale_fill_manual(values = cols) +
      labs(x = "Effect size (SI/Colon)",
           y = "",
           fill = "") +
      theme(legend.position = "none")+
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))



  return(g1)
}
