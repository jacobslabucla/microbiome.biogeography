#' Generating heatmaps to match the metabolic maps containing "GMM" - gut metabolic modules
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param path_to_significant_results_tsv filepath to Maaslin2/significant_results.tsv
#' @param ystring choose between "hierachy_annotations", "feature_annotations", "metabolic_map", or "Map"
#' @param titlestring title of plot, e.g. "Mucosal"
#' @param colorvector character vector of two colors
#' @return a ggplot2 object encoding a heatmap
#' @export 

install.packages(c(""))

generate_interregional_taxa_dotplot_TYPE <- function(path_to_significant_results_tsv, dataset,titlestring, colorvector){
  ##### test function inputs 
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
  #luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv")
  #luminal <- read.delim("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  
  #cols <- viridis::viridis(2)
  #####
  
  luminal<-read.table(path_to_significant_results_tsv, header=TRUE)
  luminal <- luminal %>% filter(metadata=="Type" & qval<0.05)
  cols=c(colorvector)
  

  if(dataset=="ucla_original"){
    annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/genus_taxonomy.csv", header=TRUE)
    luminal$feature<-gsub("X", "", luminal$feature)

  }
  
  
  if(dataset=="ucla_validation"){
    annotation <- read.csv("ImmDef-Mouse-Biogeography-Analysis/genus_Mucosal_taxonomy.csv", header=TRUE)
    annotation$feature<-annotation$X.OTU.ID
    annotation$feature<-gsub("/",".",annotation$feature)
    annotation$feature<-gsub("-",".",annotation$feature)
  }
  
  if(dataset=="cs"){
    annotation <- read.csv("CS-Facility-Analysis/Type_L6/genus_taxonomy.csv", header=TRUE)
    annotation$feature<-gsub("; ",".",annotation$feature)
    annotation$feature<-gsub(".s__.*","",annotation$feature)
    annotation$feature<-gsub("-",".",annotation$feature)
    annotation$feature<-gsub("/",".",annotation$feature)
  }
  
  if(dataset=="spf_gavage"){
    annotation <- read.csv("Humanized-Biogeography-Analysis/Genus_Taxonomy.csv", header=TRUE)
    annotation <- select(annotation, -c("X.OTU.ID", "QIIME_seqs"))
    annotation$feature<-annotation$taxonomy
    annotation$feature<-gsub("; ",".",annotation$feature)
    annotation$feature<-gsub("\\.s__.*","",annotation$feature)
    annotation$feature<-gsub("-",".",annotation$feature)
    annotation$feature<-gsub("/",".",annotation$feature)
  }
  
  if(dataset=="hum_gavage"){
    annotation <- read.csv("Humanized-Biogeography-Analysis/Genus_Taxonomy.csv", header=TRUE)
    annotation <- select(annotation, -c("X.OTU.ID", "QIIME_seqs"))
    annotation$feature<-annotation$taxonomy
    annotation$feature<-gsub("; ",".",annotation$feature)
    annotation$feature<-gsub("\\.s__.*","",annotation$feature)
    annotation$feature<-gsub("-",".",annotation$feature)
    annotation$feature<-gsub("/",".",annotation$feature)
  }
  
  data<- (merge(luminal, annotation, by = 'feature'))
  testing <- luminal$feature
  testing2<-(data$feature)
  print("Dropped features: (should say 0 )")
  print(c(setdiff(testing, testing2), setdiff(testing2, testing)))
  

  data <- data %>% select(c("feature","coef", "qval", "Phylum", "Genus"))
  data <- unique(data)
  data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
  data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")
  
  # Make plot 
  res_plot <- data %>% select(c("coef", "qval","Genus"))
  res_plot <- unique(res_plot)
  res_plot <- res_plot %>%
    mutate(site = ifelse(coef< 0, "Luminal", "Mucosal"))
  
  y = tapply(res_plot$coef, res_plot$Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  res_plot$Genus= factor(as.character(res_plot$Genus), levels = names(y))
  
  g1<- res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.05, abs(coef) > 0) %>%
    ggplot2::ggplot(aes(coef, Genus, fill = site)) +
    geom_bar(stat = "identity") +
    cowplot::theme_cowplot(12) +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_fill_manual(values = cols) +
    labs(x = "Effect size (Mucosal/Luminal)",
         y = "",
         fill = "") +
    theme(legend.position = "none")+
    ggtitle({{titlestring}}) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

return(g1)
}
