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



generate_interregional_taxa_barplot_shotgun <- function(path_to_significant_results_tsv, titlestring, colorvector){
  ##### test function inputs 
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
  #luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv")
  #luminal <- read.delim("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv")
  #cols <- viridis::viridis(2)
  #####
  luminal<-readr::read_delim(here(path_to_significant_results_tsv),delim="\t")
  significant_taxa <- luminal %>% filter(metadata=="Site" & qval<0.05)
  df <- significant_taxa$feature
  df<-as.data.frame(df)
  df$feature <- df[,1]
  
  df$Phylum <- gsub(".*\\.p__", "", df$feature)
  df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
  df$Order <- gsub(".*\\.o__", "", df$feature)
  df$Order <- gsub("\\.f__.*", "", df$Order)
  df$Order <- paste0(df$Order, " (o)")
  df$Family <- gsub(".*\\.f__", "", df$feature)
  df$Family <- gsub("\\.g__.*", "", df$Family)
  df$Family<- paste0(df$Family, " (f)")
  df$Genus <- gsub(".*\\.g__", "", df$feature)
  df$Genus <- gsub("\\.g__.*", "", df$Genus)
  df$Genus<- paste0(df$Genus, " (g)")
  df$Species <- gsub(".*\\.s__", "", df$feature)
 
  df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
  df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))
  
  df <- df %>%
    mutate(annotation = ifelse(!startsWith(Species, "GGB"), Species, Family_Species))
  
  luminal<-readr::read_delim(here(path_to_significant_results_tsv),delim="\t")
  luminal <- luminal %>% filter(metadata=="Site" & qval<0.05)
  cols=c(colorvector)
  

  data<- (merge(luminal, df, by = 'feature'))
  testing <- luminal$feature
  testing2<-(data$feature)
  print("Dropped features: (should say 0 )")
  print(c(setdiff(testing, testing2), setdiff(testing2, testing)))
  

  data <- data %>% select(c("feature","coef", "qval", "Phylum", "Genus","annotation"))
  data <- unique(data)
  data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
  data$Phylum_Genus<-paste(data$Phylum,data$annotation,sep=" : ")
  
  # Make plot 
  res_plot <- data %>% select(c("coef", "qval","Genus","Phylum_Genus","annotation","Phylum"))
  res_plot <- unique(res_plot)
  res_plot <- res_plot %>%
    filter(!grepl("Mitochondria", annotation))
  res_plot <- res_plot %>%
    mutate(site = ifelse(coef< 0, "Distal_Colon", "Jejunum"))
  
  y = tapply(res_plot$coef, res_plot$annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  res_plot$annotation= factor(as.character(res_plot$annotation), levels = names(y))
  
  
  ggplot_data <- res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.05, abs(coef) > 0) 
    
  
  g1<- res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.05, abs(coef) > 0) %>%
    ggplot2::ggplot(aes(coef, annotation, fill = site)) +
    geom_bar(stat = "identity") +
    cowplot::theme_cowplot(12) +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_fill_manual(values = cols) +
    #geom_text(aes(label = annotation, color = Phylum), hjust = -0.2) +  # Add text labels colored by Phylum
    #scale_fill_gradientn(colors = cols,breaks = bk) +
    labs(x = "Effect size (Jejunum/Distal_Colon)",
         y = "",
         fill = "") +
    theme(legend.position = "none")+
    ggtitle({{titlestring}}) +
    theme(plot.title = element_text(hjust = 0.5))
  
  newList <- list(dataframe=ggplot_data,plot=g1)
  return(newList)

}

generate_interregional_taxa_barplot_SITE <- function(path_to_significant_results_tsv, titlestring, colorvector){
  ##### test function inputs 
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
  #luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv")
  #luminal <- read.delim("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv")
  #cols <- viridis::viridis(2)
  #####
  luminal<-readr::read_delim(here(path_to_significant_results_tsv),delim="\t")
  significant_taxa <- luminal %>% filter(metadata=="Site_General" & qval<0.05)
  df <- significant_taxa$feature
  df<-as.data.frame(df)
  df$feature <- df[,1]
  
  df$Phylum <- gsub(".*\\.p__", "", df$feature)
  df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
  df$Order <- gsub(".*\\.o__", "", df$feature)
  df$Order <- gsub("\\.f__.*", "", df$Order)
  df$Order <- paste0(df$Order, " (o)")
  df$Family <- gsub(".*\\.f__", "", df$feature)
  df$Family <- gsub("\\.g__.*", "", df$Family)
  df$Family<- paste0(df$Family, " (f)")
  df$Genus <- gsub(".*\\.g__", "", df$feature)
  df$Genus <- gsub("\\.g__.*", "", df$Genus)
  df$Species <- gsub(".*\\.s__", "", df$feature)
  
  df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
  df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))
  
  df <- df %>%
    mutate(level1 = ifelse(nchar(Genus) != 0, Genus, Family))
  df <- df %>%
    mutate(annotation = ifelse(level1!= " (f)", level1, Order))
  
  luminal<-readr::read_delim(here(path_to_significant_results_tsv),delim="\t")
  luminal <- luminal %>% filter(metadata=="Site_General" & qval<0.05)
  cols=c(colorvector)
  
  
  data<- (merge(luminal, df, by = 'feature'))
  testing <- luminal$feature
  testing2<-(data$feature)
  print("Dropped features: (should say 0 )")
  print(c(setdiff(testing, testing2), setdiff(testing2, testing)))
  
  
  data <- data %>% select(c("feature","coef", "qval", "Phylum", "Genus","annotation"))
  data <- unique(data)
  data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
  data$Phylum_Genus<-paste(data$Phylum,data$annotation,sep=" : ")
  
  # Make plot 
  res_plot <- data %>% select(c("coef", "qval","Genus","Phylum_Genus","annotation","Phylum"))
  res_plot <- unique(res_plot)
  res_plot <- res_plot %>%
    filter(!grepl("Mitochondria", annotation))
  res_plot <- res_plot %>%
    mutate(site = ifelse(coef< 0, "Colon", "SI"))
  
  y = tapply(res_plot$coef, res_plot$annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  res_plot$annotation= factor(as.character(res_plot$annotation), levels = names(y))
  
  
  ggplot_data <- res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.05, abs(coef) > 0) 
  
  
  g1 <- res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.05, abs(coef) > 0) %>%
    ggplot2::ggplot(aes(coef, annotation, fill = site)) +
    geom_bar(stat = "identity") +
    cowplot::theme_cowplot(12) +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_fill_manual(values = cols) +
    #geom_text(aes(label = annotation, color = Phylum), hjust = -0.2) +  # Add text labels colored by Phylum
    #scale_fill_gradientn(colors = cols,breaks = bk) +
    labs(x = "Effect size (SI/Colon)",
         y = "",
         fill = "") +
    theme(legend.position = "none")+
    ggtitle({{titlestring}}) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(g1)
  
}
generate_interregional_taxa_barplot_SITE <- function(path_to_significant_results_tsv, dataset,titlestring, colorvector){
  ##### test function inputs 
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
  luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv")
  luminal <- read.delim("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv")
  luminal <- read.delim("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  #luminal <- read.delim("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")
  
  #cols <- viridis::viridis(2)
  #####
  
  luminal<-read.table(path_to_significant_results_tsv, header=TRUE)
  luminal <- luminal %>% filter(metadata=="Site_General" & qval<0.05)
  cols=c(colorvector)
  
  
  if(dataset=="ucla_original_lum"){
    annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/Genus_Luminal_taxonomy.csv", header=TRUE)
    #luminal$feature<-gsub("X", "", luminal$feature)
    
  }
  
  if(dataset=="ucla_original_muc"){
    annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/Genus_Mucosal_taxonomy.csv", header=TRUE)
    #luminal$feature<-gsub("X", "", luminal$feature)
    
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
  
  
  data <- data %>% select(c("feature","coef", "qval", "annotation"))
  data <- unique(data)
  #data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
  #data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")
  
  # Make plot 
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

generate_interregional_taxa_barplot_TYPE <- function(path_to_significant_results_tsv, dataset,titlestring, colorvector){
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
  
  return(g1)
}

