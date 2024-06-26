#' Generate violin plots 
#'
#' 
#' Can be used for any continuous variable (Y) against X (Site_General, Site, Type)
#' 
#' 
#' @author Julianne C. Yang
#' @author Jonathan P. Jacobs
#' @param input_data dataframe containing "SampleID", and alpha diversity metrics as columns. "SampleID" must be typed as such.
#' @param input_metadata dataframe containing "SampleID" and all relevant metadata. Can contain more SampleIDs than are in input_data.
#' @param X what you want to plot on X axis. Available options are Site, Type, or Site_General.
#' @param Y name of any alpha diversity metric that you want to plot on y-axis.
#' @param fillvariable name of any metadata column that you want to use in the fill argument of aes()
#' @param min y-axis minimum
#' @param max y-axis maximum
#' @return a ggplot2 encoding a violin plot
#' @export 
#' @examples
#' 
#'generate_adiv_plots(input_dataframe, 
#' metadata_frame,
#' Site, observed_otus, Site,0,500) +
#'ggtitle ("Mucosal") +
#'theme(plot.title = element_text(hjust = 0.5)) +
#'stat_compare_means(comparisons = list(c("DC", "PC"),
#'                                      c("DC", "Cec"),
#'                                       c("DC", "Ile"),
#'                                        c("DC", "Jej"),
#'                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
#' 


generate_adiv_plots <- function(input_data, input_metadata, X, Y, fillvariable, min, max){
  #read in files
  metadata<- as.data.frame(input_metadata)
  data<-as.data.frame(input_data)
  metadata$SampleID <-gsub("-",".",metadata$SampleID)
  
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #Shorten site names 
  data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
  data$Type= plyr::revalue(data$Type, c("Luminal" = "L", "Mucosal" = "M"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "C", "Ileum"="I", "Jejunum" = "J", "Duodenum"= "D"))
  
  #Ensure correct ordering of levels 
  data$Site <- factor(data$Site, levels = c("D","J", "I","C","PC","DC"))
 
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{fillvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

generate_adiv_plots_shotgun <- function(chao1_filepath,
                                        otus_filepath,
                                        pielou_filepath,
                                        shannon_filepath,
                                        metadata_filepath, X, Y, fillvariable, min, max){
  #chao1<- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv")
  #otus <- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv")
  #pielou<- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv")
  #shannon <- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv")
  #metadata <- read.delim("Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv")
  
  chao1<- read.delim(here(chao1_filepath))
  otus<- read.delim(here(otus_filepath))
  pielou<- read.delim(here(pielou_filepath))
  shannon<- read.delim(here(shannon_filepath))
  
  names(chao1)
  temp1<- merge(chao1,otus, by="X")
  temp2<- merge(pielou,shannon, by="X")
  data<-merge(temp1,temp2,by="X")
  data$sampleid <- data$`X`
  data$sampleid <- gsub("-",".",data$sampleid)
  
  metadata <-read.delim(here(metadata_filepath)) #mapping file
  site_metadata <- metadata %>% select("Site", "sampleid")
  intermediate<- (merge(data, site_metadata, by = 'sampleid'))
  data<- intermediate

  print(summary(data$observed_features))
  print(summary(data$pielou_evenness))
  
  #Shorten site names 
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Jejunum" = "J"))
  
  data$Site <- factor(data$Site, levels = c("J","DC"))

  plot <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{fillvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
    return(plot)
  
}

merge_adiv_shotgun <- function(chao1_filepath,
                                        otus_filepath,
                                        pielou_filepath,
                                        shannon_filepath,
                                        metadata_filepath){
  #chao1<- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv")
  #otus <- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv")
  #pielou<- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv")
  #shannon <- read.delim("Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv")
  #metadata <- read.delim("Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv")
  
  chao1<- read.delim(here(chao1_filepath))
  otus<- read.delim(here(otus_filepath))
  pielou<- read.delim(here(pielou_filepath))
  shannon<- read.delim(here(shannon_filepath))
  
  names(chao1)
  temp1<- merge(chao1,otus, by="X")
  temp2<- merge(pielou,shannon, by="X")
  data<-merge(temp1,temp2,by="X")
  data$sampleid <- data$`X`
  data$sampleid <- gsub("-",".",data$sampleid)
  
  metadata <-read.delim(here(metadata_filepath)) #mapping file
  site_metadata <- metadata
  intermediate<- (merge(data, site_metadata, by = 'sampleid'))
  data<- intermediate
  
  return(data)
  
}

generate_adiv_plots_donors <- function(chao1_filepath,
                                        otus_filepath,
                                        pielou_filepath,
                                        shannon_filepath,
                                        metadata_filepath, X, Y, fillvariable, min, max){
  chao1<- readr::read_delim(here(chao1_filepath),delim="\t")
  otus<- readr::read_delim(here(otus_filepath),delim="\t")
  pielou<- readr::read_delim(here(pielou_filepath),delim="\t")
  shannon<- readr::read_delim(here(shannon_filepath),delim="\t")
  
  names(chao1)
  temp1<- merge(chao1,otus, by="...1")
  temp2<- merge(pielou,shannon, by="...1")
  data<-merge(temp1,temp2,by="...1")
  data$SampleID <- data$`...1`
  
  metadata <-readr::read_delim(here(metadata_filepath),delim = "\t") #mapping file
  site_metadata <- metadata %>% select("Site", "SampleID")
  intermediate<- (merge(data, site_metadata, by = 'SampleID'))
  data<- intermediate
  
  print(summary(data$observed_features))
  print(summary(data$pielou_evenness))
  
  #Shorten site names 
  data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
  data$Type= plyr::revalue(data$Type, c("Luminal" = "Lum", "Mucosal" = "Muc"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  #Ensure correct ordering of levels 
  data$Site <- factor(data$Site, levels = c("Duo","Jej", "Ile","Cec","PC","DC"))
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{fillvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

