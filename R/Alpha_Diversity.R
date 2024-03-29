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
  data$Type= plyr::revalue(data$Type, c("Luminal" = "Lum", "Mucosal" = "Muc"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  #Ensure correct ordering of levels 
  data$Site <- factor(data$Site, levels = c("Jej", "DC"))
 
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{fillvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

generate_adiv_plots_shotgun <- function(chao1_filepath,
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
  data$sampleid <- data$`...1`
  
  metadata <-readr::read_delim(here(metadata_filepath),delim = "\t") #mapping file
  site_metadata <- metadata %>% select("Site", "sampleid")
  intermediate<- (merge(data, site_metadata, by = 'sampleid'))
  data<- intermediate

  print(summary(data$observed_features))
  print(summary(data$pielou_evenness))
  
  #Shorten site names 
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Jejunum" = "Jej"))
  
  data$Site <- factor(data$Site, levels = c("Jej","DC"))

  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{fillvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

