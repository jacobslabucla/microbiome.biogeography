#' Generate PCoA plots for beta diversity 
#'
#' 
#' (Site_General, Site, Type)
#' 
#' 
#' @author Julianne C. Yang
#' @author Jonathan P. Jacobs
#' @param input_data ordination file directly from QIIME2 principal_coordinates.py. read in with header=FALSE
#' @param input_metadata dataframe containing "SampleID" exactly as written
#' @param title string containing title of plot 
#' @param colorvariable string containing the variable you want to color the points by. Available options are "Site", "Type", or "Site_General".
#' @param colorvector named character vector with colors assigned to levels of colorvariable
#' @return a ggplot2 object encoding a pcoA plot
#' @export 


generate_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
  data<-as.data.frame(ordination_file)
  metadata <- as.data.frame(metadata)
  
  #store PC1 and Pc2
  PC1<-data[5,1]
  PC1 <-round(as.numeric(PC1)*100, digits=1)
  PC2<-data[5,2]
  PC2 <-round(as.numeric(PC2)*100, digits=1)
  PC1 <-as.character(PC1)
  str_PC1<-paste0("PC 1 (", PC1,"%)")
  str_PC2<-paste0("PC 2 (", PC2, "%)")
  
  #edit dataframe
  data<-data[,1:4]
  data <- dplyr::slice(data, 1:(n() - 4))     # Apply slice & n functions
  data<-data[-c(1,2,3,4,5,6,7,8,9),]
  
  #rename columns
  names(data)[names(data) == "V1"] <- "SampleID" 
  names(data)[names(data)=="V2"] <- "PC1" 
  names(data)[names(data)=="V3"] <- "PC2"
  names(data)[names(data)=="V4"] <- "PC3"
  # data$SampleID<-gsub(".","",data$SampleID)
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #declare factors
  data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))
  #data$Microbiota <-factor(data$Microbiota, levels=c("Humanized", "Cedars_SPF"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  colorvariable <- as.factor(colorvariable)
  if(colorvariable =="Site"){
    p<- ggplot2::ggplot(data, aes(x=PC1, y=PC2, colour=Site)) + 
      #geom_point(size=3) + 
      geom_point(aes(fill=Site), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      #scale_color_viridis_d()+
      xlab(str_PC1) +
      ylab(str_PC2) +
      cowplot::theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    #coord_fixed(ratio=1/2)+
    #labs(title= paste0({{title}}, " RPCA")) 
    p
  }
  else if (colorvariable =="Type"){
    data$Type = revalue(data$Type, c("Luminal"="Lum", "Mucosal"="Muc"))
    p<- ggplot2::ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
      #geom_point(size=3) + 
      geom_point(aes(fill=Type), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      #scale_color_viridis_d()+
      xlab(str_PC1) +
      ylab(str_PC2) +
      cowplot::theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    #coord_fixed(ratio=1/2)+
    #labs(title= paste0({{title}}, " RPCA")) 
    p
  }
  else if (colorvariable =="Site_General"){
    p<- ggplot2::ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
      geom_point(aes(fill=Site_General), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      xlab(str_PC1) +
      ylab(str_PC2) +
      cowplot::theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    p
  }
}
