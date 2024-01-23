#' Generate taxa summary plots by genera showing only the top 0.1% most abundant 
#'
#' 
#' Plot does not have legend - legends were graphed separately 
#' 
#' 
#' @author Julianne Yang
#' @param path_to_RDS filepath to where your L6 RDS is stored
#' @param titlestring title of plot
#' @param greppattern ".*g__" for sgenera
#' @param fillvector a named vector of colors 
#' @param graphby a string, can pass "Site" or "Type"
#' @return a ggplot2 object which you can then further customize by adding ggplot2 functions
#' @export
#' @examples
#' 
#' generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS", 
#' "Luminal ( > 0.1% Relative Abundance)",
#'  ".*g__",
#'   assign_cols,
#'    "Site")
#' 
#' 

generate_L6_taxa_plots <- function(path_to_RDS, titlestring,greppattern, fillvector, graphby){
  #L2_lum<-readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/SI_LumMuc_L6.RDS")
  #taxa <- gsub(".*g__","",taxa)
  #cols<-assign_cols
  titlestring<-c(titlestring)
  L2_lum<-readRDS(path_to_RDS)
  L2_lum<- as.matrix(L2_lum)
  L2_lum<-funrar::make_relative(L2_lum)
  L2_lum<-as.data.frame(t(L2_lum))
  toptaxa<- rowMeans(L2_lum)
  L2_lum$averageRA <-toptaxa/6
  L2_lum <- L2_lum %>% dplyr::mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
  L2_lum <-select(L2_lum,-averageRA)
  
  taxa<-L2_lum$keeptaxa
  L2_lum <- select(L2_lum,-keeptaxa)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  }
  
  else if ({{graphby}}=="SiteTypeSI"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Duodenum_L"="Duo_L", "Duodenum_M" = "Duo_M",
                                               "Jejunum_L" ="Jej_L", "Jejunum_M" ="Jej_M",
                                               "Ileum_L"="Ile_L", "Ileum_M" = "Ile_M"))
    L2_lum$Site = factor(L2_lum$Site, c("Duo_L", "Duo_M",
                                        "Jej_L", "Jej_M",
                                        "Ile_L", "Ile_M"))
  }
  else if ({{graphby}}=="SiteTypeColon"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Cecum_L"="Cec_L","Cecum_M" ="Cec_M",
                                               "Proximal_Colon_L"="PC_L","Proximal_Colon_M"="PC_M",
                                               "Distal_Colon_L"="DC_L", "Distal_Colon_M" ="DC_M"))
    L2_lum$Site = factor(L2_lum$Site, c("Cec_L","Cec_M",
                                        "PC_L","PC_M",
                                        "DC_L", "DC_M"))
  }
  
  cols <- fillvector
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme(legend.position = "none")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

#' Generate taxa summary plots by phyla showing only the top 0.1% most abundant 
#'
#' 
#' Plot does not have legend - legends were graphed separately 
#' 
#' 
#' @author Julianne Yang
#' @param input_data filepath to L2 data .csv (counts)
#' @param titlestring title of plot
#' @param greppattern ".*s__" for species 
#' @param graphby a string, can pass "Site" or "Type"
#' @return a ggplot2 object which you can then further customize by adding ggplot2 functions
#' @export
#' @examples
#' 
#' generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS", 
#' "Luminal ( > 0.1% Relative Abundance)",
#'  ".*g__",
#'   assign_cols,
#'    "Site")
#' 
#' 
generate_L2_taxa_plots <- function(input_data, titlestring,greppattern, fillvector,graphby){
  titlestring<-c(titlestring)
  L2_lum<-read.csv(input_data)
  L2_lum <- as.data.frame(t(L2_lum))
  colnames(L2_lum)<- L2_lum[1,]
  L2_lum <- L2_lum[-1,]
  taxa<-row.names(L2_lum)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  if({{graphby}} == "Site"){
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  }
  else if ({{graphby}}=="SiteTypeSI"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Duodenum_L"="Duo_L", "Duodenum_M" = "Duo_M",
                                               "Jejunum_L" ="Jej_L", "Jejunum_M" ="Jej_M",
                                               "Ileum_L"="Ile_L", "Ileum_M" = "Ile_M"))
    L2_lum$Site = factor(L2_lum$Site, c("Duo_L", "Duo_M",
                                        "Jej_L", "Jej_M",
                                        "Ile_L", "Ile_M"))
  }
  else if ({{graphby}}=="SiteTypeColon"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Cecum_L"="Cec_L","Cecum_M" ="Cec_M",
                                               "Proximal_Colon_L"="PC_L","Proximal_Colon_M"="PC_M",
                                               "Distal_Colon_L"="DC_L", "Distal_Colon_M" ="DC_M"))
    L2_lum$Site = factor(L2_lum$Site, c("Cec_L","Cec_M",
                                        "PC_L","PC_M",
                                        "DC_L", "DC_M"))
  }
  L2_lum$Value <- L2_lum$Value * 100
  
  cols <- fillvector
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values = cols)+
    theme(legend.position = "right")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="")+
    ggtitle({{titlestring}}) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
  
}

process_taxonomy_data <- function(file_path) {
  input_data <- readr::read_csv(here(file_path))
  input_data <- as.data.frame(input_data)
  row.names(input_data)<- input_data[,1]
  input_data <- input_data[,-1]
  
  
  taxa<-colnames(input_data)
  taxa <- gsub(";",".",taxa)
  colnames <- strsplit(taxa, ".o__")
  
  order=new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(colnames)) {
    order[i] <- colnames[[i]][2]
    i=i+1
  }
  
  order<-unlist(order)
  order <- strsplit(order, ".f__")
  
  family =new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(order)) {
    family[i] <- order[[i]][2]
    i=i+1
  }
  
  order<-as.list(order)
  
  i=1
  for (i in 1:length(family)) {
    if (isFALSE(family[[i]]==".g__")) {
      family[[i]] = family[[i]] 
    }
    else {
      family[[i]] <- paste0(order[[i]]," (o)")   
      family[[i]] <- family[[i]][1]
    }
    i=i+1
  }
  
  
  family<-unlist(family)
  family <- strsplit(family, ".g__")
  
  genus =new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(family)) {
    genus[i] <- family[[i]][2]
    i=i+1
  }
  
  family<-as.list(family)
  
  i=1
  for (i in 1:length(genus)) {
    if (isFALSE(genus[[i]]=="NA")) {
      genus[[i]] = genus[[i]] 
    }
    else {
      
      genus[[i]] <- paste0(family[[i]]," (f)")   
    }
    i=i+1
  }
  
  
  colnames(input_data) <-as.character(genus)
  
  input_data <- input_data[,!grepl("Mitochondria", colnames(input_data))] 
  input_data <- input_data[,!grepl("Chloroplast", colnames(input_data))] 
  return(input_data)
}
