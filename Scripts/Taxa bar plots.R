##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################

#Load libraries
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(ggforce)
library(fantaxtic)
library(readxl)
library(textshape)
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(gridExtra)
library(grid)
library(vegan)
library(ggpubr)

#Barplots using the previously generated phyloseq object

#Analyze the prevalence of the reads at phylum level and remove the low abuandant reads (less than 5 reads)
prev=apply(X = otu_table(ambysphy),
           MARGIN = ifelse(taxa_are_rows(ambysphy), yes = 1, no = 2),
           FUN = function(x){sum(x > 0)})

prev=data.frame(Prevalence = prev,
                TotalAbundance = taxa_sums(ambysphy),
                tax_table(ambysphy))

prevt=plyr::ddply(prev, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) -> prev1

kable(prev1)

filtrado=c("AncK6","Altiarchaeota","Atribacteria","Caldiserica","Calditrichaeota","Desantisbacteria",
              "PAUC34f","Schekmanbacteria","Centrohelida","Desantisbacteria")

(ambysphy=subset_taxa(ambysphy, !Phylum %in% filtrado))

ambysphy=phyloseq::prune_samples(sample_sums(ambysphy) > 10000, ambysphy)#filter samples with less that 10000 reads
ambysphy1=prune_taxa(taxa_sums(ambysphy1) > 0, ambysphy1)#filter taxa with 0 counts
taxa_names(ambysphy1=paste0("ASV", seq(ntaxa(ambysphy1)))
           
collapsedtaxa = tax_glom(ambysphy1,taxrank = "Family")#Collapse ASV at family level this step takes a lot of time(1-2H)
sample_data(collapsedtaxa)

ambys2merged2=merge_samples(collapsedtaxa, "SampleType",fun = mean)#merge samples by sample type
sample_data(ambys2merged2)$SampleType=levels(sample_data(collapsedtaxa)$SampleType)
sample_data(ambys2merged2)

barall = transform_sample_counts(ambys2merged2, function(x) 100 * x/sum(x))#Trasnform sample counts

top10all=get_top_taxa(physeq_obj = barall, n = 10, relative = T,
                         discard_other = F, other_label = "Other")#get top taxa
top10all=name_taxa(top10all, label = "", species = T, other_label = "Other")


#Cereate a fucntion for labels
make_labelstring=function(mypanels) {
  mylabels=sapply(mypanels, 
                     function(x) {LETTERS[which(mypanels == x)]})
  
  return(mylabels)
}
label_panels=ggplot2::as_labeller(make_labelstring)


#Here i made some minor modifications to the fantaxtic function of https://github.com/gmteunisse/Fantaxtic
mane_bars <- function(physeq_obj, color_by, label_by = NULL, facet_by = NULL,
                      grid_by = NULL, facet_type = "wrap", bar_width = 0.6,
                      facet_cols =  1, gen_uniq_lbls = TRUE, other_label= NULL,
                      order_alg = "aplph", color_levels = NULL,
                      base_color = "#6495ed",
                      other_color = "#f3f3f3", palette = NULL){
  
  #Check for subcoloring
  if (is.null(label_by)){
    label_by <- color_by
  }
  
  #Extract tax_tbl and add OTU names
  tax_tbl <- as.data.frame(phyloseq::tax_table(physeq_obj))
  tax_tbl$otu_name <- row.names(tax_tbl)
  
  #Replace NAs with Unknown
  tax_tbl <- as.data.frame(apply(tax_tbl, 2, function(x){
    x[is.na(x)] <- "Unknown"
    return(x)
  }))
  
  #Move Other taxa to the beginning and alter taxonomic annotations
  #of Other taxa
  if(!is.null(other_label)){
    main_ind <- which(!tax_tbl[[label_by]] %in% other_label)
    other_ind <- which(tax_tbl[[label_by]] %in% other_label)
    new_color_by <- as.character(tax_tbl[[color_by]])
    new_color_by[other_ind] <- as.character(tax_tbl[[label_by]][other_ind])
    tax_tbl[[color_by]] <- as.factor(new_color_by)
    ordr <- c(other_ind, main_ind)
    tax_tbl <- tax_tbl[ordr,]
  }
  
  #Refactor for legend ordering and order
  if (is.null(color_levels)){
    tax_levels <- unique(tax_tbl[[color_by]])
  } else {
    if (is.null(other_label)){
      tax_levels <- color_levels
    } else {
      tax_levels <- c(other_label, color_levels)
    }
  }
  tax_tbl[[label_by]] <- factor(tax_tbl[[label_by]], unique(tax_tbl[[label_by]]), ordered = T)
  tax_tbl[[color_by]] <- factor(tax_tbl[[color_by]], tax_levels, ordered = T)
  tax_tbl <- tax_tbl[order(tax_tbl[[color_by]]),]
  
  #Get the tax and OTU tables
  otu_tbl <- as.data.frame(phyloseq::otu_table(physeq_obj))
  
  #Check the orientation of the otu_tbl and change if required
  if (!taxa_are_rows(phyloseq::otu_table(physeq_obj))){
    otu_tbl <- as.data.frame(t(otu_tbl))
  }
  
  #Calculate the number of colors and color variations required
  clr_tbl <- as.data.frame(table(tax_tbl[[color_by]], useNA = "ifany"), stringsAsFactors = F)
  if(!is.null(other_label)){
    ind <- which(!clr_tbl$Var1 %in% other_label)
    clr_tbl <- clr_tbl[ind,]
  }
  
  #Generate the required color palette
  clr_pal <- gen_palette(clr_tbl = clr_tbl, clr_pal = palette, base_clr = base_color)
  names(clr_pal) <- clr_tbl$Var1
  clr_pal <- as.vector(unlist(clr_pal))
  if(!is.null(other_label)){
    n_other <- length(other_label)
    other_pal <- gen_shades_tints(n_other, clr = other_color)
    clr_pal <- c(other_pal, clr_pal)
  }
  
  #Generate unique label names if required
  if(gen_uniq_lbls){
    tax_tbl[[label_by]] <- gen_uniq_lbls(tax_tbl[[label_by]])
  }
  
  #Transform absolute taxon counts to relative values
  otu_tbl <- as.data.frame(apply(otu_tbl, 2, function(x){
    if (sum(x) > 0){x/sum(x)}
    else(x)
  }))
  
  #Match order of tax.tbl and otu.tbl
  ord <- match(tax_tbl$otu_name, row.names(otu_tbl))
  otu_tbl <- otu_tbl[ord,]
  
  #Order the samples according to the specified algorithm
  #Order according to selected taxonomies
  if (sum(order_alg %in% c("alph", "hclust", "as.is")) == 0){
    
    #Get the summed abundances
    sums <- list()
    i <- 0
    for (lvl in order_alg){
      i <- i + 1
      sums[[i]] <- round(colSums(otu_tbl[which(tax_tbl[[label_by]] == lvl),]), digits = 3)
    }
    
    #Sort
    cmd <- paste(sprintf("sums[[%d]]", 1:i), collapse = ", ")
    smpl_ord <- eval(parse(text = sprintf("order(%s)", cmd)))
    otu_tbl <- otu_tbl[,smpl_ord]
    
    #Order according to selected algorithm
  }else{
    if (order_alg == "alph"){
      otu_tbl <- otu_tbl[,order(names(otu_tbl))]
    } else {
      if(order_alg == "hclust"){
        hc <- hclust(dist(x = t(otu_tbl), method = "euclidian", upper = F))
        smpl_ord <- hc$order
        otu_tbl <- otu_tbl[,smpl_ord]
      } else {
        if (order_alg == "as.is"){
          #do nothing
        }
      }
    }
  }
  
  
  #Join labels and counts and transform to a long data format
  counts <- cbind(tax_tbl[[color_by]], tax_tbl[[label_by]], otu_tbl)
  names(counts) <- c("color_by", "label_by", colnames(otu_tbl))
  counts_long <- reshape2::melt(counts,
                                id.vars = c("color_by", "label_by"),
                                variable.name = "Sample",
                                value.name = "Abundance")
  
  #Add facet levels if needed and transform to a long data format
  if(is.null(facet_by) & !is.null(grid_by)){
    facet_by <- grid_by
    grid_by <- NULL
  }
  if (!is.null(facet_by)){
    facet <- as.data.frame(phyloseq::sample_data(physeq_obj))[[facet_by]]
    names(facet) <- row.names(phyloseq::sample_data(physeq_obj))
    ord <- match(counts_long$Sample, names(facet))
    facet <- facet[ord]
    counts_long$facet <- facet
  }
  if (!is.null(grid_by)){
    grid <- as.data.frame(phyloseq::sample_data(physeq_obj))[[grid_by]]
    names(grid) <- row.names(phyloseq::sample_data(physeq_obj))
    ord <- match(counts_long$Sample, names(grid))
    grid <- grid[ord]
    counts_long$grid <- grid
  }
  
  #Generate a plot
  p <- ggplot2::ggplot(counts_long, aes(x = Sample, y = Abundance, fill = label_by)) +
    ggplot2::geom_bar(position = "stack", stat = "identity", width = bar_width) +
    ggplot2::scale_fill_manual(values = clr_pal) +
    ggplot2::guides(fill=guide_legend(title = label_by, ncol = 2)) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_discrete(labels = function(x) {sub("\\-", "\n", x)})+
    ggplot2::theme(axis.line.x = element_line(colour = 'black'),
                   axis.line.y = element_line(colour = 'black'),
                   axis.ticks = element_line(colour = 'black'),
                   axis.text.x = element_text(angle = 0,face = "bold",
                                              size = 12, hjust = 0.5,vjust = 0.5,colour = "black"),
                   legend.background = element_rect(fill = 'transparent', colour = NA),
                   legend.key = element_rect(fill = "transparent"),
                   legend.key.size = unit(0.2, "cm"),
                   legend.title = element_text(size = 11,face = "bold"),
                   legend.position = "bottom",
                   legend.text = element_text(size = 11,face = "bold"),
                   panel.background = element_rect(fill = 'transparent', colour = NA),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_line(colour = adjustcolor('black', 0.3)),
                   panel.grid.minor = element_line(colour = NA),
                   plot.background = element_rect(fill = 'transparent', colour = NA),
                   plot.title = element_text(color = "white"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   strip.background = element_blank(),
                   strip.text = element_text(size = 12, face = "bold",hjust = -0),
                   text = element_text(size = 12))
  
  if (!is.null(facet_by)) {
    if (facet_type == "wrap"){
      if (is.null(grid_by)){
        p <- p + ggplot2::facet_wrap(~facet, scales = "free", ncol = facet_cols,labeller = label_panels)
      }else{
        p <- p + ggplot2::facet_wrap(~grid + facet, scales = "free", ncol = facet_cols,labeller = label_panels)
      }
    }else{
      if (facet_type == "grid"){
        if (is.null(grid_by)){
          p <- p + ggplot2::facet_grid(~facet, scales = "free", space = "free",labeller = label_panels)
        }else{
          p <- p + ggplot2::facet_grid(facet ~ grid, scales = "free", space = "free",labeller = label_panels)
        }
      }
    }
  }
  
  return(p)
}

#Generate a color palette
mycols<-c("#A6CEE3","#1F78B4","#FB9A99", "#E31A1C", 
          "#CAB2D6", "#6A3D9A","#B2DF8A", 
          "#33A02C","#666666", "#FF7F00","#FDBF6F")

#Plot the rel abund bars
ptop10all<-mane_bars(top10all, color_by = "Family", label_by = "Family",
                     grid_by = NULL, other_color = "Grey",gen_uniq_lbls = FALSE,palette = mycols)


