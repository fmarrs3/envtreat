



# Generate figures 1 through 4
#    maindir = main parent directory
#    resultsdir = location of results files A.txt, B.txt, and res_counterfactual.RData to read in, below main directory
#    plotdir = location to put plots below the main/parent directory
#
make_figs <- function(maindir, resultsdir="results", plotdir="figures")
{
  alpha <- .05
  options(stringsAsFactors = FALSE)
  
  resultsdir <- file.path(maindir, resultsdir)
  plotdir <- file.path(maindir, plotdir)
  datadir <- file.path(maindir, "data")
  
  
  #### Read in and file management
  if(!dir.exists(plotdir)){dir.create(plotdir)}
  
  # beta_names <- read.table(file.path(resultsdir, "beta.txt"))[,1]
  # beta_names <- beta_names[c(1:19),]
  # 
  betas <- read.table(file.path(resultsdir, "cibeta.txt"), header=TRUE)
  # # 
  # betas$covariate <- beta_names
  # # 
  # betas$signif <- ifelse(betas$pval < alpha, 1, 0)
  betas <- betas[, c("covariate", "estimate", "pval")]
  
  # imprt betas for no net fits
  betas_NN <- read.table(file.path(resultsdir, "beta_nonet.txt"))
  colnames(betas_NN) <- "Logit"

  # compile into xtable format
  covar_table <- cbind(betas, betas_NN)
  rownames(covar_table) <- covar_table$covariate


  # # Get logit results
  # sink("NewBLINFitsTex_08152018.txt")
  # xtable(round(covar_table[,c(2:9)], 2),
  #        caption = "Full BLIN Model Estimates",
  #        label = "results",
  #        auto = TRUE)
  # sink()
  # # 
  sink(file.path(plotdir, "Table1.txt"))
  print(covar_table[, c("estimate", "pval", "Logit")])
  sink()
  

  
  
  amat <- read.table(file.path( resultsdir, "ciA.txt"), header=TRUE)
  amat <- data.frame(Source = amat$cow1,
                     Target = amat$cow2,
                     Weight = amat$estimate,
                     Sig = amat$pval < alpha)
  amat <- na.omit(amat)
  
  bmat <- read.table(file.path( resultsdir, "ciB.txt"), header=TRUE)
  bmat <- data.frame(Source = bmat$treaty1,
                     Target = bmat$treaty2,
                     Weight = bmat$estimate,
                     Sig = bmat$pval < alpha)     
  bmat <- na.omit(bmat)
  
  # let's link these cow codes up to iso3 codes
  iso3 <- countrycode(amat$Source, "cown", "cowc", warn = TRUE)
  # 260 and 818 were not matched
  # 260 is german federal republic
  # 818 is viet nam
  #iso3[which(amat$Source == 260)] <- c("DFR") # order matters
  iso3[which(amat$Source == 818)] <- c("VNM")
  
  amat$Source <- iso3
  
  iso3 <- countrycode(amat$Target, "cown", "cowc", warn = TRUE)
  # 260 and 818 were not matched
  # 260 is german federal republic
  # 818 is viet nam
  #iso3[which(amat$Target == 260)] <- c("DFR") # order matters
  iso3[which(amat$Target == 818)] <- c("VNM")
  
  amat$Target <- iso3
  
  el <- amat
  
  # Write edgelists as .csv's
  # Negative Values -- A Matrix 
  neg_el <- el[el$Weight < 0,]
  neg_el$Weight <- abs(neg_el$Weight)
  
  sig_neg_el <- neg_el[neg_el$Sig == 1,]
  
  # Positive Values -- A Matrix
  pos_el <- el[el$Weight > 0,]
  
  sig_pos_el <- pos_el[pos_el$Sig == 1,]
  
  
  # nice ggplot visualizations of these networks
  
  
  pos_a_net <- network.initialize(n = length(unique(c(sig_pos_el$Source, sig_pos_el$Target))), directed = TRUE)
  network.vertex.names(pos_a_net) <- unique(c(sig_pos_el$Source, sig_pos_el$Target))
  pos_a_net[as.matrix(sig_pos_el[,c(1:2)])] <- 1
  set.edge.value(pos_a_net,"weight",sig_pos_el$Weight)
  set.vertex.attribute(pos_a_net, "outdegree1", as.character(degree(pos_a_net, cmode = "outdegree")))
  ids <- read.csv(file.path(datadir, "Updated_State_IDs.csv"))
  continent <- ids[ids$Id %in%  unique(c(sig_pos_el$Source, sig_pos_el$Target)), c("Id", "continent")]
  continent[continent$Id == "VNM","continent"] <- "Asia"
  continent[continent$Id == "GFR","continent"] <- "Europe"
  continent<-continent[match(get.vertex.attribute(pos_a_net, 'vertex.names'), continent$Id),c(1,2)]
  set.vertex.attribute(pos_a_net, 'Continent', continent$continent)
  
  
  
  posA <- ggnet2(pos_a_net, 
                 size = 14,
                 #  size.palette = c("4" = 200, "3" = 180, "2" = 160, "1" = 140, "0" = 120),
                 #mode = "circle",
                 color = "Continent",
                 color.legend = FALSE,
                 label = TRUE,
                 label.size = 3,
                 edge.color = "grey",
                 edge.size = "weight",
                 palette = "Set2",
                 arrow.size = 3, 
                 arrow.gap = 0.03) +
    guides(size=FALSE)

  ggsave(posA, width=12, height=12, filename = file.path(plotdir, "Fig2.pdf"))
  
  # tiff("Fig2.tiff", height = 800, width = 800)
  # posA
  # dev.off()
  
  
  
  neg_a_net <- network.initialize(n = length(unique(c(sig_neg_el$Source, sig_neg_el$Target))), directed = TRUE)
  network.vertex.names(neg_a_net) <- unique(c(sig_neg_el$Source, sig_neg_el$Target))
  neg_a_net[as.matrix(sig_neg_el[,c(1:2)])] <- 1
  set.edge.value(neg_a_net,"weight",sig_neg_el$Weight)
  set.vertex.attribute(neg_a_net, "outdegree1", as.character(degree(neg_a_net, cmode = "outdegree")))
  # ids <- read.csv("~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/Updated_State_IDs.csv")
  continent <- ids[ids$Id %in%  unique(c(sig_neg_el$Source, sig_neg_el$Target)), c("Id", "continent")]
  continent[continent$Id == "VNM","continent"] <- "Asia"
  continent[continent$Id == "GFR","continent"] <- "Europe"
  continent<-continent[match(get.vertex.attribute(neg_a_net, 'vertex.names'), continent$Id),c(1,2)]
  set.vertex.attribute(neg_a_net, 'Continent', continent$continent)

  
  negA <- ggnet2(neg_a_net, 
                 size = 14,
                 #  size.palette = c("4" = 200, "3" = 180, "2" = 160, "1" = 140, "0" = 120),
                 #mode = "circle",
                 color = "Continent",
                 color.legend = FALSE,
                 label = TRUE,
                 label.size = 3,
                 edge.color = "grey",
                 edge.size = "weight",
                 palette = "Set2",
                 arrow.size = 3, 
                 arrow.gap = 0.03) +
    guides(size=FALSE)
  ggsave(negA, width=12, height=12, filename = file.path(plotdir, "Fig1.pdf"))
  

  # Write edgelists as .csv's
  # Negative Values -- B Matrix 
  el <- bmat
  neg_el <- el[el$Weight < 0,]
  neg_el$Weight <- abs(neg_el$Weight)
  #write.csv(neg_el,  file = "~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/EnvTreaties_BMAT_NewResults_Neg_02132018.csv",row.names = FALSE)
  
  sig_neg_el <- neg_el[neg_el$Sig == 1,]
  
  #write.csv(sig_neg_el, file = "~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/EnvTreaties_BMAT_NewResults_SigNeg_02132018.csv", row.names = FALSE)
  
  
  # Positive Values -- B Matrix
  pos_el <- el[el$Weight > 0,]
  #write.csv(pos_el, file = "~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/EnvTreaties_BMAT_NewResults_Pos_02132018.csv", row.names = FALSE)
  
  sig_pos_el <- pos_el[pos_el$Sig == 1,]
  
  #write.csv(sig_pos_el, file = "~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/EnvTreaties_BMAT_NewResults_SigPos_02132018.csv",row.names = FALSE)
  
  
  # nice ggplot visualizations of these networks
  
  
  pos_b_net <- network.initialize(n = length(unique(c(sig_pos_el$Source, sig_pos_el$Target))), directed = TRUE)
  network.vertex.names(pos_b_net) <- as.character(unique(c(sig_pos_el$Source, sig_pos_el$Target)))
  sig_pos_el$Source <- as.character(sig_pos_el$Source)
  sig_pos_el$Target <- as.character(sig_pos_el$Target)
  pos_b_net[as.matrix(sig_pos_el[,c(1:2)])] <- 1
  set.edge.value(pos_b_net,"weight",sig_pos_el$Weight)
  set.vertex.attribute(pos_b_net, "outdegree1", as.character(degree(pos_b_net, cmode = "outdegree")))
  ids <- read.csv(file.path(datadir, "Treaty_IDs.csv"))
  topic <- ids[ids$Id %in%  unique(c(sig_pos_el$Source, sig_pos_el$Target)), c("Id", "Topic")]
  topic[is.na(topic$Topic),"Topic"] <- "Misc."
  topic<-topic[match(get.vertex.attribute(pos_b_net, 'vertex.names'), topic$Id),c(1,2)]
  set.vertex.attribute(pos_b_net, 'Topic', topic$Topic)
  
  
  posB <- ggnet2(pos_b_net, 
                 size = 14,
                 #  size.palette = c("4" = 200, "3" = 180, "2" = 160, "1" = 140, "0" = 120),
                 #mode = "circle",
                 color = "Topic",
                 color.legend = FALSE,
                 label = TRUE,
                 label.size = 3,
                 edge.color = "grey",
                 edge.size = "weight",
                 palette = "Set2",
                 arrow.size = 3, 
                 arrow.gap = 0.03) +
    guides(size=FALSE)
  
  ggsave(posB, width=12, height=12, filename = file.path( plotdir, "Fig4.pdf"))
  
  
  neg_b_net <- network.initialize(n = length(unique(c(sig_neg_el$Source, sig_neg_el$Target))), directed = TRUE)
  network.vertex.names(neg_b_net) <- as.character(unique(c(sig_neg_el$Source, sig_neg_el$Target)))
  sig_neg_el$Source <- as.character(sig_neg_el$Source)
  sig_neg_el$Target <- as.character(sig_neg_el$Target)
  neg_b_net[as.matrix(sig_neg_el[,c(1:2)])] <- 1
  set.edge.value(neg_b_net,"weight",sig_neg_el$Weight)
  set.vertex.attribute(neg_b_net, "outdegree1", as.character(degree(neg_b_net, cmode = "outdegree")))
  # ids <- read.csv("~/Dropbox/E-IGO Network Data/Data/data_envTreaties/visualizations/Gephi/results02132018/Treaty_IDs.csv")
  topic <- ids[ids$Id %in%  unique(c(sig_neg_el$Source, sig_neg_el$Target)), c("Id", "Topic")]
  topic[is.na(topic$Topic),"Topic"] <- "Misc."
  topic<-topic[match(get.vertex.attribute(neg_b_net, 'vertex.names'), topic$Id),c(1,2)]
  set.vertex.attribute(neg_b_net, 'Topic', topic$Topic)
  
  
  negB <- ggnet2(neg_b_net, 
                 size = 14,
                 #  size.palette = c("4" = 200, "3" = 180, "2" = 160, "1" = 140, "0" = 120),
                 #mode = "circle",
                 color = "Topic",
                 color.legend = FALSE,
                 label = TRUE,
                 label.size = 3,
                 edge.color = "grey",
                 edge.size = "weight",
                 palette = "Set2",
                 arrow.size = 3, 
                 arrow.gap = 0.03) +
    guides(size=FALSE)

  ggsave(negB, width=12, height=12, filename = file.path(plotdir, "Fig3.pdf"))

#   
# }
# 
# 
# 
# 
# 
# 
# # Generate figures 5 through 7
# #    maindir = main parent directory
# #    resultsdir = location of results files A.txt, B.txt, and res_counterfactual.RData to read in, below main directory
# #    plotdir = location to put plots below the main/parent directory
# #
# make_figs_5_7 <- function(maindir, resultsdir, plotdir="figures")
# {
  gpclibPermit()
  
  if(!as.logical(gpclibPermitStatus())){
    if(!as.logical(gpclibPermit())){
      stop("Need to install package gpclib and run gpclibPermit()")
    }
  }
  
  # 
  # resultsdir <- file.path(maindir, resultsdir)
  # plotdir <- file.path(maindir, plotdir)
  # datadir <- file.path(maindir, "data")
  # 
  
  #### Read in and file management
  if(!dir.exists(plotdir)){dir.create(plotdir)}
  cows <- read.csv(file.path(datadir, "COW_country_codes.csv"), header=TRUE)
  load(file.path(resultsdir, "res_counterfactual.RData"))
  A <- read.table(file.path(resultsdir, "A.txt"), header=TRUE)
  B <- read.table(file.path(resultsdir, "B.txt"), header=TRUE)
  ####
  
  
  #### Process data
  delta <- (Ymean0 - Ymean_swap)
  
  
  i <- which(dimnames(Ymean0)[[1]] == cswap)   # USA index
  j <- which(dimnames(Ymean0)[[2]] == tswap)   # UNFCCC index
  
  # cat("Changes in probabilities of the US signing treaties (percentages): \n")
  # print(sort(delta[1,] / Ymean_swap[1,])*100   )
  
  # cat("\nChanges in probabilities of countries signing UNFCCC  (percentages): \n")
  # print( sort(delta[,130] / Ymean_swap[,130])*100 ) 
  
  summary_data <- data.frame(as.matrix( cbind(c(Ymean_swap), c(Ymean0), c(delta), c(delta) / c(Ymean_swap)*100) ))
  names(summary_data) <- c("Yhat_swapped", "Yhat_unswapped", "change", "pct_change")
  summary_data$A <- rep(A[i,], times=ncol(Ymean0))
  summary_data$B <- rep(B[j,], each=nrow(Ymean0))
  summary_data$cowcode <- rep(dimnames(Ymean0)[[1]], times=ncol(Ymean0))
  summary_data$treaty <- rep(dimnames(Ymean0)[[2]], each=nrow(Ymean0))
  summary_data$country <- cows$StateNme[match(summary_data$cowcode, cows$CCode)]
  summary_data <- summary_data[, c("country", "cowcode", "treaty", "change", "pct_change", "A", "B", "Yhat_swapped", "Yhat_unswapped")]
  
  top25 <- head(summary_data[order(summary_data$change, decreasing = T), ], 25)   # top 25 largest changes in probability
  # length(unique(top25$country))
  # length(unique(top25$treaty))
  ####
  
  
  
  
  # Note US > Thailand has zero estimated influence, howver, Thailand increased probability of entry into UNFCCC
  #    triadic closure from US > Mexico > Thailand!
  
  # Cowcodes that influence Thailand positively
  # cbind(rownames(A), A[, "X800"])[A[,"X800"] > 0,]
  
  # Probability changes in these cowcodes with UNFCCC treaty
  # summary_data[summary_data$cowcode %in% rownames(A)[A[,"X800"] > 0] & summary_data$treaty ==tswap, ]
  
  # US most likely to go on to ratify (53%) the United Nations Convention To Combat Desertification In Those Countries Experiencing Serious Drought And/ Or Desertification, Particularly In Africa
  # US least likely to go on to ratify (-33%) the International Convention for the Regulation of Whaling
  # Russia is most likely to go on to ratify UNFCCC (114.32%).
  # Libya is least likely to go on to ratify the UNFCCC (-7.15%)
  
  
  
  
  # countries and treaties of interest
  countries <- c("2", "20", "70", "140", "200", "220", "260", "365", "710", "732", "740")
  treaties <- c("40793", "3273", "40827", "40824", "40828", "40869")
  
  # reduce matrix of changes by countries and treaties of interest
  reduced_delta <- delta[countries,treaties]
  
  # make it into a data frame because it'll be easer to plot and transform
  reduced_delta <- as.data.frame(reduced_delta)
  new_rows <- countrycode::countrycode(rownames(reduced_delta), "cown", "cowc")
  new_rows[is.na(new_rows)] <- c("GFR")
  rownames(reduced_delta)  <- new_rows
  colnames(reduced_delta) <- c("UNFCCC (40793)", "Kyoto (3273)", "Desertification (40827)", "Protocol to UNCOLS (40824)",
                               "Nuclear Safety (40828)", "Anti-Personnel Mines (40869)")
  
  reduced_delta$ID = rownames(reduced_delta)
  
  
  
  #### Figure 7
  # melt this for ggplot formatting
  melted <- melt(reduced_delta, id = "ID")
  melted <- melted[order(melted$ID),]
  
  p <- ggplot(melted, aes(x = variable, y =ID)) + 
    geom_tile(aes(fill = value)) +  
    scale_y_discrete(limits = unique(rev(melted$ID))) + 
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low = "white", high = "#67000D") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size=16, face="bold"),
          panel.spacing.y = grid::unit(1, "cm")) +
    # guides(fill=FALSE) +
    labs(y = "State", x = "Treaty", fill = "Change in Pr(Ratification)")
  
  # ggsave(p, filename = file.path(plotdir, "Fig7.tiff"), device = "tiff")
  
  # p
  ggsave(p, filename = file.path(plotdir, "Fig7.pdf"))
  # tiff(file.path(plotdir, "Fig7.tiff"), height = 800, width = 800)
  # p
  # dev.off()
  ####
  
  
  # #### Need this one??
  # # Plot 2: Top 25 bipartite network
  # top25EL <- data.frame(
  #   Source = top25$cowcode,
  #   Target = top25$treaty,
  #   Weight = top25$change
  # )
  # 
  # top25Ids <- data.frame(
  #   Id = c(top25$cowcode, top25$treaty),
  #   label = c(countrycode::countrycode(top25$cowcode, "cown", "cowc"), top25$treaty),
  #   type = c(rep("State", 25), rep("Treaty", 25))
  # )
  # 
  # top25Ids <- unique(top25Ids)
  # 
  # # these were then put into gephi
  # write.csv(top25EL, row.names = FALSE, file = file.path(plotdir, "top25el.csv"))
  # write.csv(top25Ids, row.names = FALSE, file = file.path(plotdir, "top25ids.csv"))
  # ####
  
  
  
  #### Fig. 5
  # install.packages('rgeos', type='source')
  # install.packages('rgdal', type='source')
  
  data(wrld_simpl)
  
  reduced_df <- summary_data[summary_data$treaty == "40793",]
  ddf <- data.frame(
    country = reduced_df$country,
    value = reduced_df$change
  )
  
  plotme <- function() {
    
    # this lets us use the contry name vs 3-letter ISO
    wrld_simpl@data$id <- wrld_simpl@data$NAME
    
    wrld <- fortify(wrld_simpl, region="id")
    wrld <- subset(wrld, id != "Antarctica") # we don't rly need Antarctica
    
    gg <- ggplot()
    
    # setup base map
    gg <- gg + geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="white", color="#7f7f7f", size=0.25)
    
    # add our colored regions
    gg <- gg + geom_map(data=ddf, map=wrld, aes(map_id=country, fill=value),  color="black", size=0.05)
    
    # this sets the scale and, hence, the legend
    pal <- colorRampPalette(brewer.pal(n = length(unique(ddf$value)), 'Reds'))(length(ddf$value))
    
    # palate of blues for negs
    #  palneg <- rev(colorRampPalette(brewer.pal(n = 6, 'Blues'))(12))
    # pal[1:12] <- palneg
    palSz <- length(pal)  # not sure what you really want/need for this range
    
    
    gg <- gg + scale_fill_gradient2(low = pal[1],
                                    mid = pal[palSz/2],
                                    high = pal[palSz],
                                    midpoint = (max(ddf$value) + min(ddf$value)) / 2,
                                    name="Change in Pr(Ratify UNFCCC)")
    
    # gg <- gg + scale_fill_gradientn(colours = pal, name = "Change in Pr(Ratify UNFCCC)")
    # this gives us proper coords. mercator proj is default
    gg <- gg + coord_map()
    gg <- gg + labs(x="", y="")
    gg <- gg + theme(plot.background = element_rect(fill = "transparent", colour = NA),
                     panel.border = element_blank(),
                     panel.background = element_rect(fill = "transparent", colour = NA),
                     panel.grid = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     legend.position = "right",
                     legend.title = element_text(face="bold"))
    gg <- gg + guides(fill = guide_colorbar(barwidth = 1, barheight = 3, nbin = 186, ticks = FALSE))
    gg
    
  }
  
  map <- plotme()
  # map
  # ggsave(map, filename = "~/Dropbox/E-IGO Network Data/Application Paper Draft/Nature Human Behavior/UNFCCMap_01142019.pdf")
  ggsave(map, filename = file.path(plotdir, "Fig5.pdf"))
  
  # tiff(file.path(plotdir, "Fig5.tiff"), height = 800, width = 800)
  # map
  # dev.off()
  ####
  
  
  
  #### Fig 6
  # get the cases for the treaty of interest
  no_direct_influence <- summary_data[summary_data$treaty == "40793",]
  no_direct_influence <- no_direct_influence[!is.na(no_direct_influence$change),]
  no_direct_influence <- no_direct_influence[order(no_direct_influence$change, decreasing = T),]
  
  
  # turn the a adj into a matrix
  A <- as.matrix(A)
  
  # rename
  colnames(A) <- rownames(A)
  inf_net <- as.network(A)
  
  # get matrices of geodist
  dists <- geodist(A)$gdist
  
  # should be asymmetric, likely
  isSymmetric(dists) # good
  
  
  colnames(dists) <- rownames(A)
  rownames(dists) <- rownames(A)
  
  # put the adj matrix of distances into an edge list i can merge with the no_direct_influence df
  el <- melt(dists)
  colnames(el) <- c("Sender", "cowcode", "dist")
  
  # only out us influence relationships
  el <- el[el$Sender == 2,]
  
  no_direct_influence_merge <- merge(no_direct_influence, el, by = "cowcode")
  
  # 3 char for plot
  # fix vietnam
  no_direct_influence_merge$cowcode[which(no_direct_influence_merge$cowcode == 818)] <- 816 
  no_direct_influence_merge$cowc <- countrycode::countrycode(no_direct_influence_merge$cowcode, "cown", "cowc")
  
  # remove us to us
  no_direct_influence_merge <- no_direct_influence_merge[no_direct_influence_merge$dist != 0,]
  
  ks.test(no_direct_influence_merge$pct_change[no_direct_influence_merge$dist == 1], no_direct_influence_merge$pct_change[no_direct_influence_merge$dist == 2], alternative = "less")
  ks.test(no_direct_influence_merge$pct_change[no_direct_influence_merge$dist == 2], no_direct_influence_merge$pct_change[no_direct_influence_merge$dist == 3], alternative = "less")
  ks.test(no_direct_influence_merge$pct_change[no_direct_influence_merge$dist == 1], no_direct_influence_merge$pct_change[no_direct_influence_merge$dist %in% c(2,3)], alternative = "less")
  
  dist_plot <- ggplot(data = no_direct_influence_merge, aes(fill = as.factor(dist), x = as.factor(dist), group = as.factor(dist), y = change))+
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    ylab("Change in Pr(Ratify UNFCCC)") +
    xlab("Geodesic Distance from US") +
    annotate("text", label = "Dist 1 > Dist 2 \n KS p = 0.000", x = 1.5, y = 0.3) +
    annotate("text", label = "Dist 2 > Dist 3 \n KS p = 0.591", x = 2.5, y = 0.25) +
    annotate("text", label = "Dist 1 > Dists 2 & 3 \n KS p = 0.000", x = 2, y = 0.4) +
    geom_hline(yintercept=0, linetype = 3) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size=16, face="bold"),
          panel.spacing.y = grid::unit(1, "cm")) +
    guides(fill=FALSE)

  ggsave(dist_plot, filename = file.path(plotdir, "Fig6.pdf"))
  
  # tiff(file.path(plotdir, "Fig6.tiff"), height = 800, width = 800)
  # dist_plot
  # dev.off()
  ####
  
}


