
requiredPackages = c('ggplot2', 'ggnewscale', 'cowplot',
                     'seqinr','stringr', 'data.table',
                     "dplyr", "tidyr", "Biostrings",
                     'grid','gridExtra','GetoptLong',
                     'ggseqlogo', 'tidyverse') 

suppressMessages(
  for (p in requiredPackages) {
    # if (!require(p, character.only = TRUE))
    #    install.packages(p)
    library(p, character.only = TRUE)
  }
)

# Functions
# map global msa to local msa (processed by FrustraEvo)
map_positions_msa = function(seq_ali, seq_subali) {
  # remove gaps from input sequences
  s1 = c2s(str_remove(seq_ali, '-'))
  s2 = c2s(str_remove(seq_subali, '-'))
  # check the position where the subalignment sequences starts
  loc = str_locate(s1, s2)
  # create a data frame containing AAs and the position of the protein sequence in the
  # global MSA
  sub_df = data.frame(aa = s2c(s2),
                      prot_pos = seq(loc[, 'start'],
                                     length(s2c(s2)) + loc[, 'start'] - 1,
                                     1))
  sub_df$sub_prot_pos = 1:nrow(sub_df)
  # create a data frame containing AAs and the position of the protein sequence in the
  # subMSA
  #submsa_pos_df = data.frame (AA = seq_subali,
  #                            subali_pos = 1:length(seq_subali))
  submsa_pos_df <- data.frame(
    AA = sapply(seq_subali, toString),
    subali_pos = 1:length(seq_subali)
  )
  submsa_pos_df = cbind(subset(submsa_pos_df, AA != '-'), sub_df)
  
  
  ddf = data.frame(aa = s2c(s1),
                   prot_pos = 1:length(s2c(s1)))
  
  #msa_pos_df = data.frame (aa = seq_ali,
  #                         ali_pos = 1:length(seq_ali))
  msa_pos_df <- data.frame(
    aa = sapply(seq_ali, toString),
    ali_pos = 1:length(seq_ali)
  )
  msa_pos_df = cbind(subset(msa_pos_df, aa != '-'), ddf)
  
  
  mapping_df = right_join(submsa_pos_df, msa_pos_df[, c(1, 2, 4)], by = c("prot_pos", "aa"))
  
  return(mapping_df)
}


map_pdb_to_ali = function(pdb_seq, prot_seq, frust_res, seqid) {
  loc = str_locate(c2s(pdb_seq), c2s(prot_seq))
  if (any(is.na(loc))) {
    PWA  = pairwiseAlignment(subject = c2s(pdb_seq),
                             pattern = c2s(prot_seq))
    loc[, 'start'] = PWA@subject@range@start
    loc[, 'end'] = PWA@subject@range@width
  }
  ff = frust_res[loc[, "start"]:loc[, "end"],]
  colnames(ff)[colnames(ff) == "Res"] = "pdb_pos"
  ff$sub_prot_pos = 1:nrow(ff)
  ff$seqid = seqid
  return(ff)
}


mapping_all = function(msas, msa_glob, frust_files, input_dirs) {
   # list to store mapping data frames
   frustic_seqic_list = list()
   seqic = list()
   ## iterate over msas to get the SRFI indexes
   for (m in 1:length(msas)) {
     msa = msas[[m]]
     # get seqic
     slog = ggseqlogo(unlist(getSequence(getFrag(msa, 1,length(msa[[1]])), as.string = T)),
                      ncol = 100,
                      font= "akrobat_bold")+
       theme(axis.text.x = element_blank(),
         axis.text.y = element_text())+
         scale_x_continuous(breaks = c(1,seq(5,2998, by = 5)), expand = c(0, 0))
     
     LOGO_DATA = ggplot_build(slog)$data[[1]]
     LOGO_DATA$sub_prot_pos = round(LOGO_DATA$x)
     seqic[[m]] = LOGO_DATA %>% group_by(sub_prot_pos) %>% summarise(Entropy = max(y, na.rm = T))
     
     
     # for each sequence in the subalignment, get the mapping of positions
     # between the global and the local alignment and the frustratometeR results
     for (seqid in names(msa)) {
       # extract the sequence from the msa and submsa
       subali_seq = msa[str_detect(names(msa), seqid)][[1]]
       ali_seq = msa_glob[str_detect(names(msa_glob), seqid)][[1]]
       # get the mapping between msa and submsa for the seqid
       mm = map_positions_msa(ali_seq, subali_seq)
       mm$seqid = seqid
       mm$ref = seqrefs[m]
       # remove gaps from submsa sequence
       prot_seq = str_remove(subali_seq, '-')
       # read frustratometer results linked to seqid
       #sel_frust_file = Sys.glob(file.path(input_dirs[m], 'Data', paste(seqid, '.done', sep=""),'FrustrationData', '*singleresidue'))
       sel_frust_file = Sys.glob(file.path(input_dirs[m], 'FrustraEvo_*', 'Data', paste(seqid, '.done', sep=""),'FrustrationData', '*singleresidue'))
       frust_res = fread(sel_frust_file)
       # get the pdb sequence
       pdb_seq = frust_res$AA
       # mapping between submsa and frustratometer results
       ff = map_pdb_to_ali(pdb_seq,  prot_seq, frust_res, seqid)
       ff$ref = seqrefs[m]
       # join results
       m2 = left_join(mm, ff, by = c("sub_prot_pos", "ref", "AA", "seqid"))
       m2$cluster = cluster_names[m]
       # add mapping data frame to list of results
       frustic_seqic_list[[paste(m,seqid, sep ="_")]] = m2
     }
   }
   return(list(rbindlist(frustic_seqic_list), seqic))
 }

# input arguments

# GetoptLong(
#   "msa_glob=s@", "List of FrustraEvo output folders.",
#   "input=s@", "Directory containing list of FrustraEvo output folders.",
#   "cluster_names=s@", "Cluster names of each group. By default 'group_1', 'group_2',..., 'group_n'.",
#   "width_plot=s@", "Width of the plot. Default = 80."
# )


# input arguments
msa_glob = "/.../global/filtered_globalmsa.fasta"
input = "/.../Inputs_FrustraEvo"
cluster_names = c("Family_target0.1", " Single_target0.1", "  Native") #as many names as sets you are analysing
width_plot = 50

if (length(cluster_names) < 1){
  cluster_names = paste('group', seq(1:length(input)), sep = "_")
}

# Get paths to necesary files
input_dirs = Sys.glob(file.path(input,'*'))
msas_files = Sys.glob(file.path(input,'*', "FrustraEvo_*", "OutPutFiles","MSA*.fasta"))
frustic_files = Sys.glob(file.path(input,'*', "FrustraEvo_*", 'OutPutFiles', 'IC_SingleRes_*'))
frustic_files = frustic_files[!str_detect(frustic_files, '.txt')]
#seqic_files = Sys.glob(file.path(input,'*', 'OutPutFiles', 'SeqIC_*.tab'))
frust_files = Sys.glob(file.path(input,'*', "FrustraEvo_*", 'Data','*.done','FrustrationData', '*.pdb_singleresidue'))

# load files
msas = lapply(msas_files, read.fasta, forceDNAtolower = F, seqtype = "AA")
# Load global msa #el codigo original de vicky olvidaba leer esto
msa_glob = read.fasta(msa_glob, forceDNAtolower = FALSE, seqtype = "AA")
frustic = lapply(frustic_files, fread)
seqic = list()#lapply(seqic_files, fread)

# get sequence of reference per protein group
seqrefs = unlist(lapply(frustic, function(x) unique(x$Prot_Ref)))

# get dataframe with positions of reference in msa and submsa
list_mapping = mapping_all(msas, msa_glob, frust_files, input_dirs)
df_mapping = list_mapping[[1]]
seqic = list_mapping[[2]]
# calculate median SRFI
median_srfi = df_mapping %>% group_by(cluster, ref, ali_pos, subali_pos) %>% summarise(median_SRFI = median(FrstIndex, na.rm=T))
# extract data for reference proteins
refs_srfi = subset(df_mapping, seqid %in% seqrefs ) 
# get consensus sequence in data frame
msas2 = lapply(msas_files, read.alignment, forceToLower = F, seqtype = "AA", format = "FASTA")
cmsa = lapply(msas2, consensus)
con = lapply(cmsa, function(c) data.frame(consensus=c, subali_pos = 1:length(c)) )
all <- do.call("rbind", con)
all$ref <- rep(seqrefs, lapply(con, nrow))
all$cluster <- rep(cluster_names, lapply(con, nrow))
# add consensus sequence to median SRFI calculations
con_median_srfi = left_join(median_srfi, all, by = c("subali_pos", "ref", "cluster"))
con_median_srfi = left_join(con_median_srfi, refs_srfi, by = c("ref", "ali_pos", "subali_pos", "cluster"))

# processing frustic data
frustic_df <- do.call("rbind", frustic)
frustic_df$cluster <- rep(cluster_names, lapply(frustic, nrow))
colnames(frustic_df)[colnames(frustic_df)=="Res"] = "sub_prot_pos"
colnames(frustic_df)[colnames(frustic_df)=="Num_Ref"] = "pdb_pos"
colnames(frustic_df)[colnames(frustic_df)=="Prot_Ref"] = "ref"


# processing of seqic data 
seqic_df <- do.call("rbind", seqic)
#le cambio el nombre a la variable entropia pq es SeqIC envd
index_entropy <- which(colnames(seqic_df) == "Entropy")
colnames(seqic_df)[index_entropy] <- "SeqIC"

seqic_df$ref <- rep(seqrefs, sapply(seqic, nrow))
seqic_df$cluster <- rep(cluster_names, lapply(seqic, nrow))

#colnames(seqic_df)[colnames(seqic_df)=="Position"] = "sub_prot_pos"

# merge frustic and seqic data
frustic_seqic_df = merge(frustic_df, seqic_df, by = c("sub_prot_pos", "ref","cluster"))

# merge median SRFI, frustic and seqic data that contains the positions in the alignemnt
# needed to generate the consensus MSFA
plot_df = left_join(con_median_srfi, frustic_seqic_df, by = c("ref", "sub_prot_pos", "pdb_pos", "cluster"))
# determine the width per row of the plot
i = seq(0,max(plot_df$ali_pos), as.numeric(width_plot))+1
j = c(seq(as.numeric(width_plot) ,max(plot_df$ali_pos), as.numeric(width_plot)), max(plot_df$ali_pos))


# Median plot
# create the plot(s)
list_plots = list()
for (k in 1:length(i)){
  # subset data referring to just a certain width
  s =subset(plot_df, ali_pos %in% i[k]:j[k] & !is.na(SeqIC))
  # generate plot
  g = ggplot(s,
             aes(x = ali_pos, y = cluster)) +
    geom_tile(
      aes(fill = median_SRFI),
      linejoin = "round",
      width = 0.96,
      height = 0.96,
      color = 'grey20'
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label  = consensus, size = SeqIC), fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_fill_gradientn(
      values = c(1, .79, .5, 0.11, 0.0909, 0),
      colours = c("green", "grey", "grey", "grey", "red", "red"),
      name = "SRLF median index",
      limits = c(-1.2, 1.1),
      oob = scales::squish,
      na.value = "white"
    ) +
    scale_color_gradient2(
      mid = "grey",
      high = "black",
      midpoint = 0.5,
      limits = c(0, log2(20)),
      oob = scales::squish
    ) +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "top")
  # store plot
  list_plots[[k]] = g + theme(legend.position = "none")
}

# combine plots
combined_plot<-plot_grid(plotlist =list_plots,ncol=1)
# extract legend to have a single legend
legend1 <- get_legend(
  g +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

# add legend to combined plot and legend using plot_grid()
p = plot_grid(NULL, legend1,NULL,   combined_plot, ncol=1,rel_heights = c(.1,.2,.1, 1))
# add x and y labels common for all the plots and not by row of the plot
y.grob <- textGrob("Cluster", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13), rot=90)
x.grob <- textGrob("Global MSA position", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13))


# save figure adjusting for height and width proportionally
height = length(cluster_names) + length(i) -1
svg("consensus_MSFA_mediant.svg", width = 15, height = height )
pp = grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

dev.off()

##########################
# Contribution plot
# # create the plot(s)
list_plots = list()
for (k in 1:length(i)) {
  # subset data referring to just a certain width
  s = subset(plot_df, ali_pos %in% i[k]:j[k] & !is.na(SeqIC))
  # generate plot
  g = ggplot(s, aes(x = ali_pos, y = cluster)) +
    geom_tile(
      aes(fill = FrustIC),
      linejoin = "round",
      width = 0.96,
      height = 0.96,
      color = 'grey20'
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label = consensus, size = SeqIC), fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_fill_manual(
      values = c("MIN" = "green", "NEU" = "grey", "MAX" = "red"),
      name = "FrustIC",
      limits = c("MIN", "NEU", "MAX"),
      na.value = "white"
    ) +
    scale_color_gradient2(
      mid = "grey",
      high = "black",
      midpoint = 0.5,
      limits = c(0, log2(20)),
      oob = scales::squish
    ) +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "top")
  
  # store plot
  list_plots[[k]] = g + theme(legend.position = "none")
}

# combine plots
combined_plot <- plot_grid(plotlist = list_plots, ncol = 1)

# extract legend to have a single legend
legend1 <- get_legend(
  g +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

# add legend to combined plot and legend using plot_grid()
p <- plot_grid(NULL, legend1, NULL, combined_plot, ncol = 1, rel_heights = c(.1, .2, .1, 1))

# add x and y labels common for all the plots and not by row of the plot
y.grob <- textGrob("Cluster", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13), rot = 90)
x.grob <- textGrob("Global MSA position", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13))

# save figure adjusting for height and width proportionally
height = length(cluster_names) + length(i) - 1
svg("consensus_MSFA_contributiont.svg", width = 15, height = height)
pp <- grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

dev.off()

######
# contribution plot with labels

x_positions <- c(6, 10, 26, 29, 39, 41, 53, 98, 123, 125, 126, 139)  
# Contribution plot
list_plots = list()
for (k in 1:length(i)) {
  # subset data referring to just a certain width
  s = subset(plot_df, ali_pos %in% i[k]:j[k] & !is.na(SeqIC))
  s$cluster <- as.factor(s$cluster)
  # filter positions
  current_x_positions <- x_positions[x_positions %in% i[k]:j[k]]
  
  # generate plot
  g = ggplot(s, aes(x = ali_pos, y = cluster)) +
    geom_tile(
      aes(fill = FrustIC),
      linejoin = "round",
      width = 0.96,
      height = 0.96,
      color = 'grey20'
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label = consensus, size = SeqIC), fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_fill_manual(
      values = c("MIN" = "green", "NEU" = "grey", "MAX" = "red"),
      name = "FrustIC",
      limits = c("MIN", "NEU", "MAX"),
      na.value = "white"
    ) +
    scale_color_gradient2(
      mid = "grey",
      high = "black",
      midpoint = 0.5,
      limits = c(0, log2(20)),
      oob = scales::squish
    ) +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "top")
  
  
  if (length(current_x_positions) > 0) {
    for (pos in current_x_positions) {
      #g <- g + annotate("text", x = pos, y = Inf, label = pos, vjust = 1.5, color = "blue", fontface = "bold", size = 4)
      g <- g + geom_tile(data = s[s$ali_pos == pos, ],
                         aes(x = ali_pos, y = cluster),
                         color = "blue", fill = NA, size = 1, width = 0.96, height = Inf) 
        #geom_text(data = subset(s, ali_pos == pos), aes(label = ali_pos),
                  #vjust = -0.5, hjust = 0.5, size = 4, color = "blue", fontface = "bold")
      g <- g +
        annotate("text", x = pos, y = Inf, label = pos, vjust = -0.6, color = "blue", fontface = "bold", size = 4)
    }
  }
  
  # store plot
  list_plots[[k]] = g + theme(legend.position = "none")
}

# combine plots
combined_plot <- plot_grid(plotlist = list_plots, ncol = 1)

for (pos in x_positions) {
  combined_plot <- combined_plot +
    annotate("text", x = pos, y = Inf, label = pos, vjust = 1.5, color = "blue", fontface = "bold", size = 4)
}

# extract legend to have a single legend
legend1 <- get_legend(
  g +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

# add legend to combined plot and legend using plot_grid()
p <- plot_grid(NULL, legend1, NULL, combined_plot, ncol = 1, rel_heights = c(.1, .2, .1, 1))

# add x and y labels common for all the plots and not by row of the plot
y.grob <- textGrob("Cluster", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13), rot = 90)
x.grob <- textGrob("Global MSA position", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13))

# save figure adjusting for height and width proportionally
height = length(cluster_names) + length(i) - 1
svg("consensus_MSFA_contribution_labelled.svg", width = 15, height = height)
pp <- grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

dev.off()

## labels aa name
x_positions <- c(6, 10, 26, 29, 39, 41, 53, 98, 123, 125, 126, 139)  
custom_titles <- c("Lys7", "Lys11", "Glu27", "Glu30", "Lys40", "Tyr42", "Gln54", "Lys99", "Ser124", "Asp126", "Lys127", "Tyr140") 


list_plots = list()
for (k in 1:length(i)) {
  # Subset data referring to just a certain width
  s = subset(plot_df, ali_pos %in% i[k]:j[k] & !is.na(SeqIC))
  s$cluster <- as.factor(s$cluster)
  current_x_positions <- x_positions[x_positions %in% i[k]:j[k]]

  g = ggplot(s, aes(x = ali_pos, y = cluster)) +
    geom_tile(
      aes(fill = ifelse(ICTot > 0.5, FrustIC, NA)),
      linejoin = "round",
      width = 0.96,
      height = 0.96,
      color = 'grey20'
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label = consensus, size = SeqIC), fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_fill_manual(
      values = c("MIN" = "green", "NEU" = "grey", "MAX" = "red"),
      name = "FrustIC",
      limits = c("MIN", "NEU", "MAX"),
      na.value = "white"
    ) +
    scale_color_gradient2(
      mid = "grey",
      high = "black",
      midpoint = 0.5,
      limits = c(0, log2(20)),
      oob = scales::squish
    ) +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "top")
  
  
  if (length(current_x_positions) > 0) {
    for (pos in current_x_positions) {
      pos_index <- which(x_positions == pos) 
      g <- g + geom_tile(data = s[s$ali_pos == pos, ],
                         aes(x = ali_pos, y = cluster),
                         color = "#3A5FCD", fill = NA, size = 1, width = 0.96, height = Inf) 
      g <- g +
        annotate("text", x = pos, y = Inf, label = custom_titles[pos_index], vjust = -0.6, hjust = 0.1, color = "#3A5FCD", fontface = "bold", size = 2.5, angle = 45)
    }
  }

  list_plots[[k]] = g + theme(legend.position = "none")
}

combined_plot <- plot_grid(plotlist = list_plots, ncol = 1)

for (pos_idx in seq_along(x_positions)) {
  pos <- x_positions[pos_idx]
  name <- custom_titles[pos_idx]
  
  combined_plot <- combined_plot +
    annotate("text", x = pos, y = Inf, label = name, vjust = 1.5, color = "blue", fontface = "bold", size = 4)
}


legend1 <- get_legend(
  g +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)


p <- plot_grid(NULL, legend1, NULL, combined_plot, ncol = 1, rel_heights = c(.1, .2, .1, 1))


y.grob <- textGrob("Cluster", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13), rot = 90)
x.grob <- textGrob("Global MSA position", 
                   gp = gpar(fontface = "bold", col = "grey20", fontsize = 13))


height = length(cluster_names) + length(i) - 1
svg("consensus_MSFA_contribution_labelled_aa.svg", width = 15, height = height)
pp <- grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

dev.off()
