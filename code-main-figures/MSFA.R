library(ggplot2)
library(dplyr)
suppressPackageStartupMessages(library("argparse"))  
parser <- ArgumentParser()
parser$add_argument("--dir", help="directory of the job")
parser$add_argument("--jobid", help="Job Id")
args <- parser$parse_args()


t <- read.table(paste(args$dir,'AuxFiles/DF_Colores', sep=""), sep = ',', header = TRUE)


df <- data.frame(
  PDBID = t$PDBID,
  pos = t$pos,
  AA = t$AA,
  FstState= t$FstState,
  FstI=t$FstI
)


colores <- c("NEU" = "gray", "MAX" = "red", "MIN" = "green", "-" = "white")


df_filtered <- df[!is.na(df$FstI), ]

filas <- list()


pdbids <- unique(df$PDBID)


for (pdbid in pdbids) {
  datos_pdbid <- df[df$PDBID == pdbid, ]
  fila <- datos_pdbid$FstI
  filas[[pdbid]] <- fila
}

df_resultado <- as.data.frame(do.call(rbind, filas))
rownames(df_resultado) <- pdbids

df_resultado[df_resultado == "N/A"] <- NA
df_resultado[is.na(df_resultado)] <- 0


num_clusters <- as.integer(length(pdbids)/3)  
set.seed(123)  
clusters <- kmeans(df_resultado, centers = num_clusters)

df_resultado$cluster <- clusters$cluster

df_ordenado <- df_resultado[order(df_resultado$cluster), ]

df_filtered$clustering <- df_resultado$cluster
cluster_labels <- as.character(clusters$cluster)
df_filtered$PDBID <- factor(df_filtered$PDBID, levels = rev(unique(df_filtered$PDBID)))
lfilas = length(filas)
laa=length(t$AA)
if ( lfilas < 20){
  laa=laa*lfilas
  lfilas = 20
}


p <- ggplot(df_filtered, aes(x = pos, y = reorder(PDBID, PDBID), fill = FstState)) +
  geom_tile(colour = NA, linewidth = 1) +
  geom_text(aes(label = AA), size = lfilas*0.06, color = "black", fontface = "bold" ) +
  scale_fill_manual(values = colores) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = lfilas * 0.2, face = "bold"),
        axis.text.y = element_text(margin = margin(r = -20)),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 8),  # Ajusta el tamaño del título de la leyenda aquí
        legend.text = element_text(size = 8),
        legend.margin = margin(t = 1, unit = "lines"),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
warnings()

p <- p + theme(legend.box.margin = margin(t = -lfilas*0.5, unit = "pt"))  
#h=lfilas*0.08
#w=laa*0.0025
ggsave(paste('MSFA_',args$jobid,'.svg', sep=""),  p, width = 9, height = 3, dpi = 300, bg = "white",limitsize = FALSE)

#### partial MSFA


posiciones_deseadas <- c(6, 10, 26, 29, 39, 41, 49, 53, 73, 98, 111, 123, 125, 126, 139, 140)
mis_labels <- c("Lys7", "Lys11", "Glu27", "Glu30", "Lys40", "Tyr42", "His50", "Gln54", 
                "Asp74", "Lys99", "His112", "Ser124", "Asp126", "Lys127", "Tyr140", "Arg141")


df_filtered_selected <- df_filtered %>%
  filter(pos %in% posiciones_deseadas)


df_filtered_selected$labels <- mis_labels[match(df_filtered_selected$pos, posiciones_deseadas)]


#unique_labels <- unique_labels[seq_along(unique(df_filtered_selected$pos))]

p <- ggplot(df_filtered_selected, aes(x = factor(pos), y = reorder(PDBID, PDBID), fill = FstState)) +
  

  geom_tile(colour = NA, linewidth = 0.8, width = 0.95, height = 0.95) +  
  
  geom_text(aes(label = AA), size = lfilas * 0.3, color = "black", fontface = "bold") +
  
 
  scale_x_discrete(breaks = posiciones_deseadas, labels = mis_labels) +
  
 
  scale_fill_manual(values = colores) +
  
  theme_minimal() +
  
  
  theme(
    panel.grid = element_blank(),  
    axis.title = element_blank(), 
    axis.text = element_text(size = lfilas * 0.6, face = "bold"),  
    axis.text.y = element_text(margin = margin(r = -1)),
    axis.text.x = element_text(color = "#3A5FCD", angle = 45, hjust = 1),
    axis.ticks = element_blank(),  
    legend.position = "top",  
    legend.title = element_text(size = lfilas * 0.7), 
    legend.text = element_text(size = lfilas * 0.7), 
    legend.margin = margin(t = 1, unit = "lines"),  
    legend.key.size = unit(0.8, "cm"),  
    legend.box.margin = margin(t = -lfilas * 0.6, unit = "pt"),  
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines")  
  )


ggsave(paste('MSFApartial_', args$jobid, '.svg', sep = ""), 
       p, 
       width = 10, 
       height = 10, 
       dpi = 300, 
       bg = "white", 
       limitsize = FALSE)
