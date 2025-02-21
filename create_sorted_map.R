### Mappa Metodo 1
data_metodo1 <- read.csv("clustersCONMAPPA.csv")
data_full <-read.csv("OECD_data.csv")

data_full <- data_full%>%
  left_join(data_metodo1, by = c("CNT" = "country"))


cluster_summary <- data_full %>%
  group_by(cluster) %>%  # Raggruppiamo per Cluster
  summarise(Mean_Y_MATH1_rate = mean(Y_MATH1_rate, na.rm = TRUE))  # Calcoliamo la media

cl_metodo1 <- data_full$cluster
cl_temp <- rep(0,length(cl_metodo1))
## Ordiniamo
for (i in 1:dim(data_full)[1]){
  if(cl_metodo1[i]==1)
    cl_temp[i] = 7
  if(cl_metodo1[i]==2)
    cl_temp[i] = 10
  if(cl_metodo1[i]==3)
    cl_temp[i] = 12
  if(cl_metodo1[i]==4)
    cl_temp[i] = 5
  if(cl_metodo1[i]==5)
    cl_temp[i] = 13
  if(cl_metodo1[i]==6)
    cl_temp[i] = 8
  if(cl_metodo1[i]==7)
    cl_temp[i] = 6
  if(cl_metodo1[i]==8)
    cl_temp[i] = 3
  if(cl_metodo1[i]==9)
    cl_temp[i] = 4
  if(cl_metodo1[i]==10)
    cl_temp[i] = 1
  if(cl_metodo1[i]==12)
    cl_temp[i] = 2
  if(cl_metodo1[i]==11)
    cl_temp[i] = 11
  if(cl_metodo1[i]==13)
    cl_temp[i] = 9
}
data_full$cluster <- as.factor(cl_temp)
cl_metodo1 <- as.factor(cl_temp)

data2_metodo1 <- read.csv("OECD_data.csv")
clusters_metodo1 <- data.frame(CNT=data2_metodo1$CNT,clust=data_full$cluster)

library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(readr)
library(rnaturalearthdata)

clusters_metodo1 <- clusters_metodo1 %>%
  mutate(CNT = case_when(
    CNT == "United States" ~ "United States of America",
    CNT == "B-S-J-Z (China)" ~ "China",
    CNT == "Chinese Taipei" ~ "Taiwan",
    CNT == "Slovak Republic" ~ "Slovakia",
    CNT == "Brunei Darussalam" ~ "Brunei",
    TRUE ~ CNT
  ))
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world %>%
  left_join(clusters_metodo1, by = c("name" = "CNT")) %>%
  filter(!sov_a3 %in% c("ATA", "GRL", "ATF"))  # Rimuove Antartide, Groenlandia e Terre Australi Francesi

n_clusters <- length(unique(na.omit(world$clust)))
distinct_colors <- c(
  "#00441B", "#1B7837","#00FF00" ,"#A6DBA0", # 3 Verdi
  "#D0FF00", "#FFFF33", "#FFE119", "#FF8C00", # Giallo, Ocra, Arancione
  "#FF6347", "#FF69B4" , "#E30B5D" ,"#D7191C","#800020"   # Rosso Chiaro, Rosso, Bordeaux
)

ggplot() +
  geom_sf(data = world, aes(fill = factor(clust)), color = "whitesmoke") +
  scale_fill_manual(values = distinct_colors[1:n_clusters], name = "Cluster Index by Country", na.value = "grey80") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Cluster Map of Countries")

