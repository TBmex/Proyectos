## Baps Valencia (R, iqtree)
- El programa utilizado es **Fastbaps** el cual parte de un alineamiento previamente obtenido por el pipeline.
- El alineamiento usado en este caso sera un alineamiento 1177 de linaje cuatro de genomas de pacientes de TB el mismo se utilizo para crear un filogenia con **iqtree**

### Generando tabla de Genotipos con Fastbaps (R)
input: **run_alignment_no_resis_1177.fas**
Utilizamos el scrip llamado **"Proyecto"**
~~~
#Libraries
library(fastbaps)
library(ape)
#Nota el archivo original para los arboles es "multi" (Ver README.md)

#Loading data OK
fasta.file.name <- "run_alignment_no_resis_1177.fas"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)

#Detalle OK
sparse.data <- optimise_prior(sparse.data, type = "baps")

#Running fastbaps "baps.hc" es el archivo a obtener.
baps.hc <- fast_baps(sparse.data)

#Bayesian hierarchical partition
Fastbaps <- multi_res_baps(sparse.data, levels = 6)
~~~
Resultados
~~~
> head(Fastbaps)
  Isolates Level 1 Level 2 Level 3 Level 4 Level 5 Level 6
1      G01       1       3       7      24      67     115
2      G02       1       3       7      22      60     108
3      G03       1       2       4      12      30      49
4      G04       2       7      24      79     173     253
5    G1000       1       2       6      21      57     102
6    G1002       2       6      18      57     138     212
~~~

### Generacion de arbol filogenetico (ITOL)
input: **run_alignment_no_resis_1177.fas.treefile**
input: **Fastbapslvl2** (Archivo Itol)

> No explico el proceso solo es crear un archivo ITOL

output: https://itol.embl.de/tree/16111121936428041608308653

### Generacion de tablas (R)
Extraemos las columnas: ID y Level 2
input: **Fastbaps**
input: **Tabla: ID_Genotipo_Spain_Cluster**

~~~
> head(Fastbaps_lvl2)
     ID Genotipo
1 G1019        1
2 G1061        1
3  G110        1
4 G1113        1
5  G111        1
6 G1154        1
~~~

~~~
## Realizar tablas de valores de Rhierbaps & Baps
## Libreria
library(dplyr)

# Data cambiamos variables valor de X dependiendo lo que quiero obtener
x <- Fastbaps_lvl2
y <- Fastbaps_lvl2$Genotipo

# Obtenemos tabla de frecuencias
#cbind( Freq=table(x), Cumul=cumsum(table(x)), relative=prop.table(table(x)))
Freq <- cbind( Freq=table(y))

# Realizamos merge de Frecuencias con españoles
ID_Genotipo_Spain <- merge (x, Spain, by = "ID", all.x=TRUE)
Genotipo_Spain <- ID_Genotipo_Spain [c(2,3)]

# Obtenemos numero de españoles por genotipo
Freq_spain <- table(Genotipo_Spain)
Freq_spain <- as.data.frame.matrix (Freq_spain)

# Juntamos frecuencias
Frecuencias <- cbind.data.frame (Genotipo = c(1:16), Freq, Freq_spain, NoSpain = c(Freq-Freq_spain))
Frecuencias <- cbind(Frecuencias, Spain_porcen = c(((Frecuencias$Spain)*1) / (Frecuencias$Freq)))

# Frecuencia total de Genotipos en transmicion (cluster)
ID_Genotipo_Spain_Cluster <- merge (ID_Genotipo_Spain, Clusters_10snps, by = "ID", all.x=TRUE)

# Transmicion total
Genotipo_Cluster <- ID_Genotipo_Spain_Cluster [c(2,4)]
Genotipo_Cluster <- cbind(Genotipo_Cluster, In_cluster = c(Genotipo_Cluster$cluster10snps==TRUE))
Genotipo_incluster <- Genotipo_Cluster [c(1,3)]
Genotipo_incluster_allfreq <- cbind.data.frame( Freq=table(Genotipo_incluster))
names(Genotipo_incluster_allfreq) = c("Genotipo", "x", "Transmicion")

Spain_Genotipo_incluster <- filter(ID_Genotipo_Spain_Cluster, Spain == "Spain")
Spain_Genotipo_incluster <- Spain_Genotipo_incluster[!is.na(Spain_Genotipo_incluster$cluster10snps),]
Spain_Genotipo_incluster <- cbind(Spain_Genotipo_incluster, In_cluster = c(Spain_Genotipo_incluster$cluster10snps==TRUE))
Spain_Genotipo_incluster <- Spain_Genotipo_incluster [c(2,5)]
Spain_Genotipo_incluster <- cbind( Freq=table(Spain_Genotipo_incluster$Genotipo))

# Renombramos columnas y agregamos Españoles en cluster y en no cluster
Frecuencias <- cbind.data.frame(Frecuencias, Spain_incluster = Spain_Genotipo_incluster)
names(Frecuencias) = c("Genotipo", "N", "Spain", "No_spain", "Spain_x", "Spain_incluster" )
Frecuencias <- cbind(Frecuencias, Spain_nocluster = c(Frecuencias$Spain - Frecuencias$Spain_incluster))
Frecuencias <- cbind(Frecuencias, Spain_incluster_x = c(((Frecuencias$Spain_incluster)*1) / (Frecuencias$Spain)))

# Agregamos transmicion
Frecuencias <- cbind(Frecuencias, Total_incluster = Genotipo_incluster_allfreq$Transmicion)
Frecuencias <- cbind(Frecuencias, Spain_incluster_xx = c(((Frecuencias$Spain_incluster)*1) / (Frecuencias$Total_incluster)))

# Adecuamos columnas
names(Frecuencias) = c("Genotipo","N","Sp","No_sp","Sp_x","Sp_incluster","Sp_nocluster","Sp_incluster_x","N_incluster","N_sp_incluster_x")

#Agregamos foreing
Frecuencias <- cbind(Frecuencias, Foreing_incluster = c(((Frecuencias$N_incluster)) - (Frecuencias$Sp_incluster)))
Frecuencias <- cbind(Frecuencias, Foreing_incluster_x = c(((Frecuencias$Foreing_incluster)*1) / (Frecuencias$N_incluster)))

#Redondeamos
Frecuencias <- round (Frecuencias, 4)

# Creamso csvs
write.csv(Frecuencias, "Frecuencias")
write.csv(ID_Genotipo_Spain_Cluster, "ID_Genotipo_Spain_Cluster")

# Borramos archivos intermedios
rm (x, y, Freq, Freq_spain, ID_Genotipo_Spain, Genotipo_Spain, Spain_Genotipo_incluster, Genotipo_incluster, Genotipo_Cluster, Genotipo_incluster_allfreq)

~~~

output: **Tabla: Frecuencias**
### Generacion de graficos (R)
> No explico los parametros de ggplot2, pero tienen que ver con colores y posiciones de las etiquetas

~~~
#Plots
library(ggplot2)
library(ggpubr)

attach(Frecuencias)
names(Frecuencias)

ggplot(Frecuencias, aes(x=N_sp_incluster_x, y=Foreing_incluster_x, label = rownames(Frecuencias), group = Genotipo)) +
  geom_point(shape = 21, size = 4, color = "gray40", fill = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  ylab("percent of foreign") + xlab("spanish in cluster/clustered") +
  ggtitle("") +
  theme(plot.title=element_text(hjust = 0.5, size = 12)) +
  geom_text(vjust = c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-3,-2,-1,-1,-1), hjust = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), size = 3)
  theme(axis.text=element_text(size=11),axis.title.x = element_text (size=14, margin = margin(t = 8, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=14, margin = margin(t = 0, r = 5, b = 0, l = 0 )))
~~~

Resultado: Se grafica

x = Ratio de españoles en cluster
y = Ratio de extranjeros en cluster


![](assets/Proyecto_BAPS-103f89da)
