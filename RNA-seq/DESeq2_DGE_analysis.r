library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(vsn)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(ggdendro)
library(apeglm)
library(ashr)

# PARTE 1: NORMALIZACION Y TRANSFORMACION DE READ COUNTS####


## PARTE 1.1: GENERACION DE LOS DATOS PARA DESeqDataSet####

# Obtener los counts
# Si se desea generar graficos con informacion de los cromosomas sexuales,
# se debe de omitir el paso de eliminar los genes del cromosoma U|V
# dlk_exonlvl.gene.txt es resultado de featurecounts
readcounts<-c("featurecounts_output/dlk_exonlvl.gene.txt")

readcounts <- read.table(readcounts, header=TRUE )

# Funcion para obtener la primer parte de un string separado por un punto
get_first_part <- function(x) {
  split_string <- strsplit(x, ".", fixed = TRUE)[[1]]
  return(split_string[1])
}

# Corregir los nombres de las columnas y filas
colnames(readcounts) <- sapply(colnames(readcounts), get_first_part)

rownames(readcounts) <- readcounts$Geneid

# Excluir columnas que no contienen counts
readcounts <- readcounts [ , -c(1:6) ]

# Revisar datos
str(readcounts)

head(readcounts)


readcounts <- readcounts [ , c(9:11,1:8) ]

# Revisar datos despues de reordenamiento de columnas
str(readcounts)

head(readcounts)

# Eliminar los genes de cromosomas sexuales
# NOTA: Omitir este paso si se desea generar figuras usando estos genes
readcounts<-readcounts[!grepl("^Mp(U|V)", rownames(readcounts)),]

# Crear los metadatos de las muestras
sample_info <- data.frame(condition = factor(gsub("_[0-9]+" , "", 
                                                  names(readcounts))),
                          row.names = names(readcounts))

# Generar el DESeqDataSet
# 1 + condition es lo mismo que  `~ condition`
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = sample_info,
                                   design = ~ 1 + condition )


# Revisar data.frame de la informacion de todos los samples 
head(colData(DESeq.ds))

# Revisar los conteos
# NOTA: Es lo mismo que usar counts(DESeq.ds)
head(assay(DESeq.ds, "counts"))

# Dado que al crear el DESeqDataSet no se le proporciono informacion de los
# genes, rowData deberia estar vacio
head(rowData(DESeq.ds))

# Revisar que valores se guardan en conteos
str(counts(DESeq.ds))


# Eliminar los genes que la suma de sus conteos entre todas las muestras sea
# igual a cero
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]

# Evaluar los diferentes tamanos de librerias
# NOTA: Se obtiene lo mismo al usar colSums(readcounts)
colSums(counts(DESeq.ds)) 


# Se debe de reordenar los datos, primero debe estar el control
# En este caso, el control ya estaba en primer lugar, pero falta indicar que
# la referencia es el WT. Si no se indica, R elige como referencia en orden
# alfabetico
str(colData(DESeq.ds)$condition)
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, ref = "WT")
str(colData(DESeq.ds)$condition)



## PARTE 1.2: NORMALIZACION DE LOS COUNTS PARA DIFERENCIAS EN PROFUNDIDAD DE SECUENCIACION ####

# Calcular el factor de tamano y agregarlo al dataset
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# Al evaluar colData, ahora contiene los sizeFactors
colData(DESeq.ds)

# Obtener los counts normalizados por los factores de tamano de las librerias
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE)



## PARTE 1.3: TRANSFORMACION DE LOS CONTEOS ####

### TRANSFORMACION LOG2 Y DISTRIBUCION DE CONTEOS####

# Transformar los counts normalizados a la escala log2 usando un pseudoconteo 
# de 1 (este valor se usa cuando se tiene un valor de 0 conteos, para que
# no de error al usar log2)
# NOTA: Tambian se puede usar normTransform
log.norm.counts <- log2(counts.sf_normalized + 1)

# Guardar los conteos transformados con log2 en un data.frame, para generar
# un boxplot con ggplot
lognormcount_df<-data.frame(lognormcount = as.vector(log.norm.counts), 
           gene_id = rep(rownames(log.norm.counts), 
                         times = dim(log.norm.counts)[2]),
           sample = rep(colnames(log.norm.counts), 
                        each = dim(log.norm.counts)[1]))

# Guardar los conteos solamente normalizados en un data.frame, para generar
# un boxplot con ggplot
normcount_df<-data.frame(normcount = as.vector(counts.sf_normalized), 
                            gene_id = rep(rownames(counts.sf_normalized), 
                                          times = dim(counts.sf_normalized)[2]),
                            sample = rep(colnames(counts.sf_normalized), 
                                         each = dim(counts.sf_normalized)[1]))

# Generar las figuras de los bloxplots para ambos data.frames
lognormcounts_boxplot<-ggplot(lognormcount_df, aes(x=sample,y=lognormcount)) +
  geom_boxplot() +
  labs(title = expression(bold("log"[2]*"-transformed read counts")),
       y = expression("log"[2]*"(normalized read counts)"),
       x = "Samples") +
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust=1),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))

normcounts_boxplot<-ggplot(normcount_df, aes(x=sample,y=normcount)) +
  geom_boxplot() +
  labs(title = "untransformed read counts",
       y = "normalized read counts",
       x = "Samples")+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust=1),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))

# Unir las graficas en una misma figura
normcounts_boxplots<-ggarrange(normcounts_boxplot, lognormcounts_boxplot,
          labels = c("A", "B"),
          font.label = list(size=18,face="bold"),
          ncol = 2, nrow = 1)

# Guardar figura de la comparacion de la distribucion de los conteos
ggsave(filename="counts_distribution_comparison.svg", 
       plot=normcounts_boxplots,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")

### TRANSFORMACION INCLUYENDO LA CONTRACCION DE LA VARIANZA (RLOG)####

# Graficar los conteos de dos replicas de una linea para comparar 
# la similitud de los valores de los conteos de cada gen entre las replicas
logncounts_plot<-ggplot(data.frame(gene=rownames(log.norm.counts),
                  rep1=log.norm.counts[,1],
                  rep2=log.norm.counts[,2]), aes(x=rep1, y=rep2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = expression(bold("log"[2]*"(normalized read counts)")),
       y = colnames(log.norm.counts)[2],
       x = colnames(log.norm.counts)[1])+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))


# Crear plot de la media y varianza de los conteos de los genes entre 
# todas las muestras

# NOTA: Observar que la varianza de los genes entre muestras con menor media
#       es mayor a los que tienen una mayor media de conteos
msd_plot_logncounts <- meanSdPlot(log.norm.counts, 
                       ranks = FALSE , # original scale,
                       plot = FALSE )
msd_plot_logncounts <- msd_plot_logncounts$gg + 
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16),
        strip.text = element_text(family="sans",size=16),
        legend.text = element_text(family="sans",size=16),
        legend.title =element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))+
  labs(title = expression(bold("log"[2]*"(normalized read counts)")),
       y="standard deviation")

# NOTA: Esto demuestra que los datos tienen un comportamiento heterocedastico, 
# y se debe corregir, dado que varios analisis y tests estadisticos asumen que 
# los datos son homocedasticos


# Obtener los conteos transformados con log2 y regularizados
# The transformation does not require that one has already estimated size 
# factors and dispersions.
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)


# Graficar los conteos de dos replicas de una linea para comparar 
# la similitud de los conteos de cada gen entre replicas despues de rlog
rlogncounts_plot<-ggplot(data.frame(gene=rownames(rlog.norm.counts),
                                   rep1=rlog.norm.counts[,1],
                                   rep2=rlog.norm.counts[,2]), 
                        aes(x=rep1, y=rep2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = "rlog(normalized read counts)",
       y = colnames(rlog.norm.counts)[2],
       x = colnames(rlog.norm.counts)[1])+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))

# Crear plot de la media y varianza de los conteos de los genes entre 
# todas las muestras
msd_plot_rlogncounts <- meanSdPlot(rlog.norm.counts,
                       ranks = FALSE,
                       plot = FALSE)

  
msd_plot_rlogncounts <- msd_plot_rlogncounts$gg +
  ggtitle("rlog(normalized read counts)") +
  ylab("standard deviation")+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16),
        strip.text = element_text(family="sans",size=16),
        legend.text = element_text(family="sans",size=16),
        legend.title =element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))



# Graficar los plots antes y despues de correccion de la dependencia de 
# la varianza sobre la media
variance_shrinkage_plots<-ggarrange(logncounts_plot, rlogncounts_plot, 
                                    msd_plot_logncounts, msd_plot_rlogncounts, 
                                    labels = c("A", "B", "C", "D"),
                                    font.label = list(size=18,face="bold"),
                                    common.legend =T,
                                    legend = "bottom",
                                    ncol = 2, nrow = 2)

# Guardar figura de la reduccion de la varianza lograda con rlog
ggsave(filename="variance_shrinkage_comparison.svg", 
       plot=variance_shrinkage_plots,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=40, units ="cm")


# Graficar los plots separados de la distribucion de los conteos y de la media
# y varianza
variance_shrinkage_plots<-ggarrange(logncounts_plot, rlogncounts_plot, 
                                    labels = c("A", "B"),
                                    font.label = list(size=18,face="bold"),
                                    common.legend =T,
                                    legend = "right",
                                    ncol = 2, nrow = 1)
variance_shrinkage_plots<-ggarrange(msd_plot_logncounts, msd_plot_rlogncounts, 
                                    labels = c("C", "D"),
                                    font.label = list(size=18,face="bold"),
                                    common.legend =T,
                                    legend = "right",
                                    ncol = 2, nrow = 1)

ggsave(filename="variance_shrinkage_comparison_read_distribution_plot.svg", 
       plot=variance_shrinkage_plots,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")

ggsave(filename="variance_shrinkage_comparison_msd_plot.svg", 
       plot=variance_shrinkage_plots,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")


## PARTE 1.4: ANALISIS EXPLORATORIOS DE LOS PATRONES DE EXPRESION ####

### GRAFICOS DE DISTRIBUCION DE CONTEOS DE CROMOSOMAS SEXUALES Y AUTOSOMALES ####

# NOTA: Para generar estas figuras, se requieren los genes de cromosomas 
# sexuales, por lo que no deben de ser eliminados si se desean volver a generar 
# estas figuras

# Generar graficos de distribucion de conteos de cromosomas sexuales y 
# autosomales.
sexual_genes<-readcounts[grepl("^Mp(U|V)", rownames(readcounts)),]
sexual_genes_df<-data.frame(count = as.vector(as.matrix(sexual_genes)),
                            gene_id = rep(rownames(sexual_genes),
                                          times = dim(sexual_genes)[2]),
                            sample = rep(colnames(sexual_genes),
                                         each = dim(sexual_genes)[1]))
sexual_genes<-log.norm.counts[grepl("^Mp(U|V)", rownames(log.norm.counts)),]
sexual_genes_df<-data.frame(count = as.vector(sexual_genes),
                            gene_id = rep(rownames(sexual_genes),
                                          times = dim(sexual_genes)[2]),
                            sample = rep(colnames(sexual_genes),
                                         each = dim(sexual_genes)[1]))

sexual_genes_df$chromosome <- ifelse(grepl("^Mp(V)", sexual_genes_df$gene_id),
                                     "ChrV", "ChrU")
autosomal_genes<-log.norm.counts[!grepl("^Mp(U|V)", rownames(log.norm.counts)),]
autosomal_genes_df<-data.frame(count = as.vector(autosomal_genes),
                               gene_id = rep(rownames(autosomal_genes),
                                             times = dim(autosomal_genes)[2]),
                               sample = rep(colnames(autosomal_genes),
                                            each = dim(autosomal_genes)[1]))

autosomal_genes_df$chromosome <- matrix(unlist(strsplit(autosomal_genes_df$gene_id, "g")),
                                        ncol=2,byrow = T)[,1]

# Generar graficos de distribucion de conteos de cromosomas sexuales y 
# autosomales.
autosomal_genes_dist_plot<-ggplot(autosomal_genes_df, aes(x=sample,y=count)) +
  geom_boxplot() +
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust = 1),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))+
  labs(title=expression("log"[2]*"(normalized read counts) of autosomal chromosomes"),
       y = expression("log"[2]*"(normalized read counts)"),
       x = "Samples") +
  facet_wrap(vars(chromosome),ncol = 3,nrow = 3)

sexual_genes_dist_plot<-ggplot(sexual_genes_df, aes(x=sample,y=count)) +
  geom_boxplot() +
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust=1),
        strip.text = element_text(family="sans",size=16),
        plot.title = element_text(family="sans",size=24,face="bold"))+
  labs(title=expression("log"[2]*"(normalized read counts) of sex chromosomes"),
       y = expression("log"[2]*"(normalized read counts)"),
       x = "Samples") +
  facet_grid(cols = vars(chromosome))

# Guardar graficos de distribucion de conteos de cromosomas sexuales y 
# autosomales.

# NOTA: Si se genera el grafico de distribucion de conteos para genes autosomales,
# sera necesario cambiar el directorio de with_sex_chromosomes, si es que
# se eliminaron previamente los genes sexuales
ggsave(filename="counts_distribution_autosomal_chrs.png",
       plot=autosomal_genes_dist_plot,
       path="figuras/with_sex_chromosomes",
       dpi=300, device = "png", width = 35,
       height=35, units ="cm")
ggsave(filename="counts_distribution_sex_chrs.png",
       plot=sexual_genes_dist_plot,
       path="figuras/with_sex_chromosomes",
       dpi=300, device = "png", width = 40,
       height=20, units ="cm")


### ANALISIS DE CALIDAD CON HEATMAP (MATRICES DE DISTANCIAS) ####
# Obtener las distancias de las muestras de la matriz transpuesta de los conteos
# Se usa la transpuesta porque la funcion dist() calcula las distancias de las 
# filas. Ademas, el metodo default para calcular las distancias es euclidean.
sampleDists <- dist(t(rlog.norm.counts))


# Generar la matriz con las distancias entre muestras, como input de pheatmap
sampleDistMatrix <- as.matrix(sampleDists)

# Generar heatmap de las distancias entre muestras usando metodo euclidiano
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
sampleDists_phm<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Calcular y generar figura de las distancias entre las muestras empleando la 
# correlacion de pearson
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson"))
distance.m_rlogmatrix <- as.matrix(distance.m_rlog)
distance_m_rlog_phm<-pheatmap(distance.m_rlogmatrix,
         clustering_distance_rows=distance.m_rlog,
         clustering_distance_cols=distance.m_rlog,
         col=colors, fontsize=20,
         treeheight_row = 100,
         treeheight_col = 100)

# Guardar figura con las distancias entre muestras
ggsave(filename="samples_distances.svg", 
       plot=distance_m_rlog_phm,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=40, units ="cm")

### CLUSTERING JERARQUICO ####

# Obtener las distancias de las muestras, mediante coeficiente de correlacion de
# Pearson
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson"))

# Generar el clustering jerarquico de las muestras y graficarlo
clustering_samples_plot<-plot(hclust(distance.m_rlog),
     labels = colnames(rlog.norm.counts),
     main = "rlog transformed read counts\ndistance: Pearson correlation")


# Obtener los datos de las distancias calculadas en la figura, para generar
# el dendograma empleando mejor ggplot
dhc <- as.dendrogram(hclust(distance.m_rlog))
data <- dendro_data(dhc, type = "rectangle")

clustering_samples_plot<-ggplot() +
  geom_segment(data = data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = data$labels, 
            aes(x = x, y = y, label = label), size=16/.pt,nudge_y=-0.0002)+
  guides(colour="none") +
  theme(
         axis.ticks.x=element_blank(),
         panel.grid.major=element_blank(), 
         panel.grid.minor=element_blank(),
         panel.background=element_rect(fill="white"),
         axis.text.x=element_blank(),
         axis.text.y=element_text(family="sans",size=12),
         axis.title.y=element_text(family="sans",size=20), 
         plot.title = element_text(family="sans",size=20,hjust = 0.5)) +
  labs(x="",y="Height", fill="Features", 
       title="rlog transformed read counts\ndistance: Pearson correlation")

clustering_samples_plot

# Guardar la figura del clustering de las muestras
ggsave(filename="samples_clustering.svg", 
       plot=clustering_samples_plot,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")

### PCA ####

# Generar PCA de las muestras
PCA_rlog<- plotPCA(DESeq.rlog)

PCA_rlog <- PCA_rlog + 
  theme_bw() +
  geom_hline(yintercept = 0 ,linewidth=1,colour="#D5D5D5") +
  geom_vline(xintercept = 0 ,linewidth=1, colour="#D5D5D5") +
  ggtitle("rlog transformed counts")+
  theme(axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16),
        text = element_text(family="sans",size=16))+
  labs(colour="Group")
PCA_rlog

# Guardar figura del PCA
ggsave(filename="samples_PCA.svg", 
       plot=PCA_rlog,
       path="figuras/without_sex_chromosomes/",
       dpi=1200, device = "svg", width = 20, 
       height=10, units ="cm")





# PARTE 2: ANALISIS DE GENES DIFERENCIALMENTE EXPRESADOS####

## PARTE 2.1: Contrasts de DESeq2, comparaciones para obtener los diferenciales ####

# Obtener la matriz del modelo del set de datos
mod_mat <- model.matrix(design(DESeq.ds), colData(DESeq.ds))
mod_mat 

# La variable condition, que es el factor de interes, es una dummy 
# variable, por tomar valores de 0 o 1.
# Este factor contiene niveles que representan los grupos de interes, las lineas
# de interes (dlk1,dlk25,OE13,WT)

# En esta matriz se representan los vectores numericos que contienen los pesos 
# de los coeficientes en cada grupo


# Obtener los vectores numericos que contengan los pesos de los coeficientes en
# cada grupo
# Obtener el vector numerico con los pesos para los coeficientes para el grupo
# de la linea silvestre
colMeans(mod_mat[DESeq.ds$condition == "WT", ])

# Obtener el vector numerico con los pesos para los coeficientes para el grupo
# de dlk1
colMeans(mod_mat[DESeq.ds$condition == "dlk1", ])
# Obtener el vector numerico con los pesos para los coeficientes para el grupo
#  de dlk25
colMeans(mod_mat[DESeq.ds$condition == "dlk25", ])
# Obtener el vector numerico con los pesos para los coeficientes para el grupo
# de OE13
colMeans(mod_mat[DESeq.ds$condition == "OE13", ])


# Para hacer las comparaciones de interes, se deben hacer las restas de los
# vectores que tengan los coeficientes de cada grupo de interes 

# Por ejemplo, el contraste de la diferencia entre dlk1 y WT se puede obtener
# mediante la resta de los vectores de grupos de dlk1 y WT
dlk1_vector<-colMeans(mod_mat[DESeq.ds$condition == "dlk1", ])
wt_vector<-colMeans(mod_mat[DESeq.ds$condition == "WT", ])
oe13_vector<-colMeans(mod_mat[DESeq.ds$condition == "OE13", ])

# Obtener el vector de contraste para permitir hacer el test para el efecto
# de la linea dlk1
dlk1_vector - wt_vector

# Al seguir esta estrategia (al guardar los vectores de cada grupo), se
# pueden hacer incluso otros contrastes mas complejos, como lineas mutantes 
# (dlk1 y dlk25) contra WT
dlk_vector <- colMeans(mod_mat[DESeq.ds$condition %in% c("dlk1","dlk25"), ])

dlk_vector - wt_vector

# En este caso, se asigna un peso de 0.4 al grupo de dlk1 y un peso de 0.6
# al grupo de dlk25. Esto se debe al numero de replicas de cada linea
dlk25_vector<-colMeans(mod_mat[DESeq.ds$condition == "dlk25", ])

# Otra forma de representar lo anterior es de la siguiente manera:
# (dlk25_vector*3+dlk1_vector*2)/5 - wt_vector

# NOTA: Otra cosa que se podria hacer es definir un nuevo grupo en la
# matriz de diseno, creando una nueva variable en column data, 
# para que aparte de la condition se tenga un nuevo factor que indique
# quien es mutante, sobreexpresora, y silvestre
# Entonces, el diseno del experimento se podria hacer ahora para la nueva 
# variable, la cual podria llamarse type: design(DESeq.ds) <- ~ 1 + type
# Sin embargo, en ese modelo, la dispersion de genes se estima en conjunto para 
# las muestras de dlk1 y dlk25 como si fueran replicas entre si, lo que puede 
# resultar en estimaciones infladas o desinfladas. 
# En cambio, el enfoque anterior estima el error dentro de cada uno de esos 
# grupos.

# Realizar comparacion entre las dos estrategias
DESeq.ds_de <- DESeq(DESeq.ds)
results1_dlkvsWT <- results(DESeq.ds_de, contrast = dlk_vector - wt_vector,
                            alpha = 0.05)


DESeq.ds$type <- factor(gsub("dlk[0-9]+" , "dlk", DESeq.ds$condition))
colData(DESeq.ds)

str(colData(DESeq.ds)$type)
colData(DESeq.ds)$type <- relevel(colData(DESeq.ds)$type, ref = "WT")
str(colData(DESeq.ds)$type)

design(DESeq.ds) <- ~ 1 + type
DESeq.ds_de <- DESeq(DESeq.ds)
resultsNames(DESeq.ds_de)
results2_dlkvsWT <- results(DESeq.ds_de, contrast = list("type_dlk_vs_WT"))


# Comparar los log2 fold-changes entre las dos estrategias
plot(results1_dlkvsWT$log2FoldChange, results2_dlkvsWT$log2FoldChange)
abline(0, 1, col = "brown", lwd = 2)

# Comparar los errores (lfcSE) entre las dos estrategias
plot(results1_dlkvsWT$lfcSE, results2_dlkvsWT$lfcSE)
abline(0, 1, col = "brown", lwd = 2)


# Al parecer, el caso mas general es dar un vector numerico. Otra
# forma es dar una lista con los nombres de las comparaciones a sumar, y los 
# nombres de las comparaciones a restar (los nombres los da resultsNames)

# Otra forma es dar un vector de caracteres con 3 elementos: el nombre del
# factor en la formula de diseno, el nombre del nivel para el numerador para
# el fold change, y el nombre del nivel para el denominador del fold change.
# Por ejemplo, para dlk1 vs WT, seria c("condition", "dlk1", "WT")



## PARTE 2.2: Analisis de expresion diferencial, y exploracion ####

# Obtener genes diferencialmente expresados
DESeq.ds_de <- DESeq(DESeq.ds)

# Devolver los efectos estimados, siendo las comparaciones realizadas
resultsNames(DESeq.ds_de) 

# Esta funcion devuelve 4 coeficientes, el primero es del intercepto (que lo
# tienen todas las muestras), mientras los otros tres son para el efecto de las 
# ldos ineas mutantes y la linea sobreexpresora.
# Esto quiere decir que los tres coeficientes restantes representan
# las diferencias entre cada linea contra la silvestre (el nivel de referencia)
# El intercepto corresponde entonces a la linea WT, el nivel de referencia

# Realizar un bloxpot de las distancias de cooks para evaluar si alguna de las 
# muestras es consistentemente mas alta que las otras
# Esto sirve para revisar en caso de que alguna muestra contega varios outliers,
# y grandes distancias de cooks indican presencia de outliers
boxplot(log10(assays(DESeq.ds_de)[["cooks"]]), range=0, las=2)


# Graficar los estimados de dispersion, donde se tienen los estimados
# de dispersion por gen
plotDispEsts(DESeq.ds_de, legend=FALSE)



# Asignar cutoff del "adjusted p-value" a 0.05
# Igualmente se puede asignar un threshold para log2foldchange, 
# usando el argumento lfcThreshold, pero esto es para hipotesis alternativas

# Devolver la tabla de resultados, para los efectos estimados (coeficientes)
# del modelo
# El argumento contrast permite especificar cual comparacion extraer
DGE.results_dlk1 <- results(DESeq.ds_de, alpha = 0.05,
                       contrast = c("condition", "dlk1", "WT"))


# Hacer revision de la tabla de resultados
summary(DGE.results_dlk1)
head(DGE.results_dlk1)
mcols(DGE.results_dlk1)$description
metadata(DGE.results_dlk1)$alpha

# Obtener informacion de todos los valores calculados
names(mcols(DESeq.ds_de))
mcols(mcols(DESeq.ds_de), use.names=TRUE)

metadata(DESeq.ds_de)[["version"]]



# Revisar genes con p-value ajustado menor a 0.05
table(DGE.results_dlk1$padj < 0.05)
rownames(subset(DGE.results_dlk1, padj < 0.05))

# Generar histogramas de p-values y p-values ajustados
hist(DGE.results_dlk1$pvalue, col = "grey", border = "white", 
     xlab = "", ylab = "",
     main = "frequencies of p-values")
hist(DGE.results_dlk1$padj, col = "grey", border = "white", 
     xlab = "", ylab = "",
     main = "frequencies of p-values")


# Generar grafico de MA
plotMA(DGE.results_dlk1, alpha = 0.05, main = "WT vs dlk1",
       ylim = c(-4 ,4))
abline(1,0)
abline(-1,0)

# Graficar la optimizacion para el filtro independentFiltering 
metadata(DGE.results_dlk1)$filterThreshold
plot(metadata(DGE.results_dlk1)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(DGE.results_dlk1)$lo.fit, col="red")
abline(v=metadata(DGE.results_dlk1)$filterTheta)


# Plot de conteos normalizados para genes de interes (lo que hace es normalizar
# por los size factors, y luego sumar 0.5 para permitir que si se desea se 
# aplique una escala logaritmica)
plotCounts(DESeq.ds, gene="Mp6g04830", intgroup="condition")
plotCounts(DESeq.ds, gene="Mp6g04830", intgroup="type")
plotCounts(DESeq.ds, gene=which.min(DGE.results_dlk1$padj), 
           intgroup="condition")


plotCounts(DESeq.ds, gene="Mp1g02580", intgroup="condition")
plotCounts(DESeq.ds, gene="MpVg00970", intgroup="condition")

# Obtener informacion de algun gen de interes
DGE.results_dlk1["Mp6g04830",]

DGE.results_dlk1["Mp1g02580",]



# COMPARACION DE USAR alpha y no usarlo

DGE.results1 <- results(DESeq.ds_de, alpha = 0.05,
                        contrast = c("condition", "dlk1", "WT"))

DGE.results2 <- results(DESeq.ds_de,
                        contrast = c("condition", "dlk1", "WT"))

DGE.results1<-as.data.frame(DGE.results1)
DGE.results2<-as.data.frame(DGE.results2)
DGE.results1<-na.omit(DGE.results1)
DGE.results2<-na.omit(DGE.results2)
DGE.results1<-DGE.results1[DGE.results1$padj < 0.05,]
DGE.results2<-DGE.results2[DGE.results2$padj < 0.05,]
DGE.results1<-DGE.results1[abs(DGE.results1$log2FoldChange) > 1,]
DGE.results2<-DGE.results2[abs(DGE.results2$log2FoldChange) > 1,]

nrow(DGE.results1)

nrow(DGE.results2)

# Numero de genes compartidos
shared_genes<-intersect(rownames(DGE.results1),rownames(DGE.results2))
length(shared_genes)
# Genes unicos de results1 (uno es un organic acid biosynthesis and exudation,
#otro es cap protein, y los demas desconocidos)
setdiff(rownames(DGE.results1),rownames(DGE.results2))
# Genes unicos de results2 (dos son desconocidos, otro es un transportador de 
# nitrito/nitrato)
setdiff(rownames(DGE.results2),rownames(DGE.results1))

# Cuando se usa alpha, se tiene que si se elige un valor menor a 0.1, aumentan
# los genes diferencialmente expresados, y si es mayor a 0.1, disminuyen.
# Esto se debe porque si se observa la grafica de los rechazos por el filtro
# de metadata(DGE.results)$filterNumRej, se ve que justo como se busca
# maximizar estos rechazos, por lo que aumentar el valor de alpha antes de
# la caida de la curva conlleva a aumentar los rechazos, pero mientras
# mas se acerque a un valor de 1 dejaría de hacer los rechazos.


# A su vez, hay genes unicos en cada set

# Evaluar si los p-values ajustados son iguales o menores
table(DGE.results1[shared_genes, "padj"]-DGE.results2[shared_genes,"padj"] > 0)

table(DGE.results1[shared_genes, "padj"] == DGE.results2[shared_genes, "padj"])


# Al parecer, practicamente todos los p-values ajustados son menores cuando se
# da el valor de alpha que sera el cutoff para elegir los diferencialmente
# expresados
# Por lo tanto, dado que el propio manual lo indica, siempre se dara a alpha el
# valor de cutoff de 0.05
# "If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha 
# should be set to that value."




## PARTE 2.3: Obtencion de genes diferencialmente expresados de todas las comparaciones ####

# GUARDAR TABLAS DE DIFERENCIALES
# Correccion: Se anadio el lfcthreshold para el valor deseado (valor 1), dado
# que el original es 0, y a partir de este es por el cual el test de wald estima
# los p-values
# It is also possible to provide thresholds for constructing Wald tests of 
# significance
# The test provides p values for the null hypothesis, the complement of the set 
# defined by the alternative.
# Esto es una prueba, y solo es para evaluar diferencias cuando se aplica este 
# threshold desde results, en vez de filtrar posteriormente usando
# abs(log2foldchange)>1

DGE.results_dlk1 <- results(DESeq.ds_de, alpha = 0.05, 
                            #lfcThreshold=1, altHypothesis="greaterAbs",
                            contrast = c("condition", "dlk1", "WT"),
                            #name="condition_dlk1_vs_WT" 
                            )

DGE.results_dlk25 <- results(DESeq.ds_de, alpha = 0.05, 
                             #lfcThreshold=1, altHypothesis="greaterAbs",
                             contrast = c("condition", "dlk25", "WT")
                             )

DGE.results_dlk <- results(DESeq.ds_de, alpha = 0.05, 
                           #lfcThreshold=1, altHypothesis="greaterAbs",
                           contrast = dlk_vector - wt_vector
                           )


DGE.results_OE13 <- results(DESeq.ds_de, alpha = 0.05, 
                            #lfcThreshold=1, altHypothesis="greaterAbs",
                            contrast = c("condition", "OE13", "WT")
                            )

# Guardar resultados originales antes de los cortes
write.csv(as.data.frame(DGE.results_dlk), 
          file="deseq2_fulldata_log2foldchange/rawresults_dlkvsWT.csv")
write.csv(as.data.frame(DGE.results_dlk1), 
          file="deseq2_fulldata_log2foldchange/rawresults_dlk1vsWT.csv")
write.csv(as.data.frame(DGE.results_dlk25), 
          file="deseq2_fulldata_log2foldchange/rawresults_dlk25vsWT.csv")
write.csv(as.data.frame(DGE.results_OE13), 
          file="deseq2_fulldata_log2foldchange/rawresults_OE13vsWT.csv")

# Visualizar como ahora solo los genes que pasan el cutoff del p-value de 0.05
# son aquellos con un abs(lfc)>1
plotMA(DGE.results_dlk, alpha = 0.05,
       ylim = c(-5 ,5)) +abline(1,0,lwd=2) +abline(-1,0,lwd=2)
plotMA(DGE.results_dlk1, alpha = 0.05,
       ylim = c(-5 ,5)) +abline(1,0,lwd=2) +abline(-1,0,lwd=2)
plotMA(DGE.results_dlk25, alpha = 0.05,
       ylim = c(-5 ,5)) +abline(1,0,lwd=2) +abline(-1,0,lwd=2)
plotMA(DGE.results_OE13, alpha = 0.05,
       ylim = c(-5 ,5)) +abline(1,0,lwd=2) +abline(-1,0,lwd=2)

# Graficos de MA 
# La funcion plotMA utiliza los valores de baseMean y log2FoldChange
# para generar los puntos en el plot

create_ma_plot<-function(dataset,ds_results,coef,path,title){
  resLFC <- lfcShrink(dataset,res=ds_results,coef=coef,
                      type="apeglm")
  resLFC<-na.omit(resLFC)
  #png(path,width = 10,height = 8,units = "in",res=300)
  svg(path,width = 10,height = 8)
  plotMA(resLFC, alpha = 0.05, main = title,
         ylim = c(-5 ,5)) +abline(1,0,lwd=2) +abline(-1,0,lwd=2)
  dev.off()
}

resultsNames(DESeq.ds_de)
create_ma_plot(DESeq.ds_de,DGE.results_dlk1,"condition_dlk1_vs_WT",
               "figuras/without_sex_chromosomes/maplot_dlk1_vs_wt_shrink.svg",
               "dlk1 vs WT")
create_ma_plot(DESeq.ds_de,DGE.results_dlk25,"condition_dlk25_vs_WT",
               "figuras/without_sex_chromosomes/maplot_dlk25_vs_wt_shrink.svg",
               "dlk25 vs WT")
create_ma_plot(DESeq.ds_de,DGE.results_OE13,"condition_OE13_vs_WT",
               "figuras/without_sex_chromosomes/maplot_OE13_vs_wt_shrink.svg",
               "OE13 vs WT")


# Crear ma plot usando ggplot
create_ma_plot<-function(resLFC,xlim,ylim,x_labels,file_name){
  resLFC<-data.frame(row.names = rownames(resLFC),
                     baseMean=resLFC@listData$baseMean,
                     log2FoldChange=resLFC@listData$log2FoldChange,
                     lfcSE=resLFC@listData$lfcSE,
                     pvalue=resLFC@listData$pvalue,
                     padj=resLFC@listData$padj
  )
  resLFC<-na.omit(resLFC)
  print(summary(resLFC$baseMean))
  resLFC$baseMean<-log10(resLFC$baseMean)
  print(summary(resLFC$baseMean))
  
  resLFC$is_DEG<-resLFC$padj<0.05 & abs(resLFC$log2FoldChange)>1
  resLFC$is_DEG<-ifelse(resLFC$is_DEG,1,0)
  resLFC$is_DEG<-ifelse(resLFC$is_DEG & resLFC$log2FoldChange > 1,2,resLFC$is_DEG)
  #print(table(resLFC$is_DEG))
  
  resLFC$flag <- resLFC$baseMean < xlim[1] | resLFC$baseMean > xlim[2] | abs(resLFC$log2FoldChange) > ylim
  resLFC$baseMean <- ifelse(resLFC$baseMean < xlim[1], -Inf, resLFC$baseMean)
  resLFC$baseMean <- ifelse(resLFC$baseMean > xlim[2], Inf, resLFC$baseMean)
  resLFC$log2FoldChange <- ifelse(abs(resLFC$log2FoldChange) > ylim, 
                                  sign(resLFC$log2FoldChange)*Inf, 
                                  resLFC$log2FoldChange)
  
  my_plot<-ggplot(resLFC,aes(x=baseMean,y=log2FoldChange,
                             colour=factor(is_DEG),shape=flag)) +
    geom_point(na.rm = T,size=2)+
    geom_hline(yintercept = c(-1, 1), col = "black", linetype = 'dashed')+
    
    theme_set(theme_classic(base_line_size=1.1) +
                theme(
                  plot.margin = margin(t=0.5,r=0.5,b=0.5,l=0.5,"cm"),
                  axis.title.y = element_text(margin = margin(0,0.5,0,0,"cm"), color = 'black'),
                  axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm"), color = 'black'),
                  axis.text.y=element_text(family="sans",size=16),
                  axis.text.x=element_text(family="sans",size=16),
                  plot.title = element_text(hjust = 0.5,face = "bold", 
                                            margin = margin(0,0.5,0.5,0,"cm"), size=18),
                  text = element_text(family="sans",size=16)
                )
              )+
    guides(shape="none")+
    
    coord_cartesian(ylim = c(-ylim, ylim), xlim = xlim,expand = F) +
    scale_x_continuous(breaks=xlim[1]:xlim[2],
                       labels=x_labels)+
    scale_colour_manual(values=c("#DDDDDD","#77AADD","#993F30"),
                        labels=c("Not Sig","Down", "Up"))+
    scale_shape_manual(values = c(16,1))+
                        #,labels=c("padj>=0.05","padj<0.05"))+
    labs(x="Mean of normalized counts",
         y=expression("log"[2]*" fold change"),
         #colour="Adjusted p-value\nthreshold",
         colour="",
         title=expression("MA plot of log"[2]*" fold changes"))
  
  my_plot
  ggsave(filename=file_name, 
         plot=my_plot,
         path="figuras/without_sex_chromosomes/",
         dpi=1200, device = "svg", 
         units ="cm",width = 30,height = 20)
  
}


create_ma_plot(DGE.results_dlk,c(0,5),10,
               expression(10^0,10^1,10^2,10^3,10^4,10^5),
               "maplot_dlk_vs_wt_ggplot.svg")
create_ma_plot(DGE.results_dlk1,c(0,5),10,
               expression(10^0,10^1,10^2,10^3,10^4,10^5),
               "maplot_dlk1_vs_wt_ggplot.svg")
create_ma_plot(DGE.results_dlk25,c(0,5),10,
               expression(10^0,10^1,10^2,10^3,10^4,10^5),
               "maplot_dlk25_vs_wt_ggplot.svg")
create_ma_plot(DGE.results_OE13,c(0,5),10,
               expression(10^0,10^1,10^2,10^3,10^4,10^5),
               "maplot_OE13_vs_wt_ggplot.svg")

# Generar funcion para eliminar genes con valores NA, y dejar solo
# aquellos con p-value ajustado menor a 0.05, y con 
# abs(log2FoldChange) mayor a 1
get_deg<-function(results){
  results<-as.data.frame(results)
  results<-na.omit(results)
  results<-results[results$padj < 0.05,]
  #s-value es solo necesario para ciertos casos de lfcshrink
  #results<-results[results$svalue < 0.05,]
  results<-results[abs(results$log2FoldChange) > 1,]
  results
}

# Realizar lfcshrink, para observar diferencias cuandos se usa y cuando no
# NOTA IMMPORTANTE: NO SE OBTIENEN LOS MISMOS RESULTADOS AL USAR
# COEF EN VEZ DE CONTRAST de lfcshrink
# Por otra parte, el estimador "normal" no funciona con vectores de contraste
# numericos, por lo que la comparacion dlk vs wt tampoco se puede obtener 
# para ese caso
DGE.results_dlk <- lfcShrink(DESeq.ds_de,res=DGE.results_dlk,
                             contrast = dlk_vector - wt_vector,
                              #lfcThreshold=1,
                              type="ashr")
DGE.results_dlk1 <- lfcShrink(DESeq.ds_de,res=DGE.results_dlk1,
                              #coef=2,
                              contrast = c("condition", "dlk1", "WT"),
                              #lfcThreshold=1,
                              type="ashr")
DGE.results_dlk25 <- lfcShrink(DESeq.ds_de,res=DGE.results_dlk25,
                               #coef=3,
                               contrast = c("condition", "dlk25", "WT"),
                              #lfcThreshold=1,
                              type="ashr")
DGE.results_OE13 <- lfcShrink(DESeq.ds_de,res=DGE.results_OE13,
                              #coef=4,
                              contrast = c("condition", "OE13", "WT"),
                              #lfcThreshold=1,
                              type="ashr")

DGE.results_dlk<-get_deg(DGE.results_dlk)
DGE.results_dlk1<-get_deg(DGE.results_dlk1)
DGE.results_dlk25<-get_deg(DGE.results_dlk25)
DGE.results_OE13<-get_deg(DGE.results_OE13)

write.csv(DGE.results_dlk, 
          file="dlk_vs_wt_results.csv")
write.csv(DGE.results_dlk1, 
          file="dlk1_vs_wt_results.csv")
write.csv(DGE.results_dlk25, 
          file="dlk25_vs_wt_results.csv")
write.csv(DGE.results_OE13, 
          file="OE13_vs_wt_results.csv")


# AGREGAR PRODUCTOS DE LOS GENES
gene_product <- c("MpTak_v6.1r2_all_gene_ids.txt")
gene_product <- read.table(gene_product, header=F)
rownames(gene_product)<-gene_product$V2

DGE.results_dlk$product <-gene_product[rownames(DGE.results_dlk),1]

write.csv(DGE.results_dlk,"dlk_vs_wt_products.csv",quote = F)

DGE.results_dlk1$product <-gene_product[rownames(DGE.results_dlk1),1]

write.csv(DGE.results_dlk1,"dlk1_vs_wt_products.csv",quote = F)

DGE.results_dlk25$product <-gene_product[rownames(DGE.results_dlk25),1]

write.csv(DGE.results_dlk25,"dlk25_vs_wt_products.csv",quote = F)

DGE.results_OE13$product <-gene_product[rownames(DGE.results_OE13),1]

write.csv(DGE.results_OE13,"OE13_vs_wt_products.csv",quote = F)


# AGREGAR CONVERSION DE IDS A DISTINTAS VERSIONES
# CARGAR TABLA DE CORRESPONDENCIA
correspondence_table<-c("Enrichment_analysis/GO_annotations/mptak_v6.1r1/mp_gene_correspondence_table_v6.1r2.tsv")

correspondence_table <- read.table(correspondence_table, sep = "\t", header = T, stringsAsFactors = FALSE)

# Anadir "-" a la columna de status para genes sin status
correspondence_table[correspondence_table$status=="", 8] <- "-"

# Agregar a la tabla de DEGs la conversion de IDs
dlk1_conversion<-DGE.results_dlk1
dlk1_conversion$id<-rownames(dlk1_conversion)
dlk1_conversion<-dlk1_conversion[,c(7,1:6)]
dlk1_conversion<-merge(dlk1_conversion, correspondence_table, by.x = "id", by.y = "MpTak_v6.1r2", all.x = TRUE)
# Evaluar si existen duplicaciones debido a que algun DEG es una fusion de dos 
# IDs de versiones anteriores
dim(dlk1_conversion)
ids<-dlk1_conversion$id
duplicated_ids <- duplicated(ids)
ids[duplicated_ids]
#"Mp1g24800" "Mp4g18590"
# Reordenar tabla, reemplazar valores NA por "-", y ordenar por lfc
dlk1_conversion<-dlk1_conversion[,c(9:14,1:7)]
dlk1_conversion[is.na(dlk1_conversion)] <- "-"
dlk1_conversion<-dlk1_conversion[order(dlk1_conversion$log2FoldChange), ]
write.csv(dlk1_conversion,quote =F, 
          file="dlk1_ids_conversion_with_locus_type.csv")
# Reemplazar los valores de locus_type por los productos de los genes
# encontrados en el gff
dlk1_conversion<-merge(dlk1_conversion, gene_product, by.x = "id", by.y = "V2", all.x = TRUE)
dlk1_conversion$locus_type<-dlk1_conversion$V1
dlk1_conversion$locus_type<-dlk1_conversion$V1
dlk1_conversion<-dlk1_conversion[,c(2:7,1,8:13)]
write.csv(dlk1_conversion,quote =F, 
          file="dlk1_ids_conversion_with_locus_type_with_GFF_products.csv")


# Agregar a la tabla de DEGs la conversion de IDs
dlk25_conversion<-DGE.results_dlk25
dlk25_conversion$id<-rownames(dlk25_conversion)
dlk25_conversion<-dlk25_conversion[,c(7,1:6)]
dlk25_conversion<-merge(dlk25_conversion, correspondence_table, by.x = "id", by.y = "MpTak_v6.1r2", all.x = TRUE)
# Evaluar si existen duplicaciones debido a que algun DEG es una fusion de dos 
# IDs de versiones anteriores
dim(dlk25_conversion)
ids<-dlk25_conversion$id
duplicated_ids <- duplicated(ids)
ids[duplicated_ids]
#"Mp1g24800" "Mp6g02170" "Mp6g02170"
# Reordenar tabla, reemplazar valores NA por "-", y ordenar por lfc
dlk25_conversion<-dlk25_conversion[,c(9:14,1:7)]
dlk25_conversion[is.na(dlk25_conversion)] <- "-"
dlk25_conversion<-dlk25_conversion[order(dlk25_conversion$log2FoldChange), ]
write.csv(dlk25_conversion,quote =F, 
          file="dlk25_ids_conversion_with_locus_type.csv")
# Reemplazar los valores de locus_type por los productos de los genes
# encontrados en el gff
dlk25_conversion<-merge(dlk25_conversion, gene_product, by.x = "id", by.y = "V2", all.x = TRUE)
dlk25_conversion$locus_type<-dlk25_conversion$V1
dlk25_conversion$locus_type<-dlk25_conversion$V1
dlk25_conversion<-dlk25_conversion[,c(2:7,1,8:13)]
write.csv(dlk25_conversion,quote =F, 
          file="dlk25_ids_conversion_with_locus_type_with_GFF_products.csv")


# Agregar a la tabla de DEGs la conversion de IDs
OE13_conversion<-DGE.results_OE13
OE13_conversion$id<-rownames(OE13_conversion)
OE13_conversion<-OE13_conversion[,c(7,1:6)]
OE13_conversion<-merge(OE13_conversion, correspondence_table, by.x = "id", by.y = "MpTak_v6.1r2", all.x = TRUE)
# Evaluar si existen duplicaciones debido a que algun DEG es una fusion de dos 
# IDs de versiones anteriores
dim(OE13_conversion)
ids<-OE13_conversion$id
duplicated_ids <- duplicated(ids)
ids[duplicated_ids]
#"Mp2g17310" "Mp5g20870"
# Reordenar tabla, reemplazar valores NA por "-", y ordenar por lfc
OE13_conversion<-OE13_conversion[,c(9:14,1:7)]
OE13_conversion[is.na(OE13_conversion)] <- "-"
OE13_conversion<-OE13_conversion[order(OE13_conversion$log2FoldChange), ]
write.csv(OE13_conversion,quote =F, 
          file="OE13_ids_conversion_with_locus_type.csv")
# Reemplazar los valores de locus_type por los productos de los genes
# encontrados en el gff
OE13_conversion<-merge(OE13_conversion, gene_product, by.x = "id", by.y = "V2", all.x = TRUE)
OE13_conversion$locus_type<-OE13_conversion$V1
OE13_conversion$locus_type<-OE13_conversion$V1
OE13_conversion<-OE13_conversion[,c(2:7,1,8:13)]
write.csv(OE13_conversion,quote =F, 
          file="OE13_ids_conversion_with_locus_type_with_GFF_products.csv")


# Graficos de volcano
# Crear volcano donde los DEGs que pasaron los cortes tengan color de
# acuerdo a si son up or down

create_volcano_plot<-function(results,x_limit,y_lim,title,file_name,
                              dpi_val,output_ft){
  results$diffexprs<-ifelse(results$log2FoldChange > 0, "UP", "DOWN")
  results$mlog10padj<- (-log10(results$padj))
  results$flag <- results$mlog10padj > y_lim | abs(results$log2FoldChange) > x_limit
  results$mlog10padj <- ifelse(results$mlog10padj > y_lim , 
                      Inf, results$mlog10padj)
  results$log2FoldChange <- ifelse(abs(results$log2FoldChange) > x_limit, 
                                sign(results$log2FoldChange)*Inf, 
                                results$log2FoldChange)
   
  
  my_plot<-ggplot(data = results, 
         aes(x = log2FoldChange, y = mlog10padj,colour=diffexprs,shape=flag)) +
    geom_vline(xintercept = c(-1, 1), col = "black",linetype = 'dashed')+
    geom_point(size=2)+
    theme_set(theme_classic(base_size = 20) +
                theme(
                  axis.title.y = element_text(face = "bold", 
                                              margin = margin(0,20,0,0), 
                                              size = rel(1.1), color = 'black'),
                  axis.title.x = element_text(hjust = 0.5, face = "bold", 
                                              margin = margin(20,0,0,0), 
                                              size = rel(1.1), color = 'black'),
                  plot.title = element_text(hjust = 0.5,face = "bold", 
                                            margin = margin(0,0,20,0), 
                                            size = rel(1.1))
                ))+

    scale_colour_manual(values=c("#77AADD","#993F30"),labels=c("Down", "Up"))+
    scale_shape_manual(values = c(16,1))+
    coord_cartesian(ylim = c(0, y_lim), xlim = c(-x_limit, x_limit),
                    expand = F) +
    scale_x_continuous(breaks = seq(-x_limit, x_limit, 2))+
    guides(shape="none")+
    labs(title = title,colour="",x = expression("log"[2]*" fold change"), 
         y = expression("-log"[10]*" adjusted p-value"))
  
  my_plot
  ggsave(filename=file_name, 
         plot=my_plot,
         path="figuras/without_sex_chromosomes/",
         dpi=dpi_val, device = output_ft, 
         units ="cm",width = 20,height = 20)
  
  
}


create_volcano_plot(DGE.results_dlk,10,300,
                    "Mutants (dlk1 + dlk25) vs WT","volcano_plot_dlk_vs_WT.svg",
                    #300,"png")
                    1200,"svg")
create_volcano_plot(DGE.results_dlk1,10,300,
                    "dlk1 vs WT","volcano_plot_dlk1_vs_WT.svg",
                    #300,"png")
                    1200,"svg")
create_volcano_plot(DGE.results_dlk25,10,300,
                    "dlk25 vs WT","volcano_plot_dlk25_vs_WT.svg",
                    #300,"png")
                    1200,"svg")
create_volcano_plot(DGE.results_OE13,10,300,
                    "OE13 vs WT","volcano_plot_OE13_vs_WT.svg",
                    #300,"png")
                    1200,"svg")


# GUARDAR LAS TABLAS ANADIENDO LOS IDS FALTANTES DE CADA COMPARACION
# Guardar todos los IDs de DEGs
degs <- c(rownames(DGE.results_dlk1),rownames(DGE.results_dlk25),
          rownames(DGE.results_OE13), rownames(DGE.results_dlk))
degs <- sort(degs)
degs <- unique(degs)

add_missing_ids<-function(degs, results){
  gene_ids<-setdiff(degs, rownames(results))
  
  num_empty_rows <- length(gene_ids)
  
  empty_rows <- data.frame(matrix(nrow = num_empty_rows, ncol = ncol(results)))
  colnames(empty_rows) <- colnames(results)
  rownames(empty_rows) <- gene_ids
  
  updated_results <- rbind(results, empty_rows)
  
  
  updated_results <- updated_results[match(degs, rownames(updated_results)), ]
  updated_results
}

results1<-add_missing_ids(degs,DGE.results_dlk)
results2<-add_missing_ids(degs,DGE.results_dlk1)
results3<-add_missing_ids(degs,DGE.results_dlk25)
results4<-add_missing_ids(degs,DGE.results_OE13)


write.csv(results1, 
          file="results1.csv")
write.csv(results2, 
          file="results2.csv")
write.csv(results3, 
          file="results3.csv")
write.csv(results4, 
          file="results4.csv")


### PARTE 2.4: Generar diagramas de Venn entre distintos los sets ####

# DIAGRAMAS DE VENN
features_palette<-c("#EEDD88","#99DDFF","#DDDDDD","#77AADD","#993F30",
                    "#FFAABB")


create_venn_diagram<-function(elements, cat_names, file_name,text_pos){
  venn_palette_col <- c("#EEBC88", "#D5D5D5","#993F30")
  venn_palette_fill <- c(alpha("#EEBC88",0.7), alpha('#D5D5D5',0.7), 
                         alpha("#993F30",0.7))
  my_plot <- venn.diagram(
    x = elements,
    disable.logging=T,
    category.names = cat_names,
    filename = NULL,
    #output = TRUE,
    lwd = 1,
    col= venn_palette_col[1:length(elements)],
    fill = venn_palette_fill[1:length(elements)],
    # Set names
    cat.cex = 3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(-25,180,25),
    cat.pos = text_pos,
    #cat.dist = c(0.055, 0.055),
    cat.fontfamily = "arial",
    
    
    fontfamily = "arial",
    cex = 4,
    margin = 0.1,
    
    imagetype = "svg",
    units = "cm"
  )
  
  size_vector<-c(20,30)
  ggsave(filename=file_name, 
         plot=my_plot,
         path="figuras/without_sex_chromosomes/",
         dpi=1200, device = "svg", 
         width = size_vector[(length(elements)-1)],
         height = size_vector[(length(elements)-1)],
         units ="cm")
  
}


# DEGs dlk1 VS DEGs dlk25
create_venn_diagram(elements = list(rownames(DGE.results_dlk1), 
                                    rownames(DGE.results_dlk25)),
                    cat_names = c("DEGs dlk1" , "DEGs dlk25"),
                    file_name ="venn_degsdlk1_vs_degsdlk25.svg",
                    c(-25,25))

# DEGs dlk1 VS Upregulated DEGs dlk25
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]), 
                                    rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]),
                                    rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,])),
                    cat_names = c("Downregulated\nDEGs dlk1", 
                                  "Upregulated\nDEGs dlk1",
                                  "Upregulated\nDEGs dlk25"),
                    file_name ="venn_degsdlk1_vs_degsdlk25_upregulated.svg",
                    c(-25,180,25))


# DEGs dlk1 VS Downregulated DEGs dlk25
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]), 
                                    rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]),
                                    rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,])),
                    cat_names = c("Downregulated\nDEGs dlk1", 
                                  "Upregulated\nDEGs dlk1",
                                  "Downregulated\nDEGs dlk25"),
                    file_name ="venn_degsdlk1_vs_degsdlk25_downregulated.svg",
                    c(-25,180,25))


# DEGs dlk1 VS DEGs dlk25 vs DEGs dlk (dlk1+dlk25)
create_venn_diagram(elements = list(rownames(DGE.results_dlk1), 
                                    rownames(DGE.results_dlk25),
                                    rownames(DGE.results_dlk)),
                    cat_names = c("DEGs dlk1" , "DEGs dlk25", "DEGs Mutants"),
                    file_name ="venn_degsdlk1_vs_degsdlk25_vs_degsdlk.svg",
                    c(-25,25,180))


# Downregulated DEGs dlk25 VS Upregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Downregulated\nDEGs dlk25" , "Upregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk25_downregulated_vs_degsOE13_upregulated.svg",
                    c(-25,25))

# Downregulated DEGs dlk1 VS Upregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Downregulated\nDEGs dlk1" , "Upregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk1_downregulated_vs_degsOE13_upregulated.svg",
                    c(-25,25))

# Downregulated DEGs dlk VS Upregulated DEGs OE13
# Nombres posibles para este set (dlk, dlk1/dlk25, Mutants)
create_venn_diagram(elements = list(rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange < 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Downregulated\nDEGs Mutants" , 
                                  "Upregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk_downregulated_vs_degsOE13_upregulated.svg",
                    c(-25,25))


# Upregulated DEGs dlk25 VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Upregulated\nDEGs dlk25" , "Downregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk25_upregulated_vs_degsOE13_downregulated.svg",
                    c(205,-205))

# Upregulated DEGs dlk1 VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Upregulated\nDEGs dlk1" , "Downregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk1_upregulated_vs_degsOE13_downregulated.svg",
                    c(205,-205))

# Upregulated DEGs dlk VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange > 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Upregulated\nDEGs Mutants" , 
                                  "Downregulated\nDEGs OE13"),
                    file_name ="venn_degsdlk_upregulated_vs_degsOE13_downregulated.svg",
                    c(205,-205))




# DEGs dlk25 VS Upregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]),
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Down in dlk25",
                                  "Up in dlk25",
                                  "Up in OE13"),
                    file_name ="venn_degsdlk25_vs_degsOE13_upregulated.svg",
                    c(-25,180,25))

# DEGs dlk1 VS Upregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Down in dlk1",
                                  "Up in dlk1",
                                  "Up in OE13"),
                    file_name ="venn_degsdlk1_vs_degsOE13_upregulated.svg",
                    c(-25,180,25))

# DEGs dlk VS Upregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange > 0,]),
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,])),
                    cat_names = c("Down in Mutants", 
                                  "Up in Mutants", 
                                  "Up in OE13"),
                    file_name ="venn_degsdlk_vs_degsOE13_upregulated.svg",
                    c(-25,180,25))


# DEGs dlk25 VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]),
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Down in dlk25",
                                  "Up in dlk25",
                                  "Down in OE13"),
                    file_name ="venn_degsdlk25_vs_degsOE13_downregulated.svg",
                    c(-25,180,25))

# DEGs dlk1 VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]), 
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Down in dlk1",
                                  "Up in dlk1",
                                  "Down in OE13"),
                    file_name ="venn_degsdlk1_vs_degsOE13_downregulated.svg",
                    c(-25,180,25))

# DEGs dlk VS Downregulated DEGs OE13
create_venn_diagram(elements = list(rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange < 0,]),
                                    rownames(DGE.results_dlk[DGE.results_dlk$log2FoldChange > 0,]),
                                    rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,])),
                    cat_names = c("Down in Mutants",
                                  "Up in Mutants",
                                  "Down in OE13"),
                    file_name ="venn_degsdlk_vs_degsOE13_downregulated.svg",
                    c(-25,180,25))

# Graficar genes diferenciales contradictorios de dlk25 y dlk1
plotCounts(DESeq.ds, gene="Mp3g09210", intgroup="condition")
plotCounts(DESeq.ds, gene="Mp5g21460", intgroup="condition")


## PARTE 2.5: Generar plots de conteos de falsos positivos por usar cromosomas sexuales ####

## FALSOS POSITIVOS POR USAR CROMOSOMAS SEXUALES ##
# Revisar a MpFGMYB
plotCounts(DESeq.ds, gene="Mp1g17210", intgroup="condition",transform = F)

# Graficar falsos positivos obtenidos por usar los genes de cromosomas sexuales
# Se toma como referencia a la comparacion de dlk25

# DEGs de dlk25 que son top4 hacia la alza y top4 hacia la baja
fp_table<-counts.sf_normalized[c("MpVg00760","MpVg00970","MpVg00440",
                                 "MpVg00980","MpUg00020","MpUg00350",
                                 "MpUg00130","MpUg00030"),]
fp_table_df<-data.frame(normcount = as.vector(fp_table),
                               gene_id = rep(rownames(fp_table),
                                             times = dim(fp_table)[2]),
                               sample = rep(colnames(fp_table),
                                            each = dim(fp_table)[1]))

fp_table_df$chr<-matrix(unlist(strsplit(fp_table_df$gene_id, "g")),
                        ncol=2,byrow = T)[,1]
fp_plot<-ggplot(fp_table_df, aes(x=sample,y=normcount, colour=gene_id)) +
  geom_point(alpha=0.8, size=2.5) +
  labs(title = "Top 8 dlk25 DEGs",
       color="Gene ID",
       y = "Normalized read counts",
       x = "Samples")+
  theme_bw()+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust = 1),
        strip.text = element_text(family="sans",size=16),
        legend.text = element_text(family="sans",size=16),
        legend.title =element_text(family="sans",size=16), 
        plot.title = element_text(family="sans",size=24,face="bold"))+
  facet_grid(cols = vars(chr))

ggsave(filename="top_dlk25_degs_with_sex_chrs.svg", 
       plot=fp_plot,
       path="figuras/with_sex_chromosomes",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")

# Generar plots de BPCU/MpUg00370 y BPCV/MpVg00350
bpcs_table<-counts.sf_normalized[c("MpUg00370","MpVg00350"),]
bpcs_table_df<-data.frame(normcount = as.vector(bpcs_table),
                        gene_id = rep(rownames(bpcs_table),
                                      times = dim(bpcs_table)[2]),
                        name = rep(c("BPCU", "BPCV"),
                                   times = dim(bpcs_table)[2]),
                        sample = rep(colnames(bpcs_table),
                                     each = dim(bpcs_table)[1]))


bpcs_plot<-ggplot(bpcs_table_df, aes(x=sample,y=normcount, colour=gene_id)) +
  geom_point(size=2.5) +
  labs(title = "",
       color="Gene ID",
       y = "Normalized read counts",
       x = "Samples")+
  theme_bw()+
  theme(axis.title = element_text(family="sans",size=20),
        axis.text.y=element_text(family="sans",size=16),
        axis.text.x=element_text(family="sans",size=16,angle=45,vjust = 1,
                                 hjust = 1),
        strip.text = element_text(family="sans",size=16),
        legend.text = element_text(family="sans",size=16),
        legend.title =element_text(family="sans",size=16), 
        plot.title = element_text(family="sans",size=24,face="bold"))+
  facet_grid(cols = vars(name))


ggsave(filename="bpcs_plot_with_sex_chrs.svg", 
       plot=bpcs_plot,
       path="figuras/with_sex_chromosomes",
       dpi=1200, device = "svg", width = 40, 
       height=20, units ="cm")


## PARTE 2.6: HEATMAPS DE LOS CONTEOS TRANSFORMADOS CON RLOG DE LOS DEGS ####

## HEATMAP ##
# Crear sets de datos, usando valores de rlog obteniendo solo los genes
# diferenciales de las intersecciones de las comparaciones

# Down Mutants vs Up OE
a<-intersect(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]),
          rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,]),
             a)

# Up Mutants vs Down OE
a<-intersect(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,]),
             a)

# Down Mutants vs Down OE
a<-intersect(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange < 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,]),
             a)

# Up Mutants vs Up OE
a<-intersect(rownames(DGE.results_dlk1[DGE.results_dlk1$log2FoldChange > 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,]),
             a)
# ALL
a<-intersect(rownames(DGE.results_dlk1),
             rownames(DGE.results_dlk25))

a<-intersect(rownames(DGE.results_OE13),
             a)

hm.mat_DGEgenes <- rlog.norm.counts[a, ]

# SIN dlk1
a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange < 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange < 0,]))

a<-intersect(rownames(DGE.results_OE13[DGE.results_OE13$log2FoldChange > 0,]),
             rownames(DGE.results_dlk25[DGE.results_dlk25$log2FoldChange > 0,]))

a<-intersect(rownames(DGE.results_OE13),
             rownames(DGE.results_dlk25))

hm.mat_DGEgenes <- rlog.norm.counts[a, -(4:5)]

# Graficar los conteos normalizados con rlog
aheatmap(hm.mat_DGEgenes, Rowv = NA , Colv = NA )
aheatmap(hm.mat_DGEgenes, Rowv = TRUE , Colv = TRUE ,distfun = "euclidean", 
         hclustfun = "average")
aheatmap(hm.mat_DGEgenes, Rowv = TRUE , Colv = TRUE ,distfun = "euclidean", 
         hclustfun = "average",scale = "row",labRow = NA,fontsize=16)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE , 
         Colv = TRUE,
         distfun = "pearson", 
         hclustfun = "complete", scale = "row",labRow = NA,fontsize=16)



# Comprobacion de como obtener distancias con correlacion de pearson, usando
# pheatmap y aheatmap
aheatmap(hm.mat_DGEgenes, 
         Rowv = list(as.dist(1 - cor(t(hm.mat_DGEgenes), method = "pearson")),
                     "complete",NA) , 
         Colv = list(as.dist(1 - cor(hm.mat_DGEgenes, method = "pearson")),
                     "complete",NA),
         revC = T,
         reorderfun = function(d, w) {return(d)},
         scale = "row",
         labRow = NA,fontsize=16)
pheatmap(hm.mat_DGEgenes,
         clustering_distance_cols=as.dist(1 - cor(hm.mat_DGEgenes, 
                                                  method = "pearson")),
         clustering_distance_rows=as.dist(1 - cor(t(hm.mat_DGEgenes), 
                                                  method = "pearson")),
         clustering_method = "complete",
         scale = "row",
         show_rownames = F)

pheatmap(hm.mat_DGEgenes,
         cluster_cols = hclust(as.dist(1 - cor(hm.mat_DGEgenes, 
                                               method = "pearson"))),
         cluster_rows = hclust(as.dist(1 - cor(t(hm.mat_DGEgenes), 
                                               method = "pearson"))),
         scale = "row",
         show_rownames = F)

# REVISAR A QUE SE DEBE DE QUE NO FUNCIONE USANDO "CORRELATION"
pheatmap(hm.mat_DGEgenes,
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method = "complete",
         scale = "row",
         show_rownames = F)

# Escalar manualmente los rlog values por row, generando z-scores
q<-t(apply(hm.mat_DGEgenes, 1, function(x){(x-mean(x))/sd(x)}))
# 2.4 (2.395803) con dlk1 y dlk25, 2.5 (2.503143) con dlk25
max(q)
# -2.13 (2.130143) con dlk1 y dlk25, -2.11 (-2.112721) con dlk25
min(q)

# COMPROBACION DE QUE REALMENTE AHEATMAP USA LA FUNCION DIST(), Y POR LO TANTO
# ESTA FUNCION TIENE UN ERROR AL NO USAR AS.DIST CUANDO SE PIDE EL METODO PEARSON
# reorderfun = function(d, w) {return(d)} ES LO MISMO A NO REORDENAR, SI SE DA
# VALOR DE NA EN COLV Y ROWV, ES LO MISMO A HACER ESTO
aheatmap(hm.mat_DGEgenes, Rowv = TRUE , 
         Colv = TRUE,
         distfun = "pearson", revC =T, reorderfun = function(d, w) {return(d)},
         hclustfun = "complete", scale = "row",labRow = NA,fontsize=16)
pheatmap(hm.mat_DGEgenes,
         clustering_distance_cols=dist(1 - cor(hm.mat_DGEgenes, 
                                               method = "pearson")),
         clustering_distance_rows=dist(1 - cor(t(hm.mat_DGEgenes), 
                                               method = "pearson")),
         clustering_method = "complete",
         scale = "row",
         show_rownames = F)
aheatmap(hm.mat_DGEgenes, 
         Rowv = list(dist(1 - cor(t(hm.mat_DGEgenes), method = "pearson")),
                     "complete",NA) , 
         Colv = list(dist(1 - cor(hm.mat_DGEgenes, method = "pearson")),
                     "complete",NA),
         revC =T,
         scale = "row",
         labRow = NA,fontsize=16)


# dist calcula las distancias, y por eso posee varios metodos
# as.dist espera una matriz de distancias, y la convierte en un objeto dist 
# (quedando una matriz triangular)
# aheatmap dice que calcula distancias usando dist(), y esto no es correcto

# POR OTRO LADO, PARECE QUE AHEATMAP Y PHEATMAP TIENEN COLORES DISTINTOS, O
# EL ESCALADO NO ES EL MISMO
# SE PUEDE COMPROBAR IMPRIMIENDO VALORES

#POR LO TANTO, SE USARA MEJOR PHEATMAP EN VEZ DE AHEATMAP


# Funcion para crear heatmap de degs
create_degs_heatmap<-function(filename, dataset){
  p<-pheatmap(dataset,
              #cluster_cols = hclust(as.dist(1 - cor(dataset, 
              #method = "pearson"))),
              cluster_cols = FALSE,
              cluster_rows = hclust(as.dist(1 - cor(t(dataset), 
                                                    method = "pearson"))),
              scale = "row",
              show_rownames = F,
              fontsize = 20,
              angle_col = 45,
              # debe tener un elemento mas que la rampa de colores
              breaks = seq(-2.5,2.5,length.out=101),
              #treeheight_col=150,
              treeheight_row=150)
  
  ggsave(filename=filename, 
         plot=p,
         path="figuras/without_sex_chromosomes/",
         dpi=1200, device = "svg", width = 40, 
         height=40, units ="cm")
  
}


create_degs_heatmap("heatmap_DEGs_UpOE_DownMutants_without_dlk1.svg",
                    hm.mat_DGEgenes)
create_degs_heatmap("heatmap_DEGs_DownOE_UpMutants_without_dlk1.svg",
                    hm.mat_DGEgenes)
create_degs_heatmap("heatmap_DEGs_DownOE_DownMutants_without_dlk1.svg",
                    hm.mat_DGEgenes)
create_degs_heatmap("heatmap_DEGs_UpOE_UpMutants_without_dlk1.svg",
                    hm.mat_DGEgenes)
create_degs_heatmap("heatmap_DEGs_OE_and_Mutants_without_dlk1.svg",
                    hm.mat_DGEgenes)

## 2.7 HEATMAP DE GENES DE FOTOSINTESIS ####
create_photosynthesis_heatmap<-function(filename, dataset,rnames,cnames,rlabels,
                                        acolors){
  p<-pheatmap(dataset,
              na_col = "purple",
              annotation_legend = F,
              color = colorRampPalette(c("#1373D1", "#DDDDDD", "#CC2408"))(100),
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              scale = "row",
              show_rownames = T,
              labels_row = rnames,
              labels_col = cnames,
              fontsize = 20,
              fontsize_col = 25,
              angle_col = 0,
              # debe tener un elemento mas que la rampa de colores
              breaks = seq(-2.5,2.5,length.out=101)
              #treeheight_col=150,
              #treeheight_row=150
              )
  
  ggsave(filename=filename, 
         plot=p,
         path="figuras/without_sex_chromosomes/",
         dpi=1200, device = "svg", width = 20, 
         height=25, units ="cm")
  
}

# PsbA y LHCB7 no pasan cortes de pvalue en ambas lineas OE13 y dlk25
# PsbaA Mp5g10950
# LHCB7 "Mp8g12010"
genes_of_interest<-c("Mp5g20680","Mp4g10720","Mp5g04200", 
                     "Mp2g15420","Mp8g13180","Mp4g00930","Mp3g07840",
                     "Mp7g05530","Mp7g06760",
                     "Mp7g05890","Mp7g06740","Mp7g09340",
                     "Mp7g05880","Mp7g05980","Mp7g06710","Mp7g06720",
                     "Mp7g06730","Mp7g06750","Mp7g06770","Mp7g06780",
                     "Mp7g06790","Mp7g09180","Mp1g29620","Mp2g13460",
                     "Mp7g08940","Mp1g16850","Mp4g10900","Mp6g01650",
                     "Mp4g09890","Mp1g08320","Mp4g21140","Mp7g02790","Mp4g07280",
                     "Mp4g15720","Mp8g11220",
                     "Mp3g03430","Mp3g17660","Mp2g20660")
genes_names<-c("VIPP1","NDHH","PsaD",
               "LHCA1","LHCA2","LHCA3","LHCA4",
               "LHCB1/2-like","LHCB1/2-like","LHCB1/2-like","LHCB1/2-like",
               "LHCB1/2-like","LHCB1/2-like","LHCB1/2-like","LHCB1/2-like",
               "LHCB1/2-like","LHCB1/2-like","LHCB1/2-like","LHCB1/2-like",
               "LHCB1/2-like","LHCB1/2-like","LHCB1/2-like","LHCB1/2-like",
               "LHCB1/2-like","LHCB3","LHCB4","LHCB5","LHCB6",
               "MpRBCS","HEMA1","CAO","HCAR", "NYC1", "NOL", "BCM2",
               "MpPR1h","MpWRKY7","MpICS")

oe_degs_ft<-DGE.results_OE13[genes_of_interest,2]
dlk25_degs_ft<-DGE.results_dlk25[genes_of_interest,2]
lth<-length(genes_of_interest)
  
annotation_matrix<-data.frame(rep("",lth),ifelse(is.na(dlk25_degs_ft), "", "*"),rep("",lth),
           rep("",lth),rep("",lth),rep("",lth),
           rep("",lth),ifelse(is.na(oe_degs_ft), "", "*"),rep("",lth))

my_labels <-as.factor(c("Thylakoid membrane biogenesis","NDH Complex","PSI",
                      rep("LHCA",4),rep("LHCB",21),"Rubisco",
                      rep("Chlorophyll biogenesis",3),
                      rep("Chlorophyll b degradation",2),
                      "Chlorophyll metabolism regulation",
                      rep("SA related genes",3)))
my_labels<-data.frame(Classification = my_labels)
rownames(my_labels)<-genes_of_interest
labels_vector<-levels(my_labels$Classification)
labels_colors<-brewer.pal(length(labels_vector), "Paired")
names(labels_colors)<-labels_vector
labels_colors<-list(Classification=labels_colors)

create_photosynthesis_heatmap("heatmap_fotosistemas_genes_fotosinteticos_rlog_LHCB_corrected_with_LHCAs_with_SA_genes.svg",
                              rlog.norm.counts[genes_of_interest, c(6:8,1:3,9:11)],
                              genes_names,
                              c("Mpdlk-25","Tak","OE13"),
                              my_labels,labels_colors)


## PARTE 2.8: HEATMAPS DE LOS GENES DE SALISILICO ####

genes_of_interest<-c("Mp5g16350","Mp2g20660","Mp1g19590","Mp8g10060", 
                     "Mp3g01570","Mp1g10150", "Mp5g19260", "Mp6g07600",
                     "Mp1g04370"
)
genes_names<-c("nahG","sid2/MpICS","MpTCP2","EFR", "EPS1","PAL4", "AIM1",
               "PBS3","EDS5")

genes_of_interest<-c("Mp2g11950","Mp3g17660","Mp7g06550","Mp2g14740",
                     "Mp3g03430","Mp7g13960","Mp3g14130","Mp5g12510"
)
genes_names<-c("MpPR1e","MpWRKY7","MpWRKY10","MpGH17.11",
               "MpPR1h","MpPR1w","MpPR1j","MpPR1m")

my_labels <-as.factor(rep("SA",8))
my_labels<-data.frame(Classification = my_labels)
rownames(my_labels)<-genes_of_interest
labels_vector<-levels(my_labels$Classification)
labels_colors<-brewer.pal(length(labels_vector), "Paired")
names(labels_colors)<-labels_vector
labels_colors<-list(Classification=labels_colors)

create_photosynthesis_heatmap("heatmap_SA_genes_rlog_response.svg",
                              rlog.norm.counts[genes_of_interest, c(6:8,1:3,9:11)],
                              genes_names,
                              c("Mpdlk-25","Tak","OE13"),
                              my_labels,labels_colors)


## PARTE 2.9: HEATMAPS DE LOS TFs en los DEGS ####

tfs_files<-c("dlk_vs_wt_nomenclature_transcription_factor.txt",
             "dlk1_vs_wt_nomenclature_transcription_factor.txt",
             "dlk25_vs_wt_nomenclature_transcription_factor.txt",
             "OE13_vs_wt_nomenclature_transcription_factor.txt"
             )

samples_columns<-list(4:8,4:5,6:8,9:11)
tf_contrast<-c("dlk_vs_wt","dlk1_vs_wt","dlk25_vs_wt","OE13_vs_wt")
for(index in 1:length(tfs_files)){
  # Generar lista de los TFs con los archivos de nomenclatura
  tfs_file<-tfs_files[index]
  results_tfs <- paste("TFs/without_sex_chromosomes/",
                     tfs_file,
                     sep="")
  results_tfs <- scan(results_tfs, what="", sep="\n")
  results_tfs<-strsplit(results_tfs, "\t")
  
  results_tfs <- lapply(results_tfs, function(x) {
    x[c(1:4,6)]
  })
  tfs_df<-results_tfs
  
  # Acortar los nombres de las familias, si tienen mas de dos comas
  tfs_df<-lapply(tfs_df, function(x) {
    fam_name<-sub("^([^,]*, [^,]*), .*$", "\\1", x[4])
    return(c(x[1:3],fam_name,x[5]))
  })
  
  # Obtener todos los nombres de cada gen
  tfs_df<-lapply(tfs_df, function(x) {
    # Asumir que siempre las primeras tres columnas contendran nombres
    gene_names<-x[1:3]
    # Separar y corregir nombres si en una columna hay dos o mas nombres anotados
    gene_names<-unlist(strsplit(gene_names,","))
    gene_names<-sub(" ", "", gene_names)
    # Eliminar nombres que esten vacios (no hay caracteres)
    gene_names<-gene_names[nchar(gene_names)>0]
    # Obtener nombres que solo inicien con Mp
    gene_names<-gene_names[grepl("^Mp", gene_names)]
    return(c(gene_names,x[4:5]))
  })
  
  
  # Obtener los IDs de cada gen y asignarlos a los vectores correspondientes
  names(tfs_df) <- sapply(tfs_df, function(x) {
    x<-sub("^.*; (.*)$", "\\1", x[length(x)])
    x<-sub("^(Mp.*)\\..*$", "\\1", x)
    x
  })
  
  # Eliminar nombres duplicados por el mismo motivo de que cada fila tenia un 
  # nombre diferente para el mismo gen, e igualmente eliminar nombres de la
  # familia del TF duplicadas
  tfs_df<-lapply(split(tfs_df, names(tfs_df)), function(sublist) {
    new_names<-unlist(lapply(sublist,function(x){x[1:(length(x)-2)]}))
    new_fams<-unlist(lapply(sublist,function(x){x[(length(x)-1)]}))
    tf_id<-names(sublist)[1]
    new_names<-unique(new_names)
    new_fams<-unique(new_fams)
  
    if(length(new_names)>1){
      #Eliminar nombres duplicados pero que esten en minusculas
      # Se asume que la columna 2 es la que los contiene
      if(toupper(new_names[1])==toupper(new_names[2])){
        new_names<-new_names[-2] 
      }
      
      new_names<-paste(new_names,collapse="/")
    }
    if(length(new_fams)>1){
      new_fams<-paste(new_fams,collapse="/")
    }
    
    c(new_names,new_fams,tf_id)
  })
  
  # Generar dataframe
  tfs_df<-data.frame(matrix(unlist(tfs_df), ncol=3,byrow = T))
  
  # Eliminar los ids de versiones anteriores a la 6.1, y corregir si estos son
  # isoformas
  tfs_df$X3<-sapply(tfs_df$X3,function(x){
    x<-sub("^.*; (.*)$", "\\1", x)
    x<-sub("^(Mp.*)\\..*$", "\\1", x)
  })
  # Asignar a cada fila los ids de los genes
  rownames(tfs_df)<-tfs_df$X3
  # Reemplazar transcription factor por family
  tfs_df$X2<-sapply(tfs_df$X2,function(x){
    x<-sub("transcription factor", "family", x)
  })
  
  # Obtener los rlog counts para las muestras, y asignar nombres a las filas
  rlogcounts_tfs<-rlog.norm.counts[tfs_df$X3,c(1:3,samples_columns[[index]])]
  rownames(rlogcounts_tfs)<-tfs_df[tfs_df$X3,1]
  
  # Generar heatmap
  tfs_heatmap<-pheatmap(rlogcounts_tfs,
                        cluster_cols = FALSE,
                        cluster_rows = hclust(as.dist(1 - cor(t(rlogcounts_tfs), 
                                                              method = "pearson"))),
                        scale = "row",
                        show_rownames = T,
                        fontsize = 20,
                        angle_col = 45,
                        # debe tener un elemento mas que la rampa de colores
                        breaks = seq(-2.5,2.5,length.out=101),
                        #treeheight_col=150,
                        treeheight_row=150)
  ggsave(filename=paste("tfs_heatmap_",tf_contrast[index],".svg",sep=""), 
         plot=tfs_heatmap,
         path="figuras/without_sex_chromosomes/",
         dpi=1200, device = "svg", width = 40, 
         height=40, units ="cm")
  
}


tfs_list<-list()
for(index in 1:length(tfs_files)){
  # Generar lista de los TFs con los archivos de nomenclatura
  tfs_file<-tfs_files[index]
  results_tfs <- paste("TFs/without_sex_chromosomes/",
                       tfs_file,
                       sep="")
  results_tfs <- scan(results_tfs, what="", sep="\n")
  results_tfs<-strsplit(results_tfs, "\t")
  
  results_tfs <- lapply(results_tfs, function(x) {
    x[c(1:4,6)]
  })
  tfs_df<-results_tfs
  
  # Acortar los nombres de las familias, si tienen mas de dos comas
  tfs_df<-lapply(tfs_df, function(x) {
    fam_name<-sub("^([^,]*, [^,]*), .*$", "\\1", x[4])
    return(c(x[1:3],fam_name,x[5]))
  })
  
  # Obtener todos los nombres de cada gen
  tfs_df<-lapply(tfs_df, function(x) {
    # Asumir que siempre las primeras tres columnas contendran nombres
    gene_names<-x[1:3]
    # Separar y corregir nombres si en una columna hay dos o mas nombres anotados
    gene_names<-unlist(strsplit(gene_names,","))
    gene_names<-sub(" ", "", gene_names)
    # Eliminar nombres que esten vacios (no hay caracteres)
    gene_names<-gene_names[nchar(gene_names)>0]
    # Obtener nombres que solo inicien con Mp
    gene_names<-gene_names[grepl("^Mp", gene_names)]
    return(c(gene_names,x[4:5]))
  })
  
  
  # Obtener los IDs de cada gen y asignarlos a los vectores correspondientes
  names(tfs_df) <- sapply(tfs_df, function(x) {
    x<-sub("^.*; (.*)$", "\\1", x[length(x)])
    x<-sub("^(Mp.*)\\..*$", "\\1", x)
    x
  })
  
  # Eliminar nombres duplicados por el mismo motivo de que cada fila tenia un 
  # nombre diferente para el mismo gen, e igualmente eliminar nombres de la
  # familia del TF duplicadas
  tfs_df<-lapply(split(tfs_df, names(tfs_df)), function(sublist) {
    new_names<-unlist(lapply(sublist,function(x){x[1:(length(x)-2)]}))
    new_fams<-unlist(lapply(sublist,function(x){x[(length(x)-1)]}))
    tf_id<-names(sublist)[1]
    new_names<-unique(new_names)
    new_fams<-unique(new_fams)
    
    if(length(new_names)>1){
      #Eliminar nombres duplicados pero que esten en minusculas
      # Se asume que la columna 2 es la que los contiene
      if(toupper(new_names[1])==toupper(new_names[2])){
        new_names<-new_names[-2] 
      }
      
      new_names<-paste(new_names,collapse="/")
    }
    if(length(new_fams)>1){
      new_fams<-paste(new_fams,collapse="/")
    }
    
    c(new_names,new_fams,tf_id)
  })
  
  # Generar dataframe
  tfs_df<-data.frame(matrix(unlist(tfs_df), ncol=3,byrow = T))
  
  # Eliminar los ids de versiones anteriores a la 6.1, y corregir si estos son
  # isoformas
  tfs_df$X3<-sapply(tfs_df$X3,function(x){
    x<-sub("^.*; (.*)$", "\\1", x)
    x<-sub("^(Mp.*)\\..*$", "\\1", x)
  })
  # Asignar a cada fila los ids de los genes
  rownames(tfs_df)<-tfs_df$X3
  # Reemplazar transcription factor por family
  tfs_df$X2<-sapply(tfs_df$X2,function(x){
    x<-sub("transcription factor", "family", x)
  })
  
  tfs_list[[index]]<-tfs_df
  
}
dlk25_families<-tfs_list[[3]]
dlk25_families$X2[1]<-c("family, C2H2-ZnF")
OE13_families<-tfs_list[[4]]
OE13_families$X2[1]<-c("family, C2H2-ZnF")
OE13_families$X2[2]<-c("family, TCP")
OE13_families$X2[4]<-c("family, C2H2-ZnF")
OE13_families$X2[20]<-c("family, ZIP")
