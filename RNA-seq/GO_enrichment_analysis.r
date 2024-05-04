library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(scales)
library(ggnewscale) #Permite agregar una nueva escala de colores en ggplot
library(httr)
library(jsonlite)
library(xml2)
library(extrafont)
font_import()

# PREPARACION DE LOS TERMS ####
# Obtener version del go ontology
requestURL <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/about"
ontology_info <- GET(requestURL, accept("application/json"))
stop_for_status(ontology_info)
ontology_info<-fromJSON(toJSON(content(ontology_info)))

ontology_info$go
# Antes de correccion y despues de correcion, son el mismo
# $version
# [1] "http://purl.obolibrary.org/obo/go/releases/2023-08-12/extensions/go-plus.owl"
#
# $timestamp
# [1] "2023-08-12"


# Obtener la informacion de los terminos
# Implementation Notes
# If possible, response fields include: 
# id, isObsolete, name, definition, ancestors, synonyms, 
# comment, aspect (for GO) and usage

get_go_terms_info <- function(IDs){
  terms_info <-lapply(IDs, function(ID){
    print(ID)
    requestURL <- paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",
                        ID, sep="")
    r <- GET(requestURL, accept("application/json"))
    stop_for_status(r)
    
    r<-fromJSON(toJSON(content(r)))
    #return(colnames(r$results))
    #return(unlist(r$numberOfHits))
    return(c(ID,
             unlist(r$results$id),
             unlist(r$results$isObsolete),
             unlist(r$results$name),
             unlist(r$results$aspect)))
  })
  return(terms_info)
}

# Evaluar nuevos terminos corregidos para GO terms secundarios
# y obsoletos (dia 23/08/2023)
go_terms_after_correction <- paste("Enrichment_analysis/GO_annotations/",
                                   "mptak_v6.1r1/",
                                   "MpTak_v6.1_GO_annotation_1line_",
                                   "first_isoform_non_sexchr_genes_",
                                   "rev2_correction_sec_corrt_",
                                   "obs_corrt_terms.txt",sep="")

#interproscan
#go_terms_after_correction <- paste("Enrichment_analysis/GO_annotations/",
#                                   "interproscan/",
#                                   "Mptakv6.1r2_go_annotation_1line_",
#                                   "ref_iso_non_sexchr_genes_terms.txt",sep="")

# GO terms antes de correccion
#go_terms_before_correction <- c("Enrichment_analysis/MpTak_v6.1_GO_annotation_1line_first_isoform_non_sexchr_genes_rev2_correction_terms.txt")
go_terms_after_correction <- scan(go_terms_after_correction, what="", sep="\n")

# Primero se hizo un analisis donde solo se devolvia el numero de hits por 
# termino, para asegurarse de que un id no obtuviera dos o mas 
# tablas de resultados
# go_terms_exit_status <- get_go_terms_info(go_terms_before_correction)
# length(go_terms_exit_status)
# [1] 2425
# table(go_terms_exit_status)
# go_terms_exit_status
# 1 
# 2425 

# De igual forma se comprobo las posibles columnas de los resultados
# para estos GO terms, para ver que coinciden con las que indican
# desde el sitio web de la API
# 
# go_terms_colnames <- get_go_terms_info(go_terms_before_correction)
# table(sapply(go_terms_colnames,length))
#   6    7    8    9 
# 128 1028 1137  132 
#
# unique(unlist(go_terms_colnames))
# [1] "id"         "isObsolete" "name"       "definition"
# [5] "synonyms"   "children"   "aspect"     "usage"     
# [9] "comment"
# En este caso tenemos todas las columnas disponibles,
# aunque la columna ancestors pasa a ser children
# Tambien la mayoria de los genes tienen 7 u 8 columnas.
#
# table(unlist(go_terms_colnames))
#     aspect   children    comment definition         id 
#       2425       1537        236       2425       2425 
# isObsolete       name   synonyms      usage 
#       2425       2425       1925       2425
#


# Como justamente las columnas de interes se encuentran presentes,
# se usaran estas para obtener la informacion faltante de los GO terms
go_terms_results <- get_go_terms_info(go_terms_after_correction)
go_terms_df<-data.frame(matrix(unlist(go_terms_results), ncol=5,byrow = T))

# Primero se evaluara que todos los GO terms posean el mismo ID en la info
#table(go_terms_df$X1==go_terms_df$X2)
# FALSE  TRUE 
#    17  2408
go_terms_df$X6 <-go_terms_df$X1==go_terms_df$X2
# Al parecer, si se buscan estos IDs en GO ontology, estan marcados como
# secundarios

# Para asegurar que todos estos IDs son secundarios, se buscaran si aparecen
# como tal
# go_terms_secids<-go_terms_df[go_terms_df$X6==FALSE,1]
# idssec_url <- paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",
#                    paste(go_terms_secids,collapse=","), "/secondaryids", sep="")
# idssec_results <- GET(idssec_url, accept("application/json"))
# 
# stop_for_status(idssec_results)
# 
# idssec_results <- fromJSON(toJSON(content(idssec_results)))

# Revisar si todos los IDs se obtuvieron, comparando los IDs que serian los 
# principales en results contra los obtenidos previamente
# setdiff(go_terms_df[go_terms_df$X6==FALSE,2],unlist(idssec_results$results$id))
# setdiff(unlist(idssec_results$results$id),go_terms_df[go_terms_df$X6==FALSE,2])

# Al parecer, todos los IDs que no coinciden son IDs secundarios, y hay dos
# IDs secundarios que comparten el mismo ID principal, GO:0006325
# Dado que este tipo de IDs son identicos en significo, seran reemplazados
# tanto su descripcion como su ID por la informacion del ID principal


# Por otra parte, los IDs obsoletos no deben de ser usados, por lo que se
# evaluara si existen en el set de datos

# Obtener los GO terms que ahora son obsoletos
# table(sapply(go_terms_results, function(x) x[3]))
# FALSE  TRUE 
#  2396    29

# Se revisara la informacion de cada termino obsoleto para evaluar 
# si seran eliminados, dado que en algunos casos se puede sugerir reemplazos

# Dado que no se encontro alguna herramienta para buscar los reemplazos,
# se hara buscando en su historia y manualmente

# go_terms_obsids<-go_terms_df[go_terms_df$X3==TRUE,1]
# idsobs_url <- paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",
#                     paste(go_terms_obsids,collapse=","), "/history", sep="")
# idsobs_results <- GET(idsobs_url, accept("application/json"))
# 
# stop_for_status(idsobs_results)
# 
# idsobs_results <- fromJSON(toJSON(content(idsobs_results)))

# Busquedas dias 19,20 del 08 del 2023
# GO:0000469  2023-04-24 replaced_by GO:0006364
# GO:0004652  2023-02-10 replaced_by GO:1990817
# GO:0005779  2022-08-26 replaced_by GO:0005778
# GO:0006073  2023-01-08 replaced_by GO:0044042
# GO:0006165  2023-01-22 The reason for obsoletion is that this 
#                        term represents a single step MF
# GO:0006379  2023-04-24 This term was obsoleted because it represents 
#                        a molecular function
# GO:0006471  2022-07-26 consider GO:1990404
# GO:0008022  2022-12-12 replaced_by GO:0005515
# GO:0015299  nnnn-nn-nn This term was obsoleted because it is an 
#                        unnecessary grouping term
# GO:0016307  2023-05-16 replaced_by GO:0052742
# GO:0018024  2023-02-28 consider GO:0140938
#             2023-02-28 consider GO:0140939
#             2023-02-28 consider GO:0140940
# GO:0018298  2022-04-24 Term describes modified products or 
#                        self-catalyzed processes
# GO:0030173  2022-08-26 replaced_by GO:0000139
# GO:0030176  2022-08-26 replaced_by GO:0005789
# GO:0031225  2022-08-26 replaced_by GO:0016020
# GO:0031307  2022-08-26 replaced_by GO:0005741
# GO:0031361  2022-08-26 replaced_by GO:0042651
# GO:0042779  2023-04-24 replaced_by GO:0042780
# GO:0043486  2022-07-18 consider GO:0006338
#             2022-07-18 consider GO:0140713
# GO:0043631  2023-05-20 The reason for obsoletion is that this 
#                        term represents a molecular function
# GO:0046658  2022-08-26 replaced_by GO:0005886
# GO:0046855  2023-01-22 consider GO:0043647
# GO:0048478  2023-01-22 consider GO:0031297 
# GO:0050072  2023-01-18 consider GO:0140933
#             2023-01-18 consider GO:0140932 
# GO:0055072  2023-01-08 consider GO:0006879
#             2023-01-08 consider GO:0060586
# GO:0070122  2022-05-07 replaced_by GO:0008233
# GO:0070940  2023-02-01 consider GO:0008420
# GO:0072321  2021-03-12 consider GO:0015031
#             2021-03-12 consider GO:0140597
# GO:0090503  2023-05-20 This term was obsoleted because it 
#                        represents a molecular function


# NOTA: CHECAR SI ESTOS IDS DE REEMPLAZO ESTAN EN EL ARCHIVO DE ENTRADA
# Esto porque algunos IDs no estan, y se debe de decidir si cambiarlos
# o no

# Dado que los IDs obsoletos sugieren reemplazos, solo seran reemplazados
# aquellos que proporcionen un unico reemplazo, y los demas seran eliminados




# Obtener los GO terms para todos los genes
# aqui se unieron todos los terminos de cada isoforma de cada gen
go_terms <- paste("Enrichment_analysis/GO_annotations/mptak_v6.1r1/",
                  "MpTak_v6.1_GO_annotation_1line_",
                  "first_isoform_non_sexchr_genes_",
                  "rev2_correction_sec_corrt_",
                  "obs_corrt.tsv",sep="")

#interproscan
#NOTA: se guardo como txt, aunque es un tsv
#go_terms <- paste("Enrichment_analysis/GO_annotations/interproscan/",
#                  "Mptakv6.1r2_go_annotation_1line_",
#                  "ref_iso_non_sexchr_genes.txt",sep="")



go_terms <- scan(go_terms, what="", sep="\n")
go_terms <- strsplit(go_terms, "\t")
# a cada vector de la lista, agregar el nombre del gen al que pertencen 
# los GO terms
names(go_terms) <- sapply(go_terms, function(x) {
  print(x)
  sub("^(Mp.*).1$", "\\1", x[1])
  })
# Remover el nombre de las isoformas de los vectores
go_terms <- lapply(go_terms, function(x) x[-1])

# Si existieran duplicaciones (de genes o de terminos), 
# fusionar y remover GO terms duplicados

dup_terms<-sapply(split(go_terms, names(go_terms)), function(sublist) {
  anyDuplicated(unlist(sublist))>0
})
# Evaluar si existen duplicaciones de terminos
table(dup_terms)

# Remover duplicaciones
go_terms <- lapply(split(go_terms, names(go_terms)), function(sublist) {
  unique(unlist(sublist))
})

dup_terms<-sapply(split(go_terms, names(go_terms)), function(sublist) {
  anyDuplicated(unlist(sublist))>0
})
# Evaluar si siguen existiendo duplicaciones de terminos
table(dup_terms)

# Evaluar distribucion del numero de GO terms de los genes
# La mayoria solo tienen 1 o 2
table(sapply(go_terms, length))

process_vector <- function(vec) {
  sub("^GO:(\\d+):.*$", "GO:\\1", vec)
}

# Guardar los GO terms con sus nombres/descripciones
term2name_list <- unique(unlist(go_terms))

# Eliminar los nombres/descripciones de los GO terms
go_terms <- lapply(split(go_terms, names(go_terms)), function(sublist) {
  unique(process_vector(unlist(sublist)))
})

# Revisar total de GO terms anotados
length(unlist(go_terms))


# Obtener los IDs que conformaran el universo al hacer el ORA analisis
# Se usaran todos los IDs de todos los genes, dado que de cualquier forma
# solo se podran usar aquellos que esten dentro de la anotacion (solo proteinas)
IDs_universe <- c("MpTak_v6.1r2_all_gene_ids.txt")
IDs_universe <- read.table(IDs_universe, header=F)
rownames(IDs_universe)<-IDs_universe$V2
# De cualquier forma, tambien se descartan los genes de cromosomas sexuales
IDs_universe<-IDs_universe[!grepl("^Mp(U|V)", rownames(IDs_universe)),2]


# Funcion para generar dataframes para interproscan enricher
create_terms_dfs<-function(go_terms_df,go_terms,go_cat){
  # Generar data.frame con los go terms y sus descripciones
  # Esto puede hacerse usando las descripciones originales
  # term2name_df<- data.frame(
  #   go_term =sub("^GO:(\\d+):.*$", "GO:\\1", term2name_list),
  #   go_name =sub("^GO:\\d+:(.*)$", "\\1", term2name_list)
  # )
  
  # O usando las descripciones nuevas obtenidas del
  # gene ontology (son sinonimos)
  # Dado que a partir de estas nuevas descripciones se hicieron
  # las correciones y se obtuvieron las categorias de los GO terms,
  # se obtara por usar estas
  
  # Generar visualizacion previa de num de GO terms para cada categoria
  print(table(go_terms_df$X5))
  
  if(go_cat==0){
    term2name_df<- data.frame(
      go_term = go_terms_df$X1,
      go_name = go_terms_df$X4
    )
  } else if(go_cat==1){
    term2name_df<- data.frame(
      go_term = go_terms_df[go_terms_df$X5=="biological_process",1],
      go_name = go_terms_df[go_terms_df$X5=="biological_process",4]
    )
  } else if(go_cat==2){
    term2name_df<- data.frame(
      go_term = go_terms_df[go_terms_df$X5=="cellular_component",1],
      go_name = go_terms_df[go_terms_df$X5=="cellular_component",4]
    )
  } else if(go_cat==3){
    term2name_df<- data.frame(
      go_term = go_terms_df[go_terms_df$X5=="molecular_function",1],
      go_name = go_terms_df[go_terms_df$X5=="molecular_function",4]
    )
  } else{
    return(0)
  }
  
  # Generar data.frame con los go terms y los genes a los cuales corresponden
  term2gene_df <- data.frame(
    go_term = unlist(go_terms),
    gene_id = rep(names(go_terms), sapply(go_terms, length))
  )
  rownames(term2gene_df) <- NULL
  
  # Si se tiene solo un tipo de categoria, ejecutar esta linea para eliminar
  # todas las filas que no esten en la categoria
  if(go_cat!=0){
    term2gene_df<-term2gene_df[term2gene_df$go_term %in% term2name_df[,1],]
  }
  return(list(term2gene_df,term2name_df))
}


# Guardar los dataframes con todas las categorias
terms_dfs<-create_terms_dfs(go_terms_df,go_terms,0)
term2gene_df<-terms_dfs[[1]]
term2name_df<-terms_dfs[[2]]
write.table(term2gene_df, "Enrichment_analysis/term2gene_df_all_cats.tsv",
            quote = F,col.names = T,row.names = F,sep = "\t")
write.table(term2name_df,"Enrichment_analysis/term2name_df_all_cats.tsv", 
            quote = F,col.names = T,row.names = F)


# Funcion para generar todos los plots
# 0 = all degs; 1 = Up; -1 = Down
create_enrichment_plots<-function(deseq2_results,term2gene_df,term2name_df,
                                  file_names,degs_set,IDs_universe,dpi_val,
                                  out_format){
  # Preparacion de genes ####
  DGE.results<-deseq2_results
  DGE.results<-DGE.results[!is.na(DGE.results$padj),]
  degs <- DGE.results[DGE.results$padj < 0.05,]
  
  if(degs_set==0){
    degs<-degs[abs(degs$log2FoldChange) > 1,]
  }else if(degs_set==(1)){
    degs<-degs[degs$log2FoldChange > 1,]
  } else if(degs_set==(-1)){
    degs<-degs[degs$log2FoldChange < -1,]
  } else{
    return(0)
  }

  degs_ids <- rownames(degs)
  
  # clusterProfiler ####
  #term2gene_df_direct_indirect<-buildGOmap(term2gene_df)
  
  ## ORA ####
  # Realizar el analisis de sobre representacion (ORA)
  ora_results <- enricher(
    # Vector con los genes de interes
    gene = degs_ids, 
    # Vector con el set de genes que seran el background
    universe = IDs_universe, 
    # NOTA: EL VALOR DE LOS CUTOFF PARECE NO INFLUIR EN LOS RESULTADOS,
    # DE IGUAL FORMA SE LES ASIGNARA EL CORTE DESEADO
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    #TERM2GENE = term2gene_df_direct_indirect,
    TERM2GENE = term2gene_df,
    TERM2NAME = term2name_df,
    # Es necesario especificar min y max de genes por categoria/termino
    # goseq no pide esto, minimo de 3 para coincidir con cluego
    # Este valor indica num de genes para la categoria ser considerada, no el
    # numero de DEGs que esten en la categoria, por lo tanto podemos seguir 
    # obteniendo hits de solo 2 genes o 1 gen
    minGSSize = 3,
    # NOTA: ORIGINALMENTE SE EMPLEO UN TAMANO MAXIMO POR CATEGORIA DE 1000000
    # GENES, PARA INCLUIR TODAS LAS CATEGORIAS, PERO AL CONTAR EL NUMERO DE 
    # GENES POR CATEGORIA PARECE QUE NINGUNO SUPERA LOS 2000 GENES, POR 
    # LO QUE SE CAMBIARA A 2000. FALTARIA REVISAR SI OCURRE ALGUN CAMBIO, 
    # PERO AL HACER PEQUENAS PRUEBAS PARECE QUE NO AFECTA PARA NADA EL PONER
    # AHORA 2000, POR LO ANTERIOR MENCIONADO
    maxGSSize = 2000
  )
  
  # Hacer una visualizacion previa de los resultados
  # dotplot(ora_results)
  
  # Ordernar y filtrar por padjust
  # Este paso se hace para el grafico de redes y el upsetplot, porque estos
  # piden como input el objeto ora_results
  ora_results@result <- ora_results@result[order(ora_results@result$p.adjust),]
  ora_results@result <- ora_results@result[ora_results@result$p.adjust<0.05,]
  
  # Obtener los resultados del ORA
  go_results_df <- data.frame(ora_results@result)
  
  # Verificar el numero de terminos que pasan el filtro
  if(nrow(go_results_df)==0){
    print("Error: No hay terminos que pasen el umbral de 0.05 o no se obtuvieron")
    return(0)
  }
  
  # Al parecer, GeneRatio es lo siguiente (NumDEGsInCat/AllDEGsInAllCat)
  # Al ser un string, es posible convertirlo a un valor numerico
  # Si se desea, se debe ejecutar lo siguiente
  # go_results_df$GeneRatio <- sapply(go_results_df$GeneRatio, function(str) {
  #   parts <- strsplit(str, "/")[[1]]
  #   as.numeric(parts[1]) / as.numeric(parts[2])
  # })
  
  
  # Obtener el top10 de terminos sobrerepresentados
  #top_10_go_terms_2<-slice_min(go_results_df, order_by=p.adjust, n=10)
  
  # Al evaluar top_10_go_terms_2 Y top_10_go_terms_1, se tienen los mismos 
  # resultados
  # Al parecer, se obtiene los mismos resultados para goseq y clusterprofiler
  
  #setdiff(rownames(go_results_df),goResults$category)
  #setdiff(goResults$category,rownames(go_results_df))
  
  # Ordenar los GO terms en caso de no estar ordenados por padjust
  #go_results_df<-go_results_df[order(go_results_df$p.adjust),]
  
  
  # Obtener num de genes en cada go term, y la proporcion
  go_results_df$numGenesInCat <- sapply(go_results_df$BgRatio, function(str) {
    parts <- strsplit(str, "/")[[1]]
    as.numeric(parts[1])
  })
  go_results_df$hitsPerc<-(go_results_df$Count*100)/go_results_df$numGenesInCat
  
  ### Guardar tabla ####
  write.table(go_results_df, file_names[5],
              quote = F,col.names = T,row.names = F,sep = "\t")
  
  
  # NOTA: La funcion de saltos de linea ya existe implementada en scales,
  # pero como hay terminos que pueden ser muy grandes, mejor se optara por
  # anadir puntos suspensivos al final
  # Generar funcion que realiza saltos de linea cada n caracteres, en caso de
  # encontrar un espacio que permita separar la linea
  addNewline_format <- function(x){
    x<-gsub("(.{40,}?)\\s", "\\1\n", x)
  }
  
  # Generar funcion que anade puntos supensivos al final de 40 caracteres, si 
  # es que encuentra espacios
  adddots_format <- function(x){
    # x es la tabla de results de ORA
    
    # Agregar puntos suspensivos
    new_description<-apply(x, MARGIN=1,function(y){
      description<-gsub("(.{40,}?)\\s", "\\1\n", y["Description"])
      is_long<-length(strsplit(description, "\n")[[1]])
      description<-strsplit(description, "\n")[[1]][1]
      if(is_long>1){
        description<-paste(description,"...",sep="")
      }
      # Anadir GO term al final
      description<-paste(description," (",y["ID"],")",sep="")
      return(description)
    }
    )
    new_description
    
  }
  
  # Agregar nuevo formato a go_results_df
  go_results_df$Description<-adddots_format(go_results_df)
  
  ### Grafico de barras ####
  # Tamano de texto default
  text_size<-12

  # Modifir las descripciones de los GO terms agregando line breaks, para
  # evitar tener nombres demasiado grandes
  #go_results_df$Description<-addNewline_format(go_results_df$Description)
  # NOTA: la funcion label_wrap permite hacer lo mismo, por lo que
  # se puede omitir lo anterior, y como a veces no quedan muy bien y se ven
  # amontonado el texto, mejor se usaran puntos suspensivos
  
  #go_results_df
  
  # Obtener el top 10, y convertir los padjust a -log10(padjust),
  # para mejorar su visualizacion con el gradiente de colores 
  # Si se usa head, se descartan todos los empates,
  # si se usa slice_min, devuelve todo incluyendo empates, por lo que
  # pueden haber mas de 10 terminos
  #top_10_go_terms <- head(go_results_df, n=10)
  top_10_go_terms <- slice_min(go_results_df, order_by=p.adjust, n=10)
  top_10_go_terms$p.adjust <- (-log10(top_10_go_terms$p.adjust))
  #print(top_10_go_terms$p.adjust)
  
  # Establecer como limite el valor de 25, 10**-25
  top_10_go_terms$p.adjust <- ifelse(top_10_go_terms$p.adjust > 25, 
                                     25, top_10_go_terms$p.adjust)
  
  # Generar grafico de barras
  barplot_top10_terms<-ggplot(top_10_go_terms, 
                              aes(x=Count,
                                  y=reorder(Description,p.adjust,decreasing=F),
                                  fill=p.adjust))+
    geom_col()+
    
    scale_y_discrete(expand = expansion(mult = c(0.075, 0.075)),
                     #labels = label_wrap(40)
                     )+
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))+
    #expansion(mult = c(0.02, 0.02)), primer valor anade 2% de 
    # espacio abajo/izquierda y el segundo 2% de espacio arriba/derecha
    theme(axis.text.x=element_text(family="Arial",size=text_size),
          axis.text.y=element_text(family="Arial",size=text_size,
                                   hjust = 1, margin = margin(l=5,r=5)),
          legend.text=element_text(family="Arial",size=text_size,
                                   hjust = 0),
          legend.title=element_text(family="Arial",size=text_size,
                                    margin = margin(b=3)),
          axis.title.x = element_text(family="Arial",size=16),
          axis.title.y = element_text(family="Arial",size=16)
    ) +
    scale_fill_gradient(breaks = seq(0, 25, 5),
                        #low="azure1", high = "royalblue4",
                        labels = expression(10^0,10^-5,10^-10,
                                            10^-15,10^-20,10^-25),
                        limits= c(0, 25))+
    # Generar los labels
    labs(x="Count", y="Terms", fill="adjust p-value")
  
  
  #barplot_top10_terms
  # Guardar barplot
  ggsave(filename=paste("barplot_top10_terms_",
                        file_names[1],sep=""), 
         plot=barplot_top10_terms,
         path="figuras/without_sex_chromosomes/ORA/",
         dpi=dpi_val, device = out_format, 
         units ="cm",width = 30, height = 15)
  
  
  
  ### dotplot ####
  # Generar dotplot
  dotplot_top10_terms<-ggplot(top_10_go_terms,
                              aes(x=hitsPerc, 
                                  y=reorder(Description,hitsPerc), 
                                  colour=p.adjust, size=Count)) +
    geom_point() +
    #scale_y_discrete(labels = label_wrap(40))+
    scale_x_continuous(breaks = seq(0,100,20))+
    coord_cartesian(xlim = c(0,100)) +
    theme(axis.text.x=element_text(family="Arial",size=text_size),
          axis.text.y=element_text(family="Arial",size=text_size,
                                   hjust = 1, margin = margin(l=5,r=5)),
          legend.text=element_text(family="Arial",size=text_size,
                                   hjust = 0),
          legend.title=element_text(family="Arial",size=text_size,
                                    margin = margin(b=3)),
          axis.title.x = element_text(family="Arial",size=16),
          axis.title.y = element_text(family="Arial",size=16))+
    scale_colour_gradient(breaks = seq(0,25,5),
                          labels = expression(10^0,10^-5,10^-10,10^-15,
                                              10^-20,10^-25),
                          limits= c(0,25))+
    # Usar scale_size_area asegura que valores de 0 se les
    # asigne tamanos de 0, y max_size es el tamano maximo que se
    # aplica al valor mas grande
    scale_size_area(max_size = 8)+
    # x es GeneRatio % (NumDEGsInCat/AllGenesInCat)
    labs(x="Gene Ratio", y="GO term", 
         colour="adjust p-value", size="Count")
  
  #dotplot_top10_terms
  # Guardar dotplot
  ggsave(filename=paste("dotplot_top10_terms_",
                        file_names[2],sep=""), 
         plot=dotplot_top10_terms,
         path="figuras/without_sex_chromosomes/ORA/",
         dpi=dpi_val, device = out_format, 
         units ="cm",width = 22, height = 12)
  
  
  ### upsetplot ####
  # Se genera un set separado de ora_results, para anadir el formato de ... a
  # este grafico, y evitar que el cnetplot lo tenga
  upsetplot_ora_results<-ora_results
  upsetplot_ora_results@result$Description<-adddots_format(upsetplot_ora_results@result)
  
  # Para usarse, requiere instalar previamente ggupset
  upsetplot_top10_terms<-enrichplot::upsetplot(upsetplot_ora_results,
                                               n=nrow(top_10_go_terms))
  
  #upsetplot_top10_terms
  # Guardar upsetplot
  ggsave(filename=paste("upsetplot_top10_terms_",
                        file_names[3],sep=""), 
         plot=upsetplot_top10_terms,
         path="figuras/without_sex_chromosomes/ORA/",
         dpi=dpi_val, device = out_format, 
         units ="cm",width = 30, height = 12)
  
  
  ### grafico de redes ####
  # Obtener los DEGs que contienen los top 10 GO terms
  # Se puede hacer de esta forma
  #col_vals<-DGE.results[intersect(degs_ids,
  #                                unique(unlist(ora_results@geneSets[top_10_go_terms$ID]))),"log2FoldChange"]
  #names(col_vals)<-intersect(degs_ids,unique(unlist(ora_results@geneSets[top_10_go_terms$ID])))
  
  # o de esta forma
  col_vals<-DGE.results[unique(unlist(strsplit(top_10_go_terms$geneID, "/"))),
                        "log2FoldChange"]
  names(col_vals)<-unique(unlist(strsplit(top_10_go_terms$geneID, "/")))
  
  
  # Revisar que ambas formas dan el mismo set de genes
  #setdiff(unique(unlist(strsplit(top_10_go_terms$geneID, "/"))),
  #        names(col_vals))
  #setdiff(names(col_vals),
  #        unique(unlist(strsplit(top_10_go_terms$geneID, "/"))))
  
  # Asignar los limites a los valores de log2foldchange
  col_vals <- ifelse(col_vals > 10, 10, col_vals)
  col_vals <- ifelse(col_vals < -10, -10, col_vals)
  
  # Dado que estos graficos no se pueden personalizar facilmente, se 
  # modificara la descripcion de ora_results agrengando line breaks a las
  # descripciones con mas de 40 caracteres
  ora_results@result$Description<-addNewline_format(ora_results@result$Description)
  # Generar grafico de redes
  cnetplot_top10_terms<-cnetplot(ora_results,
                                 node_label="category",
                                 # El valor default es dar las 5 categorias mas 
                                 # significativas
                                 showCategory = nrow(top_10_go_terms),
                                 color.params = list(foldChange = col_vals,
                                                     edge=T))
  
  # Obtener las posiciones de los GO terms en el grafico de redes
  
  pg <- ggplot_build(cnetplot_top10_terms)
  print(nrow(top_10_go_terms))
  # Guardar previamente los colores de cada termino, si es deseado
  #cat_colour<-unique(pg$data[[1]]$edge_colour)
  
  pg <- pg$data[[2]]
  pg<-pg[1:nrow(top_10_go_terms),]
  pg$size <- top_10_go_terms$Count
  # Agregar el orden en el top de los go terms, asi como sus ids y sus 
  # descripciones formateadas con line breaks
  pg$top_lvl<-1:nrow(top_10_go_terms)
  pg$id <- top_10_go_terms$ID
  pg$description <- ora_results@result$Description[1:nrow(top_10_go_terms)]
  
  
  # Tambien se pueden elegir otros colores para los terminos, pero se usara el 
  # default
  #cat_colour<-pal_unikn_dark[1:10]
  #names(cat_colour)<-pg$description
  #cat_colour<-unlist(cat_colour)
  
  # Agregar a cada color la descripcion que le corresponde
  #cat_colour<- c("#00B0F6","#D89000","#00BFC4","#A3A500","#FF62BC",
  #               "#00BF7D","#F8766D","#9590FF","#39B600","#E76BF3")
  #cat_colour<-cat_colour[1:nrow(top_10_go_terms)]
  cat_colour<-hue_pal()(nrow(top_10_go_terms))
  print(cat_colour)
  names(cat_colour)<-pg$description
  # Modificar los colores de las conexiones y los nodos, para que coincidan
  # con el termino que le corresponde
  cnetplot_top10_terms<-cnetplot_top10_terms+ guides(edge_colour="none")+
    theme_void()+
    theme(legend.text=element_text(family="Arial",size=text_size,
                                   hjust = 0),
          legend.title=element_text(family="Arial",size=text_size,
                                    margin = margin(b=3))
    ) +
    scale_colour_gradient2(breaks = seq(-10, 10, 4),
                           labels = seq(-10, 10, 4),
                           limits= c(-10, 10),
                           low="#77AADD",high="#993F30",
                           name=expression("log"[2]*" fold change"))+
    scale_colour_manual(aesthetics = "edge_colour",
                        values=cat_colour)+
    scale_size_area(max_size = 15)+
    new_scale_colour() +
    geom_point(data=pg,aes(x=x,y=y,size=size,
                           colour = reorder(factor(description),top_lvl)),
               shape=19,stroke=0.5)+
    scale_colour_manual(values=cat_colour)+
    labs(colour="GO term Description",size="Size")
  # Pasar hacia adelante los nombres de los GO terms
  cnetplot_top10_terms$layers<-cnetplot_top10_terms$layers[c(1, 2, 3, 5, 4)]
  #cnetplot_top10_terms
  #Guardar el cnetplot
  ggsave(filename=paste("cnetplot_top10_terms_",
                        file_names[4],sep=""), 
         plot=cnetplot_top10_terms,
         path="figuras/without_sex_chromosomes/ORA/",
         dpi=dpi_val, device = out_format, 
         units ="cm",width = 40, height = 30)
}

cat_names<-c("all_cats","BP_cat","CC_cat","MF_cat")
figures_formats<-c("svg","png")
dpi_formats<-c(1200,300)

degs_set_val<-(0) # Cambiar valor al que sea el deseado (0,-1,1)

# Forma automatizada para generar todas las figuras
# Si se desea cambiar a interproscan, solo anadir los datos y editar los 
# archivos de salida agregandoles el termino _ipscan_, o eliminandolo si ya 
# esta presente y se esta usando otra anotacion
for (idx_frmt in 1:2){
  print(idx_frmt)
  for (idx_cat in 1:length(cat_names)){
    # Imprimir categoria usada
    cat_name<-cat_names[idx_cat]
    print(cat_name)
    # Agregar la palabra "up" o "down", si se tiene ese set de DEGs
    if(degs_set_val==1){
      cat_name<-paste(cat_name,"_up",sep="")
    } else if(degs_set_val==(-1)){
      cat_name<-paste(cat_name,"_down",sep="")
    }
    
    # Generar dataframes a usar
    # 0 = all cats; 1 = BP; 2 = CC; 3 = MF
    terms_dfs<-create_terms_dfs(go_terms_df,go_terms,(idx_cat-1))
    term2gene_df<-terms_dfs[[1]]
    term2name_df<-terms_dfs[[2]]
    
    create_enrichment_plots(DGE.results_dlk,term2gene_df,term2name_df,
                            c(rep(paste("dlk_vs_WT_",
                                        cat_name,".",figures_formats[idx_frmt],
                                        sep=""),4),
                              paste("Enrichment_analysis/ORA_GO_terms/",
                                    "dlk_vs_WT_",
                                    cat_name, "_results.tsv", sep="")),
                            degs_set=degs_set_val,IDs_universe,
                            dpi_formats[idx_frmt],figures_formats[idx_frmt])
    
    create_enrichment_plots(DGE.results_dlk1,term2gene_df,term2name_df,
                            c(rep(paste("dlk1_vs_WT_",
                                        cat_name,".",figures_formats[idx_frmt],
                                        sep=""),4),
                              paste("Enrichment_analysis/ORA_GO_terms/",
                                    "dlk1_vs_WT_",
                                    cat_name,"_results.tsv",sep="")),
                            degs_set=degs_set_val,IDs_universe,
                            dpi_formats[idx_frmt],figures_formats[idx_frmt])
    
    create_enrichment_plots(DGE.results_dlk25,term2gene_df,term2name_df,
                            c(rep(paste("dlk25_vs_WT_",
                                        cat_name,".",figures_formats[idx_frmt],
                                        sep=""),4),
                              paste("Enrichment_analysis/ORA_GO_terms/",
                                    "dlk25_vs_WT_",
                                    cat_name,"_results.tsv",sep="")),
                            degs_set=degs_set_val,IDs_universe,
                            dpi_formats[idx_frmt],figures_formats[idx_frmt])
    
    create_enrichment_plots(DGE.results_OE13,term2gene_df,term2name_df,
                            c(rep(paste("OE13_vs_WT_",
                                        cat_name,".",figures_formats[idx_frmt],
                                        sep=""),4),
                              paste("Enrichment_analysis/ORA_GO_terms/",
                                    "OE13_vs_WT_",
                                    cat_name,"_results.tsv",sep="")),
                            degs_set=degs_set_val,IDs_universe,
                            dpi_formats[idx_frmt],figures_formats[idx_frmt])
    
  }
}





# goseq ####
# NOTA: Dado que goseq no produce diferencias con clusterprofiler,
# se puede omitir esta parte

# Crear vector binario que indique los genes que son diferenciales
# del grupo de todos los genes que pasaron los filtros de results
gene.vector <- row.names(DGE.results) %in% degs_ids %>% as.integer
names(gene.vector) <- row.names(DGE.results)

# Obtener la mediana de los tamanos para todos los genes
txdb <- makeTxDbFromGFF("MpTak_v6.1r2.gff", format = ("gff"))
txsByGene<-transcriptsBy(txdb,"gene")
lengthData<-median(width(txsByGene))


# El grafico generado tiene una tendencia a disminuir, revisar
pwf<-nullp(gene.vector,bias.data = lengthData[names(gene.vector)],
           plot.fit = F)
plotPWF(nullp)
# Non-statistical methods for determining DE, such as using a fold-change 
# cutoff, can show a decreasing trend in the proportion of DE as a function 
# of gene length 

# En el articulo de goseq, dice que hacer un cutoff del lfc puede producir 
# esa tendencia. Sin embargo, aunque no se realice ese cutoff, se
# sigue obteniendo la tendencia.


# Si hiciera contra el basemean, la tendencia es a aumentar
nullp(gene.vector,bias.data = DGE.results$baseMean)


# TEST.CATS parece no funcionar, por lo que a pesar de que goseq encuentre las
# categorias de los GO terms, no los puede separar para este caso 
# usando gene2cat
goResults <-goseq(pwf,
                  #test.cats = c("GO:MF"),
                  gene2cat = go_terms)

# Si se usa el metodo hipergeometrico, no se aplica la correccion 
# del tamano del gen
goResults <-goseq(pwf,
                  method = "Hypergeometric",
                  gene2cat = go_terms)

# Realizar corte de categorias con un numero de DEGs mayor a 0
# Estos genes deben ser descartados antes de hacer la correcion del p-value
goResults<-goResults[!goResults$numDEInCat==0,]

# Realizar correciones de los p-values
goResults$over_represented_padj<-p.adjust(goResults$over_represented_pvalue, method="BH")

goResults$under_represented_padj<-p.adjust(goResults$under_represented_pvalue, method="BH")

# Obtener el porcentaje de DEGs en un GO term con respecto al total de genes
# para ese GO term
goResults$hitsPerc<-(goResults$numDEInCat/goResults$numInCat)*100


# Obtener los terminos sobrerepresentados despues de corte de padj<0.05
goResults<-goResults[goResults$over_represented_padj<0.05,]

# Generar el TOP 10 de terminos sobrerepresentados
top_10_go_terms<-slice_min(goResults,order_by=over_represented_padj,n=10)

top_10_go_terms<-slice_min(goResults[goResults$ontology=="BP",],
                           order_by=over_represented_padj,n=10)

top_10_go_terms<-slice_min(goResults[goResults$ontology=="MF",],
                           order_by=over_represented_padj,n=10)

top_10_go_terms<-slice_min(goResults[goResults$ontology=="CC",],
                           order_by=over_represented_padj,n=10)

# Generar dotplot para el top10
ggplot(top_10_go_terms,
       aes(x=hitsPerc, y=reorder(term,hitsPerc), colour=over_represented_padj, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="adjust p-value", size="Count")


# Obtener informacion del termino
GOTERM[[top_10_go_terms$category[1]]]

GOTERM[["GO:0016021"]]

# Nota: Algunos terminos para la la version 6.1r1 pasaron a ser
# obsoletos o son IDs secundarios/sinonimos
# Parece ser tambien que los IDs secundarios son tomados de forma separada,
# por lo que a pesar de aparecer como sinonimos tienen sus propios p-values
# Tambien se obtienen p-values para terminos obsoletos, pero goseq no puede
# encontrar informacion por el status de estos terminos

# Info de GO.db
# Mappings were based on data provided by: Gene Ontology 
# http://current.geneontology.org/ontology/go-basic.obo With a date stamp 
# from the source of: 2023-01-01

# Para esta version del gene ontology se hace el goseq

