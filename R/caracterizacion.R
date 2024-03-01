#' @title Medidas de resumen
#' @description Función que devuelve las medidas de resumen de los datos de un gpkg
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param data_path  Nombre del archivo GPKG con los datos.
#' @param index_path Nombre del archivo del indice
#' Por defecto es "index.RDS"
#' @param seed Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto.
#' @param outlier valor lógico que define si sacar los outliers.  
#' Por defecto es TRUE y no saca los outliers
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param freecor Cantidad de núcleos de CPU a dejar libres al ejecutar la función. 
#' Por defecto es 1
#' @param n_col Espera un entero con el numero de la columna a resumir en caso de que tenga más de una (Datos apareados). 
#' Por defecto es 1
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera achicar el numero de capas a este valor. 
#' Por defecto es 0
#' @return Retorna un dataframe con las medidas de resumen de cada capa
#' @details Primero setea los parámetros para realizar el calculo en múltiples núcleos. 
#' Luego lista y recorre las capas para obtener la medidas de resumen. 
#' Por ultimo retorna el dataframe con las medidas de cada capa.
#' @examples summary_gpkg("D:/Facultad/Tesina/Temp/clean","Elevacion.gpkg",freecor = 2)
#' @export
summary_gpkg <- function(path, 
                         data_path, 
                         index_path="index.RDS", 
                         seed=NULL,
                         outlier=T,
                         iqr_p=3,
                         freecor=1, 
                         ncol=1,
                         debug = 0){
  #Paquetes
  require(doParallel)
  require(dplyr)
  require(sf)
  #Del total de nucleos del procesador le resta los indicados y arma el cluster
  cor <- detectCores()-freecor
  cl <- makeCluster(cor)
  registerDoParallel(cl)
  mcaffinity(1:cor)
  on.exit(stopCluster(cl))
  
  #Obtiene una lista de todas las capas del gpkg
  layers <- st_layers(file.path(path, data_path))
  #Cuando debug > 0 entramos en un entorno de prueba donde se achica el numero de layers de acuerdo al valor indicado
  if (debug > 0 && debug < length(layers$name)) {
    layers$name <- layers$name[1:debug]
  }
  
  #Obtiene el indice almancenado de manera local
  index <- readRDS(file.path(path, index_path)) |> distinct(layer,.keep_all = T)
  #Datos para crear el log
  #Crea el nombre del log de la corrida
  log_name<-"log_summary_gpkg.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  #Corre los calculos en paralelo con foreach
  ret <- foreach(layer = layers$name,
                 .combine = rbind,
                 .packages = c("sf","dplyr")) %dopar%{
                   
                   tryCatch({
                     #Levanta la capa
                     gpkg <- st_read(file.path(path, data_path),layer = layer)
                     gpkg <- gpkg |> 
                       rename(data=colnames(gpkg)[ncol])
                     gpkg$data <- gpkg$data |> as.numeric()
                     # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
                     Q1 <- quantile(gpkg$data, 0.25)
                     Q3 <- quantile(gpkg$data, 0.75)
                     # Calcular el rango intercuartílico (IQR)
                     IQR <- Q3-Q1
                     # Calcular los límites inferior y superior para identificar los outliers
                     limite_inferior <- Q1 - iqr_p * IQR
                     limite_superior <- Q3 + iqr_p * IQR
                     # Crear un campo logico que identifique si el dato es un outlier
                     gpkg <- gpkg %>%
                       mutate(outlier = if_else(
                         data<limite_inferior|data>limite_superior
                         ,T
                         ,F))
                     rm(list=c("Q1","Q3","IQR","limite_inferior","limite_superior"))
                     #Eliminando los Outliers
                     if (outlier==F) {
                       gpkg<-gpkg |> filter(outlier==F)
                     }
                     #Calculo de la distancia a todos los puntos
                     #Si tiene más de 10.000 puntos hace un sample de 10.000 puntos
                     if (nrow(gpkg)>10000) {
                       if(!is.null(seed)){
                         set.seed(seed)  # Semilla
                       }
                       gpkg_s <- sample_n(gpkg,size = 10000)
                     }else{
                       gpkg_s<-gpkg
                     }
                     #Obtiene desde el indice la superficie de la capa
                     superficie<-index$bound_sup[index$layer==layer]
                     #Obtiene una matriz con la distancia entre todos los puntos y se guarda el triangulo superior
                     distancia <- st_distance(x=gpkg_s)
                     y <- as.numeric(as.matrix(distancia)[upper.tri(distancia)])
                     x<-gpkg$data |> as.numeric()
                     #Dataframe de retorno con todas las medidas
                     data.frame(layer=layer
                                ,n=length(x)
                                ,media=mean(x)
                                ,mediana=median(x)
                                ,min=min(x)
                                ,max=max(x)
                                ,iqr=IQR(x)
                                ,sd=sd(x)
                                ,var=sd(x)**2
                                ,cv=sd(x)/mean(x)
                                ,asim=moments::skewness(x)
                                ,kurt=moments::kurtosis(x)
                                ,dist_media=mean(y)
                                ,dens_m = 10000*length(gpkg$data)/superficie
                     )
                   }
                   ,error=function(e){
                     #si tira un error abre el log, registra la capa y el error, luego lo cierra
                     sink(file.path(log_path,log_name),append = T)
                     print(paste0("Error en la capa: ",layer))
                     print(e)
                     sink()
                   })
                 }
  registerDoSEQ()
  return(ret)
}

#' @title Medidas de resumen de correlación
#' @description Función que devuelve las medidas de resumen de la correlación de un gpkg
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param data_path  Nombre del archivo GPKG con los datos.
#' @param outlier valor lógico que define si sacar los outliers.  
#' Por defecto es TRUE y no saca los outliers
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param freecor Cantidad de núcleos de CPU a dejar libres al ejecutar la función. 
#' Por defecto es 1
#' @param n_col Espera un entero con el numero de la columna a resumir en caso de que tenga más de una (Datos apareados). 
#' Por defecto es 1
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera achicar el numero de capas a este valor. 
#' Por defecto es 0
#' @return Retorna un dataframe con las medidas de resumen de cada capa
#' @details Primero setea los parámetros para realizar el calculo en múltiples núcleos. 
#' Luego lista y recorre las capas para obtener la medidas de resumen. 
#' Por ultimo retorna el dataframe con las medidas de cada capa.
#' @examples summary_gpkg("D:/Facultad/Tesina/Temp/clean", "Elevacion.gpkg",freecor = 2)
#' @export
cor_summary_gpkg <- function(path, 
                             data_path,
                             outlier=T, 
                             iqr_p=3,
                             freecor=1, 
                             ncol=1,
                             debug = 0){
  #Paquetes
  require(doParallel)
  require(dplyr)
  require(sf)
  require(gstat)
  require(tidyr)
  cor <- detectCores()-freecor
  cl <- makeCluster(cor)
  registerDoParallel(cl)
  mcaffinity(1:cor)
  on.exit(stopCluster(cl))
  
  #Obtiene una lista de todas las capas del gpkg
  layers <- st_layers(file.path(path, data_path))
  #Cuando debug > 0 entramos en un entorno de prueba donde se achica el numero de layers de acuerdo al valor indicado
  if (debug > 0 && debug < length(layers$name)) {
    layers$name <- layers$name[1:debug]
  }
  
  #Crea el nombre del log de la corrida
  log_name<-"log_cor_summary_gpkg.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Corre los calculos en paralelo con foreach
  ret <- foreach(layer = layers$name,
                 .combine = bind_rows,
                 .packages = c("sf","tidyverse","gstat")) %dopar%{
                   tryCatch({
                     #Carga el gpkg indicado en los paramentros y renombra la columna de los datos a un nombre generico
                     gpkg <- st_read(file.path(path, data_path),layer = layer)
                     gpkg <- gpkg |> 
                       rename(data=colnames(gpkg)[ncol])
                     gpkg$data <- gpkg$data |> as.numeric()
                     # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
                     Q1 <- quantile(gpkg$data, 0.25)
                     Q3 <- quantile(gpkg$data, 0.75)
                     # Calcular el rango intercuartílico (IQR)
                     IQR <- Q3-Q1
                     # Calcular los límites inferior y superior para identificar los outliers
                     limite_inferior <- Q1 - iqr_p * IQR
                     limite_superior <- Q3 + iqr_p * IQR
                     # Crear un campo logico que identifique si el dato es un outlier
                     gpkg <- gpkg %>%
                       mutate(outlier = if_else(
                         data<limite_inferior|data>limite_superior
                         ,T
                         ,F))
                     rm(list=c("Q1","Q3","IQR","limite_inferior","limite_superior"))
                     #Eliminando los Outliers
                     if (outlier==F) {
                       gpkg<-gpkg |> filter(outlier==F)
                     }
                     #Transforma el sf object a un dataframe
                     gpkg_df <- as.data.frame(gpkg) |> select(data)
                     #Obtiene el variograma 
                     x <- variogram(data ~ 1, gpkg, cloud = T, cutoff =100)
                     x_df <- as.data.frame(x)
                     x_df <- x_df |> 
                       filter(dist<=100) |> 
                       mutate(lag = cut(dist, breaks = c(0,25,50,75,100),include.lowest = T))
                     #Crea una columna id con la posicion de cada dato para poder hacer un join con el variograma
                     gpkg_df$id <- 1:nrow(gpkg_df)
                     x_df<-left_join(x_df, gpkg_df, by = c("left"= "id"))
                     x_df <- x_df |> rename(left_data = "data")
                     x_df<-left_join(x_df, gpkg_df, by = c("right"= "id"))
                     x_df <- x_df |> rename(right_data = "data")
                     #Los agrupa por los distintos lags y obtiene la correlacion para cada lag
                     y <- x_df |>
                       group_by(lag) |>
                       summarise(r = cor(left_data, right_data)) |>
                       pivot_wider(names_from = lag,values_from = r)
                   }
                   ,error=function(e){
                     #si tira un error abre el log, registra la capa y el error, luego lo cierra
                     sink(file.path(log_path,log_name),append = T)
                     print(paste0("Error en la capa: ",layer))
                     print(e)
                     sink()
                     #Crea un dataframe vacio para poder hacer el rbind
                     y<-data.frame(NA,NA,NA,NA) |> setNames(c("[0,25]", "(25,50]", "(50,75]", "(75,100]"))
                   })
                   #Agrega el nombre de la capa a los datos obtenidos ya que posteriormente vamos a necesitar hacer un merge
                   cbind(layer,y)
                 }
  #Renombra y ordena las columnas del df de retorno
  ret<-ret %>% 
    rename(lag_0_25="[0,25]", lag_25_50="(25,50]", lag_50_75="(50,75]", lag_75_100="(75,100]") %>% 
    relocate(layer,sort(colnames(.)))
  registerDoSEQ()
  return(ret)
}

#' @title Distancia promedio al vecino mas cercano
#' @description Función que devuelve las distancias promedio al vecino mas cercano de un gpkg
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param data_path  Nombre del archivo GPKG con los datos.
#' @param outlier valor lógico que define si sacar los outliers.  
#' Por defecto es TRUE y no saca los outliers
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param freecor Cantidad de núcleos de CPU a dejar libres al ejecutar la función
#' Por defecto es 1
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera achicar el numero de capas a este valor. 
#' Por defecto es 0
#' @return Retorna un dataframe con las distancias promedio al vecino mas cercano de cada capa.
#' @details Primero setea los paramentros para realizar el calculo en múltiples núcleos. 
#' Luego setea los parámetros para crear el log. 
#' Luego lista y recorre las capas para obtener las distancias. 
#' Por ultimo retorna el dataframe con las medidas de cada capa.
#' @examples dist_nn_gpkg("D:/Facultad/Tesina/Temp/clean","Elevacion.gpkg",freecor = 2)
#' @export
dist_nn_gpkg <- function(path, 
                         data_path,
                         outlier=T, 
                         iqr_p=3, 
                         freecor=1, 
                         ncol=1, 
                         debug = 0){
  #Paquetes
  require(doParallel)
  require(sf)
  #Del total de nucleos del procesador le resta los indicados y arma el cluster
  cor <- detectCores()-freecor
  cl <- makeCluster(cor)
  registerDoParallel(cl)
  mcaffinity(1:cor)
  on.exit(stopCluster(cl))
  
  #Obtiene una lista de todas las capas del gpkg
  layers <- st_layers(file.path(path, data_path))
  #Cuando debug > 0 entramos en un entorno de prueba donde se achica el numero de layers de acuerdo al valor indicado
  if (debug > 0 && debug < length(layers$name)) {
    layers$name <- layers$name[1:debug]
  }
  
  #Crea el nombre del log de la corrida
  log_name<-"log_dist_nn_gpkg.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Recorre cada capa en paralelo
  ret <- foreach(layer = layers$name,
                 .combine = rbind,
                 .packages = c("sf","dplyr")) %dopar%{
                   tryCatch({
                     #Obteniendo los datos de la capa
                     gpkg <- st_read(file.path(path, data_path),layer = layer)
                     gpkg <- gpkg |> 
                       rename(data=colnames(gpkg)[ncol])
                     gpkg$data <- gpkg$data |> as.numeric()
                     # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
                     Q1 <- quantile(gpkg$data, 0.25)
                     Q3 <- quantile(gpkg$data, 0.75)
                     # Calcular el rango intercuartílico (IQR)
                     IQR <- Q3-Q1
                     # Calcular los límites inferior y superior para identificar los outliers
                     limite_inferior <- Q1 - iqr_p * IQR
                     limite_superior <- Q3 + iqr_p * IQR
                     # Crear un campo logico que identifique si el dato es un outlier
                     gpkg <- gpkg %>%
                       mutate(outlier = if_else(
                         data<limite_inferior|data>limite_superior
                         ,T
                         ,F))
                     rm(list=c("Q1","Q3","IQR","limite_inferior","limite_superior"))
                     #Eliminando los Outliers
                     if (outlier==F) {
                       gpkg<-gpkg |> filter(outlier==F)
                     }
                     #Calculo de la distancia al vecino más cercano
                     nn <- st_nearest_feature(gpkg)
                     dist_nn <- st_distance(gpkg, gpkg[nn, ], by_element = T)
                   }
                   ,error=function(e){
                     #si tira un error abre el log, registra la capa y el error, luego lo cierra
                     sink(file.path(log_path,log_name),append = T)
                     print(paste0("Error en la capa: ",layer))
                     print(e)
                     sink()
                     #Crea un dataframe vacio para poder hacer el rbind
                     data.frame(layer=layer, dist_prom_nn = NA)
                   })
                   
                   #Dataframe de retorno con todas las medidas
                   data.frame(layer=layer, dist_prom_nn = dist_nn |> as.numeric() |> mean())
                 }
  registerDoSEQ()
  return(ret)
}