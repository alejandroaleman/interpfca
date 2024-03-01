#' @title Registrar archivos de origen
#' @description Funcion que crea un dataframe con los archivos de los datos de utilizados para crear el dataframe
#' @inheritParams get_raw_gpkg
#' @inheritParams get_raw_shp
#' @inheritParams get_raw_csv
#' @param raw_data_id Espera un string con el id del drive del archivo de origen
#' @param raw_data_path Espera un string con el id del drive de la ubicacion del archivo de origen
#' @param from Espera un dataframe con los archivos de origen agregados antes, si es el primero lo crea vacio por default
#' @return Retorna un dataframe que cada fila contiene un archivo de origen
#' @details Primero obtiene el dribble del archivo en el drive
#' Luego agrega una linea al dataframe con el nombre, id y ubicacion del archivo de origen
#' Se debe utilizar de manera recursiva
#' @examples from <- add_from_data(raw_data_id, raw_data_path)
#' from <- add_from_data(raw_data_id_2, raw_data_path_2, from)
add_from_data <- function(raw_data_id,
                          raw_data_path,
                          from = data.frame(name=character(),
                                            id=character(),
                                            path = character())){
  #Paquetes
  require(googledrive)
  #Obteniendo la metadata del archivo
  dl <- drive_get(as_id(raw_data_id))
  #Armado del df con los datos del archivo que necesitamos para el indice
  from <- from %>% add_row(name=dl$name,id=dl$id,path=raw_data_path)
  return(from)
}


#' @title Obtener datos crudos
#' @description Funcion utilizada para obtener los datos crudos desde el drive
#' @inheritParams get_raw_shp
#' @inheritParams get_raw_csv
#' @inheritParams get_raw_gpkg
#' @param raw_data_id Espera un string con el id del archivo de origen
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @return Retorna un dribble con la ruta local de los datos crudos
#' @details Obtiene el dribble del archivo de origen y lo descarga
#' @examples datos <- get_raw_data(raw_data_id,temp_dir)
get_raw_data <- function(raw_data_id, 
                         temp_path){
  #Paquetes
  require(googledrive)
  require(tibble)
  #Obtiene la metadata del archivo
  metadata <- drive_get(as_id(raw_data_id))
  #Carga el path
  path <-file.path(temp_path, metadata$name)
  #Si no está de manera local descarga el archivo
  tryCatch({
    drive_download(as_id(metadata$id), path = path, overwrite = FALSE)
  }
  ,error = function(e){ #Si el ya esta descargado solo crea el dribble
    print("El archivo ya esta descargado localmente")
  }
  )
  #Se agrega la ruta local al dribble
  dl <- as_dribble(as_id(metadata$id))
  dl<-dl |> add_column(local_path=path,.after = "name")
  return(dl)
}

#' @title Obtener archivos Shapefile
#' @description Funcion utilizada para obtener los datos crudos en formato shape desde el drive
#' @param file Espera un string con el nombre del archivo a descargar
#' @param raw_data_path Espera un string con el id de la ubicacion del archivo de origen
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @return Retorna una lista de dos entradas con los datos y los archivos de origen
#' @details Primero obtiene los dribbles del archivo de origen, 
#' Luego lo descarga y arma el dataframe con la info de origen
#' Por ultimo importa el archivo a un sf object y arma la lista con el sf object y el df de origen
#' @examples datos <- get_raw_shp(file, raw_data_path, temp_dir)
#' @export
get_raw_shp <- function(file, 
                        raw_data_path, 
                        temp_path){
  #Paquetes
  require(googledrive)
  require(dplyr)
  require(sf)
  #Obtiene una lista de todos los archivos que hay en la carpeta del drive con los datos crudos
  lista_base <- drive_ls(as_id(raw_data_path),recursive = FALSE)
  #Filtra cada archivo del shape para obtener el id, lo descarga y arma el dataframe con los datos de origen que luego se utilizara en el indice
  #.sph
  shp <- filter(lista_base, name==paste0(file,".shp"))
  shp_dl <- get_raw_data(shp$id,temp_path)
  from <- add_from_data(shp$id,raw_data_path)
  #.shx
  shx <- filter(lista_base, name==paste0(file,".shx"))
  get_raw_data(shx$id,temp_path)
  from <- add_from_data(shx$id,raw_data_path,from)
  #.dbf
  dbf <- filter(lista_base, name==paste0(file,".dbf"))
  get_raw_data(dbf$id,temp_path)
  from <- add_from_data(dbf$id,raw_data_path,from)
  #.prj
  prj <- filter(lista_base, name==paste0(file,".prj"))
  get_raw_data(prj$id,temp_path)
  from <- add_from_data(prj$id,raw_data_path,from)
  #Importa el shapefile descargado a un sf object, arma la lista de doble entrada y lo retorna
  data <- st_read(shp_dl$local_path)
  output <- list("data"=data,"from"=from)
  return(output)
}

#' @title Obtener archivos csv
#' @description Funcion utilizada para obtener los datos crudos en formato csv desde el drive
#' @param file Espera un string con el nombre del archivo a descargar
#' @param raw_data_path Espera un string con el id de la ubicacion del archivo de origen
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @return Retorna una lista de dos entradas con los datos y los archivos de origen
#' @details Primero obtiene los dribbles del archivo de origen, 
#' Luego lo descarga y arma el dataframe con la info de origen
#' Por ultimo importa el archivo a un sf object y arma la lista con el sf object y el df de origen
#' @examples datos <- get_raw_shp(file, raw_data_path, temp_dir)
#' @export
get_raw_csv <- function(file, 
                        raw_data_path, 
                        temp_path){
  #Paquetes
  require(googledrive)
  require(dplyr)
  require(sf)
  #Obtiene una lista de todos los archivos que hay en la carpeta del drive con los datos crudos
  lista_base <- drive_ls(as_id(raw_data_path),recursive = FALSE)
  #Filtra el archivo csv para obtener el id, lo descarga y arma el dataframe con los datos de origen que luego se utilizara en el indice
  #.csv
  csv <- filter(lista_base, name==paste0(file,".csv"))
  csv_dl <- get_raw_data(csv$id,temp_path)
  from <- add_from_data(csv$id,raw_data_path)
  #Importa el csv descargado a un sf object, arma la lista de doble entrada y lo retorna
  data <- st_read(csv_dl$local_path)
  output <- list("data"=data,"from"=from)
  return(output)
}

#' @title Obtener archivos geopackage
#' @description Funcion utilizada para obtener los datos crudos en formato .gpkg desde el drive
#' @param file Espera un string con el nombre del archivo a descargar
#' @param raw_data_path Espera un string con el id de la ubicacion del archivo de origen
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @return Retorna una lista de dos entradas con los datos y los archivos de origen
#' @details Primero obtiene los dribbles del archivo de origen, 
#' Luego lo descarga y arma el dframe con la info de origen
#' Por ultimo importa el archivo a un sf object y arma la lista con el sf object y el df de origen
#' @examples datos <- get_raw_shp(file, raw_data_path, temp_dir)
#' @export
get_raw_gpkg <- function(file, 
                         raw_data_path, 
                         temp_path){
  #Paquetes
  require(googledrive)
  require(dplyr)
  require(sf)
  #Obtiene una lista de todos los archivos que hay en la carpeta del drive con los datos crudos
  lista_base <- drive_ls(as_id(raw_data_path),recursive = FALSE)
  #Filtra el archivo gpkg para obtener el id, lo descarga y arma el dataframe con los datos de origen que luego se utilizara en el indice
  #.gpkg
  gpkg <- filter(lista_base, name==paste0(file,".gpkg"))
  gpkg_dl <- get_raw_data(gpkg$id,temp_path)
  from <- add_from_data(gpkg$id,raw_data_path)
  #Importa el gpkg descargado a un sf object, arma la lista de doble entrada y lo retorna
  data <- st_read(gpkg_dl$local_path)
  output <- list("data"=data,"from"=from)
  return(output)
}

#' @title Obtener los datos del perimetro del lote
#' @description Funcion utilizada para obtener los datos del perimetro del lote cargado en el drive como gpkg
#' @param file Espera un string con el nombre del archivo a descargar
#' @param raw_data_path Espera un string con el id de la ubicacion del archivo de origen
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @return Retorna un dataframe con el nombre del archivo, id de origen y la superficie
#' @details Primero obtiene los dribbles del archivo del perimetro, 
#' Luego lo descarga y importa el archivo a un sf object 
#' Por ultimo arma un dataframe con el nombre del perimetro, id en el drive y la superficie
#' @examples get_boundary_data(file, raw_data_path, temp_dir)
#' @export
get_boundary_data <- function(file, 
                              raw_data_path, 
                              temp_path){
  #Paquetes
  require(googledrive)
  require(dplyr)
  require(sf)
  #Obtiene una lista de todos los archivos que hay en la carpeta del drive con los datos crudos
  lista_base <- drive_ls(as_id(raw_data_path),recursive = FALSE)
  #Filtra el archivo gpkg para obtener el id y lo descarga
  #.gpkg
  gpkg <- filter(lista_base, name==paste0(file,".gpkg"))
  gpkg_dl <- get_raw_data(gpkg$id,temp_path)
  #Importa el archivo descargado y arma el dataframe de retorno con el nombre del archivo, si id en el drive y la superficie que abarca el perimetro
  sup <- st_read(gpkg_dl$local_path) |> st_area() |> as.numeric()
  output <- data.frame(bound_file=paste0(file,".gpkg"), bound_id=gpkg$id ,bound_sup=sup[1])
  return(output)
}


#' @title Cargar datos como capas de un gpkg
#' @description Funcion que carga el geopackage con cada dato en una capa y lo guarda localmente
#' @param data espera un sf object con el dato a guardar
#' @param name Espera un string con el nombre de la capa a guardar
#' @param from_df espera un dataframe con los datos de los archivos originales donde fué tomado
#' Se recomienda usar el dataframe de la lista obtenida con get_raw_shp
#' @param temp_path Espera un String con la ruta de la carpeta temporal donde se descargan y crean los distintos archivos
#' @param data_type espera un String con el tipo de dato a guardar para guardarlo en su correspondiente carpeta con ese nombre
#' @return Retorna un sf object pasado como parametro solo para poder usar el pipe, no tiene uso aun
#' @details Primero obtiene el directorio donde guardar los datos limpios, si no existe lo crea
#' Luego obtiene el indice, si no existe lo crea
#' Luego agrega el gpkg como una capa del tipo de dato correspondiente, si no existe el archivo lo crea
#' Por ultimo actualiza, elimina los repetidos y guarda el indice localmente
#' @examples put_gpkg_local(data, "SYM_field_01", from, "F:/Facultad/Tesina/Data processor/Temp", "mapas_rendimiento.gpkg")
#' @export
put_gpkg_local <- function(data, 
                           name, 
                           from_df, 
                           bound_df, 
                           temp_path, 
                           datatype){
  #Paquetes
  require(dplyr)
  require(sf)
  #Ruta local donde se guardan los datos ya procesados, si no existe la crea
  clean_path <- file.path(temp_path,"clean")
  dir.create(clean_path,showWarnings = F)
  #Arma la ruta con la ubicación del indice
  index_path <- file.path(clean_path,"index.RData")
  #Obtiene el indice
  index <- data.frame(layer=character(),
                      type=character(),
                      bound_file=character(),
                      bound_id=character(),
                      bound_sup=double(),
                      from=data.frame(name=character(),id=character(),path = character()))
  tryCatch({
    index <- readRDS(index_path)
  }
  ,error = function(e){ #Si el indice no está creado, lo crea
    print("Se creo un indice nuevo")
  }
  ,warning = function(w){
    print("Se creo un indice nuevo")
    saveRDS(index,file = index_path)
  }
  )
  #Guarda el archivo localmente, si ya está creado le agrega una capa
  path <- file.path(clean_path,datatype)
  st_write(data,path, layer=name, drive = "GPKG", delete_layer = T)
  #Actualiza el indice
  index <- index |> 
    add_row(layer=name,
            type=datatype,
            bound_file=bound_df$bound_file,
            bound_id=bound_df$bound_id,
            bound_sup=bound_df$bound_sup,
            from.name=from_df$name, 
            from.id=from_df$id,
            from.path=from_df$path)
  index <- index[!duplicated(index), ]
  #Guarda el indice
  saveRDS(index,file = index_path)
  return(data)
}
