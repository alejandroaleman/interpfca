#' @title Obtener los Mejores Modelos de Predicción por IDW con Validación Cruzada
#' @description Esta función calcula y devuelve los mejores modelos de predicción por el método del inverso a la distancia (IDW) para cada capa utilizando validación cruzada.
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param file_name Nombre del archivo GPKG con los datos.
#' @param outlier valor lógico que define si sacar los outliers.  
#' Por defecto es TRUE y no saca los outliers
#' @param n_col Número de columna a resumir en caso de tener múltiples columnas (e.g. datos apareados). 
#' Por defecto, es 1.
#' @param seed Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto.
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param nfold Numero de folds a utilizar en la validación cruzada. 
#' Por defecto es 5
#' @param nmax Numero máximo de vecinos utilizados. 
#' Por defecto es 10
#' @param nmin Numero mínimo de vecinos utilizados. 
#' Por defecto es 2
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera indicar que layers utilizar. 
#' Puede ser un vector o un entero, por defecto es 0.
#' @return Un dataframe con el mejor modelo de predicción por IDW y los datos de la validación cruzada para cada capa.
#' @details La función setea el archivo del log. 
#' Luego lista las capas del GPKG y arma la grilla de tuneo con diferentes valores de idp para cada capa. 
#' Luego recorre la grilla de tuneo, para cada combinación de layer-idp obtiene la validación cruzada, 
#' utilizando el método de k-folds y arma un dataframe donde cada fila tiene su combinación layer-idp y los resultados de la validación cruzada
#' Por ultimo calcula el error cuadrático medio de cada combinación layer-idp, filtra por capa las que tienen menor rmse y lo retorna
#' @examples get_idw_res(path="~/Documentos/Tesina/Data/clean", file_name="Elevacion.gpkg", seed = 701408733)
#' @export
get_idw_res <- function(path,
                        file_name,
                        sum_path,
                        outlier=T,
                        n_col=1,
                        seed=NULL,
                        iqr_p = 3,
                        nmax = 10,
                        nmin = 2,
                        nfold = 5,
                        debug = 0){
  #Paquetes
  require(sf)
  require(dplyr)
  require(tidyr)
  require(gstat)
  require(rsample)
  log_name<-"get_idw_res.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Traemos los datos del summary
  summary <- readRDS(sum_path)
  
  #Obtiene una lista de todas las capas del gpkg
  file_path<-file.path(path,file_name)
  layers <- st_layers(file_path)
  
  #Si debug es un vector solo trae los layers indicados en el vector
  #Si es un valor mayor a cero trae solo el layer indicado
  if(length(debug)>1){
    layers$name <- layers$name[debug]
  }else if (debug > 0 && debug<=length(layers$name)) {
    layers$name <- c(layers$name[debug])
  }
  
  #Seteamos los hiperparamentros, para aplicar IDW, en una grilla
  tune_grid <- expand_grid(
    layer = layers$name,#Se indica el nombre de la capa
    idp = seq(1, 3, by = 0.25)#Se setean idp desde 1 a 3 aumentando de a 0.25
  ) |> 
    group_by(layer) |> 
    nest()
  
  #Obtenemos los resultados del ajuste con los hiperparametros seteados en tune_grid
  #Recorriendo cada uno de manera paralela y se los almacena en tune_res
  columns<-c("layer", "idp", "data", "rmse","ndrop%", "time")
  tune_res <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
  rm(columns)
  imax <- nrow(tune_grid)
  for(i in 1:imax){
    #Arranca el contador de tiempo
    time_start<-Sys.time()
    #Obtenemos la combinacion layer-idp en variables separadas
    layer <- tune_grid$layer[i]
    l_idp <- tune_grid[[2]][[i]]
    #Traemos el valor de los vecinos
    layer_t<-layer
    dist_nn <- summary |> filter(layer==layer_t) %>% .$dist_prom_nn
    #Cargamos el set de datos y estandarizamos el nombre de la columna que contiene los datos
    data <- st_read(file_path,layer = layer)
    data <- data |> 
      rename(data_col=colnames(data)[n_col]) |> 
      mutate(data_col=as.numeric(data_col),
             x = st_coordinates(data)[, "X"],
             y = st_coordinates(data)[, "Y"]) |> 
      mutate(x=x-min(x),y=y-min(y)) |> 
      st_drop_geometry() |> st_as_sf(coords=c("x","y"))
    
    # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
    Q1 <- quantile(data$data_col, 0.25)
    Q3 <- quantile(data$data_col, 0.75)
    # Calcular el rango intercuartílico (IQR)
    IQR <- Q3-Q1
    # Calcular los límites inferior y superior para identificar los outliers
    limite_inferior <- Q1 - iqr_p * IQR
    limite_superior <- Q3 + iqr_p * IQR
    # Crear un campo logico que identifique si el dato es un outlier
    data <- data %>%
      mutate(outlier = if_else(
        data_col<limite_inferior|data_col>limite_superior
        ,T
        ,F))
    rm(list=c("Q1","Q3","IQR","limite_inferior","limite_superior"))
    #Eliminando los Outliers
    if (outlier==F) {
      data<-data |> filter(outlier==F)
    }
    
    # Dividiendo los datos en nfolds
    if(!is.null(seed)){
      set.seed(seed)  # Semilla
    }
    folds <- vfold_cv(data, v = nfold)
    #Sacamos de memoria las variables que no se usan
    rm(list = c("data"))
    # Realizar la validación cruzada
    columns <- c("layer","idp","data","rmse","ndrop%")
    xv <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
    rm(columns)
    maxidp <- nrow(l_idp)
    for (j in 1:maxidp){
      idp <- l_idp[[1]][[j]]
      # Inicializar una df para almacenar los resultados
      columns <- c("pred","observed","residual","geom","fold")
      ret <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
      for (k in 1:nfold) {
        
        # Crear conjuntos de entrenamiento y prueba para el fold k
        train_data <- training(folds[[1]][[k]])#Armando el subconjunto de datos para entrenar
        test_data <- testing(folds[[1]][[k]])#Armando el subconjunto de datos para realizar la validacion cruzada
        
        # Entrenando el modelo
        model <- gstat(id = layer, formula = data_col ~ 1, data = train_data,model = NULL, set = list(idp = idp),maxdist = dist_nn*10)
        
        tryCatch({
          print(paste0("Layer: ",layer," - Idp: ",idp," - Fold: ",k))
          # Realizar predicciones en el subconjunto de prueba
          if(data.class(test_data$geom)=="sfc_MULTIPOINT"){
            predictions <- predict(model, test_data |> st_cast("POINT"),debug.level=-1)
          }
          else{
            predictions <- predict(model, test_data,debug.level=-1)
          }
          names(predictions)[1]<-"pred"
          #Armando un dataframe con el predicho vs observado y su error residual, se agrega tambien el fold utilizado
          x<-cbind(pred = predictions$pred, 
                   observed = test_data$data_col, 
                   residual = test_data$data_col-predictions$pred,
                   geom = predictions$geometry |> st_as_sf(), 
                   fold = k)
        }
        ,error=function(e){
          #Si tira un error abre el log, registra la capa y el error, luego lo cierra
          sink(file.path(log_path,log_name),append = T)
          print(paste0("Error en la capa: ",layer))
          print(paste0("idp: ",idp," - Fold: ",k))
          print(e)
          sink()
        })
        #Agregamos lo calculado al df de retorno
        ret <- rbind(ret,x)
        #Sacamos de la memoria las variables que no se usan
        rm(list = c("model", "train_data", "test_data"))
      }
      lyxv <-tibble(layer=layer,idp=idp,data=list(ret)) |> 
        group_by(layer, idp) |>
        mutate(rmse = 
                 sqrt(
                   mean(
                     filter(ret,!is.na(residual))$residual^2)),
               'ndrop%'=(100*length(filter(ret,!is.na(residual))$residual)/length(ret$residual)) %>% round(digits=1))
      
      xv <- rbind(xv,lyxv)
    }
    xv <- xv |>  
      group_by(layer) |>
      filter(rmse == min(rmse))
    #Finaliza el contador de tiempo
    time_end<-Sys.time()
    xv <- cbind(xv,data.frame(time=difftime(time_end, time_start, units = "secs")))
    tune_res <- rbind(tune_res, xv)
  }
  return(tune_res)
}

#' @title Obtener los Mejores Modelos de Predicción por RF con Validación Cruzada
#' @description Esta función calcula y devuelve los mejores modelos de predicción por el método random forest (RF) para cada capa utilizando validación cruzada.
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param file_name Nombre del archivo GPKG con los datos.
#' @param outlier valor lógico que define si sacar los outliers. 
#' Por defecto es TRUE y no saca los outliers
#' @param n_col Número de columna a resumir en caso de tener múltiples columnas (e.g. datos apareados). 
#' Por defecto, es 1.
#' @param seed Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto. 
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param nfold Numero de folds a utilizar en la validación cruzada
#' Por defecto es 5
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera indicar que layers utilizar. 
#' Puede ser un vector o un entero, por defecto es 0.
#' @return Un dataframe con el mejor modelo de predicción por RF y los datos de la validación cruzada para cada capa.
#' @details La función setea el archivo del log. 
#' Luego lista las capas del GPKG y arma la grilla de tuneo con diferentes valores de numero de vecinos a utilizar para cada capa. 
#' Luego recorre la grilla de tuneo, para cada combinación de layer-knn obtiene la validación cruzada, 
#' utilizando el método de k-folds y arma un dataframe donde cada fila tiene su combinación layer-knn y los resultados de la validación cruzada
#' Por ultimo calcula el error cuadrático medio de cada combinación layer-knn, filtra por capa las que tienen menor rmse y lo retorna
#' @examples get_rf_res(path="~/Documentos/Tesina/Data/clean", file_name="Elevacion.gpkg", seed = 701408733)
#' @export
get_rf_res <- function(path,
                       file_name,
                       outlier=T,
                       n_col=1,
                       seed=NULL,
                       iqr_p=3,
                       nfold = 5,
                       debug = 0){
  #Paquetes
  require(sf)
  require(dplyr)
  require(tidyr)
  require(ranger)
  require(rsample)
  log_name <- "get_rf_res.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Obtiene una lista de todas las capas del gpkg
  file_path<-file.path(path,file_name)
  layers <- st_layers(file_path)
  
  #Si debug es un vector solo trae los layers indicados en el vector
  #Si es un valor mayor a cero trae solo el layer indicado
  if(length(debug)>1){
    layers$name <- layers$name[debug]
  }else if (debug > 0 && debug<=length(layers$name)) {
    layers$name <- c(layers$name[debug])
  }
  
  #Seteamos los hiperparamentros, para aplicar RF, en una grilla
  tune_grid <- expand_grid(
    layer = layers$name,#Se indica el nombre de la capa
    knn = seq(5, 35, by = 10)#Se setean el numero de vecinos a utilizar, segun Sekulic(pag 22) a partir de 25 practicamente no varia el rmse, aunque depende del tamaño de la muestra por lo que esto puede variar
  ) |> 
    group_by(layer) |> 
    nest()
  ntrees <- 250 #Numero de arboles para RFsi Segun sekulic pag 8 y 9
  #Obtenemos los resultados del ajuste con los hiperparametros seteados en tune_grid
  #Recorriendo cada uno de manera paralela y se los almacena en tune_res
  columns<-c("layer", "knn", "ntrees", "data", "rmse","ndrop%", "time")
  tune_res <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
  rm(columns)
  imax <- nrow(tune_grid)
  for(i in 1:imax){
    tryCatch({
      #Arranca el contador de tiempo
      time_start<-Sys.time()
      #Obtenemos la combinacion layer-knn en variables separadas
      layer <- tune_grid$layer[i]
      l_knn <- tune_grid[[2]][[i]]
      #Cargamos el set de datos y estandarizamos el nombre de la columna que contiene los datos
      data <- st_read(file_path,layer = layer)
      if(class(data$geom)[1]=="sfc_MULTIPOINT"){
        data <- data |> st_as_sf() |> st_cast("POINT")
      }
      data <- data |> 
        rename(data_col=colnames(data)[n_col]) |> 
        mutate(data_col=as.numeric(data_col),
               x = st_coordinates(data)[, "X"],
               y = st_coordinates(data)[, "Y"]) |> 
        mutate(x=x-min(x),y=y-min(y)) |> 
        st_drop_geometry() |> st_as_sf(coords=c("x","y"))
      
      # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
      Q1 <- quantile(data$data_col, 0.25)
      Q3 <- quantile(data$data_col, 0.75)
      # Calcular el rango intercuartílico (IQR)
      IQR <- Q3-Q1
      # Calcular los límites inferior y superior para identificar los outliers
      limite_inferior <- Q1 - iqr_p * IQR
      limite_superior <- Q3 + iqr_p * IQR
      # Crear un campo logico que identifique si el dato es un outlier
      data <- data %>%
        mutate(outlier = if_else(
          data_col<limite_inferior|data_col>limite_superior
          ,T
          ,F))
      rm(list=c("Q1","Q3","IQR","limite_inferior","limite_superior"))
      #Eliminando los Outliers
      if (outlier==F) {
        data<-data |> filter(outlier==F)
      }
      
      ##Nombre de la variable
      vname <- "data_col"
      # Realizar la validación cruzada
      columns<-c("layer","knn","ntrees","data","ndrop%","rmse")
      xv <-as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
      rm(columns)
      maxknn <- nrow(l_knn)
      for (j in 1:maxknn){
        knn <- l_knn[[1]][[j]]
        ## Feature names
        ft_names <- c(
          paste0("vals_nn", 1:(knn-1)),
          paste0("dist_nn", 1:(knn-1))
        )  
        ##Calculamos las distancias entre los vecinos
        nn_data <- nngeo::st_nn(data, data, k = knn, returnDist = T, sparse = T)
        
        ##Obtenemos los valores de los vecinos
        nn_data$vals <- lapply(nn_data$nn, function(x) data[[vname]][x])
        
        ##Combinamos los datos y asignamos los nombres de las columnas de los vecinos
        nn_data <- cbind(
          do.call(rbind, nn_data$vals)[, -1],
          do.call(rbind, nn_data$dist)[, -1]
        )
        colnames(nn_data) <- ft_names
        temp_data <- cbind(data, nn_data)
        #Formula del modelo
        fm <- paste0(vname, " ~ ", paste(ft_names, collapse = "+"))
        
        # Dividiendo los datos en nfolds
        if(!is.null(seed)){
          set.seed(seed)  # Semilla
        }
        folds <- vfold_cv(temp_data, v = nfold)
        #Sacamos de memoria las variables que no se usan
        rm(temp_data)
        
        # Inicializar una df para almacenar los resultados 
        columns <- c("pred","observed","residual","fold")
        ret <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
        rm(columns)
        for (k in 1:nfold) {
          
          # Crear conjuntos de entrenamiento y prueba para el fold k
          train_data <- training(folds[[1]][[k]])#Armando el subconjunto de datos para entrenar
          test_data <- testing(folds[[1]][[k]])#Armando el subconjunto de datos para realizar la validacion cruzada
          
          # Entrenando el modelo
          model <- train_data %>%  
            st_drop_geometry() %>% 
            na.omit() %>% 
            ranger::ranger(fm, data = ., num.trees = ntrees, importance = "impurity")
          
          # Realizar predicciones en el subconjunto de prueba
          tryCatch({
            x <- NULL
            predictions <- predict(model, test_data |> st_drop_geometry(),debug.level=-1)
            names(predictions)[1]<-"pred"
            #Armando un dataframe con el predicho vs observado y su error residual, se agrega tambien el fold utilizado
            x<-cbind(pred = predictions$pred,
                     observed = test_data$data_col, 
                     residual = test_data$data_col-predictions$pred,
                     fold = k)
          }
          ,error=function(e){
            #Si tira un error abre el log, registra la capa y el error, luego lo cierra
            sink(file.path(log_path,log_name),append = T)
            print(paste0("Error en el layer: ", layer," - knn: ", knn," - fold: ", k))
            print(e)
            sink()
          })
          if (is.null(x)) {
            x<-cbind(pred = NA ,observed = NA, residual = NA, fold = k)
          }
          #Agregamos lo calculado al df de retorno
          ret <- rbind(ret,x)
          #Sacamos de la memoria las variables que no se usan
          rm(list = c("model", "train_data", "test_data"))
        }
        lyxv <-tibble(layer = layer, knn = knn, ntrees = ntrees, data = list(ret)) |> 
          group_by(layer, knn, ntrees) |>
          mutate(rmse = 
                   sqrt(
                     mean(
                       filter(ret,!is.na(residual))$residual^2)),
                 'ndrop%'=(100*length(filter(ret,!is.na(residual))$residual)/length(ret$residual)) %>% round(digits=1))
        xv <- rbind(xv,lyxv)
        #Sacamos de memoria las variables que no se usan
        rm(list = c("folds"))
      }
      #Sacamos de memoria las variables que no se usan
      rm(list = c("data"))
      xv <- xv |>  
        group_by(layer) |>
        filter(rmse == min(rmse))
      #Finaliza el contador de tiempo
      time_end<-Sys.time()
      xv <- cbind(xv,data.frame(time=difftime(time_end, time_start, units = "secs")))
      # xv
    }
    ,error=function(e){
      #Si tira un error abre el log, registra la capa y el error, luego lo cierra
      sink(file.path(log_path,log_name),append = T)
      print(paste0("Error en la capa: ",layer))
      print(e)
      
      sink()
    })
    
    if(nrow(xv)==0){
      time_end<-Sys.time()
      ret <- ret |> add_row(pred=0,observed=0,residual=0,fold=0)
      xv <- data.frame(layer = layer, knn = 0, ntrees = 0, data = list(ret)) |> 
        group_by(layer, knn, ntrees) |>
        nest() |> 
        cbind(data.frame(rmse=0,time=difftime(time_end, time_start, units = "secs")))
      xv$data <- lapply(xv$data, function(df) {
        names(df) <- sub("^data\\.", "", names(df))
        return(df)
      })
    }
    tune_res<-rbind(tune_res,xv)
  }
  return(tune_res)
}

#' @title Obtener los Mejores Modelos de Predicción por OK con Validación Cruzada
#' @description Esta función calcula y devuelve los mejores modelos de predicción por el método ordinary kriging (OK) para cada capa utilizando validación cruzada.
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param file_name Nombre del archivo GPKG con los datos.
#' @param param Dataframe con los parámetros iniciales para correr el modelo
#' @param outlier valor lógico que define si sacar los outliers. 
#' Por defecto es TRUE y no saca los outliers
#' @param trainwo logico que determina si entrenar el modelo sin los outliers, Por defecto es FALSE
#' @param n_col Número de columna a resumir en caso de tener múltiples columnas (e.g. datos apareados). 
#' Por defecto, es 1.
#' @param seed Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto. 
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param nfold Numero de folds a utilizar en la validación cruzada
#' Por defecto es 5
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera indicar que layers utilizar. 
#' Puede ser un vector o un entero, por defecto es 0.
#' @return Un dataframe con el mejor modelo de predicción por OK y los datos de la validación cruzada para cada capa.
#' @details La función setea el archivo del log.
#' Luego carga en una grilla de tuneo los hiperparametros especificados en param.  
#' Luego recorre la grilla de tuneo, obtiene la validación cruzada, 
#' utilizando el método de k-folds y arma un dataframe donde cada fila tiene su layer, tiempo de ejecución y los resultados de la validacion cruzada y lo retorna
#' @examples get_idw_res(path="~/Documentos/Tesina/Data/clean", file_name="Elevacion.gpkg",param="~/Documentos/Tesina/Data/variograms/ECa30_vgm.RDS", seed = 701408733)
#' @export
get_ok_res <- function(path,
                       file_name,
                       param,
                       outlier=T,
                       trainwo=F,
                       iqr_p=3,
                       n_col=1,
                       seed=NULL,
                       nfold = 5,
                       debug = 0){
  #Paquetes
  require(sf)
  require(dplyr)
  require(tidyr)
  require(rsample)
  require(gstat)
  #Crea el nombre del log de la corrida
  log_name<-"get_ok_res.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Ruta temporal donde se guardan las salidas parciales
  temp <- file.path(path,"OK")
  dir.create(temp,showWarnings = F)
  temp <- file.path(temp,"temp")
  dir.create(temp,showWarnings = F)
  
  #Obtiene una lista de todas las capas del gpkg
  file_path<-file.path(path,file_name)
  
  #Si debug es un vector solo trae los layers indicados en el vector
  #Si es un valor mayor a cero trae solo el layer indicado
  if(length(debug)>1){
    param <- param[debug,]
  }else if (debug > 0 && debug<=length(param$layer)) {
    param <- param[debug,]
  }
  
  #Obtenemos los resultados del ajuste con los hiperparametros seteados en param
  #Recorriendo cada uno y se los almacena en ok
  columns <- c("layer", "formula", "data", "fit_param", "time", "rmse", "ndrop%")
  ok<-as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
  imax <- nrow(param)
  for(i in 1:imax){
    #Arranca el contador de tiempo
    time_start<-Sys.time()
    #Obtenemos los datos de la capa correspondiente
    layer <- param$layer[[i]]
    fm_str <- param$formula[[i]]
    fm <- formula(fm_str)
    sill <- param$sill[[i]]
    nugget <- param$nugget[[i]]
    range <- param$range[[i]]
    cutoff <-param$cutoff[[i]]
    width <- param$width[[i]]
    dist_nn <- param$dist_nn[[i]] |> as.numeric()
    
    #Cargamos el set de datos y estandarizamos el nombre de la columna que contiene los datos
    data <- st_read(file_path,layer = layer)
    data <- data %>% 
      rename(data_col=colnames(data)[n_col]) %>% 
      mutate(data_col=as.numeric(data_col),
             x = st_coordinates(data)[, "X"], #Extrae las coordenadas "x" e "y" de la columna "geom"
             y = st_coordinates(data)[, "Y"]) %>% 
      mutate(x=x-min(x),y=y-min(y)) %>% 
      st_drop_geometry() %>% st_as_sf(coords=c("x","y")) 
    data <- data %>% 
      mutate(x = st_coordinates(data)[, "X"],
             y = st_coordinates(data)[, "Y"])
    
    # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
    Q1 <- quantile(data$data_col, 0.25)
    Q3 <- quantile(data$data_col, 0.75)
    # Calcular el rango intercuartílico (IQR)
    IQR <- Q3-Q1
    # Calcular los límites inferior y superior para identificar los outliers
    limite_inferior <- Q1 - iqr_p * IQR
    limite_superior <- Q3 + iqr_p * IQR
    # Crear un campo logico que identifique si el dato es un outlier
    data <- data %>%
      mutate(outlier = if_else(
        data_col<limite_inferior|data_col>limite_superior
        ,T
        ,F))
    #Eliminando los Outliers
    if (outlier==F) {
      data<-data |> filter(outlier==F)
    }
    # Dividiendo los datos en nfolds
    if(!is.null(seed)){
      set.seed(seed)  # Semilla
    }
    folds <- vfold_cv(data, v = nfold)
    #Sacamos de memoria las variables que no se usan
    rm(list = c("data","Q1","Q3","IQR","limite_inferior","limite_superior"))
    #Realizar la validación cruzada
    #Generamos el modelo del variograma
    m<-vgm(psill = sill-nugget, nugget = nugget, range = range, model = "Mat")
    
    # Inicializar una df para almacenar los resultados 
    columns <- c("pred","observed","residual","geom","fold")
    ret <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
    columns <- c("fold","fit_psill","fit_range","fit_nugget","fit_kappa")
    fit <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
    #Inicializamos xv con valores de error para que pueda hacer el rbind de todas formas
    time_end<-Sys.time()
    xv <-tibble(layer=layer,
                formula="Error",
                data=list(data.frame()),
                fit_param=list(data.frame()),
                time=difftime(time_end, time_start, units = "secs"),
                rmse=9999,
                'ndrop%'=0)
    ##Realizamos la validacion cruzada para cada fold
    tryCatch({
      for (k in 1:nfold) {
        
        # Crear conjuntos de entrenamiento y prueba para el fold k
        train_data <- training(folds[[1]][[k]])#Armando el subconjunto de datos para entrenar
        test_data <- testing(folds[[1]][[k]])#Armando el subconjunto de datos para realizar la validacion cruzada
        
        ## Entrenando el modelo
        #Obtenemos el variograma empirico y el ajustado
        if (trainwo) {
          train_v<-filter(train_data,outlier==F)
        }else{
          train_v<-train_data
        }
        v <- variogram(fm, data=train_v,cutoff=cutoff,width=width,debug.level=-1)
        rm(train_v)
        m_fit <- fit.variogram(object = v, model = m, fit.kappa = T)
        m_fit_df <- as.data.frame(m_fit)
        #Obtenemos el modelo
        model <- gstat(formula = fm, data = train_data, model= m_fit,maxdist = dist_nn*10)
        #Removemos variables que no se usan más
        rm(list = c("train_data"))
        ## Realizar predicciones en el subconjunto de prueba
        print(paste("Layer:",layer,"- Fold:",k))
        if(data.class(test_data$geom)=="sfc_MULTIPOINT"){
          predictions <- predict(model, test_data |> st_cast("POINT"),debug.level=-1)
        }
        else{
          predictions <- predict(model, test_data,debug.level=-1)
        }
        names(predictions)[1]<-"pred"
        #Armando un dataframe con el predicho vs observado y su error residual, se agrega tambien el fold utilizado
        x<-cbind(pred = predictions$pred, 
                 observed = test_data$data_col, 
                 residual = test_data$data_col-predictions$pred,
                 geom = predictions$geometry |> st_as_sf(), 
                 fold = k)
        y <- cbind(fold = k,
                   fit_psill = m_fit_df$psill[2],
                   fit_range = m_fit_df$range[2],
                   fit_nugget = m_fit_df$psill[1],
                   fit_kappa = m_fit_df$kappa[2])
        #Agregamos lo calculado al df de retorno
        ret <- rbind(ret,x)
        fit <- rbind(fit,y)
        #Sacamos de la memoria las variables que no se usan
        rm(list = c("model", "test_data"))
      }
      #Finaliza el contador de tiempo
      time_end<-Sys.time()
      xv <-tibble(layer=layer,formula=fm_str,data=list(ret),fit_param=list(fit),time=difftime(time_end, time_start, units = "secs")) |>
        group_by(layer, formula) |>
        mutate(rmse = 
                 sqrt(
                   mean(
                     filter(ret,!is.na(residual))$residual^2)),
               'ndrop%'=(100*length(filter(ret,!is.na(residual))$residual)/length(ret$residual)) %>% round(digits=1))
    }
    ,error=function(e){
      #Si tira un error abre el log, registra la capa y el error, luego lo cierra
      sink(file.path(log_path,log_name),append = T)
      print("--------------------------------")
      print(paste("Error en:",layer))
      print(e)
      sink()
    }
    )
    ok <- rbind(ok,xv)
  }
  return(ok)
}

#' @title Obtener los Mejores Modelos de Predicción por OK con Validación Cruzada utilizando automap
#' @description Esta función calcula y devuelve los mejores modelos de predicción por el método ordinary kriging (OK) para cada capa utilizando validación cruzada.
#' @param path Ruta del directorio donde se encuentra el archivo GPKG.
#' @param file_name Nombre del archivo GPKG con los datos.
#' @param param Dataframe con los parámetros iniciales para correr el modelo
#' @param outlier valor lógico que define si sacar los outliers. 
#' Por defecto es TRUE y no saca los outliers
#' @param trainwo logico que determina si entrenar el modelo sin los outliers, Por defecto es FALSE
#' @param n_col Número de columna a resumir en caso de tener múltiples columnas (e.g. datos apareados). 
#' Por defecto, es 1.
#' @param seed Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto. 
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param nfold Numero de folds a utilizar en la validación cruzada
#' Por defecto es 5
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera indicar que layers utilizar. 
#' Puede ser un vector o un entero, por defecto es 0.
#' @return Un dataframe con el mejor modelo de predicción por OK y los datos de la validación cruzada para cada capa.
#' @details La función setea el archivo del log.
#' Luego carga en una grilla de tuneo los hiperparametros especificados en param.  
#' Luego recorre la grilla de tuneo, obtiene la validación cruzada, 
#' utilizando el método de k-folds y arma un dataframe donde cada fila tiene su layer, tiempo de ejecución y los resultados de la validacion cruzada y lo retorna
#' @examples get_idw_res(path="~/Documentos/Tesina/Data/clean", file_name="Elevacion.gpkg",param="~/Documentos/Tesina/Data/variograms/ECa30_vgm.RDS", seed = 701408733)
#' @export
get_ok_atmp <- function(path,
                        file_name,
                        param,
                        outlier=T,
                        trainwo=F,
                        iqr_p=3,
                        n_col=1,
                        seed=NULL,
                        nfold = 5,
                        debug = 0){
  #Paquetes
  require(sf)
  require(dplyr)
  require(tidyr)
  require(rsample)
  require(gstat)
  require(automap)
  #Crea el nombre del log de la corrida
  log_name<-"get_ok_res.txt"
  #Ruta local donde se guardan el log, si no existe la crea
  log_path <- file.path(path,"log")
  dir.create(log_path,showWarnings = F)
  
  #Ruta temporal donde se guardan las salidas parciales
  temp <- file.path(path,"OK")
  dir.create(temp,showWarnings = F)
  temp <- file.path(temp,"temp")
  dir.create(temp,showWarnings = F)
  
  #Obtiene una lista de todas las capas del gpkg
  file_path<-file.path(path,file_name)
  
  #Si debug es un vector solo trae los layers indicados en el vector
  #Si es un valor mayor a cero trae solo el layer indicado
  if(length(debug)>1){
    param <- param[debug,]
  }else if (debug > 0 && debug<=length(param$layer)) {
    param <- param[debug,]
  }
  
  #Obtenemos los resultados del ajuste con los hiperparametros seteados en param
  #Recorriendo cada uno y se los almacena en ok
  columns <- c("layer", "formula", "data", "fit_param", "time", "rmse", "ndrop%")
  ok<-as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
  imax <- nrow(param)
  for(i in 1:imax){
    #Arranca el contador de tiempo
    time_start<-Sys.time()
    #Obtenemos los datos de la capa correspondiente
    layer <- param$layer[[i]]
    fm_str <- param$formula[[i]]
    fm <- formula(fm_str)
    dist_nn <- param$dist_nn[[i]] |> as.numeric()
    
    #Cargamos el set de datos y estandarizamos el nombre de la columna que contiene los datos
    data <- st_read(file_path,layer = layer)
    data <- data %>% 
      rename(data_col=colnames(data)[n_col]) %>% 
      mutate(data_col=as.numeric(data_col),
             x = st_coordinates(data)[, "X"], #Extrae las coordenadas "x" e "y" de la columna "geom"
             y = st_coordinates(data)[, "Y"]) %>% 
      mutate(x=x-min(x),y=y-min(y)) %>% 
      st_drop_geometry() %>% st_as_sf(coords=c("x","y")) 
    data <- data %>% 
      mutate(x = st_coordinates(data)[, "X"],
             y = st_coordinates(data)[, "Y"])
    
    # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
    Q1 <- quantile(data$data_col, 0.25)
    Q3 <- quantile(data$data_col, 0.75)
    # Calcular el rango intercuartílico (IQR)
    IQR <- Q3-Q1
    # Calcular los límites inferior y superior para identificar los outliers
    limite_inferior <- Q1 - iqr_p * IQR
    limite_superior <- Q3 + iqr_p * IQR
    # Crear un campo logico que identifique si el dato es un outlier
    data <- data %>%
      mutate(outlier = if_else(
        data_col<limite_inferior|data_col>limite_superior
        ,T
        ,F))
    #Eliminando los Outliers
    if (outlier==F) {
      data<-data |> filter(outlier==F)
    }
    # Dividiendo los datos en nfolds
    if(!is.null(seed)){
      set.seed(seed)  # Semilla
    }
    folds <- vfold_cv(data, v = nfold)
    #Sacamos de memoria las variables que no se usan
    rm(list = c("data","Q1","Q3","IQR","limite_inferior","limite_superior"))
    #Realizar la validación cruzada
    
    # Inicializar una df para almacenar los resultados 
    columns <- c("pred","observed","residual","geom","fold")
    ret <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
    columns <- c("fold","fit_psill","fit_range","fit_nugget","fit_kappa")
    fit <- as_tibble(matrix(nrow = 0, ncol = length(columns)), .name_repair = ~ columns)
    #Inicializamos xv con valores de error para que pueda hacer el rbind de todas formas
    time_end<-Sys.time()
    xv <-tibble(layer=layer,
                formula="Error",
                data=list(data.frame()),
                fit_param=list(data.frame()),
                time=difftime(time_end, time_start, units = "secs"),
                rmse=9999,
                'ndrop%'=0)
    ##Realizamos la validacion cruzada para cada fold
    tryCatch({
      for (k in 1:nfold) {
        
        # Crear conjuntos de entrenamiento y prueba para el fold k
        train_data <- training(folds[[1]][[k]])#Armando el subconjunto de datos para entrenar
        test_data <- testing(folds[[1]][[k]])#Armando el subconjunto de datos para realizar la validacion cruzada
        
        ## Entrenando el modelo
        #Obtenemos el variograma empirico y el ajustado
        if (trainwo) {
          train_v<-filter(train_data,outlier==F)
        }else{
          train_v<-train_data
        }
        m_fit <- automap::autofitVariogram(fm,train_v)$var_model
        m_fit_df <- as.data.frame(m_fit)
        #Obtenemos el modelo
        model <- gstat(formula = fm, data = train_data, model= m_fit,maxdist = dist_nn*10)
        saveRDS
        #Removemos variables que no se usan más
        rm(list = c("train_data"))
        ## Realizar predicciones en el subconjunto de prueba
        print(paste("Layer:",layer,"- Fold:",k))
        if(data.class(test_data$geom)=="sfc_MULTIPOINT"){
          predictions <- predict(model, test_data |> st_cast("POINT"),debug.level=-1)
        }
        else{
          predictions <- predict(model, test_data,debug.level=-1)
        }
        names(predictions)[1]<-"pred"
        #Armando un dataframe con el predicho vs observado y su error residual, se agrega tambien el fold utilizado
        x<-cbind(pred = predictions$pred, 
                 observed = test_data$data_col, 
                 residual = test_data$data_col-predictions$pred,
                 geom = predictions$geometry |> st_as_sf(), 
                 fold = k)
        y <- cbind(fold = k,
                   fit_psill = m_fit_df$psill[2],
                   fit_range = m_fit_df$range[2],
                   fit_nugget = m_fit_df$psill[1],
                   fit_kappa = m_fit_df$kappa[2])
        #Agregamos lo calculado al df de retorno
        ret <- rbind(ret,x)
        fit <- rbind(fit,y)
        #Sacamos de la memoria las variables que no se usan
        rm(list = c("model", "test_data"))
      }
      ##bind----
      #Finaliza el contador de tiempo
      time_end<-Sys.time()
      xv <-tibble(layer=layer,formula=fm_str,data=list(ret),fit_param=list(fit),time=difftime(time_end, time_start, units = "secs")) |>
        group_by(layer, formula) |>
        mutate(rmse = 
                 sqrt(
                   mean(
                     filter(ret,!is.na(residual))$residual^2)),
               'ndrop%'=(100*length(filter(ret,!is.na(residual))$residual)/length(ret$residual)) %>% round(digits=1))
    }
    ,error=function(e){
      #Si tira un error abre el log, registra la capa y el error, luego lo cierra
      sink(file.path(log_path,log_name),append = T)
      print("--------------------------------")
      print(paste("Error en:",layer))
      print(e)
      sink()
    }
    )
    ok <- rbind(ok,xv)
  }
  return(ok)
}


#' @title Explorar variogramas
#' @description Esta función obtiene la tendencia que mejor R2 tenga grafica el variograma y permite setear los hiperparametros para utilizar en OK
#' @param file_path Ruta del archivo GPKG con los datos.
#' @param outlier valor lógico que define si sacar los outliers.  
#' Por defecto es TRUE y no saca los outliers
#' @param iqr_p Si outlier es FALSE especifica que outliers filtrar. 
#' Por defecto es 3
#' @param n_col Número de columna a resumir en caso de tener múltiples columnas (e.g. datos apareados). 
#' Por defecto, es 1.
#' @param debug Este parámetro se utiliza en un entorno de prueba para cuando se quiera indicar que layers utilizar. 
#' Puede ser un vector o un entero, por defecto es 0.
#' @return Un dataframe con el con la tendencia y los valores seteados para cada capa.
#' @details Primero carga las capas de los datos indicados en la ruta. 
#' Luego recorre capa por capa cargando los datos, obtiene la tendencia con mejor R2 y calcula el variograma. 
#' Luego grafica el variograma, permite ingresar datos de sill, nugget y range, los grafica y pregunta si se quiere quedar con esos valores o cargar unos nuevos.
#' Por ultimo arma y retorna un dataframe con los valores seteados para cada capa
#' @examples get_vgm_param(file_path="~/Documentos/Tesina/Data/clean/Elevacion.gpkg",debug = 1)
#' @export
get_vgm_param <- function(file_path,
                          outlier=T,
                          iqr_p=3, 
                          n_col=1,
                          debug=0){
  #Paquetes
  require(dplyr)
  require(sf)
  require(ggplot2)
  require(gstat)
  require(doParallel)
  
  #Seteamos nucleos
  cl <- makeCluster(4)
  registerDoParallel(cl)
  mcaffinity(1:4)
  on.exit(stopCluster(cl))
  #Capas
  layers <- st_layers(file_path)
  
  #Si debug es un vector solo trae los layers indicados en el vector
  #Si es un valor mayor a cero trae solo el layer indicado
  if(length(debug)>1){
    layers$name <- layers$name[debug]
  }else if (debug > 0 && debug<=length(layers$name)) {
    layers$name <- c(layers$name[debug])
  }
  # Inicializar una df para almacenar los resultados 
  columns <- c("layer","formula","sill","nugget","range","variogram","r_squared","cutoff","width","plot")
  ret <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(ret) <- columns
  for (layer in layers$name) {
    #Datos
    data <- st_read(file_path,layer = layer)
    data <- data |> 
      rename(data_col=colnames(data)[n_col]) |> 
      mutate(data_col=as.numeric(data_col),
             x = st_coordinates(data)[, "X"], #Extrae las coordenadas "x" e "y" de la columna "geom"
             y = st_coordinates(data)[, "Y"]) |> 
      mutate(x=x-min(x),y=y-min(y)) %>% 
      st_drop_geometry() %>% st_as_sf(coords=c("x","y")) 
    data <- data %>% 
      mutate(x = st_coordinates(data)[, "X"],
             y = st_coordinates(data)[, "Y"])
    
    # Calcular el primer cuartil (Q1) y el tercer cuartil (Q3)
    Q1 <- quantile(data$data_col, 0.25)
    Q3 <- quantile(data$data_col, 0.75)
    # Calcular el rango intercuartílico (IQR)
    IQR <- Q3-Q1
    # Calcular los límites inferior y superior para identificar los outliers
    limite_inferior <- Q1 - iqr_p * IQR
    limite_superior <- Q3 + iqr_p * IQR
    # Crear un campo logico que identifique si el dato es un outlier
    data <- data %>%
      mutate(outlier = if_else(data_col<limite_inferior|data_col>limite_superior,T,F))
    #Eliminando los Outliers
    if (outlier==F) {
      data<-data |> filter(outlier==F)
    }
    
    #Cutoff y width por defecto
    ##Calcula la longitud de la diagonal de la caja que abarca los datos
    diagonal_length <- sqrt(diff(range(data$x))^2 + diff(range(data$y))^2)
    cutoff <- diagonal_length / 3
    width <- cutoff/15
    #Obtenemos la tendencia con mejor R^2
    fm_rs <- data.frame(fm = c("data_col ~ 1","data_col ~ x+y","data_col ~ x + y + I(x^2) + I(y^2) + x:y"),
                        adj_r_squared = c(NA,NA,NA))
    for (i in 1:3) {
      m <-summary(lm(formula(fm_rs$fm[i]),data))  
      fm_rs$adj_r_squared[i] <- m$adj.r.squared
    }
    m2<-stats::step(lm(data_col~x+y+I(x^2)+I(y^2)+x:y,data))
    m <- summary(m2)
    m_ch<-as.character(formula(m2))
    fm_rs <- fm_rs %>% add_row(fm=paste(m_ch[2],m_ch[1],m_ch[3]),adj_r_squared=m$adj.r.squared)
    print(fm_rs)
    
    vl<-foreach(i = 1:4,.packages = c("gstat"))%dopar%{
      v <-variogram(formula(fm_rs$fm[i]), data=data,cutoff=diagonal_length,width=diagonal_length/15)
      v
    }
    v1 <- vl[[1]]
    v2 <- vl[[2]]
    v3 <- vl[[3]]
    v4 <- vl[[4]]
    print(
      ggplot()+
        geom_line(mapping = aes(x=dist,y = gamma),data = v1,color="green",show.legend = T)+
        geom_line(mapping = aes(x=dist,y = gamma),data = v2,color="purple",show.legend = T)+
        geom_line(mapping = aes(x=dist,y = gamma),data = v3,color="blue",show.legend = T)+
        geom_line(mapping = aes(x=dist,y = gamma),data = v4,color="red",show.legend = T)+
        labs(title = layer, x = "Distancia (m)", y = "Semivariograma",colour="Variogramas") +
        annotate("text", x = max(v1$dist)*0.16, y = max(v1$gamma)*0.95, label = paste("1- Constante:",fm_rs$fm[1]), color = "green", size = 3,hjust=1) +
        annotate("text", x = max(v1$dist)*0.16, y = max(v1$gamma)*0.90, label = paste("2- Linear:",fm_rs$fm[2]), color = "purple", size = 3,hjust=1) +
        annotate("text", x = max(v1$dist)*0.25, y = max(v1$gamma)*0.85, label = paste("3- Cuadratico:",fm_rs$fm[3]), color = "blue", size = 3,hjust=1) +
        annotate("text", x = max(v1$dist)*0.16, y = max(v1$gamma)*0.78, label = paste("4- Step:",fm_rs$fm[4]), color = "red", size = 3,hjust=1) +
        coord_cartesian(xlim = c(0, NA), ylim = c(0,NA)) #+
      # theme_classic()
    )
    opt<-0
    # opt<-1
    while (opt!=1&opt!=2&opt!=3&opt!=4) {
      opt<-readline(prompt = "Seleccione la tendencia a utilizar: ")
    }
    opt<-as.numeric(opt)
    fm_<-fm_rs[opt,]
    
    fm <- formula(fm_$fm[1])
    print(paste("Formula:",fm_$fm[1]))
    print(paste0("Cutoff: ",cutoff," - Width: ", width))
    # #Obtenemos el variograma
    cont <- "" #esta variable se utiliza para controlar cuando salir de los whiles
    while(cont != "N"){
      cont<-""
      #Pedimos por consola los valores
      cutoff <- readline(prompt = "Ingresa un valor de cutoff: ") |> as.numeric()
      width <- readline(prompt = "Ingresa un valor de width: ") |> as.numeric()
      if (cutoff<0 && cutoff>=-4) {
        fm_<-fm_rs[-cutoff,]
        fm <- formula(fm_$fm[1])
      }else { 
        if (is.na(cutoff) && is.na(width)) {
          cutoff <- diagonal_length / 3
          width <- cutoff/15
          v <- variogram(fm, data=data,cutoff=cutoff,width=width,debug.level=-1)
        } else if (!is.na(cutoff) && !is.na(width)) {
          #Obtenemos el variograma
          v <- variogram(fm, data=data,cutoff=cutoff,width=width,debug.level=-1)
        } else if (!is.na(cutoff)) {
          width <- cutoff/15
          #Obtenemos el variograma
          v <- variogram(fm, data=data,cutoff=cutoff,width=width,debug.level=-1)
        } else {
          cutoff <- diagonal_length / 3
          #Obtenemos el variograma
          v <- variogram(fm, data=data,cutoff=cutoff,width=width,debug.level=-1)
        }
        print(paste0("Cutoff: ",cutoff," - Width: ", width))
        print(v)
        #Graficamos el variograma con los valores cargados
        maxv <- max(v$gamma) %>% round(1)
        minv <- min(v$gamma) %>% round(1)
        print(ggplot(v)+
                aes(x = dist, y = gamma)+
                geom_line(color = "green") +
                geom_smooth(se = F, linetype = 3,fullrange=T) + 
                geom_hline(yintercept = c(maxv,minv), linetype = "dashed", color = c("red","black")) +  
                geom_vline(xintercept = 0, linetype = "solid", color = "black") +   
                labs(title = layer, x = "Distancia (m)", y = "Semivariograma") +
                annotate("text", x = max(v$dist)*0.15, y = maxv*0.95, label = paste("max:",maxv), color = "red", size = 3,hjust=1) +
                annotate("text", x = max(v$dist)*0.15, y = minv*0.95, label = paste("min:",minv), color = "black", size = 3,hjust=1) +
                coord_cartesian(xlim = c(0, NA), ylim = c(0,NA)) + 
                scale_x_continuous(breaks = round(seq(0, max(v$dist), by = 25),1))
        )
      }
      #Da la opcion de quedarse con los datos cargados o cargar unos nuevos
      #Si no responde Y/N vuelve a pedirlo
      while (cont !="N" && cont !="Y") {
        cont <- readline(prompt = "Desea cambiar los datos cargados? (Y/N): ") |> toupper()
      }
    }
    #Mientras cont sea distinto de N vuelve a pide los valores de sill, nugget y range
    cont <- ""
    while(cont != "N"){
      cont<-""
      #Pedimos por consola los valores
      sill <- readline(prompt = "Ingresa un valor de sill: ") |> as.numeric()
      nugget <- readline(prompt = "Ingresa un valor de nugget: ") |> as.numeric()
      range <- readline(prompt = "Ingresa un valor range: ") |> as.numeric()
      #Graficamos el variograma intersectado con los valores cargados
      p <- ggplot(v)+
        aes(x = dist, y = gamma)+
        geom_line(color = "green") +
        geom_smooth(se = F, linetype = 3) + 
        geom_hline(yintercept = c(sill,nugget), linetype = "dashed", color = c("red","black")) +  
        geom_vline(xintercept = c(0,range), linetype = c("solid","dashed"), color = c("black","blue")) +
        labs(title = layer, x = "Distancia (m)", y = "Semivariograma") +
        annotate("text", x = max(v$dist)*0.15, y = sill*0.95, label = paste("sill:",sill), color = "red", size = 3) +
        annotate("text", x = max(v$dist)*0.15, y = nugget*0.95, label = paste("nugget:",nugget), color = "black", size = 3) +
        annotate("text", x = range, y = nugget*0.95, label = paste("range:",range), color = "blue", size = 3) +
        coord_cartesian(xlim = c(0, NA), ylim = c(0,NA)) +
        scale_x_continuous(breaks = round(seq(0, max(v$dist), by = 25),1))
      print(p)
      #Da la opcion de quedarse con los datos cargados o cargar unos nuevos
      #Si no responde Y/N vuelve a pedirlo
      while (cont !="N" && cont !="Y") {
        cont <- readline(prompt = "Desea cambiar los datos cargados? (Y/N): ") |> toupper()
      }
    }
    #Retorna la tendencia utilizada, sill, nugget, range, variograma,r2,cutoff,width y plot
    retb<-tibble(layer=layer,formula=fm_$fm[1],sill=sill,nugget=nugget,range=range,variogram=list(v),r_squared=fm_$adj_r_squared[1],cutoff=cutoff,width=width,plot=list(p))
    ret <- ret |> rbind(retb)
    saveRDS(retb,file.path("Data/variograms/temp",paste0(layer,"_vgm.RDS")))
  }
  registerDoSEQ()
  return(ret)
}