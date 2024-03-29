% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caracterizacion.R
\name{dist_nn_gpkg}
\alias{dist_nn_gpkg}
\title{Distancia promedio al vecino mas cercano}
\usage{
dist_nn_gpkg(
  path,
  data_path,
  outlier = T,
  iqr_p = 3,
  freecor = 1,
  ncol = 1,
  debug = 0
)
}
\arguments{
\item{path}{Ruta del directorio donde se encuentra el archivo GPKG.}

\item{data_path}{Nombre del archivo GPKG con los datos.}

\item{outlier}{valor lógico que define si sacar los outliers.
Por defecto es TRUE y no saca los outliers}

\item{iqr_p}{Si outlier es FALSE especifica que outliers filtrar.
Por defecto es 3}

\item{freecor}{Cantidad de núcleos de CPU a dejar libres al ejecutar la función
Por defecto es 1}

\item{debug}{Este parámetro se utiliza en un entorno de prueba para cuando se quiera achicar el numero de capas a este valor.
Por defecto es 0}
}
\value{
Retorna un dataframe con las distancias promedio al vecino mas cercano de cada capa.
}
\description{
Función que devuelve las distancias promedio al vecino mas cercano de un gpkg
}
\details{
Primero setea los paramentros para realizar el calculo en múltiples núcleos.
Luego setea los parámetros para crear el log.
Luego lista y recorre las capas para obtener las distancias.
Por ultimo retorna el dataframe con las medidas de cada capa.
}
\examples{
dist_nn_gpkg("D:/Facultad/Tesina/Temp/clean","Elevacion.gpkg",freecor = 2)
}
