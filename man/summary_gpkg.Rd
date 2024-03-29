% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caracterizacion.R
\name{summary_gpkg}
\alias{summary_gpkg}
\title{Medidas de resumen}
\usage{
summary_gpkg(
  path,
  data_path,
  index_path = "index.RDS",
  seed = NULL,
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

\item{index_path}{Nombre del archivo del indice
Por defecto es "index.RDS"}

\item{seed}{Valor de la semilla a utilizar, si no se especifica utiliza la que queda por defecto.}

\item{outlier}{valor lógico que define si sacar los outliers.
Por defecto es TRUE y no saca los outliers}

\item{iqr_p}{Si outlier es FALSE especifica que outliers filtrar.
Por defecto es 3}

\item{freecor}{Cantidad de núcleos de CPU a dejar libres al ejecutar la función.
Por defecto es 1}

\item{debug}{Este parámetro se utiliza en un entorno de prueba para cuando se quiera achicar el numero de capas a este valor.
Por defecto es 0}

\item{n_col}{Espera un entero con el numero de la columna a resumir en caso de que tenga más de una (Datos apareados).
Por defecto es 1}
}
\value{
Retorna un dataframe con las medidas de resumen de cada capa
}
\description{
Función que devuelve las medidas de resumen de los datos de un gpkg
}
\details{
Primero setea los parámetros para realizar el calculo en múltiples núcleos.
Luego lista y recorre las capas para obtener la medidas de resumen.
Por ultimo retorna el dataframe con las medidas de cada capa.
}
\examples{
summary_gpkg("D:/Facultad/Tesina/Temp/clean","Elevacion.gpkg",freecor = 2)
}
