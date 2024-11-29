#RESERVAS PARA SEGUROS GENERALES

####################################
#            Chain Ladder          #
####################################

Cargar paquete Chain Ladder:
> library(ChainLadder)
> require(ChainLadder)
> data(package="ChainLadder")
> ## Sample triangle

Ejemplo 1:

> RAA
> plot(RAA)
> plot(RAA, lattice=TRUE)
> ?plot.triangle
> raa.inc <- cum2incr(RAA)
## Show first origin period and its incremental development
> raa.inc[1,]
> raa.cum <- incr2cum(raa.inc)
> ## Show first origin period and its cumulative development
> raa.cum[1,]


Ejemplo 2: (Importar triángulo de siniestros ocurridos desde Excel).
> file.choose() #indica la ruta de nuestro archivo .csv#
> myCSVfile <- "C:\\Users\\Riquet\\Downloads\\Libro1.csv" 
> tri <- read.csv(file=myCSVfile, header = FALSE)
> library(ChainLadder)
> tri <- as.triangle(as.matrix(tri))

Ejemplo 3: (DEMOS).
> install.packages("ROBDC")
> library(RODBC)
> demo(DatabaseExamples)
