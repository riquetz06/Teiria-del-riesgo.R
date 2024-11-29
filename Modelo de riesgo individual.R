#Cargar las librerias necesarias
library(ggplot2)
library(dplyr)
#Generar datos simulados
set.seed(123)
n<-1000
x<-rnorm(n, mean=50, sd=10) #Simulamos reclamaciones
p<-0.05 #Probabilidad de reclamación
y<-rbinom(n, size = 1,prob=p)#Indicador de reclamación (0=no reclamación, 1=reclamación)
#Crear un data frame
data<-data.frame(Reclamacion = y,Monto = x)
#Calcular el costo total de las reclamaciones
costo_total <- sum(data$Monto*data$Reclamacion)
#Calcular el costo esperado por reclamación
costo_esperado <- mean(data$Monto*data$Reclamacion)
#Visualizar los datos
ggplot(data, aes(x=Monto))+geom_histogram(bins=30,fill="blue",alpha=0.7)+
  labs(title = "Distribución de Montos de Reclamaciones", x = "Monto de Reclamación", y = "Frecuencia")
#Mostrar resultados
cat("Costo total de las reclamaciones:", costo_total,"\n")
cat("Costo esperado por reclamación:",costo_esperado,"\n")

#Utilizando una base de datos:
#Cargar base de datos de reclamaciones
datos_reclamos <- read.csv(file.choose())
datos_reclamos
x<-rnorm(n, mean=138322, sd=213283) #datos de reclamaciones
p<-0.05 #Probabilidad de reclamación
y<-rbinom(n, size = 1,prob=p)#Indicador de reclamación (0=no reclamación, 1=reclamación)
#Crear un data frame
data<-data.frame(Reclamacion = y,Monto = x)
#Calcular el costo total de las reclamaciones
costo_total <- sum(data$Monto*data$Reclamacion)
#Calcular el costo esperado por reclamación
costo_esperado <- mean(data$Monto*data$Reclamacion)
#Visualizar los datos
ggplot(data, aes(x=Monto))+geom_histogram(bins=30,fill="blue",alpha=0.7)+
  labs(title = "Distribución de Montos de Reclamaciones", x = "Monto de Reclamación", y = "Frecuencia")
#Mostrar resultados
cat("Costo total de las reclamaciones:", costo_total,"\n")
cat("Costo esperado por reclamación:",costo_esperado,"\n")

