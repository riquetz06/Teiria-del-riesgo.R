# Crear un data frame simulado
set.seed(123)  # Para reproducibilidad

auto_data <- data.frame(
  age = sample(18:80, 100, replace = TRUE),                 # Edad del conductor
  accidents = sample(0:5, 100, replace = TRUE),             # Número de accidentes
  car_type = sample(c("Sedan", "SUV", "Truck", "Coupe"), 100, replace = TRUE),  # Tipo de auto
  region = sample(c("North", "South", "East", "West"), 100, replace = TRUE),   # Región
  claim_amount = round(rnorm(100, mean = 1000, sd = 500), 2)  # Monto de la reclamación
)

# Mostrar las primeras filas del data frame
head(auto_data)
print(auto_data)# Mostrar todas las filas del data frame

# Convertir variables categóricas a factores
auto_data$age <- as.factor(auto_data$age)
auto_data$car_type <- as.factor(auto_data$car_type)
auto_data$region <- as.factor(auto_data$region)

# Verificar si hay valores negativos
sum(auto_data$claim_amount < 0)

# Eliminar o corregir valores negativos
auto_data <- auto_data[auto_data$claim_amount >= 0, ]

# Crear el modelo GLM utilizando la librería 'stats'
modelo_glm <- glm(claim_amount ~ age + accidents + car_type + region, 
                  family = poisson(link = "log"), data = auto_data)

# Resumen del modelo
summary(modelo_glm)

# Validar el modelo con predicciones
predicciones <- predict(modelo_glm, type = "response")

# Comparar predicciones con valores reales
comparacion <- data.frame(Real = auto_data$claim_amount, Prediccion = predicciones)

# Visualización de los resultados
plot(auto_data$claim)

# Crear el modelo GLM utilizando la familia Gamma
modelo_glm <- glm(claim_amount ~ age + accidents + car_type + region, 
                  family = Gamma(link = "log"), data = auto_data)

# Resumen del modelo
summary(modelo_glm)

# Validar el modelo con predicciones
predicciones <- predict(modelo_glm, type = "response")

# Comparar predicciones con valores reales
comparacion <- data.frame(Real = auto_data$claim_amount, Prediccion = predicciones)

# Visualización de los resultados
plot(auto_data$claim_amount, predicciones, main="Comparación de Predicciones vs. Valores Reales", 
     xlab="Valores Reales", ylab="Predicciones", pch=19, col="blue")
abline(0, 1, col="red")

#Exportar de R a Excel 
# Instalar y cargar la biblioteca necesaria
install.packages("writexl")
library(writexl)
# Exportar el data frame a un archivo Excel
write_xlsx(auto_data, "ruta/donde/deseas/guardar/auto_data.xlsx")

# Residuos del modelo
residuos <- residuals(modelo_glm, type = "deviance")

# Gráfico Q-Q de los residuos
qqnorm(residuos)
qqline(residuos, col = "red")

# Valores ajustados
valores_ajustados <- fitted(modelo_glm)

# Residuos estandarizados
residuos_estandarizados <- rstandard(modelo_glm)

# Gráfico de residuos estandarizados vs valores ajustados
plot(valores_ajustados, residuos_estandarizados, 
     main = "Residuos Estandarizados vs Valores Ajustados", 
     xlab = "Valores Ajustados", ylab = "Residuos Estandarizados")
abline(h = 0, col = "red")

install.packages("lmtest")
library(lmtest)

# Test de Breusch-Pagan
bptest(modelo_glm)

#Tarificación
# Crear un data frame simulado
set.seed(123)  # Para reproducibilidad

auto_data <- data.frame(
  age = sample(18:80, 100, replace = TRUE),                 # Edad del conductor
  accidents = sample(0:5, 100, replace = TRUE),             # Número de accidentes
  car_type = sample(c("Sedan", "SUV", "Truck", "Coupe"), 100, replace = TRUE),  # Tipo de auto
  region = sample(c("North", "South", "East", "West"), 100, replace = TRUE),   # Región
  claim_amount = round(rnorm(100, mean = 1000, sd = 500), 2)  # Monto de la reclamación
)

# Asegurarse de que claim_amount no tiene valores no positivos
auto_data <- auto_data[auto_data$claim_amount > 0, ]

# Convertir variables categóricas a factores
auto_data$age <- as.factor(auto_data$age)
auto_data$car_type <- as.factor(auto_data$car_type)
auto_data$region <- as.factor(auto_data$region)

# Crear el modelo GLM utilizando la familia Gamma
modelo_glm <- glm(claim_amount ~ age + accidents + car_type + region, 
                  family = Gamma(link = "log"), data = auto_data)

# Predecir el monto de la reclamación para nuevos clientes
nuevos_clientes <- data.frame(
  age = sample(18:80, 10, replace = TRUE),                 
  accidents = sample(0:5, 10, replace = TRUE),             
  car_type = sample(c("Sedan", "SUV", "Truck", "Coupe"), 10, replace = TRUE),  
  region = sample(c("North", "South", "East", "West"), 10, replace = TRUE))

# Ajustar niveles de factores para que coincidan con los datos originales
nuevos_clientes$age <- factor(nuevos_clientes$age, levels = levels(auto_data$age))
nuevos_clientes$car_type <- factor(nuevos_clientes$car_type, levels = levels(auto_data$car_type))
nuevos_clientes$region <- factor(nuevos_clientes$region, levels = levels(auto_data$region))

# Predecir el monto de la reclamación utilizando el modelo GLM
predicciones <- predict(modelo_glm, nuevos_clientes, type = "response")
nuevos_clientes$predicted_claim <- predicciones

# Añadir una tarifa fija
tarifa_fija <- 100  
nuevos_clientes$premium <- nuevos_clientes$predicted_claim + tarifa_fija

# Instalar y cargar la biblioteca writexl
install.packages("writexl")
library(writexl)

# Exportar el data frame a un archivo Excel
write_xlsx(nuevos_clientes, "ruta/donde/deseas/guardar/nuevos_clientes.xlsx")

install.packages("car")
library(car)

# Cálculo del VIF
vif(modelo_glm)

#Resultados
print(nuevos_clientes)
head(nuevos_clientes)
tail(nuevos_clientes)
# Instalar y cargar la biblioteca writexl
install.packages("writexl")
library(writexl)

# Exportar el data frame a un archivo Excel
write_xlsx(nuevos_clientes, "ruta/donde/deseas/guardar/nuevos_clientes.xlsx")

# Gráfico de predicciones vs. valores reales
plot(nuevos_clientes$predicted_claim, nuevos_clientes$premium, 
     main="Predicciones vs. Primas Calculadas", 
     xlab="Predicciones del Monto de la Reclamación", ylab="Primas Calculadas", 
     pch=19, col="blue")
abline(0, 1, col="red")
#Resumen
summary(nuevos_clientes)

# Predecir el monto de la reclamación para nuevos clientes utilizando el modelo GLM
predicciones_futuras <- predict(modelo_glm, nuevos_clientes, type = "response")

# Calcular la reserva como la suma de las predicciones futuras
reserva_estimada <- sum(predicciones_futuras)

# Predecir el monto de la reclamación para nuevos clientes
predicciones_futuras <- predict(modelo_glm, nuevos_clientes, type = "response")

# Verificar predicciones
print(predicciones_futuras)

# Manejo de valores faltantes
nuevos_clientes <- na.omit(nuevos_clientes)

# Volver a predecir después de manejar NA
predicciones_futuras <- predict(modelo_glm, nuevos_clientes, type = "response")

# Verificar las nuevas predicciones
print(predicciones_futuras)

# Calcular la reserva como la suma de las predicciones futuras
reserva_estimada <- sum(predicciones_futuras)

# Mostrar la reserva estimada
print(reserva_estimada)

#Reservas con Método Chain-Ladder
# Eliminar filas con valores NA
auto_data <- na.omit(auto_data)

# Crear un ejemplo de triángulo de desarrollo de reclamaciones sin valores NA o Inf
claims_triangle <- matrix(c(
  150, 200, 250, 300, 350, 
  160, 210, 260, 310, 360, 
  170, 220, 270, 320, 370, 
  180, 230, 280, 330, 380, 
  190, 240, 290, 340, 390), 
  nrow = 5, byrow = TRUE)

# Convertir la matriz a un triángulo de desarrollo
claims_triangle <- as.triangle(claims_triangle)
print(claims_triangle)

# Aplicar el método Chain-Ladder
chain_ladder_model <- MackChainLadder(claims_triangle)
print(chain_ladder_model)

# Calcular las reservas
reserves <- summary(chain_ladder_model)$Totals["IBNR",]

# Mostrar las reservas calculadas
print(reserves)

# Instalar y cargar la biblioteca writexl
install.packages("writexl")
library(writexl)

# Crear un data frame con las reservas estimadas
reservas_df <- data.frame(Año = 1:5,  Reservas = reserves)

# Exportar el data frame a un archivo Excel
write_xlsx(reservas_df, "ruta/donde/deseas/guardar/reservas.xlsx")

# Mostrar las primeras filas del data frame resultante
head(reservas_df)
