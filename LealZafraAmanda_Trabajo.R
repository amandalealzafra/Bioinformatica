# TRABAJO R
#
# AMANDA LEAL ZAFRA
# 
getwd() #para comprobar donde estoy y cambio el directorio de trabajo a la 
# carpeta donde están los datos que voy a utilizar

# Cargo los datos
datos = read.table ("datos-trabajoR.txt", header = T)

datos # Compruebo que se han cargado bien

# Examino los datos
head (datos)  #veo las primeras 6 filas
tail (datos)  #veo las últimas 6 filas
summary (datos)  #veo un resumen de mínimos, cuartiles, medias, medianas y
                 #máximos de cada columna
dim (datos)   #veo las dimensiones de los datos: 50 filas y 3 columnas
str (datos)   #veo los valores en línea de cada  columna
# Hay 2 variables y 5 tratamientos diferentes

# Hago boxplots, uno para cada variable
boxplot (Variable1 ~ Tratamiento, data = datos, main = "Boxplot_Variable1", 
         xlab = "Tratamiento", ylab = "Variable1", col = "yellow")  

boxplot (Variable2 ~ Tratamiento, data = datos, main = "Boxplot_Variable2", 
         xlab = "Tratamiento", ylab = "Variable2", col = "blue")

# Hago un gráfico de dispersión y le pongo su leyenda
plot (datos$Variable2 ~ datos$Variable1, col = datos$Tratamiento, 
      main = "Gráfico_Dispersión", xlab = "Variable1", ylab = "Variable2")

legend(x = "bottomright", legend = c ("Tto 1", "Tto 2", "Tto 3", "Tto 4", 
      "Tto 5"), fill = c("black", "red", "light green", "light blue", "blue"), 
       title = "Tratamiento")

# Hago histogramas, uno para cada variable
hist (datos$Variable1, main = "Histograma_Variable1", xlab = "Variable1",
	ylab = "Frecuencia", col = "yellow")

hist (datos$Variable2, main = "Histograma_Variable2", xlab = "Variable2",
	ylab = "Frecuencia", col = "blue")

# Hago un factor en tratamiento
factor_tratamiento <- factor(datos$Tratamiento)
factor_tratamiento <- factor(factor_tratamiento, levels = c("1", "2", "3", "4", "5"))

# Calculo media y desviación estandar para cada tratamiento
aggregate(Variable1~factor_tratamiento,datos, mean)
aggregate(Variable2~factor_tratamiento,datos, mean)
aggregate(Variable1~factor_tratamiento,datos, sd)
aggregate(Variable2~factor_tratamiento,datos, sd)

# Miro cuantos elementos tiene cada tratamiento
table(factor_tratamiento)
# Cada tratamiento tiene 10 elementos

# Extraigo los datos de los tratamientos 1 y 4
tratamiento1 <- datos[1:10,]
tratamiento4 <- datos[31:40,]

# Averiguo si nuestra hipótesis nula es real
shapiro.test(tratamiento1$Variable1)  #para saber si es distribución normal
shapiro.test(tratamiento4$Variable1) 
#es distribución normal porque p-values > 0,05

t.test(tratamiento1$Variable1, tratamiento4$Variable1)  #para comprobar hipótesis nula
#rechazamos hipótesis nula porque p-value < 0,05, las medias no son iguales

var.test(tratamiento1$Variable1, tratamiento4$Variable1)  # para comprobar varianzas
#varianzas diferentes, p-value < 0,05
