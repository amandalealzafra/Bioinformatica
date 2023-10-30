#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 01 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl



# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #esto mide las dimensiones 
head (data)
tail (data)

# Hacemos un primer histograma para explorar los datos
hist(data)
hist(data, col="gray", main="GSE5583 - Histogram")

# Transformamos los dos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data_log=log2(data)
hist(data_log)
#La transformación logarítmica sirve para tener una transformacion mas normal y mas bonita esteticamente, en vez de tener una forma tan sesgada 

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
#col sirve para cambiar de color
#main sirve para poner un titulo a la imagen 
#las sirve para poner los ejes en vertical 

# ¿Qué es un boxplot? Un boxplot es un diagrama de cajas y bigotes que nos permite identificar valores atípicos  y comparar distribuciones.
boxplot(data_log, col=c("blue","blue", "blue", "orange", "orange", "orange"), main="GSE5583-boxplots",las=2) 

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación? Si

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado? He generado una matrix o "array" que es una tabla con datos 
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)
head(wt)
#class(wt) nos indica el tipo de datos que tenemos 
#head (wt) nos da las primeras lineas lo que podemos ver los encabezados 


# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)
head(wt.mean)
ko.mean = apply (ko, 1, mean)
head (ko.mean)

# ¿Cuál es la media más alta? La media mas alta es en los Knock-out
max(wt.mean)
max(ko.mean)

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean~wt.mean)
plot(ko.mean~wt.mean, xlab = "W", ylab = "KO", main="GSE5583 - Scatter")

# Añadir una línea diagonal con abline
abline(0, 1, col="red")
El 0 y el 1 sirven para realizar una linea vertical (y=x)

# ¿Eres capaz de añadirle un grid?
grid()
#Grid sirve para añadirle una cuadricula al grafico

# Calculamos la diferencia entre las medias de las condiciones
diff.mean=wt.mean-ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean)
hist(diff.mean, col="pink")

# Calculamos la significancia estadística con un t-test.

# Primero crea una lista vacía para guardar los p-values
#Vamos a tener tantos p values como genes tengamos 

# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? Porque son mas fiables ya que los logaritmos sirven para ver mejor las graficas 
# ¿Cuántas valores tiene cada muestra? Dentro de cada condicion tenemos 3 muestras 
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data)) { #Para cada gen
  x=wt[i,] #gen wt numero i 
  y=ko[i,] #gen ko numero i 
  
  #Hacemos el test
  t=t.test(x,y)
  
  #Añadimos el p-value a la lista
  pvalue[i] = t$p.value
  #Añadimos las estadisticas a la lista 
  tstat[i] = t$statistic
  }

head(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue) esto es el numero total de genes  
#El pvalue nos sirve para ver si la diferencia de las dos condiciones es significativa o no. En los primeros 6 genes no tienen valor significativo 

# Hacemos un histograma de los p-values.
hist(pvalue)
# ¿Qué pasa si le ponemos con una transformación de -log10?
# Cambia las escalas y refleja mas la proporcion de la distribucion, de manera, que la mayoria estan colocados en el 0
hist(-log10(pvalue), col="gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
# El volcano plot es un diagrama de expresion que representa una variable respecto a la otra: tenemos la diferencia de medias del logaritmo
# El logaritmo en base 10 negativo nos invierte la distribucion por lo que los pvalues signifcativos estarian de la mitad del grafico hacia arriba 
plot(diff.mean, -log10(pvalue), main="GSE5583-Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff=2
pvalue_cutoff=0.01
abline(v=diff.mean_cutoff, col="blue", lwd=3)
#abline(v=-diff.mean_cutoff, col="red", lwd=3)
abline(h= -log10(pvalue_cutoff), col="green", lwd=3)
#lwd es el ancho de la linea 
# Los genes significativos se encuentran en las 2 primeras cuadriculas superiores 

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold) 
filter_by_diff.mean=abs(diff.mean) >=diff.mean_cutoff
dim(data[filter_by_diff.mean, 1])
#La diferencia de medias es 2 y -2 y lo convertimos en valor absoluto 
#abs es el valor absoluto

#Hay que filtrarlos por la diferencia de medias y por el pvalue 
#Vamos a crear 2 flitros una para cada variable y los vamos a combinar 

# Ahora el filtro de p-value
filter_by_pvalue=pvalue<=pvalue_cutoff
dim(data[filter_by_diff.mean, ])
#Filtro el pvalue y extraigo los valores que sean menores o iguales a 2 
#con el dim estoy buscando cuantos genes sobrepasan el pvalue 

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios? 426 genes cumplen los criterios
filter_combined=filter_by_diff.mean & filter_by_pvalue
filtered=data [filter_combined,]
dim (filtered)
head (filtered)
#Solo me quedo con los genes que sean comunes para los dos filtros 

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main="GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue [filter_combined]), col="red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? Porque al calcular la diferencia de medias con wt-ko, hemos visto que los sobreexpresados son los knock-out ya que sus valores son mayores que los de wild type y por ello la diferencia de medias es negativa y esta en la derecha

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap? 
heatmap(filtered)
#filteredson los datos que vamos a representar 
#cexCol es el tamaño del eje x (columnas)
#Colv y Rowv son los dendogramas (columnas y filas)
#labRow es para quitar el nombre 

# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
#En nuestro heatmap: los sobreexpresados se encuentran en la parte negativa de la tabla y los reprimidos se encuentran en la parte postiva de la tabla

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep="\t", 
             quote = FALSE)
