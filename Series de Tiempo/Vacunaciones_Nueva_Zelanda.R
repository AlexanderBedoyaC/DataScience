library(TSA)
library(tseries)
library(forecast)
library(car)
library(carData)
library(lmtest)

# Modelos EstacionalesAutorregresivos Integrados de Medias Móviles
#New_Vaccinations - New Zealand

#Carga de datos
SERIE_COMPLETA<-ts(c(scan()))

#Datos para validación cruzada (90-10)
n <- length(SERIE_COMPLETA)
n_train<-round(n*0.9)
m<-n_train+1

SERIE_TRAIN<- ts(SERIE_COMPLETA[1:n_train])
SERIE_TEST<- ts(SERIE_COMPLETA[m:n])

n_test<-length(SERIE_TEST)


#################################################
################ IDENTIFICACIÓN #################
#################################################

#Validación Gráfica
par(mfrow=c(1,1))
plot(SERIE_TRAIN)
acf(SERIE_TRAIN, lag.max = n_train/4)

#Si se identifica componente estacional
#se debe especificar S (longitud estacional) a los datos
ls<-7

SERIE_TRAIN<-ts(SERIE_TRAIN, frequency = ls)
SERIE_TEST<-ts(SERIE_TEST, frequency = ls)


#Transformaciones de potencia para problemas de heterosedasticidad
Transf=powerTransform(SERIE_TRAIN)
Transf$lambda
SERIE_BC<-SERIE_TRAIN*exp(Transf$lambda)

SERIE_LOG<-log(SERIE_TRAIN)
SERIE_RAIZ<-sqrt(SERIE_TRAIN)
SERIE_2<-SERIE_TRAIN*exp(2)

par(mfrow=c(1,1))
plot(SERIE_TRAIN)
plot(SERIE_BC)
plot(SERIE_LOG)
plot(SERIE_RAIZ)
plot(SERIE_2)

#Se escoge la serie con logaritmo
#SerieTrain vs log(SerieTrain)
par(mfrow=c(1,2))
plot(SERIE_TRAIN)
plot(SERIE_LOG)

#ACF de la serie con Logaritmo
par(mfrow=c(1,1))
acf(SERIE_LOG, lag.max = n_train/4)

#Diferenciacion para problemas de no estacionariedad en media
#Diferencia Regular

#Test Estadístico
adf.test(SERIE_LOG)

#Diferenciar la serie con tendencia
SERIE_diff<-diff(SERIE_LOG)

#Validación Gráfica
par(mfrow=c(1,1))
plot(SERIE_diff)
acf(SERIE_diff, lag.max = n_train)

#Test Estadístico
adf.test(SERIE_diff)


#Diferenciacion para problemas de estacionalidad
#Diferencia Estacional

nsdiffs(SERIE_diff)

#Diferenciar la serie con estacionalidad
SERIE_DIFF<-diff(SERIE_diff,ls)

nsdiffs(SERIE_DIFF)

#Validación Gráfica
par(mfrow=c(1,1))
plot(SERIE_DIFF)
acf(SERIE_DIFF, lag.max = n_train)
length(SERIE_DIFF) 

#IDENTIFICACIÓN DE MODELOS

#Identificar ordenes de P y Q (Estacional)
par(mfrow=c(1,1))
acf(SERIE_DIFF,lag.max = n_train/2)
abline(v=c(1,2,3,4,5,6,7,8),lty=2,col="lightgray")
Pacf(SERIE_DIFF,lag.max = n_train/2)
abline(v=c(7,14,21,28,35,42,49,56,63,70,77,84,91,98),lty=2,,col="lightgray")


#Identificar ordenes de p y q
par(mfrow=c(1,1))
acf(SERIE_DIFF,lag.max = ls)
Pacf(SERIE_DIFF,lag.max = ls)


#Modelos Identificados
#ARIMA(1,1,1)(0,1,1)
#ARIMA(1,1,0)(0,1,1)


#OJO, para la función auto.arima se utiliza la serie original (sin transformaciones)
auto.arima(SERIE_TRAIN)
#ARIMA(1,0,2)(0,1,1)


#################################################
################### ESTIMACIÓN ##################
#################################################
# K es p+q+P+Q

Modelo_1<-Arima(SERIE_TRAIN,order = c(1,1,1), seasonal = c(0,1,1),lambda = 1)
summary(Modelo_1)

Modelo_2<-Arima(SERIE_TRAIN,order = c(1,1,0), seasonal = c(0,1,1),lambda = 1)
summary(Modelo_2)

Modelo_3<-Arima(SERIE_TRAIN,order = c(1,0,2), seasonal = c(0,1,1),lambda = 1)
summary(Modelo_3)


#SIGNIFICANCIA PARÁMETROS DEL MODELO 3

k<-4

Modelo_3$coef
Modelo_3$var.coef
sqrt(diag(Modelo_3$var.coef))

estadistico3=cbind(Modelo_3$coef,sqrt(diag(Modelo_3$var.coef)))
T3=estadistico3[,1]/estadistico3[,2]
T3
qt(0.975,n-k)

valor_p3<-2*pt(-abs(T3),df=n-k)
valor_p3

LI <- estadistico3[,1]-(qt(0.975,n-k)*estadistico3[,2])
LS <- estadistico3[,1]+(qt(0.975,n-k)*estadistico3[,2])
LI
LS

#SIGNIFICANCIA PARÁMETROS DEL MODELO 1

k<-3

Modelo_1$coef
Modelo_1$var.coef
sqrt(diag(Modelo_1$var.coef))

estadistico1=cbind(Modelo_1$coef,sqrt(diag(Modelo_1$var.coef)))
T1=estadistico1[,1]/estadistico1[,2]
T1
qt(0.975,n-k)

valor_p<-2*pt(-abs(T1),df=n-k)
valor_p

LI <- estadistico1[,1]-(qt(0.975,n-k)*estadistico1[,2])
LS <- estadistico1[,1]+(qt(0.975,n-k)*estadistico1[,2])
LI
LS

#SIGNIFICANCIA PARÁMETROS DEL MODELO 2

k<-2

Modelo_2$coef
Modelo_2$var.coef
sqrt(diag(Modelo_2$var.coef))

estadistico2=cbind(Modelo_2$coef,sqrt(diag(Modelo_2$var.coef)))
T2=estadistico2[,1]/estadistico2[,2]
T2
qt(0.975,n-k)

valor_p<-2*pt(-abs(T2),df=n-k)
valor_p

LI <- estadistico2[,1]-(qt(0.975,n-k)*estadistico2[,2])
LS <- estadistico2[,1]+(qt(0.975,n-k)*estadistico2[,2])
LI
LS

#################################################
################### VALIDACIÓN ##################
#################################################

#MODELO 2
#ARIMA(1,1,0)(0,1,1)

#NORMALIDAD
par(mfrow=c(1,3))
qqPlot(Modelo_2$residuals)
hist(Modelo_2$residuals)
Boxplot(Modelo_2$residuals)
shapiro.test(Modelo_2$residuals)


#VARIANZA CONSTANTE
#Modelo_2$fitted son los valores ajustados ó Y_gorro
#Modelo_2$residuals son los errores del modelo ó e_t

par(mfrow=c(1,1))
plot(as.vector(Modelo_2$fitted),as.vector(Modelo_2$residuals))
bptest(Modelo_2$fitted~Modelo_2$residuals)

#INDEPENDENCIA
par(mfrow=c(1,2))
plot(Modelo_2$residuals)
acf(Modelo_2$residuals)
Box.test(Modelo_2$residuals)


#MODELO 3
#ARIMA(1,0,2)(0,1,1)

#NORMALIDAD
par(mfrow=c(1,3))
qqPlot(Modelo_3$residuals)
hist(Modelo_3$residuals)
Boxplot(Modelo_3$residuals)
shapiro.test(Modelo_3$residuals)


#VARIANZA CONSTANTE
#Modelo_3$fitted son los valores ajustados ó Y_gorro
#Modelo_3$residuals son los errores del modelo ó e_t

par(mfrow=c(1,1))
plot(as.vector(Modelo_3$fitted),as.vector(Modelo_3$residuals))
bptest(Modelo_3$fitted~Modelo_3$residuals)

#INDEPENDENCIA
par(mfrow=c(1,2))
plot(Modelo_3$residuals)
acf(Modelo_3$residuals)
Box.test(Modelo_3$residuals)

#################################################
################### PRONÓSTICOS #################
#################################################

##ARIMA(1,1,0)(2,1,1)
pron_m1<-ts(c(Modelo_3$fitted,forecast(Modelo_3,h=17,level = 95)$mean))
pron_m1

n_validacion<-n_train+17

validacion<- ts(SERIE_COMPLETA[1:n_validacion])

plot(validacion, col=1)
lines(pron_m1,col=2) 
abline(v=195,lty=2,col="Red")

validacion_prueba<- ts(SERIE_TEST[1:17])
pron_m1_test<-ts(forecast(Modelo_3,h=17,level = 95)$mean)

plot(validacion_prueba, col=1)
lines(pron_m1_test,col=2)

pronostico_m1_tab<-forecast(Modelo_3,h=17,level = 95)
pronostico_m1_tab

#ARIMA(1,1,0)(0,1,1)
pron_m2<-ts(c(Modelo_2$fitted,forecast(Modelo_2,h=17,level = 95)$mean))
pron_m2

n_validacion<-n_train+17

validacion<- ts(SERIE_COMPLETA[1:n_validacion])

plot(validacion, col=1)
lines(pron_m2,col=2) 
abline(v=195,lty=2,col="Red")

validacion_prueba<- ts(SERIE_TEST[1:17])
pron_m2_test<-ts(forecast(Modelo_2,h=17,level = 95)$mean)

plot(validacion_prueba, col=1)
lines(pron_m2_test,col=2)

pronostico_m2_tab<-forecast(Modelo_2,h=17,level = 95)
pronostico_m2_tab
