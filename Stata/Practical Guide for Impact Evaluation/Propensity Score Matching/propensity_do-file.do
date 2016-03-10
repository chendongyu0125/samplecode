*------------------------*
*         PSM            *
*------------------------*

*----------------------------------------------------------------------------------------------------------------*
*Algunos puntos antes de hacer el an�lisis:

*En propensity score matching utilizamos el programa "psmatch2".
*Dado que el paquete b�sico de Stata no contiene el comando "psmatch2" es necesario instalarlo:

**** Instalaci�n "psmatch2"
ssc install psmatch2, replace

*Sobre la base de datos: Hay dos observaciones para la variable talla para la edad. ha_nchs1 indica la talla 
*para la edad del individuo en el periodo 1 (pre-tratamiento) y ha_nchs2 indica la talla para la edad
*en el segundo periodo (post-tratamiento). Todas las otras variables son observaciones del periodo 1. 

*----------------------------------------------------------------------------------------------------------------*


*------------------------------------*
*1. Probabilidad de participaci�n    *
*------------------------------------*

**** Dado que el primer paso para hacer correctamente "propensity score matching"
**** es estimar las probabilidades de participaci�n en el programa, en este do-file
**** lo primero que haremos es mostrar m�todos de estimaci�n del propensity score.

**** Primero debemos determinar qu� variables incluir en el modelo.
**** Sabemos que es importante no omitir ninguna variable ni sobreespecificar el modelo, por 
**** lo tanto se debe ser cuidadoso en el momento de decidir qu� variables incluir en el modelo.
**** En este caso seguiremos la metodolog�a propuesta en el libro, incluir variable por variable
**** y dejar �nicamente aquellas que son significativas a un nivel de significancia del 5%. 
**** Resulta intuitivo suponer que a mayor ingreso, la probabilidad de participar aumenta, veamos:

dprobit D ingresos_hogar_jefe

**** Incluimos variables que consideramos razonables:


dprobit D ingresos_hogar_jefe personas

dprobit D ingresos_hogar_jefe personas orden_n

dprobit D ingresos_hogar_jefe personas orden_n educa_jefe

dprobit D ingresos_hogar_jefe personas orden_n educa_jefe ocupado_jefe

dprobit D ingresos_hogar_jefe personas orden_n educa_jefe ocupado_jefe hombre

**** El coeficiente asociado a la variable "hombre" no es significativo. Esto resulta razonable dado que, en principio,
**** no hay raz�n para focalizar el tratamiento hacia un g�nero en particular. 

**** Entonces, nuestro modelo que tiene como objetivo determinar la probabilidad de participaci�n en el programa
**** incluye las variables ingresos del jefe del hogar, personas en el hogar, orden de nacimiento del individuo,
**** la educaci�n del jefe del hogar y una dummy que indica si el jefe de hogar est� desempleado o no. Definimos
**** un vector que incluya estas variables para simplificar el an�lisis:


global X "personas orden_n ocupado_jefe educa_jefe ingresos_hogar_jefe"

dprobit D $X

**** Utilizamos las probabilidades predichas por este modelo para generar nuestro propensity score

predict pscore


**** Veamos un histograma de las probabilidades predichas:


histogram pscore, by(D)

histogram pscore if D==1, bin(100) color(blue) addplot(kdensity pscore if D==1)

kdensity pscore if D==1, epanechnikov generate(x1 y1)

histogram pscore if D==0, bin(100) color(blue) addplot(kdensity pscore if D==0)

kdensity pscore if D==0, epanechnikov generate(x0 y0)

twoway (line y1 x1) (line y0 x0, lpattern(dash)), ytitle(Densidad) xtitle(Probabilidad de ser tratado) title(Propensity Score ) legend(order(1 "Participante=1" 2 "No participante=0"))

*-----------------*
*2. Soporte com�n *
*-----------------*


**** Gr�ficamente vemos que las probabilidades predichas son similares. Sin embargo, resulta evidente que hay probabilidades de participaci�n
**** en el grupo de tratamiento superiores a la m�xima probabilidad del grupo de control y, de la misma manera, probabilidades
**** en el grupo de control inferiores a la m�nima probabilidad del grupo de tratamiento. Para solucionar este problema imponemos el 
**** soporte com�n mediante el m�ximo y el m�nimo. 

gen pscore_sc=pscore

**** Encontramos la m�xima probabilidad predicha para el grupo de control:
sum pscore_sc if D==0
scalar max_control=r(max)

**** Encontramos la m�nima probabilidad predicha para el grupo de tratamiento:
sum pscore_sc if D==1
scalar min_D=r(min)


**** Ahora, no ser�n tenidas en cuenta probabilidades del grupo de tratamiento que superen la m�xima probabilidad del grupo de control
**** ni probabilidades del grupo de control inferiores a la m�nima probabilidad del grupo de tratamiento:
replace pscore_sc=. if D==1&pscore_sc>max_control
replace pscore_sc=. if D==0&pscore_sc<min_D

**** Vemos cuantas observaciones perdemos al imponer esta restricci�n:
count if pscore!=.&pscore_sc==.

**** Veamos el resultado gr�ficamente:

drop x1 y1 x0 y0

histogram pscore_sc, by(D)

histogram pscore_sc if D==1, bin(100) color(blue) addplot(kdensity pscore_sc if D==1)

kdensity pscore_sc if D==1, epanechnikov generate(x1 y1)

histogram pscore_sc if D==0, bin(100) color(blue) addplot(kdensity pscore_sc if D==0)

kdensity pscore_sc if D==0, epanechnikov generate(x0 y0)

twoway (line y1 x1) (line y0 x0, lpattern(dash)), ytitle(Densidad) xtitle(Probabilidad de ser tratado) title(Propensity Score con Soporte Com�n ) legend(order(1 "Participante=1" 2 "No participante=0"))

**** Las probabilidades almacenadas en la variable "pscore_sc" cumplen con la propiedad del soporte com�n mediante el m�ximo y el m�nimo. 




*------------------------------*
*3. Calidad del emparejamiento *
*------------------------------*

**** Otra posibilidad para estimar la probabilidad de participaci�n es utilizar directamente el comando "pscore".
**** El comando "pscore" primero determina la probabilidad de participaci�n para cada individuo de acuerdo con el modelo que 
**** uno especifique. Posterior a esto, se dividen las observaciones en un n�mero �ptimo de bloques de manera que dentro de �stos
**** la probabilidad media del grupo de control no es estad�sticamente diferente de la probabilidad media del grupo de tratamiento. 
**** Este es el primer paso para balancear la probabilidad de participaci�n. Si se encuentra que dentro de un mismo bloque la probabilidad de 
**** participaci�n es estad�sticamente diferente, se divide el bloque en dos. Una vez se determina el n�mero de bloques mediante este procedimiento,
**** el programa prueba, bloque por bloque, que no existan diferencias estad�sticamente significativas entre los individuos de tratamiento
**** y control para las variables incluidas para predecir la probabilidad de participaci�n. Luego de esto, impone el soporte com�n.


pscore D $X, pscore(pscore_b) blockid(id) comsup det

**** La variable "id" almacena el bloque al que pertenece cada observaci�n

**** Veamos las diferencias entre el soporte com�n calculado manualmente y el estimado por "pscore":

sum pscore_sc, detail

sum pscore_b, detail

**** Gr�ficamente, veamos las diferencias: ya ten�amos la gr�fica del soporte com�n haci�ndolo manualmente.
**** Ahora veamos el soporte com�n estimado por "pscore"


drop x1 y1 x0 y0

histogram pscore_b if comsup==1, by(D)

histogram pscore_b if D==1 & comsup==1, bin(100) color(blue) addplot(kdensity pscore_b if D==1 & comsup==1)

kdensity pscore_b if D==1 & comsup==1, epanechnikov generate(x1 y1)

histogram pscore_b if D==0 & comsup==1, bin(100) color(blue) addplot(kdensity pscore_b if D==0 & comsup==1)

kdensity pscore_b if D==0 & comsup==1, epanechnikov generate(x0 y0)

twoway (line y1 x1) (line y0 x0, lpattern(dash)), ytitle(Densidad) xtitle(Probabilidad de ser tratado) title("pscore" ) legend(order(1 "Participante=1" 2 "No participante=0"))

**** Resulta evidente que hay observaciones pertenecientes al grupo de control que tienen una probabilidad predicha inferior al m�nimo del grupo de tratamiento.

sum pscore_b if D==1 & comsup==1

sum pscore_b if D==0 & comsup==1

sum pscore_sc if D==1

sum pscore_sc if D==0 

**** Otro procedimiento para garantizar la calidad del emparejamiento consiste en estimar el modelo probit con las caracter�sticas
**** especificadas controlando por la probabilidad predicha. En teor�a, los coeficientes asociados a las caracter�sticas de los individuos
**** no deben ser estad�sticamente significativos:

dprobit D pscore $X

**** Entonces, dado que ning�n coeficiente es significativo, podemos estar seguros de que el 
**** emparejamiento es adecuado. 

*-------------------------------------------------------------------*
*4. Selecci�n de un algoritmo de emparejamiento (talla para la edad)*
*-------------------------------------------------------------------*

**** 4.1 Al utilizar el emparejamiento por vecino m�s cercano, el programa "psmatch2"
**** empareja a cada individuo del grupo de tratamiento con el individuo del grupo de control
**** que tiene una probabilidad m�s cercana. Sin embargo, cuando se presentan casos en los que
**** hay varios individuos en el grupo de control a la misma distancia de un solo individuo de tratamiento
**** debemos asegurarnos de que el orden en el que se presentan los datos en nuestra base sea aleatorio.
**** Definimos la semilla en 10 arbitrariamente con el objetivo de poder replicar los resultados.

set seed 50
drawnorm orden
sort orden

**** 4.1 Estimador PSM por vecino m�s cercano. El comando "psmatch2" calcula el propensity score
**** y tiene en cuenta �nicamente a quienes est�n dentro del soporte com�n, o nosotros podemos
**** indicarle qu� variable debe utilizar como propensity score:

psmatch2 D $X, outcome(ha_nchs2) n(1) com

drop x1 y1 x0 y0

histogram pscore if _support==1, by(D)

histogram pscore if D==1 & _support==1, bin(100) color(blue) addplot(kdensity pscore_b if D==1 & _support==1)

kdensity pscore if D==1 & _support==1, epanechnikov generate(x1 y1)

histogram pscore if D==0 & _support==1, bin(100) color(blue) addplot(kdensity pscore_b if D==0 & _support==1)

kdensity pscore if D==0 & _support==1, epanechnikov generate(x0 y0)

twoway (line y1 x1) (line y0 x0, lpattern(dash)), ytitle(Densidad) xtitle(Probabilidad de ser tratado) title("psmatch2" ) legend(order(1 "Participante=1" 2 "No participante=0"))

sum pscore if D==1 & _support==1

sum pscore if D==0 & _support==1


**** Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(ha_nchs2) n(1) trim(20)

**** 4.2 Matching con 5 vecinos


**** 4.2.1 Imponiendo el soporte com�n del comando.

psmatch2 D $X, outcome(ha_nchs2) n(5) com

**** 4.2.1 Imponiendo el soporte com�n mediante trimming.

psmatch2 D $X, outcome(ha_nchs2) n(5) trim(20)


**** 4.3 Matching con 10 vecinos

**** 4.3.1 Imponiendo el soporte com�n del comando:

psmatch2 D $X, outcome(ha_nchs2) n(10) com

**** 4.3.2 Imponiendo el soporte com�n mediante trimming:

psmatch2 D $X, outcome(ha_nchs2) n(10) trim(20)

**** 4.4.1 Emparejamiento de distancia m�xima soporte com�n com

psmatch2 D $X, outcome(ha_nchs2) radius caliper(0.001) com

psmatch2 D $X, outcome(ha_nchs2) radius caliper(0.005) com		

**** 4.4.2 Emparejamiento de distancia m�xima soporte com�n trimming

psmatch2 D $X, outcome(ha_nchs2) radius caliper(0.001) trim(20)

psmatch2 D $X, outcome(ha_nchs2) radius caliper(0.005) trim(20)	

**** 4.5.1 Emparejamiento por kernel soporte com�n com

psmatch2 D $X, outcome(ha_nchs2) com kernel

bootstrap r(att) : psmatch2 D $X, out(ha_nchs2) com kernel

**** 4.5.2 Emparejamiento por kernel soporte com�n trimming

psmatch2 D $X, outcome(ha_nchs2) trim(20) kernel

bootstrap r(att) : psmatch2 D $X, out(ha_nchs2) trim(20) kernel

***** 4.6.1 Estimador por lineal local soporte com�n com

psmatch2 D $X, llr outcome(ha_nchs2) common

bootstrap r(att) : psmatch2 D $X, llr outcome(ha_nchs2) common

***** 4.6.2 Estimador por lineal local soporte com�n trimming(20)

psmatch2 D $X, llr outcome(ha_nchs2) trim(20)

bootstrap r(att) : psmatch2 D $X, llr outcome(ha_nchs2) trim(20)

*---------------------------------------------------------------------*
*5. Selecci�n de un algoritmo de emparejamiento (desnutrici�n cr�nica)*
*---------------------------------------------------------------------*



**** 5.1 Estimador PSM por vecino m�s cercano.

**** 5.1.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(desn_cr) n(1) com

**** 5.1.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(desn_cr) n(1) trim(20)

**** 5.2 Matching con 5 vecinos

**** 5.2.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(desn_cr) n(5) com

**** 5.2.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(desn_cr) n(5) trim(20)


**** 5.3 Matching con 10 vecinos

**** 5.3.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(desn_cr) n(10) com

**** 5.3.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(desn_cr) n(10) trim(20)

**** 5.4.1 Emparejamiento de distancia m�xima soporte com�n com

psmatch2 D $X, outcome(desn_cr) radius caliper(0.001) com

psmatch2 D $X, outcome(desn_cr) radius caliper(0.005) com		

**** 5.4.2 Emparejamiento de distancia m�xima soporte com�n trimming

psmatch2 D $X, outcome(desn_cr) radius caliper(0.001) trim(20)

psmatch2 D $X, outcome(desn_cr) radius caliper(0.005) trim(20)	

**** 5.5.1 Emparejamiento por kernel soporte com�n com

psmatch2 D $X, outcome(desn_cr) com kernel

bootstrap r(att) : psmatch2 D $X, out(desn_cr) com kernel

**** 5.5.2 Emparejamiento por kernel soporte com�n trimming

psmatch2 D $X, outcome(desn_cr) trim(20) kernel

bootstrap r(att) : psmatch2 D $X, out(desn_cr) trim(20) kernel

***** 5.6.1 Estimador por lineal local soporte com�n com

psmatch2 D $X, llr outcome(desn_cr) common

bootstrap r(att) : psmatch2 D $X, llr outcome(desn_cr) common

***** 5.6.2 Estimador por lineal local soporte com�n trimming

psmatch2 D $X, llr outcome(desn_cr) trim(20)

bootstrap r(att) : psmatch2 D $X, llr outcome(desn_cr) trim(20)


*-----------------------------------* 
* 6. Dobles diferencias emparejadas *
*-----------------------------------*

**** En este cap�tulo utilizaremos la metodolog�a de dobles diferencias emparejadas para evaluar
**** la evoluci�n de la variable "talla para la edad" de los individuos en dos periodos.

**** Lo primero que debemos hacer es generar la variable de diferencia entre las dos observaciones:

gen delta_ha=ha_nchs2-ha_nchs1

**** Ahora utilizamos propensity score matching para evaluar la diferencia en la evoluci�n de la talla para
**** la edad entre individuos tratados y no tratados:


**** 6.1 Estimador PSM por vecino m�s cercano.

**** 6.1.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(delta_ha) n(1) com

**** 6.1.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(delta_ha) n(1) trim(20)

**** 6.2 Matching con 5 vecinos

**** 6.2.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(delta_ha) n(5) com

**** 6.2.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(delta_ha) n(5) trim(20)


**** 6.3 Matching con 10 vecinos

**** 6.3.1 Imponiendo el soporte com�n mediante el comando:

psmatch2 D $X, outcome(delta_ha) n(10) com

**** 6.3.2 Ahora imponemos el soporte com�n mediante trimming:

psmatch2 D $X, outcome(delta_ha) n(10) trim(20)

**** 6.4.1 Emparejamiento de distancia m�xima soporte com�n com

psmatch2 D $X, outcome(delta_ha) radius caliper(0.001) com

psmatch2 D $X, outcome(delta_ha) radius caliper(0.005) com		

**** 6.4.2 Emparejamiento de distancia m�xima soporte com�n trimming

psmatch2 D $X, outcome(delta_ha) radius caliper(0.001) trim(20)

psmatch2 D $X, outcome(delta_ha) radius caliper(0.005) trim(20)	


**** 6.5.1 Emparejamiento por kernel soporte com�n com

psmatch2 D $X, outcome(delta_ha) com kernel

bootstrap r(att) : psmatch2 D $X, out(delta_ha) com kernel

**** 6.5.2 Emparejamiento por kernel soporte com�n trimming

psmatch2 D $X, outcome(delta_ha) trim(20) kernel

bootstrap r(att) : psmatch2 D $X, out(delta_ha) trim(20) kernel

***** 6.6.1 Estimador por lineal local soporte com�n com

psmatch2 D $X, llr outcome(delta_ha) common

bootstrap r(att) : psmatch2 D $X, llr outcome(delta_ha) common

***** 6.6.2 Estimador por lineal local soporte com�n trimming

psmatch2 D $X, llr outcome(delta_ha) trim(20)

bootstrap r(att) : psmatch2 D $X, llr outcome(delta_ha) trim(20)


**** 7. Utilizaci�n de pesos muestrales en matching.

**** La metodolog�a de propensity score matching tambi�n permite incluir pesos muestrales en sus estimaciones. 
**** Con el comando presentado en este archivo (psmatch2) se puede dar la opci�n de utilizar pesos muestrales
**** mediante el siguiente comando:

**** psmatch2 t x, outcome(y) w(Z)  n(1) com

**** En este caso, Z ser�a la matriz de pesos a utilizar. Otra posibilidad es utilizar el comando
**** "attnd". Su instalaci�n tambi�n se hace en la pesta�a de b�squeda de Stata. Una vez instalado,
**** el comando para hacer matching con pesos es el siguiente:

**** attnd y t x, [pweight=z]

**** En este caso z ser�a la variable que indica el peso muestral otorgado a cada observaci�n.



