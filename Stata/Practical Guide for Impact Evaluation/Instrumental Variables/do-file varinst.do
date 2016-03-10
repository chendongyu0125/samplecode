*-----------------------------------*
*    Variables instrumentales       *
*-----------------------------------*

*----------------------------------------------------------------------------------------------------------------*
*En el cap�tulo de variables instrumentales se utilizar� el comando "ivreg2". Dado que este comando no viene
*incluido en el paquete b�sico de Stata, es necesario instalarlo. Para instalarlo debemos seleccionar 
*la pesta�a "help" de Stata, en la opci�n search y seleccionamos
*search all. Buscamos "ivreg2" y seleccionamos el paquete st0030_2 para instalarlo. 
*----------------------------------------------------------------------------------------------------------------*

**** Para simplificar el an�lisis definimos un vector global con las variables que se han utilizado antes:

global X "personas orden_n ocupado_jefe educa_jefe ingresos_hogar_jefe"

**** Cuando un regresor tiene alguna relaci�n con el t�rmino de error, se viola uno de los supuestos b�sicos
**** de MCO y por lo tanto, estimaciones mediante este m�todo arrojar�an resultados sesgados. En nuestro
**** ejemplo, el tratamiento puede estar relacionado con variables no observables como el nivel de preocupaci�n
**** de una madre frente a su hijo o el conocimiento que �sta tiene sobre el desarrollo adecuado de los ni�os.

**** Veamos los resultados de una estimaci�n simple por MCO:

reg ha_nchs D $X

**** Los resultados nos indican que el tratamiento tiene un efecto de 0.248 en la talla para la edad
**** de los individuos. Sin embargo, dado el problema de autoselecci�n, este resultado puede estar sesgado. 
**** Para proceder a utilizar la metodolog�a de variables instrumentales debemos seleccionar instrumentos
**** que sean relevantes y ex�genos. En principio, el n�mero de oficinas operadoras en el municipio determina
**** en cierta medida la probabilidad de participaci�n sin tener relaci�n con las variables no observables. 
**** Lo mismo sucede con la distancia del hogar a la oficina operadora m�s cercana. Antes de utilizar
**** estos instrumentos realizamos pruebas de relevancia y exogeneidad:

*------------------*
* 1. Relevancia    *
*------------------*

**** Para evaluar la relevancia de un instrumento se pueden utilizar diferentes pruebas. En este caso vamos a utilizar
**** una estimaci�n por MCO, la prueba can�nica de Anderson, la prueba de Cragg-Donald y la prueba de Stock y Yoko. 

**** 1.1 Distancia

**** 1.1.1 Regresi�n MCO 

reg D distancia

reg D $X distancia

**** Vemos que la distancia s� puede ser un buen predictor de la variable "D". 
**** El valor del estad�stico F es 19,23. Seg�n Stock y Watson pg. 443, (2.007), si el estad�stico
**** F es superior a diez, se puede asegurar que el instrumento es relevante. Adem�s, todos los coeficientes
**** son estad�sticamente significativos.

**** 1.1.2 Prueba can�nica de Anderson

ivreg2 ha_nchs $X (D=distancia)

**** La hip�tesis nula de esta prueba es que la ecuaci�n est� subidentificada. Dado que la rechazamos, la
**** prueba can�nica de Anderson (Underidentification test) nos permite estar m�s seguros de la relevancia de nuestro instrumento.

**** 1.1.3 Prueba Cragg-Donald

ivreg2 ha_nchs $X (D=distancia), first

**** La prueba Cragg-Donald tiene la misma hip�tesis nula que la prueba can�nica de Anderson pero 
**** el estad�stico de prueba se construye de manera diferente. En este caso, tambi�n rechazamos
**** la hip�tesis nula. 

**** 1.1.4 Prueba de Stock y Yoko

ivreg2 ha_nchs $X (D=distancia), first

**** Para esta prueba utilizamos el estad�stico Cragg-Donald y nos remitimos al valor reportado en la tabla de Stock y Yoko (2002).

*Realizamos las mismas pruebas para el instrumento "of_op":

**** 1.2 of_op

**** 1.2.1 Regresi�n MCO 

reg D of_op

reg D $X of_op

**** 1.2.2 Prueba can�nica de Anderson

ivreg2 ha_nchs $X (D=of_op)

**** 1.2.3 Prueba Cragg-Donald

ivreg2 ha_nchs $X (D=of_op), first

**** 1.2.4 Prueba de Stock y Yoko

ivreg2 ha_nchs $X (D=of_op), first

**** Ahora realizamos las mismas pruebas utilizando los dos instrumentos:

**** 1.3 Distancia y oficinas operadoras:

**** 1.3.1 Regresi�n MCO 

reg D of_op distancia

reg D $X of_op distancia

test distancia=of_op=0

**** 1.3.2 Prueba can�nica de Anderson

ivreg2 ha_nchs $X (D=of_op distancia)

**** 1.3.3 Prueba Cragg-Donald

ivreg2 ha_nchs $X (D=of_op distancia), first

**** 1.3.4 Prueba de Stock y Yoko

ivreg2 ha_nchs $X (D=of_op distancia), first

**** Podemos concluir, por lo tanto, que los instrumentos son relevantes.


*------------------*
* 2. Exogeneidad   *
*------------------*

**** Al igual que probamos la relevancia de los instrumentos, la exogeneidad la probaremos
**** en cada uno de los instrumentos y despu�s utilizando los dos. 

**** 2.1 Distancia

**** 2.1.1 Regresi�n de errores en funci�n del instrumento.

**** Lo que vamos a hacer es predecir los errores de MCO y verificar que estos no est�n relacionados
**** con el instrumento

reg ha_nchs D $X

predict uhat, residuals

reg uhat distancia

**** Vemos que la variable distancia no es significativa y tampoco lo es el modelo.
**** Para utilizar las otras pruebas de exogeneidad es necesario que el n�mero
**** de instrumentso exceda el n�mero de variables a instrumentar. Estas pruebas
**** ser�n llevadas a cabo cuando evaluemos la exogeneidad de la distancia y las 
**** oficinas operadoras de manera conjunta.

**** 2.2 Oficinas operadoras

**** 2.2.1 Regresi�n de errores en funci�n del instrumento

reg uhat of_op

**** Vemos que of_op aparentemente no predice bien los errores.

**** 2.3 Distancia y Oficinas operadoras

**** 2.2.1 Regresi�n de errores en funci�n de los instrumentos

reg uhat of_op distancia

**** 2.2.2 Prueba de Sargan

ivreg2 ha_nchs $X (D=of_op distancia)

**** La hip�tesis nula de esta prueba es que los instrumentos son ex�genos. Rechazarla nos llevar�a a dudar
**** sobre la validez de los instrumentos. En este caso no podemos rechazar esta prueba.

**** 2.2.3 Sargan y Basmann

overid

**** Utilizando el comando "overid" luego de llevar a cabo una estimaci�n por variables instrumentales
**** se lleva a cabo la prueba de Sargan y la de Basmann. La prueba de Sargan ya la vimos, la de Basmann
**** es una prueba an�loga. Sin embargo, el estad�stico de prueba se construye de manera diferente. 

****2.2.4 Prueba "J".

**** La prueba "J" es an�loga a la prueba de Sargan pero se utiliza cuando la estimaci�n se hace mediante
**** el m�todo generalizado de momentos:

ivreg2 ha_nchs $X (D=of_op distancia), gmm2s robust


*----------------------------*
*  3. Estimaci�n por MC2E    *
*----------------------------*

**** Dado que ya probamos que los instrumentos que tenemos pueden ser utilizados, procedemos
**** a realizar la estimaci�n mediante MCO en dos etapas. 

**** 3.1 Primera etapa, predicci�n de la variable end�gena (D) utilizando los 
**** regresores ex�genos y los instrumentos (distancia y of_op):

drop uhat

probit D $X distancia of_op

predict trat_hat


**** 3.2 Segunda etapa, regresi�n donde la variable dependiente es la talla para la edad y 
**** se utilizan las variables predichas por la anterior regresi�n (trata_hat) para instrumentar
**** la participaci�n en el programa:

ivreg2 ha_nchs  $X (D=trat_hat), first

**** Para llevar a cabo la estimaci�n de manera directa con el comando "ivreg2" instrumentando el tratamiento con la distancia y el
**** n�mero de oficinas operadoras en el municipio:

ivreg2 ha_nchs $X (D=distancia of_op), first

**** Utilizando m�nimos cuadrados en dos etapas vemos que el efecto del tratamiento se disminuye significativamente. 
**** Mediante MCO el efecto del tratamiento en la talla para la edad era de 0.24 desviaciones est�ndar. Sin embargo,
**** mediante MC2E pudimos confirmar que en realidad el efecto es 0.21. 
	