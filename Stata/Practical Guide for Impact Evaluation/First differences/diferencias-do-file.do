*----------------------------------*
*      Modelo de diferencias       *
*----------------------------------*



*----------------------------------------------------------------------------------------------------------------*
*Algunos puntos antes de hacer el an�lisis:

*El s�mbolo asterisco (*) reconoce una linea del do-file como comentarios
*y cuando se corre el do-file las l�neas precedidas por este s�mbolo no son 
*reconocidas como comandos.

*Lo primero que debemos hacer es abrir la base de datos del modelo de diferencias

*Para simplificar los comandos definimos un vector X compuesto de las principales caracter�sticas del hogar

global X "personas hombre orden_n baja ocupado_jefe educa_jefe ingresos_hogar_jefe"
*----------------------------------------------------------------------------------------------------------------*


				***************************
				**** Talla para la edad****
				***************************


*--------------------------------------*
*1. Estad�sticas descriptivas
*--------------------------------------*


sum ha_nchs, detail

**** Podemos ver algunas caracter�sticas de la variable ha_nchs. Entre estas, el m�ximo, m�nimo, desviaci�n est�ndar y promedio.

*--------------------------------------*
*2. Estad�sticas descriptivas por grupo
*--------------------------------------*

sum ha_nchs if D==1, detail

sum ha_nchs if D==0, detail

*-------------------------------------------------------------------*
*3.Regresi�n con "D" siendo la �nica variable explicativa *
*-------------------------------------------------------------------*

reg ha_nchs D

**** Mediante la estimaci�n de MCO podemos ver una prueba de diferencia de medias entre los dos grupos (tratamiento
**** y control), sin controles adicionales. 

*------------------------------------------------*
*4. Regresi�n adicionando variables explicativas
*------------------------------------------------*


reg ha_nchs D $X

**** Dado que en realidad la variable talla para la edad depende de factores adicionales
**** como el ingreso del jefe de hogar, el orden de nacimiento, la educaci�n del jefe de hogar
**** y el n�mero de personas en el hogar, entre otras, realizamos la prueba de diferencia de medias
**** entre los dos grupos pero ahora controlando por estas variables.



*---------------------------------------------------------------------*
*5. Estimador de diferencias con efectos heterog�neos			    *
*---------------------------------------------------------------------*


reg ha_nchs D $X D_baja

**** Al incorporar la variable "D_baja"="D"*"baja" podemos encontrar que el 
**** tratamiento tiene efectos heterogeneos sobre los individuos dependiendo de la 
**** raza a la que pertenecen.


*--------------------------------------------------*
*6. Prueba estad�stica "D-D_baja=0".
*--------------------------------------------------*

reg ha_nchs D $X D_baja

test D+D_baja=0

**** Al correr la anterior regresi�n pudimos ver que el impacto del tratamiento
**** disminuye si el individuo pertenece a la raza "baja". Sin embargo, al probar que
**** "D+D_baja" es diferente de cero, y al encontrar que "D"
**** y "D_baja" son significativas de manera individual, encontramos que pese a 
**** que el programa tiene un efecto reducido sobre la talla para la edad en los individuos
**** de raza baja, este efecto es positivo.



				*****************************
				**** Desnutrici�n Cr�nica****
				*****************************

*----------------------------*
*1. Estad�sticas descriptivas
*----------------------------*


sum desn_cr

**** Por medio de este comando podemos verificar que 
**** aproximadamente 11,75% de los individuos en la muestra sufren
**** de desnutrici�n cr�nica.


*----------------------*
*2. Promedios por grupo
*----------------------*


sum desn_cr if D==0

sum desn_cr if D==1

**** Encontramos los datos para desnutrici�n cr�nica para el grupo de tratamiento y control.


*--------------------------------------------------------------------*
*3. Regresi�n con "D" siendo la �nica variable explicativa
*--------------------------------------------------------------------*

dprobit desn_cr D

logit desn_cr D

mfx

**** Al igual que con la talla para la edad, lo primero que hacemos es una prueba de diferencia de medias sin controlar por otras variables.
**** Dado que la variable dependiente es dicot�mica, utilizamos un modelo probit y un modelo logit, respectivamente. 
**** Es importante recordar que la diferencia entre estos dos modelos es la distribuci�n que suponen para
**** la procedencia de los datos. Mientras que el modelo probit supone una distribuci�n normal, logit supone una log�stica.

**** El comando "dprobit" nos permite identificar los efectos marginales en probabilidad. "mfx" hace lo mismo para el modelo logit.

*-----------------------------------------------*
*4. Regresi�n adicionando variables explicativas
*-----------------------------------------------*


dprobit desn_cr D $X

logit desn_cr D $X

mfx


**** Ahora hacemos la diferencia de medias pero controlando por m�s variables. Entre estas la educaci�n del jefe del hogar,
**** el orden de nacimiento del individuo, los ingresos del jefe del hogar, g�nero del individuo, raza del individuo etc.

				****************************************
				**** Verificaci�n de aleatorizaci�n ****
				****************************************

dprobit D $X

logit D $X

mfx

**** En este �ltimo paso verificamos que en realidad el tratamiento no se asigna de acuerdo
**** a caracter�sticas observables de los individuos. Ninguna de las variables es significativa
**** en ninguno de los dos modelos y ning�n modelo es significativo. Esto nos garantiza la aleatorizaci�n
**** en la asignaci�n del tratamiento.

