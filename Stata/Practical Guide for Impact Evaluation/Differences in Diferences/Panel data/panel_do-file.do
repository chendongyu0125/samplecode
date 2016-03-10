                          *------------------------------------------------------------------------------------*
				  *-------------------------    Diferencias en diferencias    -------------------------*
				  *------------------------------------------------------------------------------------*

								  	    *    Datos Panel      *

*--------------------------------------------------------------------------------------------------------------------------------------------*
*Lo primero que debemos hacer es abrir la base de datos del modelo de diferencias:

*Esta base de datos cuenta con informaci�n para 4.000 individuos. Para cada individuo, tenemos informaci�n de dos periodos. El sub�ndice de cada
*variable corresponde al periodo de observaci�n. Por ejemplo "ha_nchs2" indica el ingreso mensual, en decenas de miles de pesos
*del jefe del hogar en el periodo 2. 
*--------------------------------------------------------------------------------------------------------------------------------------------*



*-------------------------------------------*
*1. Veamos algunas estad�sticas descriptivas*
*-------------------------------------------*


*Talla para la edad de los individuos del grupo de tratamiento en el primer periodo
sum ha_nchs1 if D==1

*Talla para la edad de los individuos del grupo de control en el primer periodo
sum ha_nchs1 if D==0

*Talla para la edad de los individuos del grupo de tratamiento en el segundo periodo
sum ha_nchs2 if D==1

*Talla para la edad de los individuos del grupo de control en el segundo periodo
sum ha_nchs2 if D==0

**** Aparentemente, existe una diferencia entre los individuos de tratamiento y de control en el primer periodo en lo referente
**** a la variable talla para la edad "ha_nchs. Adem�s, esta diferencia se incrementa en el tiempo.

**** Haciendo una prueba de diferencia de medias, podemos ver que en los dos periodos hay una diferencia estad�sticamente significativa
**** que es mayor en el segundo periodo:

ttest ha_nchs1, by(D)

ttest ha_nchs2, by(D)

**** Mientras que en el primer periodo la diferencia es de 0.1, para el segundo se incrementa hasta 0.3. 


*--------------------------------------------------------*
*2. Modelo de diferencias en diferencias con datos panel *
*--------------------------------------------------------*


**** Al ver las estad�sticas descriptivas de la variable talla para la edad y desnutrici�n cr�nica
**** vemos que en realidad hay una diferencia preexistente entre los individuos pertenecientes al grupo de tratamiento
**** y al grupo de control. En el periodo despu�s de aplicar el tratamiento es posible que la diferencia entre los dos grupos
**** se deba a la diferencia preexistente o a la aplicaci�n del tratamiento. El modelo de diferencias en diferencias nos
**** permite controlar por diferencias existentes antes de la aplicaci�n del programa. 

**** Para correr el modelo b�sico de diferencias en diferencias, utilizando una base de datos panel, 
**** utilizamos como variable dependiente el cambio en la talla para la edad en funci�n del tratamiento. 
**** Primero debemos generar la variable "delta_ha_nchs" y despu�s si podemos correr el programa:

gen delta_ha_nchs=ha_nchs2-ha_nchs1


reg delta_ha_nchs D


**** El coeficiente asociado a la variable "D" tiene una magnitud de 0.2 y es significativo. 
**** Esto nos indica que la diferencia de talla para la edad entre los individuos de tratamiento y
**** control se incrementa en 0.2 por la aplicaci�n del tratamiento. 

*------------------------------------------------------------------------------------------*
*3. Modelo de diferencias en diferencias con datos panel utilizando regresores adicionales *
*------------------------------------------------------------------------------------------*


**** Podemos suponer que el crecimiento en la variable talla para la edad, adem�s de deberse al tratamiento y 
**** al crecimiento natural en el tiempo, puede tener su origen en caracter�sticas de los individuos. 
**** Por ejemplo, resulta razonable suponer que el crecimiento en la talla para la edad de un individuo
**** tambi�n est� asociado al ingreso en el primer periodo del jefe del hogar. En este punto verificamos esta hip�tesis
**** con los ingresos del jefe de hogar y la educaci�n del jefe de hogar. 



reg delta_ha_nchs D ingresos_hogar_jefe1 educa_jefe1 


