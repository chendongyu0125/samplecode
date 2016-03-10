*-------------------------------------*
*      Regression discontinuity       *
*-------------------------------------*

*----------------------------------------------------------------------------------------------------------------*
*En regression discontinuity  el programa "rd".
*Dado que el paquete b�sico de Stata no contiene el comando "rd" es necesario instalarlo:

**** Instalaci�n "rd"
ssc install rd, replace

cd "C:\Documents and Settings\r-azuero\Mis documentos\Rodrigo\EdeI\Archivos para entregar\Gu�a pr�ctica para la evaluaci�n de impacto\Regresi�n discontinua"
use discontinuity_base.dta, clear

*Sobre la base de datos: Hay dos observaciones para la variable talla para la edad. Tenemos el puntaje de sisben
* que va de 0 a 30. Se supone que los individuos con un puntaje sisben inferior a 10 tienen una mayor probabilidad
* de participar en el programa. 
*----------------------------------------------------------------------------------------------------------------*

**** En este do-file presentaremos las cuatro metodolog�as presentadas en el texto. Estas son, Imbens Lemieux, Kernel Triangular,
**** Metodolog�a de Imbens y Lemieux con Kernel triangular y finalmente la metodolog�a propuesta por el programa "RD".
**** Como se explica en el documento, los resultados de las estimaciones dependen del ancho de banda utilizado. Por
**** esta raz�n se har�n estimativos para todas las metodolog�as y se almacenar�n en una matriz de manera que
**** al final tendremos resultados de todas las metodolog�as para todos los anchos de banda posibles. 


**** I. Metodolog�a de Imbens y Lemieux

**** Los resultados de la metodolog�a Imbens y Lemieux se almacenar�n en el vector IL.

mat IL=J(10,1,0)

**** Lo que hacemos es generar un ciclo para el ancho de banda. Hacemos estimativos para anchos de banda 
**** desde uno hasta diez. Debemos hacer cuatro regresiones para cada ancho de banda. Para la talla para la edad
**** corremos una regresi�n a la izquierda y otra a la derecha. Lo mismo hacemos para la probabilidad de participar
**** en el programa. Luego de correr cada regresi�n, el resultado de la constante queda guardado como un escalar. 
**** El sub�ndice utilizado para identificar la constante de la regresi�n del tratamiento por la izquierda es alfa_`h'l
**** donde `h' es el ancho de banda. Para las regresiones que se corren a la derecha, las constantes quedan guardadas como 
**** alfa_`h'r. Cuando hacemos las regresiones por izquierda y derecha para estimar la probabilidad de participaci�n en el programa
**** las constantes se guardan con el nombre alfa_d`h'l y alfa_d`h'r respectivamente. Finalmente, se construye el estimador
**** del efecto del programa mediante la relaci�n de las diferencias de las constantes:
**** Efecto=[alfa_`h'l-alfa_`h'r]/[alfa_d`h'l -alfa_d`h'r].

**** Debemos generar una variable que normalice el puntaje de sisben:

gen sisbenalin=sisben-10

local h=1
while `h'<=10 {
	display " Regresi�n talla para la edad por izquierda. Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben<10 & sisben>=10-`h'
	mat define A_10=e(b)
	scalar define alfa_`h'l=el(A_10,1,2)
	mat drop A_10
	
	display " Regresi�n talla para la edad por derecha. Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben>=10 & sisben<=10+`h'
	mat define A_10=e(b)
	scalar define alfa_`h'r=el(A_10,1,2)
	mat drop A_10

	display "Regresi�n de la probabilidad de participaci�n por izquierda. Ancho de banda `h'"
	reg D sisbenalin if  sisben<10 & sisben>=10-`h'
	mat define A_10=e(b)
	scalar define alfa_d`h'l=el(A_10,1,2)
	mat drop A_10

	display "Regresi�n de la probabilidad de participaci�n por derecha. Ancho de banda `h'"
	reg D sisbenalin if  sisben>=10 & sisben<=10+`h'
	mat define A_10=e(b)
	scalar define alfa_d`h'r=el(A_10,1,2)
	mat drop A_10

	mat IL[`h',1]=(alfa_`h'l-alfa_`h'r)/(alfa_d`h'l-alfa_d`h'r)

local h=`h'+1
}
mat coln IL=IL
mat list IL

**** Podemos ver entonces que los resultados quedan almacenados en el vector IL. Cada fila es un
**** ancho de banda diferente. Si se desea ver los resultados de las constantes estimadas 
**** utilizando el comando "scalar list" se pueden ver todas las estimaciones de las constante. 


**** II. Kernel triangular con la metodolog�a Imbens Lemieux

**** El Kernel triangular con metodolog�a Imebns y Lemieux es una combinaci�n de la metodolog�a
**** propuesta por Imbens y Lemieux pero con una ponderaci�n dada por un kernel triangular. 
**** Los resultados de las estimaciones para el efecto del programa ser�n almacenados en el vector IL_KT,
**** los resultados de la discontinuidad en la talla para la edad quedan en la matriz Y y del tratamiento
**** en la matriz T. La metodolog�a es id�ntica a la propuesta en el primer cap�tulo, lo �nico que cambia
**** es que cada observaci�n ser� ponderada de acuerdo a su distancia frente al umbral. Los pesos
**** de cada observaci�n ser�n almacenados en la variable k. Note que en el caso donde h=1 o h=2 
**** la estimaci�n es la misma que bajo la metodolog�a de IL dado que los pesos para observaciones
**** diferentes a aquellas con puntaje sisben 9, 10 u 11 es cero. 


mat IL_KT=J(10,1,0)
mat Y=J(1,10,0)
mat T=J(1,10,0)


gen norm1=.
gen k=.

local h=2
	while `h'<=10 {
	replace norm1=1-abs(sisbenalin/`h')
	replace k=0 if norm1<0
	replace k=norm1 if norm1>=0
	
	display " Regresi�n talla para la edad por izquierda. Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben<10 & sisben>=10-`h' [aw=k], robust
	mat define A_10=e(b)
	scalar define alfa_`h'l=el(A_10,1,2)
	mat drop A_10
	
	display " Regresi�n talla para la edad por derecha. Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben>=10 & sisben<=10+`h' [aw=k], robust
	mat define A_10=e(b)
	scalar define alfa_`h'r=el(A_10,1,2)
	mat drop A_10

	display " Regresi�n probabilidad participaci�n por izquierda. Ancho de banda `h'"
	reg D sisbenalin if  sisben<10 & sisben>=10-`h' [aw=k], robust
	mat define A_10=e(b)
	scalar define alfa_d`h'l=el(A_10,1,2)
	mat drop A_10

	display " Regresi�n probabilidad participaci�n por derecha. Ancho de banda `h'"
	reg D sisbenalin if  sisben>=10 & sisben<=10+`h' [aw=k], robust
	mat define A_10=e(b)
	scalar define alfa_d`h'r=el(A_10,1,2)
	mat drop A_10

	mat IL_KT[`h',1]=(alfa_`h'l-alfa_`h'r)/(alfa_d`h'l-alfa_d`h'r)

local h=`h'+1
}

mat IL_KT[1,1]=IL[1,1]
mat coln IL_KT=IL_KT
mat list IL_KT


**** III. Kernel triangular

**** A diferencia de la metodolog�a de Imbens y Lemieux, el kernel triangular nos permite
**** estimar la discontinuidad en una variable haciendo uso de una sola regresi�n. 
**** Los resultados quedan almacenados en el vector KT, los resultados de la discontinuidad
**** en la talla para la edad son guardados en el vector T y los resultados de la discontinuidad
**** en la probabilidad de participaci�n quedan almacenados en el vector T. 

mat KT=J(10,1,0)
mat Y=J(1,10,0)
mat T=J(1,10,0)

**** Debemos generar una variable que nos indique si el individuo est� encima o debajo del umbral:

gen W=0 if sisben>=10
replace W=1 if sisben<10

**** As� mismo, generamos la interacci�n entre esta variable reci�n creada y el puntaje de sisben alineado. 

gen W_sisbenalin=W*sisbenalin

**** La variable k ser� la que contenga los pesos de cada observaci�n.


local h=1
while `h'<=10 {
replace norm1=1-abs(sisbenalin/`h')
replace k=0 if norm1<0
replace k=norm1 if norm1>=0
display "Ancho de banda `h'"
reg ha_nchs W sisbenalin W_sisbenalin [aw=k], robust
mat define A_`h'=e(b)
mat Y[1,`h']=el(A_`h',1,1)
mat drop A_`h'
display "Ancho de banda `h'"
reg D W sisbenalin W_sisbenalin [aw=k], robust 
mat define A_T_`h'=e(b)
mat T[1,`h']=el(A_T_`h',1,1)
mat drop A_T_`h'
mat KT[`h',1]=el(Y,1,`h')/el(T,1,`h')
local h=`h'+1
}

mat coln KT=KT
mat list KT





**** IV.Comando RD:

**** El comando rd es bastante f�cil de manejar. Utiliza una estimaci�n mediante regresi�n local polin�mica
**** de grado uno y presenta su gr�fica si uno desea observarla. Veamos un ejemplo con un ancho de banda de 4:

rd ha_nchs D sisben, bw(2) mbw(100) k(gau) z0(10) gr 

**** Se especifican, en respectivo orden, la variable de resultado (ha_nchs), la variable
**** que indica si el individuo fue o no tratado (D) y la variable en donde que indica
**** el umbral de la discontinuidad (sisben). En caso de una regresi�n discontinua n�tida no ser� necesario
**** especificar una variable de tratamiento ya que todos los individuos a la derecha del umbral
**** son tratados y a la izquierda no lo ser�n. 

**** En este caso estamos utilizando un kernel gaussiano. Sin embargo, la estimaci�n se 
**** puede hacer con otro tipo de kernel. La opci�n gr nos permite ver los dos gr�ficos
**** (en la talla para la edad y en la probabilidad de participaci�n) con su respectiva
**** continuidad. z0(10) nos indica que la discontinuidad se presenta en el nivel de sisben 10. 
**** Dado que estamos utilizando un m�todo no param�trico de estimaci�n, no hay un m�todo sencillo
**** para estimar los errores est�ndar. Por esta raz�n debemos hacer bootstrap.

**** En este caso guardaremos los resultados de la estimaci�n del efecto del programa en el vector
**** RD. Los errores est�ndar ser�n almacenados en el vector EE y los estad�sticos t de cada estimaci�n en
**** el vector ET.


mat RD=J(10,1,0)
mat EE=J(10,1,0)
mat ET=J(10,1,0)

local h=1
	while `h'<=10 {
	rd ha_nchs D sisben, bw(`h') mbw(100) k(gau) z0(10) 
	mat define A_10=e(b)
	mat RD[`h',1]=el(A_10,1,3)
	bs: rd ha_nchs D sisben, bw(`h') mbw(100) k(gau) z0(10)
	mat define B_10=e(se)
	mat EE[`h',1]=el(B_10,1,3)
	mat ET[`h',1]=el(A_10,1,3)/el(B_10,1,3)



local h=`h'+1
}

mat coln RD=RD
mat list RD
mat list EE
mat list ET

**** Finalmente creamos la variable EFECTO donde las columnas indican el estimador para cada programa
**** y las filas el ancho de banda.

mat EFECTO=(IL,KT,IL_KT,RD)

mat list EFECTO




**** V. Selecci�n de h


**** Como hemos visto, todas las metodolog�as de estimaci�n necesitan un h especificado. 
**** En este punto presentamos la metodolog�a propuesta por Imbens y Lemieux para la selecci�n
**** del ancho de banda. Los valores de los Cross-validation Criteria quedar�n almacenados en 
**** el vector CR_V para la talla para la edad y en el vector CR_V_T para el tratamiento. 

mat define C=J(1,10,0)
mat define CR_V=J(10,1,0)

**** La gr�fica realizada a continuaci�n corresponde a la gr�fica 6.18 del libro en donde se presentan las observaciones tenidas en cuenta
**** en caso que el ancho de banda fuera 2 y el umbral fuera 7:

twoway (scatter ha_nchs sisben if sisben>=4 & sisben<=9)

**** En el siguiente programa comparamos la distancia de las observaciones frente al l�mite de la discontinuidad
**** estimado y seleccionamos el que minimice la distancia promedio. 

local h=1
while `h'<=10 {
gen cv_`h'=.

	display "Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben<10 & sisben>=10-`h'
	mat define A_`h'=e(b)
	mat C[1,`h']=el(A_`h',1,2)
	mat drop A_`h'
	replace cv_`h'=(ha_nchs-el(C,1,`h'))*(ha_nchs-el(C,1,`h')) if sisben<10 & sisben>=10-`h'


	display "Ancho de banda `h'"
	reg ha_nchs sisbenalin if  sisben<=10+`h' & sisben>=10
	mat define A_`h'=e(b)
	mat C[1,`h']=el(A_`h',1,2)
	mat drop A_`h'
	replace cv_`h'=(ha_nchs-el(C,1,`h'))*(ha_nchs-el(C,1,`h')) if sisben<=10+`h' & sisben>=10


sum cv_`h'
mat CR_V[`h',1]=r(mean)
drop cv_`h'
local h=`h'+1
}

mat list CR_V


**** Ahora hacemos este mismo procedimiento para la probabilidad de tratamiento y evaluamos
**** qu� nivel de ancho de banda reporta el menor Cross Validation Criteria. 


mat define C_T=J(1,10,0)
mat define CR_V_T=J(10,1,0)


local h=1
while `h'<=10 {
gen cv_t_`h'=.

	display "Ancho de banda `h'"
	reg D sisbenalin if  sisben<10 & sisben>=10-`h'
	mat define A_`h'=e(b)
	mat C_T[1,`h']=el(A_`h',1,2)
	mat drop A_`h'
	replace cv_t_`h'=(D-el(C_T,1,`h'))*(D-el(C_T,1,`h')) if sisben<10 & sisben>=10-`h'


	display "Ancho de banda `h'"
	reg D sisbenalin if  sisben<=10+`h' & sisben>=10
	mat define A_`h'=e(b)
	mat C_T[1,`h']=el(A_`h',1,2)
	mat drop A_`h'
	replace cv_t_`h'=(D-el(C_T,1,`h'))*(D-el(C_T,1,`h')) if sisben<=10+`h' & sisben>=10


sum cv_t_`h'
mat CR_V_T[`h',1]=r(mean)
drop cv_t_`h'
local h=`h'+1
}

mat list CR_V_T
mat list CR_V


**** Podemos ver que en ambos casos el "Cross Validation Criteria"
**** nos indica que lo �ptimo es utilizar un ancho de banda de 2. 



**** VI. Ejemplo con m�todo Delta.

**** En esta secci�n presentaremos c�mo estimar el error est�ndar del estimador presentado en (III) 
**** con el m�todo delta para un ancho de banda de dos (2).  

local h=2

	replace norm1=1-abs(sisbenalin/`h')
	replace k=0 if norm1<0
	replace k=norm1 if norm1>=0
	reg ha_nchs W sisbenalin W_sisbenalin [aw=k], robust
	predict er_`h' if e(sample), residuals 
	mat define A_`h'=e(b)
	mat define V_`h'=e(V)
	mat Y[1,`h']=el(A_`h',1,1)
	scalar define bet_`h'=el(A_`h',1,1)
	scalar define var_bet_`h'=el(V_`h',1,1)
	mat drop A_`h'
	mat drop V_`h'

	reg D W sisbenalin W_sisbenalin [aw=k], robust 
	predict er_D_`h' if e(sample), residuals
	mat define A_T_`h'=e(b)
	mat define V_T_`h'=e(V)
	mat T[1,`h']=el(A_T_`h',1,1)
	scalar define gam_`h'=el(A_T_`h',1,1)
	scalar define var_gam_`h'=el(V_T_`h',1,1)
	mat drop A_T_`h'
	mat drop V_T_`h'
	
	mat KT[`h',1]=el(Y,1,`h')/el(T,1,`h')
	
	gen int_err_`h'=er_D_`h'*er_`h'
	sum int_err_`h'
	scalar define cov_`h'=r(mean)
	
	local varianza_t_2=((var_bet_`h')/(gam_`h'^2)) + (((bet_`h'^2)*(var_gam_`h'))/(gam_`h'^4))-((2*bet_`h'*cov_`h')/(gam_`h'^3))
	drop int_err_`h' er_D_`h' er_`h'

display "La varianza del estimador presentado en III por el m�todo delta es `varianza_t_2'"

**** VII. Continuidad de los covariados. 
**** En este punto se verificar� el supuesto de continuidad en las otras variables en el 
**** umbral sisben=10


foreach var of varlist  personas educa_jefe ocupado_jefe hombre ingresos_hogar_jefe{
     rd `var' D sisben, bw(2) mbw(100) k(gau) z0(10) gr
     bs: rd `var' D sisben, bw(2) mbw(100) k(gau) z0(10) 
     }

rd ingresos_hogar_jefe D sisben, bw(2) mbw(100) k(gau) z0(10) gr
bs: rd ingresos_hogar_jefe D sisben, bw(2) mbw(100) k(gau) z0(10) 
