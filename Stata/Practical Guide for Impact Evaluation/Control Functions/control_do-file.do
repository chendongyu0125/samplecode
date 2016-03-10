*-------------------------------------*
*      Funciones de Control           *
*-------------------------------------*

*----------------------------------------------------------------------------------------------------------------*
*El programa que se utiliza para funciones de control es "treatreg"
*----------------------------------------------------------------------------------------------------------------*

**** En este caso se quiere identificar el impacto del programa "canasta" sobre la talla para la edad de los individuos.
**** Sin embargo, se sospecha que la decisi�n de asistir depende del nivel de capital social de los principales respondables
**** de los ni�os. Dado que no es posible observar este indicador directamente, lo �ptimo es utilizar funciones de control.

**** 1. Estimaci�n mediante M�nimos Cuadrados Ordinarios

ttest ha_nchs, by(D)

**** Una prueba de diferencia de medias nos muestra que los individuos que participan en el programa tienen un promedio de talla para 
**** la edad debajo de los individuos que no participan en el programa. Veamos qu� nos permite afirmar una regresi�n por m�nimos cuadrados ordinarios. 

global X "educa_jefe ocupado_jefe ingresos_hogar_jefe "


reg ha_nchs D $X 

**** Seg�n los resultados de la estimaci�n por m�nimos cuadrados ordinarios, el programa canasta tiene un efecto negativo sobre la
**** talla para la edad de los participantes. Sin embargo, no estamos controlando por el hecho de que la decisi�n de participaci�n
**** depende del nivel de capital social con el que cuentan los principales responsables de los ni�os, variable que tambi�n
**** afecta el nivel de talla para la edad de los individuos. En este punto utilizaremos el comando "treatreg".

**** 2. Estimaci�n mediante el comando Treatreg

**** El comando treatreg incorpora una serie de estimaciones para evaluaci�n de programas cuando la decisi�n de participar en el programa 
**** est� determinada por una variable no observable que a su vez determina la variable de inter�s del programa. Treatreg estima las ecuaciones:

**** Y=b0D+b2X+u
**** D=aZ+v

**** Donde Y es la variable de inter�s para la evaluaci�n del programa, X es un vector de caracter�sticas observables de los individuos,
**** D es una variable dicot�mica que nos indica si el individuo participa o no en el programa, Z es un vector de caracter�sticas observables
**** que determinan la probabilidad de participar en el programa y tanto u como v son t�rminos de error. La estructura del programa es la siguiente:

**** treatreg Y X, treat(D = Z) 

**** Al utilizar esta versi�n b�sica del programa con las variables de la base de datos, tenemos:

treatreg ha_nchs  $X, treat(D=$X of_op distancia)

**** El output obtenido es el siguiente:



**** Iteration 0:   log likelihood = -2637.3161  
**** Iteration 1:   log likelihood = -2637.1498  
**** Iteration 2:   log likelihood = -2637.1475  
**** Iteration 3:   log likelihood = -2637.1475  
**** 
**** Treatment-effects model -- MLE                    Number of obs   =       4000
****                                                   Wald chi2(4)    =      43.51
**** Log likelihood = -2637.1475                       Prob > chi2     =     0.0000
**** 
**** ------------------------------------------------------------------------------
****              |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
**** -------------+----------------------------------------------------------------
**** ha_nchs      |
****   educa_jefe |   .0026954   .0010125     2.66   0.008     .0007109    .0046799
**** ocupado_jefe |   -.015154   .0105939    -1.43   0.153    -.0359176    .0056096
**** ingresos_h~e |  -6.47e-06   .0000449    -0.14   0.885    -.0000945    .0000816
****            D |   .2306015   .0376373     6.13   0.000     .1568338    .3043693
****        _cons |  -.4602233   .0189447   -24.29   0.000    -.4973543   -.4230923
**** -------------+----------------------------------------------------------------
**** D            |
****   educa_jefe |  -.0076233   .0052429    -1.45   0.146    -.0178991    .0026525
**** ocupado_jefe |   .2212016    .052308     4.23   0.000     .1186798    .3237235
**** ingresos_h~e |    .000508   .0002309     2.20   0.028     .0000554    .0009606
****        of_op |   .0862736    .009768     8.83   0.000     .0671286    .1054186
****    distancia |  -.0001259   .0000177    -7.11   0.000    -.0001607   -.0000912
****        _cons |  -.1668645   .0697852    -2.39   0.017    -.3036411    -.030088
**** -------------+----------------------------------------------------------------
****      /athrho |   .2687221   .0996296     2.70   0.007     .0734516    .4639926
****     /lnsigma |   -1.41145   .0192509   -73.32   0.000    -1.449181   -1.373719
**** -------------+----------------------------------------------------------------
****          rho |   .2624353   .0927679                      .0733198    .4333327
****        sigma |   .2437896   .0046932                      .2347625    .2531638
****       lambda |    .063979   .0236295                      .0176661    .1102919
**** ------------------------------------------------------------------------------
**** LR test of indep. eqns. (rho = 0):   chi2(1) =     6.89   Prob > chi2 = 0.0087
**** ------------------------------------------------------------------------------


**** A continuaci�n haremos una breve explicaci�n sobre los resultados arrojados por el programa:

**** 3. Explicaci�n del programa

**** 3.1 N�mero de iteraciones

**** Como podemos observar en la parte superior, encontramos un resultado de la forma "Iteration 3:   log likelihood = -2637.1475"
**** Esto nos indica que se hicieron tres iteraciones para la estimaci�n y en la final, la funci�n de verosimilitud evaluada en los par�metros
**** toma el valor de 2637.1475.

**** 3.2 Metodolog�a de estimaci�n.

**** Inmediatamente despu�s del n�mero de iteraciones podemos ver la informaci�n "Treatment-effects model -- MLE". 
**** Esto nos indica que la estimaci�n se hizo por m�xima verosimilitud (Maximum Likelihood). Sin embargo, si se desea
**** tambi�n se puede llevar a cabo la estimaci�n en dos etapas utilizando la opci�n "twostep". A continuaci�n presentamos un ejemplo:


treatreg ha_nchs  $X, treat(D=$X of_op distancia) twostep

**** 3.3 Significancia global del modelo

**** En la parte superior derecha del resultado de la estimaci�n podemos ver "Wald chi2(4)    =      54.20" y "Prob > chi2     =     0.0000".
**** Esta es una prueba de significancia global del modelo. Lo que se hace es comparar 
**** el modelo donde se incluyen todas las variables con un modelo en el que s�lo hay una constante. Esta prueba es asint�ticamente
**** un test de Wald que se distribuye chi-cuadrado. En este caso el valor del estad�stico es 54.20 y por lo tanto 
**** el modelo es significativo a un nivel de 100% de confianza. Si se desea ver el modelo donde s�lo se incluye una constante
**** se debe seleccionar la opci�n "noskip".

treatreg ha_nchs  $X, treat(D=$X of_op distancia) noskip

**** Cuando se hace una estimaci�n mediante m�xima verosimilitud se comparan los valores de la funci�n de verosimilitud,
**** la metodolog�a es diferente cuando se hace en dos pasos. 

**** 3.4 Prueba de correlaci�n entre "u" y "v"

**** En la parte inferior del output podemos ver: "LR test of indep. eqns. (rho = 0):   chi2(1) =    6.89   Prob > chi2 = 0.0087"
**** Esta es una prueba de hip�tesis donde la hip�tesis nula es que la covarianza entre los t�rminos de error "u" y "v" es cero.
**** En este caso la prueba se rechaza, es decir hay una correlaci�n negativa entre los errores. Dado que el estad�stico de prueba
**** se construye utilizando los valores de la funci�n de verosimilitud en las dos ecuaciones, esta prueba solo se presenta cuando
**** la estimaci�n se hace por m�xima verosimilitud. 

**** 3.5 Sigma, rho y lambda

**** Sigma es el valor estimado de la varianza del t�rmino de error "u". rho es el coeficiente de correlaci�n entre 
**** los errores "u" y "v" y lambda es el producto entre sigma y rho.

**** 3.6 La raz�n de mills. 

**** Dado que, seg�n se presenta en la teor�a de funciones de control, una parte fundamental de la ecuaci�n a estimar es 
**** el inverso de la raz�n de mills para cada observaci�n, tambi�n llamada "hazard ratio", existe una opci�n que nos permite
**** predecir lo valores de esta funci�n para cada observaci�n:

treatreg ha_nchs  $X, treat(D=$X of_op distancia) hazard(var)











