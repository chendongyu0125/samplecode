#Packages installed:
  #1.Distributions (normal distribution)
  #2. Optim (optimization package)
  #3. NLopt (ooptimization package)

using Distributions
#using Optim
using NLopt
using JuMP
using Debug
using GLPKMathProgInterface


##Before everything, I will define the frac function



#Fraction function. It returns 0.5 in case both non-labor incomes are zero
function frac(Yf,Ym)
    if Yf==0 && Ym==0
        l=0.5
    else
        l=Yf./(Yf+Ym)
    end
  return(l)
end


#And now to the likelihood function
count = 0
function likelihood(x::Vector,Ef,Em,A,D,Sb0,Pga,Age,S,Hf,Hm,Inv,Wf,Wm,Yf,Ym,Barg,Hhchores,Fage,Mage,Fed,Med)
  #Block of parameter definition
  #_________________________________________________
  #1.a Utility of the father
  α1f=x[1]  #Consumption father
  α2f=x[2]  #Utility from child quality
  α3f=x[3]  #Disutility from work
  α40f=x[4] #Disutility from effort
  α41f=x[5] #Disutility from effort * someone helps

  #A.1.b Utility of the mother
  α1m=x[6] #Consumption father
  α2m=x[7]  #Utility from child quality
  α3m=x[8]  #Disutility from work
  α40m=x[9] #Disutility from effort
  α41m=x[10] #Disutility from effort * someone helps

  #2 Production of skills

  γ0=x[11]  #Investments
  γ1=x[12]   #Productivity of the father effort
  γ2=x[13]   #Productivity of the mother effort
  γ3=x[14]   #Productivity of childcare
  γ4=x[15]   #Discipline
  γ5=x[16]  #PG characteristics
  γ6=x[17]  #Birth conditions
  γ7=x[18]  #Age of child (in months)
  γ8=x[19]    #Intercept

  #3 Mincer equation
  #(this preliminary estimates come from ols estimates)

  #3.a Mincer for father
  β0f=x[20] #Intercept
  β1f=x[21] #Education
  β2f=x[22]  #Age
  β3f=x[23]  #Age^2


  #3.b Mincer for Mother
  β0m=x[24]  #Intercept
  β1m=x[25]  #Education
  β2m=x[26]  #Age
  β3m=x[27] #Age^2
  #4 Price of childcare
  Pa=x[28]




  #5 Variance terms
  VarepsS=exp(x[29])  #Skills
  VarepsEF=exp(x[30]) #Effort father
  VarepsEM=exp(x[31]) #Effort mother
  VarepsWF=exp(x[32]) #Wages father
  VarepsWM=exp(x[33]) #Wages mother
  VarepsF1=exp(x[34])
  VarepsF2=exp(x[35])
  VarepsF3=exp(x[36])
  VarepsF4=exp(x[37])
  VarepsF5=exp(x[38])
  VarepsF6=exp(x[39])
  VarepsM1=exp(x[40])
  VarepsM2=exp(x[41])
  VarepsM3=exp(x[42])
  VarepsM4=exp(x[43])
  VarepsM5=exp(x[44])
  VarepsM6=exp(x[45])


  #6 Bargaining power
  λ0=x[46]
  λ1=x[47]
  λ2=x[48]
  λ3=x[49]
  λ4=x[50]

  #7 Last variance terms
  VarMMU=exp(x[51])  #Pareto weight
  VarINV=exp(x[52])  #Investment decisions



  #3. Defining intermediate functions for the total likelihood
  #___________________________________________________________


  #3.1 Intermediate functions
  Κ1(μ)=μ.*α1f+(1-μ).*α1m
  Κ2(μ)=μ.*α2f+(1-μ).*α2m
  Κ3(μ)=μ.*α3f+(1-μ).*α3m
  IT(hf,hm,wf,wm,Yf,Ym)=(Yf+Ym+wf.*hf+wm.*hm)
  R(μ,α1f,α1m,γ0)=(μ.*α1f+Κ2(μ).*γ0+(1-μ).*α1m)
  α4f(hchores)=α40f-α41f.*hchores
  α4m(hchores)=α40m-α41m.*hchores

  #3.2 Endogenous functions
  Cf(μ,hf,hm,a,wf,wm,Yf,Ym)=max(((IT(hf,hm,wf,wm,Yf,Ym)-Pa.*a).*(μ.*α1f./R(μ,α1f,α1m,γ0))),1.0000e-100)
  Cm(μ,hf,hm,a,wf,wm,Yf,Ym)=max(((IT(hf,hm,wf,wm,Yf,Ym)-Pa.*a).*((1-μ).*α1m./R(μ,α1f,α1m,γ0))),1.0000e-100)
  I (μ,hf,hm,a,wf,wm,Yf,Ym)=max(((IT(hf,hm,wf,wm,Yf,Ym)-Pa.*a).*(Κ2(μ).*γ0./R(μ,α1f,α1m,γ0))),0)
  Eff(μ,hf,hchores)=max((Κ2(μ).*γ1./(μ.*(1+hf)))-α4f(hchores),0)
  Emf(μ,hm,hchores)=max((Κ2(μ).*γ2./(μ.*(1+hm)))-α4m(hchores),0)

  #3.3 Skills
  s(μ,hf,hm,hchores,a,d,s0,pg,agemonths,wf,wm,Yf,Ym)=(γ0.*log(max(I(μ,hf,hm,a,wf,wm,Yf,Ym),1.0e-15))+γ1.*Eff(μ,hf,hchores)+γ2.*Emf(μ,hf,hchores)+γ3.*a+γ4.*d+γ5.*s0+γ6.*pg+γ7.*agemonths+γ8)

  #3.4 Welfare function
  function WELF(μ,hf,hm,hchores,a,d,s0,pg,agemonths,wf,wm,Yf,Ym)
    μ.*(α1f.*log(Cf(μ,hf,hm,a,wf,wm,Yf,Ym))-
                α3f.*hf-(1+hf).*((Eff(μ,hf,hchores).^2/2)+α4f(hchores).*Eff(μ,hf,hchores)))+
      (1-μ).*(α1m.*log(Cm(μ,hf,hm,a,wf,wm,Yf,Ym))-
                α3m.*hm-(1+hm).*(((Emf(μ,hm,hchores).^2/2)+α4m(hchores).*Emf(μ,hm,hchores))))+
      Κ2(μ).*s(μ,hf,hm,hchores,a,d,s0,pg,agemonths,wf,wm,Yf,Ym)
  end


  #3.5
  #Predicted wages
  WMPRED(Mage,Medu)=β0m+β1m.*Medu+β2m.*Mage+β3m.*(Mage.^2)
  WFPRED(Fage,Fedu)=β0f+β1f.*Fedu+β2f.*Fage+β3f.*(Fage.^2)

  #3.6
  #Lambda predicted
  lpred(wf,wm,Yf,Ym,Fage,Mage,Fed,Med)=λ0+λ1.*(wf./wm)+λ2.*frac(Yf,Ym)+λ3.*(Fage-Mage)+λ4.*(Fed-Med)


  #_____________________________________
  #4. Final intermediate inputs

  #4.1 Size of the loop
  NN=size(Wf,1)

  #4.2 Vector of labor supply

  #Father's labor supply
  HM=[0;0;0.5;0.5;1;1;0;0;0.5;0.5;1;1;0;0;0.5;0.5;1;1]
  HF=[0;0;0;0;0;0;0.5;0.5;0.5;0.5;0.5;0.5;1;1;1;1;1;1]
  CHILD=[0;1;0;1;0;1;0;1;0;1;0;1;0;1;0;1;0;1]

  HCH=[HF HM CHILD]

  #4.3 The matrix of variances

  VHCT=[VarepsF1 VarepsM1;
      VarepsF4 VarepsM4;
      VarepsF1 VarepsM2;
      VarepsF4 VarepsM5;
      VarepsF1 VarepsM3;
      VarepsF4 VarepsM6;
      VarepsF2 VarepsM1;
      VarepsF5 VarepsM4;
      VarepsF2 VarepsM2;
      VarepsF5 VarepsM5;
      VarepsF2 VarepsM3;
      VarepsF5 VarepsM6;
      VarepsF3 VarepsM1;
      VarepsF6 VarepsM4;
      VarepsF3 VarepsM2;
      VarepsF6 VarepsM5;
      VarepsF3 VarepsM3;
      VarepsF6 VarepsM6]

  #4.4 Vector of equalities across observations
  Obsec=zeros(18)

  #4.5 And imposing bargaining maximum of zero or 0.001
  Barr=max(Barg,0.001)

  #5. Starting the iteration for the likelihood function
  #__________________________________________________


  loglik=0 #Starting it at zero
  for ii=1:NN

    #5.1 By default the likelihood of wages will be set to zero:
    llWF=0
    llWM=0
    #5.2 Identifying the decision taken by observation ii:
    obs=[Hf[ii] Hm[ii] A[ii]]

    for jj=1:18
      Obsec[jj]=1-float(obs==HCH[jj,:])
    end
    Obsec
    (zMin,l)=findmin(Obsec)

    #5.3. Taking the variance of the observed decision:
    VF=VHCT[l,1]
    VM=VHCT[l,2]

    #And now I have to consider the cases if spouses provide labor supply or not


    #--------------------------------------
    #5.4 Case 1: Both spouses supply labor:
    #--------------------------------------

    if Hf[ii]>0 && Hm[ii]>0
      #No need to predict the wages. The WPREDICTED values are observed
      WFPREDICTED=Wf[ii]
      WMPREDICTED=Wm[ii]
      temp1=WELF(Barr[ii],HCH[:,1],HCH[:,2],Hhchores[ii],HCH[:,3],D[ii],Sb0[ii],Pga[ii],Age[ii],Wf[ii],Wm[ii],Yf[ii],Ym[ii])
      temp2=WELF(Barr[ii],Hf[ii]  , Hm[ii] ,Hhchores[ii],A[ii]   ,D[ii],Sb0[ii],Pga[ii],Age[ii],Wf[ii],Wm[ii],Yf[ii],Ym[ii])
      var=(Barr[ii].^2).*(VF+VHCT[:,1])+((1-Barr[ii]).^2).*(VM+VHCT[:,2])
      dif=(temp1-temp2)./sqrt(var)
      arg=1-cdf(Normal(),dif)
      llH=log(Obsec.*arg)
      llH2=(1-float(isfinite (llH)))

      #Changing the inf to large negative numbers
      for kk=1:18
        if llH2[kk]==1
        llH[kk]=-100000000
        end
      end

      #And computing the likelihood of the behavioral model
      llH=sum(llH)

      #Likelihood of the observed wages

      llWF=pdf(Normal(),(log(Wf[ii])-WFPRED(Fage[ii],Fed[ii]))./(VarepsWF))
      llWM=pdf(Normal(),(log(Wm[ii])-WMPRED(Mage[ii],Med[ii]))./(VarepsWM))


      #---------------------------------------------------
      #Case 2. Now, if father doesn't work and mother does
      #---------------------------------------------------

    elseif Hf[ii]==0 && Hm[ii]>0

      #I have to predict the wages of the father
      WFPREDICTED=WFPRED(Fage[ii],Fed[ii])

      temp1=WELF(Barr[ii],HCH[:,1],HCH[:,2],Hhchores[ii],HCH[:,3],D[ii],Sb0[ii],Pga[ii],Age[ii],WFPREDICTED,Wm[ii],Yf[ii],Ym[ii])
      temp2=WELF(Barr[ii],Hf[ii]  , Hm[ii] ,Hhchores[ii],A[ii]   ,D[ii],Sb0[ii],Pga[ii],Age[ii],WFPREDICTED,Wm[ii],Yf[ii],Ym[ii])
      var=(Barr[ii].^2).*(VF+VHCT[:,1])+((1-Barr[ii]).^2).*(VM+VHCT[:,2])
      dif=(temp1-temp2)./sqrt(var)
      arg=1-cdf(Normal(),dif)
      llH=log(Obsec.*arg)
      llH2=(1-float(isfinite (llH)))

      #Changing the inf to large negative numbers
      for kk=1:18
        if llH2[kk]==1
        llH[kk]=-100000000
        end
      end

      #And computing the likelihood of the behavioral model
      llH=sum(llH)

      #Likelihood of the observed wages (ONLY ONE OBSERVED IS THE MOTHER)
      llWM=pdf(Normal(),(log(Wm[ii])-WMPRED(Mage[ii],Med[ii]))./(VarepsWM))


      #---------------------------------------------------
      #Case 3. Now, if mother doesn't work and father does
      #---------------------------------------------------

    elseif Hf[ii]>0 && Hm[ii]==0

      #I have to predict the wages of the father
      WMPREDICTED=WMPRED(Mage[ii],Med[ii])

      temp1=WELF(Barr[ii],HCH[:,1],HCH[:,2],Hhchores[ii],HCH[:,3],D[ii],Sb0[ii],Pga[ii],Age[ii],Wf[ii],WMPREDICTED,Yf[ii],Ym[ii])
      temp2=WELF(Barr[ii],Hf[ii]  , Hm[ii] ,Hhchores[ii],A[ii]   ,D[ii],Sb0[ii],Pga[ii],Age[ii],Wf[ii],WMPREDICTED,Yf[ii],Ym[ii])
      var=(Barr[ii].^2).*(VF+VHCT[:,1])+((1-Barr[ii]).^2).*(VM+VHCT[:,2])
      dif=(temp1-temp2)./sqrt(var)
      arg=1-cdf(Normal(),dif)
      llH=log(Obsec.*arg)
      llH2=(1-float(isfinite (llH)))

      #Changing the inf to large negative numbers
      for kk=1:18
        if llH2[kk]==1
        llH[kk]=-100000000
        end
      end

      #And computing the likelihood of the behavioral model
      llH=sum(llH)

      #Likelihood of the observed wages (ONLY ONE OBSERVED IS THE MOTHER)
      llWF=pdf(Normal(),(log(Wf[ii])-WFPRED(Fage[ii],Fed[ii]))./(VarepsWF))


      #---------------------------------------------------
      #Case 4. Now, if mother doesn't work and father does
      #---------------------------------------------------
    elseif Hf[ii]==0 && Hm[ii]==0

      #I have to predict the wages of the father and the mother
      WFPREDICTED=WFPRED(Fage[ii],Fed[ii])
      WMPREDICTED=WMPRED(Mage[ii],Med[ii])

      temp1=WELF(Barr[ii],HCH[:,1],HCH[:,2],Hhchores[ii],HCH[:,3],D[ii],Sb0[ii],Pga[ii],Age[ii],WFPREDICTED,WMPREDICTED,Yf[ii],Ym[ii])
      temp2=WELF(Barr[ii],Hf[ii]  , Hm[ii] ,Hhchores[ii],A[ii]   ,D[ii],Sb0[ii],Pga[ii],Age[ii],WFPREDICTED,WMPREDICTED,Yf[ii],Ym[ii])
      var=(Barr[ii].^2).*(VF+VHCT[:,1])+((1-Barr[ii]).^2).*(VM+VHCT[:,2])
      dif=(temp1-temp2)./sqrt(var)
      arg=1-cdf(Normal(),dif)
      llH=log(Obsec.*arg)
      llH2=(1-float(isfinite (llH)))

      #Changing the inf to large negative numbers
      for kk=1:18
        if llH2[kk]==1
        llH[kk]=-100000000
        end
      end


      #And computing the likelihood of the behavioral model
      llH=sum(llH)
    end #End conditionals on labor supply of people

    #------------------------
    #5.5 Likelihood of skills
    #------------------------
    llS=pdf(Normal(),(s(Barr[ii],Hf[ii],Hm[ii],Hhchores[ii],A[ii],D[ii],Sb0[ii],Pga[ii],Age[ii],WFPREDICTED,WMPREDICTED,Yf[ii],Ym[ii])-S[ii])./(VarepsS))
    llS=max(log((1./VarepsS).*llS),-100)

    #------------------------
    #5.6 Likelihood of effort
    #------------------------
    llEF=pdf(Normal(),(Eff(Barr[ii],Hf[ii],Hhchores[ii])-Ef[ii])./VarepsEF)
    llEF=log((1./VarepsEF).*llEF)

    llEM=pdf(Normal(),(Emf(Barr[ii],Hm[ii],Hhchores[ii])-Em[ii])./VarepsEM)
    llEM=log((1./VarepsEM).*llEM)

    #----------------------------------------------
    #5.7 Likelihood of investment measurement error
    #----------------------------------------------

    llINV=pdf(Normal(),(I(Barr[ii],Hf[ii],Hm[ii],A[ii],WFPREDICTED,WMPREDICTED,Yf[ii],Ym[ii])-Inv[ii])./VarINV)
    llINV=log((1./VarINV).*llINV)


    #-------------------------------
    #5.8 Likelihood of Pareto weight
    #-------------------------------
    llMMU=pdf(Normal(),(log(Barr[ii]./(1-Barr[ii]))-lpred(WFPREDICTED,WMPREDICTED,Yf[ii],Ym[ii],Fage[ii],Mage[ii],Fed[ii],Med[ii]))/VarMMU)
    llMMU=max(log((1./VarMMU).*llMMU),-100)

    #-------------------------------------------
    #5.9 Computing the total likelihood function
    #-------------------------------------------
    loglik=llMMU+llINV+llEM+llEF+llS+llH+llWM+llWF
    loglik=-loglik



  end #End loop over observations
  global count
  count::Int += 1
  println("f_$count($x)")
  print(loglik)
  return(loglik)
end #End the function


#_____________________________________________________
#Load the data, install packages and organize the data

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# I. First step: Getting the data:
#1. See current directory
pwd()

#2. Set current directory to the julia in Chile
cd("/Users/rodrigoazuero/Documents/Research/Chile/julia")

#3. Read the data
M = readcsv("matlabexport.csv")

#4. Name the vectors correspondingly
A=M[:,1]
Sb0=M[:,2]
Inv=M[:,3]
D=M[:,4]
Hhchores=M[:,5]
S=M[:,6]
Barg=M[:,7]
Age=M[:,8]
Mpart=M[:,9]
Wm=M[:,10]
Wf=M[:,11]
Yf=M[:,12]
Ym=M[:,13]
Em=M[:,14]
Ef=M[:,15]
Pga=M[:,16]
Hf=M[:,17]
Hm=M[:,18]
Fage=M[:,19]
Mage=M[:,20]
Med=M[:,21]
Fed=M[:,22]
Folio=M[:,23]

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#Impute wages on the non workers and get rid of N

β0f= 2.481863 #Intercept
β1f=.1156016  #Education
β2f=.0785408  #Age
β3f=-.0010014 #%Age^2
β0m=1.831828  #Intercept
β1m=.1156016   #Education
β2m=.0785408   #Age
β3m=-.0010014  #Age^2


Fwpred=exp(β0f+β1f.*Fed+β2f.*Fage+β3f*Fage.^2)
Mwpred=exp(β0m+β1m.*Med+β2m.*Mage+β3m*Mage.^2)


#Wherever I got a NaN I will change it to zeros.
Fwpred[isnan(Fwpred)]=0
Mwpred[isnan(Mwpred)]=0
Wf[isnan(Wf)]=0
Wm[isnan(Wm)]=0

#And defining offered wage as the predicted for those that are not working
#or as the observed for those who work.
endLoop=size(Wf,1)


Wofff=zeros(endLoop)
Woffm=zeros(endLoop)

#I will do it as a loop
for ii=1:endLoop
  #For men
  if Hf[ii]==0
    Wofff[ii]=Fwpred[ii]
  elseif Hf[ii]>0
    Wofff[ii]=Wf[ii]
  end

  #For women
  if Hm[ii]==0
    Woffm[ii]=Mwpred[ii]
  elseif Hm[ii]>0
    Woffm[ii]=Wm[ii]
  end

  #In this section of the loop replace zeros by NaN
  if Wofff[ii]==0
    Wofff[ii]=NaN
  end
  if Woffm[ii]==0
    Woffm[ii]=NaN
  end
end

DAT=hcat(Wofff,Woffm,Yf,Ym,Sb0,D,Hhchores,S,Barg,Age,Hf,Hm,A,Em,Ef,Fage,Mage,Pga,Fed,Med,Inv,Folio)
#Removing all the NaN:

#First, identifying which rows have NaN:
Naninfo=!any(isnan(DAT),2)

#Generating array with the not NaN info:
endLoop=size(DAT,1)
colN=size(DAT,2)
DATNONA=zeros(endLoop,colN)
jj=1
ii=1
for ii=1:endLoop
  if Naninfo[ii]==true
    DATNONA[jj,:]=DAT[ii,:]
    jj=jj+1
  end
end

#And storing the final data ii
DATOTAL=zeros(jj-1,colN)
for ii=1:jj-1
  DATOTAL[ii,:]=DATNONA[ii,:]
end

DATOTAL
Dat=DATOTAL

#Re generating the variables
Wf=Dat[:,1]
Wm=Dat[:,2]
Yf=Dat[:,3]
Ym=Dat[:,4]
Sb0=Dat[:,5]
D=Dat[:,6]
Hhchores=Dat[:,7]
S=Dat[:,8]
Barg=Dat[:,9]
Age=Dat[:,10]
Hf=Dat[:,11]
Hm=Dat[:,12]
A=Dat[:,13]
Em=Dat[:,14]
Ef=Dat[:,15]
Fage=Dat[:,16]
Mage=Dat[:,17]
Pga=Dat[:,18]
Fed=Dat[:,19]
Med=Dat[:,20]
Inv=Dat[:,21]
Folio=Dat[:,22]

#And then turning back zeros to NaN
endLoop=size(Wf,1)
for ii=1:endLoop
  if Wf[ii]==0
    Wf[ii]=NaN
  end
  if Wm[ii]==0
    Wm[ii]=NaN
  end
end

#__________________________________________





#_________________________________________________________
# Block of parameter definitions for the initial value:
#_________________________________________________

#III .1.a Utility of the father
α1f=0.5  #Consumption father
α2f=2.6  #Utility from child quality
α3f=1.1  #Disutility from work
α40f=0.35 #Disutility from effort
α41f=0.01 #Disutility from effort * someone helps

#A.1.b Utility of the mother
α1m=0.5 #Consumption father
α2m=3.1;  #Utility from child quality
α3m=2.0;  #Disutility from work
α40m=0.21; #Disutility from effort
α41m=0.01; #Disutility from effort * someone helps

# III.2 Production of skills

γ0=0.05  #Investments
γ1=0.1   #Productivity of the father effort
γ2=0.1   #Productivity of the mother effort
γ3=0.2   #Productivity of childcare
γ4=0.0   #Discipline
γ5=0.07  #PG characteristics
γ6=0.01  #Birth conditions
γ7=0.07  #Age of child (in months)
γ8=-4    #Intercept

#III.3 Mincer equation
#(this preliminary estimates come from ols estimates)

#III.3.a Mincer for father
β0f= 2.481863 #Intercept
β1f=.1156016 #Education
β2f=.0785408  #Age
β3f=-.0010014  #Age^2


#III.3.b Mincer for Mother
β0m=1.974835  #Intercept
β1m=.1156016  #Education
β2m=.0785408  #Age
β3m=-.0010014 #Age^2

#III.4 Price of childcare
Pa=0.1

#III.5 Bargaining power
λ0=-5.0
λ1=0.1
λ2=0.01
λ3=0.01
λ4=0.01


#III.6 Variance terms
VarepsS=1.0  #Skills
VarepsEF=1.0 #Effort father
VarepsEM=1.0 #Effort mother
VarepsWF=1.0 #Wages father
VarepsWM=1.0 #Wages mother
VarepsF1=1.0
VarepsF2=1.0
VarepsF3=1.0
VarepsF4=1.0
VarepsF5=1.0
VarepsF6=1.0
VarepsM1=1.0
VarepsM2=1.0
VarepsM3=1.0
VarepsM4=1.0
VarepsM5=1.0
VarepsM6=1.0
VarMMU=1.0  #1.0areto weight
VarINV=1.0  #Investment decisions
#And setting everything into one vector
xo=[α1f α2f α3f α40f α41f α1m α2m α3m α40m α41m γ0 γ1 γ2 γ3 γ4 γ5 γ6 γ7 γ8 β0f β1f β2f β3f β0m β1m β2m β3m Pa VarepsS VarepsEF VarepsEM VarepsWF VarepsWM VarepsF1 VarepsF2 VarepsF3 VarepsF4 VarepsF5 VarepsF6 VarepsM1 VarepsM2 VarepsM3 VarepsM4 VarepsM5 VarepsM6 λ0 λ1 λ2 λ3 λ4 VarMMU VarINV]

xo=vec(xo)

likelihood(xo,Ef,Em,A,D,Sb0,Pga,Age,S,Hf,Hm,Inv,Wf,Wm,Yf,Ym,Barg,Hhchores,Fage,Mage,Fed,Med)

function evlik(x::Vector,m)
  a=likelihood(x,Ef,Em,A,D,Sb0,Pga,Age,S,Hf,Hm,Inv,Wf,Wm,Yf,Ym,Barg,Hhchores,Fage,Mage,Fed,Med)
end
 evlik(xo)
evlik
#Attempt with NlOpt
count = 0 # keep track of # function evaluations

#Can't use ld_MMA algorithm because is gradient based. Error in evlik as it will demand the second vector
opt = Opt(:LN_PRAXIS, 52)
xtol_rel!(opt,1e-4)
xo=vec(xo)
evlik(vec(xo))
min_objective!(opt, evlik)
(minf,minx,ret) = optimize(opt,xo)
(minf,minx,ret) = optimize(opt,[α1f, α2f, α3f, α40f, α41f, α1m, α2m, α3m, α40m, α41m, γ0, γ1, γ2, γ3, γ4, γ5, γ6, γ7, γ8, β0f, β1f, β2f, β3f, β0m, β1m, β2m, β3m, Pa, VarepsS, VarepsEF, VarepsEM, VarepsWF, VarepsWM, VarepsF1, VarepsF2, VarepsF3, VarepsF4, VarepsF5, VarepsF6, VarepsM1, VarepsM2, VarepsM3, VarepsM4, VarepsM5, VarepsM6, λ0, λ1, λ2, λ3, λ4, VarMMU, VarINV] )



opt
min_objective!
println("got $minf at $minx after $count iterations (returned $ret)")


opt



#Attempt with JuMP

m=Model()
m=Model(solver =NLoptSolver(algorithm=:LN_PRAXIS))
m=Model(solver =GLPKSolverLP(method=:Exact, presolve=true))
GLPKSolverLP(method=:Exact, presolve=true)
m=Model(solver =SCS())


#Setting the initial values
@defVar(m,x1,start=α1f)
@defVar(m,x2,start=α2f)
@defVar(m,x3,start=α3f)
@defVar(m,x4,start=α40f)
@defVar(m,x5,start=α41f)
@defVar(m,x6,start=α1m)
@defVar(m,x7,start=α2m)
@defVar(m,x8,start=α3m)
@defVar(m,x9,start=α40m)
@defVar(m,x10,start=α41m)
@defVar(m,x11,start=γ0)
@defVar(m,x12,start=γ1)
@defVar(m,x13,start=γ2)
@defVar(m,x14,start=γ3)
@defVar(m,x15,start=γ4)
@defVar(m,x16,start=γ5)
@defVar(m,x17,start=γ6)
@defVar(m,x18,start=γ7)
@defVar(m,x19,start=γ8)
@defVar(m,x20,start=β0f)
@defVar(m,x21,start=β1f)
@defVar(m,x22,start=β2f)
@defVar(m,x23,start=β3f)
@defVar(m,x24,start=β0m)
@defVar(m,x25,start=β1m)
@defVar(m,x26,start=β2m)
@defVar(m,x27,start=β3m)
@defVar(m,x28,start=Pa)
@defVar(m,x29,start=VarepsS)
@defVar(m,x30,start=VarepsEF)
@defVar(m,x31,start=VarepsEM)
@defVar(m,x32,start=VarepsWF)
@defVar(m,x33,start=VarepsWM)
@defVar(m,x34,start=VarepsF1)
@defVar(m,x35,start=VarepsF2)
@defVar(m,x36,start=VarepsF3)
@defVar(m,x37,start=VarepsF4)
@defVar(m,x38,start=VarepsF5)
@defVar(m,x39,start=VarepsF6)
@defVar(m,x40,start=VarepsM1)
@defVar(m,x41,start=VarepsM2)
@defVar(m,x42,start=VarepsM3)
@defVar(m,x43,start=VarepsM4)
@defVar(m,x44,start=VarepsM5)
@defVar(m,x45,start=VarepsM6)
@defVar(m,x46,start=λ0)
@defVar(m,x47,start=λ1)
@defVar(m,x48,start=λ2)
@defVar(m,x49,start=λ3)
@defVar(m,x50,start=VarMMU)
@defVar(m,x51,start=VarINV)



@setNLObjective(m, Min, evlik)
@setNLObjective(m, Min, evlik(x1,x2,x4,x5,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51))



@setNLObjective(m, Min, )
using MathProgBase

linprog(c, A, sense, b, solver) = linprog(c, A, sense, b, 0, Inf, solver)

using(KNITRO)
print("a")
evlik(xo)

m
## Setting everything for the solver
using(SCS)
@debug
status=solve(m) #User defined functions cannot be used with JuMP!

MathProgBase.numvar(m::Model)
MathProgBase.numlinconstr(m::Model)
MathProgBase.getsolvetime(m::Model)





xo
(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)=optimize(evlik,e,method=:cg)
optimize(evlik,xo)
e
optimize(evlik,e,method=:cg)
Results=optimize(evlik,e,method=:simulated_annealing)


optimize(evlik,[α1f α2f α3f α40f α41f α1m α2m α3m α40m α41m γ0 γ1 γ2 γ3 γ4 γ5 γ6 γ7 γ8 β0f β1f β2f β3f β0m β1m β2m β3m Pa VarepsS VarepsEF VarepsEM VarepsWF VarepsWM VarepsF1 VarepsF2 VarepsF3 VarepsF4 VarepsF5 VarepsF6 VarepsM1 VarepsM2 VarepsM3 VarepsM4 VarepsM5 VarepsM6 λ0 λ1 λ2 λ3 λ4 VarMMU VarINV])


optimize(evlik,[α1f, α2f, α3f, α40f, α41f, α1m, α2m, α3m, α40m, α41m, γ0, γ1, γ2, γ3, γ4, γ5, γ6, γ7, γ8, β0f, β1f, β2f, β3f, β0m, β1m, β2m, β3m, Pa, VarepsS, VarepsEF, VarepsEM, VarepsWF, VarepsWM, VarepsF1, VarepsF2, VarepsF3, VarepsF4, VarepsF5, VarepsF6, VarepsM1, VarepsM2, VarepsM3, VarepsM4, VarepsM5, VarepsM6, λ0, λ1, λ2, λ3, λ4, VarMMU, VarINV], method = :l_bfgs)

Results=optimize(evlik,e,method = :nelder_mead,iterations=500000000)
Results

use(JuMP)

?
  optimize
using(Optim)


function sqerror(betas)
    err = 0.0
    for i in 1:length(x)
        pred_i = betas[1] + betas[2] * x[i]
        err += (y[i] - pred_i)^2
    end
    return err
end











count = 0 # keep track of # function evaluations

function myfunc(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end

    global count
    count::Int += 1
    println("f_$count($x)")

    sqrt(x[2])
end

function myconstraint(x::Vector, grad::Vector, a, b)
    if length(grad) > 0
        grad[1] = 3a * (a*x[1] + b)^2
        grad[2] = -1
    end
    (a*x[1] + b)^3 - x[2]
end

opt = Opt(:LD_MMA, 2)
lower_bounds!(opt, [-Inf, 0.])
xtol_rel!(opt,1e-4)

min_objective!(opt, myfunc)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,2,0), 1e-8)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)

xa=[1.234, 3.3]
(minf,minx,ret) = optimize(opt, xa)
println("got $minf at $minx after $count iterations (returned $ret)")







count = 0 # keep track of # function evaluations

function myfunc(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end

    global count
    count::Int += 1
    println("f_$count($x)")

    sqrt(x[2])
end


opt = Opt(:LD_MMA, 2)
lower_bounds!(opt, [-Inf, 0.])
xtol_rel!(opt,1e-4)

min_objective!(opt, myfunc)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,2,0), 1e-8)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)

xa=[1.234, 3.3]
(minf,minx,ret) = optimize(opt, xa)
println("got $minf at $minx after $count iterations (returned $ret)")
