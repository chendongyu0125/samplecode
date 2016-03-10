//
//  main.cpp
//  BEHAVIORAL38
//
//  Created by Rodrigo Azuero on 2/24/16.
//  Copyright (c) 2016 Rodrigo Azuero Melo. All rights reserved.
//

#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <unistd.h>
#include <nlopt.hpp>
using std::vector;
using namespace std;

// Functions: Start with `F_' and then underscore continues
// Data inputs: All start with `in'
// Greek letters: start with double letters. e.g. aalpha, bbeta and ggamma.
// Likelihood functions: After the usual `F_' for functions we use `like_' to denote them.

//The structure of the program is as follows: The first block will contain the predicting equations. This are endogenous variables
//determined within the model or data that depends on parameters of the model (Wages, Skills, Effort, Consumption). The second part
//will include the likelihood computation for each part (i.e. likelihood of observing such wages, likelihood of observing skills, not
//necessarily of consumption but indeed of effort levels given the observed level). The third part includes the computation of all the likelihood
//function

//This file computes the likelihood as specified in the folder "Likelihoodfunction"


//===================================================================================//
//-1. Before starting with the function definitions I will define the transformations
//that are necessary for each parameter. I will not export them to R and they will all
//start with FT to denote Function and Transformation.
//==================================================================================//

//This function generates a discrete distribution.
//It is used exclusively as an example for how to draw from a
//discrete distribution so that later we get draws for the
//bootstrap filter.

double DISCDIST(double a){
    double resto=a*a;
    //NOW ATEMPT WITH DOUBLE VECTOR
    std::vector<double> probs(5);
    probs[0]=1;
    probs[1]=1;
    probs[2]=1;
    probs[3]=1;
    probs[4]=10;
    boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
    //boost::random::discrete_distribution<> dist(probs);
    boost::mt19937 eng(std::time(0));
    std::cout << dist(eng);
    return(resto);
}

//Define template for random number generator
template<class T>
double gen_normal_3(T &generator)
{
    return generator();
}



//-1.0 Theta0 transformation
//--------------------------
// [[Rcpp::export]]
long double FT_ttheta0(double ttheta0, double ttheta1){
    //If ttheta0 too large we need to fix it
    if (ttheta0>700){
        ttheta0=700;
    }
    double result=exp(ttheta0)/(1+exp(ttheta0)+exp(ttheta1));
    if (result>0.999){
        result=0.99;
    }
    if (result<1.0e-5){
        result=1.0e-5;
    }
    return(result);
}

//-1.1 Ttheta1 transformation
//-------------------------
// [[Rcpp::export]]
long double FT_ttheta1(double ttheta0, double ttheta1){
    if (ttheta1>700){
        ttheta1=700;
    }
    
    double result=exp(ttheta1)/(1+exp(ttheta0)+exp(ttheta1));
    if (result>0.999){
        result=0.99;
    }
    if (result<1.0e-5){
        result=1.0e-5;
    }
    
    return(result);
}

//-1.2 Pphi transformation
//Function checked!
//-------------------------
// [[Rcpp::export]]
long double FT_pphi(double pphi){
    double result=1-exp(-pphi);
    // Fix pphi if it is zero
    if (result==0){
        result=1.0e-20;
    }
    
    // Fix pphi if it is too large
    if (pphi>4){
        result=1-exp(-4);
    }
    
    //Fix pphi if it is too small
    else if (pphi<-2.5){
        result = 1-exp(2.5);
    }
    return(result);
}


//-1.3 Exp transformation
//-------------------------
// [[Rcpp::export]]
long double FT_exp(double logstd){
    double R=exp(logstd);
    if (R==0){
        R=1.0e-17;
    }
    //0.0.1 We also need to take into account the fact that R can be too large
    if (logstd>700){
        R=exp(700);
    }
    return(R);
}


//-1.4 Gf transformation
//-----------------------------------------------------------------
// [[Rcpp::export]]
long double FT_gf(double gf){
    double result=gf;
    //We need to arrange this in case ggammaf or ggammam are small
    if (gf>0.999){
        result=0.99;
    }
    if (gf<1.0e-5){
        result=1.0e-5;
    }
    return(result);
}
//-1.5// aalphas. We need constraints in aalphas. They should be positive and sum up to one
//moreover, we need aalpha4=aalpha40+aalpha41 to be positive because otherwise we might get
//into trouble when trying to find the adequate level of effort for single mothers. The normalization of making them add
//to one needs to be done in place as we need all three arguments to be normalized and divide them by the sum of
//everyhing.
//-1.5.0 aalpha40
//------------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha40(double aalpha40){
    double a4m=1/(1+exp(aalpha40));
    double result=a4m;
    if (a4m<1.0e-5){
        result=1.0e-5;
    }
    if (a4m>0.9999){
        result=0.999;
    }
    return(result);
}
//-1.5.1 aalpha41
//------------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha41(double aalpha40, double aalpha41){
    double a40=FT_aalpha40(aalpha40);
    double a41=1/(1+exp(aalpha41));
    double result=a41*a40;
    return(result);
}

//-1.5.2 aalpha1
//-----------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha1(double aalpha1m){
    double a1m=1/(1+exp(aalpha1m));
    double result=a1m;
    if (a1m<1.0e-5){
        result=1.0e-5;
    }
    if (a1m>0.9999){
        result=0.999;
    }
    return(result);
}


//-1.5.2 aalpha2
//---------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha2(double aalpha2m){
    double a2m=1/(1+exp(aalpha2m));
    double result=a2m;
    if (a2m<1.0e-5){
        result=1.0e-5;
    }
    if (a2m>0.9999){
        result=0.999;
    }
    return(result);
}

//-1.5. aalpha3
//---------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha3(double aalpha3m){
    double a3m=1/(1+exp(aalpha3m));
    double result=a3m;
    if (a3m<1.0e-5){
        result=1.0e-5;
    }
    if (a3m>0.9999){
        result=0.999;
    }
    return(result);
}
//-1.6. Bounding transformations. This will serve
//the purpose of putting close to zero
//parameters that are really small or putting
//close to one other ones. Such is the case of the aalpha parameters
//----------------------------------------------------------------------
// [[Rcpp::export]]
long double FT_approach(double aalpha){
    double result=aalpha;
    if (aalpha<1.0e-5){
        result=1.0e-5;
    }
    if (aalpha>0.9999){
        result=0.999;
    }
    return(result);
}

//===========================================
//Second transformations of the aalphas:
//==========================================
//

long double FT_aalpha_trans(double aalpha1, double aalpha2, double aalpha3, double aalpha40, double aalpha41){
    if (aalpha1>700){
        aalpha1=700;
    }
    double result=exp(aalpha1)/(exp(aalpha1)+exp(aalpha2)+exp(aalpha3)+exp(aalpha40)+exp(aalpha41));
    return(result);
}




//===================================================================================//
//0. Block of function definition
//===================================================================================//




//----------
//0.0. Wages
//------------
//===================================================================================//
// [[Rcpp::export]]
long double F_predwage(double bbeta0, double bbeta1, double bbeta2, double bbeta3,  double Schooling, double Age){
    //1. Defining the output
    double wage=bbeta0+bbeta1*Schooling+bbeta2*Age+bbeta3*Age*Age;
    //1.1 If wage>700 we need to impose the max numerically possible
    if (wage>700){
        wage=exp(700);
    }
    //1.2
    //If such is not the case, we can define it as exp of the predicted value
    else{
        wage=exp(wage);
    }
    return(wage);
}





//---------------------------
//0.1 Production of Skills
//---------------------------
//===================================================================================//
// [[Rcpp::export]]
long double F_predskills(double ddelta0, double ddelta1, double ddelta2, double ddelta3,
                         double ddelta4,
                         double Age,  double ttheta0, double ttheta1, double ttheta2,
                         double pphi, double gf, double gm, double Fe, double Me,
                         double I, double S0,
                         double childcare, double PG, double Hmembers){
    
    //0 Compute the level of effort in the CES production of Effort
    //0.0 If effort levels are zero, replace them by something close to zero.
    if (Fe<=0){
        Fe=0.01;
    }
    if (Me<=0){
        Me=0.01;
    }
    
    
    //And compute the effort by the CES production function
    double e=gf*pow(Fe,pphi)+gm*pow(Me,pphi);
    e=pow(e,1/pphi);
    if (e<=0){
        e=1.0e-5;
    }
    
    double R=exp(ddelta0+ddelta1*Age+ddelta2*childcare+ddelta3*PG+ddelta4*Hmembers);
    
    
    //1 If investments are zero or negative, put them close to zero
    if (I<=0){
        I=0.01;
    }
    //2 The production of skills
    
    double skills=R*pow(S0,ttheta0)*pow(I,ttheta1)*pow(e,ttheta2);
    double isf=std::isfinite(skills);
    if (isf==0){
        skills=5.0e+10;
    }
    return(skills);
}

//----------------------------------------------
//0.2 Effort
// This function is defined in 09-02-2015, pg7.
//------------------------------------------------

//We will have to define several intermediate functions


//---------------------------------------------------
//0.2.0: Tau_j. This is an intermediate function defined
//in the same page as mentioned before
//---------------------------------------------------
//===================================================================================//
//0.2.0.0: tau_f
//===================================================================================//
// [[Rcpp::export]]
long double F_tau_f(double ggammaf, double aalpha4, double mmu, double hf){
    double tau_f=aalpha4*(mmu)*(1+hf);
    if (tau_f==0){
        tau_f=1.0e-10;
    }
    tau_f=ggammaf/tau_f;
    return(tau_f);
}
//===================================================================================//
//0.2.0.1: tau_m
//===================================================================================//
// [[Rcpp::export]]
long double F_tau_m(double ggammam, double aalpha4, double mmu, double hm){
    double tau_m=aalpha4*(1-mmu)*(1+hm);
    if (tau_m==0){
        tau_m=1.0e-10;
    }
    tau_m=ggammam/tau_m;
    return(tau_m);
}

//---------------------------------------------------
//0.2.1. ces_tau as defined in page 7
//----------------------------------------------------

//===================================================================================//
//0.2.1.0 ces_tau_f
//===================================================================================//
// [[Rcpp::export]]
long double F_ces_tau_f(double tau_f, double tau_m, double ggammaf,  double pphi){
    //0. Define the denominator of the  exponent
    double denexp=1-pphi;
    double ggammam=1-ggammaf;
    //0.0 If the denominator is infinite we need to Fix
    if (pphi==1){
        denexp=1.0e-10;
    }
    //0.1 Defining the exponent
    double exponent=pphi/denexp;
    if (exponent==0){
        exponent=1.0e-5;
    }
    //1. Now define the denominator
    double denominator=ggammaf*pow(tau_f,exponent)+ggammam*pow(tau_m,exponent);
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.5e+10){
        denominator=1.5e+10;
    }
    //2. Defining the numerator
    double numerator=pow(tau_f,exponent);
    if (numerator>=1.5e+10){
        numerator=1.5e+10;
    }
    //3. Final definition of the ces_tau_f Functions
    double ces_tau_f=numerator/denominator;
    return(ces_tau_f);
}
//===================================================================================//
//0.2.1.1 ces_tau_m
//===================================================================================//
// [[Rcpp::export]]
long double F_ces_tau_m(double tau_f, double tau_m, double ggammaf,  double pphi){
    //0. Define the denominator of the  exponent
    double denexp=1-pphi;
    double ggammam=1-ggammaf;
    //0.0 If the denominator is infinite we need to Fix
    if (pphi==1){
        denexp=1.0e-10;
    }
    //0.1 Defining the exponent
    double exponent=pphi/denexp;
    if (exponent==0){
        exponent=1.0e-5;
    }
    
    //1. Now define the denominator
    double denominator=ggammaf*pow(tau_f,exponent)+ggammam*pow(tau_m,exponent);
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.5e+10){
        denominator=1.5e+10;
    }
    //2. Defining the numerator
    double numerator=pow(tau_m,exponent);
    if (numerator>=1.5e+10){
        numerator=1.5e+10;
    }
    //3. Final definition of the ces_tau_f Functions
    double ces_tau_m=numerator/denominator;
    return(ces_tau_m);
}


//---------------------------------------------------
//0.2.2. kkappa2 weighted utility of skills
//----------------------------------------------------
//===================================================================================//
// [[Rcpp::export]]
long double F_kkappa2(double mmu, double aalpha2f, double aalpha2m){
    double kkappa2=aalpha2f*mmu+aalpha2m*(1-mmu);
    return(kkappa2);
}

//--------------------------------------------
//0.2.3 Effort levels
//--------------------------------------------
//===================================================================================//
//0.2.3.0 Effort of father
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_f(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,  double ttheta2){
    //0. Defining the tau inputs
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_m( tau_f,  tau_m,  ggammaf,  pphi);
    
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    
    //3. Finally we have the total effort
    double effort_f=tau_f*kkappa2*ttheta2*ces_tau_m;
    if (effort_f>5.0e+10){
        effort_f=5.0e+10;
    }
    if (effort_f<=0){
        effort_f=1.0e-10;
    }
    return(effort_f);
}
//===================================================================================//
//0.2.3.1 Effort of mother
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_m(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,  double ttheta2){
    //0. Defining the tau inputs
    
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_f( tau_f,  tau_m,  ggammaf,  pphi);
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    //3. Finally we have the total effort
    double effort_m=tau_m*kkappa2*ttheta2*ces_tau_m;
    
    if (effort_m>5.0e+10){
        effort_m=5.0e+10;
    }
    if (effort_m<=0){
        effort_m=1.0e-10;
    }
    return(effort_m);
    
    
}
//======================================
//Effort of father in the first period
//===================================================================================//
//0.2.3.0 Effort of father
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_f10(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,
                         double ttheta0,  double ttheta2, double BBETA){
    //0. Defining the tau inputs
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_m( tau_f,  tau_m,  ggammaf,  pphi);
    
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    
    //3. Finally we have the total effort
    double effort_f=tau_f*(kkappa2*ttheta2+BBETA*kkappa2*ttheta2*ttheta0)*ces_tau_m;
    if (effort_f>5.0e+10){
        effort_f=5.0e+10;
    }
    if (effort_f<=0){
        effort_f=1.0e-10;
    }
    return(effort_f);
}
//===================================================================================//
//0.2.3.1 Effort of mother
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_m10(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,
                         double ttheta0,  double ttheta2, double BBETA){
    //0. Defining the tau inputs
    
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_f( tau_f,  tau_m,  ggammaf,  pphi);
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    //3. Finally we have the total effort
    double effort_m=tau_m*((kkappa2*ttheta2+BBETA*kkappa2*ttheta2*ttheta0))*ces_tau_m;
    
    if (effort_m>5.0e+10){
        effort_m=5.0e+10;
    }
    if (effort_m<=0){
        effort_m=1.0e-10;
    }
    return(effort_m);
    
    
}



//----------------------------------------------------------------
//0.2.3.2 Effort of person if single parent in 2012
//----------------------------------------------------------------
//
//
//===============================================================
//Note that in this case I will consider for single mothers simply
//to react to the action of the father. In this point
//I needed to install the gsl library.
//==============================================================
//-----------------------------------------------------

//===================================================================================//
//Defining the structure of parameters for the effort solver
//===================================================================================//
struct em_my_params {double EM_aalpha2m;  double EM_aalpha4m; double EM_ggammaf;  double EM_ef; double EM_pphi; double EM_ttheta; double EM_hm;};
//Defining the function to which we will define the roots
double F_root_em (double em, void *params){
    //0. Assigning the parameters
    struct em_my_params *p=(struct em_my_params *) params;
    double aalpha2m=p->EM_aalpha2m;
    double aalpha4m=p->EM_aalpha4m;
    double ggammaf=p->EM_ggammaf;
    double ggammam=1-ggammaf;
    double ef=p->EM_ef;
    double pphi=p->EM_pphi ;
    double ttheta=p->EM_ttheta ;
    double hm=p->EM_hm;
    //1. Getting the function
    double residual=aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+ggammam*pow(em,pphi))-ttheta*aalpha2m*ggammam*pow(em,pphi-1);
    
    return(residual);
}

//===================================================================================//
//Defining the structure of parameters for the effort solver in 2010
//===================================================================================//
struct em_my_params10 {double EM_aalpha2m;  double EM_aalpha2m10; double EM_aalpha4m; double EM_ggammaf;  double EM_ef; double EM_pphi; double EM_ttheta0; double EM_ttheta2; double EM_hm;};
//Defining the function to which we will define the roots
double F_root_em10 (double em, void *params){
    //0. Assigning the parameters
    struct em_my_params10 *p=(struct em_my_params10 *) params;
    double aalpha2m=p->EM_aalpha2m;
    double aalpha2m10=p->EM_aalpha2m10;
    double aalpha4m=p->EM_aalpha4m;
    double ggammaf=p->EM_ggammaf;
    double ggammam=1-ggammaf;
    double ef=p->EM_ef;
    double pphi=p->EM_pphi ;
    double ttheta0=p->EM_ttheta0 ;
    double ttheta2=p->EM_ttheta2 ;
    double hm=p->EM_hm;
    //1. Getting the function
    double residual=aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+ggammam*pow(em,pphi))-(ttheta2*aalpha2m10+aalpha2m*0.92*ttheta0*ttheta2)*ggammam*pow(em,pphi-1);
    
    return(residual);
}



//===================================================================================//
// Solver if the parent is single. Solver will be used in case any of the parents is single
// but I call it Meffortsolver because it is larger the number of women single than fathers
//===================================================================================//
// [[Rcpp::export]]
long double Meffortsolver( double aalpha2m, double aalpha4m, double ggammaf,  double ef, double pphi, double ttheta, double hm) {
    int status;
    int iter=0, max_iter=100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r=0, r_expected=22;
    double x_lo=1.0e-20, x_hi=100000;
    //First we need to make sure that there is not a corner solution at zero:
    double emtry=1.0e-10;
    double costzero=aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+(1-ggammaf)*(pow(emtry,pphi)))-ttheta*aalpha2m*(1-ggammaf)*(pow(emtry,pphi-1));
    //It will only make sense to do the solution if the profits are zero are positive. otherwise
    //corner solution at em=0;
    if (costzero<0){
        gsl_function F;       struct em_my_params params={aalpha2m, aalpha4m, ggammaf,  ef, pphi, ttheta, hm};
        F.function = &F_root_em;
        F.params = &params;
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        
        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        //printf ("using %s method /n", gsl_root_fsolver_name (s));
        //printf ("%5s [%9s, %9s] %9s %10s %9s/n",
        //   "iter" , "lower" , "upper" , "root" , "err" , "err(est)");
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            r = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
            
            //if (status == GSL_SUCCESS)
            // printf ("Converged:\n");
            //printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
            //iter, x_lo, x_hi,r, r - r_expected, x_hi - x_lo);
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fsolver_free (s);
    }
    //If cost is excessive then we set effort of mom at zero:
    else if (costzero>=0){
        r=1.0e-5;
    }
    return(r);
}




//2010

//===================================================================================//
// Solver if the parent is single. Solver will be used in case any of the parents is single
// but I call it Meffortsolver because it is larger the number of women single than fathers
//===================================================================================//
// [[Rcpp::export]]
long double Meffortsolver10( double aalpha2m, double aalpha2m10,
                            double aalpha4m, double ggammaf,  double ef, double pphi, double ttheta0, double ttheta2, double hm) {
    //cout << aalpha2m << "aalpha2m "<< endl;
    //cout << aalpha2m10 << " aalpha2m10 " << endl;
    //cout << aalpha4m << " aalpha4m" << endl;
    //cout << ggammaf << "ggammaf " << endl;
    //cout << ef << " ef " << endl;
    //cout << pphi << " phi"<< endl;
    //cout << ttheta0 << " ttheta0 "<< endl;
    //cout << ttheta2 << " ttheta2"<< endl;
    //cout << hm << "hm "<< endl;
    int status;
    int iter=0, max_iter=100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r=0, r_expected=22;
    double x_lo=1.0e-20, x_hi=100000;
    //First we need to make sure that there is not a corner solution at zero:
    //cout << aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+(1-ggammaf)*(pow(x_lo,pphi)))-(ttheta2*aalpha2m10+0.92*ttheta0*ttheta2*aalpha2m)*(1-ggammaf)*(pow(x_lo,pphi-1)) << " lo" << endl;
    //cout << aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+(1-ggammaf)*(pow(x_hi,pphi)))-(ttheta2*aalpha2m10+0.92*ttheta0*ttheta2*aalpha2m)*(1-ggammaf)*(pow(x_hi,pphi-1)) << " high" << endl;
    double emtry=1.0e-20;
    double costzero=aalpha4m*(1+hm)*(ggammaf*pow(ef,pphi)+(1-ggammaf)*(pow(emtry,pphi)))-(ttheta2*aalpha2m10+0.92*ttheta0*ttheta2*aalpha2m)*(1-ggammaf)*(pow(emtry,pphi-1));
    //It will only make sense to do the solution if the profits are zero are positive. otherwise
    //corner solution at em=0;
    //cout << costzero << " costzero  " << endl;
    if (costzero<0){
        gsl_function F;       struct em_my_params10 params={aalpha2m, aalpha2m10,
            aalpha4m, ggammaf,  ef, pphi, ttheta0, ttheta2, hm};
        F.function = &F_root_em10;
        F.params = &params;
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        
        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        //printf ("using %s method /n", gsl_root_fsolver_name (s));
        //printf ("%5s [%9s, %9s] %9s %10s %9s/n",
        //   "iter" , "lower" , "upper" , "root" , "err" , "err(est)");
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            r = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
            
            //if (status == GSL_SUCCESS)
            // printf ("Converged:\n");
            //printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
            //iter, x_lo, x_hi,r, r - r_expected, x_hi - x_lo);
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fsolver_free (s);
    }
    //If cost is excessive then we set effort of mom at zero:
    else if (costzero>=0){
        r=1.0e-5;
    }
    //cout << r << " r " << endl;
    
    return(r);
}





//===================================================================================//
//0.3 Utility Functions
//===================================================================================//
//---------------------------
// [[Rcpp::export]]
long double F_utility(double aalpha1j, double aalpha2j, double aalpha3j, double aalpha4j,  double cj, double ej, double hj, double S){
    if (S<=0){
        S=1.0e-5;
    }
    
    if (S>=1.0e+10){
        S=5.0e+10;
    }
    
    
    double utility=aalpha1j*log(cj+0.1)+aalpha2j*log(S)-aalpha3j*hj-aalpha4j*(1+hj)*ej;
    if (utility>5.0e+10){
        utility=5.0e+10;
    }
    
    if (utility<-5.0e+10){
        utility=-5.0e+10;
    }
    
    return(utility);
}

//Utility 10
// [[Rcpp::export]]
long double F_utility10(double aalpha1j, double aalpha2j, double aalpha3j, double aalpha4j,
                        double aalpha5j, double cj, double ej, double hj, double S,
                        double a){
    
    
    if (S<=0){
        S=1.0e-5;
    }
    
    if (S>=1.0e+10){
        S=5.0e+10;
    }
    
    double utility=aalpha1j*log(cj+0.1)+aalpha2j*log(S)-aalpha3j*hj-
    aalpha4j*(1+hj)*ej-aalpha5j*hj*(1-a);
    if (utility>5.0e+10){
        utility=5.0e+10;
    }
    if (utility<-5.0e+10){
        utility=-5.0e+10;
    }
    
    return(utility);
}

//===================================================================================//
//0.4 Consumption level
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption(double aalpha1j, double aalpha2j, double ttheta1, double wj, double yj, double hj){
    double num=aalpha1j*(yj+wj*hj);
    //double isf=std::isfinite(num);
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=aalpha2j*ttheta1+aalpha1j;
    if (den<=0){
        den=1.0e-5;
    }
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Father
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption_TOG(double aalpha1f, double aalpha2f, double aalpha2m,
                              double Investment,double ttheta1, double mmu,
                              double priceINV){
    
    double num=aalpha1f*mmu*Investment*priceINV;
    if (num<=0){
        num=1.0e-10;
    }
    if (num>=1.0e+10){
        num=1.0e+10;
    }
    double den=aalpha2f*mmu+aalpha2m*(1-mmu);
    den=den*ttheta1;
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    double consumption=num/den;
    
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Mother
//===================================================================================//
// [[Rcpp::export]]
long double M_consumption_TOG(double aalpha1m, double aalpha2f, double aalpha2m,
                              double Investment,double ttheta1, double mmu,
                              double priceINV){
    
    double num=aalpha1m*(1-mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=aalpha2f*mmu+aalpha2m*(1-mmu);
    
    den=den*ttheta1;
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT father 10
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption_TOG10(double aalpha1f10, double aalpha2f, double aalpha2m,
                                double aalpha2f10, double aalpha2m10,
                                double Investment, double ttheta0,
                                double ttheta1,
                                double mmu, double priceINV){
    
    double num=aalpha1f10*(mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=ttheta1*(aalpha2f10*mmu+aalpha2m10*(1-mmu))+
    0.92*ttheta0*ttheta1*(aalpha2f*mmu+aalpha2m*(1-mmu));
    
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Mother 10
//===================================================================================//
// [[Rcpp::export]]
long double M_consumption_TOG10(double aalpha1m10, double aalpha2f, double aalpha2m,
                                double aalpha2f10, double aalpha2m10,
                                double Investment, double ttheta0,
                                double ttheta1,
                                double mmu, double priceINV){
    
    double num=aalpha1m10*(1-mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=ttheta1*(aalpha2f10*mmu+aalpha2m10*(1-mmu))+
    0.92*ttheta0*ttheta1*(aalpha2f*mmu+aalpha2m*(1-mmu));
    
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}
//===================================================================================//
//0.5 Investment decision if single
//===================================================================================//
// [[Rcpp::export]]
long double F_investment(double aalpha1j, double aalpha2j, double ttheta1, double wj, double yj, double hj,double priceINV){
    double num=aalpha2j*ttheta1*(yj+wj*hj);
    
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=priceINV*(aalpha1j+aalpha2j*ttheta1);
    if (den==0){
        den=1.0e-5;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}



//===================================================================================//
//0.5 Investment decision if single 2010
//===================================================================================//
// [[Rcpp::export]]
long double F_investment10(double aalpha1j10, double aalpha2j,
                           double aalpha2j10,
                           double ttheta1,
                           double ttheta0,
                           double wj, double yj, double hj,double priceINV){
    
    
    double num=(aalpha2j10*ttheta1+0.92*aalpha2j*ttheta0*ttheta1)*(yj+wj*hj);
    
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=priceINV*(aalpha1j10+aalpha2j10*ttheta1+aalpha2j*ttheta0*ttheta1*0.92);
    if (den==0){
        den=1.0e-5;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}


//===================================================================================//
//0.5.1 Investment decision if joint
//===================================================================================//
//Verified
// [[Rcpp::export]]
long double F_invcouple(double aalpha1m, double aalpha1f, double aalpha2m, double aalpha2f,
                        double mmu, double hf,double hm,double wf,
                        double wm, double Yf, double Ym, double ttheta1, double priceINV){
    
    
    
    
    double TI=hf*wf+hm*wm+Yf+Ym;
    double k1=aalpha1f*mmu+aalpha1m*(1-mmu);
    double k2=aalpha2f*mmu+aalpha2m*(1-mmu);
    double denominator=(k1+k2*ttheta1)*priceINV;
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.0e+10){
        denominator=1.0e+10;
    }
    double Investment=TI*k2*ttheta1/denominator;
    if (Investment>5.0e+10){
        Investment=5.0e+10;
    }
    
    if (Investment<=0){
        Investment=1.0e-5;
    }
    return(Investment);
}

//===================================================================================//
//0.5.1 Investment decision if joint 2010
//===================================================================================//
// [[Rcpp::export]]
long double F_invcouple10(double aalpha1f10, double aalpha1m10, double aalpha2f,
                          double aalpha2m,
                          double aalpha2f10,
                          double aalpha2m10,
                          double mmu,
                          double hf,
                          double hm,
                          double wf,
                          double wm,
                          double Yf, double Ym,
                          double ttheta0,
                          double ttheta1, double priceINV){
    
    double TI=hf*wf+hm*wm+Yf+Ym;
    double k1=aalpha1f10*mmu+aalpha1m10*(1-mmu);
    double k2=aalpha2f*mmu+aalpha2m*(1-mmu);
    double denominator=(k1+k2*ttheta1+
                        0.92*ttheta0*ttheta1*(aalpha2f10*mmu+(1-mmu)*aalpha2m10))*priceINV;
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.0e+10){
        denominator=1.0e+10;
    }
    double Investment=TI*(k2*ttheta1+k2*ttheta0*ttheta1*0.92)/denominator;
    if (Investment>5.0e+10){
        Investment=5.0e+10;
    }
    
    if (Investment<=1.0e-5){
        Investment=1.0e-5;
    }
    
    return(Investment);
}



//================================
//0.6 Bargaining power
//================================
// [[Rcpp::export]]
long double F_mmu(double llambda0, double llambda1, double llambda2, double llambda3,
                  double llambda4, double llambda5, double llambda6, double llambda7,
                  double llambda8,
                  double wf, double wm,
                  double Yf, double Ym, double kkappa, double mmuLB,
                  double mmuUB, double Fage12, double Mage12,
                  double Fyrschool12, double Myrschool12,
                  double FMRATIO, double Unemployment, double Wageratio, double Distance){
    //0 Defining ratio of nonlabor income
    
    double ymratio=0;
    if (Yf+Ym==0){
        ymratio=0.5;
    }
    else {
        ymratio=Yf/(Yf+Ym);
    }
    if (wm==0){
        wm=1.0e-5;
    }
    double wmratio=wf/wm;
    //1 Defining what goes inside of the exponent
    double inexp=llambda0+llambda1*(wmratio)+llambda2*ymratio+
    llambda3*(Fage12-Mage12)+llambda4*(Fyrschool12-Myrschool12)+
    llambda5*FMRATIO+llambda6*Unemployment+llambda7*Wageratio+llambda8*Distance
    +kkappa;
    if (inexp>100){
        inexp=100;
    }
    double mmu=mmuLB+(mmuUB-mmuLB)*((exp(inexp))/(1+exp(inexp)));
    return(mmu);
}




//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================







//--------------------------------------------
//1. Block of Likelihood Functions definition
//--------------------------------------------
//=================================================================================
//1.0. Likelihood of Wages
//Checked conditional on F_predwage. Works fine
//=================================================================================
// [[Rcpp::export]]
long double F_likelihood_wage(double bbeta0, double bbeta1, double bbeta2, double bbeta3, double stdwage, double Schooling, double Age, double Wage){
    
    //1. Obtaining the predicted Wage:
    double predwage=F_predwage(bbeta0,bbeta1,bbeta2,bbeta3,Schooling,Age);
    //2. Once we obtain the predicted wage we proceed to
    //computing the likelihood Function
    //2.0. The first step is to define the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    
    //2.1 Now I will compute the actual likelihood of the observed Wage
    double likelihood_wage=log(Wage)-log(predwage);
    likelihood_wage=likelihood_wage/stdwage;
    likelihood_wage=pdf(normdist,likelihood_wage);
    
    likelihood_wage=likelihood_wage/stdwage;
    //2.2 If the pdf given is equal to zero we will approximate it to very close to zero
    if (likelihood_wage==0){
        likelihood_wage=1.0e-100;
    }
    likelihood_wage=log(likelihood_wage);
    return(likelihood_wage);
}



//=================================================================================
//1.1 Likelihood of skills --Checked conditional on F_predskills
//=================================================================================
// [[Rcpp::export]]
long double F_likelihood_skills(double ddelta0, double ddelta1,
                                double ddelta2, double ddelta3,
                                double ddelta4,
                                double ttheta0, double ttheta1,
                                double ttheta2, double pphi,
                                double gf, double gm,
                                double Fe, double Me,
                                double I, double S0,
                                double obS, double stdskills,
                                double Age, double CHILDCARE,
                                double PG,
                                double Hmembers){
    
    
    //1. Getting predicted level of skills
    double predskills=F_predskills(ddelta0, ddelta1, ddelta2, ddelta3, ddelta4,Age, ttheta0, ttheta1, ttheta2, pphi, gf, gm, Fe, Me, I, S0, CHILDCARE, PG,Hmembers);
    
    
    //1. Computing the likelihood of Skills
    boost::math::normal_distribution<> normdist (0,1);
    
    if (predskills==0){
        predskills=1.0e-5;
    }
    double loglik=log(obS)-log(predskills);
    
    loglik=loglik/stdskills;
    loglik=pdf(normdist,loglik);
    loglik=loglik/stdskills;
    if(loglik==0){
        loglik=-28;
    }
    
    else{
        loglik=log(loglik);
    }
    
    return(loglik);
}


//=================================================================================
//1.2 Likelihood of effort observed
//=================================================================================

//------------------------------------------
//1.2.0 Effort observed by father-not single -- Checked, conditional on F_effort_f
//------------------------------------------
// [[Rcpp::export]]
long double F_likelihood_feffort(double Feffort_observed, double mmu, double ggammaf,  double aalpha2f, double aalpha2m, double aalpha4m,  double aalpha4f,   double hf, double hm, double ttheta0, double ttheta1, double ttheta2,double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of father
    //-----------------------------------------
    
    double effort_f=F_effort_f(mmu,ggammaf,aalpha2f,aalpha2m,aalpha4f,aalpha4m,hf,hm,pphi,ttheta2);
    
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_feffort=Feffort_observed-effort_f;
    
    loglik_feffort=loglik_feffort/stdeffort;
    loglik_feffort=pdf(normdist,loglik_feffort);
    loglik_feffort=loglik_feffort/stdeffort;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_feffort==0){
        loglik_feffort=1.0e-100;
    }
    loglik_feffort=log(loglik_feffort);
    return(loglik_feffort);
}

//--------------------------------------------------------
// 1.2.1 Likelihood of effort observed by mother-not single -- checked conditional on F_effort_m
// -------------------------------------------------------
// [[Rcpp::export]]
long double F_likelihood_meffort(double Meffort_observed, double mmu, double ggammaf,  double aalpha2f, double aalpha2m,  double aalpha4m,   double aalpha4f,  double hf, double hm, double ttheta0, double ttheta1, double ttheta2, double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of mother
    //-----------------------------------------
    double effort_m=F_effort_m(mmu,ggammaf,aalpha2f,aalpha2m,aalpha4f,aalpha4m,hf,hm,pphi,ttheta2);
    
    
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_meffort=Meffort_observed-effort_m;
    
    loglik_meffort=loglik_meffort/stdeffort;
    
    loglik_meffort=pdf(normdist,loglik_meffort);
    
    loglik_meffort=loglik_meffort/stdeffort;
    
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_meffort==0){
        loglik_meffort=1.0e-100;
    }
    
    loglik_meffort=log(loglik_meffort);
    return(loglik_meffort);
}



//------------------------------------------------------
//Effort observed by mother-single -- checked conditional on Meffortsolver 2012
//-----------------------------------------------------

// ----------------------------------
// [[Rcpp::export]]
long double F_likelihood_meffort_single(double Meffort_observed, double ef, double ggammaf,double ggammam,   double aalpha2m, double aalpha4m,  double hm,  double ttheta2, double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of mother
    //-----------------------------------------
    double effort_m=Meffortsolver(  aalpha2m,  aalpha4m, ggammaf,  ef, pphi, ttheta2, hm);
    
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_meffort=Meffort_observed-effort_m;
    loglik_meffort=loglik_meffort/stdeffort;
    loglik_meffort=pdf(normdist,loglik_meffort);
    loglik_meffort=loglik_meffort/stdeffort;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_meffort==0){
        loglik_meffort=1.0e-100;
    }
    loglik_meffort=log(loglik_meffort);
    return(loglik_meffort);
}





//------------------------------------------------------
//Effort observed by father-single: checked conditional on Meffortsolver
//-----------------------------------------------------

// ----------------------------------
// [[Rcpp::export]]
long double F_likelihood_feffort_single(double Feffort_observed, double em, double ggammaf, double ggammam, double  aalpha2f,  double aalpha4f,  double hf, double ttheta2, double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of mother
    //-----------------------------------------
    double effort_m=Meffortsolver(  aalpha2f,  aalpha4f, ggammam,  em, pphi, ttheta2, hf);
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_meffort=Feffort_observed-effort_m;
    loglik_meffort=loglik_meffort/stdeffort;
    loglik_meffort=pdf(normdist,loglik_meffort);
    loglik_meffort=loglik_meffort/stdeffort;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_meffort==0){
        loglik_meffort=1.0e-100;
    }
    loglik_meffort=log(loglik_meffort);
    return(loglik_meffort);
}



//------------------------------------------------------------
// 1.3 Likelihood of utility observed for discrete decisions
// we will do the four combinations: mother and father if married
// or single--checked conditional on intermediate functions
// ---------------------------------------------------------
// [[Rcpp::export]]
long double F_like_util_momsingle(double aalpha1m,  double aalpha3m, double aalpha4m,  double hm, double wm, double ym, double hchores, double Investment, double em, double S ,  double stdh1 ){
    //1. Defining the non-random components of the utilities
    //1.0 Consumption levels
    double cons_h0=ym+0.1-Investment;
    //1.0.0 If consumption level given by ym-I is smaller or equal to zero we need to fix it:
    if (cons_h0<=0){
        cons_h0=1.0e-5;
    }
    //1.0.1 Consumption level if works defined in reference to non-work Consumption
    double cons_h1=cons_h0+wm;
    
    
    
    
    //2. Defining the non-random components of the utilities:
    
    
    double UNR0=aalpha1m*log(cons_h0)-aalpha4m*em;
    double UNR1=aalpha1m*log(cons_h1)-aalpha3m-aalpha4m*em*2;
    //3.0. And now computing the likelihood. Defining numerator
    double nume=UNR1-UNR0;
    //3.1 denominator
    //double deno=pow(pow(stdh1,2)+1,0.5); //Remember one is normalized to one
    //before I put normalized one sd to be one. Now the problem is that it cannot be as small as we want
    //the whole variance which might be an issue. In this point I will normalize the sum of both
    double deno=stdh1;
    if (deno==0){
        deno=1.0e-10;
    }
    //1.3 Defining what goes inside of the cdf
    double incdf=nume/deno;
    //1.4 Computing the corresponding cdf
    boost::math::normal_distribution<> normdist (0 ,1);
    double cdfT=cdf(normdist,incdf);
    //If one or the other we need to take that into account
    if (hm==0){
        cdfT=1-cdfT;
    }
    if (cdfT==0){
        cdfT=1.0e-10;
    }
    cdfT=log(cdfT);
    return(cdfT);
}

//1.4 Likelihood of Investment levels without using predicted levels function
// [[Rcpp::export]]
long double F_likelihood_investment_single(double OBsinvestment,double priceINV, double aalpha1, double aalpha2, double ttheta0, double  Cfactorinv12, double  Mwage12, double  Mnly12, double  Mfraclabor12, double stdinvestment){
    double Totalincome=Mwage12*Mfraclabor12+Mnly12;
    double den=aalpha2*ttheta0+aalpha1;
    //0. If denominator close to zero re-arrange it:
    if (den==0){
        den=1.0e-5;
    }
    double num=aalpha2*ttheta0;
    double Investment=Totalincome*(num/den);
    Investment=Investment/priceINV;
    //1.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //1.1 And computing the likelihood
    double loglik_investment=OBsinvestment-Investment;
    loglik_investment=loglik_investment/stdinvestment;
    loglik_investment=pdf(normdist,loglik_investment);
    loglik_investment=loglik_investment/stdinvestment;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_investment==0){
        loglik_investment=1.0e-100;
    }
    loglik_investment=log(loglik_investment);
    return(loglik_investment);
    
}
//1.4.1 Likelihood of Investment levels  using predicted levels function
// [[Rcpp::export]]
long double F_likelihood_investment(double OBsinvestment,
                                    double Predinvestment,double MEASINV){
    //Getting a normal distribution
    boost::math::normal_distribution<> normdist (0,1);
    
    //Fixing not to take logs of zero
    if (OBsinvestment<=1.0e-5){
        OBsinvestment=1.0e-5;
    }
    if (Predinvestment<=1.0e-5){
        Predinvestment=1.0e-5;
    }
    
    //If variance is too small we need to fix it
    if (MEASINV<=1.0e-5){
        MEASINV=1.0e-5;
    }
    
    double loglikeinv=log(OBsinvestment)-log(Predinvestment);
    loglikeinv=loglikeinv/MEASINV;
    loglikeinv=pdf(normdist,loglikeinv);
    loglikeinv=loglikeinv/MEASINV;
    if (loglikeinv<=1.0e-5){
        loglikeinv=1.0e-5;
    }
    loglikeinv=log(loglikeinv);
    return(loglikeinv);
}



//1.5 Likelihood of mmu

// [[Rcpp::export]]
long double F_likelihood_mmu(double llambda0, double llambda1, double llambda2,
                             double llambda3, double llambda4,
                             double llambda5, double llambda6, double llambda7, double llambda8,
                             double Fwageinlik, double Mwageinlik,
                             double Mnly12, double Fnly12,
                             double kkappa, double mmuLB, double mmuUB,
                             double Fage12, double Mage12,
                             double Fyrschool12, double Myrschool12,
                             double FMRATIO, double Unemployment, double Wageratio, double Distance,
                             double Hbarg,
                             double MEASMMu){
    
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    
    //First I will find the predicted level of mmu
    double loglikmmu=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                           llambda5,llambda6,llambda7,llambda8,
                           Fwageinlik,Mwageinlik,
                           Mnly12,Fnly12,kkappa,mmuLB,mmuUB,Fage12,Mage12,
                           Fyrschool12,Myrschool12,FMRATIO,Unemployment,Wageratio,Distance);
    
    if(loglikmmu<=0){
        loglikmmu=1.0e-5;
    }
    loglikmmu=log(Hbarg)-log(loglikmmu);
    loglikmmu=loglikmmu/MEASMMu;
    loglikmmu=pdf(normdist,loglikmmu);
    loglikmmu=loglikmmu/MEASMMu;
    if (loglikmmu==0){
        loglikmmu=1.0e-10;
    }
    loglikmmu=log(loglikmmu);
    return(loglikmmu);
}


//1.6 I will define a generic loglikelihood for an observed variable
//And an error term that follows a log-normal distribution


// [[Rcpp::export]]
long double F_loglikelihood_generic(double observed, double predicted, double variance){
    
    //Fixing all the parameters
    if (observed<=1.0e-5){
        observed=1.0e-5;
    }
    if (predicted<=1.0e-5){
        predicted=1.0e-5;
    }
    if (variance<=1.0e-5){
        variance=1.0e-5;
    }
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    //Getting the estimate
    double likelihood=log(observed)-log(predicted);
    likelihood=likelihood/variance;
    likelihood=pdf(normdist,likelihood);
    likelihood=likelihood/variance;
    if (likelihood<=1.0e-5){
        likelihood=1.0e-5;
    }
    return(likelihood);
}



//This one is not for log-normal but simply normal
// [[Rcpp::export]]
long double F_likgeneric(double observed, double predicted, double variance){
    
    
    if (variance<=1.0e-5){
        variance=1.0e-5;
    }
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    //Getting the estimate
    double likelihood=(observed)-(predicted);
    likelihood=likelihood/variance;
    likelihood=pdf(normdist,likelihood);
    likelihood=likelihood/variance;
    if (likelihood<=1.0e-5){
        likelihood=1.0e-5;
    }
    return(likelihood);
}


//----------------------------------------------------
//2. Block of defining likelihood for all the dataset
//---------------------------------------------------
// [[Rcpp::export]]
long double F_likelihood(vector<double> par_likelihood,
                         vector<double>  Cliveswithfather12,
                         vector<double> Cliveswithmother12,
                         vector<double> Hhchores12,
                         vector<double> Meffort12,
                         vector<double> Feffort12,
                         vector<double> Mage12,
                         vector<double> Fage12,
                         vector<double> Cedad_meses12,
                         vector<double> Ctestsfactorsss2012,
                         vector<double> Ctestsfactor2ss_10,
                         vector<double> Myrschool12,
                         vector<double> Fyrschool12,
                         vector<double> Ffraclabor12,
                         vector<double> Mfraclabor12,
                         vector<double> Mwage12,
                         vector<double> Fwage12,
                         vector<double> Mnly12,
                         vector<double> Fnly12,
                         vector<double> Cfactorinv12,
                         vector<double> Cchildcare12,
                         vector<double> Ccareskills12,
                         vector<double> Hbarg,
                         vector<double> Feffort10,
                         vector<double> Meffort10,
                         vector<double> Cfactorinv10,
                         vector<double> Cbirthfactor,
                         vector<double> Mwage10,
                         vector<double> Fwage10,
                         vector<double> Mnly10,
                         vector<double> Fnly10,
                         vector<double> Cedad_meses10,
                         vector<double> Cliveswithfather10,
                         vector<double> Cliveswithmother10,
                         vector<double> Mage10,
                         vector<double> Fage10,
                         vector<double> Mfraclabor10,
                         vector<double> Ffraclabor10,
                         vector<double> Cchildcare10,
                         vector<double> Hchildcareobs,
                         vector<double> PG,
                         vector<double> MTJH,
                         vector<double> MRATIO,
                         vector<double> Unemployment,
                         vector<double> Wageratio,
                         vector<double>  Distance,
                         vector<double> Magegroup10,
                         vector<double>  Magegroup12,
                         vector<double>  Hmemberstotal10,
                         vector<double>  Hmemberstotal12,
                         int SIZE) {
    
    //------------------------------------------
    //0. Performing the right transformations to
    //guarantee that parameters lie within the
    //specified set
    //-----------------------------------------
    
    //0.1 Aalphas
    //0.1 Aalphas
    double aalpha1m=FT_approach(FT_exp(par_likelihood[0])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha2m=FT_approach(FT_exp(par_likelihood[1])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha3m=FT_approach(FT_exp(par_likelihood[2])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha40m=FT_approach(FT_exp(par_likelihood[3])/
                                 (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    
    
    
    double aalpha41m=FT_approach(FT_exp(par_likelihood[4])/
                                 (1+FT_exp(par_likelihood[4])))*aalpha40m;
    
    double aalpha4m=0;
    
    
    
    
    
    
    //0.2 Skills
    double ttheta0=FT_exp(par_likelihood[5])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double ttheta1=FT_exp(par_likelihood[6])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double ttheta2=FT_exp(par_likelihood[7])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double pphi=FT_pphi(par_likelihood[8]);
    double stdskills=FT_exp(par_likelihood[9]);
    
    
    double ggammaf=FT_gf(FT_exp(par_likelihood[10])/
                         (FT_exp(par_likelihood[10])+FT_exp(par_likelihood[11])));
    
    
    double ggammam=FT_gf(FT_exp(par_likelihood[11])/
                         (FT_exp(par_likelihood[10])+FT_exp(par_likelihood[11])));
    
    
    //0.3 Wages for mother
    double bbeta0m=par_likelihood[12];
    double bbeta1m=par_likelihood[13];
    double bbeta2m=par_likelihood[14];
    double bbeta3m=par_likelihood[15];
    double stdwm=FT_exp(par_likelihood[16]);
    
    //0.4 Measurement system
    double stdeffortmom=FT_exp(par_likelihood[17]);
    double stdinvestment=FT_exp(par_likelihood[18]);
    
    //0.5 Price of investment
    double priceINV=FT_exp(par_likelihood[19]);
    
    //0.6 Childcare situation
    double ddelta0=par_likelihood[20]; //Intercept
    double ddelta1=par_likelihood[21]; //Age shifter
    double ddelta2=par_likelihood[22]; //Childcare production
    double ddelta3=par_likelihood[23]; //PG in 2010
    double ddelta3_12=par_likelihood[24]; //Childcare production PG in 2012
    
    
    
    
    //Preferences for the father
    //0.1 Aalphas
    double aalpha1f=FT_approach(FT_exp(par_likelihood[25])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha2f=FT_approach(FT_exp(par_likelihood[26])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha3f=FT_approach(FT_exp(par_likelihood[27])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha40f=FT_approach(FT_exp(par_likelihood[28])/
                                 (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                  FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha41f=FT_approach(FT_exp(par_likelihood[29])/
                                 (1+FT_exp(par_likelihood[29])))*aalpha40f;
    double aalpha4f=0;
    
    
    
    //0.3 Wages
    double bbeta0f=par_likelihood[30];
    double bbeta1f=par_likelihood[31];
    double bbeta2f=par_likelihood[32];
    double bbeta3f=par_likelihood[33];
    double stdwf=FT_exp(par_likelihood[34]);
    
    //0.4 Measurement system
    double stdeffortfat=FT_exp(par_likelihood[35]);
    
    
    //0.6 Bargaining power parameters
    double llambda0=par_likelihood[36]; //Intercept
    double llambda1=par_likelihood[37]; //Wage ratio
    double llambda2=par_likelihood[38]; //Non-labor income ratio
    double llambda3=par_likelihood[39]; //Age difference
    double llambda4=par_likelihood[40]; //Education difference
    double stdmmu=FT_exp(par_likelihood[41]); //STD of barg shock
    double mmuLB=0.2; //Lower bound for bargaining power
    double mmuUB=0.8; //Upper bound for bargaining power
    double mupred=0; //Predicted level of bargaining power
    double mupred10=0; //Predicted level of bargaining power 2010
    //0.7 Measurement system parameters
    
    double MEASSkills=FT_exp(par_likelihood[42]);
    double MEASFeffort=FT_exp(par_likelihood[43]);
    double MEASMeffort=FT_exp(par_likelihood[44]);
    
    double MEASMMu=FT_exp(par_likelihood[45]);
    double MEASINV=FT_exp(par_likelihood[46]);
    
    
    //0.8 Preferences in the first period
    //Mother
    double aalpha1m10=FT_approach(FT_exp(par_likelihood[47])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha2m10=FT_approach(FT_exp(par_likelihood[48])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha3m10=FT_approach(FT_exp(par_likelihood[49])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha40m10=FT_approach(FT_exp(par_likelihood[50])/
                                   (FT_exp(par_likelihood[47])+
                                    FT_exp(par_likelihood[48])+
                                    FT_exp(par_likelihood[49])+
                                    FT_exp(par_likelihood[50])+
                                    FT_exp(par_likelihood[51])));
    
    double aalpha5m10=FT_approach(FT_exp(par_likelihood[51])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha41m10=FT_approach(FT_exp(par_likelihood[52])/
                                   (1+FT_exp(par_likelihood[52])))*aalpha40m10;
    double aalpha4m10=0;
    //Father
    double aalpha1f10=FT_approach(FT_exp(par_likelihood[53])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    double aalpha2f10=FT_approach(FT_exp(par_likelihood[54])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    double aalpha3f10=FT_approach(FT_exp(par_likelihood[55])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    
    double aalpha40f10=FT_approach(FT_exp(par_likelihood[56])/
                                   (FT_exp(par_likelihood[53])+
                                    FT_exp(par_likelihood[54])+
                                    FT_exp(par_likelihood[55])+
                                    FT_exp(par_likelihood[56])+
                                    FT_exp(par_likelihood[57])));
    
    double aalpha5f10=FT_approach(FT_exp(par_likelihood[57])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    
    double aalpha41f10=FT_approach(FT_exp(par_likelihood[58])/
                                   (1+FT_exp(par_likelihood[58])))*aalpha40m10;
    double aalpha4f10=0;
    
    //Price of childcare
    double pchildcare0=FT_exp(par_likelihood[59]);
    double pchildcare1=FT_exp(par_likelihood[60]); //Regarding the Distance
    
    //Preference shocks
    double MshockWA=FT_exp(par_likelihood[61]);
    double MshockNWA=FT_exp(par_likelihood[62]);
    double MshockWNA=FT_exp(par_likelihood[63]);
    double MshockNWNA=FT_exp(par_likelihood[64]);
    double FshockWA=FT_exp(par_likelihood[65]);
    double FshockNWA=FT_exp(par_likelihood[66]);
    double FshockWNA=FT_exp(par_likelihood[67]);
    double FshockNWNA=FT_exp(par_likelihood[68]);
    
    
    //Final elements to take into account:
    //Distribution factors
    double llambda5=par_likelihood[69];//Female to male ratio
    double llambda6=par_likelihood[70];//Unemployment ratio
    double llambda7=par_likelihood[71];//Wage ratio
    double llambda8=par_likelihood[72];//Distance a centros de protección
    
    //Preference heterogeneity
    
    //Preference heterogeneity terms will have the same restrictions as aalpha41.
    //Between zero and one.
    
    //Mujer trabajadora y jefa de hogar 2010
    double aalpha3m10_mtjh=FT_approach(FT_exp(par_likelihood[73])/
                                       (1+FT_exp(par_likelihood[73])))*aalpha3m10;
    
    //Mujer trabajadora y jefa de hogar 2012
    double aalpha3m12_mtjh=FT_approach(FT_exp(par_likelihood[74])/
                                       (1+FT_exp(par_likelihood[74])))*aalpha3m;
    
    //Age preferences for leisure
    double aalpha3mage10=FT_approach(FT_exp(par_likelihood[75])/
                                     (1+FT_exp(par_likelihood[75])))*aalpha3m10;
    
    double aalpha3mage12=FT_approach(FT_exp(par_likelihood[76])/
                                     (1+FT_exp(par_likelihood[76])))*aalpha3m;
    
    
    //People with children
    double ddelta4=par_likelihood[77]; //How having one additional person affects exp(R).
    
    
    
    //Measurement systems for the bootstrap filter
    double STDS0=FT_exp(0); //Leave it out of the estimation parameters from now
    double STDPG=FT_exp(0);
    
    //==============================
    
    //1. Initializing Likelihood
    double loglik=0;
    double logcontribobs_ii=0; //Contribution of observation ii
    double logcontribobs_rr=0; //Contribution of simulation rr
    double logcontribobs_rr10=0;
    double logcontribdecision=0; // Sum of contributions of decisions correct
    double logcontribdecision10=0;
    double logcontribdecision_rr=0;// Contribution of decision in simulation rr
    double logcontribdecision_rr10=0;
    double likwageM=0;
    double likwageM10=0;
    double likwageF=0;
    double likwageF10=0;
    double likbehavior=0;
    double likbehavior10=0;
    double loglikskills=0;
    double loglikskills10=0;
    double loglikskills_rr=0;
    double loglikskills_rr10=0;
    double loglikeffm=0;
    double loglikeffm10=0;
    double loglikeffm_rr=0;
    double loglikeffm_rr10=0;
    double loglikefff=0;
    double loglikefff_rr=0;
    double loglikefff10=0;
    double loglikefff_rr10=0;
    double loglikmmu=0;
    double loglikmmu_rr=0;
    double loglikeinv=0;
    double loglikeinv_rr=0;
    double loglikeinv10=0;
    double loglikeinv_rr10=0;
    
    //2. All the initializations that take place in for and while loop are done here
    //NW->No work; W-> Work
    int rr=0; //Initializing the number of times for the integration
    //2.1 We need to initialize all factor variables in the possible alternatives.
    //work-no work, childcare-no childcare, etc.
    // Effort of mother
    double Meffort_factor12shock=0; //Shock of effort for mother
    double MeffortMean12NW=0; //Mean of effort of mother initializer if single
    double MeffortMean12NWNW_TOG=0; //Mean of effort of mother initializer if nw with spouse
    double MeffortMean12WNW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12NWW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12WW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12W=0;  //Mean of effort if mother works
    double Meffort_factor12NW=0; // Effort level of mother initializer
    double Meffort_factor12W=0; // Effort level of mother initializer
    double Meffort_factor10shock=0; //Shock of effort for mother
    double MeffortMean10NW=0; //Mean of effort of mother initializer if single
    double MeffortMean10NWNW_TOG=0; //Mean of effort of mother initializer if nw with spouse
    double MeffortMean10WNW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10NWW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10WW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10W=0;  //Mean of effort if mother works
    double Meffort_factor10NW=0; // Effort level of mother initializer
    double Meffort_factor10W=0; // Effort level of mother initializer
    // Effort of father
    double Feffort_factor12shock=0; //Shock of effort for father
    
    double FeffortMean12NWNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12WNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12NWW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12WW_TOG=0; // mean of effort if father nw and with spouse
    
    double Feffort_factor10shock=0; //Shock of effort for father
    
    double FeffortMean10NWNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10WNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10NWW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10WW_TOG=0; // mean of effort if father nw and with spouse
    
    // Skills
    double CSkills_factor12NW=0; //Skills level initializer
    double CSkills_factor12W=0; //Skills level initializer
    double CSkills_factor12NWNW_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor12NWW_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor12WNW_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor12WW_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor12shock=0;  //Skills shock
    double CSkillsMean12NW=0; // Mean of skills level initializer
    
    double CSkillsMean12W=0; // Mean of skills level initializer
    
    double CSkills_factor10NWA=0; //Skills level initializer
    double CSkills_factor10NWNA=0; //Skills level initializer
    double CSkills_factor10WA=0; //Skills level initializer
    double CSkills_factor10WNA=0; //Skills level initializer
    double CSkills_factor10NWNWA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10NWNWNA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10NWWA_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor10NWWNA_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor10WNWA_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor10WNWNA_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor10WWA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10WWNA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10shock=0;  //Skills shock
    double CSkillsMean10NWA=0; // Mean of skills level initializer
    double CSkillsMean10NWNA=0; // Mean of skills level initializer
    double CSkillsMean10WA=0; // Mean of skills level initializer
    double CSkillsMean10WNA=0; // Mean of skills level initializer
    
    //Investment
    double CInvestmentMean12NW=0;  // Investment level initializer
    double CInvestmentMean12W=0;  // Investment level initializer
    double CInvestmentMean12NWNW_TOG=0; //Investment nwnw tog
    double CInvestmentMean12WNW_TOG=0; //Investment nw tog
    double CInvestmentMean12NWW_TOG=0; //Investment nww tog
    double CInvestmentMean12WW_TOG=0; //Investment ww tog
    double CInv_factor12shock=0; // Shock of investment
    double CInv_factor12NW=0;//Mean of investment level
    double CInv_factor12W=0;//Mean of investment level
    
    double CInvestmentMean10NWA=0;  // Investment level initializer
    double CInvestmentMean10NWNA=0;  // Investment level initializer
    double CInvestmentMean10WA=0;  // Investment level initializer
    double CInvestmentMean10WNA=0;  // Investment level initializer
    double CInvestmentMean10NWNWA_TOG=0; //Investment nwnw tog
    double CInvestmentMean10NWNWNA_TOG=0; //Investment nwnw tog
    double CInvestmentMean10WNWA_TOG=0; //Investment nw tog
    double CInvestmentMean10WNWNA_TOG=0; //Investment nw tog
    double CInvestmentMean10NWWA_TOG=0; //Investment nww tog
    double CInvestmentMean10NWWNA_TOG=0; //Investment nww tog
    double CInvestmentMean10WWA_TOG=0; //Investment ww tog
    double CInvestmentMean10WWNA_TOG=0; //Investment ww tog
    double CInv_factor10shock=0; // Shock of investment
    double CInv_factor10NWA=0;//Mean of investment level
    double CInv_factor10NWNA=0;//Mean of investment level
    double CInv_factor10WA=0;//Mean of investment level
    double CInv_factor10WNA=0;//Mean of investment level
    
    //Consumption
    //Mother
    double MconsumptionNW=0; //nw single
    double MconsumptionW; //w single
    double MconsumptionNWNW_TOG=0; //nwnw tog
    double MconsumptionWNW_TOG=0; //wnw tog
    double MconsumptionNWW_TOG=0; //nww tog
    double MconsumptionWW_TOG=0; //ww tog
    
    
    double MconsumptionNWA10=0; //nw single
    double MconsumptionNWNA10=0; //nw single
    double MconsumptionWA10=0; //w single
    double MconsumptionWNA10=0; //w single
    double MconsumptionNWNWA_TOG10=0; //nwnw tog
    double MconsumptionNWNWNA_TOG10=0; //nwnw tog
    double MconsumptionWNWNA_TOG10=0; //wnw tog
    double MconsumptionWNWA_TOG10=0; //wnw tog
    double MconsumptionNWWNA_TOG10=0; //nww tog
    double MconsumptionNWWA_TOG10=0; //nww tog
    double MconsumptionWWA_TOG10=0; //ww tog
    double MconsumptionWWNA_TOG10=0; //ww tog
    
    //Father
    
    double FconsumptionNWNW_TOG=0; //nwnw tog
    
    double FconsumptionWNW_TOG=0; //wnw tog
    double FconsumptionNWW_TOG=0; //nww tog
    double FconsumptionNWWNA_TOG=0; //nww tog
    double FconsumptionWW_TOG=0; //ww tog
    
    
    
    double FconsumptionNWNWA_TOG10=0; //nwnw tog
    double FconsumptionNWNWNA_TOG10=0; //nwnw tog
    double FconsumptionWNWA_TOG10=0; //wnw tog
    double FconsumptionWNWNA_TOG10=0; //wnw tog
    double FconsumptionNWWA_TOG10=0; //nww tog
    double FconsumptionNWWNA_TOG10=0; //nww tog
    double FconsumptionWWA_TOG10=0; //ww tog
    double FconsumptionWWNA_TOG10=0; //ww tog
    
    
    //Preference shocks initializer
    //Mother
    double MprefshockW=0; //w single
    double MprefshockNW=0; //nw single
    double MprefshockWA10=0; //w single
    double MprefshockWNA10=0; //w single
    double MprefshockNWA10=0; //nw single
    double MprefshockNWNA10=0; //nw single
    
    //Father
    double FprefshockW=0; //w single
    double FprefshockNW=0; //nw single
    double FprefshockWA10=0; //w single
    double FprefshockWNA10=0; //w single
    double FprefshockNWA10=0; //nw single
    double FprefshockNWNA10=0; //nw single
    
    
    //Utility levels
    //Mother
    double MutilityNW=0; //Utility for mother if she works sing
    double MutilityW=0; //Utility for mother if she doesn't work single
    double MutilityNWNW_TOG=0; //Utility nwnw tog
    double MutilityWNW_TOG=0; //Utility wnw tog
    double MutilityNWW_TOG=0; //Utility nww tog
    double MutilityWW_TOG=0; //Utility ww tog
    
    double MutilityNWA10=0; //Utility for mother if she works sing
    double MutilityNWNA10=0; //Utility for mother if she works sing
    double MutilityWA10=0; //Utility for mother if she doesn't work single
    double MutilityWNA10=0; //Utility for mother if she doesn't work single
    double MutilityNWNWA_TOG10=0; //Utility nwnw tog
    double MutilityNWNWNA_TOG10=0; //Utility nwnw tog
    double MutilityWNWA_TOG10=0; //Utility wnw tog
    double MutilityWNWNA_TOG10=0; //Utility wnw tog
    double MutilityNWWA_TOG10=0; //Utility nww tog
    double MutilityNWWNA_TOG10=0; //Utility nww tog
    double MutilityWWA_TOG10=0; //Utility ww tog
    double MutilityWWNA_TOG10=0; //Utility ww tog
    //Father
    
    double FutilityNWNW_TOG=0; //Utility nwnw tog
    double FutilityWNW_TOG=0; //Utility wnw tog
    double FutilityNWW_TOG=0; //Utility nww tog
    double FutilityWW_TOG=0; //Utility ww tog
    
    
    double FutilityNWNWA_TOG10=0; //Utility nwnw tog
    double FutilityNWNWNA_TOG10=0; //Utility nwnw tog
    double FutilityWNWA_TOG10=0; //Utility wnw tog
    double FutilityWNWNA_TOG10=0; //Utility wnw tog
    double FutilityNWWA_TOG10=0; //Utility nww tog
    double FutilityNWWNA_TOG10=0; //Utility nww tog
    double FutilityWWA_TOG10=0; //Utility ww tog
    double FutilityWWNA_TOG10=0; //Utility ww tog
    
    //Welfare of both
    double WelfareNWNW=0;
    double WelfareWNW=0;
    double WelfareNWW=0;
    double WelfareWW=0;
    
    double WelfareNWNWA10=0;
    double WelfareNWNWNA10=0;
    double WelfareWNWA10=0;
    double WelfareWNWNA10=0;
    double WelfareNWWA10=0;
    double WelfareNWWNA10=0;
    double WelfareWWA10=0;
    double WelfareWWNA10=0;
    
    //Indicator for the decision taken
    double Decision=0; //Predicted decision
    double obDecision=0; //Observed decision
    
    double Decision10=0; //Predicted decision
    double obDecision10=0; //Observed decision
    
    //Wage initializer
    double Mwageinlik=0; //Wage of mother initializer
    double Mwageinlik10=0; //Wage of mother in 2010 initializer
    double Fwageinlik=0; //Wage of mother initializer
    double Fwageinlik10=0;//Wage of father in 2010 initializer
    
    //Bargaining power
    double MMushock=0; //Shock for unobserved heterogeneity of Pareto weight
    double MMushock10=0; //Shock for unobserved heterogeneity of Pareto weight
    
    //Price of childcare initializer
    double pricechildcare=0;
    
    //Set number of simulations to perform
    double RR=1000;
    
    //Particles definitions
    //0 Skills at birth
    double S0rr=0; //Skills at birth
    vector<double> S0Vector; //Vector storking skills at birth
    S0Vector.resize(RR);
    double loglik_S0rr=0; //Likelihood contribution of skills rr
    double loglik_S0=0; //Total likelihood contribution of sum of rr
    
    //1. Skills of primarycaregiver
    double PGrr=0; //Skills of PG
    vector<double> PGVector; //Vector storking skills at birth
    PGVector.resize(RR);
    double loglik_PGrr=0; //Likelihood contribution of particle rr
    double loglik_PG=0; //likelihood contribution of the sum of all particles.
    
    
    //2. Skills at t=1
    vector<double> S1Vector; //Vector storking skills at t=1
    S1Vector.resize(RR);
    double loglik_S1rr=0; //Likelihood contribution of skills rr
    
    
    //2. Skills at t=2
    vector<double> S2Vector; //Vector storking skills at t=1
    S2Vector.resize(RR);
    double loglik_S2rr=0; //Likelihood contribution of skills rr
    
    
    
    //Defining the remaining parameters for the bootstrap filter
    vector<double> w1_rr_vector; //Vector storing weights
    w1_rr_vector.resize(RR);
    
    vector<double> what1_rr_vector; //Vector storing weights
    what1_rr_vector.resize(RR);
    
    vector<double> wtilde1_rr_vector; //Vector storing weights
    wtilde1_rr_vector.resize(RR);
    
    
    double wsum1=0;
    vector<double> wsum1_rr_vector; //Vector storing weights
    wsum1_rr_vector.resize(RR);
    
    vector<double> w2_rr_vector; //Vector storing weights
    w2_rr_vector.resize(RR);
    
    vector<double> what2_rr_vector; //Vector storing weights
    what2_rr_vector.resize(RR);
    
    vector<double> wtilde2_rr_vector; //Vector storing weights
    wtilde2_rr_vector.resize(RR);
    
    
    double wsum2=0;
    vector<double> wsum2_rr_vector; //Vector storing weights
    wsum2_rr_vector.resize(RR);
    
    
    
    
    double randomaux=0; //Auxiliary random number just in case
    double randomchosen=0; //Auxiliary random number also.
    double condition=0; //auxiliary condition
    
    
    
    //Setting seed to get draws.
    //This will generate the same random number every time.
    typedef boost::mt19937 RNGType;
    RNGType rng;
    
    boost::variate_generator<RNGType, boost::normal_distribution<> >
    generator(rng,
              boost::normal_distribution<>());
    
    //Boost to define uniform draws between 0 and 1
    boost::random::uniform_real_distribution<> unidistrib(0,1);
    
    //If I want different random numbers each time:
    //iboost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    //    generator(boost::mt19937(time(0)),
    //        boost::normal_distribution<>());
    loglik=0;
    for (int ii=0;ii<SIZE;ii=ii+1){
        //Initializing the likelihood contributions for observation ii
        logcontribobs_rr=0; //Contribution of simulation rr
        likwageM=0;
        likwageM10=0;
        likwageF=0;
        likbehavior=0;
        likbehavior10=0;
        loglikskills=0;
        loglikskills10=0;
        loglikeffm=0;
        loglikeffm10=0;
        loglikefff=0;
        loglikefff10=0;
        loglikmmu=0;
        loglikeinv=0;
        loglikeinv10=0;
        logcontribdecision=0;
        logcontribdecision10=0;
        
        //=============================
        //Definitions needed within the household:
        //=============================
        
        //Declaring the aalpha4
        aalpha4m=aalpha40m+aalpha41m*Hhchores12[ii]*1;
        aalpha4f=aalpha40f+aalpha41f*Hhchores12[ii]*1;
        aalpha4m10=aalpha40m10+aalpha41m10*Hhchores12[ii]*1;
        aalpha4f10=aalpha40f10+aalpha41f10*Hhchores12[ii]*1;
        //The price of childcare
        pricechildcare=pchildcare0+(1/(1+Hchildcareobs[ii]))*pchildcare1;
        
        //Declaring aalpha3:
        
        
        aalpha3m10=aalpha3m10-aalpha3m10_mtjh*MTJH[ii]+aalpha3mage10*Magegroup10[ii];
        aalpha3m=aalpha3m-aalpha3m12_mtjh*MTJH[ii]+aalpha3mage12*Magegroup12[ii];
        //cout << aalpha3m10 << " aalpha3m10 post" << endl;
        
        
        //========================
        //1.1 Wages likelihood
        //========================
        //1.1.1 2010
        //=======================
        //1.1.1.1  For mothers
        if (Mfraclabor10[ii]==0){
            likwageM10=0;
            Mwageinlik10=F_predwage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                    Myrschool12[ii],Mage10[ii]);
        }
        if (Mfraclabor10[ii]==1){
            likwageM10=F_likelihood_wage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                         stdwm,Myrschool12[ii],Mage10[ii],Mwage10[ii]);
            Mwageinlik10=Mwage10[ii];
        }
        
        //1.1.1.2 For fathers
        if (Ffraclabor10[ii]==0){
            likwageM10=0;
            Fwageinlik10=F_predwage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                    Fyrschool12[ii],Fage10[ii]);
        }
        if (Ffraclabor10[ii]==1){
            likwageF10=F_likelihood_wage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                         stdwf,Fyrschool12[ii],Fage10[ii],Fwage10[ii]);
            Fwageinlik10=Fwage10[ii];
        }
        
        //=======================
        //2012
        //=======================
        //1.1.1 For mother
        if (Mfraclabor12[ii]==0){
            likwageM=0;
            Mwageinlik=F_predwage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,Myrschool12[ii],Mage12[ii]);
        }
        else if (Mfraclabor12[ii]==1){
            likwageM=F_likelihood_wage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                       stdwm,Myrschool12[ii],Mage12[ii],Mwage12[ii]);
            Mwageinlik=Mwage12[ii];
        }
        //1.1.2 For father
        
        if (Ffraclabor12[ii]==0){
            likwageF=0;
            Fwageinlik=F_predwage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,Fyrschool12[ii],Fage12[ii]);
        }
        else if (Ffraclabor12[ii]==1){
            likwageF=F_likelihood_wage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                       stdwf,Fyrschool12[ii],Fage12[ii],Fwage12[ii]);
            Fwageinlik=Fwage12[ii];
        }
        
        
        //Likelihood that requires simulation
        rr=0;
        
        //Defining the mean of everything that can be set out of the loop for the
        //RR simulations
        
        //0. Single mothers 2010
        //========================================================
        
        //1.1 Effort levels
        //1.1.1 If doesn't work
        MeffortMean10NW=Meffortsolver10(aalpha2m,aalpha2m10,aalpha4m10,ggammaf,
                                        Feffort10[ii],pphi,ttheta0,ttheta2,0);
        
        //1.1.2 If works
        MeffortMean10W=Meffortsolver10(aalpha2m,aalpha2m10,aalpha4m,ggammaf,
                                       Feffort10[ii],pphi,ttheta0,ttheta2,1);
        //1.2 Investment levels
        //1.2.1 If doesn't work
        CInvestmentMean10NWNA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                             ttheta1,ttheta0,
                                             Mwageinlik10,Mnly10[ii],0,priceINV);
        
        CInvestmentMean10NWA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                            ttheta1,ttheta0,
                                            Mwageinlik10,Mnly10[ii]-pricechildcare,0,priceINV);
        
        
        //1.2.2 If works
        CInvestmentMean10WA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                           ttheta1,ttheta0,
                                           Mwageinlik10,Mnly10[ii],1,priceINV);
        CInvestmentMean10WNA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                            ttheta1,ttheta0,
                                            Mwageinlik10,Mnly10[ii]-pricechildcare,1,priceINV);
        //1. Single mothers 2012
        //========================================================
        
        //1.1 Effort levels
        //1.1.1 If doesn't work
        MeffortMean12NW=Meffortsolver(aalpha2m,aalpha4m,ggammaf,
                                      Feffort12[ii],pphi,ttheta2,0);
        //1.1.2 If works
        MeffortMean12W=Meffortsolver(aalpha2m,aalpha4m,ggammaf,
                                     Feffort12[ii],pphi,ttheta2,1);
        //1.2 Investment levels
        //1.2.1 If doesn't work
        CInvestmentMean12NW=F_investment(aalpha1m,aalpha2m,ttheta1,
                                         Mwageinlik,Mnly12[ii],0,priceINV);
        
        
        //1.2.2 If works
        CInvestmentMean12W=F_investment(aalpha1m,aalpha2m,ttheta1,
                                        Mwageinlik,Mnly12[ii],1,priceINV);
        
        
        
        //======================================================
        //Identifying the observed decision in each year
        //======================================================
        
        //2010
        //2010 Single mothers
        if (Cliveswithfather10[ii]==0 & Cliveswithmother10[ii]==1){
            //Decisions:
            //1. Mother does not work and no childcare
            //2. Mother does not work and childcare
            //3. Mother works and no childcare
            //4. Mother works and childcare
            obDecision10=(1-Mfraclabor10[ii])*(1-Cchildcare10[ii])*1+
            (1-Mfraclabor10[ii])*Cchildcare10[ii]*2+
            Mfraclabor10[ii]*(1-Cchildcare10[ii])*3+
            Mfraclabor10[ii]*Cchildcare10[ii];
        }
        
        //2010 Married couples
        if (Cliveswithfather10[ii]==1 & Cliveswithmother10[ii]==1){
            //Decisions:
            //1. Both work and childcare
            //2. Both work and no childcare
            //3. Father works, mother doesn't and childcare
            //4. Father works, mother doesn't and no childcare
            //5. Father doesn't work, mother works and childcare
            //6. Father doesn't work, mother works and no childcare
            //7. Father doesn't work, mother doesn't work and childcare
            //8. Father doesn't work, mother doesn't work and no childcare
            
            
            obDecision10=Ffraclabor10[ii]*Mfraclabor10[ii]*Cchildcare10[ii]*1+
            Ffraclabor10[ii]*Mfraclabor10[ii]*(1-Cchildcare10[ii])*2+
            Ffraclabor10[ii]*(1-Mfraclabor10[ii])*Cchildcare10[ii]*3+
            Ffraclabor10[ii]*(1-Mfraclabor10[ii])*(1-Cchildcare10[ii])*4+
            (1-Ffraclabor10[ii])*Mfraclabor10[ii]*Cchildcare10[ii]*5+
            (1-Ffraclabor10[ii])*Mfraclabor10[ii]*(1-Cchildcare10[ii])*6+
            (1-Ffraclabor10[ii])*(1-Mfraclabor10[ii])*(Cchildcare10[ii])*7+
            (1-Ffraclabor10[ii])*(1-Mfraclabor10[ii])*(1-Cchildcare10[ii])*8;
            
            
        }
        
        //2012
        //2012 Single mothers
        if (Cliveswithfather12[ii]==0 & Cliveswithmother12[ii]==1){
            //Decisions:
            //1. Mother works
            //0. Mother does not work
            obDecision=(Mfraclabor12[ii]==1);
        }
        
        //2012 Married couples
        if (Cliveswithfather12[ii]==1 & Cliveswithmother12[ii]==1){
            //Decisions:
            //1. Both work
            //2. Father works, mother doesn't
            //3. Mother works, father doesn't
            //4. Neither father or mother work.
            
            obDecision=Ffraclabor12[ii]*Mfraclabor12[ii]*1+
            Ffraclabor12[ii]*(1-Mfraclabor12[ii])*2+
            (1-Ffraclabor12[ii])*Mfraclabor12[ii]*3+
            (1-Ffraclabor12[ii])*(1-Mfraclabor12[ii])*4;
        }
        //cout << " here1 "<< endl;
        
        
        
        //Start of simulations. I will start with the particles that need
        // to be drawn in t=0. These are the skills at period zero or at birth,
        //given by S0, and the skills of the primary caregiver. In first stance
        //I will assume just normal distributions for both. Analyzing the distribution of
        //both variables in STATA, the assumption of standard normal does not seem unreasonable.
        
        
        //==========================================
        //Draws of the birth outcomes and PG skills
        //==========================================
        
        rr=0;
        loglik_S0rr=0;
        loglik_PGrr=0;
        loglik_S1rr=0;
        loglik_S2rr=0;
        wsum1=0;
        wsum2=0;
        
        while (rr<RR){
            //Skills at birth
            S0rr=gen_normal_3(generator);
            S0Vector[rr]=S0rr;
            
            //Skills of primary caregiver
            PGrr=gen_normal_3(generator);
            PGVector[rr]=PGrr;
            
            //Likelihood of skills at birth and skills of primary caregiver
            loglik_S0rr+=F_likgeneric(S0rr,Cbirthfactor[ii],STDS0);
            loglik_PGrr+=F_likgeneric(PGrr,Ccareskills12[ii],STDPG);
            
            rr=rr+1;
        }
        
        
        //========================================================
        //2010 simulations
        //========================================================
        rr=0;
        rr=0;
        loglik_S0rr=0;
        loglik_PGrr=0;
        loglik_S1rr=0;
        loglik_S2rr=0;
        wsum1=0;
        wsum2=0;
        
        
        while (rr<RR){
            logcontribobs_rr10=0;
            loglikskills_rr10=0;
            loglikeffm_rr10=0;
            loglikefff_rr10=0;
            loglikeinv_rr10=0;
            Decision10=0;
            logcontribdecision_rr10=0;
            
            //1. Draw shocks from the corresponding distributions
            //===========================================================
            //1.1 Single mothers:
            if (Cliveswithfather10[ii]==0 & Cliveswithmother10[ii]==1){
                //cout << " here10 "<< endl;
                //========================================================
                //1. Effort levels
                //========================================================
                //Draw a shock for the effort levels
                Meffort_factor10shock=gen_normal_3(generator);
                
                //========================================================
                //1.1 If mom does'nt work
                //========================================================
                // Getting a draw from a standard normal distribution
                
                // Now transforming it as a draw from the given distribution
                Meffort_factor10NW=Meffort_factor10shock*stdeffortmom;
                
                Meffort_factor10NW=exp(Meffort_factor10NW)*MeffortMean10NW;
                
                //Effort can't be negative, we are drawing from a truncated normal distribution
                if (Meffort_factor10NW<0){
                    Meffort_factor10NW=0;
                }
                //cout << " here11 "<< endl;
                //We also bound the effort from above
                if (Meffort_factor10NW>=1.0e+6){
                    Meffort_factor10NW=1.0e+6;
                }
                //========================================================
                //1.2 If mom works
                // Getting a draw from a standard normal distribution
                
                // Now transforming it as a draw from the given distribution
                
                Meffort_factor10W=Meffort_factor10shock*stdeffortmom;
                
                Meffort_factor10W=exp(Meffort_factor10W)*MeffortMean10W;
                //cout << " here12 "<< endl;
                //Effort can't be negative, we are drawing from a truncated normal distribution
                if (Meffort_factor10W<0){
                    Meffort_factor10W=0;
                }
                if (Meffort_factor10W>=1.0e+6){
                    Meffort_factor10W=1.0e+6;
                }
                //========================================================
                
                
                //========================================================
                //2. Investment levels
                //========================================================
                //Getting the draw for the investment shock
                CInv_factor10shock=gen_normal_3(generator);
                //========================================================
                //2.1 If mom doesn't work
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                //CHILDCARE
                CInv_factor10NWA=CInv_factor10shock*stdinvestment;
                CInv_factor10NWA=CInv_factor10NWA+CInvestmentMean10NWA;
                //cout << " here13 "<< endl;
                //Investment can't be negative
                if (CInv_factor10NWA<0){
                    CInv_factor10NWA=0;
                }
                
                //NO CHILDCARE
                CInv_factor10NWNA=CInv_factor10shock*stdinvestment;
                CInv_factor10NWNA=CInv_factor10NWNA+CInvestmentMean10NWNA;
                //cout << " here13 "<< endl;
                //Investment can't be negative
                if (CInv_factor10NWNA<0){
                    CInv_factor10NWNA=0;
                }
                //========================================================
                //2.2 If mom works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                //CHILDCARE
                CInv_factor10WA=CInv_factor10shock*stdinvestment;
                CInv_factor10WA=CInv_factor10WA+CInvestmentMean10WA;
                //cout << " here14 "<< endl;
                //Investment can't be negative
                if (CInv_factor10WA<0){
                    CInv_factor10WA=0;
                }
                
                
                //NO CHILDCARE
                CInv_factor10WNA=CInv_factor10shock*stdinvestment;
                CInv_factor10WNA=CInv_factor10WNA+CInvestmentMean10WNA;
                //cout << " here14 "<< endl;
                //Investment can't be negative
                if (CInv_factor10WNA<0){
                    CInv_factor10WNA=0;
                }
                //==========================================================
                //3. Skills levels
                //==========================================================
                //Shock of skills
                CSkills_factor10shock=gen_normal_3(generator);
                //========================================================
                //1. If mom doesn't works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                //No childcare
                CSkillsMean10NWNA=F_predskills(ddelta0,ddelta1,
                                               ddelta2,ddelta3,ddelta4,
                                               Cedad_meses10[ii],ttheta0,
                                               ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                               Feffort10[ii],Meffort_factor10NW,
                                               CInv_factor10NWNA,
                                               S0Vector[rr],
                                               0,
                                               PGVector[rr],Hmemberstotal10[ii]);
                
                //Notice that it is no longer
                //Cbirthfactor[ii] and PG[ii] as we
                //are doing the bootstrap filter.
                
                if (CSkills_factor10shock*stdskills<=500){
                    CSkills_factor10NWNA=CSkills_factor10shock*stdskills;
                }
                if (CSkills_factor10shock*stdskills>=500){
                    CSkills_factor10NWNA=500;
                }
                
                CSkills_factor10NWNA=exp(CSkills_factor10NWNA)*CSkillsMean10NWNA;
                //Childcare
                CSkillsMean10NWA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                              Cedad_meses10[ii],ttheta0,
                                              ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                              Feffort10[ii],Meffort_factor10NW,
                                              CInv_factor10NWA,S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii]);
                //cout << " here15 "<< endl;
                if (CSkills_factor10shock*stdskills<=500){
                    CSkills_factor10NWA=CSkills_factor10shock*stdskills;
                }
                if (CSkills_factor10shock*stdskills>=500){
                    CSkills_factor10NWA=500;
                }
                
                CSkills_factor10NWA=exp(CSkills_factor10NWA)*CSkillsMean10NWA;
                //========================================================
                //2. If mom works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                //no childcare
                CSkillsMean10WNA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                              Cedad_meses10[ii],ttheta0,
                                              ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                              Feffort10[ii],Meffort_factor10W,
                                              CInv_factor10WNA,S0Vector[rr],0,PGVector[rr],Hmemberstotal10[ii]);
                //cout << " here16 "<< endl;
                if (CSkills_factor10shock*stdskills<=500){
                    CSkills_factor10WNA=CSkills_factor10shock*stdskills;
                }
                if (CSkills_factor10shock*stdskills>=500){
                    CSkills_factor10WNA=500;
                }
                
                
                CSkills_factor10WNA=exp(CSkills_factor10WNA)*CSkillsMean10WNA;
                //childcare
                CSkillsMean10WA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                             Cedad_meses10[ii],ttheta0,
                                             ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                             Feffort10[ii],Meffort_factor10W,
                                             CInv_factor10WA,S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii]);
                //cout << " here16 "<< endl;
                if (CSkills_factor10shock*stdskills<=500){
                    CSkills_factor10WA=CSkills_factor10shock*stdskills;
                }
                if (CSkills_factor10shock*stdskills>=500){
                    CSkills_factor10WA=500;
                }
                
                
                CSkills_factor10WA=exp(CSkills_factor10WA)*CSkillsMean10WA;
                //cout << " here17 "<< endl;
                //===========================================================
                //===========================================================
                //4 Preference shocks
                //===========================================================
                //1. If mom doesn't work
                //===========================================================
                //No childcare
                MprefshockNWNA10=gen_normal_3(generator);
                MprefshockNWNA10=MprefshockWNA10*MshockNWNA;
                //Childcare
                MprefshockNWA10=gen_normal_3(generator);
                MprefshockNWA10=MprefshockWA10*MshockNWA;
                //===========================================================
                //2. If mom works
                //===========================================================
                //No childcare
                MprefshockWNA10=gen_normal_3(generator);
                MprefshockWNA10=MprefshockWNA10*MshockWNA;
                //childcare
                MprefshockWA10=gen_normal_3(generator);
                MprefshockWA10=MprefshockWA10*MshockWA;
                //===========================================================
                //===========================================================
                //5 Consumption
                //===========================================================
                //1. If mom doesn't work
                //===========================================================
                //No childcare
                MconsumptionNWNA10=Mnly10[ii]-priceINV*CInv_factor10NWNA;
                if (MconsumptionNWNA10<0){
                    MconsumptionNWNA10=1.0e-10;
                }
                //childcare
                MconsumptionNWA10=Mnly10[ii]-priceINV*CInv_factor10NWA-pricechildcare;
                if (MconsumptionNWA10<0){
                    MconsumptionNWA10=1.0e-10;
                }
                //===========================================================
                //2. If mom works
                //===========================================================
                //no childcare
                MconsumptionWNA10=Mnly10[ii]-priceINV*CInv_factor10NWNA+Mwageinlik10;
                if (MconsumptionWNA10<0){
                    MconsumptionWNA10=1.0e-10;
                }
                //childcare
                MconsumptionWA10=Mnly10[ii]-priceINV*CInv_factor10NWA+Mwageinlik10-pricechildcare;
                if (MconsumptionWA10<0){
                    MconsumptionWA10=1.0e-10;
                }
                
                
                
                //cout << " here18 "<< endl;
                //===========================================================
                //===========================================================
                //6. Behavioral model
                //==========================================================
                //1. If mom doesn't work
                //==========================================================
                //no childcare
                MutilityNWNA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                           aalpha4m10,aalpha5m10,
                                           MconsumptionNWNA10,
                                           Meffort_factor10NW,
                                           0,CSkills_factor10NWNA,0)+MprefshockNWNA10;
                //childcare
                MutilityNWA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                          aalpha4m10,aalpha5m10,
                                          MconsumptionNWA10,
                                          Meffort_factor10NW,
                                          0,CSkills_factor10NWA,1)+MprefshockNWA10;
                //===========================================================
                //2. If mom works
                //===========================================================
                
                //no childcare
                MutilityWNA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                          aalpha4m,aalpha5m10,
                                          MconsumptionWNA10,
                                          Meffort_factor10W,
                                          1,CSkills_factor10WNA,0)+MprefshockWNA10;
                //CHILDCARE
                MutilityWA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                         aalpha4m,aalpha5m10,
                                         MconsumptionWA10,
                                         Meffort_factor10W,1,CSkills_factor10WA,
                                         1)+MprefshockWA10;
                //cout << " here19 "<< endl;
                //===========================================================
                //===========================================================
                //cout << " here20 "<< endl
                //Identifying the decision
                if ((MutilityNWNA10> MutilityNWA10) &&
                    (MutilityNWNA10> MutilityWNA10) &&
                    (MutilityNWNA10> MutilityWA10)){
                    Decision10=1;
                }
                else if ((MutilityNWA10> MutilityWNA10) &&
                         (MutilityNWA10> MutilityWA10)){
                    Decision10=2;
                }
                else if (MutilityWNA10> MutilityWA10){
                    Decision10=3;
                }
                else{
                    Decision10=4;
                }
                
                Decision10=(Decision10==obDecision10);
                
                //Likelihood contributions of the measurement systems
                if (Mfraclabor10[ii]==1 && Cchildcare10[ii]==0){
                    
                    
                    
                    
                    
                    
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WNA,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WNA;
                    
                    
                    //2. Effort
                    loglikeffm_rr10=
                    F_loglikelihood_generic(Meffort10[ii],
                                            Meffort_factor10W,MEASMeffort);
                    
                    loglikefff_rr10=0;
                    
                    //3. Investment
                    loglikeinv_rr10=
                    F_loglikelihood_generic(Cfactorinv10[ii],CInv_factor10WNA,
                                            MEASINV);
                    
                    
                }
                
                if (Mfraclabor10[ii]==1 && Cchildcare10[ii]==1){
                    
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WA,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WA;
                    
                    
                    
                    //2. Effort
                    loglikeffm_rr10=
                    F_loglikelihood_generic(Meffort10[ii],
                                            Meffort_factor10W,MEASMeffort);
                    
                    loglikefff_rr10=0;
                    
                    //3. Investment
                    loglikeinv_rr10=
                    F_loglikelihood_generic(Cfactorinv10[ii],CInv_factor10WA,
                                            MEASINV);
                    
                    
                }
                if (Mfraclabor10[ii]==0 && Cchildcare10[ii]==0){
                    
                    //1. Skills
                    
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                            CSkills_factor10NWNA,MEASSkills);
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWNA;
                    
                    
                    
                    
                    
                    //Storing the particle weight
                    what1_rr_vector[rr]=loglikskills_rr10;
                    
                    
                    //2. Effort
                    loglikeffm_rr=
                    F_loglikelihood_generic(Meffort10[ii],
                                            Meffort_factor10NW,MEASMeffort);
                    
                    
                    loglikefff_rr=0;
                    
                    //3. Investment
                    loglikeinv_rr=
                    F_loglikelihood_generic(Cfactorinv10[ii],CInv_factor10NWNA,
                                            MEASINV);
                    
                    
                }
                if (Mfraclabor10[ii]==0 && Cchildcare10[ii]==1){
                    
                    //1. Skills
                    
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                            CSkills_factor10NWA,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWA;
                    
                    
                    
                    //Storing the particle weight
                    what1_rr_vector[rr]=loglikskills_rr10;
                    
                    
                    //2. Effort
                    loglikeffm_rr=
                    F_loglikelihood_generic(Meffort10[ii],
                                            Meffort_factor10NW,MEASMeffort);
                    //cout << " here28 "<< endl;
                    
                    loglikefff_rr=0;
                    
                    //3. Investment
                    loglikeinv_rr=
                    F_loglikelihood_generic(Cfactorinv10[ii],CInv_factor10NWA,
                                            MEASINV);
                    //cout << " here29 "<< endl;
                    
                }
                
                
                
            }//End of single mothers
            //If both parents are present.
            
            
            if (Cliveswithfather10[ii]==1 & Cliveswithmother10[ii]==1){
                
                //======================
                //0. Mmu
                //======================
                //Getting the draw for the unobserved heterogeneity of mmy
                MMushock10=gen_normal_3(generator);
                MMushock10=MMushock10*stdmmu;
                //Getting the predicted level of bargaining power
                mupred10=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                               llambda5,llambda6,llambda7,llambda8,
                               Fwageinlik10,Mwageinlik10,Fnly10[ii],Mnly10[ii],MMushock,mmuLB,mmuUB,
                               Fage10[ii],Mage10[ii],Fyrschool12[ii],Myrschool12[ii],
                               MRATIO[ii],Unemployment[ii],Wageratio[ii],Distance[ii]);
                
                
                
                //========================================================
                //1. Effort levels
                //========================================================
                
                //Draw a shock for the effort levels
                Meffort_factor10shock=exp(gen_normal_3(generator)*stdeffortmom);
                Feffort_factor10shock=exp(gen_normal_3(generator)*stdeffortfat);
                //Getting the means
                //1.1 Effort levels
                //1.1.1 NWNW
                MeffortMean10NWNW_TOG=F_effort_m10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f10,aalpha4m10,0,0,pphi,ttheta0,
                 ttheta2,0.92)*Meffort_factor10shock;
                FeffortMean10NWNW_TOG=F_effort_f10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f,aalpha4m,0,0,pphi,ttheta0,
                 ttheta2,0.92)*Feffort_factor10shock;
                
                
                //1.1.2 Father works mother not
                MeffortMean10WNW_TOG=F_effort_m10
                (mupred10,ggammaf,aalpha2f,aalpha2m,
                 aalpha4f,aalpha4m,1,0,pphi,ttheta0,
                 ttheta2,0.92)*Meffort_factor10shock;
                FeffortMean10WNW_TOG=F_effort_f10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f,aalpha4m,1,0,pphi,ttheta0,
                 ttheta2,0.92)*Feffort_factor10shock;
                
                
                //1.1.3 Mother works father doesn't
                MeffortMean10NWW_TOG=F_effort_m10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
                 ttheta2,0.92)*Meffort_factor10shock;
                
                FeffortMean10NWW_TOG=F_effort_f10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
                 ttheta2,0.92)*Feffort_factor10shock;
                
                //1.1.4 Both work
                MeffortMean10WW_TOG=F_effort_m10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
                 ttheta2,0.92)*Meffort_factor10shock;
                FeffortMean10WW_TOG=F_effort_f10
                (mupred10,ggammaf,aalpha2f10,aalpha2m10,
                 aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
                 ttheta2,0.92)*Feffort_factor10shock;
                
                //Need to fix the effort levels either if they are negative or
                //if they are too large
                //Effort can't be negative, we are drawing from a truncated normal distribution
                //Fix the effort levels
                //Mother fixing
                if (MeffortMean10NWNW_TOG<0){
                    MeffortMean10NWNW_TOG=0;
                }
                if (MeffortMean10NWNW_TOG>=1.0e+6){
                    MeffortMean10NWNW_TOG=1.0e+6;
                }
                
                
                if (MeffortMean10WNW_TOG<0){
                    MeffortMean10WNW_TOG=0;
                }
                if (MeffortMean10WNW_TOG>=1.0e+6){
                    MeffortMean10WNW_TOG=1.0e+6;
                }
                
                
                if (MeffortMean10NWW_TOG<0){
                    MeffortMean10NWW_TOG=0;
                }
                if (MeffortMean10NWW_TOG>=1.0e+6){
                    MeffortMean10NWW_TOG=1.0e+6;
                }
                
                
                if (MeffortMean10WW_TOG<0){
                    MeffortMean10WW_TOG=0;
                }
                if (MeffortMean10WW_TOG>=1.0e+6){
                    MeffortMean10WW_TOG=1.0e+6;
                }
                //Father fixing
                if (FeffortMean10NWNW_TOG<0){
                    FeffortMean10NWNW_TOG=0;
                }
                if (FeffortMean10NWNW_TOG>=1.0e+6){
                    FeffortMean10NWNW_TOG=1.0e+6;
                }
                
                
                if (FeffortMean10WNW_TOG<0){
                    FeffortMean10WNW_TOG=0;
                }
                if (FeffortMean10WNW_TOG>=1.0e+6){
                    FeffortMean10WNW_TOG=1.0e+6;
                }
                
                
                if (FeffortMean10NWW_TOG<0){
                    FeffortMean10NWW_TOG=0;
                }
                if (FeffortMean10NWW_TOG>=1.0e+6){
                    FeffortMean10NWW_TOG=1.0e+6;
                }
                
                
                if (FeffortMean10WW_TOG<0){
                    FeffortMean10WW_TOG=0;
                }
                if (FeffortMean10WW_TOG>=1.0e+6){
                    FeffortMean10WW_TOG=1.0e+6;
                }
                
                
                //========================================================
                //2. Investment levels
                //========================================================
                //Getting the draw for the investment shock
                CInv_factor10shock=gen_normal_3(generator);
                //========================================================
                //2.1 If mom doesn't work
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CInv_factor10shock=CInv_factor10shock*stdinvestment;
                //4. Now get the mean of the investment as a function of mmu
                //for the possible combinations
                //1. NWNWNA
                CInvestmentMean10NWNWNA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,0,0,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii],ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10NWNWNA_TOG<=0){
                    CInvestmentMean10NWNWNA_TOG=1.0e-5;
                }
                
                //2. NWNWA
                CInvestmentMean10NWNWA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,0,0,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii]-pricechildcare,ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10NWNWA_TOG<=0){
                    CInvestmentMean10NWNWA_TOG=1.0e-5;
                }
                
                //3. WNWNA
                CInvestmentMean10WNWNA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,1,0,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii],ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10WNWNA_TOG<=0){
                    CInvestmentMean10WNWNA_TOG=1.0e-5;
                }
                
                //4. WNWA
                CInvestmentMean10WNWNA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,1,0,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii]-pricechildcare,ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10WNWA_TOG<=0){
                    CInvestmentMean10WNWA_TOG=1.0e-5;
                }
                
                //5. NWWNA
                CInvestmentMean10NWWNA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,0,1,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii],ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10NWWNA_TOG<=0){
                    CInvestmentMean10NWWNA_TOG=1.0e-5;
                }
                
                //6. NWWNA
                CInvestmentMean10NWWA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,0,1,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii]-pricechildcare,ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10NWWA_TOG<=0){
                    CInvestmentMean10NWWA_TOG=1.0e-5;
                }
                
                //7. WWNA
                CInvestmentMean10WWNA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,1,1,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii],ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10WWNA_TOG<=0){
                    CInvestmentMean10WWNA_TOG=1.0e-5;
                }
                
                
                //8. WWA
                CInvestmentMean10WWA_TOG=F_invcouple10
                (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
                 aalpha2f10,aalpha2m10,
                 mupred10,1,1,Fwageinlik10,Mwageinlik10,
                 Fnly10[ii],Mnly10[ii]-pricechildcare,ttheta0,
                 ttheta1,priceINV)+CInv_factor10shock;
                if (CInvestmentMean10WWA_TOG<=0){
                    CInvestmentMean10WWA_TOG=1.0e-5;
                }
                //==========================================================
                //3. Skills levels
                //==========================================================
                //Shock of skills
                CSkills_factor10shock=gen_normal_3(generator);
                if (CSkills_factor10shock*stdskills<500){
                    CSkills_factor10shock=exp(CSkills_factor10shock*stdskills);
                }
                if (CSkills_factor10shock*stdskills>=500){
                    CSkills_factor10shock=1000;
                }
                
                //Now getting the different levels of skills for all possi-
                //bilities
                //3.1 NWNWNA
                CSkills_factor10NWNWNA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10NWNW_TOG,
                 MeffortMean10NWNW_TOG,
                 CInvestmentMean10NWNWNA_TOG,
                 S0Vector[rr],0,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                
                
                //3.1 NWNWA
                CSkills_factor10NWNWA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10NWNW_TOG,
                 MeffortMean10NWNW_TOG,
                 CInvestmentMean10NWNWA_TOG,
                 S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                
                //3.2 WNWNA
                CSkills_factor10WNWNA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10WNW_TOG,
                 MeffortMean10WNW_TOG,
                 CInvestmentMean10WNWNA_TOG,
                 S0Vector[rr],0,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                //3.2 WNWA
                CSkills_factor10WNWA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10WNW_TOG,
                 MeffortMean10WNW_TOG,
                 CInvestmentMean10WNWA_TOG,
                 S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                
                //3.3 NWWNA
                CSkills_factor10NWWNA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10NWW_TOG,
                 MeffortMean10NWW_TOG,
                 CInvestmentMean10NWWNA_TOG,
                 S0Vector[rr],0,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                //3.3 NWWA
                CSkills_factor10NWWA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10NWW_TOG,
                 MeffortMean10NWW_TOG,
                 CInvestmentMean10NWWA_TOG,
                 S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                
                //3.4 WWNA
                CSkills_factor10WWNA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10WW_TOG,
                 MeffortMean10WW_TOG,
                 CInvestmentMean10WWNA_TOG,
                 S0Vector[rr],0,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                //3.4 WWNA
                
                CSkills_factor10WWA_TOG=F_predskills
                (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10[ii],
                 ttheta0,ttheta1,ttheta2,pphi,
                 ggammaf,ggammam,FeffortMean10WW_TOG,
                 MeffortMean10WW_TOG,
                 CInvestmentMean10WWA_TOG,
                 S0Vector[rr],1,PGVector[rr],Hmemberstotal10[ii])*CSkills_factor10shock;
                
                //===========================================================
                //===========================================================
                //4 Preference shocks
                
                //1. If mom doesn't work
                MprefshockNWNA10=gen_normal_3(generator);
                MprefshockNWNA10=MprefshockNWNA10*MshockNWNA;
                
                MprefshockNWA10=gen_normal_3(generator);
                MprefshockNWA10=MprefshockNWA10*MshockNWA;
                
                
                //2. If mom works
                MprefshockWNA10=gen_normal_3(generator);
                MprefshockWNA10=MprefshockWNA10*MshockWNA;
                
                MprefshockWA10=gen_normal_3(generator);
                MprefshockWA10=MprefshockWA10*MshockWA;
                
                //3. If FATHER doesn't work
                FprefshockNWNA10=gen_normal_3(generator);
                FprefshockNWNA10=FprefshockWNA10*FshockNWNA;
                
                FprefshockWNA10=gen_normal_3(generator);
                FprefshockWNA10=FprefshockWNA10*FshockWNA;
                
                //4. If FATHER works
                FprefshockWNA10=gen_normal_3(generator);
                FprefshockWNA10=FprefshockWNA10*FshockWNA;
                
                FprefshockWA10=gen_normal_3(generator);
                FprefshockWA10=FprefshockWA10*FshockWA;
                
                
                //===========================================================
                //===========================================================
                //5 Consumption
                //===========================================================
                //1. Mother
                //===========================================================
                
                //1. NWNW
                MconsumptionNWNWNA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                MconsumptionNWNWA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                //2. WNW
                MconsumptionWNWNA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                MconsumptionWNWA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                //3. NWW
                MconsumptionNWWNA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                MconsumptionNWWA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                //4. WW
                MconsumptionWWNA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                MconsumptionWWA_TOG10=M_consumption_TOG10
                (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                //===========================================================
                //1. Father
                //===========================================================
                
                //1. NWNW
                FconsumptionNWNWNA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                FconsumptionNWNWA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                //2. WNW
                FconsumptionWNWNA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                FconsumptionWNWA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                //3. NWW
                FconsumptionNWWNA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                FconsumptionNWWA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10NWWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                //4 WW
                FconsumptionWWNA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                FconsumptionWWA_TOG10=F_consumption_TOG10
                (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
                 CInvestmentMean10WWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
                
                //===========================================================
                //===========================================================
                //6. Behavioral model
                //==========================================================
                //1. Mother
                //==========================================================
                //1. NWNW
                
                MutilityNWNWNA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,
                 aalpha4m10,aalpha5m10,MconsumptionNWNWNA_TOG10,MeffortMean10NWNW_TOG,0,CSkills_factor10NWNWNA_TOG,0)
                +MprefshockNWNA10;
                
                MutilityNWNWA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,
                 aalpha4m10,aalpha5m10,MconsumptionNWNWA_TOG10,MeffortMean10NWNW_TOG,0,CSkills_factor10NWNWA_TOG,1)
                +MprefshockNWA10;
                
                //2. WNW
                
                MutilityWNWNA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionWNWNA_TOG10,MeffortMean10WNW_TOG,0,
                 CSkills_factor10WNWNA_TOG,0)+MprefshockNWNA10;
                
                MutilityWNWA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionWNWA_TOG10,MeffortMean10WNW_TOG,0,
                 CSkills_factor10WNWA_TOG,1)+MprefshockNWA10;
                
                
                //3. NWW
                
                MutilityNWWNA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionNWWNA_TOG10,MeffortMean10NWW_TOG,1,
                 CSkills_factor10NWWNA_TOG,0)+MprefshockWNA10;
                
                MutilityNWWA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionNWWA_TOG10,MeffortMean10NWW_TOG,1,
                 CSkills_factor10NWWA_TOG,1)+MprefshockWA10;
                
                //4. WW
                MutilityWWNA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionWWNA_TOG10, MeffortMean10WW_TOG,1,
                 CSkills_factor10WWNA_TOG,0)+MprefshockWNA10;
                
                MutilityWWA_TOG10=F_utility10
                (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
                 MconsumptionWWA_TOG10, MeffortMean10WW_TOG,1,
                 CSkills_factor10WWA_TOG,1)+MprefshockWA10;
                //==========================================================
                //2. Father
                //==========================================================
                
                //1. NWNW
                FutilityNWNWNA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionNWNWNA_TOG10, FeffortMean10NWNW_TOG,0,
                 CSkills_factor10NWNWNA_TOG,0)+ FprefshockNWNA10;
                
                FutilityNWNWA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionNWNWA_TOG10, FeffortMean10NWNW_TOG,0,
                 CSkills_factor10NWNWA_TOG,1)+ FprefshockNWA10;
                
                //2. WNW
                
                FutilityWNWNA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionWNWNA_TOG10,FeffortMean10WNW_TOG,1,
                 CSkills_factor10WNWNA_TOG,0)+FprefshockWNA10;
                
                FutilityWNWA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionWNWA_TOG10,FeffortMean10WNW_TOG,1,
                 CSkills_factor10WNWA_TOG,1)+FprefshockWA10;
                //3.NWW
                FutilityNWWNA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionNWWNA_TOG,FeffortMean10NWW_TOG,0,
                 CSkills_factor10NWWNA_TOG,0)+FprefshockNWNA10;
                
                FutilityNWWA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionNWWA_TOG10,FeffortMean10NWW_TOG,0,
                 CSkills_factor10NWWA_TOG,1)+FprefshockNWA10;
                
                //4. WW
                FutilityWWNA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionWWNA_TOG10,FeffortMean12WW_TOG,1,
                 CSkills_factor10WWNA_TOG,0)+FprefshockWNA10;
                
                FutilityWWA_TOG10=F_utility10
                (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
                 FconsumptionWWA_TOG10,FeffortMean12WW_TOG,1,
                 CSkills_factor10WWA_TOG,1)+FprefshockWA10;
                
                //Seeing the utility levels for each case:
                
                
                //==========================================================
                // Defining the utilities for the possible combinations
                //==========================================================
                
                WelfareNWNWNA10=mupred10*FutilityNWNWNA_TOG10+
                (1-mupred10)*MutilityNWNWNA_TOG10;
                
                WelfareNWNWA10=mupred10*FutilityNWNWA_TOG10+
                (1-mupred10)*MutilityNWNWA_TOG10;
                
                
                WelfareWNWNA10=mupred10*FutilityWNWNA_TOG10+
                (1-mupred10)*MutilityWNWNA_TOG10;
                
                WelfareWNWA10=mupred10*FutilityWNWA_TOG10+
                (1-mupred10)*MutilityWNWA_TOG10;
                
                WelfareNWWNA10=mupred10*FutilityNWWNA_TOG10+
                (1-mupred10)*MutilityNWWNA_TOG10;
                
                WelfareNWWA10=mupred10*FutilityNWWA_TOG10+
                (1-mupred10)*MutilityNWWA_TOG10;
                
                WelfareWWNA10=mupred10*FutilityWWNA_TOG10+
                (1-mupred10)*MutilityWWNA_TOG10;
                
                WelfareWWA10=mupred10*FutilityWWA_TOG10+
                (1-mupred10)*MutilityWWA_TOG10;
                
                
                
                //==========================================================
                //Identifying the correct decision
                //==========================================================
                if ( (WelfareNWNWNA10>WelfareNWNWA10) &&
                    ( WelfareNWNWNA10>WelfareNWWNA10)  &&
                    ( WelfareNWNWNA10>WelfareNWWA10)  &&
                    ( WelfareNWNWNA10>WelfareWNWNA10)  &&
                    ( WelfareNWNWNA10>WelfareWNWA10)  &&
                    ( WelfareNWNWNA10>WelfareWWNA10)  &&
                    ( WelfareNWNWNA10>WelfareWWA10)) {
                    Decision10=8;
                }
                else if ((WelfareNWNWA10>WelfareNWWNA10) &&
                         (WelfareNWNWA10>WelfareNWWA10) &&
                         (WelfareNWNWA10>WelfareWNWNA10) &&
                         (WelfareNWNWA10>WelfareWNWA10) &&
                         (WelfareNWNWA10>WelfareWWNA10) &&
                         (WelfareNWNWA10>WelfareWWA10)) {
                    Decision10=7;
                }
                
                else if ((WelfareNWWNA10>WelfareNWWA10) &&
                         (WelfareNWWNA10>WelfareWNWNA10) &&
                         (WelfareNWWNA10>WelfareWNWA10) &&
                         (WelfareNWWNA10>WelfareWWNA10) &&
                         (WelfareNWWNA10>WelfareWWA10)) {
                    Decision10=6;
                }
                
                else if ((WelfareNWWA10>WelfareWNWNA10) &&
                         (WelfareNWWA10>WelfareWNWA10) &&
                         (WelfareNWWA10>WelfareWWNA10) &&
                         (WelfareNWWA10>WelfareWWA10)) {
                    Decision10=5;
                }
                
                else if ((WelfareWNWNA10>WelfareWNWA10) &&
                         (WelfareWNWNA10>WelfareWWNA10) &&
                         (WelfareWNWNA10>WelfareWWA10)) {
                    Decision10=4;
                }
                
                else if ((WelfareWNWA10>WelfareWWNA10) &&
                         (WelfareWNWA10>WelfareWWA10)) {
                    Decision10=3;
                }
                
                else if ((WelfareWWNA10>WelfareWWA10)) {
                    Decision10=2;
                }
                
                else{
                    Decision10=1;
                }
                
                
                
                
                Decision10=(Decision10==obDecision10);
                //7. Likelihood of measurement system
                //===============================
                
                //7.1 NWNWNA
                if (Mfraclabor10[ii]==0 && Ffraclabor10[ii]==0 && Cchildcare10[ii]==0){
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10NWNWNA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWNWNA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10NWNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10NWNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10NWNWNA_TOG,MEASINV);
                }
                
                //7.2 NWNWA
                if (Mfraclabor10[ii]==0 && Ffraclabor10[ii]==0 && Cchildcare10[ii]==1){
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10NWNWA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWNWA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10NWNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10NWNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10NWNWA_TOG,MEASINV);
                }
                
                
                //7.3 WNWNA
                if (Mfraclabor10[ii]==0 && Ffraclabor10[ii]==1 && Cchildcare10[ii]==0){
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WNWNA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WNWNA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10WNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10WNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10WNWNA_TOG,MEASINV);
                }
                
                
                //7.4 WNWA
                if (Mfraclabor10[ii]==0 && Ffraclabor10[ii]==1 && Cchildcare10[ii]==1){
                    
                    
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WNWA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WNWA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10WNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10WNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10WNWA_TOG,MEASINV);
                }
                
                
                //7.5 NWWNA
                if (Mfraclabor10[ii]==1 && Ffraclabor10[ii]==0 && Cchildcare10[ii]==0){
                    //1. Skills
                    
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10NWWNA_TOG,MEASSkills);
                    
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWWNA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10NWW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10NWW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10NWWNA_TOG,MEASINV);
                }
                
                
                //7.6 NWWA
                if (Mfraclabor10[ii]==1 && Ffraclabor10[ii]==0 && Cchildcare10[ii]==1){
                    //1. Skills
                    
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10NWWA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10NWWA_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10NWW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10NWW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10NWWA_TOG,MEASINV);
                }
                //7.7 WWNA
                if (Mfraclabor10[ii]==1 && Ffraclabor10[ii]==1 && Cchildcare10[ii]==0){
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WWNA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WWNA_TOG;
                    
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10WW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10WW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10WWNA_TOG,MEASINV);
                    
                }
                
                //7.7 WWA
                if (Mfraclabor10[ii]==1 && Ffraclabor10[ii]==1 && Cchildcare10[ii]==1){
                    //1. Skills
                    loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10[ii],
                                                              CSkills_factor10WWA_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what1_rr_vector[rr]=loglikskills_rr10;
                    S1Vector[rr]=CSkills_factor10WWA_TOG;
                    
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr10=F_loglikelihood_generic(Meffort10[ii],
                                                            MeffortMean10WW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr10=F_loglikelihood_generic(Feffort10[ii],
                                                            FeffortMean10WW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr10=F_loglikelihood_generic(Cfactorinv10[ii],CInvestmentMean10WWA_TOG,MEASINV);
                    
                }
                
                
                //5. Measurement system for mmu. No in 2010
                
            }//If living together
            
            //Computing ww: weights for resampling of period 1.
            w1_rr_vector[rr]=what1_rr_vector[rr]*(1/RR);
            
            //Computing the sum of the weights
            
            wsum1+=w1_rr_vector[rr];
            //Adding up likelihood contributions of simulations
            
            
            //0 Decisions (behavioral system)
            logcontribdecision10+=Decision10; //Storing the number of times correct thing
            //1. Skills
            
            loglikskills10+=loglikskills_rr10;
            
            
            //2. Mother's effort
            loglikeffm10+=loglikeffm_rr10;
            //3. Father's effort
            loglikefff10+=loglikefff_rr10;
            //4. Investment
            loglikeinv10+=loglikeinv_rr10;
            
            //5. Mmu
            
            rr=rr+1;
        } //Finish the RR simulations
        
        //We need to re-draw the particles with the given weights specified previously.
        //First, draw uniform number between 0 and the sum of the weights
        randomaux=unidistrib(rng)*wsum1;
        
        
        //Getting RR draws
        for (int mm=0; mm<RR; mm=mm+1){
            randomaux=unidistrib(rng)*wsum1;
            rr=0;
            while (condition==0){
                if (randomaux<w1_rr_vector[rr]){
                    randomchosen=S1Vector[rr];
                    condition=1;
                    S1Vector[rr]=randomchosen;
                }
                randomaux-=w1_rr_vector[rr];
                rr=rr+1;
            }
            condition=0;
        }
        
        
        //=====
        //Test of the previous stuff
        //=====
        //randomaux=unidistrib(rng);
        //cout << wsum1 << " wsum1 "<< endl;
        //cout << randomaux << " randomaux " << endl;
        
        //w1_rr_vector[0]=0.0;
        //w1_rr_vector[1]=0.0;
        //w1_rr_vector[2]=1.0;
        //S1Vector[0]=80;
        //S1Vector[1]=25;
        //S1Vector[2]=10;
        
        
        //Getting one draw:
        //for (int mm=0; mm<50; mm=mm+1){
        //    randomaux=unidistrib(rng);
        //    rr=0;
        
        //    while (condition==0){
        //        if (randomaux<w1_rr_vector[rr]){
        //            randomchosen=S1Vector[rr];
        //            condition=1;
        //        }
        //        randomaux-=w1_rr_vector[rr];
        //        rr=rr+1;
        //
        //    }
        //    cout << randomchosen << " random number " << endl;
        //    condition=0;
        //}
        
        
        
        //========================================================
        //2012 simulations
        //========================================================
        rr=0;
        while (rr<RR){
            logcontribobs_rr=0;
            loglikskills_rr=0;
            loglikeffm_rr=0;
            loglikefff_rr=0;
            loglikmmu_rr=0;
            loglikeinv_rr=0;
            Decision=0;
            logcontribdecision_rr=0;
            
            //1. Draw shocks from the corresponding distributions
            //===========================================================
            //1.1 Single mothers:
            if (Cliveswithfather12[ii]==0 & Cliveswithmother12[ii]==1){
                
                //========================================================
                //1. Effort levels
                //========================================================
                //Draw a shock for the effort levels
                Meffort_factor12shock=gen_normal_3(generator);
                //========================================================
                //1.1 If mom does'nt work
                //========================================================
                // Getting a draw from a standard normal distribution
                
                // Now transforming it as a draw from the given distribution
                Meffort_factor12NW=Meffort_factor12shock*stdeffortmom;
                Meffort_factor12NW=exp(Meffort_factor12NW)*MeffortMean12NW;
                //Effort can't be negative, we are drawing from a truncated normal distribution
                if (Meffort_factor12NW<0){
                    Meffort_factor12NW=0;
                }
                //========================================================
                //1.2 If mom works
                // Getting a draw from a standard normal distribution
                
                // Now transforming it as a draw from the given distribution
                Meffort_factor12W=Meffort_factor12shock*stdeffortmom;
                Meffort_factor12W=exp(Meffort_factor12W)*MeffortMean12W;
                
                //Effort can't be negative, we are drawing from a truncated normal distribution
                if (Meffort_factor12W<0){
                    Meffort_factor12W=0;
                }
                //========================================================
                
                
                //========================================================
                //2. Investment levels
                //========================================================
                //Getting the draw for the investment shock
                CInv_factor12shock=gen_normal_3(generator);
                //========================================================
                //2.1 If mom doesn't work
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CInv_factor12NW=CInv_factor12shock*stdinvestment;
                CInv_factor12NW=CInv_factor12NW+CInvestmentMean12NW;
                
                //Investment can't be negative
                if (CInv_factor12NW<0){
                    CInv_factor12NW=0;
                }
                //========================================================
                //2.2 If mom works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CInv_factor12W=CInv_factor12shock*stdinvestment;
                CInv_factor12W=CInv_factor12W+CInvestmentMean12W;
                //Investment can't be negative
                if (CInv_factor12W<0){
                    CInv_factor12W=0;
                }
                
                //==========================================================
                //3. Skills levels
                //==========================================================
                //Shock of skills
                CSkills_factor12shock=gen_normal_3(generator);
                //========================================================
                //1. If mom doesn't works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CSkillsMean12NW=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                             Cedad_meses12[ii],ttheta0,
                                             ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                             Feffort12[ii],Meffort_factor12NW,
                                             CInv_factor12NW,S1Vector[rr],
                                             Cchildcare12[ii],PGVector[rr],Hmemberstotal12[ii]);
                
                if (CSkills_factor12shock*stdskills<=500){
                    CSkills_factor12NW=CSkills_factor12shock*stdskills;
                }
                if (CSkills_factor12shock*stdskills>=500){
                    CSkills_factor12NW=500;
                }
                
                
                CSkills_factor12NW=exp(CSkills_factor12NW)*CSkillsMean12NW;
                //========================================================
                //2. If mom works
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CSkillsMean12W=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                            Cedad_meses12[ii],ttheta0,
                                            ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                            Feffort12[ii],Meffort_factor12W,
                                            CInv_factor12W,
                                            S1Vector[rr],
                                            Cchildcare12[ii],
                                            PGVector[rr],Hmemberstotal12[ii]);
                
                if (CSkills_factor12shock*stdskills<=500){
                    CSkills_factor12W=CSkills_factor12shock*stdskills;
                }
                if (CSkills_factor12shock*stdskills>=500){
                    CSkills_factor12W=500;
                }
                
                
                CSkills_factor12W=exp(CSkills_factor12W)*CSkillsMean12W;
                
                //===========================================================
                //===========================================================
                //4 Preference shocks
                //===========================================================
                //1. If mom doesn't work
                //===========================================================
                MprefshockNW=gen_normal_3(generator);
                MprefshockNW=MprefshockW*MshockNWA;
                //===========================================================
                //2. If mom works
                //===========================================================
                MprefshockW=gen_normal_3(generator);
                MprefshockW=MprefshockW*MshockWA;
                //===========================================================
                //===========================================================
                //5 Consumption
                //===========================================================
                //1. If mom doesn't work
                //===========================================================
                MconsumptionNW=Mnly12[ii]-priceINV*CInv_factor12NW;
                if (MconsumptionNW<0){
                    MconsumptionNW=1.0e-10;
                }
                //===========================================================
                //2. If mom works
                //===========================================================
                MconsumptionW=Mnly12[ii]-priceINV*CInv_factor12NW+Mwageinlik;
                if (MconsumptionW<0){
                    MconsumptionW=1.0e-10;
                }
                
                
                
                //===========================================================
                //===========================================================
                //6. Behavioral model
                //==========================================================
                //1. If mom doesn't work
                //==========================================================
                MutilityNW=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,MconsumptionNW,
                                     Meffort_factor12NW,0,CSkills_factor12NW)+MprefshockNW;
                
                //===========================================================
                //2. If mom works
                //===========================================================
                MutilityW=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,MconsumptionW,
                                    Meffort_factor12W,1,CSkills_factor12W)+MprefshockW;
                
                //===========================================================
                //===========================================================
                //Computing the actual likelihood
                
                Decision=(MutilityW>MutilityNW);
                //We want to call the decision=1 if the decission is correct. Then if person
                //decides not to work the thing will be 1-that.
                Decision=(obDecision==Decision);
                
                //Likelihood contributions of the measurement systems
                if (Mfraclabor12[ii]==1){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12W,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12W;
                    
                    //2. Effort
                    loglikeffm_rr=
                    F_loglikelihood_generic(Meffort12[ii],
                                            Meffort_factor12W,MEASMeffort);
                    loglikefff_rr=0;
                    
                    //3. Investment
                    loglikeinv_rr=
                    F_loglikelihood_generic(Cfactorinv12[ii],CInv_factor12W,
                                            MEASINV);
                    
                }
                if (Mfraclabor12[ii]==0){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12NW,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12W;
                    
                    
                    //2. Effort
                    loglikeffm_rr=
                    F_loglikelihood_generic(Meffort12[ii],
                                            Meffort_factor12NW,MEASMeffort);
                    loglikefff_rr=0;
                    
                    //3. Investment
                    loglikeinv_rr=
                    F_loglikelihood_generic(Cfactorinv12[ii],CInv_factor12W,
                                            MEASINV);
                    
                }
                
                
                
            }//End of single mothers
            //If both parents are present.
            
            if (Cliveswithfather12[ii]==1 & Cliveswithmother12[ii]==1){
                
                //======================
                //0. Mmu
                //======================
                //Getting the draw for the unobserved heterogeneity of mmy
                MMushock=gen_normal_3(generator);
                MMushock=MMushock*stdmmu;
                //Getting the predicted level of bargaining power
                mupred=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                             llambda5,llambda6,llambda7,llambda8,
                             Fwageinlik,Mwageinlik,
                             Fnly12[ii],Mnly12[ii],MMushock,mmuLB,mmuUB,
                             Fage12[ii],Mage12[ii],Fyrschool12[ii],Myrschool12[ii],
                             MRATIO[ii],Unemployment[ii],Wageratio[ii],Distance[ii]);
                
                
                
                //========================================================
                //1. Effort levels
                //========================================================
                
                //Draw a shock for the effort levels
                Meffort_factor12shock=exp(gen_normal_3(generator)*stdeffortmom);
                Feffort_factor12shock=exp(gen_normal_3(generator)*stdeffortfat);
                
                //Getting the means
                //1.1 Effort levels
                //1.1.1 NWNW
                MeffortMean12NWNW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                                 aalpha4f,aalpha4m,0,0,pphi,
                                                 ttheta2)*Meffort_factor12shock;
                
                FeffortMean12NWNW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                                 aalpha4f,aalpha4m,0,0,pphi,
                                                 ttheta2)*Feffort_factor12shock;
                
                //1.1.2 Father works mother not
                MeffortMean12WNW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                                aalpha4f,aalpha4m,1,0,pphi,
                                                ttheta2)*Meffort_factor12shock;
                
                FeffortMean12WNW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                                aalpha4f,aalpha4m,1,0,pphi,
                                                ttheta2)*Feffort_factor12shock;
                
                //1.1.3 Mother works father doesn't
                MeffortMean12NWW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                                aalpha4f,aalpha4m,0,1,pphi,
                                                ttheta2)*Meffort_factor12shock;
                
                FeffortMean12NWW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                                aalpha4f,aalpha4m,0,1,pphi,
                                                ttheta2)*Feffort_factor12shock;
                //1.1.4 Both work
                MeffortMean12WW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                               aalpha4f,aalpha4m,1,1,pphi,
                                               ttheta2)*Meffort_factor12shock;
                
                FeffortMean12WW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                               aalpha4f,aalpha4m,1,1,pphi,
                                               ttheta2)*Feffort_factor12shock;
                
                
                //========================================================
                //2. Investment levels
                //========================================================
                //Getting the draw for the investment shock
                CInv_factor12shock=gen_normal_3(generator);
                //========================================================
                //2.1 If mom doesn't work
                //========================================================
                //3. Now transforming it as a draw from the given distribution
                CInv_factor12shock=CInv_factor12shock*stdinvestment;
                //4. Now get the mean of the investment as a function of mmu
                //for the possible combinations
                
                CInvestmentMean12NWNW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                      aalpha2m,aalpha2f,
                                                      mupred,0,0,
                                                      Fwageinlik,Mwageinlik,
                                                      Fnly12[ii],Mnly12[ii],
                                                      ttheta1,priceINV)+CInv_factor12shock;
                
                if (CInvestmentMean12WW_TOG<=0){
                    CInvestmentMean12WW_TOG=1.0e-3;
                }
                CInvestmentMean12WNW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                     aalpha2m,aalpha2f,
                                                     mupred,1,0,
                                                     Fwageinlik,Mwageinlik,
                                                     Fnly12[ii],Mnly12[ii],
                                                     ttheta1,priceINV)+CInv_factor12shock;
                if (CInvestmentMean12WNW_TOG<=0){
                    CInvestmentMean12WNW_TOG=1.0e-3;
                }
                CInvestmentMean12NWW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                     aalpha2m,aalpha2f,
                                                     mupred,0,1,
                                                     Fwageinlik,Mwageinlik,
                                                     Fnly12[ii],Mnly12[ii],
                                                     ttheta1,priceINV)+CInv_factor12shock;
                if (CInvestmentMean12NWW_TOG<=0){
                    CInvestmentMean12NWW_TOG=1.0e-3;
                }
                CInvestmentMean12WW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                    aalpha2m,aalpha2f,
                                                    mupred,1,1,
                                                    Fwageinlik,Mwageinlik,
                                                    Fnly12[ii],Mnly12[ii],
                                                    ttheta1,priceINV)+CInv_factor12shock;
                if (CInvestmentMean12WW_TOG<=0){
                    CInvestmentMean12WW_TOG=1.0e-3;
                }
                
                //==========================================================
                //3. Skills levels
                //==========================================================
                //Shock of skills
                CSkills_factor12shock=gen_normal_3(generator);
                if (CSkills_factor12shock*stdskills<=500){
                    CSkills_factor12shock=exp(CSkills_factor12shock*stdskills);
                }
                if (CSkills_factor12shock*stdskills>500){
                    CSkills_factor12shock=exp(500);
                }
                
                
                //Now getting the different levels of skills for all possi-
                //bilities
                //3.1 NWNW
                CSkills_factor12NWNW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                      Cedad_meses12[ii],
                                                      ttheta0,ttheta1,ttheta2,pphi,
                                                      ggammaf,ggammam,FeffortMean12NWNW_TOG,
                                                      MeffortMean12NWNW_TOG,
                                                      CInvestmentMean12NWNW_TOG,
                                                      S1Vector[rr],
                                                      Cchildcare12[ii],
                                                      PGVector[rr],Hmemberstotal12[ii])*CSkills_factor12shock;
                
                //3.2 WNW
                CSkills_factor12WNW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                     Cedad_meses12[ii],
                                                     ttheta0,ttheta1,ttheta2,pphi,
                                                     ggammaf,ggammam,FeffortMean12WNW_TOG,
                                                     MeffortMean12WNW_TOG,
                                                     CInvestmentMean12WNW_TOG,
                                                     S1Vector[rr],
                                                     Cchildcare12[ii],
                                                     PGVector[rr],Hmemberstotal12[ii])*CSkills_factor12shock;
                
                //3.3 NWW
                CSkills_factor12NWW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                     Cedad_meses12[ii],
                                                     ttheta0,ttheta1,ttheta2,pphi,
                                                     ggammaf,ggammam,FeffortMean12NWW_TOG,
                                                     MeffortMean12NWW_TOG,
                                                     CInvestmentMean12NWW_TOG,
                                                     S1Vector[rr],Cchildcare12[ii],
                                                     PGVector[rr],Hmemberstotal12[ii])*CSkills_factor12shock;
                
                //3.4 WW
                CSkills_factor12WW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                    Cedad_meses12[ii],
                                                    ttheta0,ttheta1,ttheta2,pphi,
                                                    ggammaf,ggammam,FeffortMean12WW_TOG,
                                                    MeffortMean12WW_TOG,
                                                    CInvestmentMean12WW_TOG,
                                                    S1Vector[rr],Cchildcare12[ii],
                                                    PGVector[rr],Hmemberstotal12[ii])*CSkills_factor12shock;
                
                //===========================================================
                //===========================================================
                //4 Preference shocks
                
                //1. If mom doesn't work
                MprefshockNW=gen_normal_3(generator);
                MprefshockNW=MprefshockW*MshockNWA;
                //2. If mom works
                MprefshockW=gen_normal_3(generator);
                MprefshockW=MprefshockW*MshockWA;
                //3. If FATHER doesn't work
                FprefshockNW=gen_normal_3(generator);
                FprefshockNW=FprefshockW*FshockNWA;
                //4. If FATHER works
                FprefshockW=gen_normal_3(generator);
                FprefshockW=FprefshockW*FshockWA;
                
                
                
                //===========================================================
                //===========================================================
                //5 Consumption
                //===========================================================
                //1. Mother
                //===========================================================
                MconsumptionNWNW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                       CInvestmentMean12NWNW_TOG,
                                                       ttheta1,mupred,priceINV);
                
                
                MconsumptionWNW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                      CInvestmentMean12WNW_TOG,
                                                      ttheta1,mupred,priceINV);
                
                MconsumptionNWW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                      CInvestmentMean12NWW_TOG,
                                                      ttheta1,mupred,priceINV);
                
                MconsumptionWW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                     CInvestmentMean12WW_TOG,
                                                     ttheta1,mupred,priceINV);
                
                //===========================================================
                //1. Father
                //===========================================================
                FconsumptionNWNW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                       CInvestmentMean12NWNW_TOG,
                                                       ttheta1,mupred,priceINV);
                
                FconsumptionWNW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                      CInvestmentMean12WNW_TOG,
                                                      ttheta1,mupred,priceINV);
                
                FconsumptionNWW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                      CInvestmentMean12NWW_TOG,
                                                      ttheta1,mupred,priceINV);
                
                FconsumptionWW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                     CInvestmentMean12WW_TOG,
                                                     ttheta1,mupred,priceINV);
                
                
                
                
                //===========================================================
                //===========================================================
                //6. Behavioral model
                //==========================================================
                //1. Mother
                //==========================================================
                
                MutilityNWNW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                           MconsumptionNWNW_TOG,MeffortMean12NWNW_TOG,0,CSkills_factor12NWNW_TOG)+
                MprefshockNW;
                
                MutilityWNW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                          MconsumptionWNW_TOG,MeffortMean12WNW_TOG,0,CSkills_factor12WNW_TOG)+MprefshockNW;
                
                MutilityNWW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                          MconsumptionNWW_TOG,MeffortMean12NWW_TOG,1,CSkills_factor12NWW_TOG)+MprefshockW;
                
                MutilityWW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                         MconsumptionWW_TOG,MeffortMean12WW_TOG,1,CSkills_factor12WW_TOG)+
                MprefshockW;
                
                //==========================================================
                //2. Father
                //==========================================================
                
                FutilityNWNW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f
                                           ,FconsumptionNWNW_TOG,
                                           FeffortMean12NWNW_TOG,0,CSkills_factor12NWNW_TOG)+
                FprefshockNW;
                
                FutilityWNW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                          FconsumptionWNW_TOG,
                                          FeffortMean12WNW_TOG,1,CSkills_factor12WNW_TOG)+
                FprefshockW;
                
                FutilityNWW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                          FconsumptionNWW_TOG,
                                          FeffortMean12NWW_TOG,0,CSkills_factor12NWW_TOG)+
                FprefshockNW;
                
                FutilityWW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                         FconsumptionWW_TOG,
                                         FeffortMean12WW_TOG,1,CSkills_factor12WW_TOG)+
                FprefshockW;
                
                //Seeing the utility levels for each case:
                
                
                //==========================================================
                // Defining the utilities for the possible combinations
                //==========================================================
                
                WelfareNWNW=mupred*FutilityNWNW_TOG+(1-mupred)*MutilityNWNW_TOG;
                WelfareWNW=mupred*FutilityWNW_TOG+(1-mupred)*MutilityWNW_TOG;
                WelfareNWW=mupred*FutilityNWW_TOG+(1-mupred)*MutilityNWW_TOG;
                WelfareWW=mupred*FutilityWW_TOG+(1-mupred)*MutilityWW_TOG;
                
                //==========================================================
                //Identifying the correct decision
                //==========================================================
                if ( (WelfareNWNW>WelfareWNW) && ( WelfareNWNW>WelfareNWW) && ( WelfareNWNW>WelfareWW)) {
                    Decision=4;
                }
                else if ( (WelfareWNW>WelfareNWW) && ( WelfareNWNW>WelfareWW) ) {
                    Decision=2;
                }
                else if ( (WelfareNWW>WelfareWW)){
                    Decision=3;
                }
                else{
                    Decision=1;
                }
                
                
                Decision=(Decision==obDecision);
                //7. Likelihood of measurement system
                //===============================
                
                //7.1 NWNW
                if (Mfraclabor12[ii]==0 & Ffraclabor12[ii]==0){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12NWNW_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12NWNW_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr=F_loglikelihood_generic(Meffort12[ii],
                                                          MeffortMean12NWNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr=F_loglikelihood_generic(Feffort12[ii],
                                                          FeffortMean12NWNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr=F_loglikelihood_generic(Cfactorinv12[ii],CInvestmentMean12NWNW_TOG,MEASINV);
                }
                
                
                //7.2 WNW
                if (Mfraclabor12[ii]==0 & Ffraclabor12[ii]==1){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12WNW_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12WNW_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr=F_loglikelihood_generic(Meffort12[ii],
                                                          MeffortMean12WNW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr=F_loglikelihood_generic(Feffort12[ii],
                                                          FeffortMean12WNW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr=F_loglikelihood_generic(Cfactorinv12[ii],CInvestmentMean12WNW_TOG,MEASINV);
                }
                //7.3 NWW
                if (Mfraclabor12[ii]==1 & Ffraclabor12[ii]==0){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12NWW_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12NWW_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr=F_loglikelihood_generic(Meffort12[ii],
                                                          MeffortMean12NWW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr=F_loglikelihood_generic(Feffort12[ii],
                                                          FeffortMean12NWW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    loglikeinv_rr=F_loglikelihood_generic(Cfactorinv12[ii],CInvestmentMean12NWW_TOG,MEASINV);
                }
                //7.4 WW
                if (Mfraclabor12[ii]==1 & Ffraclabor12[ii]==1){
                    
                    //1. Skills
                    loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012[ii],
                                                            CSkills_factor12WW_TOG,MEASSkills);
                    
                    //Storing the particle weight and the particle itself
                    
                    what2_rr_vector[rr]=loglikskills_rr;
                    S2Vector[rr]=CSkills_factor12WW_TOG;
                    
                    
                    //2. Mother's effort
                    loglikeffm_rr=F_loglikelihood_generic(Meffort12[ii],
                                                          MeffortMean12WW_TOG,MEASMeffort);
                    
                    //3. Father's effort
                    loglikefff_rr=F_loglikelihood_generic(Feffort12[ii],
                                                          FeffortMean12WW_TOG,MEASFeffort);
                    
                    //4. Investment level
                    
                    loglikeinv_rr=F_loglikelihood_generic(Cfactorinv12[ii],CInvestmentMean12WW_TOG,MEASINV);
                }
                
                //5. Measurement system for mmu, doesn't depend on labor supply
                loglikmmu_rr=F_loglikelihood_generic(Hbarg[ii],mupred,MEASMMu);
                
            }//If living together
            
            //Adding up likelihood contributions of simulations
            
            
            //0 Decisions (behavioral system)
            logcontribdecision+=Decision; //Storing the number of times correct thing
            //1. Skills
            
            loglikskills+=loglikskills_rr;
            
            
            
            //2. Mother's effort
            loglikeffm+=loglikeffm_rr;
            
            //3. Father's effort
            loglikefff+=loglikefff_rr;
            
            //4. Investment
            loglikeinv+=loglikeinv_rr;
            
            //5. Mmu
            loglikmmu+=loglikmmu_rr;
            
            rr=rr+1;
        } //Finish the RR simulations
        //At the end of the R computations we will perform the fraction of times person
        //decided to work.
        
        //0. Behavioral decisions
        
        logcontribdecision=logcontribdecision/RR;
        logcontribdecision10=logcontribdecision10/RR;
        
        //We need to fix in those cases where we didn't catch a single time. The likelihood
        //will then be zero.
        
        
        //==================
        //Decisions
        //==================
        //2010
        if (logcontribdecision10==0){
            logcontribdecision10=1.0e-20;
        }
        logcontribdecision10=log(logcontribdecision10);
        //2012
        if (logcontribdecision==0){
            logcontribdecision=1.0e-20;
        }
        logcontribdecision=log(logcontribdecision);
        
        //=============
        //Skills
        //=============
        //1. Skills 2010
        loglikskills10=loglikskills10/RR;
        if (loglikskills10==0){
            loglikskills10=1.0e-20;
        }
        loglikskills10=log(loglikskills10);
        
        //1. Skills 2012
        loglikskills=loglikskills/RR;
        if (loglikskills==0){
            loglikskills=1.0e-20;
        }
        loglikskills=log(loglikskills);
        
        //==================
        //2. Mother's effort
        //==================
        //2010
        if (Cliveswithmother10[ii]==1){
            loglikeffm10=loglikeffm10/RR;
            if (loglikeffm10==0){
                loglikeffm10=1.0e-20;
            }
            loglikeffm10=log(loglikeffm10);
        }
        
        //2012
        if (Cliveswithmother12[ii]==1){
            loglikeffm=loglikeffm/RR;
            if (loglikeffm==0){
                loglikeffm=1.0e-20;
            }
            loglikeffm=log(loglikeffm);
        }
        
        
        
        //3. Father's effort
        
        //2010
        if (Cliveswithfather10[ii]==1){
            loglikefff10=loglikefff10/RR;
            if (loglikefff10==0){
                loglikefff10=1.0e-20;
            }
            loglikefff10=log(loglikefff10);
        }
        
        //2012
        if (Cliveswithfather12[ii]==1){
            loglikefff=loglikefff/RR;
            if (loglikefff==0){
                loglikefff=1.0e-20;
            }
            loglikefff=log(loglikefff);
        }
        //==================
        //4. Investment
        //=================
        
        //2010
        loglikeinv10=loglikeinv10/RR;
        if (loglikeinv10==0){
            loglikeinv10=1.0e-20;
        }
        loglikeinv10=log(loglikeinv10);
        
        //2012
        loglikeinv=loglikeinv/RR;
        if (loglikeinv==0){
            loglikeinv=1.0e-20;
        }
        loglikeinv=log(loglikeinv);
        
        //5. Mmu
        //2012
        if (Cliveswithfather12[ii]==1){
            loglikmmu=loglikmmu/RR;
            if (loglikmmu==0){
                loglikmmu=1.0e-20;
            }
            loglikmmu=log(loglikmmu);
        }
        
        //Skills at birth
        loglik_S0=loglik_PGrr/RR;
        if (loglik_S0==0){
            loglik_S0=1.0e-20;
        }
        loglik_S0=log(loglik_S0);
        
        //Primary caregiver skills
        loglik_PG=loglik_PGrr/RR;
        if (loglik_PG==0){
            loglik_PG=1.0e-20;
        }
        loglik_PG=log(loglik_PG);
        
        
        
        logcontribobs_ii+=likwageF10+likwageM10+likwageF+likwageM+
        loglikefff+loglikefff10+
        loglikeffm10+loglikeffm+
        logcontribdecision+logcontribdecision10+
        loglikeinv+loglikeinv10+
        loglikskills+loglikskills10+
        loglikmmu;
        
        
        
        
        
    } //Finish the ii computations
    //To the
    loglik=logcontribobs_ii;
    loglik=-loglik;
    return(loglik);
}//End of function



//====================
//Trying to define likelihood function exclusively as function of parameters
long double F_likelihood_FIN(const double *PARDOUBLE){
    //0. Define the current directory
    //chdir("/Users/rodrigoazuero/Documents/Research/Chile/RR/BEHAVIORAL38/BEHAVIORAL38/BEHAVIORAL38");
    
    int NVAR=63; //Number of variables in file
    int SIZEOBS=1278; //Number of observations in file
    
    //First step making the vector<double> from the const(double *par)
    vector<double> PAR;
    PAR.resize(78);
    
    for (int it=0; it<77; it=it+1){
        PAR[it]=PARDOUBLE[it];
    }
    
    //========================================
    //BLOCK OF GETTING THE DATA INTO VECTORS
    //=========================================
    std::ifstream theFile ("BEHAVIORAL38excel1.csv");
    double MYARRAY[SIZEOBS+1][NVAR+1];
    // ...
    
    std::string line;
    std::vector<std::vector<std::string> > values;
    int it=0;
    int it2=0;
    std::string line_value;
    std::vector<std::string> line_values;
    std::stringstream ss;
    while(std::getline(theFile, line))
    {
        
        ss<<line;
        //std::stringstream ss(line);
        //std::string item;
        //cout <<  << "linevalprev"<<endl;
        while(std::getline(ss, line_value, ','))
        {
            line_values.push_back(line_value);
            MYARRAY[it][it2] = ::atof(line_value.c_str());
            //MYARRAY[it][it2]=std::stod (line_value); only for c++11 compi
            it2=it2+1;
            if (it2==NVAR){ //later change 4 for
                it2=0;
            }
        }
        values.push_back(line_values);
        
        //For c++11 used values.emplace_back(line_values);
        
        //cout << line_value<< "line_value2"<< endl;
        it=it+1;
        //Free the string types
        line_value.clear();
        line_values.clear();
        ss.clear();
        
    }
    
    //Defining the vectors
    vector<double> Hmemberstotal12;
    Hmemberstotal12.resize(SIZEOBS);
    
    vector<double> Hmemberstotal10;
    Hmemberstotal10.resize(SIZEOBS);
    
    vector<double> Cedad_meses10;
    Cedad_meses10.resize(SIZEOBS);
    
    vector<double> Ctestsfactor1_10;
    Ctestsfactor1_10.resize(SIZEOBS);
    
    vector<double> Mfraclabor10;
    Mfraclabor10.resize(SIZEOBS);
    
    vector<double> Mwage10;
    Mwage10.resize(SIZEOBS);
    
    vector<double> Mnlincome10;
    Mnlincome10.resize(SIZEOBS);
    
    vector<double> Ffraclabor10;
    Ffraclabor10.resize(SIZEOBS);
    
    vector<double> Fwage10;
    Fwage10.resize(SIZEOBS);
    
    vector<double> Fnlincome10;
    Fnlincome10.resize(SIZEOBS);
    
    vector<double> Cfactorbirth;
    Cfactorbirth.resize(SIZEOBS);
    
    vector<double> Meffort10;
    Meffort10.resize(SIZEOBS);
    
    vector<double> Feffort10;
    Feffort10.resize(SIZEOBS);
    
    vector<double> CfactorInv10;
    CfactorInv10.resize(SIZEOBS);
    
    vector<double> Cchildcare10;
    Cchildcare10.resize(SIZEOBS);
    
    vector<double> Cliveswithfather10;
    Cliveswithfather10.resize(SIZEOBS);
    
    vector<double> Cliveswithmother10;
    Cliveswithmother10.resize(SIZEOBS);
    
    vector<double> Cedad_meses12;
    Cedad_meses12.resize(SIZEOBS);
    
    vector<double> Ctestsfactor4_2012;
    Ctestsfactor4_2012.resize(SIZEOBS);
    
    vector<double> Cliveswithfather12;
    Cliveswithfather12.resize(SIZEOBS);
    
    vector<double> Cliveswithmother12;
    Cliveswithmother12.resize(SIZEOBS);
    
    vector<double> Mage12;
    Mage12.resize(SIZEOBS);
    
    vector<double> Myrschool12;
    Myrschool12.resize(SIZEOBS);
    
    vector<double> Mfraclabor12;
    Mfraclabor12.resize(SIZEOBS);
    
    vector<double> Fage12;
    Fage12.resize(SIZEOBS);
    
    vector<double> Fyrschool12;
    Fyrschool12.resize(SIZEOBS);
    
    vector<double> Ffraclabor12;
    Ffraclabor12.resize(SIZEOBS);
    
    vector<double> Mwage12;
    Mwage12.resize(SIZEOBS);
    
    vector<double> Mnlincome12;
    Mnlincome12.resize(SIZEOBS);
    
    vector<double> Fwage12;
    Fwage12.resize(SIZEOBS);
    
    vector<double> Fnlincome12;
    Fnlincome12.resize(SIZEOBS);
    
    vector<double> Hbarg;
    Hbarg.resize(SIZEOBS);
    
    vector<double> Hbarg1;
    Hbarg1.resize(SIZEOBS);
    
    vector<double> Hbarg2;
    Hbarg2.resize(SIZEOBS);
    
    vector<double> Hbarg3;
    Hbarg3.resize(SIZEOBS);
    
    vector<double> Hbarg4;
    Hbarg4.resize(SIZEOBS);
    
    vector<double> Meffort12;
    Meffort12.resize(SIZEOBS);
    
    vector<double> Feffort12;
    Feffort12.resize(SIZEOBS);
    
    vector<double> CfactorInv12;
    CfactorInv12.resize(SIZEOBS);
    
    vector<double> Hchores12;
    Hchores12.resize(SIZEOBS);
    
    vector<double> Cchildcare12;
    Cchildcare12.resize(SIZEOBS);
    
    vector<double> Ccareskills;
    Ccareskills.resize(SIZEOBS);
    
    vector<double> Hchildcareobs;
    Hchildcareobs.resize(SIZEOBS);
    
    vector<double> FMRATIO;
    FMRATIO.resize(SIZEOBS);
    
    vector<double> Unemployment;
    Unemployment.resize(SIZEOBS);
    
    vector<double> Wageratio;
    Wageratio.resize(SIZEOBS);
    
    vector<double> Distance;
    Distance.resize(SIZEOBS);
    
    vector<double> MTJH;
    MTJH.resize(SIZEOBS);
    
    vector<double> Magegroup10;
    Magegroup10.resize(SIZEOBS);
    
    vector<double> Magegroup12;
    Magegroup12.resize(SIZEOBS);
    
    vector<double> Fage10;
    Fage10.resize(SIZEOBS);
    
    vector<double> Mage10;
    Mage10.resize(SIZEOBS);
    
    vector<double> Mage2_12;
    Mage2_12.resize(SIZEOBS);
    
    vector<double> Fage2_12;
    Fage2_12.resize(SIZEOBS);
    
    vector<double> Mage2_10;
    Mage2_10.resize(SIZEOBS);
    
    vector<double> Fage2_10;
    Fage2_10.resize(SIZEOBS);
    
    vector<double> Mwagepred;
    Mwagepred.resize(SIZEOBS);
    
    vector<double> Fwagepred;
    Fwagepred.resize(SIZEOBS);
    
    vector<double> Agedif;
    Agedif.resize(SIZEOBS);
    
    vector<double> Edudif;
    Edudif.resize(SIZEOBS);
    
    vector<double> ymratio;
    ymratio.resize(SIZEOBS);
    
    vector<double> ymdif;
    ymdif.resize(SIZEOBS);
    
    vector<double> wageratio;
    wageratio.resize(SIZEOBS);
    
    for (int it=0; it<SIZEOBS;it++){
        Hmemberstotal12[it]=MYARRAY[it][0];
        Hmemberstotal10[it]=MYARRAY[it][1];
        Cedad_meses10[it]=MYARRAY[it][2];
        Ctestsfactor1_10[it]=exp(MYARRAY[it][3]);
        Mfraclabor10[it]=MYARRAY[it][4];
        Mwage10[it]=MYARRAY[it][5];
        Mnlincome10[it]=MYARRAY[it][6];
        Ffraclabor10[it]=MYARRAY[it][7];
        Fwage10[it]=MYARRAY[it][8];
        Fnlincome10[it]=MYARRAY[it][9];
        Cfactorbirth[it]=exp(MYARRAY[it][10]);
        Meffort10[it]=MYARRAY[it][11];
        Feffort10[it]=MYARRAY[it][12];
        CfactorInv10[it]=exp(MYARRAY[it][13]);
        Cchildcare10[it]=MYARRAY[it][14];
        Cliveswithfather10[it]=MYARRAY[it][15];
        Cliveswithmother10[it]=MYARRAY[it][16];
        Cedad_meses12[it]=MYARRAY[it][17];
        Ctestsfactor4_2012[it]=exp(MYARRAY[it][18]);
        Cliveswithfather12[it]=MYARRAY[it][19];
        Cliveswithmother12[it]=MYARRAY[it][20];
        Mage12[it]=MYARRAY[it][21];
        Myrschool12[it]=MYARRAY[it][22];
        Mfraclabor12[it]=MYARRAY[it][23];
        Fage12[it]=MYARRAY[it][24];
        Fyrschool12[it]=MYARRAY[it][25];
        Ffraclabor12[it]=MYARRAY[it][26];
        Mwage12[it]=MYARRAY[it][27];
        Mnlincome12[it]=MYARRAY[it][28];
        Fwage12[it]=MYARRAY[it][29];
        Fnlincome12[it]=MYARRAY[it][30];
        Hbarg[it]=MYARRAY[it][31];
        Hbarg1[it]=MYARRAY[it][32];
        Hbarg2[it]=MYARRAY[it][33];
        Hbarg3[it]=MYARRAY[it][34];
        Hbarg4[it]=MYARRAY[it][35];
        Meffort12[it]=MYARRAY[it][36];
        Feffort12[it]=MYARRAY[it][37];
        CfactorInv12[it]=exp(MYARRAY[it][38]);
        Hchores12[it]=MYARRAY[it][39];
        Cchildcare12[it]=MYARRAY[it][40];
        Ccareskills[it]=MYARRAY[it][41];
        Hchildcareobs[it]=MYARRAY[it][42];
        FMRATIO[it]=MYARRAY[it][43];
        Unemployment[it]=MYARRAY[it][44];
        Wageratio[it]=MYARRAY[it][45];
        Distance[it]=MYARRAY[it][46];
        MTJH[it]=MYARRAY[it][47];
        Magegroup10[it]=MYARRAY[it][48];
        Magegroup12[it]=MYARRAY[it][49];
        Fage10[it]=MYARRAY[it][50];
        Mage10[it]=MYARRAY[it][51];
        Mage2_12[it]=MYARRAY[it][52];
        Fage2_12[it]=MYARRAY[it][53];
        Mage2_10[it]=MYARRAY[it][54];
        Fage2_10[it]=MYARRAY[it][55];
        Mwagepred[it]=MYARRAY[it][56];
        Fwagepred[it]=MYARRAY[it][57];
        Agedif[it]=MYARRAY[it][58];
        Edudif[it]=MYARRAY[it][59];
        ymratio[it]=MYARRAY[it][60];
        ymdif[it]=MYARRAY[it][61];
        wageratio[it]=MYARRAY[it][62];
    }//Finish loading data
    
    
    //===============================================
    //Block to get the initial parameters to optimize
    //===============================================
    
    long double test2=F_likelihood(PAR, Cliveswithfather12, Cliveswithmother12, Hchores12, Meffort12, Feffort12, Mage12, Fage12, Cedad_meses12, Ctestsfactor4_2012, Ctestsfactor1_10, Myrschool12, Fyrschool12, Ffraclabor12, Mfraclabor12, Mwage12, Fwage12, Mnlincome12, Fnlincome12, CfactorInv12, Cchildcare12, Ccareskills, Hbarg4, Feffort10, Meffort10, CfactorInv10, Cfactorbirth, Mwage10, Fwage10, Mnlincome10, Fnlincome10, Cedad_meses10, Cliveswithfather10, Cliveswithmother10, Mage10, Fage10, Mfraclabor10, Ffraclabor10, Cchildcare10, Hchildcareobs, Ccareskills, MTJH,FMRATIO,Unemployment,
                                   Wageratio,Distance,Magegroup10,Magegroup12,Hmemberstotal10,
                                   Hmemberstotal12,1278);
    
    return(test2);
}

//Defining the likelihood function that is going to be minimized finally
int iterat=0;
double FLIKMINIMIZED(unsigned n, const double *x, double *grad, void *FLIKMNIMZED_data){
    double result=F_likelihood_FIN(x);
    cout << result << " evalF"<< endl;
    iterat=iterat+1;
    printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(result);
}

//Definitions of functions for the example of nlopt
// http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial

int main(int argc, const char * argv[])
{
    //0. Define the current directory
    //chdir("/Users/rodrigoazuero/Documents/Research/Chile/RR/BEHAVIORAL38/BEHAVIORAL38/BEHAVIORAL38");
    
    //===============================================
    //Block to get the initial parameters to optimize
    //===============================================
    
    //I will generate two vectors. One of type vector and
    //another one double because it seems that nlopt only allows
    //me to use double types objects for the minimization.
    std::string line;
    std::vector<std::vector<std::string> > values;
    vector<double> PAR;
    PAR.resize(78);
    
    double PARDOUBLE[78]={};
    
    
    std::ifstream PARAMETERS ("OPTIMIZEDPARAMETERS1.csv");
    
    std::string linePARAM;
    int itpar=0;
    while(std::getline(PARAMETERS, line))
    {
        PAR[itpar]=::atof (line.c_str());
        //PAR[itpar]=std::stod (line); only C++11s
        PARDOUBLE[itpar]=PAR[itpar];
        itpar=itpar+1;
    }
    //And the ones that were not stored
    PAR[69]=0;
    PAR[70]=0;
    PAR[71]=0;
    PAR[72]=0;
    PAR[73]=-10;
    PAR[74]=-10;
    PAR[75]=-10;
    PAR[76]=-10;
    PAR[77]=0;
    
    PARDOUBLE[69]=0;
    PARDOUBLE[70]=0;
    PARDOUBLE[71]=0;
    PARDOUBLE[72]=0;
    PARDOUBLE[73]=-10;
    PARDOUBLE[74]=-10;
    PARDOUBLE[75]=-10;
    PARDOUBLE[76]=-10;
    PARDOUBLE[77]=0;
    
    //And trying with the other function
    cout << F_likelihood_FIN(PARDOUBLE) << "final attempt " << endl;
    std::cout << "Hello, World!\n";
    //Now doing for the likelihood function
    cout <<"---------------"<< endl;
    nlopt_opt opt3;
    opt3 = nlopt_create(NLOPT_LN_COBYLA, 78); /* algorithm and dimensionality */
    nlopt_set_min_objective(opt3, FLIKMINIMIZED, NULL);
    nlopt_set_xtol_rel(opt3, 1.0e-4);
    double x3[78]={0};
    for(int it=0;it<78;it=it+1){
        x3[it]=PARDOUBLE[it];
    }
    
    
    
    
    double minf3; /* the minimum objective value, upon return */
    nlopt_optimize(opt3, x3, &minf3);
    cout << minf3 << " minf3 " << endl;
    cout << nlopt_optimize(opt3, x3, &minf3) << " code " << endl;
    if (nlopt_optimize(opt3, x3, &minf3) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
    }
    
    //Saving the vector of optimal parameters in csv file
    ofstream optparam("PAROPTFOUND.csv");
    for (int it=0;it<=78;it=it+1){
        optparam<<x3[it] << endl;
    }
    optparam.close();
    
    
    nlopt_destroy(opt3);
    
    return 0;
}
