//
//  main.cpp
//  2
//
//  Created by Rodrigo Azuero on 11/13/15.
//  Copyright (c) 2015 Rodrigo Azuero Melo. All rights reserved.
//

#include <iostream>
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
#include <algorithm>
using std::vector;
using namespace std;




vector<vector<vector<double> > > valueo(vector<double> ECON, vector<double> PAR,
                                        vector<double> Bgrid, vector<double> Agrid, vector<double> Thgrid, double tolerance){
    
    
    
    //0. Parameter definition
    //0.0 Economy
    
    double wh=ECON[0];
    double wl=ECON[1];
    double Ph=ECON[2];
    double r=ECON[3];
    
    //0.1 Utility function
    double ppsi=PAR[0];
    double eeta=PAR[1];
    double bbeta=PAR[2];
    double ddelta=PAR[3];
    
    //0.2 Grids
    long double bbmin=Bgrid[0];
    long double bbmax=Bgrid[1];
    int bbnumber=Bgrid[2];
    long double bbstep=(bbmax-bbmin)/(bbnumber-1);
    
    
    long double aamin=Agrid[0];
    long double aamax=Agrid[1];
    int aanumber=Agrid[2];
    long double aastep=(aamax-aamin)/(aanumber-1);
    
    
    
    long double tthetamin=Thgrid[0];
    long double tthetamax=Thgrid[1];
    int tthetanumber=Thgrid[2];
    long double tthetastep=(tthetamax-tthetamin)/(tthetanumber-1);
    
    //0.3 Initializing answer
    vector<vector<vector<double> > > ans;
    ans.resize(2);
    for(int k=0; k<2; k++){
        ans[k].resize(aanumber);
        for(int j=0; j<aanumber; j++){
            ans[k][j].resize(tthetanumber);
            for(int l=0; l<tthetanumber; l++){
                ans[k][j][l] = 0;
            }
        }
    }
    //VV1 initializtaion
    vector<vector<vector<vector<double> > > > vv1;
    vv1.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        vv1[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            vv1[n][k].resize(aanumber);
            for(int j=0; j<aanumber; j++){
                vv1[n][k][j].resize(bbnumber);
                for(int l=0; l<bbnumber; l++){
                    vv1[n][k][j][l] = 0;
                }
            }
        }
    }
    
    //VV0 initializtaion
    vector<vector<vector<vector<double> > > > vv0;
    vv0.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        vv0[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            vv0[n][k].resize(aanumber);
            for(int j=0; j<aanumber; j++){
                vv0[n][k][j].resize(bbnumber);
                for(int l=0; l<bbnumber; l++){
                    vv0[n][k][j][l] = 0;
                }
            }
        }
    }
    
    
    //VVhm initialization
    vector<vector<double> > VV1NEW;
    VV1NEW.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        VV1NEW[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            VV1NEW[n][k] = 0;
        }
    }
    
    //VV1old
    vector<vector<double> > VV0NEW;
    VV0NEW.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        VV0NEW[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            VV0NEW[n][k] = 0;
        }
    }
    
    //VVhm initialization
    vector<vector<double> > VV1OLD;
    VV1OLD.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        VV1OLD[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            VV1OLD[n][k] = 0;
        }
    }
    
    //VV1old
    vector<vector<double> > VV0OLD;
    VV0OLD.resize(aanumber);
    for (int n=0; n<aanumber; n++){
        VV0OLD[n].resize(tthetanumber);
        for(int k=0; k<tthetanumber; k++){
            VV0OLD[n][k] = 0;
        }
    }


    
    
    //1. Start value function iteration
    double error=10;
    
    

    

    
    
    int bbcoordinate=0;
    int aacoordinate=0;
    int aaSTATEcoordinate=0;
    int ttcoordinate=0;
    int step=1;
    double u1=0; //If education
    double u0=0; //If no education
    double incomeold1=0; //Income when old IF STUDY
    double incomeold0=0; //Income when old IF NO STUDY
    double incomenext0=0; //Income next generation if no education
    double incomenext1=0; //Income next generation if education

    double v0candidate=-10000;
    double v1candidate=-10000;
    double temp=0;
    double tempt2=0;
    
    //vv1[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=0;
    //vv0[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=0;
    
    
    //while (error>tolerance){
    for (double total=0; total<=3; total=total+1){

        error=0;
        aaSTATEcoordinate=0;
        for (long double aaSTATE=aamin; aaSTATE<=aamax; aaSTATE=aaSTATE+aastep){
            ttcoordinate=0;
            for (long double tt=tthetamin; tt<=tthetamax; tt=tt+tthetastep){
                aacoordinate=0;
                for (long double aa=aamin; aa<=aamax; aa=aa+aastep){
                    bbcoordinate=0;
                    for (long double bb=bbmin; bb<=bbmax; bb=bb+bbstep){
                        incomeold1=wh+aaSTATE*(1+r)-bb;
                        incomeold0=wl+aaSTATE*(1+r)-bb;
                        incomenext0=wl+bb-aa;
                        incomenext1=bb-aa-Ph;
                        u0=-1000;
                        u1=-1000;

                        //First if father studied
                        //===================================
                        //Utility if h==1
                        if (incomeold1>0 && incomenext1>0){
                            u1=pow(incomeold1,1-ppsi)/(1-ppsi)+
                            ddelta*pow(incomenext1,1-ppsi)/(1-ppsi)+
                            bbeta*ddelta*
                            (pow(tt,eeta)*VV1OLD[aacoordinate][ttcoordinate] +
                             (1-pow(tt,eeta))*VV0OLD[aacoordinate][ttcoordinate]);
                            
                        }
                        if (incomeold1<=0 || incomenext1<=0){
                            u1=-1.0e+15;
                        }
                        
                        //Utility for hh==0
                        if (incomeold1>0 && incomenext0>0){
                            u0=pow(incomeold1,1-ppsi)/(1-ppsi)+
                            ddelta*pow(incomenext0,1-ppsi)/(1-ppsi)+
                            bbeta*ddelta*VV0OLD[aacoordinate][ttcoordinate];
                        }
                        if (incomeold1<=0 || incomenext0<=0){
                            u0=-1.0e+15;
                        }
                        //Identify if person wants to study or not
                        vv1[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=u0;
                        if (u1>=u0){
                            vv1[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=u1;
                        }

                        //===================================
                        //If father did not studie
                        //===================================
                        //Utility if h==1
                        u0=-1000;
                        u1=-1000;
                        if (incomeold0>0 && incomenext1>0){
                            u1=pow(incomeold0,1-ppsi)/(1-ppsi)+
                            ddelta*pow(incomenext1,1-ppsi)/(1-ppsi)+
                            bbeta*ddelta*
                            (pow(tt,eeta)*VV1OLD[aacoordinate][ttcoordinate] +
                             (1-pow(tt,eeta))*VV0OLD[aacoordinate][ttcoordinate]);
                        }
                        if (incomeold0<=0 || incomenext1<=0){
                            u1=-1.0e+15;
                        }
                        
                        //Utility for hh==0
                        if (incomeold0>0 && incomenext0>0){
                            u0=pow(incomeold0,1-ppsi)/(1-ppsi)+
                            ddelta*pow(incomenext0,1-ppsi)/(1-ppsi)+
                            bbeta*ddelta*VV0OLD[aacoordinate][ttcoordinate];
                            
                        }
                        if (incomeold0<0 || incomenext0<=0){
                            u0=-1.0e+15;
                            
                        }
                        //Identify if person wants to study or not
                        
                        vv0[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=u0;
                        if (u1>=u0){
                            vv0[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate]=u1;
                        }
                        bbcoordinate+=step;
                    }//end bb
                    aacoordinate+=step;
                    
                    
                    
                }//end aa
                //Now we cand find the maximum. Initializing it at coordinate 00
                
                VV1NEW[aaSTATEcoordinate][ttcoordinate]=
                -1000;
                
                VV0NEW[aaSTATEcoordinate][ttcoordinate]=
                -10000;
                
                aacoordinate=0;
                
                for (long double aaNEXT=aamin; aaNEXT<=aamax; aaNEXT=aaNEXT+aastep){
                    bbcoordinate=0;
                    for (long double bbNEXT=bbmin; bbNEXT<=bbmax; bbNEXT=bbNEXT+bbstep){
                        
                        v1candidate=vv1[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate];
                        if (v1candidate>VV1NEW[aaSTATEcoordinate][ttcoordinate]){
                            tempt2=v1candidate;
                            VV1NEW[aaSTATEcoordinate][ttcoordinate]=v1candidate;
                        }
                        
                        
                        v0candidate=vv0[aaSTATEcoordinate][ttcoordinate][aacoordinate][bbcoordinate];
                        if (v0candidate>VV0NEW[aaSTATEcoordinate][ttcoordinate]){
                            VV0NEW[aaSTATEcoordinate][ttcoordinate]=v0candidate;
                        }
                        
                        
                        bbcoordinate+=step;
                    }//end bb
                    aacoordinate+=step;
                }//end aa
                
                //We have computed the maximum:
                
                error+=
                pow(VV1OLD[aaSTATEcoordinate][ttcoordinate]-
                    VV1NEW[aaSTATEcoordinate][ttcoordinate],2)+
                pow(VV0OLD[aaSTATEcoordinate][ttcoordinate]-
                    VV0NEW[aaSTATEcoordinate][ttcoordinate],2);
                //if (aaSTATEcoordinate==0 && ttcoordinate==0){
                    //cout << " ==a=== " << endl;
                    //cout << VV0NEW[0][0] << " new " << endl;
                    //cout << VV0OLD[0][0] << " old " << endl;
                    //cout << " ==b=== " << endl;
                //}
                tempt2=VV1OLD[aaSTATEcoordinate][ttcoordinate];
                temp=VV1NEW[aaSTATEcoordinate][ttcoordinate];
                VV1OLD[aaSTATEcoordinate][ttcoordinate]=temp;
                
                tempt2=VV0OLD[aaSTATEcoordinate][ttcoordinate];
                temp=VV0NEW[aaSTATEcoordinate][ttcoordinate];
                VV0OLD[aaSTATEcoordinate][ttcoordinate]=temp;

               
                ttcoordinate+=step;
            }//end ttheta
            aaSTATEcoordinate+=step;
            
        }//end aaSTATEcoordinate
        
    }//end while
    aacoordinate=0;
    for (long double aaFINAL=aamin; aaFINAL<=aamax; aaFINAL=aaFINAL+aastep){
        ttcoordinate=0;
        for (long double ttFINAL=tthetamin; ttFINAL<=tthetamax; ttFINAL=ttFINAL+tthetastep){
            temp=VV1NEW[aacoordinate][ttcoordinate];
            temp=VV0NEW[aacoordinate][ttcoordinate];
            ans[0][aacoordinate][ttcoordinate]=VV1NEW[aacoordinate][ttcoordinate];
            ans[1][aacoordinate][ttcoordinate]=VV0NEW[aacoordinate][ttcoordinate];
            ttcoordinate+=1;
        }
        aacoordinate+=1;
    }
    
    
    return ans;
    
    
}



int main(int argc, const char * argv[])
{
    //Parameter definitions
    double ppsi=2;
    double eeta=0.5;
    double bbeta=0.971;
    double ddelta=0.1;
    vector<double> PAR;
    PAR.resize(4);
    PAR[0]=ppsi;
    PAR[1]=eeta;
    PAR[2]=bbeta;
    PAR[3]=ddelta;
    
    //Economy parameters
    double wh=7;
    double wl=1;
    double Ph=0.3;
    double r=0.02;
    vector<double> ECON;
    ECON.resize(4);
    ECON[0]=wh;
    ECON[1]=wl;
    ECON[2]=Ph;
    ECON[3]=r;
    
    //Grid parameters
    long double bbmin=0;
    long double bbmax=1;
    long double bbnumber=20;
    vector<double>Bgrid;
    Bgrid.resize(3);
    Bgrid[0]=bbmin;
    Bgrid[1]=bbmax;
    Bgrid[2]=bbnumber;
    
    double aamin=0;
    double aamax=1;
    double aanumber=20;
    vector<double>Agrid;
    Agrid.resize(3);
    Agrid[0]=aamin;
    Agrid[1]=aamax;
    Agrid[2]=aanumber;
    
    double tthetamin=0.3;
    double tthetamax=0.5;
    double tthetanumber=20;
    vector<double>Thgrid;
    Thgrid.resize(3);
    Thgrid[0]=tthetamin;
    Thgrid[1]=tthetamax;
    Thgrid[2]=tthetanumber;
    
    double tolerance=1.0e-6;
    vector<vector<vector<double> > > Vy;
    Vy=valueo(ECON,PAR,Bgrid,Agrid,Thgrid,tolerance);
    //Printing the corresponding values of the function
    int aa=0;
    int tt=0;
    int step=1;
    for (aa=0; aa<=19; aa=aa+step){
        for (tt=0; tt<=19; tt=tt+step){
            cout << Vy[0][aa][tt]<< " vy" << endl;
            //cout << Vy[1][aa][tt]<< " vy" << endl;
        }
    }
    return 0;
}

