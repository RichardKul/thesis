/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on June 27, 2019, 10:13 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double ninf(double V, double k, double Vh);
double init(double x[],int n); 

int main(int argc, char** argv) {
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    ofstream myfile;
    myfile.open ("countI2.txt");
    
    double I0=0; // 0
    double gL=0.0029; // 8
    double EL=-80; // -80  
    double gNa=0.01; // 20
    double ENa=60; // 60
    double gK=0.004; // 9
    double EK=-90; // -90
   
    //Na current
    double km=14; // 15 
    double vm=-18; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-25; // -20 
    double tau=400; // 0.152  
    int N0;
    N0=5000000;
    int points=10000;
    int runs=5;
    int sampling=(N0*runs)/points;
    double count[runs],vout[points],hout[points];
    double v,vs,nf,nfs,D;
    D=0.01;
    double th=-25;
    vs=-30;
    nfs=0.2;
    int j,f,p;
    p=0;
    double dt;
    dt=0.01; 
    init(count,runs);
for(int a=0;a<runs;a++){
for(j=0;j<N0;j++){
    v=I0*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
if(v > th && vs < th){
        count[a]=count[a]+1;
    }
    if(j+a*N0==p*sampling){vout[p]=v;
    hout[p]=nf;
    p=p+1;
    }
vs=v;
nfs=nf;
}
}
for(f=0;f<points;f++){
myfile << dt*sampling*f << " " << vout[f] << " "<< hout[f] << "\n";
}
 myfile.close();   
    return 0;
}

double ninf(double V, double k, double Vh){
    double f=1/(1+exp((Vh-V)/k));
    return f;
}

double init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}
