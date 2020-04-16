/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on June 27, 2019, 9:20 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double ninf(double V, double k, double Vh);
void init(double x[],int n); 

int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("countInew2.txt");
    
    double I0=0; // 0
    double gL=0.3; // 8
    double EL=-80; // -80  
    double gNa=1; // 20
    double ENa=60; // 60
    double gK=0.4; // 9
    double EK=-90; // -90
   
    //Na current
    double km=14; // 15 
    double vm=-18; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-25; // -20 
    double tau=3; // 0.152  
    int N0;
    N0=50000000;
    int runs=50;
    double count[runs];
    double v,vs,nf,nfs;
    double th=-25;
    vs=-30;
    nfs=0.2;
    int j,f;
    double dt,dt0;
    dt0=0.0005; 
    init(count,runs);
for(int a=0;a<runs;a++){
    dt=dt0*(a+1);
int N=round(N0/(a+1));
for(j=0;j<N;j++){
    v=I0*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
if(v > th && vs < th){
        count[a]=count[a]+1;
    }
vs=v;
nfs=nf;
}
}
for(f=0;f<runs;f++){
myfile << dt0*(f+1) << " " << count[f] << "\n";
}
 myfile.close();   
    return 0;
}

double ninf(double V, double k, double Vh){
    double f=1/(1+exp((Vh-V)/k));
    return f;
}

void init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}

