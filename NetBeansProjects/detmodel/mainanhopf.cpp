/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on January 8, 2019, 12:48 PM
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
    myfile.open ("countanhopf.txt");
    
    double I0=43; // 0
    double gL=1; // 8
    double EL=-78; // -80  
    double gNa=4; // 20
    double ENa=60; // 60
    double gK=4; // 9
    double EK=-90; // -90
    
    //Na current
    double km=7; // 15 
    double vm=-30; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-45; // -20 
    double tau=1; // 0.152  
    
    int N0;
    N0=50000000;
    int runs=51;
    double istep=0.1;
    double count[runs];
    double v,vs,nf,nfs;
    double th=-60;
    vs=-50;
    nfs=0.7;
    int j,f;
    double dt;
    dt=0.005; 
    init(count,runs);
for(int a=0;a<runs;a++){
for(j=0;j<N0;j++){
    v=(I0+istep*a)*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
if(v > th && vs < th){
        count[a]=count[a]+1;
    }
vs=v;
nfs=nf;
}
}
for(f=0;f<runs;f++){
myfile << I0+istep*f << " " << count[f] << "\n";
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
