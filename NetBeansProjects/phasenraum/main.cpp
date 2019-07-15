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
double init(double x[],int n); 

int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("phasenraumstab.txt");
    
    double I=0; // 0
    double gL=8; // 8
    double EL=-80; // -80  
    double gNa=20; // 20
    double ENa=60; // 60
    double gK=9; // 9
    double EK=-90; // -90
    double gM=5; // 5
    //slow K+ current
    double kinf=5; // 5  
    double vinf=-20; // -20  
    double tauM=30; // 20 
    //Na current
    double km=15; // 15 
    double vm=-20; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-25; // -20 
    double tau=0.152; // 0.152  
    int N=1000000;
    int nvalues=11;
    int vvalues=11;
    int runs=nvalues*vvalues;
    double count[runs],vendvalue[runs],nendvalue[runs];
    double v,vs,nf,nfs;
    double th=-25;
    double vstart=-80;
    double vend=0;
    double dv=(vend-vstart)/(vvalues-1);
    double nstart=0.1;
    double nend=0.9;
    double dn=(nend-nstart)/(nvalues-1);
    nfs=0;
    int j,f,f2;
    double dt=0.00001; 
    init(count,runs);
    for(int b=0;b<vvalues;b++){
for(int a=0;a<nvalues;a++){
    vs=vstart+b*dv;
    nfs=nstart+a*dn;
for(j=0;j<N;j++){
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
if(v > th && vs < th){
        count[b*nvalues+a]=count[b*nvalues+a]+1;
    }
vs=v;
nfs=nf;
}
    vendvalue[b*nvalues+a]=v;
    nendvalue[b*nvalues+a]=nf;
}
    }
    for(f=0;f<nvalues;f++){
        for(f2=0;f2<vvalues;f2++){
myfile << vstart+f2*dv << " " << nstart+f*dn << " " << count[f2*nvalues+f] << " " 
        << nendvalue[f2*nvalues+f] << " " << vendvalue[f2*nvalues+f] << "\n";
}
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

