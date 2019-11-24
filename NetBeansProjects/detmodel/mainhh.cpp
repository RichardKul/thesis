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
double an(double V);
double bn(double V);
double am(double V);
double bm(double V);
double ah(double V);
double bh(double V);

int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("counthh.txt");
    
    double mu=6.5; // 0
    double dI=0.06;
    double gL=0.3; // 8
    double EL=10; // -80  
    double gNa=120; // 20
    double ENa=115; // 60
    double gK=36; // 9
    double EK=-12; // -90
    int N0;
    N0=50000000;
    int runs=51;
    double count[runs];
    double v,vs,nf,nfs,mf,mfs,hf,hfs;
    double th=50;
    vs=0;
    nfs=0.35;
    mfs=0.6;
    hfs=0.6;
    int j,f;
    double dt;
    dt=0.0005; 
    init(count,runs);
for(int a=0;a<runs;a++){
for(j=0;j<N0;j++){
  v=(mu+a*dI)*dt+vs+gL*(EL-vs)*dt+gNa*pow(mfs,3)*hfs*(ENa-vs)*dt+gK*pow(nfs,4)*(EK-vs)*dt;
  nf=nfs+an(vs)*(1-nfs)*dt-bn(vs)*nfs*dt;
  mf=mfs+am(vs)*(1-mfs)*dt-bm(vs)*mfs*dt;
  hf=hfs+ah(vs)*(1-hfs)*dt-bh(vs)*hfs*dt;
if(v > th && vs < th){
        count[a]=count[a]+1;
    }
 
vs=v;
nfs=nf;
mfs=mf;
hfs=hf;
}
}
for(f=0;f<runs;f++){
myfile << mu+dI*f << " " << count[f] << "\n";
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

double an(double V){
  double f=(10-V)/(100*(exp((10-V)/10)-1));
  return f;
}

double am(double V){
  double f=(25-V)/(10*(exp((25-V)/10)-1));
  return f;
}

double ah(double V){
  double f=7*exp(-V/20)/100;
  return f;
}

double bn(double V){
  double f=exp(-V/80)/8;
  return f;
}

double bm(double V){
  double f=4*exp(-V/18);
  return f;
}

double bh(double V){
  double f=1/(exp((30-V)/10)+1);
  return f;
}