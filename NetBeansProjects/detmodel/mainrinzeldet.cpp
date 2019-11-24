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
#include <math.h>

using namespace std;

double init(double x[],int n); 
double Ninf(double V);
double minf(double V);
double hinf(double V);

int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("countrinzel.txt");
    
    double I0=-21; // 0
    double S=1.27;
    double phi=3.82;
    double dI=0.3;
    double gL=0.3; // 8
    double EL=10; // -80  
    double gNa=120; // 20
    double ENa=115; // 60
    double gK=36; // 9
    double EK=12; // -90
    int N0;
    N0=50000000;
    int runs=51;
    double count[runs];
    double v,vs,wf,wfs;
    double th=40;
    vs=30;
    wfs=0.2;
    int j,f;
    double dt;
    dt=0.0005; 
    init(count,runs);
for(int a=0;a<runs;a++){
for(j=0;j<N0;j++){
  v =(I0+a*dI)*dt+vs-gL*(vs-EL)*dt-gNa*pow(minf(vs),3)*(1-wfs)*(vs-ENa)*dt-gK*pow(wfs/S,4)*(vs-EK)*dt;
  wf=wfs+phi*((S*(Ninf(vs)+S*(1-hinf(vs)))/(1+(S*S)))-wfs)*dt/(5*exp(-pow((vs+10)/55,2))+1);
  if(wf >1){
    wf=1;
  }
if(v > th && vs < th){
        count[a]=count[a]+1;
    }
vs=v;
wfs=wf;
}
}
for(f=0;f<runs;f++){
myfile << I0+dI*f << " " << count[f] << "\n";
}
 myfile.close();   
    return 0;
}

double init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}

double Ninf(double V){
  double f=((10-V)/(100*(exp((10-V)/10)-1)))/((10-V)/(100*(exp((10-V)/10)-1))+exp(-V/80)/8);
  return f;
}

double minf(double V){
  double f=((25-V)/(10*(exp((25-V)/10)-1)))/((25-V)/(10*(exp((25-V)/10)-1))+4*exp(-V/18));
  return f;
  
}

double hinf(double V){
  double f=(7*exp(-V/20)/100)/(7*exp(-V/20)/100+1/(exp((30-V)/10)+1));
  return f;
}