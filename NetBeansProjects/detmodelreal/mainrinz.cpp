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
double minf(double V);
double ninf(double V);
double hinf(double V);

int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("countrinz2.txt");
    
    double I=-15; // -15
    double S=1.27; // 1.27
    double phi=3.82; //3.82
    double dI=0.75; // 0.75 
    double gL=0.3; // 0.3
    double EL=10; // 10 
    double gNa=120; // 120
    double ENa=115; // 115
    double gK=36; // 36
    double EK=12; // -12
    int N0;
    N0=500000000;
    int runs=50;
    double count[runs];
    double v,vs,wf,wfs;
    double th=40;
    vs=30;
    wfs=0.2;
    int j,f;
    double dt,dt0;
    dt0=0.00005; 
    init(count,runs);
for(int a=0;a<runs;a++){
    dt=dt0*(a+1);
int N=round(N0/(a+1));
for(j=0;j<N;j++){
  v =I*dt+vs-gL*(vs-EL)*dt-gNa*pow(minf(vs),3)*(1-wfs)*(vs-ENa)*dt-gK*pow(wfs/S,4)*(vs-EK)*dt;
  wf=wfs+phi*((S*(ninf(vs)+S*(1-hinf(vs)))/(1+(S*S)))-wfs)*dt/(5*exp(-pow((vs+10)/55,2))+1);
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
myfile << dt0*f << " " << count[f] << "\n";
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

double ninf(double V){
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
