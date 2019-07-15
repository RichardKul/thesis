/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on December 27, 2018, 1:54 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <map>
#include <string>
#include <iomanip>

using namespace std;
double average(double x[],int n);
double square(double x[],int n);
double init(double x[],int n);

int main(int argc, char** argv) {
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    myfile.open (argv[2]);
    myfile2.open (argv[3]);
    myfile3.open (argv[4]);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
 
int T=1200;//Simulationszeit
int O=1;//Anzahl der Kraftschritte
int P=200;//Zeit, ab der gemittelt wird
double df=0.02;


double v[100],vs[100];
double x[100],xs[100]; 
double deff[O];
deff[0]=0;
double deff2[O];
deff2[0]=0;
double vav[O];

double dt=0.01;
int N=T/dt;
double F0=atof(argv[1]);//kleinste Kraft 0.55
double kT=0.033;
double gam=0.4;
int Neq=P/dt;

int j,l,a2,f,g,g2,s;
double F;
double xav[O];
xav[0]=0;
double Fano[O];
Fano[0]=0;
double vavt,xavt,x2avt,defft;

//Schleife über versch Kräfte
for(s=0;s<O;s++){
F=F0+s*df;
//Initialisierung von x
init(x,100);
init(xs,100);
//Initialisierung von v
init(v,100);
init(vs,100);
//Schleife über die ersten Zeitschritte
for(j=0;j<Neq;j++){ 
//Schleife über alle Teilchen
    for(l=0;l<100;l++){
	v[l]=vs[l]-gam*vs[l]*dt+F*dt-sin(xs[l])*dt+sqrt(2*kT*gam*dt)*n(gen);
	x[l]=xs[l]+vs[l]*dt;
        xs[l]=x[l];
        vs[l]=v[l];
    }
}
int e=0;
double av=0;
double avv=0;
//Schleife über die letzten Zeitschritte
for(a2=Neq;a2<N;a2++){
    //Schleife über alle Teilchen
    for(l=0;l<100;l++){
        v[l]=vs[l]-gam*vs[l]*dt+F*dt-sin(xs[l])*dt+sqrt(2*kT*gam*dt)*n(gen);
        x[l]=xs[l]+vs[l]*dt;
        vs[l]=v[l];
        xs[l]=x[l];
        }
    
    vavt=average(v,100);//mittl Geschwindigkeit
    xavt=average(x,100);//mittl Position
    square(x,100);
    x2avt=average(x,100);//Mittel der Positionsquadrate
    
    defft=(x2avt-xavt*xavt)/(2*(a2+1)*dt);//Diffusionskoeffizient
        
//Finden des Mittels der letzten Werte
av=(e*av+defft)/(e+1);
avv=(e*avv+vavt)/(e+1);
e=e+1;   
deff2[s]=av;//mittl Diffusionskoeffizient
}
vav[s]=avv;
Fano[s]=deff2[s]/vav[s];
}
    
for(f=0;f<O;f++){
myfile << F0+f*df << " " << deff2[f] <<"\n";
}

for(g=0;g<O;g++){
    myfile2 << F0+g*df << " " << vav[g] << "\n";
}

for(g2=0;g2<O;g2++){
    myfile3 << F0+g2*df << " " <<Fano[g2] << "\n";
}

    myfile.close();
    myfile2.close();
    myfile3.close();
    
    return 0;
}

double average(double x[],int n){ 
    int i;
    double sum=0;
    int e2=0;
    for(i=0;i<n;i++){
        sum=(sum*e2+x[i])/(e2+1);
        e2=e2+1;
    }
    return sum;
}

double square(double x[],int n){
    int j;
    for(j=0;j<n;j++){
        x[j]*=x[j];
    }
}

double init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}
