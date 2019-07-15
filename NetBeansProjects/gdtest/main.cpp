/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on December 28, 2018, 12:03 AM
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

int main(int argc, char** argv) {
    ofstream myfile;
    ofstream myfile2;
    myfile.open ("difff.txt");
    myfile2.open ("geschf.txt");
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
 
int N=50000000;//Anzahl Zeitschritte

double v[100],vs[100];
double x[100],xs[100]; 
double deff;
double deff2;

double dt=0.01;
double F=0.7;//kleinste Kraft
double kT=0.056;
double gam=0.4;

int h,i,j,l,a,b,c,a2,b2,c2,d2,f,g;
double vavt[1000],xavt[1000],x2avt[1000],defft[1000];
//Initialisierung von x
for(h=0;h<100;h++){
x[h]=0;
}
//Initialisierung von v
for(i=0;i<100;i++){
v[i]=0;
}
d2=0;
//Schleife über alle Zeitschritte
for(j=0;j<N;j++){
    if(j==(d2+1)*50000-1){
        d2=d2+1;
         for(l=0;l<100;l++){
             x[l]=xs[l]+vs[l]*dt;
	v[l]=vs[l]-gam*vs[l]*dt+F*dt-sin(xs[l])*dt+sqrt(2*kT*gam*dt)*n(gen);
        vs[l]=v[l];
        xs[l]=x[l];
    }
   double sum4=0;
    for(a2=0;a2<100;a2++){
    sum4=sum4+v[a2];
    }
    vavt[d2]=sum4/100;//mittl Geschwindigkeit
    
    double sum5=0;
    for(b2=0;b2<100;b2++){
    sum5=sum5+x[b2];
    }
    xavt[d2]=sum5/100;//mittl Position
    
    double sum6=0;
    for(c2=0;c2<100;c2++){
    sum6=sum6+x[c2]*x[c2];
    }
    x2avt[d2]=sum6/100;//Mittel der Positionsquadrate
    
    defft[d2]=(x2avt[d2]-xavt[d2]*xavt[d2])/(2*(j+1)*dt);//Diffusionskoeffizient
        }
    else{
        
    for(l=0;l<100;l++){
        x[l]=xs[l]+vs[l]*dt;
	v[l]=vs[l]-gam*vs[l]*dt+F*dt-sin(xs[l])*dt+sqrt(2*kT*gam*dt)*n(gen);
        xs[l]=x[l];
        vs[l]=v[l];
    }
    }
}
vavt[0]=0;xavt[0]=0;x2avt[0]=0;defft[0]=0;
for(f=0;f<1000;f++){
myfile << f*500 << ' ' << defft[f] <<"\n";
}

for(g=0;g<1000;g++){
    myfile2 << g*500 << ' ' << vavt[g] << "\n";
}
    myfile.close();
    myfile2.close();
  
    return 0;
}


