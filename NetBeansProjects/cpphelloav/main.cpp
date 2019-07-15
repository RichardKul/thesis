/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on February 13, 2019, 3:58 PM
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
    int filename=atof(argv[1]);
    ostringstream filename_convert;
    filename_convert << filename;
    string filename_str = filename_convert.str();
    string filename_str1="deffaism"+filename_str+".txt";
    string filename_str2="vavaism"+filename_str+".txt";
    string filename_str3="fanoaism"+filename_str+".txt";
    string filename_str4="deffais"+filename_str+".txt";
    string filename_str5="vavais"+filename_str+".txt";
    string filename_str6="fanoais"+filename_str+".txt";
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    ofstream myfile6;
    myfile.open (filename_str1);
    myfile2.open (filename_str2);
    myfile3.open (filename_str3);
    myfile4.open (filename_str4);
    myfile5.open (filename_str5);
    myfile6.open (filename_str6);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
 
int T=5000;//Simulationszeit
//int O=1;//Anzahl der Kraftschritte
int P=1000;//Zeit, ab der gemittelt wird
//double df=0.02;

double dvalues[5]={0.033,0.043,0.056,0.072,0.094};
double v[100],vs[100];
double x[100],xs[100]; 

double dt=0.01;
int N=T/dt;
double F=0.7;//kleinste Kraft 0.55
double kT=dvalues[filename-1];
double gam=0.4;
int Neq=P/dt;
int points=1000;
int sampling=(N-Neq)/points;
int j,l,a2,f,g,g2;
double vavt,xavt,x2avt,defft;
double Deff[points],vav[points],Fano[points],Deffm[points],vavm[points],Fanom[points];

//Schleife über versch Kräfte
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
int p=0;
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
        if(a2==Neq+p*sampling){
    Deff[p]=defft;
    Deffm[p]=av;
    vav[p]=vavt;
    vavm[p]=avv;
    Fano[p]=2*Deff[p]/vav[p];
    Fanom[p]=2*Deffm[p]/vavm[p];
    p++;
    }
}

    
for(f=0;f<points;f++){
myfile << Neq*dt+f*sampling*dt << " " << Deffm[f] <<"\n";
}

for(g=0;g<points;g++){
    myfile2 << Neq*dt+g*sampling*dt<< " " << vavm[g] << "\n";
}

for(g2=0;g2<points;g2++){
    myfile3 << Neq*dt+g2*sampling*dt<< " " <<Fanom[g2] << "\n";
}

for(f=0;f<points;f++){
myfile4 << Neq*dt+f*sampling*dt << " " << Deff[f] <<"\n";
}

for(g=0;g<points;g++){
    myfile5 << Neq*dt+g*sampling*dt<< " " << vav[g] << "\n";
}

for(g2=0;g2<points;g2++){
    myfile6 << Neq*dt+g2*sampling*dt<< " " <<Fano[g2] << "\n";
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
