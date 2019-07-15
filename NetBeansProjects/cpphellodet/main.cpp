/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on March 25, 2019, 4:02 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

double average(long long int x[],int n);
double square(long long int x[],int n);
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
    
    ofstream myfile;
    myfile.open ("mechdet2.txt");

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
int O=1;//Anzahl der Kraftschritte
int P=20;//Zeit, ab der gemittelt wird

double dt=0.0001;
double F=0.7;//kleinste Kraft 0.55
double kT=0;
double gam=0.4;
int Neq=P/dt;

int j,l,a2,f,g,g2,s;

double vavt,xavt,x2avt,defft;
    int points=10000;
    int sampling=Neq/points;
    double v,x,xs,vs,vout[points],xout[points];
    int p=0;
    vs=1.4;

//Schleife über die ersten Zeitschritte
for(j=0;j<Neq;j++){ 
//Schleife über alle Teilchen
	v=vs-gam*vs*dt+F*dt-sin(xs)*dt+sqrt(2*kT*gam*dt)*n(gen);
	x=xs+vs*dt;

   
    if(j==p*sampling){vout[p]=v;
    xout[p]=x;
    p=p+1;
    }
    xs=x;
    vs=v;
    }

    for(f=0;f<points;f++){
myfile << f*sampling*dt << " " << vout[f] << " "<< xout[f] << "\n";
}
 myfile.close(); 

    return 0;
}

double ninf(double V, double k, double Vh){
    double f=1/(1+exp((Vh-V)/k));
    return f;
}

double average(long long int x[],int n){ 
    int i;
    double sum=0;
    int e2=0;
    for(i=0;i<n;i++){
        sum=(sum*e2+x[i])/(e2+1);
        e2=e2+1;
    }
    return sum;
}

double square(long long int x[],int n){
    int j;
    for(j=0;j<n;j++){
        x[j]=x[j]*x[j];
    }
}

double init(long long int x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}

double initd(double x[],int n){
    for(int k2=0;k2<n;k2++){
        x[k2]=0;
    }
}

