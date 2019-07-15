/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on March 26, 2019, 9:43 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

double average(long long int x[],int n);
double square(long long int x[],int n);
double ninf(double V, double k, double Vh);
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
    
    ofstream myfile;

    myfile.open ("chaymodelvm62.txt");
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=0; // 0
    double dI=0.5;
    double gm=8; // 8
    double Vm=-62; // -81  
    double gI=20; // 20
    double VI=60; // 60
    double gK=11; // 9
    double VK=-90; // -90
        //slow K+ current
    double kn=7; // 12  
    double vn=-21; // -26  
    double tau=0.38; // 3 
    //Na current
    double km=10; // 15 
    double vm=-20; // -20 
    //fast K+ current
    double kh=-7; // -7 
    double vh=-3; // 7 
    
    int Neq=3000000;
    int points=10000;
    int sampling=Neq/points;
    int j,f;
    double dt=0.0001;
    double D=0;
    int runs=1;
    double v[runs],nf[runs],nfs[runs],vs[runs],vout[points],hout[points];
    int p=0;

    nfs[0]=0.08;
    vs[0]=-20;
for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=I0*dt+vs[a]-gm*(vs[a]-Vm)*dt-gI*ninf(vs[a],km,vm)*ninf(vs[a],kh,vh)*(vs[a]-VI)*dt-gK*(vs[a]-VK)*nfs[a]*dt-sqrt(2*D*dt)*n(gen);
    nf[a]=nfs[a]+(ninf(vs[a],kn,vn)-nfs[a])*dt/tau;
    if(j==p*sampling){vout[p]=v[a];
    hout[p]=nf[a];
    p=p+1;
    }
    vs[a]=v[a];
    nfs[a]=nf[a];
    }
}

    for(f=0;f<points;f++){
myfile << f*sampling*dt << " " << vout[f] << " "<< hout[f] << "\n";
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



