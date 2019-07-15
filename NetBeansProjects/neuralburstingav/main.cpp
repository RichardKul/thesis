/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on February 8, 2019, 1:40 PM
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
    ofstream myfile2;
    ofstream myfile3;
    myfile.open ("defft.txt");
    myfile2.open ("vavt.txt");
    myfile3.open ("fanot.txt");
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=0; // 0
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
    int N=500000;
    int Neq=100000;
    int j,f;
    double dt=0.00001;
    double D=5;
    //int spikes=100;
    //double* spike=new double[spikes];
    double th=-25;
    int runs=100;
    long long int count[runs],sqcount[runs];
    double Nav,N2av,T;
    int points=1000;
    int sampling=(N-Neq)/points;
    double v[runs],nf[runs],nfs[runs],vs[runs],Deff,Deff2,Deffv[points],vav[points],Fano[points];
    initd(nfs,runs);
    initd(vs,runs);
    init(count,runs);
for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=I0*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*(vs[a]-ENa)*dt-gK*nfs[a]*(vs[a]-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf[a]=nfs[a]+(ninf(vs[a],kn,vn)-nfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        count[a]=count[a]+1;
    }
    vs[a]=v[a];
    nfs[a]=nf[a];
    }
}
    int e=0;
    double av=0;
    double avv=0;
    double avd=0;
    double uG=0;
    int p=0;
    for(int j2=Neq;j2<N;j2++){
            for(int a=0;a<runs;a++){ 
    v[a]=I0*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*(vs[a]-ENa)*dt-gK*nfs[a]*(vs[a]-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf[a]=nfs[a]+(ninf(vs[a],kn,vn)-nfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        count[a]=count[a]+1;
    }

    vs[a]=v[a];
    nfs[a]=nf[a];
    }
    Nav=average(count,runs);
    for(int j4=0;j4<runs;j4++){
        sqcount[j4]=count[j4]*count[j4];
    }
    N2av=average(sqcount,runs);
    T=2*j2*dt;
    Deff=(N2av-Nav*Nav)/T;
    Deff2=Deff*Deff;
  av=(e*av+Deff)/(e+1);
avv=(e*avv+Nav/(T/2))/(e+1);
uG=(e*uG+4*Nav*Nav/(T*T))/(e+1);
avd=(e*avd+Deff2)/(e+1);
e=e+1;
        if(j2==Neq+p*sampling){
    Deffv[p]=Deff;
    vav[p]=avv;
    Fano[p]=2*Deffv[p]/vav[p];
    p++;
    }
}
Deffv[999]=Deff;
vav[999]=avv;
Fano[999]=2*Deff/avv;
    for(f=0;f<points;f++){
myfile << Neq*dt+f*sampling*dt << " " << Deffv[f] << "\n";
}
 myfile.close(); 
    for(f=0;f<points;f++){
myfile2 << Neq*dt+f*sampling*dt << " " << vav[f] <<"\n";
}
 myfile2.close();     
 for(f=0;f<points;f++){
myfile3 << Neq*dt+f*sampling*dt << " " << Fano[f] << "\n";
}
 myfile3.close(); 
 
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
