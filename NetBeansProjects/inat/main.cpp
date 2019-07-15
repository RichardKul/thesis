/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on March 16, 2019, 1:03 PM
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
    myfile.open ("timed.txt");
    myfile2.open ("timev.txt");
    myfile3.open ("timef.txt");
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=0; // 0
    double dI=0.5;
    double gL=1; // 1
    double EL=-70; // -70  
    double gNa=15; // 15
    double ENa=60; // 60
    //double gM=5; // 5
    /*//slow K+ current
    double kinf=5; // 5  
    double vinf=-20; // -20  
    double tauM=30; // 20*/ 
    //Na current
    double km=15; // 15 
    double vm=-40; // -20 
    //fast K+ current
    double kh=-7; // 5 
    double vh=-62; // -20 
    double tau=5; // 0.152  
    int N=500000;
    int Neq=100000;
    int j,f;
    double dt=0.00001;
    double D=3;
    //int spikes=100;
    //double* spike=new double[spikes];
    double th=-25;
    int runs=100;
    int Ivalues=1;
    long long int count[runs],sqcount[runs];
    double Nav,N2av,T;
    double v[runs],hf[runs],hfs[runs],vs[runs],uGe[Ivalues],uFano[Ivalues],uDeff[Ivalues],Deff,Deff2,deff[Ivalues],vav[Ivalues],Fano[Ivalues];
    int s;
for(s=0;s<Ivalues;s++){
    initd(hfs,runs);
    initd(vs,runs);
    init(count,runs);
for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=(I0+s*dI)*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*ninf(vs[a],km,vm)*ninf(vs[a],km,vm)*hfs[a]*(vs[a]-ENa)*dt-sqrt(2*D*dt)*n(gen);
    hf[a]=hfs[a]+(ninf(vs[a],kh,vh)-hfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        count[a]=count[a]+1;
    }
    vs[a]=v[a];
    hfs[a]=hf[a];
    }
}
    int e=0;
    double av=0;
    double avv=0;
    double avd=0;
    double uG=0;
    for(int j2=Neq;j2<N;j2++){
    for(int a=0;a<runs;a++){ 
    v[a]=(I0+s*dI)*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*ninf(vs[a],km,vm)*ninf(vs[a],km,vm)*hfs[a]*(vs[a]-ENa)*dt-sqrt(2*D*dt)*n(gen);
    hf[a]=hfs[a]+(ninf(vs[a],kh,vh)-hfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        count[a]=count[a]+1;
    }
    vs[a]=v[a];
    hfs[a]=hf[a];
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
}
deff[s]=av;//mittl Diffusionskoeffizient
uDeff[s]=sqrt(avd-av*av);
vav[s]=avv;
uGe[s]=sqrt(uG-avv*avv);
Fano[s]=2*deff[s]/avv;   
uFano[s]=2*uDeff[s]/avv+2*av*uGe[s]/(avv*avv);
    } 
    for(f=0;f<Ivalues;f++){
myfile << I0+f*dI << " " << deff[f] << " "<< uDeff[f] <<"\n";
}
 myfile.close(); 
    for(f=0;f<Ivalues;f++){
myfile2 << I0+f*dI << " " << vav[f] << " "<< uGe[f] <<"\n";
}
 myfile2.close();     
 for(f=0;f<Ivalues;f++){
myfile3 << I0+f*dI << " " << Fano[f] << " "<< uFano[f] << "\n";
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
