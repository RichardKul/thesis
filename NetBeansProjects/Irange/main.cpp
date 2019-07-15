/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on January 9, 2019, 2:23 PM
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
    
    double I0=-1; // 0
    double dI=0.2;
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
    int N=100000;
    double v,nf;
    double vs=0,nfs=0;
    int j,f;
    double dt=0.0004;
    double D=atof(argv[1]);
    double th=-25;
    int runs=50;
    long long int count[runs];
    int Ivalues=20;
    double Fano[Ivalues],Deff[Ivalues],avN2[Ivalues];
    double Nav,N2av;
    
for(int b=0;b<Ivalues;b++){ 
    init(count,runs);
for(int a=0;a<runs;a++){ 
for(j=0;j<N;j++){ 
    v=(I0+b*dI)*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    if(v > th && vs < th){
        count[a]=count[a]+1;
    }
    vs=v;
    nfs=nf;
    }
}
    Nav=average(count,runs);
    square(count,runs);
    N2av=average(count,runs);
    Fano[b]=(N2av-Nav*Nav)/Nav;
    Deff[b]=(N2av-Nav*Nav)/;
    avN2[b]=N2av;
}    
    
    
    for(f=0;f<Ivalues;f++){
myfile << I0+f*dI << " " << Fano[f] << "\n";
}
 myfile.close(); 

for(int x=0;x<Ivalues;x++){
    myfile2 << I0+x*dI << " " << avN[x] << "\n";
} 
 myfile2.close();
 
for(int y=0;y<Ivalues;y++){
    myfile3 << I0+y*dI << " " << avN2[y] << "\n";
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
    for(i=0;i<n;i++){
        sum=sum+x[i];
    }
    double av=sum/n;
    return av;
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
