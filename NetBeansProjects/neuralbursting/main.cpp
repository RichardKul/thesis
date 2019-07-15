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
    myfile.open ("deffD5.txt");
    myfile2.open ("vavD5.txt");
    myfile3.open ("fanoD5.txt");
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=0; // 0
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
    long int N=5000000;
    int points=1000;
    int sample=N/points;
    int j,f;
    double dt=0.00001;
    double D=5;
    //int spikes=100;
    //double* spike=new double[spikes];
    double th=-25;
    int runs=50;
    long long int count[runs],sqcount[runs];
    double Nav[points],N2av[points],vav[points],Fano[points];
    double v[runs],nf[runs],nfs[runs],vs[runs],Deff[points];
    int p=1;
    Deff[0]=0;
    initd(nfs,runs);
    initd(vs,runs);
    init(count,runs);
for(j=0;j<N;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=I0*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*(vs[a]-ENa)*dt-gK*nfs[a]*(vs[a]-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf[a]=nfs[a]+(ninf(vs[a],kn,vn)-nfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        count[a]=count[a]+1;
    }
    vs[a]=v[a];
    nfs[a]=nf[a];
    }
    if(j==p*sample){
    Nav[p]=average(count,runs);
    for(int j2=0;j2<runs;j2++){
        sqcount[j2]=count[j2]*count[j2];
    }
    N2av[p]=average(sqcount,runs);
    Deff[p]=(N2av[p]-Nav[p]*Nav[p])/(2*j*dt);
    vav[p]=Nav[p]/(j*dt);
    Fano[p]=2*Deff[p]/vav[p];
    p++;
    }    
}    
    
    for(f=0;f<points;f++){
myfile << f*sample*dt << " " << Deff[f] << "\n";
}
 myfile.close(); 
    for(f=0;f<points;f++){
myfile2 << f*sample*dt << " " << vav[f] << "\n";
}
 myfile2.close();     
 for(f=0;f<points;f++){
myfile3 << f*sample*dt << " " << Fano[f] << "\n";
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

double initd(double x[],int n){
    for(int k2=0;k2<n;k2++){
        x[k2]=0;
    }
}
