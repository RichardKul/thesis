/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on July 2, 2019, 5:19 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

double average(long long int x[],int n);
double square(long long int x[],int n);
double ninf(double V, double k, double Vh);
double minf(double V);
double ninf(double V);
double hinf(double V);
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
    
    ofstream myfile;
    //ofstream myfile2;
    //ofstream myfile3;
    myfile.open ("realstaterinzel.txt");
    //myfile2.open ("statechange.txt");
    //myfile3.open (argv[4]
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I=-15; // 0
    double S=1.27;
    double phi=3.82;
    double dI=0.5;
    double gL=0.3; // 8
    double EL=10; // -80  
    double gNa=120; // 20
    double ENa=115; // 60
    double gK=36; // 9
    double EK=12; // -90
    //double gM=5; // 5
    /*//slow K+ current
    double kinf=5; // 5  
    double vinf=-20; // -20  
    double tauM=30; // 20*/ 
    //Na current
    //int N=500000;
    int Neq=30000000;
    int points=100000;
    int sampling=Neq/points;
    int j,f;
    double dt=0.05;
    double D=35;
    //int spikes=100;
    //double* spike=new double[spikes];
    int runs=1;
    int Ivalues=1;
    long long int sqcount[runs];
    double Nav,N2av,T;
    double v[runs],wf[runs],wfs[runs],vs[runs],uGe[Ivalues],uFano[Ivalues],uDeff[Ivalues],Deff,Deff2,deff[Ivalues],vav[Ivalues],Fano[Ivalues],vout[points],wout[points],lastcount,prelastchange,lastchange;
    int p=0;
    int bv,bvp,sv;
    lastcount=0;
    prelastchange=0;lastchange=0;
    sv=0;
//for(s=0;s<Ivalues;s++){
    //initd(nfs,runs);
    wfs[0]=0.2;
    vs[0]=20;
for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=I*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*pow(minf(v[a]),3)*(1-wfs[a])*(vs[a]-ENa)*dt-gK*pow(wfs[a]/S,4)*(vs[a]-EK)*dt+sqrt(2*D*dt)*n(gen);
    wf[a]=wfs[a]+phi*((S*(ninf(vs[a])+S*(1-hinf(vs[a])))/(1+(S*S)))-wfs[a])*dt/(5*exp(-pow((vs[a]+10)/55,2))+1);
    if(wf[a]>1){
      wf[a]=1;
    }
    if(j==p*sampling){vout[p]=v[a];
    wout[p]=wf[a];
    p=p+1;
    }
    vs[a]=v[a];
    wfs[a]=wf[a];
    }
}

    for(f=0;f<points;f++){
myfile << f*sampling*dt << " " << vout[f] << " "<< wout[f] <<" " <<  "\n";
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

