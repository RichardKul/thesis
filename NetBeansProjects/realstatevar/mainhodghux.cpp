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
double an(double V);
double bn(double V);
double am(double V);
double bm(double V);
double ah(double V);
double bh(double V);
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
  
  ofstream myfile;
  //ofstream myfile2;
  //ofstream myfile3;
  myfile.open ("realstatehh.txt");
  //myfile2.open ("statechange.txt");
  //myfile3.open (argv[4]
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> n(0,1);
  
  double mu=6.5; // 0
  double dI=0.5;
  double gL=0.3; // 8
  double EL=10; // -80  
  double gNa=120; // 20
  double ENa=115; // 60
  double gK=36; // 9
  double EK=-12; // -90
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
  double dt=0.0005;
  double D=0.5;
  //int spikes=100;
  //double* spike=new double[spikes];
  int runs=1;
  int Ivalues=1;
  long long int sqcount[runs];
  double Nav,N2av,T;
  double v[runs],nf[runs],nfs[runs],vs[runs],mf[runs],mfs[runs],hf[runs],hfs[runs],uGe[Ivalues],uFano[Ivalues],uDeff[Ivalues],Deff,Deff2,deff[Ivalues],vav[Ivalues],Fano[Ivalues],vout[points],nout[points],mout[points],hout[points],lastcount,prelastchange,lastchange;
  int p=0;
  int bv,bvp,sv;
  lastcount=0;
  prelastchange=0;lastchange=0;
  sv=0;
  //for(s=0;s<Ivalues;s++){
  //initd(nfs,runs);
  mfs[0]=0.06;
  hfs[0]=0.6;
  nfs[0]=0.35;
  vs[0]=0;
  for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
      v[a]=mu*dt+vs[a]+gL*(EL-vs[a])*dt+gNa*pow(mfs[a],3)*hfs[a]*(ENa-vs[a])*dt+gK*pow(nfs[a],4)*(EK-vs[a])*dt+sqrt(2*D*dt)*n(gen);
      nf[a]=nfs[a]+an(vs[a])*(1-nfs[a])*dt-bn(vs[a])*nfs[a]*dt;
      mf[a]=mfs[a]+am(vs[a])*(1-mfs[a])*dt-bm(vs[a])*mfs[a]*dt;
      hf[a]=hfs[a]+ah(vs[a])*(1-hfs[a])*dt-bh(vs[a])*hfs[a]*dt;
      if(j==p*sampling){vout[p]=v[a];
        nout[p]=nf[a];
        mout[p]=mf[a];
        hout[p]=hf[a];
        p=p+1;
      }
      vs[a]=v[a];
      nfs[a]=nf[a];
      hfs[a]=hf[a];
      mfs[a]=mf[a];
    }
  }
  
  for(f=0;f<points;f++){
    myfile << f*sampling*dt << " " << vout[f] << " "<< nout[f] <<" " << mout[f] << " "<< hout[f] << "\n";
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

double an(double V){
  double f=(10-V)/(100*(exp((10-V)/10)-1));
  return f;
}

double am(double V){
  double f=(25-V)/(10*(exp((25-V)/10)-1));
  return f;
}

double ah(double V){
  double f=7*exp(-V/20)/100;
  return f;
}

double bn(double V){
  double f=exp(-V/80)/8;
  return f;
}

double bm(double V){
  double f=4*exp(-V/18);
  return f;
}

double bh(double V){
  double f=1/(exp((30-V)/10)+1);
  return f;
}