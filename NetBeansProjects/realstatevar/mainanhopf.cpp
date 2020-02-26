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
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
    
    ofstream myfile;
    ofstream myfile2;
    //ofstream myfile3;
    myfile.open ("realstateanhopf5.txt");
    myfile2.open ("paramsanhopf5.txt");
    //myfile3.open (argv[4]
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=45; // 0
    double dI=0.5;
    double gL=1; // 8
    double EL=-78; // -80  
    double gNa=4; // 20
    double ENa=60; // 60
    double gK=4; // 9
    double EK=-90; // -90
    //double gM=5; // 5
    /*//slow K+ current
    double kinf=5; // 5  
    double vinf=-20; // -20  
    double tauM=30; // 20*/ 
    //Na current
    double km=7; // 15 
    double vm=-30; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-45; // -20 
    double tau=1; // 0.152  
    //int N=500000;
    int Neq=10000000;
    int points=200000;
    int sampling=Neq/points;
    int j,f;
    double dt=0.005;
    double D=0.3;
    //int spikes=100;
    //double* spike=new double[spikes];
    // simulate equilibrium point
    int Nsim=5000000;
    double vsim,nfsim,nfsims,vsims,neq,veq;
    nfsims=0.2;
    vsims=-50;
    for(int jsim=0;jsim<Nsim;jsim++){ 
      vsim=I0*dt+vsims-gL*(vsims-EL)*dt-gNa*ninf(vsims,km,vm)*(vsims-ENa)*dt-gK*nfsims*(vsims-EK)*dt;
      nfsim=nfsims+(ninf(vsims,kn,vn)-nfsims)*dt/tau;
      vsims=vsim;
      nfsims=nfsim;
    }
    veq=vsims;
    neq=nfsims;
    /*//simulate instable focus
    int Nfoc=20000000;
    double dtfoc,vfoc,nffoc,nffocs,vfocs,th,thn;
    dtfoc=-dt;
    nffocs=0.5;
    vfocs=-25;
    for(int jfoc=0;jfoc<Nfoc;jfoc++){
      vfoc=I0*dtfoc+vfocs-gL*(vfocs-EL)*dtfoc-gNa*ninf(vfocs,km,vm)*(vfocs-ENa)*dtfoc-gK*nffocs*(vfocs-EK)*dtfoc;
      nffoc=nffocs+(ninf(vfocs,kn,vn)-nffocs)*dtfoc/tau;
      vfocs=vfoc;
      nffocs=nffoc;
    }
    if(vfocs<=0 && vfocs >=-100 && nffocs <= 1 && nffocs >=0){
      th=vfocs;
      thn=nffocs;
    }
    else{
      th=-25;
      thn=0.7;
    }*/
    double th=-30;
    double thn=0.7;
    int runs=1;
    int Ivalues=1;
    long long int count;
    double Nav,N2av,T;
    double mx,nx,nx2,v,nf,nfs,vs,uGe[Ivalues],uFano[Ivalues],uDeff[Ivalues],Deff,Deff2,deff[Ivalues],vav[Ivalues],Fano[Ivalues],vout[points],nout[points],bvout[points];
    int cout[points];
    int p=0;
    int bv2,sv,svlocked;
    sv=0;
//for(s=0;s<Ivalues;s++){
    //initd(nfs,runs);
    nfs=0.3;
    vs=-50;
    sv=0;svlocked=0;
count=0;
for(j=0;j<Neq;j++){ 
  v=I0*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    if(v > th && vs < th){
        sv=1;
    }
    if(sv==1 && nf > thn && nfs < thn){
        count=count+1;
        sv=0;
        bv2=1;
        svlocked=0;
    }
    
    if(bv2==1 && v<veq+5 && nf<neq+0.1){
      if(svlocked==0){
        if(v < veq && vs > veq){
          svlocked=1;
        }}
      if(svlocked==0){
        if(nf > neq && nfs < neq){
          svlocked=-1;
        }
      }
      if(svlocked==1){
        if(nf < neq && nfs > neq){
          bv2=0;
          svlocked=2;
        }
      }
      if(svlocked==-1){
        if(v < veq && vs > veq){
          bv2=0;
          svlocked=2;
        }
      }
    }
   
    if(j==p*sampling){vout[p]=v;
    nout[p]=nf;
    bvout[p]=bv2;
    cout[p]=count;
    p=p+1;
    }
    vs=v;
    nfs=nf;
  }

    for(f=0;f<points;f++){
myfile << f*sampling*dt << " " << vout[f] << " "<< nout[f] << " "<< bvout[f]<<" " << cout[f] << "\n";
}
    myfile2 << th<<" "<<thn<<" "<<veq<<" "<<neq;
 myfile.close(); 
myfile2.close();
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

