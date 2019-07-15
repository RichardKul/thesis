/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on June 27, 2019, 9:20 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <valarray>

using namespace std;

double ninf(double V, double k, double Vh);
void init(double x[],int n); 

int main(int argc, char** argv) {
    
//ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
//    myfile.open ("countInew.txt");
    myfile2.open("voltageplusstate.txt");
    myfile3.open("param.txt");
    myfile4.open ("time.txt");
    myfile5.open ("timenew.txt");
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I=0; // 0
    double gL=0.3; // 8
    double EL=-80; // -80  
    double gNa=1; // 20
    double ENa=60; // 60
    double gK=0.4; // 9
    double EK=-90; // -90
   
    //Na current
    double km=14; // 15 
    double vm=-18; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-25; // -20 
    double tau=3; // 0.152  
    int N0;
    N0=5000000;
    int runs=50;
    int points=100000;
    int sampling=N0*runs/points;
    double th=-25;
    double thn=0.7;
  double D=1;
    int j,f;
    double dt;
    dt=0.0001;
    int timepoints=100000;
    int Nsim=5000000;
    double vsim,nfsim,nfsims,vsims,neq,veq;
    nfsims=0.4;
    vsims=-70;
    for(int jsim=0;jsim<Nsim;jsim++){ 
      vsim=I*dt+vsims-gL*(vsims-EL)*dt-gNa*ninf(vsims,km,vm)*(vsims-ENa)*dt-gK*nfsims*(vsims-EK)*dt;
      nfsim=nfsims+(ninf(vsims,kn,vn)-nfsims)*dt/tau;
      vsims=vsim;
      nfsims=nfsim;
    }
    veq=vsims;
    neq=nfsims;
    
    long long int count;
    std::valarray<double> state(0.0,timepoints),state2(0.0,timepoints),times(0.0,timepoints),times2(0.0,timepoints);
    double Nav,Nava,Nava2,N2ava,T,vout[points],nout[points],bvout[points],bv2out[points],cout[points];
    double mx,nx,nx2,v,nf,nfs,vs,Deff,vav,vavalt,sigmav,Fano,lastcount;
    std::valarray<int> avalue(timepoints),avalue2(timepoints),j2value(timepoints),j2value2(timepoints),jrvalue(timepoints),jrvalue2(timepoints);
    avalue[0]=0;
    j2value[0]=0;
    jrvalue[0]=0;
    avalue2[0]=0;
    j2value2[0]=0;
    jrvalue2[0]=0;
    count=0;
    int bv,bv2,lastbv,lastbv2,sv,svlocked,t,t2,tv,j22,p;
    times[0]=0;times2[0]=0;lastcount=0;Nava=0;Nava2=0;sv=0;t=0;t2=0;tv=0;svlocked=0;
    mx=0.0309-0.0047*I;
    nx=1.654+0.01*I;
    nx2=mx*63;
    p=0;
    vs=-30;
    nfs=0.2;
    if(nfs>mx*vs+nx2){
      bv=0;
      bv2=0;
    }
    else{
      bv=1;
      bv2=1;
    }
lastbv=bv;
lastbv2=bv2;
for(int a=0;a<runs;a++){
for(j=0;j<N0;j++){
  j22=a*N0+j;
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen);
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    if(v > th && vs < th){
      sv=1;
    }
    if(sv==1 && nf > thn && nfs < thn){
      count=count+1;
      sv=0;
      bv2=1;
      svlocked=0;
      lastcount=j22;
    }
    if(bv2==1 && v<veq+5 && nf<neq+0.1){
      if(svlocked==0){
        if(v < veq && vs > veq){
          svlocked=1;
        }}
      if(svlocked==0){
        if(nf < neq && nfs > neq){
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
    if((bv==1 && nf>mx*v+nx2) || (bv==1 && (j22-lastcount)>20/dt)){
      bv=0;
    }
    if(bv==0 && nf<mx*v+nx){
      bv=1;  
    }
    if(j22==4000*tv){
      if((t<timepoints) && (0.5<bv+lastbv) && (1.5>bv+lastbv)){
        times[t]=j22*dt;
        state[t]=bv;
        avalue[t]=a;
        jrvalue[t]=j;
        t=t+1;
      }
      tv=tv+1;
      lastbv=bv;
    }
    if((t2<timepoints) && (0.5<bv2+lastbv2) && (1.5>bv2+lastbv2)){
      times2[t2]=j22*dt;
      state2[t2]=bv2;
      avalue2[t2]=a;
      jrvalue2[t2]=j;
      t2=t2+1;
      lastbv2=bv2;
    }
    if(j22==p*sampling){vout[p]=v;
      nout[p]=nf;
      bvout[p]=bv;
      bv2out[p]=bv2;
      cout[p]=count;
      p=p+1;
    }
    vs=v;
    nfs=nf;
}
}
/*for(f=0;f<runs;f++){
myfile << dt0*f << " " << count[f] << "\n";
}
 myfile.close(); */

for(int rv=0;rv<timepoints;rv++){
  myfile4 << times[rv] <<" " << avalue[rv] << " " << jrvalue[rv] <<" " <<state[rv] <<"\n";
}
myfile4.close();

for(int rv=0;rv<timepoints;rv++){
  myfile5 << times2[rv] <<" " << avalue2[rv] << " " << jrvalue2[rv] <<" " <<state2[rv] << "\n";
}
myfile5.close();

for(f=0;f<points;f++){
  myfile2 << f*sampling*dt << " " << vout[f] << " "<< nout[f] << " "<< bvout[f]<<" " << cout[f] << "\n";
}
myfile2.close(); 
myfile3 << veq << " "<< neq << "\n";
myfile3.close();

    return 0;
}

double ninf(double V, double k, double Vh){
    double f=1/(1+exp((Vh-V)/k));
    return f;
}

void init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}

