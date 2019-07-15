/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on April 9, 2019, 9:48 AM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
using namespace std;

double average(long long int x[],int n);
double square(long long int x[],int n);
double ninf(double V, double k, double Vh);
double init(long long int x[],int n); 
double initd(double x[],int n);
int main(int argc, char** argv) {
    
    double D=atof(argv[1]);
    double d=D*10;
    double nr=atof(argv[2]);
    double I=-1+0.25*nr;
    ostringstream d_convert;
    ostringstream nr_convert;
    d_convert << d;
    nr_convert << nr;
    string d_str = d_convert.str();
    string nr_str = nr_convert.str();
    string filenamef="fa"+d_str+nr_str+".txt";
    string filenamed="da"+d_str+nr_str+".txt";
    string filenameg="ga"+d_str+nr_str+".txt";
    string filenametime="timea"+d_str+nr_str+".txt";
    string filenamevalues="valuesa"+d_str+nr_str+".txt";
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    myfile.open (filenamed);
    myfile2.open (filenameg);
    myfile3.open (filenamef);
    myfile4.open (filenametime);
    myfile5.open (filenamevalues);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    //double I0=-1; // 0
    //double dI=0.2;
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
    int N=50000000; //5*10⁷
    int Neq=10000000; //10⁷
    int points=10000;
    
    int j;
    double dt=0.00001;
    //int spikes=100;
    //double* spike=new double[spikes];
    double th=-25;
    double thn=0.5;
    int runs=2; //5000
    int sampling=runs*(N-Neq)/points;
    //int Ivalues=25;
long long int count;
double Nav,N2av,T;
double mx,nx,nx2,v,nf,nfs,vs,Deff,vav,Fano,state[10000],times[10000],lastcount,prelastchange,lastchange,lasttime;
double vout[points],nout[points],bvout[points],cout[points];
int bv,lastbv,sv,t,p,j22,tv;
tv=0;
lastcount=0;
prelastchange=0;
lastchange=0;
lasttime=0;
sv=0;
t=1;
times[0]=0;
p=0;
//int s;
//for(s=0;s<Ivalues;s++){
nfs=0;
vs=-60;
mx=0.0285-0.0006*I;
nx=1.5-I/80;
nx2=1.65-I/80;
if(nfs>mx*vs+nx2){
  bv=0;
}
else{
  bv=1;
}
lastbv=bv;
state[0]=bv;
count=0;
for(int a=0;a<runs;a++){ 
  for(j=0;j<Neq;j++){ 
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt/*+sqrt(2*D*dt)*n(gen)*/;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    vs=v;
    nfs=nf;
  }
}
long long int e=0;
for(int a=0;a<runs;a++){
  count=0;
  for(int j2=0;j2<N-Neq;j2++){
      j22=a*(N-Neq)+j2;
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt/*+sqrt(2*D*dt)*n(gen)*/;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    if(v > th && vs < th){
      sv=1;
    }
    if(sv==1 && nf > thn && nfs < thn){
      count=count+1;
      sv=0;
      lastcount=j22;
    }
    if((bv==1 && nf>mx*v+nx2) || (bv==1 && (j22-lastcount)*dt>2)){
      //if((j22-lastchange)*dt<0.4){
        //bv=0;
        //lastchange=prelastchange;
      //}
      //else{
        bv=0; 
        //prelastchange=lastchange;
        //lastchange=j22;
             // if(t<10000){
         //         if(lastchange==times[t-1]*100000){}
         //         else{
        //times[t]=lastchange*dt;
        //t=t+1;
        //}
    }
if(bv==0 && nf<mx*v+nx){
      //if((j22-lastchange)*dt<0.4){
       // bv=1;
        //lastchange=prelastchange;
      //}
      //else{
        bv=1;
        //prelastchange=lastchange;
        //lastchange=j22;
        //if(lastchange==times[t-1]*100000){}
        //          else{
        //  if(t<10000){
        //times[t]=lastchange*dt;
        //t=t+1;
        //}   
    }          
    if(j22==40000*tv){
       if((t<10000) && (0.5 < bv+lastbv) && (1.5 > bv+lastbv)){
        times[t]=j22*dt;
        state[t]=bv;
        t=t+1;
        }
       tv=tv+1;
       lastbv=bv;
    }
    if(j22==p*sampling){vout[p]=v;
    nout[p]=nf;
    bvout[p]=bv;
    cout[p]=count;
    p=p+1;
    }
    vs=v;
    nfs=nf;
    Nav=(Nav*e+count/runs)/(e+1);
    N2av=(N2av*e+(count*count)/runs)/(e+1);
    e=e+1;
  }
} 
T=2*N*dt;
Deff=(N2av-Nav*Nav)/T;
vav=2*Nav/T;
Fano=2*Deff/vav;   
//} 
myfile << I << " " << Deff << "\n";

myfile.close(); 

myfile2 << I << " " << vav <<"\n";

myfile2.close();     

myfile3 << I << " " << Fano << "\n";
myfile3.close(); 

for(int rv=0;rv<10000;rv++){
  myfile4 << times[rv] << " "<< state[rv]<< "\n";
}
myfile4.close();

for(int k=0;k<points;k++){
    myfile5 << k*sampling*dt << " " << vout[k] << " "<< nout[k] << " " << bvout[k] << " " << cout[k]<< "\n";
}
myfile5.close();
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

