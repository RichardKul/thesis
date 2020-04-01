/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on March 18, 2019, 6:14 PM
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
    //ofstream myfile2;
    //ofstream myfile3;
    myfile.open ("inaprealanhopfunstable2.txt");
    //myfile2.open ("statechange.txt");
    //myfile3.open (argv[4]
     
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I0=46; // 0
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
    int Neq=3000000; //30000000
    int points=100000;//100000
    int sampling=Neq/points;
    int j,f;
    double dt=0.005;
    double D=0;
    //int spikes=100;
    //double* spike=new double[spikes];
    double th=-25;
    double thn=0.5;
    int runs=1;
    int Ivalues=1;
    long long int count[runs],sqcount[runs];
    double Nav,N2av,T;
    double mx,nx,nx2,v[runs],nf[runs],nfs[runs],vs[runs],uGe[Ivalues],uFano[Ivalues],uDeff[Ivalues],Deff,Deff2,deff[Ivalues],vav[Ivalues],Fano[Ivalues],vout[points],nout[points],bvout[points],cout[points],lastcount,prelastchange,lastchange;
    int p=0;
    int bv,bvp,sv;
    lastcount=0;
    prelastchange=0;lastchange=0;
    sv=0;
//for(s=0;s<Ivalues;s++){
    //initd(nfs,runs);
    nfs[0]=0.433;
    vs[0]=-50;
    //vs[0]=-50.5567; 
    //nfs[0]=0.177142;
    mx=0.0285-0.0006*I0;
    nx=1.5-I0/80;
    nx2=1.6-I0/80;
    if(nfs[0]>mx*vs[0]+nx2){
        bv=0;
    }
    else{
        bv=1;
    }
    init(count,runs);
for(j=0;j<Neq;j++){ 
    for(int a=0;a<runs;a++){ 
    v[a]=I0*dt+vs[a]-gL*(vs[a]-EL)*dt-gNa*ninf(vs[a],km,vm)*(vs[a]-ENa)*dt-gK*nfs[a]*(vs[a]-EK)*dt;
    nf[a]=nfs[a]+(ninf(vs[a],kn,vn)-nfs[a])*dt/tau;
    if(v[a] > th && vs[a] < th){
        sv=1;
    }
    if(sv==1 && nf[a] > thn && nfs[a] < thn){
        count[a]=count[a]+1;
        sv=0;
        lastcount=j;
    }
    if(( bv==1 && nf[a]>mx*v[a]+nx2) || (j-lastcount)*dt>2){
        if((j-lastchange)*dt<0.4){
                bv=0;
                lastchange=prelastchange;
        }
        else{
            bv=0;
            prelastchange=lastchange;
            lastchange=j;
        }
    }
    if(bv==0 && nf[a]<mx*v[a]+nx){
        if((j-lastchange)*dt<0.4){
        bv=1;
        lastchange=prelastchange;
        }
        else{
            bv=1;
            prelastchange=lastchange;
            lastchange=j;
        }  
    }
    if(j==p*sampling){vout[p]=v[a];
    nout[p]=nf[a];
    bvout[p]=bv;
    cout[p]=count[a];
    p=p+1;
    }
    vs[a]=v[a];
    nfs[a]=nf[a];
    }
}

    for(f=0;f<points;f++){
myfile << f*sampling*dt << " " << vout[f] << " "<< nout[f] << " "<< bvout[f]<<" " << cout[f] << "\n";
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



