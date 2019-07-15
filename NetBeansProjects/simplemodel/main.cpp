/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on March 27, 2019, 2:51 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    ofstream myfile;
    myfile.open("simplemodel5d2.txt");
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> n(0,1);
    
    double I,D,dt,vr,vs,vt,vpeak,k,b,a,c,d,vnew,vold,unew,uold;
    int p,N,T,points,sampling,j,f;
    I=0;
    D=2;
    vr=-60;        
    vs=-57;
    vt=-55;
    k=0.001;
    b=0.2;
    a=0.03;
    c=-70;
    d=-4;
    vpeak=35;
    T=10000;
    dt=0.1;
    N=T/dt;
    points=10000;
    double vout[points],uout[points];
    sampling=N/points;
    vold=0;
    uold=0;
    p=0;
    for(j=0;j<N;j++){
        vnew=vold+k*(vold-vr)*(vold-vt)*dt-uold*dt+I*dt+sqrt(2*D*dt)*n(gen);
        unew=uold+a*(b*(vold-vs)-uold)*dt;
        if(vnew>=vpeak){
            vnew=c;
            unew=d;
        }
        if(j==p*sampling){
            vout[p]=vnew;
            uout[p]=unew;
            p=p+1;
        }
        uold=unew;
        vold=vnew;
    }
    for(f=0;f<points;f++){
    myfile << f*sampling*dt << " " << vout[f]<< " "<< uout[f]<<"\n";
    }
    myfile.close();  
    return 0;
}

