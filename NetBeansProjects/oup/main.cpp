/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on May 23, 2019, 9:24 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <fftw3.h>
#include <math.h>
#include <complex>

#define REAL 0
#define IMAG 1

void fft(fftw_complex *in, fftw_complex *out, int N)
{
	// create a DFT plan
	fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	// execute the plan
	fftw_execute(plan);
	// do some cleaning
	fftw_destroy_plan(plan);
	fftw_cleanup();
}
void initd(double x[],int n){
    for(int k2=0;k2<n;k2++){
        x[k2]=0;
    }
}

int main(int argc, char** argv) {
std::ofstream myfile;
    myfile.open("xtraje2j.txt");

const double pi = 3.14159265358979323846;    
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> n(0,1);
std::uniform_real_distribution<> dis(0, 2*pi);

int points=1000000;
int runs=2000;
int N=100000000;
int sampling=N/points;
double dt=0.001;
double dtint=sampling*dt;
double T=N*dt;
double tau=1;
double sigma=1;
double epsilon=0.01;
double omega=0.01;
double t,phi;
//double xv[points];
//double xvout[500][2];
double sqout[points];
initd(sqout,points);
double x;
int p;

//fftw_complex in[points], out[points];
fftw_complex* in = new fftw_complex[points];
fftw_complex* out = new fftw_complex[points];

for(int a=0;a<runs;a++){
    p=0;
    x=0;
    phi=dis(gen);
  for(int j2=0;j2<N;j2++){
      t=j2*dt;
    x=x*(1-dt/tau)+epsilon*cos(omega*t*pi+phi)*dt+sqrt((2*sigma*sigma*dt)/tau)*n(gen);
    if(j2==p*sampling){
    in[p][REAL]=x;
    in[p][IMAG]=0;
    p=p+1;
    }
  }
fft(in, out,points);  
for(int ind=0;ind<points;ind++){
    sqout[ind]=(a*sqout[ind]+(out[ind][REAL]*out[ind][REAL]+out[ind][IMAG]*out[ind][IMAG])*dtint*dtint)/(a+1);
}  
}

for(int rv=0;rv<points;rv++){
  myfile << rv/T << " " << sqout[rv] <<"\n";
}
myfile.close();
return 0;
}
