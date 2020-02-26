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
  std::ofstream myfile2;
  std::ofstream myfile3;
    myfile.open("xtraje18ftfteneq2long.txt");
myfile2.open("xtraje18ftfteneq2longshort.txt");
myfile3.open("param18ftften2.txt");
const double pi = 3.14159265358979323846;    
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> n(0,1);
std::uniform_real_distribution<> dis(0, 2*pi);

int points=100000;
int runs=200;
int N=10000000;
int sampling=N/points;
double dt=0.001;
double dtint=sampling*dt;
double T=N*dt;


//int fpointscheck=2*sampling*10;
int fpointsshort,shortsamp;
double deltapeakshort,dtshort,Tfshort;
/*if(fpointscheck<points){
  fpointsshort=fpointscheck;
  shortsamp=1;
  dtshort=dt;
}
else{
  fpointsshort=points;
  shortsamp=floor(fpointscheck/points);
  dtshort=shortsamp*dt;
}
deltapeakshort=1/dtshort;
Tfshort=dtshort*fpointsshort;
*/

fpointsshort=100000;
shortsamp=1;
dtshort=dt;
deltapeakshort=1/dtshort;
Tfshort=dtshort*fpointsshort;
 
 
int foutpoints=1000;
int foutvalues[foutpoints];
foutvalues[0]=1;
double logstep=exp(log(points/2)/foutpoints);
int wints=round(1/(logstep-1));
int wnew,wold,wlength;
for(int w=1;w<wints;w++){
  foutvalues[w]=w+1;
}
for(int w=wints;w<foutpoints;w++){
  wold=foutvalues[w-1];
  wnew=round(wold*logstep);
  if(wnew>points/2){
    wlength=w;
    break;
  }
  else{
    foutvalues[w]=wnew;
  }
}
int foutfinal[wlength];
for (int w2=0;w2<wlength;w2++){
  foutfinal[w2]=foutvalues[w2];
  
} 

int foutvaluesshort[foutpoints];
foutvaluesshort[0]=1;
double logstepshort=exp(log(fpointsshort/2)/foutpoints);
int wintshort=round(1/(logstepshort-1));
int wnewshort,woldshort,wlengthshort;
for(int w=1;w<wintshort;w++){
  foutvaluesshort[w]=w+1;
}
for(int w=wintshort;w<foutpoints;w++){
  woldshort=foutvaluesshort[w-1];
  wnewshort=round(woldshort*logstepshort);
  if(wnewshort>fpointsshort/2){
    wlengthshort=w;
    break;
  }
  else{
    foutvaluesshort[w]=wnewshort;
  }
}

int foutfinalshort[wlengthshort];
for (int w2=0;w2<wlengthshort;w2++){
  foutfinalshort[w2]=foutvaluesshort[w2];
}

double tau=3;
double sigma=5;
double epsilon=0.01;
double omega=0.01;
double t,phi;
//double xv[points];
//double xvout[500][2];
double sqout[points];
initd(sqout,points);
double sqoutshort[fpointsshort];
initd(sqoutshort,fpointsshort);
double x;
int p,qshort,qvarshort;

//fftw_complex in[points], out[points];
fftw_complex* in = new fftw_complex[points];
fftw_complex* out = new fftw_complex[points];

fftw_complex* inshort = new fftw_complex[fpointsshort];
fftw_complex* outshort = new fftw_complex[fpointsshort];

for(int a=0;a<runs;a++){
    p=0;
    x=0;
    qshort=0;
    qvarshort=0;
    phi=dis(gen);
  for(int j2=0;j2<N;j2++){
    qvarshort=qvarshort+1;
      t=j2*dt;
    x=x*(1-dt/tau)+epsilon*cos(omega*t*pi+phi)*dt+sqrt((2*sigma*sigma*dt)/tau)*n(gen);
    if(j2==p*sampling){
    in[p][REAL]=x;
    in[p][IMAG]=0;
    p=p+1;
    }
    if(qshort<fpointsshort){
      if(qvarshort==shortsamp){
        qvarshort=0;
        inshort[qshort][REAL]=x;
        inshort[qshort][IMAG]=0;
        qshort=qshort+1;
      }
    }
  }
fft(in, out,points);  
for(int ind=0;ind<points;ind++){
    sqout[ind]=(a*sqout[ind]+(out[ind][REAL]*out[ind][REAL]+out[ind][IMAG]*out[ind][IMAG])*dtint*dtint)/(a+1);
}
fft(inshort, outshort,fpointsshort);  
for(int ind=0;ind<fpointsshort;ind++){
  sqoutshort[ind]=(a*sqoutshort[ind]+(outshort[ind][REAL]*outshort[ind][REAL]+outshort[ind][IMAG]*outshort[ind][IMAG])*dtshort*dtshort)/(a+1);
}
}

for(int rv=0;rv<wlength;rv++){
  myfile << foutvalues[rv]/T << " " << sqout[foutvalues[rv]] <<"\n";
}
myfile.close();
for(int rv=0;rv<wlengthshort;rv++){
  myfile2 << foutvaluesshort[rv]/Tfshort << " " << sqoutshort[foutvaluesshort[rv]] <<"\n";
}
myfile2.close();
myfile3<<T<<" "<<Tfshort<<" "<<points<<" "<<fpointsshort<<"\n";
myfile3.close();

return 0;
}
