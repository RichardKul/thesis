/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <fftw3.h>

using namespace std;

int main(int argc, char** argv) {
    ofstream myfile;
    myfile.open("xtraj.txt");
    
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> n(0,1);

const int d=1;
const int * zeiger;
zeiger=&d;
int runs=10;
int N=1000000;
int points=1000000;
int sampling=N*runs/points;
int j22;
double dt=0.001;
double tau=1;
double sigma=1;
double xv[points];
double xvout[points][2];
int p=0;
double x=0;

for(int a=0;a<runs;a++){
  for(int j2=0;j2<N;j2++){
    j22=j2+a*N;
    x=x*(1-dt/tau)+sqrt((2*sigma*sigma*dt)/tau)*n(gen);
    if(j22==p*sampling){
    xv[p]=x;
    p=p+1;
    }
  }
}

fftw_plan q;
q=fftw_plan_dft_r2c(N,zeiger,xv,xvout,FFTW_ESTIMATE);
fftw_execute(q);
fftw_destroy_plan(q); 
for(int rv=0;rv<points;rv++){
  myfile << rv*sampling*dt << " " << xvout[rv] <<"\n";
}
myfile.close();
return 0;
}