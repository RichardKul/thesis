/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on May 26, 2019, 6:06 PM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <complex>
#include <iostream>
#include <valarray>
 
const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}
 
// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}

int main(int argc, char** argv) {
std::ofstream myfile;
    myfile.open("xtraj.txt");
    
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> n(0,1);

int runs=10;
int N=100000;
int points=100000;
int sampling=N*runs/points;
int j22;
double dt=0.001;
double T0=N*runs*dt;
double tau=1;
double sigma=1;
Complex xv[points];
double xvout[points];
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

CArray data(xv, points);
 
    // forward fft
    fft(data);
    
    for(int b=0;b<points;b++){
        xvout[b]=std::abs(data[b])*std::abs(data[b]);
    }
for(int rv=0;rv<points;rv++){
  myfile << rv/T0 << " " << xvout[rv] <<"\n";
}
myfile.close();
return 0;
}

